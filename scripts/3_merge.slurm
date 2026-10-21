#!/bin/bash
#SBATCH --job-name=func_merge
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%j.log

# =============================================================================
#  Step 3 of 3: sample x function tables (gene depth from mag x eggNOG function).
#  Submit from the script folder once all 2_annotate.slurm tasks are done:
#     sbatch 3_merge.slurm
#  Stops without merging if any sample is not annotated yet.
# =============================================================================

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ./func_config.sh
require_tools python3

echo "== Checking annotations"
missing=()
row=0
while IFS=$'\t' read -r SAMPLE _rest; do
    row=$((row + 1))
    [[ -s $OUT/$SAMPLE/$SAMPLE.emapper.annotations ]] || missing+=("$row:$SAMPLE")
done < <(tail -n +2 func_manifest.tsv)
if (( ${#missing[@]} )); then
    echo "ERROR: not annotated yet (task:sample): ${missing[*]}"
    echo "Rerun those tasks, e.g.  sbatch --array=$(printf '%s\n' "${missing[@]}" | cut -d: -f1 | paste -sd,) 2_annotate.slurm"
    exit 1
fi
echo "   all $row samples annotated"

echo "== Merging"
python3 - --func-dir "$OUT" --manifest func_manifest.tsv --out-dir "$OUT/tables" <<'PY'
#!/usr/bin/env python3
"""
Sample x function abundance tables from eggNOG-mapper annotations and
nf-core/mag contig depths. Python 3 standard library only (no pandas).

Per sample
  genes      mag Prodigal GFF: gene ID = <contig>_<gene no.> (matches the protein FASTA
             header and therefore the eggNOG-mapper 'query'); gene length = end - start + 1
  depth      mag contig depth file GenomeBinning/depths/contigs/<asm>-<sample>-depth.txt.gz
             (MetaBAT2 jgi_summarize_bam_contig_depths). Only the column of the sample
             itself is used: '<asm>-<sample>-<sample>.bam' = the sample's reads on its own
             assembly. Each gene takes the depth of the contig it lies on.
  function   eggNOG-mapper .emapper.annotations

Per category (KO, KEGG_Pathway [ko IDs only], KEGG_Module, EC, COG [COG IDs from eggNOG_OGs], CAZy, PFAMs)
a gene with n terms gives 1/n of its value to each term; values are summed over genes.

Three tables per category:
  depth_<cat>.tsv      summed gene depth (coverage units)
  dpm_<cat>.tsv        depth per million: gene depth / sum of depth over ALL genes of the
                       sample (annotated or not) x 1e6; proportional to TPM
  estreads_<cat>.tsv   estimated reads: gene depth x gene length / mean read length
                       (mean read length = fastp after-filtering total_bases / total_reads);
                       an approximation of read counts, not observed counts

annotation_summary.tsv per sample: genes, genes without depth, genes annotated, mean read
length, and % of total gene depth carried by annotated genes in each category.

Usage:
  (run by 3_merge.slurm)
"""
import argparse
import csv
import gzip
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

CATEGORIES = {
    "KO": "KEGG_ko",
    "KEGG_Pathway": "KEGG_Pathway",
    "KEGG_Module": "KEGG_Module",
    "EC": "EC",
    "COG": "eggNOG_OGs",          # COG IDs (e.g. COG0513) from the eggNOG orthologous groups
    "CAZy": "CAZy",
    "PFAMs": "PFAMs",
}
VALUES = ("depth", "dpm", "estreads")
GFF_ID = re.compile(r"(?:^|;)ID=(\d+)_(\d+)(?:;|$)")


class SampleError(Exception):
    pass


def opener(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def split_terms(value, category):
    if not value or value == "-":
        return []
    if category == "COG":                                   # "COG0513@1|root,COG0513@2|Bacteria" -> COG0513
        return list(dict.fromkeys(re.findall(r"(?:^|,)\s*(COG\d{4})@", value)))
    terms = [t.strip() for t in value.split(",") if t.strip()]
    if category == "KO":
        terms = [t.replace("ko:", "") for t in terms]
    if category == "KEGG_Pathway":                          # keep ko00010, drop map00010 duplicates
        terms = [t for t in terms if t.startswith("ko")]
    return terms


def read_genes(gff):
    """gene ID -> (contig, length)"""
    genes = {}
    with opener(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            m = GFF_ID.search(f[8])
            if not m:
                raise SampleError(f"{gff}: CDS without Prodigal ID=<n>_<m>: {f[8][:60]}")
            gid = f"{f[0]}_{m.group(2)}"
            if gid in genes:
                raise SampleError(f"{gff}: duplicated gene ID {gid}")
            genes[gid] = (f[0], int(f[4]) - int(f[3]) + 1)
    return genes


def read_own_depth(path, sample):
    """Contig depth of `sample` on its own assembly: column '<asm>-<sample>-<sample>.bam'."""
    prefix = Path(path).name[: -len("-depth.txt.gz")]          # e.g. MEGAHIT-B1
    col = f"{prefix}-{sample}.bam"
    depth = {}
    with opener(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if col not in header:
            bams = [h for h in header if h.endswith(".bam")]
            raise SampleError(f"{Path(path).name}: no column '{col}'; columns are {bams[:5]}")
        j = header.index(col)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            depth[f[0]] = float(f[j])
    return depth


def mean_read_length(json_list):
    bases = reads = 0
    for j in json_list.split(","):
        with open(j) as fh:
            af = json.load(fh)["summary"]["after_filtering"]
        bases += af["total_bases"]
        reads += af["total_reads"]
    if reads == 0:
        raise SampleError(f"zero reads after filtering in {json_list}")
    return bases / reads


def read_emapper(path):
    """eggNOG-mapper .annotations -> {query: {column: value}}"""
    cols, ann = None, {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#query"):
                cols = line.lstrip("#").rstrip("\n").split("\t")
                continue
            if line.startswith("#") or not line.strip():
                continue
            if cols is None:
                raise SampleError(f"no '#query' header line before data in {path}")
            f = line.rstrip("\n").split("\t")
            ann[f[0]] = dict(zip(cols, f))
    if cols is None:
        raise SampleError(f"no '#query' header line in {path}")
    missing = [c for c in CATEGORIES.values() if c not in cols]
    if missing:
        raise SampleError(f"{path}: missing columns {missing}")
    return ann


def process_sample(row, func_dir):
    s = row["sample"]
    ann_f = func_dir / s / f"{s}.emapper.annotations"
    if not ann_f.exists():
        raise SampleError(f"missing {ann_f}")

    genes = read_genes(row["gff"])
    depth = read_own_depth(row["depth"], s)
    rl = mean_read_length(row["fastp_json"])

    gdepth, glen = {}, {}
    n_nodepth = 0
    for g, (contig, length) in genes.items():
        if contig in depth:
            gdepth[g] = depth[contig]
            glen[g] = length
        else:
            n_nodepth += 1
    total = sum(gdepth.values())

    ann = read_emapper(ann_f)
    shared = [q for q in ann if q in gdepth]
    if not shared:
        raise SampleError("no gene IDs shared between eggNOG-mapper output and the GFF")

    summary = {"sample": s, "genes_in_gff": len(genes), "genes_without_depth": n_nodepth,
               "genes_annotated": len(shared),
               "annotated_ids_not_in_gff": sum(1 for q in ann if q not in genes),
               "mean_read_length": round(rl, 2)}
    per_cat = {}
    for cat, col in CATEGORIES.items():
        agg = defaultdict(lambda: [0.0, 0.0, 0.0])      # term -> [depth, dpm, estreads]
        cat_depth = 0.0
        for q in shared:
            terms = split_terms(ann[q].get(col, "-"), cat)
            if not terms:
                continue
            d = gdepth[q]
            cat_depth += d
            w = 1 / len(terms)
            vals = (d * w, (d / total * 1e6 if total else 0.0) * w, d * glen[q] / rl * w)
            for t in terms:
                a = agg[t]
                a[0] += vals[0]
                a[1] += vals[1]
                a[2] += vals[2]
        per_cat[cat] = agg
        summary[f"pct_depth_{cat}"] = round(100 * cat_depth / total, 2) if total else 0
    return per_cat, summary


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--func-dir", required=True, type=Path, help="OUT folder of func_annot_array.slurm")
    ap.add_argument("--manifest", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    a = ap.parse_args()

    with open(a.manifest, newline="") as fh:
        manifest = list(csv.DictReader(fh, delimiter="\t"))
    a.out_dir.mkdir(parents=True, exist_ok=True)

    tables = {cat: {} for cat in CATEGORIES}            # cat -> sample -> {term: [d, dpm, est]}
    summaries, failed = [], []
    for row in manifest:
        try:
            per_cat, summ = process_sample(row, a.func_dir)
        except (SampleError, OSError, KeyError, ValueError) as e:
            failed.append(f"{row['sample']}: {e}")
            continue
        summaries.append(summ)
        for cat, agg in per_cat.items():
            tables[cat][row["sample"]] = agg
        print(f"  {summ['sample']}: {summ['genes_in_gff']} genes, {summ['genes_without_depth']} without depth, "
              f"{summ['genes_annotated']} annotated, {summ['pct_depth_KO']}% of depth on KO genes")

    if failed:
        print("Samples skipped:\n  " + "\n  ".join(failed))
    if not summaries:
        sys.exit("ERROR: no sample could be processed")

    done = [x["sample"] for x in summaries]
    for cat in CATEGORIES:
        terms = sorted(set().union(*(tables[cat][s].keys() for s in done)))
        for i, v in enumerate(VALUES):
            with open(a.out_dir / f"{v}_{cat}.tsv", "w") as fh:
                fh.write(cat + "\t" + "\t".join(done) + "\n")
                for t in terms:
                    fh.write(t + "\t" + "\t".join(
                        str(round(tables[cat][s][t][i], 4)) if t in tables[cat][s] else "0.0"
                        for s in done) + "\n")
        print(f"{cat}: {len(terms)} features")

    with open(a.out_dir / "annotation_summary.tsv", "w") as fh:
        cols = list(summaries[0].keys())
        fh.write("\t".join(cols) + "\n")
        for summ in summaries:
            fh.write("\t".join(str(summ[c]) for c in cols) + "\n")
    print(f"Tables written to {a.out_dir}")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
PY

echo "== Checks (annotation_summary.tsv)"
awk -F'\t' 'NR==1 { for (i=1;i<=NF;i++) c[$i]=i; next }
    { n++; nd += $c["genes_without_depth"]; ni += $c["annotated_ids_not_in_gff"] }
    END { printf "   samples: %d\n   genes without depth (should be 0): %d\n   annotated IDs not in GFF (should be 0): %d\n", n, nd, ni }' \
    "$OUT/tables/annotation_summary.tsv"
echo "== Done: tables in $OUT/tables"
