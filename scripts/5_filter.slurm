#!/bin/bash
#SBATCH --job-name=func_bacteria
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%j.log

# =============================================================================
#  Step 5: bacteria-only functional abundance tables.
#  Needs all samples annotated (2_annotate.slurm) and excluded_features.tsv
#  (4_build_exclusions.sh). Submit from the script folder:
#     sbatch 5_filter.slurm
#
#  Layer 1  gene origin   (gene removed from all categories)
#     tiara     contig = eukarya / archaea / organelle (mitochondrion, plastid) in mag's Tiara output
#     origin    eggNOG lineage must include Bacteria and must not include Archaea, Eukaryota
#               (protists, fungi, plants, animals), Viruses, Cyanobacteria or Chlorobi;
#               genes without an eggNOG hit have no known origin and are removed
#  Layer 2  gene function (gene removed from all categories)
#     function  all KOs (or, without a KO, all Pfams) are photosynthesis / eukaryote / organelle terms
#  Layer 3  terms of kept genes
#     terms     every term in excluded_features.tsv is dropped in KO, KEGG pathway, KEGG module,
#               EC, COG (COG IDs) and Pfam before the abundance is split among the gene's terms
#
#  Output: $OUT/tables_bacteria/  depth_/dpm_/estreads_<category>.tsv, annotation_summary.tsv,
#          filter_summary.tsv, excluded_genes.tsv.gz. $OUT/tables is not changed.
# =============================================================================

RULES=tiara,origin,function,terms      # remove a word to switch that rule off

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ./func_config.sh
REF_DIR=${REF_DIR:-$OUT/reference}
require_tools python3

EXCL=$REF_DIR/excluded_features.tsv
[[ -s $EXCL ]] || { echo "ERROR: $EXCL missing; run 'bash 4_build_exclusions.sh' first"; exit 1; }

missing=()
while IFS=$'\t' read -r SAMPLE _rest; do
    [[ -s $OUT/$SAMPLE/$SAMPLE.emapper.annotations ]] || missing+=("$SAMPLE")
done < <(tail -n +2 func_manifest.tsv)
(( ${#missing[@]} == 0 )) || { echo "ERROR: not annotated yet: ${missing[*]}"; exit 1; }

python3 - --func-dir "$OUT" --manifest func_manifest.tsv --mag-outdir "$MAG_OUTDIR" \
          --exclusions "$EXCL" --out-dir "$OUT/tables_bacteria" --rules "$RULES" <<'PY'
#!/usr/bin/env python3
"""
Bacteria-only sample x function tables. Same abundance method as 3_merge.slurm
(gene depth = depth of its contig in the sample's own-assembly column of mag's depth file;
a gene with n terms in a category gives 1/n of its value to each term).

Layer 1 - gene ORIGIN (gene removed from every category):
  tiara     contig labelled eukarya, archaea or organelle (mitochondrion / plastid) in mag's
            per-contig Tiara output (contigs below mag's tiara_min_length, default 3000 bp,
            are not classified and are kept)
  origin    eggNOG ortholog lineage (eggNOG_OGs, 'OG@taxid|name', broadest -> narrowest) must
            include Bacteria (2) and must not include Archaea (2157), Eukaryota (2759: protists,
            fungi, plants, animals), Viruses (10239), Cyanobacteria (1117 or name
            'Cyanobacteri*') or Chlorobi (1090 or name 'Chlorobi*', green sulfur bacteria).
            Genes without an eggNOG assignment have no known origin and are removed.
Layer 2 - gene FUNCTION (gene removed from every category):
  function  all of the gene's KOs are photosynthesis / eukaryote / organelle KOs by name, or
            (no KO) all of its Pfams are, per excluded_features.tsv (removes_gene = yes / if_no_KO)
Layer 3 - TERMS (kept genes): every term listed in excluded_features.tsv is dropped from the
  gene's term list before the 1/n split, in KO, KEGG_Pathway, KEGG_Module, EC, COG
  and PFAMs. A kept gene left with no term in a category does not count in that category.

DPM is computed over the kept genes only.

Outputs (--out-dir): depth_/dpm_/estreads_<category>.tsv, annotation_summary.tsv,
filter_summary.tsv (per sample and rule: genes and % of depth removed),
excluded_genes.tsv.gz (each removed gene and why).
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
NONBACTERIAL_TAXIDS = {"2157": "Archaea", "2759": "Eukaryota", "10239": "Viruses",
                      "1117": "Cyanobacteria", "1090": "Chlorobi"}
NONBACTERIAL_NAMES = re.compile(r"^(cyanobacteri|chlorobi)", re.I)
TIARA_FST = {"eukarya", "archaea", "organelle"}
TIARA_SND = {"mitochondrion", "plastid"}


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
    if category == "KEGG_Pathway":
        terms = [t for t in terms if t.startswith("ko")]
    return terms


def read_genes(gff):
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
    prefix = Path(path).name[: -len("-depth.txt.gz")]
    col = f"{prefix}-{sample}.bam"
    depth = {}
    with opener(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if col not in header:
            raise SampleError(f"{Path(path).name}: no column '{col}'")
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
    need = list(CATEGORIES.values()) + ["eggNOG_OGs"]
    missing = [c for c in need if c not in cols]
    if missing:
        raise SampleError(f"{path}: missing columns {missing}")
    return ann


def origin_problem(value):
    """None if the eggNOG lineage is bacterial and not an excluded lineage, else a reason."""
    levels = []
    for og in (value or "").split(","):
        m = re.search(r"@(\d+)\|(.+)$", og.strip())
        if m:
            levels.append((m.group(1), m.group(2).strip()))
    if not levels:
        return "origin:no eggNOG lineage"
    for taxid, name in levels:
        if taxid in NONBACTERIAL_TAXIDS:
            return f"origin:{NONBACTERIAL_TAXIDS[taxid]}"
        if NONBACTERIAL_NAMES.match(name):
            return f"origin:{name}"
    if "2" not in {t for t, _ in levels}:
        return "origin:no Bacteria level"
    return None


def read_exclusions(path):
    """excluded_features.tsv -> ({category: set(ids)}, KOs that remove genes, Pfams that remove genes if no KO)"""
    terms, ko_gene, pf_gene = defaultdict(set), set(), set()
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            terms[r["category"]].add(r["id"])
            if r["category"] == "KO" and r["removes_gene"] == "yes":
                ko_gene.add(r["id"])
            if r["category"] == "PFAMs" and r["removes_gene"] == "if_no_KO":
                pf_gene.add(r["id"])
    if not terms:
        raise SystemExit(f"ERROR: no rows in {path}")
    return terms, ko_gene, pf_gene


def assign(name, samples):
    hits = [s for s in samples if re.search(rf"(?<![A-Za-z0-9]){re.escape(s)}(?![A-Za-z0-9])", name)]
    return max(hits, key=len) if hits else None


def find_tiara(mag_outdir, samples):
    """sample -> path of mag's per-contig Tiara table (header has sequence_id, class_fst_stage)"""
    found = defaultdict(list)
    root = Path(mag_outdir) / "Taxonomy"
    for p in sorted(root.rglob("*")) if root.exists() else []:
        if not p.is_file() or not re.search(r"\.(txt|tsv)(\.gz)?$", p.name):
            continue
        try:
            with opener(p) as fh:
                header = fh.readline().rstrip("\n").split("\t")
        except (OSError, UnicodeDecodeError, EOFError):
            continue
        if {"sequence_id", "class_fst_stage"} <= set(header):
            s = assign(p.name, samples)
            if s:
                found[s].append(p)
    return found


def read_tiara_excluded(path):
    """contig -> Tiara label, for contigs labelled eukarya/archaea/organelle (or mito/plastid)"""
    excl = {}
    with opener(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            fst = (r.get("class_fst_stage") or "").lower()
            snd = (r.get("class_snd_stage") or "").lower()
            if fst in TIARA_FST or snd in TIARA_SND:
                excl[r["sequence_id"].split()[0]] = snd if snd in TIARA_SND else fst
    return excl


def process_sample(row, func_dir, tiara_file, rules, excl):
    excl_terms, ko_gene, pf_gene = excl
    s = row["sample"]
    ann_f = func_dir / s / f"{s}.emapper.annotations"
    if not ann_f.exists():
        raise SampleError(f"missing {ann_f}")

    genes = read_genes(row["gff"])
    depth = read_own_depth(row["depth"], s)
    rl = mean_read_length(row["fastp_json"])
    ann = read_emapper(ann_f)
    tiara_excl = read_tiara_excluded(tiara_file) if ("tiara" in rules and tiara_file) else {}

    gdepth, glen = {}, {}
    n_nodepth = 0
    for g, (contig, length) in genes.items():
        if contig in depth:
            gdepth[g] = depth[contig]
            glen[g] = length
        else:
            n_nodepth += 1

    # ---- layers 1 and 2: genes removed ----
    removed = {}
    for g in gdepth:
        reasons = []
        contig = genes[g][0]
        if contig in tiara_excl:
            reasons.append(f"tiara:{tiara_excl[contig]}")
        a = ann.get(g)
        if "origin" in rules:
            why = origin_problem(a.get("eggNOG_OGs")) if a else "origin:no eggNOG annotation"
            if why:
                reasons.append(why)
        if "function" in rules and a is not None:
            kos = split_terms(a.get("KEGG_ko"), "KO")
            pfs = split_terms(a.get("PFAMs"), "PFAMs")
            if kos and all(k in ko_gene for k in kos):
                reasons.append("function:KO " + "+".join(kos))
            elif not kos and pfs and all(p in pf_gene for p in pfs):
                reasons.append("function:Pfam " + "+".join(pfs))
        if reasons:
            removed[g] = reasons

    total_all = sum(gdepth.values())
    kept = {g: d for g, d in gdepth.items() if g not in removed}
    total = sum(kept.values())

    fsum = {"sample": s, "genes_with_depth": len(gdepth), "genes_removed": len(removed),
            "genes_kept": len(kept),
            "pct_depth_removed": round(100 * (total_all - total) / total_all, 3) if total_all else 0,
            "tiara_file": str(tiara_file) if tiara_file else "none"}
    for rule in ("tiara", "origin", "function"):
        gs = [g for g, rs in removed.items() if any(r.startswith(rule + ":") for r in rs)]
        fsum[f"genes_{rule}"] = len(gs)
        fsum[f"pct_depth_{rule}"] = round(100 * sum(gdepth[g] for g in gs) / total_all, 3) if total_all else 0

    # ---- layer 3: terms ----
    shared = [q for q in ann if q in kept]
    summary = {"sample": s, "genes_in_gff": len(genes), "genes_without_depth": n_nodepth,
               "genes_removed_by_filter": len(removed), "genes_kept": len(kept),
               "genes_annotated_kept": len(shared),
               "annotated_ids_not_in_gff": sum(1 for q in ann if q not in genes),
               "mean_read_length": round(rl, 2)}
    per_cat = {}
    for cat, col in CATEGORIES.items():
        agg = defaultdict(lambda: [0.0, 0.0, 0.0])
        cat_depth = dropped_depth = 0.0
        drop = excl_terms.get(cat, set()) if "terms" in rules else set()
        for q in shared:
            terms0 = split_terms(ann[q].get(col, "-"), cat)
            terms = [t for t in terms0 if t not in drop]
            if terms0 and not terms:
                dropped_depth += kept[q]
            if not terms:
                continue
            d = kept[q]
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
        fsum[f"pct_kept_depth_all_terms_dropped_{cat}"] = round(100 * dropped_depth / total, 3) if total else 0
    excluded_rows = [(s, g, genes[g][0], round(gdepth[g], 4), ";".join(rs)) for g, rs in removed.items()]
    return per_cat, summary, fsum, excluded_rows


def write_tsv(path, rows):
    with open(path, "w") as fh:
        cols = list(rows[0].keys())
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--func-dir", required=True, type=Path)
    ap.add_argument("--manifest", required=True, type=Path)
    ap.add_argument("--mag-outdir", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--exclusions", required=True, type=Path, help="excluded_features.tsv from 4_build_exclusions.sh")
    ap.add_argument("--rules", default="tiara,origin,function,terms",
                    help="comma list of: tiara, origin, function, terms")
    a = ap.parse_args()
    rules = {r.strip() for r in a.rules.split(",") if r.strip()}
    bad = rules - {"tiara", "origin", "function", "terms"}
    if bad:
        sys.exit(f"ERROR: unknown rules {bad}")

    with open(a.manifest, newline="") as fh:
        manifest = list(csv.DictReader(fh, delimiter="\t"))
    samples = [r["sample"] for r in manifest]

    tiara = {}
    if "tiara" in rules:
        found = find_tiara(a.mag_outdir, samples)
        multi = {s: p for s, p in found.items() if len(p) > 1}
        if multi:
            sys.exit(f"ERROR: several Tiara tables for one sample, e.g. {list(multi.items())[:1]}")
        tiara = {s: p[0] for s, p in found.items()}
        missing = [s for s in samples if s not in tiara]
        if missing:
            sys.exit(f"ERROR: no Tiara per-contig table under {a.mag_outdir}/Taxonomy for "
                     f"{len(missing)} samples, e.g. {missing[:3]}. Use RULES without 'tiara' to skip this rule.")
    exclusions = read_exclusions(a.exclusions)
    print(f"Rules: {', '.join(sorted(rules))}")
    print("Excluded terms loaded: " + ", ".join(f"{k} {len(v)}" for k, v in sorted(exclusions[0].items())))

    a.out_dir.mkdir(parents=True, exist_ok=True)
    tables = {cat: {} for cat in CATEGORIES}
    summaries, fsums, failed = [], [], []
    with gzip.open(a.out_dir / "excluded_genes.tsv.gz", "wt") as ex:
        ex.write("sample\tgene\tcontig\tdepth\treason\n")
        for row in manifest:
            try:
                per_cat, summ, fsum, excl = process_sample(row, a.func_dir, tiara.get(row["sample"]), rules, exclusions)
            except (SampleError, OSError, KeyError, ValueError) as e:
                failed.append(f"{row['sample']}: {e}")
                continue
            for r in excl:
                ex.write("\t".join(map(str, r)) + "\n")
            summaries.append(summ)
            fsums.append(fsum)
            for cat, agg in per_cat.items():
                tables[cat][row["sample"]] = agg
            print(f"  {fsum['sample']}: removed {fsum['genes_removed']} of {fsum['genes_with_depth']} genes "
                  f"({fsum['pct_depth_removed']}% of depth; tiara {fsum['genes_tiara']}, "
                  f"origin {fsum['genes_origin']}, function {fsum['genes_function']})")

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
    write_tsv(a.out_dir / "annotation_summary.tsv", summaries)
    write_tsv(a.out_dir / "filter_summary.tsv", fsums)
    print(f"Tables written to {a.out_dir}")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
PY

cp "$EXCL" "$OUT/tables_bacteria/"
[[ -s $REF_DIR/reference_versions.txt ]] && cp "$REF_DIR/reference_versions.txt" "$OUT/tables_bacteria/"
echo "== Summary (filter_summary.tsv; a gene can match several rules)"
awk -F'\t' 'NR==1 { for (i=1;i<=NF;i++) c[$i]=i; next }
    { n++; t += $c["pct_depth_removed"]; ti += $c["pct_depth_tiara"]; o += $c["pct_depth_origin"]; f += $c["pct_depth_function"] }
    END { printf "   samples: %d\n   mean %% of depth removed: %.2f (tiara %.2f, origin %.2f, function %.2f)\n", n, t/n, ti/n, o/n, f/n }' \
    "$OUT/tables_bacteria/filter_summary.tsv"
echo "== Done: $OUT/tables_bacteria"
