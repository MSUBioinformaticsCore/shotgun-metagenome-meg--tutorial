#!/bin/bash
# =============================================================================
#  Step 1 of 3: check everything and build the sample list. Submits nothing.
#  Run on a dev node, from the folder with the scripts:
#     bash 1_prepare.sh
# =============================================================================

set -euo pipefail
cd "$(dirname "$0")"
source ./func_config.sh

echo "== Software"
require_tools python3 emapper.py
for t in python3 emapper.py; do echo "   $t: $(command -v "$t")"; done
emapper.py --data_dir "$EGGNOG_DB" --version 2>&1 | sed 's/\x1b\[[0-9;]*m//g'

echo "== eggNOG database: $EGGNOG_DB"
for f in eggnog.db eggnog_proteins.dmnd; do
    [[ -s $EGGNOG_DB/$f ]] || { echo "ERROR: $EGGNOG_DB/$f missing"; exit 1; }
    echo "   $f: $(du -h "$EGGNOG_DB/$f" | cut -f1)"
done

echo "== Inputs from $MAG_OUTDIR"
[[ -s $SAMPLESHEET ]] || { echo "ERROR: samplesheet not found: $SAMPLESHEET"; exit 1; }
python3 - --mag-outdir "$MAG_OUTDIR" --samplesheet "$SAMPLESHEET" --out func_manifest.tsv <<'PY'
#!/usr/bin/env python3
"""
Build the input manifest for functional annotation (option B: contig depths).
One row per sample, all files taken from the nf-core/mag output:

  sample   faa                         gff                         depth                                   fastp_json
           Annotation/Prodigal/...     Annotation/Prodigal/...     GenomeBinning/depths/contigs/           QC_shortreads/fastp/...
           *_prodigal.faa.gz           *_prodigal.gff.gz           MEGAHIT-<sample>-depth.txt.gz           (comma-separated if several runs)

Files are matched to samples by name: the sample name must appear in the file
name with no letter or digit directly before or after it, and the longest
matching sample name wins (so 'B1' never claims 'B11' files).

Usage:
  (run by 1_prepare.sh)
"""
import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


def assign(name, samples):
    hits = [s for s in samples if re.search(rf"(?<![A-Za-z0-9]){re.escape(s)}(?![A-Za-z0-9])", name)]
    return max(hits, key=len) if hits else None


def collect(root, samples, keep):
    out = defaultdict(list)
    if root.exists():
        for p in sorted(root.rglob("*")):
            if p.is_file() and keep(p):
                s = assign(p.name, samples)
                if s:
                    out[s].append(p.resolve())
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mag-outdir", required=True, type=Path)
    ap.add_argument("--samplesheet", required=True, type=Path)
    ap.add_argument("--out", default="func_manifest.tsv", type=Path)
    a = ap.parse_args()

    with open(a.samplesheet, newline="") as fh:
        samples = sorted({r["sample"] for r in csv.DictReader(fh)})
    o = a.mag_outdir

    prodigal = lambda p: "prodigal" in str(p).lower()
    faa = collect(o / "Annotation", samples, lambda p: prodigal(p) and re.search(r"\.faa(\.gz)?$", p.name))
    gff = collect(o / "Annotation", samples, lambda p: prodigal(p) and re.search(r"\.gff(\.gz)?$", p.name))
    depth = collect(o / "GenomeBinning", samples,
                    lambda p: "depths" in p.parts and "contigs" in p.parts and p.name.endswith("-depth.txt.gz"))
    fastp = collect(o / "QC_shortreads", samples,
                    lambda p: "fastp" in str(p.relative_to(o)).lower() and p.name.endswith(".json"))

    rows, problems = [], []
    for s in samples:
        f, g, d, j = faa.get(s, []), gff.get(s, []), depth.get(s, []), fastp.get(s, [])
        for label, files in (("Prodigal .faa", f), ("Prodigal .gff", g), ("contig depth file", d)):
            if len(files) != 1:
                problems.append(f"{s}: {len(files)} {label} files")
        if not j:
            problems.append(f"{s}: no fastp .json (needed for mean read length)")
        rows.append([s,
                     f[0] if len(f) == 1 else "NA",
                     g[0] if len(g) == 1 else "NA",
                     d[0] if len(d) == 1 else "NA",
                     ",".join(map(str, j)) or "NA"])

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "faa", "gff", "depth", "fastp_json"])
        w.writerows(rows)

    print(f"Wrote {a.out}: {len(rows)} samples")
    if problems:
        print("PROBLEMS (fix before submitting):")
        for p in problems:
            print("  " + p)
        sys.exit(1)
    print(f"OK: {len(rows)} samples ready.")


if __name__ == "__main__":
    main()
PY

mkdir -p logs "$OUT"
N=$(( $(wc -l < func_manifest.tsv) - 1 ))

cat <<EOF

== Ready ($N samples). Next:
   sbatch --array=1 2_annotate.slurm          # test sample 1 first
   sbatch --array=1-$N%20 2_annotate.slurm      # then all (finished samples are skipped)
   sbatch 3_merge.slurm                       # after all annotations are done
EOF
