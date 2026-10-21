#!/usr/bin/env python3
"""
Build an nf-core/mag samplesheet (sample,run,group,short_reads_1,short_reads_2,long_reads)
from a folder of paired FASTQ files and a metadata table.

metadata.tsv (tab-separated, header required):
    run_accession   sample   group
    ERR1234567      dog01_B  cohortA
    ERR1234568      dog01_B  cohortA      <- second run of the same sample (gets run=1)
    ...

FASTQ names expected: <run_accession>_1.fastq.gz / <run_accession>_2.fastq.gz
(ENA / fasterq-dump naming); _R1/_R2 variants are also recognised.

Usage:
    python3 make_samplesheet.py --fastq-dir raw/ --metadata metadata.tsv --out samplesheet.csv
"""
import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

R1_PATTERNS = ["{r}_1.fastq.gz", "{r}_R1.fastq.gz", "{r}_R1_001.fastq.gz", "{r}_1.fq.gz", "{r}_R1.fq.gz"]


def find_pair(fastq_dir: Path, run: str):
    for pat in R1_PATTERNS:
        r1 = fastq_dir / pat.format(r=run)
        if r1.exists():
            r2 = Path(str(r1).replace("_R1", "_R2") if "_R1" in r1.name else re.sub(r"_1(\.f(ast)?q\.gz)$", r"_2\1", str(r1)))
            if r2.exists():
                return r1.resolve(), r2.resolve()
            sys.exit(f"ERROR: found {r1.name} but not its mate {r2.name}")
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fastq-dir", required=True, type=Path)
    ap.add_argument("--metadata", required=True, type=Path)
    ap.add_argument("--out", default="samplesheet.csv", type=Path)
    a = ap.parse_args()

    with open(a.metadata, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    need = {"run_accession", "sample", "group"}
    if not rows or not need.issubset(rows[0].keys()):
        sys.exit(f"ERROR: metadata must have columns {sorted(need)}")

    run_idx = defaultdict(int)
    out, missing = [], []
    for row in rows:
        run, sample, group = row["run_accession"].strip(), row["sample"].strip(), row["group"].strip()
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", sample):
            sys.exit(f"ERROR: sample name '{sample}' has characters other than letters, digits, _ . -")
        pair = find_pair(a.fastq_dir, run)
        if pair is None:
            missing.append(run)
            continue
        out.append([sample, run_idx[sample], group, str(pair[0]), str(pair[1]), ""])
        run_idx[sample] += 1

    if missing:
        sys.exit(f"ERROR: no FASTQ pair found for {len(missing)} runs, e.g. {missing[:5]}")

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "run", "group", "short_reads_1", "short_reads_2", "long_reads"])
        w.writerows(out)

    groups = defaultdict(set)
    for s, _, g, *_ in out:
        groups[g].add(s)
    print(f"Wrote {a.out}: {len(out)} runs, {len(run_idx)} samples, {len(groups)} groups")
    for g, s in sorted(groups.items()):
        print(f"  group {g}: {len(s)} samples")
    big = [g for g, s in groups.items() if len(s) > 20]
    if big:
        print(f"WARNING: groups with >20 samples ({', '.join(big)}) make binning_map_mode=group "
              f"costly (samples^2 mappings). Consider splitting them.")


if __name__ == "__main__":
    main()
