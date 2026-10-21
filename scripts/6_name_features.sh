#!/bin/bash
# =============================================================================
#  Step 6: add the official name to every functional ID in the bacteria-only tables.
#  Rows whose ID has no exact current name in the official source, and any
#  excluded (non-bacterial) term, are removed and listed in unmapped_features.tsv.
#  Run on a dev node (downloads name files once), from the script folder:
#     bash 6_name_features.sh                    (FORCE=1 ... re-downloads)
#     TABLES=tables bash 6_name_features.sh      (name the unfiltered tables instead)
#
#  Name sources:
#    KO            KEGG rest.kegg.jp/list/ko
#    KEGG_Pathway  KEGG rest.kegg.jp/list/pathway/ko
#    KEGG_Module   KEGG rest.kegg.jp/list/module
#    EC            ExPASy ENZYME enzyme.dat (IUBMB recommended name); enzclass.txt for x.x.x.-
#                  transferred / deleted EC numbers are removed
#    COG           NCBI COG 2024 cog-24.def.tab (COG name); COG IDs not in COG 2024 are removed
#    CAZy          CAZy class + family number (Glycoside Hydrolase family 13, ...)
#    PFAMs         Pfam-A.clans.tsv family description: current release, then archived
#                  releases 37.0 ... 31.0 (eggNOG may use names from an older release)
#
#  Output: $OUT/<TABLES>_named/  depth_/dpm_/estreads_<category>.tsv with a 'name' column,
#          feature_names.tsv (id -> name -> source), unmapped_features.tsv (removed rows + reason),
#          naming_coverage.tsv (per table: rows in = rows named + rows removed)
# =============================================================================

set -euo pipefail
cd "$(dirname "$0")"
source ./func_config.sh
REF_DIR=${REF_DIR:-$OUT/reference}
TABLES=${TABLES:-tables_bacteria}
require_tools python3 curl
mkdir -p "$REF_DIR"
[[ -d $OUT/$TABLES ]] || { echo "ERROR: $OUT/$TABLES not found (run 5_filter.slurm first)"; exit 1; }

get () {   # get <url> <file>
    if [[ -s $2 && -z ${FORCE:-} ]]; then echo "   have $(basename "$2")"; return; fi
    echo "   downloading $(basename "$2")"
    curl -fsSL --retry 3 -o "$2.tmp" "$1" && mv "$2.tmp" "$2" || { rm -f "$2.tmp"; return 1; }
}

echo "== Name files -> $REF_DIR"
get https://rest.kegg.jp/list/ko                 "$REF_DIR/kegg_list_ko.tsv"
get https://rest.kegg.jp/list/pathway/ko         "$REF_DIR/kegg_list_pathway_ko.tsv"
get https://rest.kegg.jp/list/module             "$REF_DIR/kegg_list_module.tsv"
get https://ftp.expasy.org/databases/enzyme/enzyme.dat   "$REF_DIR/enzyme.dat"
get https://ftp.expasy.org/databases/enzyme/enzclass.txt "$REF_DIR/enzclass.txt"
get https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/cog-24.def.tab "$REF_DIR/cog-24.def.tab"
get https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz "$REF_DIR/Pfam-A.clans.tsv.gz"
# Archived Pfam releases: eggNOG may report domain names from an older release.
# A release that cannot be downloaded is skipped (its names are then simply not available).
for v in 37.0 36.0 35.0 34.0 33.1 32.0 31.0; do
    get "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam$v/Pfam-A.clans.tsv.gz" \
        "$REF_DIR/Pfam-A.$v.clans.tsv.gz" || echo "   (Pfam $v not available; skipped)"
done
echo "names downloaded: $(date)" >> "$REF_DIR/reference_versions.txt"

echo "== Naming $OUT/$TABLES"
python3 - --tables "$OUT/$TABLES" --ref "$REF_DIR" --out "$OUT/${TABLES}_named" <<'PY'
#!/usr/bin/env python3
"""
Add the official name of every functional ID to the bacteria-only tables.
A row is kept only if its ID has an exact current name in the official source and is not
an excluded (non-bacterial) term; every other row is removed and listed in
unmapped_features.tsv with the reason.

Name sources (exact text from the source, nothing edited except a trailing period):
  KO            KEGG  rest.kegg.jp/list/ko                 (KEGG definition: symbols; name [EC])
  KEGG_Pathway  KEGG  rest.kegg.jp/list/pathway/ko
  KEGG_Module   KEGG  rest.kegg.jp/list/module
  EC            ExPASy ENZYME enzyme.dat DE line (IUBMB recommended name); partial ECs
                (e.g. 1.1.1.-) from ExPASy enzclass.txt. Transferred / deleted ECs -> removed.
                ECs whose official name matches the photosynthesis / eukaryote / organelle
                name rules of 4_build_exclusions.sh -> removed.
  COG           NCBI COG 2024 cog-24.def.tab (column 3: COG name), COG IDs taken by
                3_merge / 5_filter from eggNOG_OGs. COG IDs not in COG 2024 -> removed.
  CAZy          CAZy class name + family number (GH, GT, PL, CE, AA, CBM; cazy.org class names)
  PFAMs         Pfam-A.clans.tsv (column 5: family description; eggNOG writes Pfam IDs, e.g.
                'Toprim'). Current release first, then archived releases
                (Pfam-A.<version>.clans.tsv.gz), because eggNOG may use names from an older release.
Every row of every input table is either named or removed (checked; naming_coverage.tsv).
"""
import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

# Same name rules as 4_build_exclusions.sh, applied here to the official EC names
PHOTO = re.compile(r"photosynth|photosystem|chlorophyll|light-harvesting|phycobili|phycocyan|"
                   r"phycoerythr|reaction cent(?:er|re)|thylakoid|cyanobacteri|"
                   r"magnesium chelatase|mg-chelatase|mg-protoporphyrin|magnesium-protoporphyrin", re.I)
EUK = re.compile(r"eukaryot|\bplants?\b|\byeasts?\b|\banimals?\b|\bfly\b|\binsects?\b|\bfung(?:i|al|us)\b|"
                 r"mitochondri|chloroplast|plastid|nucleosome|nuclear pore|spliceos|meiosis|oocyte", re.I)
KEEP = re.compile(r"\bprokaryot\w*|\bbacteria\b|\bbacterial\b", re.I)


def name_rule(text):
    if PHOTO.search(text):
        return "name: photosynthesis"
    if EUK.search(text) and not KEEP.search(text):
        return "name: eukaryote/organelle"
    return None


CATS = ("KO", "KEGG_Pathway", "KEGG_Module", "EC", "COG", "CAZy", "PFAMs")
CAZY_CLASS = {"GH": "Glycoside Hydrolase", "GT": "GlycosylTransferase", "PL": "Polysaccharide Lyase",
              "CE": "Carbohydrate Esterase", "AA": "Auxiliary Activity", "CBM": "Carbohydrate-Binding Module"}


def opener(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def read_kegg_list(path, prefix):
    """KEGG REST list: '<id>\\t<name>' with optional 'ko:'/'path:'/'md:' prefix."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if "\t" not in line:
                continue
            i, name = line.rstrip("\n").split("\t", 1)
            i = i.split(":", 1)[-1]
            if re.match(prefix, i):
                out[i] = name.strip()
    return out


def read_enzyme(dat, enzclass):
    """EC -> (name or None, reason)"""
    out = {}
    ec, de = None, []
    with open(dat) as fh:
        for line in fh:
            tag, val = line[:2], line[5:].rstrip("\n")
            if tag == "ID":
                ec, de = val.strip(), []
            elif tag == "DE":
                de.append(val.strip())
            elif tag == "//" and ec:
                name = " ".join(de).strip()
                if name.startswith("Transferred entry"):
                    out[ec] = (None, "EC transferred (" + name.rstrip(".") + ")")
                elif name.startswith("Deleted entry"):
                    out[ec] = (None, "EC deleted")
                else:
                    out[ec] = (name[:-1] if name.endswith(".") else name, "")
                ec = None
    n_class = 0
    if enzclass and Path(enzclass).exists():
        with open(enzclass, errors="replace") as fh:
            for line in fh:
                m = re.match(r"^\s*(\d+)\.\s*([\d-]+)\.\s*([\d-]+)\.\s*([\d-]+)\s+(.+?)\s*$", line)
                if m:
                    key = ".".join(m.group(i) for i in range(1, 5))
                    name = m.group(5)
                    out.setdefault(key, (name[:-1] if name.endswith(".") else name, ""))
                    n_class += 1
    return out, n_class


def read_cog(path):
    """cog-24.def.tab: COG ID, category letters, COG name, gene, pathway, PubMed, PDB"""
    raw = open(path, "rb").read()
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        text = raw.decode("cp1252")
    out = {}
    for f in csv.reader(text.splitlines(), delimiter="\t", quoting=csv.QUOTE_NONE):
        if len(f) >= 3 and re.fullmatch(r"COG\d{4}", f[0]) and f[2].strip():
            out[f[0]] = f[2].strip()
    return out


def read_pfam(paths):
    """Pfam ID -> (description, release). Earlier files (newer releases) take precedence."""
    out = {}
    for path, release in paths:
        with opener(path) as fh:
            for f in csv.reader(fh, delimiter="\t"):
                if len(f) >= 5 and f[0].startswith("PF") and f[3] not in out:
                    out[f[3]] = (f[4].strip(), release)
    return out


def cazy_name(i):
    m = re.fullmatch(r"(GH|GT|PL|CE|AA|CBM)(\d+)(?:_(\d+))?", i)
    if not m:
        return None
    name = f"{CAZY_CLASS[m.group(1)]} family {int(m.group(2))}"
    return name + (f", subfamily {int(m.group(3))}" if m.group(3) else "")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tables", required=True, type=Path, help="tables_bacteria folder")
    ap.add_argument("--ref", required=True, type=Path, help="REF_DIR with downloaded name files")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--no-size-check", action="store_true", help=argparse.SUPPRESS)  # tests only
    a = ap.parse_args()
    r = a.ref

    enz, n_class = read_enzyme(r / "enzyme.dat", r / "enzclass.txt")
    names = {
        "KO": read_kegg_list(r / "kegg_list_ko.tsv", r"^K\d{5}$"),
        "KEGG_Pathway": read_kegg_list(r / "kegg_list_pathway_ko.tsv", r"^ko\d{5}$"),
        "KEGG_Module": read_kegg_list(r / "kegg_list_module.tsv", r"^M\d{5}$"),
        "COG": read_cog(r / "cog-24.def.tab"),
    }
    pfam_files = [(r / "Pfam-A.clans.tsv.gz", "Pfam current release")]
    pfam_files += sorted(((p, "Pfam " + p.name.split(".", 1)[1].rsplit(".clans", 1)[0])
                          for p in r.glob("Pfam-A.*.clans.tsv.gz")),
                         key=lambda x: [int(n) for n in re.findall(r"\d+", x[1])], reverse=True)
    pfam = read_pfam(pfam_files)
    names["PFAMs"] = {k: v[0] for k, v in pfam.items()}
    sizes = {"KO": 20000, "KEGG_Pathway": 300, "KEGG_Module": 300, "COG": 4000, "PFAMs": 10000}
    for k, n in sizes.items():
        if len(names[k]) < n and not a.no_size_check:
            sys.exit(f"ERROR: name source for {k} looks incomplete ({len(names[k])} entries)")
    if len(enz) < 5000 and not a.no_size_check:
        sys.exit(f"ERROR: enzyme.dat looks incomplete ({len(enz)} entries)")
    source = {"KO": "KEGG list/ko", "KEGG_Pathway": "KEGG list/pathway/ko", "KEGG_Module": "KEGG list/module",
              "EC": "ExPASy ENZYME", "COG": "NCBI COG 2024 cog-24.def.tab",
              "CAZy": "CAZy class + family", "PFAMs": "Pfam-A.clans.tsv"}

    excluded = defaultdict(set)
    ef = a.tables / "excluded_features.tsv"
    if ef.exists():
        with open(ef, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                excluded[row["category"]].add(row["id"])
    else:
        print(f"NOTE: {ef} not found; excluded-term check skipped")

    def lookup(cat, i):
        if i in excluded[cat]:
            return None, "excluded non-bacterial term"
        if cat == "EC":
            n, why = enz.get(i, (None, "EC not in ExPASy ENZYME"))
            if n and name_rule(n):
                return None, f"official EC name matches {name_rule(n)} ({n})"
            return n, why
        if cat == "CAZy":
            n = cazy_name(i)
            return n, "" if n else "not a CAZy family ID (GH/GT/PL/CE/AA/CBM + number)"
        n = names[cat].get(i)
        if cat == "PFAMs":
            return n, "" if n else "ID not in any downloaded Pfam release"
        return n, "" if n else f"ID not in {source[cat]}"

    a.out.mkdir(parents=True, exist_ok=True)
    unmapped, named, coverage = {}, {}, []
    total_depth = defaultdict(dict)
    for cat in CATS:
        for v in ("depth", "dpm", "estreads"):
            src = a.tables / f"{v}_{cat}.tsv"
            if not src.exists():
                sys.exit(f"ERROR: missing {src}")
            with open(src) as fh, open(a.out / f"{v}_{cat}.tsv", "w") as out:
                header = fh.readline().rstrip("\n").split("\t")
                out.write("\t".join([header[0], "name"] + header[1:]) + "\n")
                n_in = n_named = n_removed = 0
                for line in fh:
                    n_in += 1
                    f = line.rstrip("\n").split("\t")
                    i = f[0]
                    name, why = lookup(cat, i)
                    if v == "depth":
                        total_depth[cat][i] = sum(float(x) for x in f[1:])
                    if name is None:
                        unmapped[(cat, i)] = why
                        n_removed += 1
                        continue
                    named[(cat, i)] = name
                    n_named += 1
                    out.write("\t".join([i, name.replace("\t", " ")] + f[1:]) + "\n")
                if n_named + n_removed != n_in:
                    sys.exit(f"ERROR: {src.name}: {n_in} rows read but {n_named} named + {n_removed} removed")
                coverage.append((v, cat, n_in, n_named, n_removed))

    with open(a.out / "feature_names.tsv", "w") as fh:
        fh.write("category\tid\tname\tsource\n")
        for (cat, i), n in sorted(named.items()):
            src = pfam[i][1] if cat == "PFAMs" else source[cat]
            fh.write(f"{cat}\t{i}\t{n}\t{src}\n")
    with open(a.out / "unmapped_features.tsv", "w") as fh:
        fh.write("category\tid\treason\ttotal_depth_all_samples\n")
        for (cat, i), why in sorted(unmapped.items()):
            fh.write(f"{cat}\t{i}\t{why}\t{round(total_depth[cat].get(i, 0), 4)}\n")

    with open(a.out / "naming_coverage.tsv", "w") as fh:
        fh.write("table\tcategory\trows_in\trows_named\trows_removed\n")
        for v, cat, n_in, n_named, n_removed in coverage:
            fh.write(f"{v}_{cat}.tsv\t{cat}\t{n_in}\t{n_named}\t{n_removed}\n")
    print("Pfam releases used: " + ", ".join(rel for _, rel in pfam_files))
    print(f"Name sources: KO {len(names['KO'])}, pathways {len(names['KEGG_Pathway'])}, "
          f"modules {len(names['KEGG_Module'])}, EC {len(enz)} (incl. {n_class} classes), "
          f"COG {len(names['COG'])}, Pfam {len(names['PFAMs'])}")
    bad_excl = sum(1 for why in unmapped.values() if why == "excluded non-bacterial term")
    for cat in CATS:
        ids = total_depth[cat]
        tot = sum(ids.values())
        lost = [i for i in ids if (cat, i) in unmapped]
        pct = 100 * sum(ids[i] for i in lost) / tot if tot else 0
        print(f"  {cat}: {len(ids) - len(lost)} named, {len(lost)} removed ({pct:.2f}% of depth)")
    if bad_excl:
        print(f"WARNING: {bad_excl} excluded terms were present in the input tables and were removed")
    print(f"Wrote {a.out}")


if __name__ == "__main__":
    main()
PY

echo
echo "== Check what was removed:  less -S $OUT/${TABLES}_named/unmapped_features.tsv"
