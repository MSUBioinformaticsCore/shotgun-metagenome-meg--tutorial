#!/bin/bash
# =============================================================================
#  Step 4: download KEGG, Pfam and COG reference files and build the list of
#  non-bacterial functional terms (eukaryote, archaea, mitochondria, chloroplast,
#  photosynthesis) for every functional category. Nothing is filtered here.
#  Run on a dev node (needs internet), from the script folder:
#     bash 4_build_exclusions.sh            (FORCE=1 bash 4_build_exclusions.sh re-downloads)
#
#  Output in $REF_DIR:
#     kegg_ko00001.keg, kegg_ko00002.keg, Pfam-A.clans.tsv.gz,
#     cog-24.def.tab                                            reference files
#     reference_versions.txt                                    KEGG / Pfam release info
#     excluded_features.tsv   category, id, name, rule, removes_gene  -> REVIEW THIS FILE
# =============================================================================

set -euo pipefail
cd "$(dirname "$0")"
source ./func_config.sh
REF_DIR=${REF_DIR:-$OUT/reference}
require_tools python3 curl
mkdir -p "$REF_DIR"

get () {   # get <url> <file>
    if [[ -s $2 && -z ${FORCE:-} ]]; then echo "   have $(basename "$2")"; return; fi
    echo "   downloading $(basename "$2")"
    curl -fsSL --retry 3 -o "$2.tmp" "$1" && mv "$2.tmp" "$2"
}

echo "== Reference files -> $REF_DIR"
get https://rest.kegg.jp/get/br:ko00001 "$REF_DIR/kegg_ko00001.keg"
get https://rest.kegg.jp/get/br:ko00002 "$REF_DIR/kegg_ko00002.keg"
get https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz "$REF_DIR/Pfam-A.clans.tsv.gz"
get https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/cog-24.def.tab "$REF_DIR/cog-24.def.tab"
{
    echo "downloaded: $(date)"
    echo "KEGG:"; curl -fsSL https://rest.kegg.jp/info/kegg | head -3
    echo "COG: NCBI COG 2024 (cog-24.def.tab)"
    echo "Pfam:"; curl -fsSL https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/relnotes.txt | grep -m2 -i "release"
} > "$REF_DIR/reference_versions.txt" 2>/dev/null || true

echo "== Building excluded_features.tsv"
python3 - --ko00001 "$REF_DIR/kegg_ko00001.keg" --ko00002 "$REF_DIR/kegg_ko00002.keg" \
          --pfam-clans "$REF_DIR/Pfam-A.clans.tsv.gz" \
          --cog-def "$REF_DIR/cog-24.def.tab" --out "$REF_DIR/excluded_features.tsv" <<'PY'
#!/usr/bin/env python3
"""
Build the list of functional terms that are NOT bacterial functions, for every
functional category, from KEGG and Pfam reference files. Output is a reviewable table
(excluded_features.tsv); nothing is filtered here.

Matching rules (applied to names / definitions):
  PHOTO  photosynthesis-related, always excluded (also in bacteria, e.g. purple/green
         bacteria, cyanobacteria):
         photosynth | photosystem | chlorophyll (incl. bacteriochlorophyll) | light-harvesting |
         phycobili/phycocyan/phycoerythr | reaction center | thylakoid | cyanobacteri |
         magnesium chelatase | Mg-protoporphyrin (chlorophyll-branch enzymes)
         Not photosynthesis-specific, therefore kept: RuBisCO / Calvin cycle (also in
         chemoautotrophs), carotenoids, F-type ATP synthase.
  EUK    eukaryote / organelle related, excluded UNLESS the same name also says prokaryotes
         or bacteria (e.g. M00157 "F-type ATPase, prokaryotes and chloroplasts" is kept):
         eukaryot | plant(s) | yeast | animal(s) | fly | insect | fungi/fungal | mitochondri |
         chloroplast | plastid | nucleosome | nuclear pore | spliceos | meiosis | oocyte
  KEEP   prokaryote(s) | bacteria | bacterial  (whole words; 'cyanobacteria' does not count)

Per category
  KEGG_Pathway  KEGG BRITE ko00001 classes: all of Organismal Systems and Human Diseases;
                subclasses Cellular community - eukaryotes, Transport and catabolism,
                Signaling molecules and interaction, Chromosome, Information processing in
                viruses; Signal transduction except Two-component system; Cell growth and death
                except Cell cycle - Caulobacter; Cell motility except Bacterial chemotaxis and
                Flagellar assembly; name rules; plus EXTRA_PATHWAYS below.
  KEGG_Module   BRITE ko00002 'Photosynthesis' group, name rules, C4/CAM modules.
  KO            KO definition matches the name rules  -> term removed AND the gene removed
                (if all of the gene's KOs are such KOs).
                KO whose KEGG maps (overview maps ko011xx/ko012xx ignored) are ALL excluded
                maps -> term removed, gene kept (its other functions stay).
  EC            EC numbers carried only by name-rule-excluded KOs.
  COG           NCBI COG 2024 cog-24.def.tab (term removed, gene kept):
                - COG 2024 pathway (column 5) is Photosystem I, Photosystem II, Archaeal ribosomal
                  proteins or Archaella;
                - COG name matches the name rules, except the COGs in COG_NAME_KEEP (bacterial COGs
                  whose name mentions a eukaryote, organelle or Mg-chelatase only as a homolog,
                  an alternative function or a protein architecture);
                - functional category letters include Y (Nuclear structure).
                No name rule for 'archaeal': many COG names are shared bacterial/archaeal families
                (e.g. COG0472 MraY, COG0636 ATP synthase subunit c); archaeal genes are removed by
                the gene-origin rule of 5_filter.slurm.
  PFAMs         Pfam ID or description matches the name rules -> term removed; gene removed
                if it has no KO and all its Pfams are excluded.
  CAZy          no term-level rule (CAZy families are not domain-specific); non-bacterial
                CAZymes are removed by the gene-origin rules in 5_filter.slurm.
"""
import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

PHOTO = re.compile(r"photosynth|photosystem|chlorophyll|light-harvesting|phycobili|phycocyan|"
                   r"phycoerythr|reaction cent(?:er|re)|thylakoid|cyanobacteri|"
                   r"magnesium chelatase|mg-chelatase|mg-protoporphyrin|magnesium-protoporphyrin", re.I)
EUK = re.compile(r"eukaryot|\bplants?\b|\byeasts?\b|\banimals?\b|\bfly\b|\binsects?\b|\bfung(?:i|al|us)\b|"
                 r"mitochondri|chloroplast|plastid|nucleosome|nuclear pore|spliceos|meiosis|oocyte", re.I)
KEEP = re.compile(r"\bprokaryot\w*|\bbacteria\b|\bbacterial\b", re.I)

EXCLUDED_A = {"Organismal Systems", "Human Diseases"}
EXCLUDED_B = {"Cellular community - eukaryotes", "Transport and catabolism",
              "Signaling molecules and interaction", "Chromosome", "Information processing in viruses"}
KEEP_IN_B = {"Signal transduction": {"Two-component system"},
             "Cell growth and death": {"Cell cycle - Caulobacter"},
             "Cell motility": {"Bacterial chemotaxis", "Flagellar assembly"}}
# Eukaryote/plant/fungal-specific maps whose names carry no keyword (names as in KEGG PATHWAY)
EXTRA_PATHWAYS = {
    "ko03040": "Spliceosome", "ko03013": "Nucleocytoplasmic transport",
    "ko03015": "mRNA surveillance pathway", "ko04141": "Protein processing in endoplasmic reticulum",
    "ko04130": "SNARE interactions in vesicular transport", "ko04120": "Ubiquitin mediated proteolysis",
    "ko03022": "Basal transcription factors", "ko03460": "Fanconi anemia pathway",
    "ko00073": "Cutin, suberine and wax biosynthesis", "ko00100": "Steroid biosynthesis",
    "ko00120": "Primary bile acid biosynthesis", "ko00140": "Steroid hormone biosynthesis",
    "ko00062": "Fatty acid elongation", "ko00590": "Arachidonic acid metabolism",
    "ko00592": "alpha-Linolenic acid metabolism", "ko00510": "N-Glycan biosynthesis",
    "ko00513": "Various types of N-glycan biosynthesis", "ko00512": "Mucin type O-glycan biosynthesis",
    "ko00515": "Mannose type O-glycan biosynthesis", "ko00514": "Other types of O-glycan biosynthesis",
    "ko00532": "Glycosaminoglycan biosynthesis - chondroitin sulfate / dermatan sulfate",
    "ko00534": "Glycosaminoglycan biosynthesis - heparan sulfate / heparin",
    "ko00533": "Glycosaminoglycan biosynthesis - keratan sulfate",
    "ko00563": "Glycosylphosphatidylinositol (GPI)-anchor biosynthesis",
    "ko00601": "Glycosphingolipid biosynthesis - lacto and neolacto series",
    "ko00603": "Glycosphingolipid biosynthesis - globo and isoglobo series",
    "ko00604": "Glycosphingolipid biosynthesis - ganglio series", "ko00830": "Retinol metabolism",
    "ko00905": "Brassinosteroid biosynthesis", "ko00908": "Zeatin biosynthesis",
    "ko00940": "Phenylpropanoid biosynthesis", "ko00945": "Stilbenoid, diarylheptanoid and gingerol biosynthesis",
    "ko00941": "Flavonoid biosynthesis", "ko00944": "Flavone and flavonol biosynthesis",
    "ko00942": "Anthocyanin biosynthesis", "ko00943": "Isoflavonoid biosynthesis",
    "ko00901": "Indole alkaloid biosynthesis", "ko00403": "Indole diterpene alkaloid biosynthesis",
    "ko00950": "Isoquinoline alkaloid biosynthesis", "ko00960": "Tropane, piperidine and pyridine alkaloid biosynthesis",
    "ko00232": "Caffeine metabolism", "ko00965": "Betalain biosynthesis", "ko00966": "Glucosinolate biosynthesis",
    "ko00402": "Benzoxazinoid biosynthesis", "ko00254": "Aflatoxin biosynthesis",
}
EXTRA_MODULES = {"M00168": "CAM (Crassulacean acid metabolism), dark",
                 "M00169": "CAM (Crassulacean acid metabolism), light",
                 "M00170": "C4-dicarboxylic acid cycle, phosphoenolpyruvate carboxykinase type",
                 "M00171": "C4-dicarboxylic acid cycle, NAD - malic enzyme type",
                 "M00172": "C4-dicarboxylic acid cycle, NADP - malic enzyme type"}
COG_PATHWAYS = {"Photosystem I", "Photosystem II", "Archaeal ribosomal proteins", "Archaella"}
COG_NAME_KEEP = {   # name matches a rule word, but the COG is a bacterial function (reason)
    "COG0270": "Dcm DNA-cytosine methylase; eukaryotic tRNA methylase named as alternative",
    "COG0293": "RlmE/FtsJ 23S rRNA methylase; eukaryotic tRNA methylase named as alternative",
    "COG0398": "YdjX membrane protein; fungal transporter named as a relative",
    "COG0606": "YifB, 'predicted Mg-chelatase' ATPase family, not chlorophyll synthesis",
    "COG1399": "YceD 23S rRNA accumulation protein; 'essential in plants' is a remark",
    "COG1429": "CobN, COG 2024 pathway Cobalamine/B12 biosynthesis",
    "COG3476": "TspO/CrtK sensory protein; mitochondrial receptor named as homolog",
    "COG4981": "'yeast-type' fatty acid synthase describes the protein architecture",
    "COG4982": "'yeast-type' fatty acid synthase describes the protein architecture",
    "COG5034": "'plant homeodomain' (PHD) is a domain name",
    "COG1467": "'eukaryotic-type' DNA primase names a protein family (AEP)",
    "COG2219": "'eukaryotic-type' DNA primase names a protein family (AEP)",
    "COG6116": "'eukaryotic-type' DNA primase names a protein family (AEP)",
    "COG6117": "'eukaryotic-type' DNA primase names a protein family (AEP)",
    "COG6118": "'eukaryotic-type' DNA primase names a protein family (AEP), LigD-like domain",
}
OVERVIEW = re.compile(r"^ko01[12]\d\d$")


def name_rule(text):
    """'photo' / 'eukaryote-organelle' / None"""
    if PHOTO.search(text):
        return "name: photosynthesis"
    if EUK.search(text) and not KEEP.search(text):
        return "name: eukaryote/organelle"
    return None


def read_cog_def(path):
    """cog-24.def.tab rows: COG ID, category letters, COG name, gene, pathway, PubMed, PDB"""
    raw = open(path, "rb").read()
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        text = raw.decode("cp1252")
    for f in csv.reader(text.splitlines(), delimiter="\t", quoting=csv.QUOTE_NONE):
        if len(f) >= 3 and re.fullmatch(r"COG\d{4}", f[0]):
            yield f


def parse_brite(path):
    """KEGG BRITE text -> list of (A, B, C, D) lines with level letter"""
    rows = []
    A = B = C = ""
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line[0] not in "ABCD":
                continue
            lvl, text = line[0], line[1:].strip()
            if lvl == "A":
                A, B, C = re.sub(r"^\d+\s+", "", text), "", ""
            elif lvl == "B":
                B, C = re.sub(r"^\d+\s+", "", text), ""
            elif lvl == "C":
                C = text
            rows.append((lvl, A, B, C, text))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ko00001", required=True, type=Path)
    ap.add_argument("--ko00002", required=True, type=Path)
    ap.add_argument("--pfam-clans", required=True, type=Path)
    ap.add_argument("--cog-def", required=True, type=Path, help="NCBI COG 2024 cog-24.def.tab")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--no-size-check", action="store_true", help=argparse.SUPPRESS)  # tests only
    a = ap.parse_args()
    out = []          # (category, id, name, rule, removes_gene)

    # ---------------- KEGG pathways + KO definitions ----------------
    pw_info, ko_def, ko_pw, ko_ec = {}, {}, defaultdict(set), defaultdict(set)
    for lvl, A, B, C, text in parse_brite(a.ko00001):
        if lvl == "C":
            m = re.match(r"(\d{5})\s+(.*?)\s*\[PATH:(ko\d{5})\]", text)
            if m:
                pw_info[m.group(3)] = (m.group(2), A, B)
        elif lvl == "D":
            m = re.match(r"(K\d{5})\s+(.*)$", text)
            if not m:
                continue
            ko, d = m.group(1), m.group(2)
            ko_def.setdefault(ko, d)
            for ec in re.findall(r"\[EC:([^\]]+)\]", d):
                ko_ec[ko].update(ec.split())
            pm = re.search(r"\[PATH:(ko\d{5})\]", C)
            if pm:
                ko_pw[ko].add(pm.group(1))
    if not a.no_size_check and (len(pw_info) < 300 or len(ko_def) < 10000):
        sys.exit(f"ERROR: {a.ko00001} looks incomplete ({len(pw_info)} maps, {len(ko_def)} KOs)")

    excl_pw = {}
    for pid, (name, A, B) in pw_info.items():
        rule = None
        if A in EXCLUDED_A:
            rule = f"KEGG class: {A}"
        elif B in EXCLUDED_B:
            rule = f"KEGG subclass: {B}"
        elif B in KEEP_IN_B and name not in KEEP_IN_B[B]:
            rule = f"KEGG subclass: {B} (non-bacterial)"
        elif pid in EXTRA_PATHWAYS:
            rule = "eukaryote/plant/fungal-specific map"
        else:
            rule = name_rule(name)
        if rule:
            excl_pw[pid] = (name, rule)
    missing_extra = [p for p in EXTRA_PATHWAYS if p not in pw_info]
    for pid, (name, rule) in sorted(excl_pw.items()):
        out.append(("KEGG_Pathway", pid, name, rule, "no"))

    # ---------------- KEGG modules ----------------
    n_mod = 0
    for lvl, A, B, C, text in parse_brite(a.ko00002):
        if lvl != "D":
            continue
        m = re.match(r"(M\d{5})\s+(.*?)(?:\s*\[PATH:.*\])?$", text)
        if not m:
            continue
        n_mod += 1
        mid, name = m.group(1), m.group(2)
        rule = ("KEGG module group: Photosynthesis" if C == "Photosynthesis"
                else "C4/CAM plant carbon fixation" if mid in EXTRA_MODULES else name_rule(name))
        if rule:
            out.append(("KEGG_Module", mid, name, rule, "no"))
    if not a.no_size_check and n_mod < 300:
        sys.exit(f"ERROR: {a.ko00002} looks incomplete ({n_mod} modules)")

    # ---------------- KOs ----------------
    ko_named = set()
    for ko, d in sorted(ko_def.items()):
        rule = name_rule(d)
        if rule:
            ko_named.add(ko)
            out.append(("KO", ko, d, rule, "yes"))
            continue
        pws = {p for p in ko_pw.get(ko, ()) if not OVERVIEW.match(p)}
        if pws and pws <= set(excl_pw):
            out.append(("KO", ko, d, "all KEGG maps of this KO are excluded maps", "no"))

    # ---------------- EC ----------------
    ec_kos = defaultdict(set)
    for ko, ecs in ko_ec.items():
        for ec in ecs:
            ec_kos[ec].add(ko)
    for ec, kos in sorted(ec_kos.items()):
        if kos <= ko_named:
            out.append(("EC", ec, "carried by: " + ",".join(sorted(kos)),
                        "all KOs with this EC are excluded by name", "no"))

    # ---------------- COG ----------------
    n_cog = 0
    for f in read_cog_def(a.cog_def):
        n_cog += 1
        cog, letters, name = f[0], f[1], f[2]
        pathway = f[4].strip() if len(f) > 4 else ""
        rule = f"COG 2024 pathway: {pathway}" if pathway in COG_PATHWAYS else None
        if rule is None and cog not in COG_NAME_KEEP:
            rule = name_rule(name)
        if rule is None and "Y" in letters:
            rule = "COG category Y (Nuclear structure)"
        if rule:
            out.append(("COG", cog, name, rule, "no"))
    if not a.no_size_check and n_cog < 4000:
        sys.exit(f"ERROR: {a.cog_def} looks incomplete ({n_cog} COGs)")

    # ---------------- Pfam ----------------
    n_pf = 0
    with gzip.open(a.pfam_clans, "rt") as fh:
        for f in csv.reader(fh, delimiter="\t"):
            if len(f) < 5 or not f[0].startswith("PF"):
                continue
            n_pf += 1
            pid, desc = f[3], f[4]
            rule = name_rule(f"{pid} {desc}")
            if rule:
                out.append(("PFAMs", pid, desc, rule, "if_no_KO"))
    if not a.no_size_check and n_pf < 10000:
        sys.exit(f"ERROR: {a.pfam_clans} looks incomplete ({n_pf} families)")

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["category", "id", "name", "rule", "removes_gene"])
        w.writerows(out)

    counts = defaultdict(int)
    for r in out:
        counts[r[0]] += 1
    print(f"Reference: {len(pw_info)} KEGG maps, {n_mod} modules, {len(ko_def)} KOs, "
          f"{len(ec_kos)} EC numbers, {n_cog} COGs, {n_pf} Pfam families")
    for cat in ("KEGG_Pathway", "KEGG_Module", "KO", "EC", "COG", "PFAMs"):
        print(f"  excluded {cat}: {counts[cat]}")
    print("  CAZy: no term-level rule (handled by gene origin)")
    if missing_extra:
        print(f"NOTE: {len(missing_extra)} listed maps not in this KEGG release: {missing_extra}")
    print(f"Wrote {a.out}")


if __name__ == "__main__":
    main()
PY

echo
echo "== Review before filtering:  less -S $REF_DIR/excluded_features.tsv"
echo "   (columns: category, id, name, rule, removes_gene)"
echo "== Next:  sbatch 5_filter.slurm"
