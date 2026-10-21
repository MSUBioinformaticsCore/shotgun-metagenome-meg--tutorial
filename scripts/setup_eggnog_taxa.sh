#!/bin/bash
# =============================================================================
#  One-time setup: add the eggNOG-mapper NCBI taxonomy files (needed by
#  --target_taxa / --excluded_taxa in 2_annotate.slurm) next to the shared eggNOG
#  database, in a folder you can write to.
#
#  Creates $NEW_DB with:
#     eggnog.db, eggnog_proteins.dmnd      links to the shared database (not copied)
#     eggnog.taxa.db, eggnog.taxa.db.traverse.pkl   unpacked from eggnog.taxa.tar(.gz) in the
#                                  shared folder, or downloaded (emapperdb-5.0.2) if it has none
#  Run on a dev node, from the script folder:
#     bash setup_eggnog_taxa.sh
#  It then sets EGGNOG_DB=$NEW_DB in func_config.sh (old copy kept as func_config.sh.bak).
#
#  Download host: eggnog5.embl.de. The old host in eggNOG-mapper's
#  download_eggnog_data.py (eggnogdb.embl.de) returns 404; it is tried second.
# =============================================================================

set -euo pipefail
cd "$(dirname "$0")"
source ./func_config.sh

SHARED_DB=/mnt/research/bioinformaticsCore/shared/eggnog_db   # shared eggNOG database (read only)
NEW_DB=/mnt/scratch/$USER/eggnog_data             # writable folder (edit if you like)
URLS=(http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz
      http://eggnogdb.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz)

mkdir -p "$NEW_DB"
for f in eggnog.db eggnog_proteins.dmnd; do
    src=$(readlink -f "$SHARED_DB/$f")
    [[ -s $src ]] || { echo "ERROR: $SHARED_DB/$f not found"; exit 1; }
    [[ $src == "$(readlink -f "$NEW_DB/$f" 2>/dev/null)" ]] || ln -sfn "$src" "$NEW_DB/$f"
    echo "   linked $f -> $src"
done

if [[ -s $NEW_DB/eggnog.taxa.db && -s $NEW_DB/eggnog.taxa.db.traverse.pkl ]]; then
    echo "   taxonomy files already present"
else
    archive=
    for a in "$SHARED_DB/eggnog.taxa.tar" "$SHARED_DB/eggnog.taxa.tar.gz"; do
        [[ -s $a ]] && { archive=$a; echo "   using $a"; break; }
    done
    if [[ -z $archive ]]; then
        archive=$NEW_DB/eggnog.taxa.tar.gz
        for u in "${URLS[@]}"; do
            echo "   downloading $u"
            curl -fL --retry 3 -o "$archive" "$u" && break || rm -f "$archive"
        done
        [[ -s $archive ]] || { echo "ERROR: no eggnog.taxa archive in $SHARED_DB and download failed"; exit 1; }
    fi
    tmp=$(mktemp -d "$NEW_DB/taxa_tmp.XXXX")
    tar -xf "$archive" -C "$tmp"              # tar detects gzip by itself
    for f in eggnog.taxa.db eggnog.taxa.db.traverse.pkl; do
        hit=$(find "$tmp" -type f -name "$f" | head -1)
        [[ -n $hit ]] || { echo "ERROR: $f not inside $archive:"; tar -tf "$archive"; rm -rf "$tmp"; exit 1; }
        mv "$hit" "$NEW_DB/$f"
    done
    rm -rf "$tmp"
    [[ $archive == "$NEW_DB/"* ]] && rm -f "$archive"
fi

for f in eggnog.taxa.db eggnog.taxa.db.traverse.pkl; do
    [[ -s $NEW_DB/$f ]] || { echo "ERROR: $f missing after extraction:"; ls -l "$NEW_DB"; exit 1; }
done
ls -lh "$NEW_DB"

if [[ $(readlink -f "$EGGNOG_DB") != "$(readlink -f "$NEW_DB")" ]]; then
    cp func_config.sh func_config.sh.bak
    sed -i "s#^EGGNOG_DB=.*#EGGNOG_DB=$NEW_DB#" func_config.sh
fi
grep -q "^EGGNOG_DB=$NEW_DB\$" func_config.sh || { echo "ERROR: could not set EGGNOG_DB in func_config.sh"; exit 1; }
echo
echo "== Done. func_config.sh now has:  $(grep '^EGGNOG_DB=' func_config.sh)"
echo "   Next: bash 1_prepare.sh, then sbatch --array=1 2_annotate.slurm"
