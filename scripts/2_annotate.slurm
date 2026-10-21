#!/bin/bash
#SBATCH --job-name=func_annot
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G                      # eggNOG-mapper --dbmem loads ~45 GB
#SBATCH --output=logs/%x_%A_%a.log

# =============================================================================
#  Step 2: eggNOG-mapper on mag's Prodigal proteins, one sample per task.
#  Submit from the script folder after 1_prepare.sh:
#     sbatch --array=1 2_annotate.slurm         (test one sample)
#     sbatch --array=1-60%20 2_annotate.slurm   (all; 20 at a time)
#  Task N = row N of func_manifest.tsv. A sample is skipped only if it was already
#  annotated with the same settings (emapper_settings.txt); otherwise it is re-run.
#
#  Bacteria-only annotation (eggNOG-mapper 2.1 options):
#    --tax_scope bacteria        a protein is annotated only if its best hit (seed ortholog)
#                                belongs to a bacterial orthologous group; proteins whose
#                                best hit is eukaryotic, archaeal or viral get no annotation
#    --target_taxa 2             functions are transferred only from bacterial orthologs
#    --excluded_taxa 1117,1090   ... and not from Cyanobacteria or Chlorobi orthologs
#  The search itself still runs against the whole eggNOG database, so a host or diet
#  protein finds its true (non-bacterial) best hit and is left unannotated rather than
#  being forced onto a bacterial protein. eggNOG_OGs still lists the full lineage of the
#  best hit, which 5_filter.slurm uses.
# =============================================================================

TAX_OPTS=(--tax_scope bacteria --target_taxa 2 --excluded_taxa 1117,1090)
# TAX_OPTS=()     # <- use this line instead for the eggNOG-mapper defaults (all taxa)

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ./func_config.sh
require_tools emapper.py
if (( ${#TAX_OPTS[@]} )) && [[ ! -s $EGGNOG_DB/eggnog.taxa.db ]]; then
    echo "ERROR: $EGGNOG_DB/eggnog.taxa.db (NCBI taxonomy) is needed for --target_taxa / --excluded_taxa"
    echo "       Run once on a dev node: bash setup_eggnog_taxa.sh  (it also updates EGGNOG_DB in func_config.sh)"
    exit 1
fi

IFS=$'\t' read -r SAMPLE FAA _GFF _DEPTH _FASTP \
    < <(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" func_manifest.tsv)
[[ -n ${SAMPLE:-} ]] || { echo "No manifest row for task $SLURM_ARRAY_TASK_ID"; exit 1; }

D=$OUT/$SAMPLE
W=$TMP/$SAMPLE
EMAPPER_VERSION=$(emapper.py --data_dir "$EGGNOG_DB" --version 2>&1 | sed 's/\x1b\[[0-9;]*m//g' | grep -m1 '^emapper' || true)
[[ -n $EMAPPER_VERSION ]] || { echo "ERROR: emapper.py --version gave no version line"; exit 1; }
SETTINGS="$EMAPPER_VERSION | ${TAX_OPTS[*]:-default taxa}"
echo "== $SAMPLE  (task $SLURM_ARRAY_TASK_ID, $(date))"
echo "   settings: $SETTINGS"

if [[ -s $D/$SAMPLE.emapper.annotations && -s $D/emapper_settings.txt ]] \
   && [[ $(cat "$D/emapper_settings.txt") == "$SETTINGS" ]]; then
    echo "already annotated with these settings; skipping"
    exit 0
fi
mkdir -p "$D" "$W"
rm -f "$D/emapper_settings.txt"

# proteins without Prodigal's trailing stop '*'
zcat -f "$FAA" | sed '/^>/!s/\*$//' > "$W/proteins.faa"
echo "   proteins: $(grep -c '^>' "$W/proteins.faa")"

emapper.py \
    -i "$W/proteins.faa" --itype proteins \
    -m diamond \
    --data_dir "$EGGNOG_DB" \
    --dbmem \
    "${TAX_OPTS[@]}" \
    --cpu "$SLURM_CPUS_PER_TASK" \
    --output "$SAMPLE" \
    --output_dir "$D" \
    --temp_dir "$W" \
    --override

echo "$SETTINGS" > "$D/emapper_settings.txt"
rm -rf "$W"
echo "   annotated genes: $(grep -vc '^#' "$D/$SAMPLE.emapper.annotations")"
echo "== $SAMPLE done ($(date))"
