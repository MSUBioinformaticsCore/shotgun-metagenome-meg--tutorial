# =============================================================================
#  Settings for the functional annotation scripts. Edit ONLY this file.
#  Read by 1_prepare.sh, 2_annotate.slurm and 3_merge.slurm.
# =============================================================================

# nf-core/mag output folder and the samplesheet given to mag
MAG_OUTDIR=/mnt/scratch/$USER/mag/results
SAMPLESHEET=/mnt/scratch/$USER/mag/samplesheet.csv

# eggNOG-mapper database (must contain eggnog.db and eggnog_proteins.dmnd)
EGGNOG_DB=/mnt/research/bioinformaticsCore/shared/eggnog_db

# Results: per-sample annotations in $OUT/<sample>/, final tables in $OUT/tables/
OUT=/mnt/scratch/$USER/mag/function
TMP=/mnt/scratch/$USER/func_tmp

# Software: no environment is activated. The scripts use whatever emapper.py and
# python3 are on your PATH when you run 1_prepare.sh / submit with sbatch
# (sbatch passes your current environment to the job).
# -----------------------------------------------------------------------------
require_tools () {            # usage: require_tools python3 emapper.py
    for t in "$@"; do
        command -v "$t" >/dev/null || { echo "ERROR: $t not on PATH; submit from a shell where it works"; exit 1; }
    done
}

# Reference files for the bacteria-only filter (4_build_exclusions.sh, 5_filter.slurm)
REF_DIR=/mnt/scratch/$USER/mag/function/reference
