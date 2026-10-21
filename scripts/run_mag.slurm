#!/bin/bash --login
#SBATCH --job-name=nfcore_mag
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.log
#SBATCH --mail-type=END,FAIL

# =============================================================================
#  nf-core/mag 5.5.0 on the MSU HPCC
#
#  This job runs only the Nextflow head process (2 CPUs, 16 GB). Nextflow submits
#  every pipeline step as its own SLURM job.
#  Submit from the run folder:   sbatch run_mag.slurm
#  To resume after a failure or the 7-day limit, submit again (-resume reuses
#  every finished step, as long as WORKDIR is the same).
# =============================================================================

set -euo pipefail

# ---- 1. Settings (edit) --------------------------------------------------------
RUN_DIR=/mnt/scratch/$USER/mag                       # run folder (logs, .nextflow/)
SAMPLESHEET=$RUN_DIR/samplesheet.csv                  # from make_samplesheet.py
OUTDIR=$RUN_DIR/results
WORKDIR=$RUN_DIR/work                                 # keep on scratch; needed for -resume
HOST_FASTA=/mnt/research/YOUR_LAB/refs/host/dog_human_host.fna   # from prep_references.sh
PIPELINE_DIR=/mnt/scratch/$USER/pipelines/nf-core-mag-5.5.0      # pipeline code
ACCOUNT=""                                            # SLURM account for the pipeline jobs,
                                                      # e.g. "data-machine"; "" = your default

# ---- 2. Checks ---------------------------------------------------------------------
mkdir -p "$RUN_DIR/logs"
cd "$RUN_DIR"
[[ -s $SAMPLESHEET ]] || { echo "ERROR: samplesheet not found: $SAMPLESHEET"; exit 1; }
[[ -s $HOST_FASTA ]]  || { echo "ERROR: host reference not found: $HOST_FASTA"; exit 1; }
export NXF_OPTS="-Xms1g -Xmx8g"

if [[ ! -d $PIPELINE_DIR/bin ]]; then
    git clone --branch 5.5.0 --depth 1 https://github.com/nf-core/mag.git "$PIPELINE_DIR"
fi
chmod +x "$PIPELINE_DIR"/bin/*

# Use ICER's system Apptainer/Singularity, not a conda copy
NEXTFLOW=$(command -v nextflow) || { echo "ERROR: nextflow not found"; exit 1; }
export JAVA_CMD=$(command -v java) || { echo "ERROR: java not found"; exit 1; }
export PATH=$(echo "$PATH" | tr ':' '\n' | grep -v -E 'miniforge|miniconda|anaconda' | paste -sd:)
case "$(command -v singularity)" in
    *conda*|*miniforge*|"") echo "ERROR: system singularity/apptainer not found"; exit 1 ;;
esac

# ---- 3. Nextflow settings for the HPCC ----------------------------------------------
cat > hpcc.config <<EOF
process {
    executor       = 'slurm'
    ${ACCOUNT:+clusterOptions = '--account=$ACCOUNT'}
    resourceLimits = [ cpus: 128, memory: 1900.GB, time: 168.h ]
    errorStrategy  = { task.exitStatus in [137,140,143] ? 'retry' : 'finish' }
    maxRetries     = 2
    withName: '.*:BOWTIE2_HOST_REMOVAL_BUILD' { cpus = 16; memory = 64.GB; time = 24.h }
    withName: '.*:MEGAHIT' { cpus = 32; memory = { 120.GB * task.attempt }; time = { 48.h * task.attempt } }
}
executor {
    queueSize       = 100
    submitRateLimit = '10/1min'
}
singularity {
    enabled    = true
    autoMounts = true
}
EOF

# ---- 4. Run nf-core/mag ---------------------------------------------------------------
"$NEXTFLOW" run "$PIPELINE_DIR" \
    -profile singularity \
    -c hpcc.config \
    -work-dir "$WORKDIR" \
    -resume \
    --input "$SAMPLESHEET" \
    --outdir "$OUTDIR" \
    --host_fasta "$HOST_FASTA" \
    --reads_minlength 50 \
    --skip_spades \
    --binning_map_mode group \
    --skip_concoct --skip_comebin --skip_semibin \
    --bin_domain_classification \
    --exclude_unbins_from_postbinning \
    --refine_bins_dastool \
    --postbinning_input refined_bins_only \
    --run_busco false \
    --run_checkm2

# Job information and resource usage
scontrol show job "$SLURM_JOB_ID"
module load powertools 2>/dev/null && js -j "$SLURM_JOB_ID" || true
