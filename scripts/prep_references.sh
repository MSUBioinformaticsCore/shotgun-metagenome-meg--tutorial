#!/bin/bash --login
#SBATCH --job-name=mag_host_ref
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%j.log

# =============================================================================
#  Build ONE host reference for nf-core/mag (--host_fasta takes a single file):
#     dog   ROS_Cfam_1.0    (GCF_014441545.1) + dog mitogenome   NC_002008.4
#     human T2T-CHM13v2.0   (GCF_009914755.1) + human mitogenome NC_012920.1
#  Sequence names get a dog_ / human_ prefix, so each removed read can be
#  traced to its host.
#
#  Run once:  sbatch prep_references.sh   (or run directly on a dev node)
# =============================================================================

set -euo pipefail

REF_DIR=/mnt/research/YOUR_LAB/refs/host
OUT=$REF_DIR/dog_human_host.fna

NCBI=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF
EFETCH="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text"

mkdir -p "$REF_DIR" logs
cd "$REF_DIR"

# ---- helper: download a genome, add mitogenome if absent, prefix names -------
add_host () {
    local prefix=$1 url=$2 mito=$3
    local gz; gz=$(basename "$url")

    [[ -s $gz ]] || wget -c "$url"
    zcat "$gz" > "$prefix.fna"

    if grep -q "^>${mito%.*}" "$prefix.fna"; then
        echo "$prefix: mitogenome $mito already present"
    else
        echo "$prefix: adding mitogenome $mito"
        wget -qO "$mito.fna" "$EFETCH&id=$mito"
        grep -q '^>' "$mito.fna" || { echo "Download of $mito failed"; exit 1; }
        cat "$mito.fna" >> "$prefix.fna"
    fi

    sed "s/^>/>${prefix}_/" "$prefix.fna" >> "$OUT"
    rm "$prefix.fna"
}

# ---- build -------------------------------------------------------------------------
rm -f "$OUT"
add_host dog   "$NCBI/014/441/545/GCF_014441545.1_ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna.gz"   NC_002008.4
add_host human "$NCBI/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz" NC_012920.1

# ---- summary -----------------------------------------------------------------------
echo "dog sequences:   $(grep -c '^>dog_'   "$OUT")"
echo "human sequences: $(grep -c '^>human_' "$OUT")"
echo "mitogenomes:     $(grep -E '^>(dog|human)_(NC_002008|NC_012920)' "$OUT" | cut -c1-40)"
echo "Host reference:  $OUT"
