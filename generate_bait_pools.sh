#!/bin/bash
#SBATCH --job-name=generate_bait_pools
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/ry00555/RNASeqPaper2026/Proteome/SingleBaitPools/logs/generate_pools_%j.out
#SBATCH --error=/scratch/ry00555/RNASeqPaper2026/Proteome/SingleBaitPools/logs/generate_pools_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ry00555@uga.edu

# ══════════════════════════════════════════════════════════════════════
# generate_bait_pools.sh
#
# Generates AF3 pool JSONs for 4 individual bait screens:
#   SET7 / SUZ12 / EED / EAF3 — each as chain A
#   Candidates = full N. crassa proteome on chains B, C, D...
#
# Output structure:
#   /scratch/ry00555/RNASeqPaper2026/Proteome/SingleBaitPools/
#     SET7/
#       SET7_PooledPPI_InputJSONs/
#         set7_pool0001.json ...
#         SET7_pool_key.csv
#         SET7_pool_key.xlsx
#     SUZ12/
#       SUZ12_PooledPPI_InputJSONs/
#     EED/
#       EED_PooledPPI_InputJSONs/
#     EAF3/
#       EAF3_PooledPPI_InputJSONs/
#
# Upload to cluster before running:
#   scp generate_individual_bait_pools.py ry00555@xfer.gacrc.uga.edu:/home/ry00555/Research/FungiCutAndRun/
#   scp baits_SET7_SUZ12_EED_EAF3.fasta   ry00555@xfer.gacrc.uga.edu:/scratch/ry00555/RNASeqPaper2026/Proteome/SingleBaitPools/
#   scp StringDB_Ncrassa_Sequences.txt     ry00555@xfer.gacrc.uga.edu:/scratch/ry00555/RNASeqPaper2026/Proteome/SingleBaitPools/
# ══════════════════════════════════════════════════════════════════════

set -euo pipefail

SCRATCH="/scratch/ry00555/RNASeqPaper2026/Proteome/SingleBaitPools"
SCRIPTS="/home/ry00555/Research/FungiCutAndRun"

PROTEOME="${SCRATCH}/StringDB_Ncrassa_Sequences.txt"
BAIT_FASTA="${SCRATCH}/baits_SET7_SUZ12_EED_EAF3.fasta"

mkdir -p "${SCRATCH}/logs"

# Create output subdirectories
for BAIT in SET7 SUZ12 EED EAF3; do
    mkdir -p "${SCRATCH}/${BAIT}/${BAIT}_PooledPPI_InputJSONs"
done

echo "════════════════════════════════════════════════════════════"
echo " Individual Bait Pool Generator"
echo " Job ID  : ${SLURM_JOB_ID}"
echo " Started : $(date)"
echo " Proteome: ${PROTEOME}"
echo " Baits   : SET7 SUZ12 EED EAF3"
echo "════════════════════════════════════════════════════════════"

# Check input files exist
[[ ! -f "${PROTEOME}" ]]   && echo "ERROR: proteome FASTA not found: ${PROTEOME}"   && exit 1
[[ ! -f "${BAIT_FASTA}" ]] && echo "ERROR: bait FASTA not found: ${BAIT_FASTA}"     && exit 1

echo ""
echo "Input files confirmed."
echo "Proteome: $(grep -c '^>' ${PROTEOME}) sequences"
echo "Baits:    $(grep -c '^>' ${BAIT_FASTA}) sequences"
echo ""

# Run pool generator
python3 "${SCRIPTS}/generate_individual_bait_pools.py" \
    --proteome_fasta "${PROTEOME}" \
    --bait_fasta     "${BAIT_FASTA}" \
    --baits          SET7 SUZ12 EED EAF3 \
    --outdir         "${SCRATCH}" \
    --max_aa         5000 \
    --max_candidates 6 \
    --seeds          42 137

echo ""
echo "Pool generation complete: $(date)"
echo ""

# Move output into correct subdirectory structure
# The script writes to ${SCRATCH}/SET7_pools/ etc
# We need to move to ${SCRATCH}/SET7/SET7_PooledPPI_InputJSONs/
echo "Reorganising output directories..."

for BAIT in SET7 SUZ12 EED EAF3; do
    SRC="${SCRATCH}/${BAIT}_pools"
    DST="${SCRATCH}/${BAIT}/${BAIT}_PooledPPI_InputJSONs"

    if [[ -d "${SRC}" ]]; then
        mv "${SRC}"/* "${DST}/"
        rmdir "${SRC}"
        N=$(ls "${DST}"/*.json 2>/dev/null | wc -l)
        echo "  ${BAIT}: ${N} JSONs -> ${DST}"
    else
        echo "  WARNING: ${SRC} not found"
    fi
done

echo ""
echo "════════════════════════════════════════════════════════════"
echo " Done: $(date)"
echo ""
echo " Output directories:"
for BAIT in SET7 SUZ12 EED EAF3; do
    N=$(ls "${SCRATCH}/${BAIT}/${BAIT}_PooledPPI_InputJSONs"/*.json 2>/dev/null | wc -l)
    echo "   ${BAIT}: ${N} pool JSONs"
    echo "         ${SCRATCH}/${BAIT}/${BAIT}_PooledPPI_InputJSONs/"
done
echo ""
echo " Next steps:"
echo "   1. Submit MSA + INF jobs for each bait using PRC2_AF3_Pools.sh pattern"
echo "   2. Run AF3pool_scraper_v1.py on each bait output directory"
echo "════════════════════════════════════════════════════════════"
