#!/bin/bash
#SBATCH --job-name=EAF3_MSA
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=150gb
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/ry00555/RNASeqPaper2026/Proteome/SingleBaitPools/logs/EAF3_MSA.%A_%a.out
#SBATCH --error=/scratch/ry00555/RNASeqPaper2026/Proteome/SingleBaitPools/logs/EAF3_MSA.%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ry00555@uga.edu

# ══════════════════════════════════════════════════════════════════════
# EAF3_AF3_MSA_Pool.sh — MSA (data pipeline only) for one EAF3 pool
# per array task. Submitted by EAF3_AF3_Pools.sh with --array=N-M%THROTTLE
# ══════════════════════════════════════════════════════════════════════

PROTEIN="EAF3"
BASE_DIR="/scratch/ry00555/RNASeqPaper2026/Proteome/SingleBaitPools"
INPUT_DIR="${BASE_DIR}/${PROTEIN}/${PROTEIN}_PooledPPI_InputJSONs"
OUTPUT_DIR="${BASE_DIR}/${PROTEIN}/${PROTEIN}_PooledPPI_output"
MODEL_DIR="/home/ry00555/Research/AlphaFold3ModelParameters"
PUBLIC_DB="/db/AlphaFold3/20241114"

mkdir -p "$OUTPUT_DIR"

file=$(ls $INPUT_DIR/*.json | awk "NR==${SLURM_ARRAY_TASK_ID}")
[ -z "$file" ] && echo "ERROR: no file at array index ${SLURM_ARRAY_TASK_ID}" && exit 1

echo "Running MSA for: $(basename $file)"

singularity exec \
     --bind ${INPUT_DIR}:/root/af_input \
     --bind ${OUTPUT_DIR}:/root/af_output \
     --bind ${MODEL_DIR}:/root/models \
     --bind ${PUBLIC_DB}:/root/public_databases \
     /apps/singularity-images/alphafold-3.0.1.sif \
     python /app/alphafold/run_alphafold.py \
     --json_path=/root/af_input/$(basename $file) \
     --model_dir=/root/models \
     --db_dir=/root/public_databases \
     --output_dir=/root/af_output \
     --run_data_pipeline=true \
     --run_inference=false
