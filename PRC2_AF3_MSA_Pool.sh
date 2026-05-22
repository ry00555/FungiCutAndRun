#!/bin/bash
#SBATCH --job-name=PRC2_MSA
#SBATCH --partition=batch #All MSAs (batch partition, no GPU cost) → All INFs (gpu_p, batched)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120gb
#SBATCH --time=16:00:00
#SBATCH --output=MSA.%j.out
#SBATCH --error=MSA.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu

# ── Only change this line when swapping proteins ──────────────────────────────
PROTEIN="PRC2"
# ─────────────────────────────────────────────────────────────────────────────

cd $SLURM_SUBMIT_DIR

BASE_DIR="/scratch/ry00555/RNASeqPaper2026/Proteome/PRC2_Proteome_pools/PooledPPI"
INPUT_DIR="${BASE_DIR}/${PROTEIN}_AF3_PooledJSONs"
OUTPUT_DIR="${BASE_DIR}/${PROTEIN}_AF3_PooledJSON_output"
MODEL_DIR="/home/ry00555/Research/AlphaFold3ModelParameters"
PUBLIC_DB="/db/AlphaFold3/20241114"

mkdir -p $OUTPUT_DIR

file=$(ls $INPUT_DIR/*.json | awk "NR==${SLURM_ARRAY_TASK_ID}")

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
     --run_inference=false  # this is essentially Ankur's run_msa.py script
