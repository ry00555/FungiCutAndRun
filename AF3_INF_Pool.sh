#!/bin/bash
#SBATCH --job-name=EAF3_INF
#SBATCH --partition=gpu_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:A100:1
#SBATCH --constraint=Milan|SapphireRapids
#SBATCH --mem=60gb
#SBATCH --time=4:00:00
#SBATCH --output=INF.%j.out
#SBATCH --error=INF.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu

# ── Only change this line when swapping proteins ──────────────────────────────
PROTEIN="EAF3"
# ─────────────────────────────────────────────────────────────────────────────

cd $SLURM_SUBMIT_DIR

BASE_DIR="/scratch/ry00555/EpigeneticMemoryPaper2026/AlphaFold3/PooledPPI"
INPUT_DIR="${BASE_DIR}/${PROTEIN}_AF3_PooledJSONs"
OUTPUT_DIR="${BASE_DIR}/${PROTEIN}_AF3_PooledJSON_output"
MODEL_DIR="/home/ry00555/Research/AlphaFold3ModelParameters"
PUBLIC_DB="/db/AlphaFold3/20241114"

# Get the json file for this array task
file=$(ls $INPUT_DIR/*.json | awk "NR==${SLURM_ARRAY_TASK_ID}")

# Derive the output folder name — AF3 lowercases the job name
job_name=$(basename $file .json | tr '[:upper:]' '[:lower:]')
MSA_JSON="${OUTPUT_DIR}/${job_name}/${job_name}_data.json"

echo "Running INF for: $job_name"
echo "Using MSA JSON: $MSA_JSON"

# Check MSA exists before running
if [ ! -f "$MSA_JSON" ]; then
    echo "ERROR: Missing MSA JSON for $job_name at $MSA_JSON"
    echo "Available folders:"
    ls $OUTPUT_DIR/
    exit 1
fi

export XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"

singularity exec \
     --nv \
     --bind ${INPUT_DIR}:/root/af_input \
     --bind ${OUTPUT_DIR}:/root/af_output \
     --bind ${MODEL_DIR}:/root/models \
     --bind ${PUBLIC_DB}:/root/public_databases \
     /apps/singularity-images/alphafold-3.0.1.sif \
     python /app/alphafold/run_alphafold.py \
     --json_path=/root/af_output/${job_name}/${job_name}_data.json \
     --model_dir=/root/models \
     --db_dir=/root/public_databases \
     --output_dir=/root/af_output \
     --run_data_pipeline=false \
     --run_inference=true
