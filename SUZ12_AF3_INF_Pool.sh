#!/bin/bash
#SBATCH --job-name=SUZ12_INF
#SBATCH --partition=gpu_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:A100:1
#SBATCH --constraint=Milan|SapphireRapids
#SBATCH --mem=60gb
#SBATCH --time=4:00:00
#SBATCH --output=/scratch/ry00555/RNASeqPaper2026/Proteome/SingleBaitPools/logs/SUZ12_INF.%A_%a.out
#SBATCH --error=/scratch/ry00555/RNASeqPaper2026/Proteome/SingleBaitPools/logs/SUZ12_INF.%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ry00555@uga.edu

# ══════════════════════════════════════════════════════════════════════
# SUZ12_AF3_INF_Pool.sh — inference for one SUZ12 pool per array task.
# Submitted by SUZ12_AF3_Pools.sh (inf_scan/recheck modes) with --array=...
# NOTE: this pipeline's GPU jobs share the account-wide 8-running /
# 20-submitted gpu_qos cap with the other 3 bait pipelines -- see the
# INF_THROTTLE / INF_BATCH comments in SUZ12_AF3_Pools.sh before changing.
# ══════════════════════════════════════════════════════════════════════

PROTEIN="SUZ12"
BASE_DIR="/scratch/ry00555/RNASeqPaper2026/Proteome/SingleBaitPools"
INPUT_DIR="${BASE_DIR}/${PROTEIN}/${PROTEIN}_PooledPPI_InputJSONs"
OUTPUT_DIR="${BASE_DIR}/${PROTEIN}/${PROTEIN}_PooledPPI_output"
MODEL_DIR="/home/ry00555/Research/AlphaFold3ModelParameters"
PUBLIC_DB="/db/AlphaFold3/20241114"

file=$(ls $INPUT_DIR/*.json | awk "NR==${SLURM_ARRAY_TASK_ID}")
[ -z "$file" ] && echo "ERROR: no file at array index ${SLURM_ARRAY_TASK_ID}" && exit 1

job_name=$(basename $file .json | tr '[:upper:]' '[:lower:]')
MSA_JSON="${OUTPUT_DIR}/${job_name}/${job_name}_data.json"

echo "Running INF for: $job_name"
echo "Using MSA JSON: $MSA_JSON"

if [ ! -f "$MSA_JSON" ]; then
    echo "ERROR: Missing MSA JSON for $job_name at $MSA_JSON"
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

# ── Move timestamped output into the clean MSA folder ─────────────────────────
timestamped=$(ls -dt ${OUTPUT_DIR}/${job_name}_*/ 2>/dev/null | head -1)
if [ -n "$timestamped" ]; then
    echo "Moving output from $timestamped into ${OUTPUT_DIR}/${job_name}/"
    mv ${timestamped}/* ${OUTPUT_DIR}/${job_name}/
    rmdir ${timestamped}
    echo "Done — output consolidated into ${OUTPUT_DIR}/${job_name}/"
else
    echo "No timestamped folder found for $job_name"
fi
