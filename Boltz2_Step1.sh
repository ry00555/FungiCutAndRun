#!/bin/bash
#SBATCH --job-name=PRC2_boltz
#SBATCH --partition=gpu_p
#SBATCH --gres=gpu:H100:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=06:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --output=/scratch/ry00555/Boltz2/logs/boltz2_%A_%a.out
#SBATCH --error=/scratch/ry00555/Boltz2/logs/boltz2_%A_%a.err

module load Boltz/2.2.0

# ── Only change these lines when swapping proteins ────────────────────────────
PROJECT="PRC2"
BOLTZDIR="/scratch/ry00555/RNASeqPaper2026/Proteome/PRC2_Proteome_pools/Boltz2_Inputs"
OUTDIR="/scratch/ry00555/RNASeqPaper2026/Proteome/PRC2_Proteome_pools/Boltz2_Output"
# ─────────────────────────────────────────────────────────────────────────────

mkdir -p $OUTDIR /scratch/ry00555/Boltz2/logs

YAML_FILES=("$BOLTZDIR"/*.yaml)
INPUT="${YAML_FILES[$SLURM_ARRAY_TASK_ID]}"
BASENAME=$(basename "$INPUT" .yaml)

echo "[$(date)] Running Boltz-2 on: $INPUT"

boltz predict "$INPUT" \
    --out_dir  $OUTDIR \
    --cache    "/scratch/ry00555/Boltz2/cache" \
    --use_msa_server \
    --diffusion_samples 3 \
    --diffusion_samples_affinity 3 \
    --accelerator gpu \
    --devices 1

echo "[$(date)] Done: $BASENAME"
