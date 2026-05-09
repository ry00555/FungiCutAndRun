#!/bin/bash
#SBATCH --job-name=EAF3
#SBATCH --partition=gpu_p
#SBATCH --gres=gpu:H100:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --output=/scratch/ry00555/Boltz2/logs/boltz2_%A_%a.out
#SBATCH --error=/scratch/ry00555/Boltz2/logs/boltz2_%A_%a.err

module load Boltz/2.2.0

PROJECT="EAF3" #change this and job name to dock protein

BOLTZDIR="/scratch/ry00555/Boltz2/$PROJECT/boltz2_inputs"
OUTDIR="/scratch/ry00555/Boltz2/$PROJECT/boltz2_outputs"

mkdir -p $OUTDIR

YAML_FILES=("$BOLTZDIR"/*.yaml)

    INPUT="${YAML_FILES[$SLURM_ARRAY_TASK_ID]}"
   BASENAME=$(basename "$INPUT" .yaml | sed 's/binder_\([0-9]*\)_\([^-]*\).*/\1_\2/')

    echo "[$(date)] Running Boltz-2 on: $INPUT"

    boltz predict "$INPUT" \
        --out_dir  $OUTDIR \
--cache "/scratch/ry00555/Boltz2/" \
        --use_msa_server  \
        --diffusion_samples 5 \
        --diffusion_samples_affinity 5 \
        --accelerator  gpu \
        --devices      1

    echo "[$(date)] Done: $BASENAME"
