#!/bin/bash
#SBATCH --job-name=EAF3
#SBATCH --partition=gpu_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gres=gpu:A100:1
#SBATCH --mem=60gb
#SBATCH --time=10:00:00
#SBATCH --output=AF3Pool.%A_%a.out
#SBATCH --error=AF3Pool.%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --array=1-20%1
 
#change this to the number of input jsons the %10 means it will run only 10 at a time but it will go through all 113, otherwise you can do 1-20, 21-40 etc



cd $SLURM_SUBMIT_DIR

# Directories
INPUT_DIR="/scratch/ry00555/EpigeneticMemoryPaper2026/AlphaFold3/Input/EAF3_MS_JsonInputs"
OUTPUT_DIR="/scratch/ry00555/EpigeneticMemoryPaper2026/AlphaFold3/Output/EAF3"
MODEL_DIR="/home/ry00555/Research/FungiCutAndRun/AlphaFold3"
PUBLIC_DB="/db/AlphaFold3/20241114"

# Pick the JSON file for this array task
file=$(ls $INPUT_DIR/*.json | awk "NR==${SLURM_ARRAY_TASK_ID}")

singularity exec \
     --nv \
     --bind ${INPUT_DIR}:/root/af_input \
     --bind ${OUTPUT_DIR}:/root/af_output \
     --bind ${MODEL_DIR}:/root/models \
     --bind ${PUBLIC_DB}:/root/public_databases \
     /apps/singularity-images/alphafold-3.0.1.sif \
     python /app/alphafold/run_alphafold.py \
     --json_path=/root/af_input/$(basename $file) \
     --model_dir=/root/models \
     --db_dir=/root/public_databases \
     --output_dir=/root/af_output