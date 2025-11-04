#!/bin/bash
#SBATCH --job-name=alphafold3_rtt109_naf2			#Name your job something original
#SBATCH --partition=gpu_p			#Use the GPU partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32			#If you use the default options, AlphaFold3 will run four simutaneous Jackhmmer processes with 8 CPUs each
#SBATCH --gres=gpu:1				#If you don’t care whether your job uses an A100 node or an H100 node (and there isn’t much difference in run time)…
#SBATCH --constraint=Milan|SapphireRapids	#…this is the easiest way to specify either one without accidentally using a P100 or L4, which lack sufficient device memory
#SBATCH --mem=60gb
#SBATCH --time=10:00:00
#SBATCH --output=../alphafold3_rtt109_naf2.%j.out
#SBATCH --error=../alphafold3_rtt109_naf2.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu


cd $SLURM_SUBMIT_DIR
# Directories
INPUT_DIR="/scratch/ry00555/AlphaFold/input"
OUTPUT_DIR="/scratch/ry00555/AlphaFold/output"
MODEL_DIR="/home/ry00555/Research/FungiCutAndRun/AlphaFold3"        # where your AlphaFold3 parameters live
PUBLIC_DB="/db/AlphaFold3/20241114"                # Sapelo2 path for public database

singularity exec \
     --nv \
     --bind ${INPUT_DIR}:/root/af_input \
     --bind ${OUTPUT_DIR}:/root/af_output \
     --bind ${MODEL_DIR}:/root/models \
     --bind ${PUBLIC_DB}:/root/public_databases \
     /apps/singularity-images/alphafold-3.0.0.sif \
     python /app/alphafold/run_alphafold.py \
     --json_path=/root/af_input/RTT109_D145A_DD304_5AA_H3H4heterodimer.json \
     --model_dir=/root/models \
     --db_dir=/root/public_databases \
     --output_dir=/root/af_output

     singularity exec \
          --nv \
          --bind ${INPUT_DIR}:/root/af_input \
          --bind ${OUTPUT_DIR}:/root/af_output \
          --bind ${MODEL_DIR}:/root/models \
          --bind ${PUBLIC_DB}:/root/public_databases \
          /apps/singularity-images/alphafold-3.0.0.sif \
          python /app/alphafold/run_alphafold.py \
          --json_path=/root/af_input/RTT109D145A_H3H4heterodimer.json \
          --model_dir=/root/models \
          --db_dir=/root/public_databases \
          --output_dir=/root/af_output

     singularity exec \
          --nv \
          --bind ${INPUT_DIR}:/root/af_input \
          --bind ${OUTPUT_DIR}:/root/af_output \
          --bind ${MODEL_DIR}:/root/models \
          --bind ${PUBLIC_DB}:/root/public_databases \
          /apps/singularity-images/alphafold-3.0.0.sif \
          python /app/alphafold/run_alphafold.py \
          --json_path=/root/af_input/AF_3_RTT109_H3H3HeteroDimer.json \
          --model_dir=/root/models \
          --db_dir=/root/public_databases \
          --output_dir=/root/af_output

          singularity exec \
               --nv \
               --bind ${INPUT_DIR}:/root/af_input \
               --bind ${OUTPUT_DIR}:/root/af_output \
               --bind ${MODEL_DIR}:/root/models \
               --bind ${PUBLIC_DB}:/root/public_databases \
               /apps/singularity-images/alphafold-3.0.0.sif \
               python /app/alphafold/run_alphafold.py \
               --json_path=/root/af_input/AF3_RTT109_NAF2_H3H4tetramer.json \
               --model_dir=/root/models \
               --db_dir=/root/public_databases \
               --output_dir=/root/af_output
