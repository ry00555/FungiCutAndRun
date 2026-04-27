#!/bin/bash
#SBATCH --job-name=alphafold_step2
#SBATCH --partition=gpu_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:A100:1 #You can specify other GPU devices, like 'H100', here
#SBATCH --mem=32gb                    
#SBATCH --time=12:00:00           
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --output=AlphaFold2_Step2_StructurePred.%j.out
#SBATCH --error=AlphaFold2_Step2_StructurePred.%j.err
#SBATCH --array=8,16,17,24,25,26,30,31,37,39,48,49,50,54,55,59,67,70,75,82,94,96,112

cd $SLURM_SUBMIT_DIR

ml purge
ml AlphaFold/2.3.2-foss-2023a-CUDA-12.1.1
export ALPHAFOLD_DATA_DIR=/db/AlphaFold/2.3.2
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_MEM_FRACTION=8

WORKDIR="/scratch/ry00555/AlphaFold"
file=$(awk "NR==${SLURM_ARRAY_TASK_ID}" $WORKDIR/FastaforAlphaFold/input.lst)

alphafold \
--models_to_relax=none \
--data_dir=$ALPHAFOLD_DATA_DIR \
--model_preset=multimer \
--use_precomputed_msas=true \
--num_multimer_predictions_per_model=1 \
--max_template_date=2023-10-01 \
--db_preset=full_dbs \
--output_dir=$WORKDIR/output/AF2/$(basename $file .fa) \
--fasta_paths=$WORKDIR/FastaforAlphaFold/$file