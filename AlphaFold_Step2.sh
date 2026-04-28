#!/bin/bash
#SBATCH --job-name=alphafold_step2
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --gres=gpu:V100:1 
#SBATCH --time=4:00:00
#SBATCH --output=StructurePred.%j.out
#SBATCH --error=StructurePred.%j.err
#SBATCH --array=1-113

cd $SLURM_SUBMIT_DIR

#check if the msas were properly made from step1
for dir in /scratch/ry00555/AlphaFold/output/AF2/*/; do
    protein=$(basename $dir)
    if [ ! -d "$dir/$protein/msas/A" ] || [ ! -d "$dir/$protein/msas/B" ]; then
        echo "INCOMPLETE: $protein"
    fi
done
    
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
