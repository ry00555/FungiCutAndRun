#!/bin/bash
#SBATCH --job-name=alphaMSAs_step1
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --time=8:00:00 #this is the longest step according to EMBL EBI; adjust accordingly Small proteins (<300 aa): 1-3 hours Medium proteins (300-800 aa): 3-8 hours Large proteins (>800 aa): 8-12 hours
#SBATCH --output=../AF2_MSABuilder.%j.out
#SBATCH --error=../AF2_MSABuilder.%j.err
#SBATCH --array=1-2 #change to the actual number of interactions otherwise x-200 will fail, 113 predictions took 4 hours
#Substep1 before running check and stop elements script
cd $SLURM_SUBMIT_DIR

ml purge
ml AlphaFold/2.3.2-foss-2023a-CUDA-12.1.1
export ALPHAFOLD_DATA_DIR=/db/AlphaFold/2.3.2
WORKDIR="/scratch/ry00555/AlphaFold"

file=$(awk "NR==${SLURM_ARRAY_TASK_ID}" $WORKDIR/FastaforAlphaFold/input.lst)
#test
alphafold \
--models_to_relax=none \
--data_dir=$ALPHAFOLD_DATA_DIR \
--model_preset=multimer \
--num_multimer_predictions_per_model=1 \
--max_template_date=2023-10-01 \
--db_preset=full_dbs \
--output_dir=$WORKDIR/output/AF2/$(basename $file .fa) \
--fasta_paths=$WORKDIR/FastaforAlphaFold/$file

#in the command line you can run ls /scratch/ry00555/AlphaFold/output/AF2/ | wc -l to see how many predictions you ran. the number doesn't mean how many were successfully made, run step 2 to check '
