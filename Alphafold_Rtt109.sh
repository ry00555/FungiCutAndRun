#!/bin/bash
#SBATCH --job-name=AlphaFold
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --time=08:00:00
#SBATCH --output=AlphaFold.%j.out
#SBATCH --error=AlphaFold.%j.err



#!/bin/bash

accession_file="/scratch/ry00555/AlphaFold/Rtt109_FilteredHits_Accessions.txt"
output_dir="/scratch/ry00555/AlphaFold/Rtt109_Accessions_Fastas"

for i in $(cat $accession_file)
do
  curl -s https://www.uniprot.org/uniprotkb?query=$i.fasta.gz | gunzip -c > $output_dir/$i.fa
#grep -v study $i.txt - this would remove any line that has the word study - not robust
#separate tabs by comlumn NR - take only the second row print NF only the last column
done


# cd $SLURM_SUBMIT_DIR
# ml purge
# ml AlphaFold/2.3.1-foss-2022a-CUDA-11.7.0
# export ALPHAFOLD_DATA_DIR=/apps/db/AlphaFold/2.3.1
#
# alphafold --data_dir=$ALPHAFOLD_DATA_DIR  \
# --model_preset=multimer \
# --max_template_date=2023-1-1 \
# --db_preset=full_dbs \
# --output_dir=./ \
# --fasta_paths=/scratch/zlewis/alphafold/MULTI/mus30_parp.fa
