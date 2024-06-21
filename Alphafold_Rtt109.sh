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

# Load the BLAST+ module
module load BLAST+/2.7.1-foss-2016b-Python-2.7.14

# Check if the accession file exists
if [ ! -f "$accession_file" ]; then
  echo "Accession file not found: $accession_file"
  exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each accession and extract amino acid sequence to a separate FASTA file
while read -r accession; do
  output_file="$output_dir/$accession.fa"
  blastdbcmd -db nr -entry "$accession" -dbtype prot -out "$output_file"
done < "$accession_file"

echo "Amino acid sequences extracted and saved to $output_dir"
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
