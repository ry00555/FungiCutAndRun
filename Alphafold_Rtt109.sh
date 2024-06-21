#!/bin/bash
#SBATCH --job-name=j_GATK
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --time=08:00:00
#SBATCH --output=AlphaFold.%j.out
#SBATCH --error=AlphaFold.%j.err



#!/bin/bash

# Define paths
INPUT_LIST="/scratch/ry00555/AlphaFold/Rtt109_FilteredHits_Accessions.txt"
OUTPUT_DIR="/scratch/ry00555/AlphaFold/Rtt109_Accessions_Fastas"

# Check if the output directory exists, if not, create it
[ ! -d "$OUTPUT_DIR" ] && mkdir -p "$OUTPUT_DIR"

# Loop through each accession in the list
while IFS= read -r ACCESSION
do
    # Construct the output file path based on the accession
    OUTPUT_FILE="$OUTPUT_DIR/$ACCESSION.fa"

    # Fetch sequence from NCBI and write to the output file
    efetch -db protein -id "$ACCESSION" -format fasta > "$OUTPUT_FILE"

    # Check if sequence retrieval was successful
    if [ -s "$OUTPUT_FILE" ]; then
        echo "Sequence fetched for accession $ACCESSION and saved to $OUTPUT_FILE"
    else
        echo "WARNING: Sequence not found for accession $ACCESSION"
    fi

    # Break the loop after fetching the first sequence (for rtt109)
    if [ "$ACCESSION" = "NCU09825-t26_1-p1" ]; then
        break
    fi
done < "$INPUT_LIST"



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
