#!/bin/bash
#SBATCH --job-name=extract_ipTM
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --mem=8gb
#SBATCH --time=1:00:00
#SBATCH --output=%AlphaFold2_Step3_PullipTMs.%j.out
#SBATCH --error=%AlphaFold2_Step3_PullipTMs.%j.err



WORKDIR="/lustre2/scratch/ry00555/AlphaFold"
AF2_DIR="$WORKDIR/output/AF2"
FASTA_DIR="$WORKDIR/FastaforAlphaFold"
OUTPUT_CSV="$AF2_DIR/RTT109_MSHits_ipTM_scores.csv"

echo "protein,interactor1,interactor2,ipTM" > $OUTPUT_CSV

for dir in $AF2_DIR/*/; do
    protein=$(basename $dir)
    json="$dir/$protein/summary_confidences_0.json"
    fasta="$FASTA_DIR/${protein}.fa"

    # Extract the two sequence names from the fasta file
    interactors=($(grep "^>" $fasta | sed 's/>//'))
    interactor1=${interactors[0]}
    interactor2=${interactors[1]}

    if [ -f "$json" ]; then
        python3 -c "
import json
with open('$json') as f:
    data = json.load(f)
print('$protein,$interactor1,$interactor2,' + str(data.get('iptm', 'NA')))
" >> $OUTPUT_CSV
    else
        echo "$protein,$interactor1,$interactor2,MISSING" >> $OUTPUT_CSV
    fi
done

echo "Done! Scores written to $OUTPUT_CSV"