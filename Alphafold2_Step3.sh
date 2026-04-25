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

echo "interactor2,ipTM_model0,mean_ipTM,mean_PAE,ipTM_sd" > $OUTPUT_CSV

for dir in $AF2_DIR/*/; do
    protein=$(basename $dir)
    fasta="$FASTA_DIR/${protein}.fa"
    interactor2=$(grep "^>" $fasta | sed 's/>//' | sed -n '2p')

    python3 -c "
import json, statistics

workdir = '$dir/$protein'
interactor2 = '$interactor2'

# Get ipTM from model 0
try:
    with open(workdir + '/summary_confidences_0.json') as f:
        data0 = json.load(f)
    iptm_model0 = data0.get('iptm', 'NA')
except:
    iptm_model0 = 'NA'

# Read all 5 models for ipTM and SD
iptm_vals = []
iptm_sd_vals = []
pae_vals = []

for i in range(5):
    # ipTM and SD from summary_confidences
    try:
        with open(workdir + f'/summary_confidences_{i}.json') as f:
            data = json.load(f)
        iptm = data.get('iptm', None)
        if iptm is not None:
            iptm_vals.append(iptm)
        samples = data.get('diffusion_samples', {})
        sd_iptm = samples.get('iptm', [])
        if sd_iptm:
            iptm_sd_vals.extend(sd_iptm)
    except:
        pass

    # PAE matrix from full_data files
    try:
        with open(workdir + f'/_full_data_{i}.json') as f:
            full_data = json.load(f)
        pae_matrix = full_data.get('pae', None)
        if pae_matrix is not None:
            # Flatten and average the entire matrix
            flat = [val for row in pae_matrix for val in row]
            pae_vals.append(statistics.mean(flat))
    except:
        pass

mean_iptm = round(statistics.mean(iptm_vals), 4) if iptm_vals else 'NA'
mean_pae  = round(statistics.mean(pae_vals), 4)  if pae_vals  else 'NA'
iptm_sd   = round(statistics.stdev(iptm_sd_vals), 4) if len(iptm_sd_vals) > 1 else 'NA'

print(f'{interactor2},{iptm_model0},{mean_iptm},{mean_pae},{iptm_sd}')
" >> $OUTPUT_CSV

done

echo "Done! Scores written to $OUTPUT_CSV"