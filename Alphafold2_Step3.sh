#!/bin/bash
#SBATCH --job-name=extract_iptm
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --time=2:00:00
#SBATCH --output=AF2_Step3.%j.out
#SBATCH --error=AF2_Step3.%j.err


WORKDIR="/scratch/ry00555/AlphaFold"
AF2_DIR="$WORKDIR/output/AF2"
FASTA_DIR="$WORKDIR/FastaforAlphaFold"
OUTPUT_CSV="$WORKDIR/output/AF2/RTT109_MSHits_ipTM_scores.csv"

ml purge
ml AlphaFold/2.3.2-foss-2023a-CUDA-12.1.1

echo "interactor2,best_model,iptm_best,ptm_best,ranking_confidence_best,mean_iptm,iptm_sd,mean_contacts_5A,mean_contacts_10A" > $OUTPUT_CSV

for dir in $AF2_DIR/*/; do
    protein=$(basename $dir)
    inner_dir="$dir/$protein"
    fasta="$FASTA_DIR/${protein}.fa"

    # Skip if fasta doesn't exist
    if [ ! -f "$fasta" ]; then
        echo "SKIPPING $protein - no fasta found"
        continue
    fi

    # Skip if ranking_debug.json doesn't exist
    if [ ! -f "$inner_dir/ranking_debug.json" ]; then
        echo "SKIPPING $protein - no ranking_debug.json"
        continue
    fi

    interactor2=$(grep "^>" $fasta | sed 's/>//' | sed -n '2p')

    python3 << PYEOF
import json, os, pickle
import numpy as np

inner_dir = '$inner_dir'
interactor2 = '$interactor2'
fasta_path = '$fasta'

# Chain lengths from fasta
chain_lengths = []
current_seq = ''
with open(fasta_path) as f:
    for line in f:
        if line.startswith('>'):
            if current_seq:
                chain_lengths.append(len(current_seq))
            current_seq = ''
        else:
            current_seq += line.strip()
    if current_seq:
        chain_lengths.append(len(current_seq))

len_A = chain_lengths[0]
len_B = chain_lengths[1]

# Best model from ranking_debug.json
with open(os.path.join(inner_dir, 'ranking_debug.json')) as f:
    ranking = json.load(f)
best_model = ranking.get('order', [None])[0]
best_idx   = int(best_model.split('_')[1])

# Extract from all 5 pkl files
iptm_vals      = []
pae_5A_counts  = []
pae_10A_counts = []
iptm_best      = 'NA'
ptm_best       = 'NA'
ranking_best   = 'NA'

for i in range(1, 6):
    pkl_path = os.path.join(inner_dir, f'result_model_{i}_multimer_v3_pred_0.pkl')
    with open(pkl_path, 'rb') as f:
        result = pickle.load(f)

    iptm = float(result.get('iptm'))
    ptm  = float(result.get('ptm'))
    rc   = float(result.get('ranking_confidence'))
    iptm_vals.append(iptm)

    if i == best_idx:
        iptm_best    = round(iptm, 4)
        ptm_best     = round(ptm, 4)
        ranking_best = round(rc, 4)

    pae = np.array(result.get('predicted_aligned_error'))
    AB  = pae[0:len_A, len_A:len_A+len_B]
    BA  = pae[len_A:len_A+len_B, 0:len_A]
    pae_5A_counts.append(int(np.sum(AB < 5)  + np.sum(BA < 5)))
    pae_10A_counts.append(int(np.sum(AB < 10) + np.sum(BA < 10)))

mean_iptm = round(float(np.mean(iptm_vals)), 4)
iptm_sd   = round(float(np.std(iptm_vals)), 4)
mean_5A   = round(float(np.mean(pae_5A_counts)), 1)
mean_10A  = round(float(np.mean(pae_10A_counts)), 1)

print(f'{interactor2},{best_model},{iptm_best},{ptm_best},{ranking_best},{mean_iptm},{iptm_sd},{mean_5A},{mean_10A}')
PYEOF

done >> $OUTPUT_CSV

echo "Done! Scores written to $OUTPUT_CSV"
