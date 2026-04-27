#!/bin/bash
#SBATCH --job-name=extract_AF2_scores
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --time=2:00:00
#SBATCH --output=extract_AF2_scores.%j.out
#SBATCH --error=extract_AF2_scores.%j.err

WORKDIR="/scratch/ry00555/AlphaFold"
AF2_DIR="$WORKDIR/output/AF2"
FASTA_DIR="$WORKDIR/FastaforAlphaFold"
OUTPUT_CSV="$WORKDIR/ipTM_scores.csv"

ml purge
ml AlphaFold/2.3.2-foss-2023a-CUDA-12.1.1

echo "interactor2,best_model,iptm_best,ptm_best,ranking_confidence_best,mean_iptm,iptm_sd,mean_contacts_5A,mean_contacts_10A" > $OUTPUT_CSV

for dir in $AF2_DIR/*/; do
    protein=$(basename $dir)
    inner_dir="$dir/$protein"
    fasta="$FASTA_DIR/${protein}.fa"

    interactor2=$(grep "^>" $fasta | sed 's/>//' | sed -n '2p')

    python3 -c "
import json, os, pickle
import numpy as np

inner_dir = '$inner_dir'
interactor2 = '$interactor2'
fasta_path = '$fasta'

# ── Chain lengths from fasta ──────────────────────────────────────────────────
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

# ── Best model from ranking_debug.json ───────────────────────────────────────
try:
    with open(os.path.join(inner_dir, 'ranking_debug.json')) as f:
        ranking = json.load(f)
    best_model = ranking.get('order', [None])[0]
    best_idx   = int(best_model.split('_')[1])  # e.g. model_3 -> 3
except:
    best_model, best_idx = 'NA', None

# ── Extract from all 5 pkl files ─────────────────────────────────────────────
iptm_vals      = []
pae_5A_counts  = []
pae_10A_counts = []
iptm_best      = 'NA'
ptm_best       = 'NA'
ranking_best   = 'NA'

for i in range(1, 6):
    pkl_path = os.path.join(inner_dir, f'result_model_{i}_multimer_v3_pred_0.pkl')
    try:
        with open(pkl_path, 'rb') as f:
            result = pickle.load(f)

        iptm = float(result.get('iptm', float('nan')))
        ptm  = float(result.get('ptm',  float('nan')))
        rc   = float(result.get('ranking_confidence', float('nan')))
        iptm_vals.append(iptm)

        # Save best model scores
        if best_idx is not None and i == best_idx:
            iptm_best    = round(iptm, 4)
            ptm_best     = round(ptm, 4)
            ranking_best = round(rc, 4)

        # Inter-chain PAE contact counts
        pae = result.get('predicted_aligned_error', None)
        if pae is not None:
            pae = np.array(pae)
            AB = pae[0:len_A, len_A:len_A+len_B]
            BA = pae[len_A:len_A+len_B, 0:len_A]
            pae_5A_counts.append(int(np.sum(AB < 5)  + np.sum(BA < 5)))
            pae_10A_counts.append(int(np.sum(AB < 10) + np.sum(BA < 10)))
    except:
        continue

mean_iptm  = round(float(np.mean(iptm_vals)), 4)      if iptm_vals      else 'NA'
iptm_sd    = round(float(np.std(iptm_vals)), 4)       if len(iptm_vals) > 1 else 'NA'
mean_5A    = round(float(np.mean(pae_5A_counts)), 1)  if pae_5A_counts  else 'NA'
mean_10A   = round(float(np.mean(pae_10A_counts)), 1) if pae_10A_counts else 'NA'

print(f'{interactor2},{best_model},{iptm_best},{ptm_best},{ranking_best},{mean_iptm},{iptm_sd},{mean_5A},{mean_10A}')
" >> $OUTPUT_CSV

done

echo "Done! Scores written to $OUTPUT_CSV"