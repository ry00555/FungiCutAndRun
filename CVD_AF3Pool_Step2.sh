#!/bin/bash
#SBATCH --job-name=PA14_52800
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --mem=8gb
#SBATCH --time=1:00:00
#SBATCH --output=AF3_Step2_PullipTMs.%j.out
#SBATCH --error=AF3_Step2_PullipTMs.%j.err

WORKDIR="/scratch/ry00555/CVD_AF3"
AF3_DIR="$WORKDIR/output"
OUTPUT_CSV="$AF3_DIR/CVD_AF3_ipTM_scores.csv"

echo "dock_protein,interactor2,iptm_sample0,ranking_score_sample0,mean_iptm,std_iptm,mean_ranking_score,std_ranking_score,chain_pair_iptm_AB,chain_pair_pae_min_AB" > $OUTPUT_CSV

for dir in $AF3_DIR/*/; do
    job_name=$(basename $dir)

    # Extract interactor2 name from job name
    # job names are like pa14_52800pa14_41980biocycmanual
    # interactor2 is everything after pa14_52800
    interactor2=$(echo $job_name | sed 's/pa14_52800//')

    python3 << PYEOF
import json, os
import statistics

job_dir = '$dir'
job_name = '$job_name'
interactor2 = '$interactor2'
dock_protein = 'PA14_52800'
    
# Read sample 0 (best model) from top-level summary_confidences.json
try:
    with open(os.path.join(job_dir, f'{job_name}_summary_confidences.json')) as f:
        top = json.load(f)
    iptm_sample0        = top.get('iptm', 'NA')
    ranking_score_sample0 = top.get('ranking_score', 'NA')
    # chain_pair_iptm[0][1] = A->B inter-chain ipTM
    cpp = top.get('chain_pair_iptm', None)
    chain_pair_iptm_AB  = cpp[0][1] if cpp else 'NA'
    # chain_pair_pae_min[0][1] = A->B minimum PAE
    cpm = top.get('chain_pair_pae_min', None)
    chain_pair_pae_AB   = cpm[0][1] if cpm else 'NA'
except Exception as e:
    iptm_sample0, ranking_score_sample0, chain_pair_iptm_AB, chain_pair_pae_AB = 'NA', 'NA', 'NA', 'NA'

# Read all 5 samples for mean ipTM and ranking score
iptm_vals         = []
ranking_score_vals = []

for i in range(5):
    sample_dir = os.path.join(job_dir, f'seed-1_sample-{i}')
    try:
        with open(os.path.join(sample_dir, 'summary_confidences.json')) as f:
            data = json.load(f)
        iptm = data.get('iptm', None)
        rs   = data.get('ranking_score', None)
        if iptm is not None:
            iptm_vals.append(iptm)
        if rs is not None:
            ranking_score_vals.append(rs)
    except:
        continue

mean_iptm          = round(statistics.mean(iptm_vals), 4)          if iptm_vals          else 'NA'
mean_ranking_score = round(statistics.mean(ranking_score_vals), 4) if ranking_score_vals else 'NA'
std_iptm = round(statistics.stdev(iptm_vals), 4) if len(iptm_vals) > 1 else 'NA'
std_ranking_score = round(statistics.stdev(ranking_score_vals), 4) if len(ranking_score_vals) > 1 else 'NA'
print(f'{dock_protein},{interactor2},{iptm_sample0},{ranking_score_sample0},{mean_iptm},{std_iptm},{mean_ranking_score},{std_ranking_score},{chain_pair_iptm_AB},{chain_pair_pae_AB}')
PYEOF

done >> $OUTPUT_CSV

echo "Done! Scores written to $OUTPUT_CSV"