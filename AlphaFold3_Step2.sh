#!/bin/bash
#SBATCH --job-name=EAF3
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --mem=8gb
#SBATCH --time=1:00:00
#SBATCH --output=AF3_Extraction.%j.out
#SBATCH --error=AF3_Extraction.%j.err

# ── Only change these lines when swapping proteins ────────────────────────────
PROTEIN="EAF3"
SCREEN_TYPE="pooled"    # "pooled" or "pairwise"
ROW_INDEX=1             # 1 = chain A (bait), 2 = chain B, etc.
# ─────────────────────────────────────────────────────────────────────────────

SCRIPTS_DIR="/home/ry00555/Research/FungiCutAndRun"
BASE_DIR="/scratch/ry00555/EpigeneticMemoryPaper2026/AlphaFold3/PooledPPI"

if [ "$SCREEN_TYPE" == "pooled" ]; then
    # Pooled screen — uses AF3pool_scraper_v1.py and key CSV made from generate_AF3_pooled_jsons.py 
    KEY_CSV="${BASE_DIR}/${PROTEIN}_AF3_PooledJSONs/${PROTEIN}_pool_key_scraper.csv"
    OUTPUT_DIR="${BASE_DIR}/${PROTEIN}_AF3_PooledJSON_output"
    OUTPUT_CSV="${OUTPUT_DIR}/${PROTEIN}_AF3_ipTM_scores_Pool.csv"

    python3 ${SCRIPTS_DIR}/AF3pool_scraper_v1.py \
        $KEY_CSV \
        $OUTPUT_DIR \
        $OUTPUT_CSV \
        --row-index $ROW_INDEX

elif [ "$SCREEN_TYPE" == "pairwise" ]; then
    # Pairwise one-vs-all screen — uses inline python
    AF3_DIR="/scratch/ry00555/AlphaFold/output/AF3"
    OUTPUT_CSV="$AF3_DIR/${PROTEIN}_AF3_ipTM_scores_Pairwise.csv"
    dock_lower=$(echo $PROTEIN | tr '[:upper:]' '[:lower:]')

    echo "dock_protein,interactor2,iptm_best,ranking_score_best,mean_iptm,std_iptm,mean_ranking_score,std_ranking_score,chain_pair_iptm_AB,chain_pair_pae_min_AB" > $OUTPUT_CSV

    for dir in $AF3_DIR/*/; do
        job_name=$(basename $dir)
        interactor2=$(echo $job_name | sed "s/${dock_lower}//")

        python3 << PYEOF
import json, os, statistics

job_dir      = '$dir'
job_name     = '$job_name'
interactor2  = '$interactor2'
dock_protein = '$PROTEIN'

try:
    with open(os.path.join(job_dir, f'{job_name}_summary_confidences.json')) as f:
        top = json.load(f)
    iptm_best          = top.get('iptm', 'NA')
    ranking_score_best = top.get('ranking_score', 'NA')
    cpp = top.get('chain_pair_iptm', None)
    chain_pair_iptm_AB = cpp[0][1] if cpp else 'NA'
    cpm = top.get('chain_pair_pae_min', None)
    chain_pair_pae_AB  = cpm[0][1] if cpm else 'NA'
except Exception:
    iptm_best, ranking_score_best, chain_pair_iptm_AB, chain_pair_pae_AB = 'NA', 'NA', 'NA', 'NA'

iptm_vals, ranking_score_vals = [], []
for i in range(5):
    sample_dir = os.path.join(job_dir, f'seed-1_sample-{i}')
    try:
        with open(os.path.join(sample_dir, 'summary_confidences.json')) as f:
            data = json.load(f)
        iptm = data.get('iptm', None)
        rs   = data.get('ranking_score', None)
        if iptm is not None: iptm_vals.append(iptm)
        if rs   is not None: ranking_score_vals.append(rs)
    except:
        continue

mean_iptm          = round(statistics.mean(iptm_vals), 4)           if iptm_vals                   else 'NA'
std_iptm           = round(statistics.stdev(iptm_vals), 4)          if len(iptm_vals) > 1          else 'NA'
mean_ranking_score = round(statistics.mean(ranking_score_vals), 4)  if ranking_score_vals          else 'NA'
std_ranking_score  = round(statistics.stdev(ranking_score_vals), 4) if len(ranking_score_vals) > 1 else 'NA'

print(f'{dock_protein},{interactor2},{iptm_best},{ranking_score_best},{mean_iptm},{std_iptm},{mean_ranking_score},{std_ranking_score},{chain_pair_iptm_AB},{chain_pair_pae_AB}')
PYEOF

    done >> $OUTPUT_CSV
fi

echo "Done! Scores written to $OUTPUT_CSV"
