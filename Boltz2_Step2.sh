#!/bin/bash
#SBATCH --job-name=Boltz2_scores
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --time=2:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --output=/scratch/ry00555/Boltz2/logs/Boltz2_Step2.%j.out
#SBATCH --error=/scratch/ry00555/Boltz2/logs/Boltz2_Step2.%j.err

# ── Only change these lines when swapping proteins ────────────────────────────
PROJECT="PRC2"
N_ANCHORS=3        # 3 for PRC2 trimer, 1 for EAF3/RTT109, 2 for dimer etc
OUTDIR="/scratch/ry00555/RNASeqPaper2026/Proteome/PRC2_Proteome_pools/Boltz2_Inputs"
KEY_CSV="/scratch/ry00555/RNASeqPaper2026/Proteome/PRC2_Proteome_pools/Boltz2_Inputs/PRC2_pool_key.csv"
# ─────────────────────────────────────────────────────────────────────────────

OUTPUT_CSV="$OUTDIR/${PROJECT}_Boltz2_confidence_scores.csv"

# Remove existing CSV so header gets written fresh
rm -f $OUTPUT_CSV

echo "Extracting Boltz2 scores for $PROJECT"
echo "N_ANCHORS : $N_ANCHORS"
echo "Output    : $OUTPUT_CSV"

for pred_dir in "$OUTDIR"/boltz_results_*/predictions/*/; do

    [ -d "$pred_dir" ] || continue

    folder=$(basename "$pred_dir")
    pool_name=$(echo "$folder" | tr '[:upper:]' '[:lower:]' | grep -oP "[a-z0-9]+_pool[0-9]+|[a-z0-9]+_trimer" | head -1)
    [ -z "$pool_name" ] && pool_name="$folder"

    python3 << PYEOF
import json, os, statistics, csv

pred_dir   = '$pred_dir'
folder     = os.path.basename(pred_dir.rstrip('/'))
pool_name  = '$pool_name'
dock       = '$PROJECT'
n_anchors  = $N_ANCHORS
output_csv = '$OUTPUT_CSV'
key_csv    = '$KEY_CSV'

# ── Read chain names from key CSV ─────────────────────────────────────────────
anchor_names    = []
candidate_names = []

try:
    with open(key_csv) as f:
        for line in f:
            parts = line.strip().split(',')
            if parts[0].lower() == pool_name.lower():
                anchor_names    = parts[1:1 + n_anchors]
                candidate_names = parts[1 + n_anchors:]
                break
except Exception as e:
    print(f'Warning: could not read key CSV: {e}')

# ── Read all 5 models ─────────────────────────────────────────────────────────
vals    = {m: {} for m in range(5)}
n_chains = 0

for m in range(5):
    cf = os.path.join(pred_dir, f'confidence_{folder}_model_{m}.json')
    try:
        with open(cf) as f:
            d = json.load(f)
        vals[m]['confidence_score'] = d.get('confidence_score', 'NA')
        vals[m]['ptm']              = d.get('ptm',              'NA')
        vals[m]['iptm']             = d.get('iptm',             'NA')
        vals[m]['protein_iptm']     = d.get('protein_iptm',     'NA')
        vals[m]['complex_plddt']    = d.get('complex_plddt',    'NA')
        vals[m]['pair_chains_iptm'] = d.get('pair_chains_iptm', {})
        if n_chains == 0 and vals[m]['pair_chains_iptm']:
            n_chains = len(vals[m]['pair_chains_iptm'])
    except:
        for k in ['confidence_score', 'ptm', 'iptm', 'protein_iptm', 'complex_plddt']:
            vals[m][k] = 'NA'
        vals[m]['pair_chains_iptm'] = {}

# Fallback chain names if key CSV lookup failed
if not anchor_names:
    anchor_names    = [f'anchor_{i}' for i in range(n_anchors)]
if not candidate_names:
    n_cands = max(0, n_chains - n_anchors)
    candidate_names = [f'candidate_{i}' for i in range(n_cands)]

anchor_indices    = list(range(n_anchors))
candidate_indices = list(range(n_anchors, n_chains)) if n_chains > n_anchors else []

# ── Helper functions ──────────────────────────────────────────────────────────
def col(metric):
    return [vals[m][metric] for m in range(5)]

def mean_std(values):
    v = [x for x in values if isinstance(x, (int, float))]
    mean = round(statistics.mean(v), 4)  if v          else 'NA'
    std  = round(statistics.stdev(v), 4) if len(v) > 1 else 'NA'
    return mean, std

def pair_iptm_scores(chain_i, chain_j):
    scores = []
    for m in range(5):
        try:
            score = vals[m]['pair_chains_iptm'].get(str(chain_i), {}).get(str(chain_j), 'NA')
        except:
            score = 'NA'
        scores.append(score)
    return scores

# ── Build header ──────────────────────────────────────────────────────────────
write_header = not os.path.exists(output_csv)

header = [
    'dock', 'pool_name', 'candidate_name', 'candidate_chain_idx',
    'confidence_score_m0','confidence_score_m1','confidence_score_m2',
    'confidence_score_m3','confidence_score_m4',
    'mean_confidence_score','std_confidence_score',
    'iptm_m0','iptm_m1','iptm_m2','iptm_m3','iptm_m4',
    'mean_iptm','std_iptm',
    'protein_iptm_m0','protein_iptm_m1','protein_iptm_m2',
    'protein_iptm_m3','protein_iptm_m4',
    'mean_protein_iptm','std_protein_iptm',
    'complex_plddt_m0','complex_plddt_m1','complex_plddt_m2',
    'complex_plddt_m3','complex_plddt_m4',
    'mean_complex_plddt','std_complex_plddt',
]

# One set of pair_iptm columns per anchor
for anchor_name in anchor_names:
    for suffix in ['m0','m1','m2','m3','m4','mean','std']:
        header.append(f'pair_iptm_{anchor_name}_vs_cand_{suffix}')

# ── Write rows ────────────────────────────────────────────────────────────────
with open(output_csv, 'a', newline='') as csvfile:
    writer = csv.writer(csvfile)

    if write_header:
        writer.writerow(header)

    # If no candidates (trimer-only control) write one summary row
    if not candidate_indices:
        cs     = col('confidence_score')
        iptm   = col('iptm')
        piptm  = col('protein_iptm')
        cplddt = col('complex_plddt')
        m_cs,  s_cs    = mean_std(cs)
        m_i,   s_i     = mean_std(iptm)
        m_pi,  s_pi    = mean_std(piptm)
        m_cp,  s_cp    = mean_std(cplddt)

        row = ([dock, pool_name, 'trimer_only', 'NA'] +
               cs + [m_cs, s_cs] +
               iptm + [m_i, s_i] +
               piptm + [m_pi, s_pi] +
               cplddt + [m_cp, s_cp])

        # No candidates so pair_iptm all NA
        for anchor_name in anchor_names:
            row += ['NA']*5 + ['NA', 'NA']

        writer.writerow(row)

    else:
        # One row per candidate chain
        for cand_idx in candidate_indices:
            cand_pos  = cand_idx - n_anchors
            cand_name = (candidate_names[cand_pos]
                         if cand_pos < len(candidate_names)
                         else f'candidate_{cand_pos}')

            cs     = col('confidence_score')
            iptm   = col('iptm')
            piptm  = col('protein_iptm')
            cplddt = col('complex_plddt')
            m_cs,  s_cs = mean_std(cs)
            m_i,   s_i  = mean_std(iptm)
            m_pi,  s_pi = mean_std(piptm)
            m_cp,  s_cp = mean_std(cplddt)

            row = ([dock, pool_name, cand_name, cand_idx] +
                   cs + [m_cs, s_cs] +
                   iptm + [m_i, s_i] +
                   piptm + [m_pi, s_pi] +
                   cplddt + [m_cp, s_cp])

            # pair_iptm for each anchor vs this candidate
            for anchor_idx, anchor_name in zip(anchor_indices, anchor_names):
                scores          = pair_iptm_scores(anchor_idx, cand_idx)
                mean_p, std_p   = mean_std(scores)
                row += scores + [mean_p, std_p]

            writer.writerow(row)

print(f'Written: {pool_name} — {len(candidate_indices) or 1} row(s)')
PYEOF

done

echo "Done! Scores written to $OUTPUT_CSV"
