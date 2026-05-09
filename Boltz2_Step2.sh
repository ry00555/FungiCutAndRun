#!/bin/bash
#SBATCH --job-name=ASH1_scores
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --output=/scratch/ry00555/Boltz2/ASH1/logs/Boltz2_Step2.%j.out
#SBATCH --error=/scratch/ry00555/Boltz2/ASH1/logs/Boltz2_Step2.%j.err

OUTDIR="/scratch/ry00555/Boltz2/ASH1/boltz2_outputs"
OUTPUT_CSV="$OUTDIR/ASH1_Boltz2_confidence_scores.csv"

# One row per binder, one column per model for each metric
echo "dock_protein,binder,\
confidence_score_model0,confidence_score_model1,confidence_score_model2,confidence_score_model3,confidence_score_model4,\
ptm_model0,ptm_model1,ptm_model2,ptm_model3,ptm_model4,\
iptm_model0,iptm_model1,iptm_model2,iptm_model3,iptm_model4,\
protein_iptm_model0,protein_iptm_model1,protein_iptm_model2,protein_iptm_model3,protein_iptm_model4,\
complex_plddt_model0,complex_plddt_model1,complex_plddt_model2,complex_plddt_model3,complex_plddt_model4,\
pair_iptm_AB_model0,pair_iptm_AB_model1,pair_iptm_AB_model2,pair_iptm_AB_model3,pair_iptm_AB_model4,\
mean_confidence_score,std_confidence_score,\
mean_iptm,std_iptm,\
mean_protein_iptm,std_protein_iptm,\
mean_complex_plddt,std_complex_plddt" > $OUTPUT_CSV

for pred_dir in "$OUTDIR"/boltz_results_*/predictions/*/; do

    folder=$(basename "$pred_dir")
    binder=$(echo "$folder" | sed 's/binder_[0-9]*_//' | sed 's/-t26_[0-9]*//')

    python3 << PYEOF
import json, os, statistics

pred_dir = '$pred_dir'
binder   = '$binder'
dock     = 'ASH1'
folder   = os.path.basename(pred_dir.rstrip('/'))

vals = {m: {} for m in range(5)}

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
        pair = d.get('pair_chains_iptm', {})
        vals[m]['pair_iptm_AB']     = pair.get('0', {}).get('1', 'NA')
    except:
        for k in ['confidence_score','ptm','iptm','protein_iptm','complex_plddt','pair_iptm_AB']:
            vals[m][k] = 'NA'

def col(metric):
    return ','.join(str(vals[m][metric]) for m in range(5))

def stats(metric):
    v = [vals[m][metric] for m in range(5) if isinstance(vals[m][metric], float)]
    mean = round(statistics.mean(v), 4)  if v         else 'NA'
    std  = round(statistics.stdev(v), 4) if len(v) > 1 else 'NA'
    return f'{mean},{std}'

row = (f'{dock},{binder},'
       f'{col("confidence_score")},'
       f'{col("ptm")},'
       f'{col("iptm")},'
       f'{col("protein_iptm")},'
       f'{col("complex_plddt")},'
       f'{col("pair_iptm_AB")},'
       f'{stats("confidence_score")},'
       f'{stats("iptm")},'
       f'{stats("protein_iptm")},'
       f'{stats("complex_plddt")}')
print(row)
PYEOF

done >> $OUTPUT_CSV

echo "Done! Scores written to $OUTPUT_CSV"
