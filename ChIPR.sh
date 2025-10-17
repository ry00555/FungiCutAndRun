#!/bin/bash
#SBATCH --job-name=MacsPeakCalling
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=90gb
#SBATCH --time=8:00:00
#SBATCH --output=../ChIP-R.%j.out
#SBATCH --error=../ChIP-R.%j.err


META="/scratch/ry00555/RNASeqPaper/Oct2025/BAM_File_Metadata_with_index_merged_V2.csv"
MACSDIR="/scratch/ry00555/RNASeqPaper/Oct2025/MACSPeaks"
CHIPR_OUT="/scratch/ry00555/RNASeqPaper/Oct2025/ChIPR"


module load ChIP-R/1.1.0-foss-2023a-Python-3.11.3
dos2unix "$META" 2>/dev/null || true

header=$(head -n 1 "$META")

# Find column numbers for Tissue and DesiredPeakName
tissue_col=$(awk -F',' '{for(i=1;i<=NF;i++){if($i=="Tissue") print i}}' <(echo "$header"))
desired_col=$(awk -F',' '{for(i=1;i<=NF;i++){if($i=="DesiredPeakName") print i}}' <(echo "$header"))

if [[ -z "$tissue_col" || -z "$desired_col" ]]; then
    echo "❌ Could not find 'Tissue' or 'DesiredPeakName' columns in $META"
    exit 1
fi

echo "📊 Tissue column: $tissue_col | DesiredPeakName column: $desired_col"

# --- Get unique tissue names ---
tissues=$(tail -n +2 "$META" | cut -d',' -f"$tissue_col" | sort -u)

# --- Loop over each tissue group ---
for tissue in $tissues; do
    echo "➡️ Processing group: $tissue"

    # Collect all DesiredPeakNames for this tissue
    peak_names=$(awk -F',' -v t="$tissue" -v tc="$tissue_col" -v dc="$desired_col" \
        'NR>1 && $tc==t {print $dc}' "$META")

    peak_files=()
    for name in $peak_names; do
        file="${MACSDIR}/${name}_peaks.broadPeak"
        if [[ -f "$file" ]]; then
            peak_files+=("$file")
        else
            echo "⚠️ Missing file: $file"
        fi
    done

    n=${#peak_files[@]}
    if (( n < 2 )); then
        echo "❌ Not enough replicates for $tissue (found $n)"
        continue
    fi

    # --- Determine minentries (m) based on 55% of replicates ---
    m=$(awk -v n="$n" 'BEGIN{printf "%d\n", int(n*0.55 + 0.999)}')  # round up
    if (( m < 1 )); then m=1; fi
    echo "   🧮 $n replicates → using --minentries $m (55% of replicates)"

    # --- Run ChIP-R ---
    echo "   🧬 Running ChIP-R..."
    chipr -i "${peak_files[@]}" -m "$m" -o "${CHIPR_OUT}/${tissue}_consensus"

done

echo "✅ All ChIP-R consensus peak sets written to: $CHIPR_OUT"
