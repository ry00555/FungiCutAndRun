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

set -euo pipefail
META="/scratch/ry00555/RNASeqPaper/Oct2025/BAM_File_Metadata_with_index_merged_V2.csv"
MACSDIR="/scratch/ry00555/RNASeqPaper/Oct2025/MACSPeaks"
CHIPR_OUT="/scratch/ry00555/RNASeqPaper/Oct2025/ChIPR"
MIN_FRAC=0.55               # require peaks to appear in ≥55% of replicates

dos2unix "$META" 2>/dev/null || true

# --- Get column indices for Tissue and Peaks ---
header=$(head -n 1 "$META")
tissue_col=$(awk -F',' '{for(i=1;i<=NF;i++) if($i=="Tissue") print i}' <(echo "$header"))
peaks_col=$(awk -F',' '{for(i=1;i<=NF;i++) if($i=="Peaks") print i}' <(echo "$header"))

if [[ -z "$tissue_col" || -z "$peaks_col" ]]; then
  echo "❌ ERROR: Could not detect 'Tissue' or 'Peaks' columns in $META"
  exit 1
fi

# --- Loop through unique tissue groups ---
tissues=$(tail -n +2 "$META" | cut -d',' -f"$tissue_col" | sort -u)

for tissue in $tissues; do
  echo "🧬 Processing tissue: $tissue"

  # Get all peak filenames for this tissue
  mapfile -t peaks_list < <(awk -F',' -v t="$tissue" -v tc="$tissue_col" -v pc="$peaks_col" \
    'NR>1 && $tc==t {print $pc}' "$META")

  # Resolve full paths
  peak_files=()
  for peak in "${peaks_list[@]}"; do
    peak_path="${MACSDIR}/${peak}"
    if [[ -f "$peak_path" ]]; then
      peak_files+=("$peak_path")
    else
      echo "⚠️ Missing peak file: $peak_path"
    fi
  done

  num_reps=${#peak_files[@]}
  if (( num_reps < 2 )); then
    echo "⚠️ Not enough replicates ($num_reps) for $tissue — skipping."
    continue
  fi

  # --- Determine m dynamically (≥55% of replicates, rounded) ---
  m=$(awk -v n="$num_reps" 'BEGIN {
      prop = int(n * 0.55 + 0.5);
      if (prop < 2) prop = 2;
      if (prop > n) prop = n;
      print prop;
  }')

  echo "   ✅ Found $num_reps replicates → using m=$m"
  echo "   Peak files: ${peak_files[*]}"

  # --- Run ChIP-R ---
  chipr -i "${peak_files[@]}" -m "$m" -o "${CHIPR_OUT}/${tissue}_consensus"

  echo "   ✅ Consensus peaks saved to: ${CHIPR_OUT}/${tissue}_consensus.bed"
  echo
done

echo "🎯 All consensus peaks complete! Results in: $CHIPR_OUT/"
