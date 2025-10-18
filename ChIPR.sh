#!/bin/bash
#SBATCH --job-name=MacsPeakCalling
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=8:00:00
#SBATCH --output=../ChIP-R.%j.out
#SBATCH --error=../ChIP-R.%j.err

set -euo pipefail
META="/scratch/ry00555/RNASeqPaper/Oct2025/BAM_File_Metadata_with_index_merged_V2.csv"
MACSDIR="/scratch/ry00555/RNASeqPaper/Oct2025/MACSPeaks"
CHIPR_OUT="/scratch/ry00555/RNASeqPaper/Oct2025/ChIPR"
MIN_FRAC=0.55               # require peaks to appear in ≥55% of replicates
SUMMARY="${CHIPR_OUT}/consensus_summary.tsv"

dos2unix "$META" 2>/dev/null || true
ml ChIP-R/1.1.0-foss-2023a-Python-3.11.3

echo -e "Tissue\tNumReps\tm_value\tOutputFile\tPeakFiles" > "$SUMMARY"

# === STEP 1: Collect Tissue + Peaks columns ===
# Expecting columns: RunID,...,Tissue,...,Peaks,...
# Adjust field positions if needed (Tissue=$?, Peaks=$?)
tail -n +2 "$META" | awk -F',' '{print $6 "\t" $12}' > tmp_peaks_by_tissue.tsv

# === STEP 2: Loop over each unique tissue ===
for tissue in $(cut -f1 tmp_peaks_by_tissue.tsv | sort | uniq); do
    echo "🔍 Processing tissue: $tissue"
    peak_files=()

    # Gather all .broadPeak files for this tissue
    while IFS=$'\t' read -r t peak; do
        [[ "$t" != "$tissue" ]] && continue

        # Add full path if only filename is provided
        if [[ ! "$peak" =~ / ]]; then
            peak="${MACSDIR}/${peak}"
        fi

        # Confirm existence
        if [[ -f "$peak" ]]; then
            peak_files+=("$peak")
        else
            echo "⚠️ Missing peak file: $peak"
        fi
    done < tmp_peaks_by_tissue.tsv

    n=${#peak_files[@]}
    if (( n < 2 )); then
        echo "⚠️ Skipping $tissue (only $n valid replicates)"
        continue
    fi

    # === Compute m = max(2, floor(0.55 * n)) ===
    m=$(echo "scale=0; n=$n; val=(n*0.55)/1; if (val<2) val=2; val" | bc)

    echo "🧬 Running ChIPR for $tissue ($n replicates, m=$m)"
    echo "   Peak files: ${peak_files[*]}"

    prefix="${CHIPR_OUT}/${tissue}_consensus"

    # === Run ChIPR ===
    chipr -i "${peak_files[@]}" -m "$m" -o "$prefix"

    echo -e "${tissue}\t${n}\t${m}\t${prefix}.bed\t${peak_files[*]}" >> "$SUMMARY"
done

rm -f tmp_peaks_by_tissue.tsv

echo "✅ ChIPR consensus generation complete!"
echo "📄 Summary written to: $SUMMARY"
