#!/bin/bash
#SBATCH --job-name=ChIP-R
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=30gb
#SBATCH --time=2:00:00
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

echo -e "Tissue\tNumReps\tm_value\tNumPeaks\tOutputFile\tValidPeakFiles\tSkippedPeakFiles" > "$SUMMARY"

# --- STEP 1: Extract Tissue + Peaks ---
# Adjust field positions if needed (Tissue=$6, Peaks=$12)
tail -n +2 "$META" | awk -F',' '{print $6 "\t" $12}' > tmp_peaks_by_tissue.tsv

# --- STEP 2: Loop over tissues ---
for tissue in $(cut -f1 tmp_peaks_by_tissue.tsv | sort | uniq); do
    echo "🔍 Processing tissue: $tissue"

    peak_files=()
    skipped_files=()

    # Collect valid .broadPeak files
    while IFS=$'\t' read -r t peak; do
        [[ "$t" != "$tissue" ]] && continue

        # Prepend full path if needed
        [[ ! "$peak" =~ / ]] && peak="${MACSDIR}/${peak}"

        if [[ -f "$peak" && -s "$peak" ]]; then
            peak_files+=("$peak")
        else
            echo "⚠️ Skipping missing or empty peak file: $peak"
            skipped_files+=("$peak")
        fi
    done < tmp_peaks_by_tissue.tsv

    n=${#peak_files[@]}
    if (( n < 2 )); then
        echo "⚠️ Skipping $tissue: fewer than 2 valid replicates ($n)"
        continue
    fi

    # --- Compute m = max(2, floor(MIN_FRAC * n)) ---
    m=$(echo "scale=0; val=int($n * $MIN_FRAC); if(val<2) val=2; val" | bc)

    prefix="${CHIPR_OUT}/${tissue}_consensus"
    outfile="${prefix}_optimal.bed"

    # --- Skip if consensus already exists ---
    if [[ -s "$outfile" ]]; then
        num_peaks=$(wc -l < "$outfile")
        echo "✅ Skipping $tissue (consensus exists, $num_peaks peaks)"
        echo -e "${tissue}\t${n}\t${m}\t${num_peaks}\t${outfile}\t${peak_files[*]}\t${skipped_files[*]}" >> "$SUMMARY"
        continue
    fi

    # --- Run ChIPR ---
    echo "🧬 Running ChIPR for $tissue ($n replicates, m=$m)"
    echo "   Valid peak files: ${peak_files[*]}"
    chipr -i "${peak_files[@]}" -m "$m" -o "$prefix"

    # --- Record results ---
    if [[ -f "$outfile" && -s "$outfile" ]]; then
        num_peaks=$(wc -l < "$outfile")
        echo "   📊 $num_peaks consensus peaks retained for $tissue"
    else
        num_peaks=0
        echo "   ⚠️ No output found for $tissue"
    fi

    echo -e "${tissue}\t${n}\t${m}\t${num_peaks}\t${outfile}\t${peak_files[*]}\t${skipped_files[*]}" >> "$SUMMARY"
done

rm -f tmp_peaks_by_tissue.tsv
echo "✅ ChIPR consensus generation complete!"
echo "📄 Summary written to: $SUMMARY"
