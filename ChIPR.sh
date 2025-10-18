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

dos2unix "$META" 2>/dev/null || true
ml ChIP-R/1.1.0-foss-2023a-Python-3.11.3

# --- Get column indices for Tissue and Peaks ---
header=$(head -n 1 "$META")
tissue_col=$(awk -F',' '{for(i=1;i<=NF;i++) if($i=="Tissue") print i}' <(echo "$header"))
peaks_col=$(awk -F',' '{for(i=1;i<=NF;i++) if($i=="Peaks") print i}' <(echo "$header"))

SUMMARY="${OUTDIR}/consensus_summary.tsv"
echo -e "Tissue\tNumReps\tm_value\tOutputFile\tPeakFiles" > "$SUMMARY"

# === LOOP OVER TISSUES ===
# Skip header, group by Tissue, and process each group
tail -n +2 "$META" | awk -F',' '{print $0}' | while IFS=',' read -r SampleID Tissue Peaks; do
    echo -e "${Tissue}\t${Peaks}" >> tmp_peaks_by_tissue.tsv
done

# Get unique tissues
for tissue in $(cut -f1 tmp_peaks_by_tissue.tsv | sort | uniq); do
    # Get peak files for this tissue
    peak_files=($(awk -v t="$tissue" '$1 == t {print $2}' tmp_peaks_by_tissue.tsv))
    n=${#peak_files[@]}

    if (( n < 2 )); then
        echo "Skipping $tissue (only $n replicate(s))"
        continue
    fi

    # Calculate m = max(2, floor(0.55 * n))
    m=$(echo "scale=0; n=$n; val=(n*0.55)/1; if (val<2) val=2; val" | bc)
    echo "Processing $tissue with $n replicates (m=$m)"

    # Output file prefix
    prefix="${OUTDIR}/${tissue}_consensus"

    # Run ChIPR
    chipr -i "${peak_files[@]}" -m "$m" -o "$prefix"

    # Append summary info
    echo -e "${tissue}\t${n}\t${m}\t${prefix}.bed\t${peak_files[*]}" >> "$SUMMARY"
done

# Cleanup temp file
rm -f tmp_peaks_by_tissue.tsv

echo "✅ Consensus generation complete. Summary written to: $SUMMARY"
