#!/bin/bash
#SBATCH --job-name=IDR
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=8:00:00
#SBATCH --output=../IDR.%j.out
#SBATCH --error=../IDR.%j.err

#ml R/4.4.2-gfbf-2024a
ml IDR/2.0.3-foss-2023a BEDTools/2.31.1-GCC-13.3.0

META="/scratch/ry00555/RNASeqPaper/Oct2025/BAM_File_Metadata_with_index_merged_V2.csv"
MACSDIR="/scratch/ry00555/RNASeqPaper/Oct2025/MACSPeaks"
CHIPR_OUT="/scratch/ry00555/RNASeqPaper/Oct2025/ChIPR"
MIN_FRAC=0.55               # require peaks to appear in ≥55% of replicates
SUMMARY="${CHIPR_OUT}/consensus_summary.tsv"
IDR="/scratch/ry00555/RNASeqPaper/Oct2025/IDR"
# Remove carriage returns
dos2unix "$META" 2>/dev/null || true

# Output TSV
SUMMARY="${IDR}/replicate_vs_CHIPRconsensus.tsv"
echo -e "Tissue\tReplicate\tPeakFile\tNumPeaks\tNumOverlapWithConsensus\tFractionOverlap" > "$SUMMARY"

# --- Process each tissue ---
tail -n +2 "$META" | while IFS=, read -r RunID bamReads BamIndex SampleID Factor Tissue Condition Replicate bamControl bamInputIndex ControlID Peaks PeakCaller DesiredPeakName; do
    [[ -z "$Tissue" ]] && continue

    # Paths
    consensus="${CHIPR_DIR}/${Tissue}_consensus_optimal.bed"
    rep_peak="${MACSDIR}/${Peaks}"

    if [[ ! -f "$rep_peak" ]]; then
        echo "⚠️ Missing replicate peak: $rep_peak"
        continue
    fi
    if [[ ! -f "$consensus" ]]; then
        echo "⚠️ Missing consensus peak: $consensus"
        continue
    fi

    # Count total peaks in replicate
    num_peaks=$(wc -l < "$rep_peak")

    # Intersect with consensus
    overlap_file="${IDR}/${SampleID}_vs_consensus_optimal.bed"
    bedtools intersect -u -a "$rep_peak" -b "$consensus" > "$overlap_file"
    num_overlap=$(wc -l < "$overlap_file")

    # Fraction overlapping
    if (( num_peaks > 0 )); then
        frac=$(awk "BEGIN {printf \"%.3f\", $num_overlap/$num_peaks}")
    else
        frac=0
    fi

    echo -e "${Tissue}\t${Replicate}\t${rep_peak}\t${num_peaks}\t${num_overlap}\t${frac}" >> "$SUMMARY"
done

echo "✅ Ranking complete! See results in: $SUMMARY"
