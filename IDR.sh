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
CHIPR_DIR="/scratch/ry00555/RNASeqPaper/Oct2025/ChIPR"
MIN_FRAC=0.55               # require peaks to appear in ≥55% of replicates
IDR="/scratch/ry00555/RNASeqPaper/Oct2025/IDR"
SUMMARY="${IDR}/replicate_vs_CHIPRconsensus.tsv"

# Remove carriage returns
dos2unix "$META" 2>/dev/null || true
echo -e "Tissue\tReplicate\tPeakFile\tNumPeaks\tNumOverlapWithConsensus\tFractionOverlap\tIDR_Passing" > "$SUMMARY"

# Output TSV
tail -n +2 "$META" | while IFS=, read -r RunID bamReads BamIndex SampleID Factor Tissue Condition Replicate bamControl bamInputIndex ControlID Peaks PeakCaller DesiredPeakName; do
    [[ -z "$Tissue" ]] && continue

    rep_peak="${MACSDIR}/${Peaks}"
    consensus="${CHIPR_OUT}/${Tissue}_consensus_optimal.bed"

    if [[ ! -s "$rep_peak" ]]; then
        echo "⚠️ Missing replicate peak: $rep_peak"
        continue
    fi
    if [[ ! -s "$consensus" ]]; then
        echo "⚠️ Missing consensus peak: $consensus"
        continue
    fi

    num_peaks=$(wc -l < "$rep_peak")

    # === BEDTOOLS overlap fraction ===
    overlap_file="${IDR}/${SampleID}_vs_consensus_overlap.bed"
    bedtools intersect -u -a "$rep_peak" -b "$consensus" > "$overlap_file"
    num_overlap=$(wc -l < "$overlap_file")
    frac=$(awk "BEGIN {printf \"%.3f\", ($num_peaks>0)?$num_overlap/$num_peaks:0}")

    # === True IDR Analysis ===
    idr_output="${IDR}/${SampleID}_vs_consensus.idr.txt"
    idr_plot="${IDR}/${SampleID}_vs_consensus.png"

    idr --samples "$rep_peak" "$consensus" \
        --input-file-type broadPeak \
        --rank p.value \
        --output-file "$idr_output" \
        --plot \
        --log-output-file "${IDR}/${SampleID}_vs_consensus.log" 2>/dev/null

    idr_pass=$(awk '($12 < 0.05){c++} END{print c+0}' "$idr_output" 2>/dev/null || echo 0)

    echo -e "${Tissue}\t${Replicate}\t${rep_peak}\t${num_peaks}\t${num_overlap}\t${frac}\t${idr_pass}" >> "$SUMMARY"

done
