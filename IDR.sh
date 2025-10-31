#!/bin/bash
#SBATCH --job-name=IDR
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=400gb
#SBATCH --time=4:00:00
#SBATCH --output=../IDR.%j.out
#SBATCH --error=../IDR.%j.err

set -euo pipefail

# ================================
# Load modules
# ================================
ml BEDTools/2.30.0-GCC-11.3.0 deepTools SAMtools/1.16.1-GCC-11.3.0 BamTools/2.5.2-GCC-11.3.0

# ================================
# Paths
# ================================
META="/scratch/ry00555/RNASeqPaper/Oct2025/BAM_File_Metadata_with_index_merged_V2.csv"
MACSDIR="/scratch/ry00555/RNASeqPaper/Oct2025/MACSPeaks"
CHIPR_DIR="/scratch/ry00555/RNASeqPaper/Oct2025/ChIPR"
BAMDIR="/scratch/ry00555/RNASeqPaper/Oct2025/SortedBamFiles"
OUTDIR="/scratch/ry00555/RNASeqPaper/Oct2025/IDR"

MASTER_SUMMARY="${OUTDIR}/master_summary.tsv"
FRIP_TSV="${OUTDIR}/frip_metrics.tsv"
OVERLAP_TSV="${OUTDIR}/replicate_overlap.tsv"
JACCARD_TSV="${OUTDIR}/pairwise_jaccard.tsv"
BAM_CORR_NPZ="${OUTDIR}/bam_corr.npz"
CORR_HEAT="${OUTDIR}/bam_correlation_heatmap.pdf"

mkdir -p "$OUTDIR"

# ================================
# Initialize master summary
# ================================
#echo -e "SampleID\tTissue\tFactor\tTotalReads\tReadsInPeaks\tFRiP\tPeakFile\tConsensusPeak\tNumPeaks\tNumOverlap\tFracOverlap" > "$MASTER_SUMMARY"
#echo -e "SampleID\tTissue\tFactor\tTotalReads\tReadsInPeaks\tFRiP" > "$FRIP_TSV"
#echo -e "Tissue\tSampleID\tPeakFile\tNumPeaks\tNumOverlap\tFracOverlap" > "$OVERLAP_TSV"

# ================================
# Index BAMs only if needed
# ================================
# echo "---- Checking and reindexing BAMs ----"
# for bam in "${BAMDIR}"/*.bam; do
#     bai="${bam}.bai"
#     if [[ ! -f "$bai" || "$bam" -nt "$bai" ]]; then
#         echo "Indexing BAM: $bam"
#         samtools index -@ 4 "$bam"
#     else
#         echo "BAM index up-to-date: $bai"
#     fi
# done

# ================================
# FRiP and peak overlap
# ================================
echo "---- Calculating FRiP and overlap ----"
#############################################
# Reprocess only the last N processed samples
#############################################

REPROCESS_LAST_N=3

# If MASTER_SUMMARY doesn't exist or is empty, create empty lists
if [[ ! -s "$MASTER_SUMMARY" ]]; then
    echo "" > /tmp/processed_ids.txt
    echo "" > /tmp/lastN_ids.txt
else
    # All processed SampleIDs
    awk '{print $1}' "$MASTER_SUMMARY" > /tmp/processed_ids.txt

    # Last N SampleIDs
    tail -n "$REPROCESS_LAST_N" "$MASTER_SUMMARY" | awk '{print $1}' > /tmp/lastN_ids.txt
fi

#############################################
# Main metadata loop (with skip logic)
#############################################

tail -n +2 "$META" | while IFS=, read -r RunID bamReads BamIndex SampleID Factor Tissue Condition Replicate bamControl bamInputIndex ControlID Peaks PeakCaller DesiredPeakName MACS3minlength MACS3maxgap; do
    [[ -z "$SampleID" || -z "$Tissue" || -z "$bamReads" ]] && continue

    # Skip if already processed AND not one of the re-run entries
    if grep -qx "$SampleID" /tmp/processed_ids.txt && ! grep -qx "$SampleID" /tmp/lastN_ids.txt; then
        echo "Skipping (already processed): $SampleID"
        continue
    fi

    echo "Reprocessing: $SampleID"

    bam="${BAMDIR}/${bamReads}"
    peak="${MACSDIR}/${DesiredPeakName}_peaks.broadPeak"
    consensus="${CHIPR_DIR}/${Tissue}_consensus_optimal_peaks.broadPeak"

    [[ ! -s "$bam" ]] && { echo "Missing BAM → skipping $SampleID"; continue; }
    [[ ! -s "$peak" ]] && { echo "Missing peak file → skipping $SampleID"; continue; }

    # FRiP calculation
    total_reads=$(samtools view -c -F 260 "$bam" 2>/dev/null || echo 0)
    reads_in_peaks=$(bedtools intersect -a "$peak" -b "$bam" -c | awk '{sum+=$NF} END{print sum+0}')
    frip=$(awk "BEGIN{printf \"%.4f\", ($reads_in_peaks/$total_reads)}")

    # Peak overlap
    num_peaks=$(wc -l < "$peak" | tr -d '[:space:]')
    if [[ -s "$consensus" && "$num_peaks" -gt 0 ]]; then
        num_overlap=$(bedtools intersect -u -a "$peak" -b "$consensus" | wc -l | tr -d '[:space:]')
        frac_overlap=$(awk "BEGIN{printf \"%.4f\", ($num_overlap/$num_peaks)}")
    else
        num_overlap=0
        frac_overlap=0
    fi

    # Append to outputs
    echo -e "${SampleID}\t${Tissue}\t${Factor}\t${total_reads}\t${reads_in_peaks}\t${frip}\t${peak}\t${consensus}\t${num_peaks}\t${num_overlap}\t${frac_overlap}" >> "$MASTER_SUMMARY"
    echo -e "${SampleID}\t${Tissue}\t${Factor}\t${total_reads}\t${reads_in_peaks}\t${frip}" >> "$FRIP_TSV"
    echo -e "${Tissue}\t${SampleID}\t${peak}\t${num_peaks}\t${num_overlap}\t${frac_overlap}" >> "$OVERLAP_TSV"
done



# ================================
# JACCARD SIMILARITY
# ================================
# echo "---- Calculating Jaccard similarity ----"
# echo -e "fileA\tfileB\tjaccard" > "$JACCARD_TSV"
#
# PEAK_FILES=( $(awk -F'\t' 'NR>1{print $7}' "$MASTER_SUMMARY" | sort -u) )
# for ((i=0;i<${#PEAK_FILES[@]};i++)); do
#     for ((j=i+1;j<${#PEAK_FILES[@]};j++)); do
#         f1="${PEAK_FILES[i]}"
#         f2="${PEAK_FILES[j]}"
#         if [[ -s "$f1" && -s "$f2" ]]; then
#             jacc=$(bedtools jaccard -a "$f1" -b "$f2" | awk 'NR==2{print $3}')
#             echo -e "$(basename "$f1")\t$(basename "$f2")\t${jacc:-0}" >> "$JACCARD_TSV"
#         fi
#     done
# done

# ================================
# MULTIBAMSUMMARY + PLOT CORRELATION
# ================================
BAMLIST="${OUTDIR}/bamlist.txt"
> "$BAMLIST"

tail -n +2 "$META" | while IFS=, read -r RunID bamReads BamIndex SampleID Factor Tissue Condition Replicate bamControl bamInputIndex ControlID Peaks PeakCaller DesiredPeakName MACS3minlength MACS3maxgap; do
    bam="${BAMDIR}/${bamReads}"
    [[ -s "$bam" ]] && echo "$bam" >> "$BAMLIST"
done

if [[ -s "$BAMLIST" ]]; then
    echo "---- Running deeptools correlation ----"
    multiBamSummary bins --bamfiles $(paste -sd ' ' "$BAMLIST") -o "$BAM_CORR_NPZ" --binSize 5000
    plotCorrelation -in "$BAM_CORR_NPZ" -c pearson --corMethod pearson \
        --plotFileName "$CORR_HEAT" --plotNumbers
    echo "✅ Correlation heatmap saved: $CORR_HEAT"
else
    echo "⚠️ No BAMs found for correlation step"
fi

echo "🎉 All steps complete! Master summary is: $MASTER_SUMMARY"
