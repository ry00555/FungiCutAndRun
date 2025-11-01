#!/bin/bash
#SBATCH --job-name=BamQC
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=200gb
#SBATCH --time=6:00:00
#SBATCH --output=../BamQC.%j.out
#SBATCH --error=../BamQC.%j.err

set -euo pipefail

 #================================
 #Load modules
 #================================
ml BEDTools/2.30.0-GCC-11.3.0 deepTools SAMtools/1.16.1-GCC-11.3.0 BamTools/2.5.2-GCC-11.3.0

 #================================
 #Paths
 #================================
META="/scratch/ry00555/RNASeqPaper/Oct2025/BAM_File_Metadata_with_index_merged_V2.csv"
MACSDIR="/scratch/ry00555/RNASeqPaper/Oct2025/MACSPeaks"
CHIPR_DIR="/scratch/ry00555/RNASeqPaper/Oct2025/ChIPR"
BAMDIR="/scratch/ry00555/RNASeqPaper/Oct2025/SortedBamFiles"
OUTDIR="/scratch/ry00555/RNASeqPaper/Oct2025/IDR"
dos2unix "$META" 2>/dev/null || true

MASTER_SUMMARY="${OUTDIR}/master_summary.tsv"
FRIP_TSV="${OUTDIR}/frip_metrics.tsv"
OVERLAP_TSV="${OUTDIR}/replicate_overlap.tsv"
JACCARD_TSV="${OUTDIR}/pairwise_jaccard.tsv"
BAM_CORR_NPZ="${OUTDIR}/bam_corr.npz"
CORR_HEAT="${OUTDIR}/bam_correlation_heatmap.pdf"

mkdir -p "$OUTDIR"

 #================================
 #Initialize master summary
 #================================
echo -e "SampleID\tTissue\tFactor\tTotalReads\tReadsInPeaks\tFRiP\tPeakFile\tConsensusPeak\tNumPeaks\tNumOverlap\tFracOverlap" > "$MASTER_SUMMARY"
echo -e "SampleID\tTissue\tFactor\tTotalReads\tReadsInPeaks\tFRiP" > "$FRIP_TSV"
echo -e "Tissue\tSampleID\tPeakFile\tNumPeaks\tNumOverlap\tFracOverlap" > "$OVERLAP_TSV"

 #================================
 #Index BAMs only if needed
 #================================
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

# If MASTER_SUMMARY doesn't exist or is empty, create processed_ids.txt
if [[ ! -s "$MASTER_SUMMARY" ]]; then
    echo "" > /tmp/processed_ids.txt
else
    awk '{print $1}' "$MASTER_SUMMARY" > /tmp/processed_ids.txt
fi

# Main metadata loop
tail -n +2 "$META" | while IFS=, read -r RunID bamReads BamIndex SampleID Factor Tissue Condition Replicate bamControl bamInputIndex ControlID Peaks PeakCaller DesiredPeakName MACS3minlength MACS3maxgap; do

    [[ -z "$SampleID" || -z "$Tissue" || -z "$bamReads" ]] && continue

    # Skip if processed already
    if grep -qx "$SampleID" /tmp/processed_ids.txt; then
        echo "Skipping (already processed): $SampleID"
        continue
    fi

    echo "Processing: $SampleID"

    bam="${BAMDIR}/${bamReads}"
    peak="${MACSDIR}/${DesiredPeakName}_peaks.broadPeak"
    consensus="${CHIPR_DIR}/${Tissue}_consensus_optimal_peaks.broadPeak"

    [[ ! -s "$bam" ]] && { echo "Missing BAM → skipping $SampleID"; continue; }
    [[ ! -s "$peak" ]] && { echo "Missing peak file → skipping $SampleID"; continue; }

    # FRiP calculation
    total_reads=$(samtools view -c -F 260 "$bam" 2>/dev/null || echo 0)
    reads_in_peaks=$(bedtools intersect -a "$peak" -b "$bam" -c | awk '{sum+=$NF} END{print sum+0}')
    frip=$(awk "BEGIN{printf \"%.4f\", ($reads_in_peaks/$total_reads)}")

    # Peak overlap with consensus
    num_peaks=$(wc -l < "$peak" | tr -d '[:space:]')
    if [[ -s "$consensus" && "$num_peaks" -gt 0 ]]; then
        num_overlap=$(bedtools intersect -u -a "$peak" -b "$consensus" | wc -l | tr -d '[:space:]')
        frac_overlap=$(awk "BEGIN{printf \"%.4f\", ($num_overlap/$num_peaks)}")
    else
        num_overlap=0
        frac_overlap=0
    fi

    # Append outputs
    echo -e "${SampleID}\t${Tissue}\t${Factor}\t${total_reads}\t${reads_in_peaks}\t${frip}\t${peak}\t${consensus}\t${num_peaks}\t${num_overlap}\t${frac_overlap}" >> "$MASTER_SUMMARY"
    echo -e "${SampleID}\t${Tissue}\t${Factor}\t${total_reads}\t${reads_in_peaks}\t${frip}" >> "$FRIP_TSV"
    echo -e "${Tissue}\t${SampleID}\t${peak}\t${num_peaks}\t${num_overlap}\t${frac_overlap}" >> "$OVERLAP_TSV"

done


# ================================
# JACCARD SIMILARITY
# ================================
echo "---- Calculating Jaccard similarity ----"
echo -e "fileA\tfileB\tjaccard" > "$JACCARD_TSV"

# Collect unique peak files
mapfile -t PEAK_FILES < <(awk -F'\t' 'NR>1{print $7}' "$MASTER_SUMMARY" | sort -u)

for ((i=0; i<${#PEAK_FILES[@]}; i++)); do
    for ((j=i+1; j<${#PEAK_FILES[@]}; j++)); do

        f1="${PEAK_FILES[i]}"
        f2="${PEAK_FILES[j]}"

        if [[ -s "$f1" && -s "$f2" ]]; then
            jacc=$(bedtools jaccard -a "$f1" -b "$f2" | awk 'NR==2{print $3}')
            echo -e "$(basename "$f1")\t$(basename "$f2")\t${jacc:-0}" >> "$JACCARD_TSV"
        fi
    done
done


# ===============================
# multiBamSummary per Tissue
# ===============================
SKIPPED_BAMS="/tmp/skipped_bams.txt"
> "$SKIPPED_BAMS"

declare -A groups   # associative array

echo "Grouping BAMs by Tissue from META file..."

# Build tissue → bam list
tail -n +2 "$META" | while IFS=, read -r RunID bamReads BamIndex SampleID Factor Tissue Condition Replicate bamControl bamInputIndex ControlID Peaks PeakCaller DesiredPeakName MACS3minlength MACS3maxgap; do

    bam="$BAMDIR/${SampleID}.bam"

    if [[ ! -f "$bam" ]]; then
        echo "⚠️ Missing BAM: $bam" | tee -a "$SKIPPED_BAMS"
        continue
    fi

    if ! samtools quickcheck "$bam" 2>/dev/null; then
        echo "⚠️ Unreadable BAM: $bam" | tee -a "$SKIPPED_BAMS"
        continue
    fi

    groups["$Tissue"]+="$bam "

done

echo "✅ Finished grouping BAMs by Tissue."

# Run multiBamSummary per group
for tissue in "${!groups[@]}"; do
    echo "-----------------------------"
    echo "Processing Tissue: $tissue"
    echo "-----------------------------"

    bams="${groups[$tissue]}"
    bam_count=$(echo $bams | wc -w)

    echo "Found $bam_count BAMs"

    if [[ $bam_count -lt 2 ]]; then
        echo "⚠️ Skipping $tissue (needs ≥2 BAMs)"
        continue
    fi

    tissue_npz="${OUTDIR}/${tissue}_multiBamSummary.npz"
    tissue_tab="${OUTDIR}/${tissue}_readCounts.tab"
    tissue_heat="${OUTDIR}/${tissue}_correlation_heatmap.png"
    tissue_matrix="${OUTDIR}/${tissue}_PearsonCorrMatrix.tab"

    multiBamSummary bins \
        --bamfiles $bams \
        --binSize 50000 \
        -p max \
        --smartLabels \
        -o "$tissue_npz" \
        --outRawCounts "$tissue_tab"

    plotCorrelation \
        -in "$tissue_npz" \
        -c pearson \
        --whatToPlot heatmap \
        --plotNumbers \
        --skipZeros \
        --plotFile "$tissue_heat" \
        --outFileCorMatrix "$tissue_matrix"

    echo "✅ Finished Tissue group: $tissue"

done

if [[ -s "$SKIPPED_BAMS" ]]; then
    echo "⚠️ Skipped BAMs:"
    cat "$SKIPPED_BAMS"
fi
