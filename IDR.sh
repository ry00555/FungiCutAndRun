#!/bin/bash
#SBATCH --job-name=IDR
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=4:00:00
#SBATCH --output=../IDR.%j.out
#SBATCH --error=../IDR.%j.err

set -euo pipefail

ml BEDTools/2.30.0-GCC-11.3.0 deepTools SAMtools/1.16.1-GCC-11.3.0 BamTools/2.5.2-GCC-11.3.0

META="/scratch/ry00555/RNASeqPaper/Oct2025/BAM_File_Metadata_with_index_merged_V2.csv"
MACSDIR="/scratch/ry00555/RNASeqPaper/Oct2025/MACSPeaks"
CHIPR_DIR="/scratch/ry00555/RNASeqPaper/Oct2025/ChIPR"
BAMDIR="/scratch/ry00555/RNASeqPaper/Oct2025/SortedBamFiles"
OUTDIR="/scratch/ry00555/RNASeqPaper/Oct2025/IDR"

MASTER_SUMMARY="${OUTDIR}/master_summary.tsv"
JACCARD_TSV="${OUTDIR}/pairwise_jaccard.tsv"
BAM_CORR_NPZ="${OUTDIR}/bam_corr.npz"
CORR_HEAT="${OUTDIR}/bam_correlation_heatmap.pdf"

mkdir -p "$OUTDIR"

echo -e "SampleID\tTissue\tFactor\tReplicate\tTotalReads\tReadsInPeaks\tFRiP\tPeakFile\tNumPeaks\tNumOverlap\tFracOverlap" > "$MASTER_SUMMARY"

echo "---- Checking and reindexing BAMs if needed ----"
for b in "$BAMDIR"/*.bam; do
    bai="${b}.bai"
    if [[ ! -s "$bai" || "$b" -nt "$bai" ]]; then
        echo "Indexing BAM: $b"
        samtools index -@ 4 "$b"
    else
        echo "BAM index up-to-date: $bai"
    fi
done

echo "---- Calculating FRiP and overlap ----"
tail -n +2 "$META" | while IFS=, read -r RunID bamReads BamIndex SampleID Factor Tissue Condition Replicate bamControl bamInputIndex ControlID Peaks PeakCaller DesiredPeakName; do
    [[ -z "$SampleID" || -z "$Tissue" || -z "$bamReads" ]] && continue

    bam="${BAMDIR}/${bamReads}"
    peak="${MACSDIR}/${SampleID}_${Factor}_${Replicate}_peaks.broadPeak"
    consensus="${CHIPR_DIR}/${Tissue}_consensus_optimal_peaks.broadPeak"

    [[ ! -s "$bam" || ! -s "$peak" ]] && continue

    echo "Processing sample: $SampleID ($Tissue, $Factor, replicate $Replicate)"

    total_reads=$(samtools view -c -F 260 "$bam" 2>/dev/null || echo 0)
    reads_in_peaks=$(bedtools intersect -a "$peak" -b "$bam" -c | awk '{sum+=$NF} END{print sum+0}')
    num_peaks=$(wc -l < "$peak" | tr -d '[:space:]')

    if [[ -s "$consensus" && "$num_peaks" -gt 0 ]]; then
        num_overlap=$(bedtools intersect -u -a "$peak" -b "$consensus" | wc -l | tr -d '[:space:]')
        frac_overlap=$(awk "BEGIN{printf \"%.4f\", ($num_overlap/$num_peaks)}")
    else
        num_overlap=0; frac_overlap=0
    fi

    frip=$(awk "BEGIN{printf \"%.4f\", ($reads_in_peaks/$total_reads)}")

    # Write to master summary
    echo -e "${SampleID}\t${Tissue}\t${Factor}\t${Replicate}\t${total_reads}\t${reads_in_peaks}\t${frip}\t${peak}\t${num_peaks}\t${num_overlap}\t${frac_overlap}" >> "$MASTER_SUMMARY"
done

# ================================
# JACCARD SIMILARITY
# ================================
echo "---- Calculating Jaccard similarity ----"
echo -e "fileA\tfileB\tjaccard" > "$JACCARD_TSV"

PEAK_FILES=( $(awk -F'\t' 'NR>1{print $8}' "$MASTER_SUMMARY" | sort -u) )
for ((i=0;i<${#PEAK_FILES[@]};i++)); do
    for ((j=i+1;j<${#PEAK_FILES[@]};j++)); do
        f1="${PEAK_FILES[i]}"
        f2="${PEAK_FILES[j]}"
        if [[ -s "$f1" && -s "$f2" ]]; then
            echo "Comparing peaks: $(basename "$f1") vs $(basename "$f2")"
            jacc=$(bedtools jaccard -a "$f1" -b "$f2" | awk 'NR==2{print $3}')
            echo -e "$(basename "$f1")\t$(basename "$f2")\t${jacc:-0}" >> "$JACCARD_TSV"
        fi
    done
done

# ================================
# MULTIBAMSUMMARY + PLOT CORRELATION
# ================================
BAMLIST="${OUTDIR}/bamlist.txt"
awk -F'\t' 'NR>1{print $5}' "$MASTER_SUMMARY" | while read -r b; do
    [[ -s "$b" ]] && echo "$b"
done > "$BAMLIST"

if [[ -s "$BAMLIST" ]]; then
    echo "---- Running multiBamSummary ----"
    while read -r bamfile; do
        echo "Including BAM: $bamfile"
    done < "$BAMLIST"

    multiBamSummary bins --bamfiles $(paste -sd ' ' "$BAMLIST") -o "$BAM_CORR_NPZ" --binSize 5000
    echo "---- Plotting correlation heatmap ----"
    plotCorrelation -in "$BAM_CORR_NPZ" -c pearson --corMethod pearson \
        --plotFileName "$CORR_HEAT" --plotNumbers
    echo "✅ Correlation heatmap saved: $CORR_HEAT"
else
    echo "⚠️ No BAMs found for correlation step"
fi

echo "🎉 All steps complete! Master summary is: $MASTER_SUMMARY"
