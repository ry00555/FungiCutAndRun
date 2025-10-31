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
module load GCC/12.3.0
ml BEDTools/2.30.0-GCC-11.3.0 deepTools SAMtools/0.1.20-GCC-11.3.0 BamTools/2.5.2-GCC-11.3.0

META="/scratch/ry00555/RNASeqPaper/Oct2025/BAM_File_Metadata_with_index_merged_V2.csv"
MACSDIR="/scratch/ry00555/RNASeqPaper/Oct2025/MACSPeaks"
CHIPR_DIR="/scratch/ry00555/RNASeqPaper/Oct2025/ChIPR"
BAMDIR="/scratch/ry00555/RNASeqPaper/Oct2025/SortedBamFiles"
OUTDIR="/scratch/ry00555/RNASeqPaper/Oct2025/IDR"
SUMMARY="${OUTDIR}/replicate_vs_CHIPRconsensus.tsv"
COMBINED="${OUTDIR}/all_combined.tsv"


FRIP_TSV="${OUTDIR}/frip_metrics.tsv"
OVERLAP_TSV="${OUTDIR}/replicate_overlap.tsv"
JACCARD_TSV="${OUTDIR}/pairwise_jaccard.tsv"
BAM_CORR_NPZ="${OUTDIR}/bam_corr.npz"
CORR_HEAT="${OUTDIR}/bam_correlation_heatmap.pdf"


echo "---- Calculating FRiP and overlap ----"
echo -e "SampleID\tTissue\tFactor\tTotalReads\tReadsInPeaks\tFRiP" > "$FRIP_TSV"
echo -e "Tissue\tSampleID\tPeakFile\tNumPeaks\tNumOverlap\tFracOverlap" > "$OVERLAP_TSV"

tail -n +2 "$META" | while IFS=, read -r RunID bamReads BamIndex SampleID Factor Tissue Condition Replicate bamControl bamInputIndex ControlID Peaks PeakCaller DesiredPeakName; do
  [[ -z "$SampleID" || -z "$Tissue" ]] && continue


  bam="${BAMDIR}/${bamReads}"
  for bam in /scratch/ry00555/RNASeqPaper/Oct2025/SortedBamFiles/*.bam; do
      samtools index -@ 4 "$bam"
  done

  peak="${MACSDIR}/${SampleID}_${Factor}_${Replicate}_peaks.broadPeak"
  consensus="${CHIPR_DIR}/${Tissue}_consensus_optimal_peaks.broadPeak"

  [[ ! -s "$bam" || ! -s "$peak" ]] && continue

  total_reads=$(samtools view -c -F 260 "$bam" 2>/dev/null || echo 0)
  reads_in_peaks=$(bedtools coverage -a "$peak" -b "$bam" -counts | awk '{sum += $7} END{print sum+0}')
  num_peaks=$(wc -l < "$peak" | tr -d '[:space:]')

  if [[ -s "$consensus" && "$num_peaks" -gt 0 ]]; then
    num_overlap=$(bedtools intersect -u -a "$peak" -b "$consensus" | wc -l | tr -d '[:space:]')
    frac_overlap=$(awk "BEGIN{printf \"%.4f\", ($num_overlap/$num_peaks)}")
  else
    num_overlap=0; frac_overlap=0
  fi

  frip=$(awk "BEGIN{printf \"%.4f\", ($reads_in_peaks/$total_reads)}")
  echo -e "${SampleID}\t${Tissue}\t${Factor}\t${total_reads}\t${reads_in_peaks}\t${frip}" >> "$FRIP_TSV"
  echo -e "${Tissue}\t${SampleID}\t${peak}\t${num_peaks}\t${num_overlap}\t${frac_overlap}" >> "$OVERLAP_TSV"
done

# ================================
# 5. JACCARD SIMILARITY
# ================================
echo "---- Calculating Jaccard similarity ----"
echo -e "fileA\tfileB\tjaccard" > "$JACCARD_TSV"

PEAK_FILES=( $(awk -F'\t' 'NR>1{print $3}' "$OVERLAP_TSV" | sort -u) )
for ((i=0;i<${#PEAK_FILES[@]};i++)); do
  for ((j=i+1;j<${#PEAK_FILES[@]};j++)); do
    f1="${PEAK_FILES[i]}"
    f2="${PEAK_FILES[j]}"
    if [[ -s "$f1" && -s "$f2" ]]; then
      jacc=$(bedtools jaccard -a "$f1" -b "$f2" | awk 'NR==2{print $3}')
      echo -e "$(basename "$f1")\t$(basename "$f2")\t${jacc:-0}" >> "$JACCARD_TSV"
    fi
  done
done

# ================================
# 6. MULTIBAMSUMMARY + PLOT CORRELATION
# ================================
BAMLIST="${OUTDIR}/bamlist.txt"
awk -F, 'NR>1{print $2}' "$META" | while read -r b; do
  bampath="${BAMDIR}/${b}"
  [[ -s "$bampath" ]] && echo "$bampath"
done > "$BAMLIST"

if [[ -s "$BAMLIST" ]]; then
  echo "---- Running deeptools correlation ----"
  multiBamSummary bins --bamfiles $(paste -sd ' ' "$BAMLIST") \
    -o "$BAM_CORR_NPZ" --binSize 5000
  plotCorrelation -in "$BAM_CORR_NPZ" -c pearson --corMethod pearson \
    --plotFileName "$CORR_HEAT" --plotNumbers
  echo "✅ Correlation heatmap saved: $CORR_HEAT"
else
  echo "⚠️ No BAMs found for correlation step"
fi

echo "🎉 All steps complete!"
