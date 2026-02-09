#!/bin/bash
#SBATCH --job-name=multiBamSummary
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100gb
#SBATCH --time=5:00:00
#SBATCH --output=../multiBamSummary.%j.out
#SBATCH --error=../multiBamSummary.%j.err

set -euo pipefail

 #================================
 #Load modules
 #================================
ml BEDTools/2.30.0-GCC-11.3.0 deepTools SAMtools/1.16.1-GCC-11.3.0 BamTools/2.5.2-GCC-11.3.0

 #================================
 #Paths
 #================================
META="/lustre2/scratch/ry00555/RTT109PaperFigures/RTT109Paper_AllMeta_CorrPlot.csv"
BAMDIR="/lustre2/scratch/ry00555/RTT109PaperFigures/SortedBamFiles"
OUTDIR="/lustre2/scratch/ry00555/RTT109PaperFigures/IDR"
dos2unix "$META" 2>/dev/null || true

BAM_CORR_NPZ="${OUTDIR}/bam_corr.npz"
CORR_HEAT="${OUTDIR}/bam_correlation_heatmap.pdf"


# ===============================
# multiBamSummary per Tissue
# ===============================
SKIPPED_BAMS="/tmp/skipped_bams.txt"
> "$SKIPPED_BAMS"

# Build tissue → bam list
declare -A groups
echo "Grouping BAMs by Tissue from META file..."

while IFS=, read -r RunID bamReads BamIndex SampleID Factor Tissue Condition Replicate bamControl bamInputIndex ControlID Peaks PeakCaller DesiredPeakName MACS3minlength MACS3maxgap; do

    bam="$BAMDIR/${bamReads}"

    if [[ ! -f "$bam" ]]; then
        echo "⚠️ Missing BAM: $bam" | tee -a "$SKIPPED_BAMS"
        continue
    fi

    if ! samtools quickcheck "$bam" 2>/dev/null; then
        echo "⚠️ Unreadable BAM: $bam" | tee -a "$SKIPPED_BAMS"
        continue
    fi

    groups["$Tissue"]+="$bam "

done < <(tail -n +2 "$META")


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
        --numberOfProcessors max \
        -p max \
        --smartLabels \
        -o "$tissue_npz" \
        --outRawCounts "$tissue_tab"

    plotCorrelation \
        -in "$tissue_npz" \
        -c spearman \
        --whatToPlot heatmap \
        --plotNumbers \
        --skipZeros \
        --plotHeight 3.5 \
        --plotWidth 3.5 \
        --plotFile "$tissue_heat" \
        --outFileCorMatrix "$tissue_matrix"

    echo "✅ Finished Tissue group: $tissue"

done

if [[ -s "$SKIPPED_BAMS" ]]; then
    echo "⚠️ Skipped BAMs:"
    cat "$SKIPPED_BAMS"
fi

# source code:https://deeptools.readthedocs.io/en/develop/content/tools/plotCorrelation.html
# Performing correlation analysis between bam files (selecting the pertinent bigwig files)
# multiBigwigSummary bins \
# -b IP-1_1.bigwig IP-1_2.bigwig \
# --labels IP-1_1 IP-1_2 \
# -out scores_per_bin.npz --outRawCounts scores_per_bin.tab
# plotCorrelation \
# -in scores_per_bin.npz \
# --corMethod pearson --skipZeros \
# --plotTitle "Pearson Correlation of Average Scores Per Transcript" \
# --whatToPlot scatterplot \
# -o scatterplot_PearsonCorr_bigwigScores.png \
# --outFileCorMatrix PearsonCorr_bigwigScores.tab
