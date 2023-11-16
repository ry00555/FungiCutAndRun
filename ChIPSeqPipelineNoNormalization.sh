#!/bin/bash
#SBATCH --job-name=ChIPSeqPipeline
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=../ChIPSeqPipelinePooled.%j.out
#SBATCH --error=../ChIPSeqPipelinePooled.%j.err

OUTDIR=/scratch/ry00555/PoolChIPSeq

#Load all the modules that are needed for the entire pipeline
#ml BWA/0.7.17-GCCcore-11.3.0
#ml SAMtools/1.16.1-GCC-11.3.0
#ml Homer/4.11-foss-2022a
ml deepTools/3.5.2-foss-2022a
# ml Perl/5.34.1-GCCcore-11.3.0
#ml Trim_Galore/0.6.5-GCCcore-8.3.0-Java-11-Python-3.7.4
ml BEDTools/2.30.0-GCC-11.3.0


#mkdir -p "${OUTDIR}/SortedBamFiles"
#mkdir -p "${OUTDIR}/BigWigs"
#mkdir -p "${OUTDIR}/Peaks"
#mkdir -p "$OUTDIR/HomerTagDirectories"
#mkdir -p "$OUTDIR/TdfFiles"
#mkdir -p "${OUTDIR}/Heatmaps"
#mkdir -p "${OUTDIR}/Matrices"
#mkdir -p "${OUTDIR}/NormalizedBigWigs"
#mkdir -p "${OUTDIR}/Beds"
#mkdir -p "${OUTDIR}/Counts"


#PEAKDIR="${OUTDIR}/Peaks"
#TAGDIR="${OUTDIR}/HomerTagDirectories"
#BAMDIR="${OUTDIR}/SortedBamFiles"
#BEDDIR="${OUTDIR}/Beds"

#for bam_file in "${BAMDIR}"/*.bam; do
  # Get the sample ID from the BAM file name
#  sample_id=$(basename "${bam_file}" .bam)
  # Remove everything after "Rep_1" in the sample ID
#  sample_id="${sample_id%%_Rep_1*}"

# Make tag directory
#makeTagDirectory "${TAGDIR}/${sample_id}" "${bam_file}"
# Check if the sample_id contains "input" and set the -i argument accordingly
# if [[ "${sample_id}" == *input* ]]; then
   #findPeaks "${TAGDIR}/${sample_id}" -style histone -region -size 150 -minDist 530 -o "${PEAKDIR}/${sample_id}_peaks.txt" -i "${TAGDIR}/${sample_id}"
# else
#   findPeaks "${TAGDIR}/${sample_id}" -style histone -region -size 150 -minDist 530 -o "${PEAKDIR}/${sample_id}_peaks.txt"
# fi

 # for file_path in "${OUTDIR}/BigWigs"/*.bw; do
 #   # Get the base name of the file
 #   BW_name=$(basename "${file_path}" .bw)
 #
 #   # Remove the file extension to create the sample ID
 #   BW_id="${BW_name%.*}"
 #   # Replace special characters with underscores in the sample ID
 #  BW_id=${BW_id//[^a-zA-Z0-9]/_}
 #
 #  # Limit the length of the sample ID to avoid long filenames
 #  BW_id=${BW_id:0:50}
  # Compute matrix for the reference-point TSS
#done
#  computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S "${OUTDIR}/BigWigs/133-90_ChIP_WT_H3K27me3_Rep1_S87_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw" "${OUTDIR}/BigWigs/133-78_ChIP_NCU00423_H3K27me3_Rep1_S75_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw" -R "/scratch/ry00555/neurospora.bed" --skipZeros -o "${OUTDIR}/Matrices/NCU00423_H3K27me3.gz"
#  plotHeatmap --matrixFile "${OUTDIR}/Matrices/NCU00423_H3K27me3.gz" --outFileName "${OUTDIR}/Heatmaps/NCU00423_H3K27me3_hclust.png" --samplesLabel WT NCU00423-KO --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1

#computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S ${OUTDIR}/BigWigs/133-91_ChIP_WT_H3K36me3_Rep1_S88_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/133-79_ChIP_NCU00423_H3K36me3_Rep1_S76_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw -R "/scratch/ry00555/neurospora.bed" --skipZeros -o "${OUTDIR}/Matrices/NCU00423_H3K36me3.gz"
#  plotHeatmap --matrixFile "${OUTDIR}/Matrices/NCU00423_H3K36me3.gz" --outFileName "${OUTDIR}/Heatmaps/NCU00423_H3K36me3_hclust.png" --samplesLabel WT NCU00423-KO --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1


 computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S "${OUTDIR}/BigWigs/133-90_ChIP_WT_H3K27me3_Rep1_S87_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw" "${OUTDIR}/BigWigs/133-78_ChIP_NCU00423_H3K27me3_Rep1_S75_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw" -R "/scratch/ry00555/heatmapPRC2genes.bed" --skipZeros -o "${OUTDIR}/Matrices/NCU00423_H3K27me3_133only_1stsetPRC2target.gz"
  plotHeatmap --matrixFile "${OUTDIR}/Matrices/NCU00423_H3K27me3_133only_1stsetPRC2target.gz" --outFileName "${OUTDIR}/Heatmaps/NCU00423_H3K27me3_hclust2.png" --samplesLabel WT NCU00423-KO --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1

computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S ${OUTDIR}/BigWigs/133-91_ChIP_WT_H3K36me3_Rep1_S88_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/133-79_ChIP_NCU00423_H3K36me3_Rep1_S76_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw -R "/scratch/ry00555/heatmapPRC2genes.bed" --skipZeros -o "${OUTDIR}/Matrices/NCU00423_H3K36me3_133only_1stsetPRC2target.gz"
  plotHeatmap --matrixFile "${OUTDIR}/Matrices/NCU00423_H3K36me3_133only_1stsetPRC2target.gz" --outFileName "${OUTDIR}/Heatmaps/NCU00423_H3K36me3_hclust2.png" --samplesLabel WT NCU00423-KO --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1



#Pool all chip seq together 11/12/23
#NCU00423 H3K27me3
computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S "${OUTDIR}/BigWigs/133-90_ChIP_WT_H3K27me3_Rep1_S87_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw" ${OUTDIR}/BigWigs/128-53_ChIP_-NCU00423__hph_H3K27me3_Rep1_S32_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/132-30_ChIP_ncu00423_H3K27me3_Rep_1_S27_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/132-34_ChIP_ncu00423_Input_Rep_1_S31_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/132-35_ChIP_ncu00423_H3K27me3_Rep_1_S32_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/133-78_ChIP_NCU00423_H3K27me3_Rep1_S75_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw -R "/scratch/ry00555/neurospora.bed" --skipZeros -o "${OUTDIR}/Matrices/NCU00423_H3K27me3_Pool2.gz"
plotHeatmap --matrixFile "${OUTDIR}/Matrices/NCU00423_H3K27me3_Pool2.gz" --outFileName "${OUTDIR}/Heatmaps/NCU00423_H3K27me3_Pool2.png" --samplesLabel WT 128-53 131-52 132-30 132-34 133-78 --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1


computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S "${OUTDIR}/BigWigs/133-90_ChIP_WT_H3K27me3_Rep1_S87_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw" ${OUTDIR}/BigWigs/128-53_ChIP_-NCU00423__hph_H3K27me3_Rep1_S32_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/132-30_ChIP_ncu00423_H3K27me3_Rep_1_S27_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/132-34_ChIP_ncu00423_Input_Rep_1_S31_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/132-35_ChIP_ncu00423_H3K27me3_Rep_1_S32_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw  -R "/scratch/ry00555/heatmapPRC2genes.bed" --skipZeros -o "${OUTDIR}/Matrices/NCU00423_H3K27me3_Pool_1stsetPRC2target3.gz"
plotHeatmap --matrixFile"${OUTDIR}/Matrices/NCU00423_H3K27me3_Pool_1stsetPRC2target3.gz" --outFileName "${OUTDIR}/Heatmaps/NCU00423_H3K27me3_Pool_1stsetPRC2target3.png" --samplesLabel WT 128-53 131-52 132-30 132-34 133-78 --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1



#NCU00423 NCU00423_H3K36me3
#computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S ${OUTDIR}/BigWigs/133-91_ChIP_WT_H3K36me3_Rep1_S88_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw  ${OUTDIR}/BigWigs/132-32_ChIP_ncu00423_H3K36me3_Rep_1_S29_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/132-37_ChIP_ncu00423_H3K36me3_Rep_1_S34_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/133-79_ChIP_NCU00423_H3K36me3_Rep1_S76_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/128-51_ChIP_-NCU00423__hph_H3K36me3_Rep1_S30_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw  -R "/scratch/ry00555/neurospora.bed" --skipZeros -o "${OUTDIR}/Matrices/NCU00423_H3K36me3_Pool2.gz"
#plotHeatmap --matrixFile "${OUTDIR}/Matrices/NCU00423_H3K36me3_Pool2.gz" --outFileName "${OUTDIR}/Heatmaps/NCU00423_H3K36me3_Pool2.png" --samplesLabel WT 131-54 132-37 133-79 128-51 --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1

#computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S ${OUTDIR}/BigWigs/133-91_ChIP_WT_H3K36me3_Rep1_S88_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/131-54_ChIP_WT_H3K36me3_Rep1_S42_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/132-32_ChIP_ncu00423_H3K36me3_Rep_1_S29_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/132-37_ChIP_ncu00423_H3K36me3_Rep_1_S34_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/133-79_ChIP_NCU00423_H3K36me3_Rep1_S76_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/128-51_ChIP_-NCU00423__hph_H3K36me3_Rep1_S30_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw  -R "/scratch/ry00555/heatmapPRC2genes.bed" --skipZeros -o "${OUTDIR}/Matrices/NCU00423_H3K36me3.gz"
#plotHeatmap --matrixFile "${OUTDIR}/Matrices/NCU00423_H3K36me3_1stsetPRC2target.gz" --outFileName "${OUTDIR}/Heatmaps/NCU00423_H3K36me3_1stsetPRC2target2.png" --samplesLabel WT 131-54 132-37 133-79 128-51 --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1


#NCU00423 NCU00423_H3K27ac
#computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S ${OUTDIR}/BigWigs/131-53_ChIP_WT_H3K27ac_Rep1_S41_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/132-31_ChIP_ncu00423_H3K27ac_Rep_1_S28_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/132-36_ChIP_ncu00423_H3K27ac_Rep_1_S33_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw -R "/scratch/ry00555/heatmapPRC2genes.bed" --skipZeros -o "${OUTDIR}/Matrices/NCU00423_H3K27ac_Pool_1stsetPRC2target.gz"
#  plotHeatmap --matrixFile "${OUTDIR}/Matrices/NCU00423_H3K27ac_Pool_1stsetPRC2target.gz" --outFileName "${OUTDIR}/Heatmaps/NCU00423_H3K27ac_Pool_1stsetPRC2target2.png" --samplesLabel WT 132-31 132-36 --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1



#rtt109 by Run
#H3K27me3
computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S ${OUTDIR}/BigWigs/133-90_ChIP_WT_H3K27me3_Rep1_S87_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/131-64_ChIP_ncu00548_H3K27me3_Rep1_S52_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw -R "/scratch/ry00555/neurospora.bed" --skipZeros -o "${OUTDIR}/Matrices/rtt109_H3K27me3_Run131Only2.gz"
  plotHeatmap --matrixFile "${OUTDIR}/Matrices/rtt109_H3K27me3_Run131Only2.gz" --outFileName "${OUTDIR}/Heatmaps/rtt109_H3K27me3_Run131Only2.png" --samplesLabel WT 131-64 --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1


computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S ${OUTDIR}/BigWigs/133-90_ChIP_WT_H3K27me3_Rep1_S87_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/133-81_ChIP_rtt109_hph_H3K27me3_Rep1_S78_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw -R "/scratch/ry00555/neurospora.bed" --skipZeros -o "${OUTDIR}/Matrices/rtt109_H3K27me3_Run133Only2.gz"
plotHeatmap --matrixFile "${OUTDIR}/Matrices/rtt109_H3K27me3_Run133Only2.gz" --outFileName "${OUTDIR}/Heatmaps/rtt109_H3K27me3_Run133Only2.png" --samplesLabel WT 133-90 133-81 --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1

#H3K27ac
#computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S ${OUTDIR}/BigWigs/131-53_ChIP_WT_H3K27ac_Rep1_S41_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/131-65_ChIP_ncu00548_H3K27ac_Rep1_S53_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw -R "/scratch/ry00555/neurospora.bed" --skipZeros -o "${OUTDIR}/Matrices/rtt109_H3K27ac_Run131Only.gz"
#plotHeatmap --matrixFile "${OUTDIR}/Matrices/rtt109_H3K27ac_Run131Only.gz" --outFileName "${OUTDIR}/Heatmaps/rtt109_H3K27ac_Run131Only_hclust2.png" --samplesLabel WT 131-65 --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1

#rtt109_H3K36me3
#computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S ${OUTDIR}/BigWigs/133-91_ChIP_WT_H3K36me3_Rep1_S88_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/131-66_ChIP_ncu00548_H3K36me3_Rep1_S54_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw -R "/scratch/ry00555/neurospora.bed" --skipZeros -o "${OUTDIR}/Matrices/rtt109_H3K36me3_Run131Only.gz"
#plotHeatmap --matrixFile "${OUTDIR}/Matrices/rtt109_H3K36me3_Run131Only.gz" --outFileName "${OUTDIR}/Heatmaps/ncu00548_H3K36me3_Run131Only_hclust.png" --samplesLabel WT 131-66 --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1

computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S ${OUTDIR}/BigWigs/133-91_ChIP_WT_H3K36me3_Rep1_S88_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/133-82_ChIP_rtt109_hph_H3K36me3_Rep1_S79_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw -R "/scratch/ry00555/neurospora.bed" --skipZeros -o "${OUTDIR}/Matrices/rtt109_H3K36me3_Run133Only3.gz"
plotHeatmap --matrixFile "${OUTDIR}/Matrices/rtt109_H3K36me3_Run133Only3.gz" --outFileName "${OUTDIR}/Heatmaps/rtt109_H3K36me3_Run133Only3.png" --samplesLabel WT 133-82 --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1

computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S ${OUTDIR}/BigWigs/133-91_ChIP_WT_H3K36me3_Rep1_S88_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/133-82_ChIP_rtt109_hph_H3K36me3_Rep1_S79_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw -R "/scratch/ry00555/heatmapPRC2genes.bed" --skipZeros -o "${OUTDIR}/Matrices/rtt109_H3K36me3_Run133Only_PRC2Genes3.gz"
plotHeatmap --matrixFile "${OUTDIR}/Matrices/rtt109_H3K36me3_Run133Only_PRC2Genes3.gz" --outFileName "${OUTDIR}/Heatmaps/rtt109_H3K36me3_Run133Only_PRC2Genes3.png" --samplesLabel WT 133-82 --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1


#computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S ${OUTDIR}/BigWigs/132-17_ChIP_WT_H4K16ac_Rep_1_S14_R1_001_val_1.fq.gz.bin_25.smooth_75_MNase.bw ${OUTDIR}/BigWigs/132-40_ChIP_rtt109_H4K16ac_Rep_1_S37_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw "/scratch/ry00555/neurospora.bed" --skipZeros -o "${OUTDIR}/Matrices/rtt109_H4K16ac.gz"
#plotHeatmap --matrixFile "${OUTDIR}/Matrices/rtt109_H4K16ac.gz" --outFileName "${OUTDIR}/Heatmaps/rtt109_H4K16ac_Run132Only_hclust2.png" --samplesLabel WT rtt109 --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1

#ncu06788

computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S ${OUTDIR}/BigWigs/133-90_ChIP_WT_H3K27me3_Rep1_S87_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw  ${OUTDIR}/BigWigs/134-14_ChIP_ncu06788_H3K27me3_Rep2.bin_25.smooth_75Bulk.bw -R "/scratch/ry00555/neurospora.bed" --skipZeros -o "${OUTDIR}/Matrices/ncu06788_H3K27me3_13414_2.gz"
plotHeatmap --matrixFile "${OUTDIR}/Matrices/ncu06788_H3K27me3_13414_2.gz" --outFileName "${OUTDIR}/Heatmaps/ncu06788_H3K27me3_13414_2.png" --samplesLabel WT 134-14  --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1

computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S ${OUTDIR}/BigWigs/133-90_ChIP_WT_H3K27me3_Rep1_S87_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/132-25_ChIP_ncu06787_H3K27me3_Rep_1_S22_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw -R "/scratch/ry00555/neurospora.bed" --skipZeros -o "${OUTDIR}/Matrices/ncu06788_H3K27me3_13225_3.gz"
plotHeatmap --matrixFile "${OUTDIR}/Matrices/ncu06788_H3K27me3_13225_3.gz" --outFileName "${OUTDIR}/Heatmaps/ncu06788_H3K27me3_13225_3.png" --samplesLabel WT 132-25   --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1

computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S ${OUTDIR}/BigWigs/133-90_ChIP_WT_H3K27me3_Rep1_S87_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/134-14_ChIP_ncu06788_H3K27me3_Rep2.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/132-25_ChIP_ncu06787_H3K27me3_Rep_1_S22_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw -R "/scratch/ry00555/neurospora.bed" --skipZeros -o "${OUTDIR}/Matrices/ncu06788_H3K27me3_Pool_3.gz"
plotHeatmap --matrixFile "${OUTDIR}/Matrices/ncu06788_H3K27me3_Pool_3.gz" --outFileName "${OUTDIR}/Heatmaps/ncu06788_H3K27me3_Pool_3.png" --samplesLabel WT 134-14 132-25  --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1

#K27me3 pooled
computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S ${OUTDIR}/BigWigs/133-90_ChIP_WT_H3K27me3_Rep1_S87_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/132-25_ChIP_ncu06787_H3K27me3_Rep_1_S22_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/134-14_ChIP_ncu06788_H3K27me3_Rep2.bin_25.smooth_75Bulk.bw -R "/scratch/ry00555/neurospora.bed" --skipZeros -o "${OUTDIR}/Matrices/ncu06788_H3K27me3_Pool2.gz"
plotHeatmap --matrixFile "${OUTDIR}/Matrices/ncu06788_H3K27me3_Pool2.gz" --outFileName "${OUTDIR}/Heatmaps/ncu06788_H3K27me3_Pool2.png" --samplesLabel WT 132-25 134-14  --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1

computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S ${OUTDIR}/BigWigs/133-90_ChIP_WT_H3K27me3_Rep1_S87_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/132-25_ChIP_ncu06787_H3K27me3_Rep_1_S22_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/134-14_ChIP_ncu06788_H3K27me3_Rep2.bin_25.smooth_75Bulk.bw -R "/scratch/ry00555/heatmapPRC2genes.bed" --skipZeros -o "${OUTDIR}/Matrices/ncu06788_H3K27me3_Pool_1stsetPRC2target_2.gz"
plotHeatmap --matrixFile "${OUTDIR}/Matrices/ncu06788_H3K27me3_Pool_1stsetPRC2target_2.gz" --outFileName "${OUTDIR}/Heatmaps/ncu06788_H3K27me3_Pool_1stsetPRC2target_2.png" --samplesLabel WT 132-25 134-14 --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1


#K36me3
computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S ${OUTDIR}/BigWigs/133-91_ChIP_WT_H3K36me3_Rep1_S88_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/132-27_ChIP_ncu06787_H3K36me3_Rep_1_S24_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/134-12_ChIP_ncu06787_H3K36me3_Rep2_S9_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw  -R "/scratch/ry00555/neurospora.bed" --skipZeros -o "${OUTDIR}/Matrices/ncu06788_H3K36me3_Pool.gz"
plotHeatmap --matrixFile "${OUTDIR}/Matrices/ncu06788_H3K36me3_Pool.gz" --outFileName "${OUTDIR}/Heatmaps/ncu06788_H3K27me3_Pool.png" --samplesLabel WT 132-27 134-12  --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1

computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S ${OUTDIR}/BigWigs/133-91_ChIP_WT_H3K36me3_Rep1_S88_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/132-27_ChIP_ncu06787_H3K36me3_Rep_1_S24_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/134-12_ChIP_ncu06787_H3K36me3_Rep2_S9_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw  -R "/scratch/ry00555/heatmapPRC2genes.bed" --skipZeros -o "${OUTDIR}/Matrices/ncu06788_H3K36me3_Pool_1stsetPRC2target_2.gz"
plotHeatmap --matrixFile "${OUTDIR}/Matrices/ncu06788_H3K36me3_Pool_1stsetPRC2target_2.gz" --outFileName "${OUTDIR}/Heatmaps/ncu06788_H3K36me3_Pool_1stsetPRC2target_2.png" --samplesLabel WT 132-27 134-12  --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1


#K27ac
#computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S ${OUTDIR}/BigWigs/131-53_ChIP_WT_H3K27ac_Rep1_S41_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw ${OUTDIR}/BigWigs/132-26_ChIP_ncu06787_H3K27ac_Rep_1_S23_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw -R "/scratch/ry00555/neurospora.bed" --skipZeros -o "${OUTDIR}/Matrices/ncu06788_H3K27ac_pool.gz"
#plotHeatmap --matrixFile "${OUTDIR}/Matrices/ncu06788_H3K27ac_pool.gz" --outFileName "${OUTDIR}/Heatmaps/ncu06788_H3K27ac_pool.png" --samplesLabel WT 132-26  --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1

# K416ac
#computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S ${OUTDIR}/BigWigs/132-17_ChIP_WT_H4K16ac_Rep_1_S14_R1_001_val_1.fq.gz.bin_25.smooth_75_MNase.bw ${OUTDIR}/BigWigs/132-28_ChIP_ncu06787_H4K16ac_Rep_1_S25_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw -R "/scratch/ry00555/neurospora.bed" --skipZeros -o "${OUTDIR}/Matrices/ncu06788_H4K16ac_pool.gz"
#plotHeatmap --matrixFile "${OUTDIR}/Matrices/ncu06788_H4K16ac_pool.gz" --outFileName "${OUTDIR}/Heatmaps/ncu06788_H4K16ac_pool.png" --samplesLabel WT 132-28  --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1
