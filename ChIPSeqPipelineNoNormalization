#!/bin/bash
#SBATCH --job-name=ChIPSeqPipeline
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=500gb
#SBATCH --time=48:00:00
#SBATCH --output=../ChIPSeqPipelinePooled.%j.out
#SBATCH --error=../ChIPSeqPipelinePooled.%j.err

OUTDIR=/scratch/ry00555/PoolChipSeq

#Load all the modules that are needed for the entire pipeline
ml SAMtools/1.16.1-GCC-10.2.0
ml BWA/0.7.17-GCC-10.2.0
ml Homer/4.11-foss-2020b
ml deepTools/3.5.1-foss-2020b-Python-3.8.6
# ml Perl/5.30.0-GCCcore-8.3.0
#ml Trim_Galore/0.6.5-GCCcore-8.3.0-Java-11-Python-3.7.4
ml BEDTools/2.30.0-GCC-10.2.0


#mkdir -p "${OUTDIR}/SortedBamFiles"
#mkdir -p "${OUTDIR}/BigWigs"
mkdir -p "${OUTDIR}/Peaks"
mkdir -p "$OUTDIR/HomerTagDirectories"
#mkdir -p "$OUTDIR/TdfFiles"
mkdir -p "${OUTDIR}/Heatmaps"
mkdir -p "${OUTDIR}/Matrices"
#mkdir -p "${OUTDIR}/NormalizedBigWigs"
mkdir -p "${OUTDIR}/Beds"
mkdir -p "${OUTDIR}/Counts"


PEAKDIR="${OUTDIR}/Peaks"
TAGDIR="${OUTDIR}/HomerTagDirectories"
BAMDIR="${OUTDIR}/SortedBamFiles"
BEDDIR="${OUTDIR}/Beds"


# Check if the sample_id contains "input" and set the -i argument accordingly
 if [[ "${sample_id}" == *input* ]]; then
   findPeaks "${TAGDIR}/${sample_id}" -style histone -region -size 150 -minDist 530 -o "${TAGDIR}/${sample_id}_peaks.txt" -i "${TAGDIR}/${sample_id}"
 else
   findPeaks "${TAGDIR}/${sample_id}" -style histone -region -size 150 -minDist 530 -o "${TAGDIR}/${sample_id}_peaks.txt"
 fi

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

  computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S 133-90_ChIP_WT_H3K27me3_Rep1_S87_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw 133-78_ChIP_NCU00423_H3K27me3_Rep1_S75_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw -R "/scratch/ry00555/neurospora.bed" --skipZeros -o "${OUTDIR}/Matrices/NCU00423_H3K27me3.gz"
  plotHeatmap --matrixFile "${OUTDIR}/Matrices/NCU00423_H3K27me3.gz" --outFileName "${OUTDIR}/Heatmaps/NCU00423_H3K27me3_hclust.png" --samplesLabel WT NCU00423-KO --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1

computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S 133-91_ChIP_WT_H3K36me3_Rep1_S88_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw 133-79_ChIP_NCU00423_H3K36me3_Rep1_S76_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw -R "/scratch/ry00555/neurospora.bed" --skipZeros -o "${OUTDIR}/Matrices/NCU00423_H3K36me3.gz"


  plotHeatmap --matrixFile "${OUTDIR}/Matrices/NCU00423_H3K36me3.gz" --outFileName "${OUTDIR}/Heatmaps/NCU00423_H3K36me3_hclust.png" --samplesLabel WT NCU00423-KO --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1


  # computeMatrix reference-point --referencePoint TSS -b 1000 -a 1000 -R H3K27me3_domains_5_boundary.bed -S 119_41_CUT_RUN_WT_H3K27me3_Rep1_R1_val_1.fq.gz.bin_25.smooth_75Bulk.bw 106_50_ChIP_hda1_H3K27me23_Rep1_R1_merged_trimmed.fq.gz.bin_25.smooth_75Bulk.bw 107_4_ChIP_hda1_H3K27me23_Rep2_R1_merged_val_1.fq.gz.bin_25.smooth_75Bulk.bw  129_20_ChIP_hda1_H3K27me3_Rep3_R1_val_1.fq.gz.bin_25.smooth_75Bulk.bw ChIP_dim5_H3K27me3_Rep1_R1_trimmed.fq.gz.bin_25.smooth_75Bulk.bw ChIP_dim5_H3K27me3_Rep2_R1_trimmed.fq.gz.bin_25.smooth_75Bulk.bw ChIP_dim5_H3K27me3_Rep3_R1_trimmed.fq.gz.bin_25.smooth_75Bulk.bw -o H3K27me3_domains_5matrix2.gz --outFileSortedRegions H3K27me3_domains_5-2.bed --skipZeros --verbose -p max/2  --missingDataAsZero
  #
  # plotHeatmap --matrixFile H3K27me3_domains_5matrix2.gz  -o H3K27me3_domains_5.png  --outFileSortedRegions sorted_clustered_H3K27me3_domains_5matrix2.gz --samplesLabel WT-Cu-1  hda1-Ch-1  hda-Cu-1 hda1-Ch-2 dim5-1 dim5-2 dim5-3 --colorMap Greens --sortRegions descend  --heatmapHeight 20  --missingDataColor white --sortUsingSamples 1 --kmeans 3 --zMax 15
