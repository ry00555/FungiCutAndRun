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
ml BWA/0.7.17-GCCcore-11.3.0
ml SAMtools/1.16.1-GCC-11.3.0
ml Homer/4.11-foss-2022a
ml deepTools/3.5.2-foss-2022a
# ml Perl/5.34.1-GCCcore-11.3.0
#ml Trim_Galore/0.6.5-GCCcore-8.3.0-Java-11-Python-3.7.4
ml BEDTools/2.30.0-GCC-11.3.0


#mkdir -p "${OUTDIR}/SortedBamFiles"
#mkdir -p "${OUTDIR}/BigWigs"
mkdir -p "${OUTDIR}/Peaks"
mkdir -p "$OUTDIR/HomerTagDirectories"
#mkdir -p "$OUTDIR/TdfFiles"
#mkdir -p "${OUTDIR}/Heatmaps"
#mkdir -p "${OUTDIR}/Matrices"
#mkdir -p "${OUTDIR}/NormalizedBigWigs"
#mkdir -p "${OUTDIR}/Beds"
#mkdir -p "${OUTDIR}/Counts"


PEAKDIR="${OUTDIR}/Peaks"
TAGDIR="${OUTDIR}/HomerTagDirectories"
BAMDIR="${OUTDIR}/SortedBamFiles"
BEDDIR="${OUTDIR}/Beds"

for bam_file in "${BAMDIR}"/*.bam; do
  # Get the sample ID from the BAM file name
  sample_id=$(basename "${bam_file}" .bam)
  # Remove everything after "Rep_1" in the sample ID
  sample_id="${sample_id%%_Rep_1*}"
# Check if the sample_id contains "input" and set the -i argument accordingly
 if [[ "${sample_id}" == *input* ]]; then
   findPeaks "${TAGDIR}/${sample_id}" -style histone -region -size 150 -minDist 530 -o "${PEAKDIR}/${sample_id}_peaks.txt" -i "${TAGDIR}/${sample_id}"
 else
   findPeaks "${TAGDIR}/${sample_id}" -style histone -region -size 150 -minDist 530 -o "${PEAKDIR}/${sample_id}_peaks.txt"
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
done
  computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S 133-90_ChIP_WT_H3K27me3_Rep1_S87_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw 133-78_ChIP_NCU00423_H3K27me3_Rep1_S75_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw -R "/scratch/ry00555/neurospora.bed" --skipZeros -o "${OUTDIR}/Matrices/NCU00423_H3K27me3.gz"
  plotHeatmap --matrixFile "${OUTDIR}/Matrices/NCU00423_H3K27me3.gz" --outFileName "${OUTDIR}/Heatmaps/NCU00423_H3K27me3_hclust.png" --samplesLabel WT NCU00423-KO --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1

computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S 133-91_ChIP_WT_H3K36me3_Rep1_S88_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw 133-79_ChIP_NCU00423_H3K36me3_Rep1_S76_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw -R "/scratch/ry00555/neurospora.bed" --skipZeros -o "${OUTDIR}/Matrices/NCU00423_H3K36me3.gz"


  plotHeatmap --matrixFile "${OUTDIR}/Matrices/NCU00423_H3K36me3.gz" --outFileName "${OUTDIR}/Heatmaps/NCU00423_H3K36me3_hclust.png" --samplesLabel WT NCU00423-KO --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1
