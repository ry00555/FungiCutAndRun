#!/bin/bash
#SBATCH --job-name=PeakCalling
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=24:00:00
#SBATCH --output=../callpeaks.%j.out
#SBATCH --error=../callpeaks.%j.err

#load modules
module load Homer/4.11-foss-2019b SAMtools/1.16.1-GCC-11.3.0

#Make directory where you want your called peaks file to go at the end, I recommend the "Peaks" directory created after running MapCutAndRun script
OUTDIR="/scratch/ry00555/Run133/Peaks"

#Make directory for tag files, will need to change Run# for your set of files
TAGDIR="/scratch/ry00555/Run133/homer/tag"

if [ ! -d $TAGDIR ]
then
    mkdir -p $TAGDIR
fi

#Make directory for bam files, This should exist under whichever Run directory you are using if you ran the MapCutAndRun script
BAMDIR="/scratch/ry00555/Run133/SortedBamFiles"

#Make tag directories for each of your bam files, you will need to change the bam file used for each line of code, and the directory if you want to separate your tag files
makeTagDirectory ${TAGDIR}/129-33_ChIP_WT_K27me3_AM_Rep_2 ${BAMDIR}/129-33_ChIP_WT_K27me3_AM_Rep_2_S32_L001_R1_001_val_1.fq.gz.bam
makeTagDirectory ${TAGDIR}/129-34_ChIP_cac-1_K27me3_AM_Rep_2 ${BAMDIR}/129-34_ChIP_cac-1_K27me3_AM_Rep_2_S33_L001_R1_001_val_1.fq.gz.bam
makeTagDirectory ${TAGDIR}/129-35_ChIP_cac-2_K27me3_AM_Rep_2 ${BAMDIR}/129-35_ChIP_cac-2_K27me3_AM_Rep_2_S34_L001_R1_001_val_1.fq.gz.bam
makeTagDirectory ${TAGDIR}/129-36_ChIP_cac-3_K27me3_AM_Rep_2 ${BAMDIR}/129-36_ChIP_cac-3_K27me3_AM_Rep_2_S35_L001_R1_001_val_1.fq.gz.bam
makeTagDirectory ${TAGDIR}/129-37_ChIP_set-7_K27me3_AM_Rep_2 ${BAMDIR}/129-37_ChIP_set-7_K27me3_AM_Rep_2_S36_L001_R1_001_val_1.fq.gz.bam
makeTagDirectory ${TAGDIR}/129-38_ChIP_WT_K27me3_AbC_Rep_1 ${BAMDIR}/129-38_ChIP_WT_K27me3_AbC_Rep_1_S37_L001_R1_001_val_1.fq.gz.bam
makeTagDirectory ${TAGDIR}/129-39_ChIP_cac-1_K27me3_AbC_Rep_1 ${BAMDIR}/129-39_ChIP_cac-1_K27me3_AbC_Rep_1_S38_L001_R1_001_val_1.fq.gz.bam
makeTagDirectory ${TAGDIR}/129-40_ChIP_cac-2_K27me3_AbC_Rep_1 ${BAMDIR}/129-40_ChIP_cac-2_K27me3_AbC_Rep_1_S39_L001_R1_001_val_1.fq.gz.bam
makeTagDirectory ${TAGDIR}/129-41_ChIP_cac-3_K27me3_AbC_Rep_1 ${BAMDIR}/129-41_ChIP_cac-3_K27me3_AbC_Rep_1_S40_L001_R1_001_val_1.fq.gz.bam
makeTagDirectory ${TAGDIR}/129-42_ChIP_set-7_K27me3_AbC_Rep_1 ${BAMDIR}/129-42_ChIP_set-7_K27me3_AbC_Rep_1_S41_L001_R1_001_val_1.fq.gz.bam
makeTagDirectory ${TAGDIR}/129-43_ChIP_WT_input ${BAMDIR}/129-43_ChIP_WT_input_S42_L001_R1_001_val_1.fq.gz.bam
makeTagDirectory ${TAGDIR}/129-44_ChIP_cac-1_input ${BAMDIR}/129-44_ChIP_cac-1_input_S43_L001_R1_001_val_1.fq.gz.bam
makeTagDirectory ${TAGDIR}/129-45_ChIP_cac-2_input ${BAMDIR}/129-45_ChIP_cac-2_input_S44_L001_R1_001_val_1.fq.gz.bam
makeTagDirectory ${TAGDIR}/129-46_ChIP_cac-3_input ${BAMDIR}/129-46_ChIP_cac-3_input_S45_L001_R1_001_val_1.fq.gz.bam
makeTagDirectory ${TAGDIR}/129-47_ChIP_set-7_input ${BAMDIR}/129-47_ChIP_set-7_input_S46_L001_R1_001_val_1.fq.gz.bam

#Run findPeaks command on tag files to call peaks, Use respective input files for each sample (i.e. WT_input -> WT_K27)
findPeaks ${TAGDIR}/129-33_ChIP_WT_K27me3_AM_Rep_2 -style histone -region -size 150 -minDist 530 -o ${OUTDIR}/129-33_ChIP_WT_K27me3_AM_Rep_2_peaks.txt -i ${TAGDIR}/129-43_ChIP_WT_input
findPeaks ${TAGDIR}/129-34_ChIP_cac-1_K27me3_AM_Rep_2 -style histone -region -size 150 -minDist 530 -o ${OUTDIR}/129-34_ChIP_cac-1_K27me3_AM_Rep_2_peaks.txt -i ${TAGDIR}/129-44_ChIP_cac-1_input
findPeaks ${TAGDIR}/129-35_ChIP_cac-2_K27me3_AM_Rep_2 -style histone -region -size 150 -minDist 530 -o ${OUTDIR}/129-35_ChIP_cac-2_K27me3_AM_Rep_2_peaks.txt -i ${TAGDIR}/129-45_ChIP_cac-2_input
findPeaks ${TAGDIR}/129-36_ChIP_cac-3_K27me3_AM_Rep_2 -style histone -region -size 150 -minDist 530 -o ${OUTDIR}/129-36_ChIP_cac-3_K27me3_AM_Rep_2_peaks.txt -i ${TAGDIR}/129-46_ChIP_cac-3_input
findPeaks ${TAGDIR}/129-37_ChIP_set-7_K27me3_AM_Rep_2 -style histone -region -size 150 -minDist 530 -o ${OUTDIR}/129-37_ChIP_set-7_K27me3_AM_Rep_2_peaks.txt -i ${TAGDIR}/129-47_ChIP_set-7_input
findPeaks ${TAGDIR}/129-38_ChIP_WT_K27me3_AbC_Rep_1 -style histone -region -size 150 -minDist 530 -o ${OUTDIR}/129-38_ChIP_WT_K27me3_AbC_Rep_1_peaks.txt -i ${TAGDIR}/129-43_ChIP_WT_input
findPeaks ${TAGDIR}/129-39_ChIP_cac-1_K27me3_AbC_Rep_1 -style histone -region -size 150 -minDist 530 -o ${OUTDIR}/129-39_ChIP_cac-1_K27me3_AbC_Rep_1_peaks.txt -i ${TAGDIR}/129-44_ChIP_cac-1_input
findPeaks ${TAGDIR}/129-40_ChIP_cac-2_K27me3_AbC_Rep_1 -style histone -region -size 150 -minDist 530 -o ${OUTDIR}/129-40_ChIP_cac-2_K27me3_AbC_Rep_1_peaks.txt -i ${TAGDIR}/129-45_ChIP_cac-2_input
findPeaks ${TAGDIR}/129-41_ChIP_cac-3_K27me3_AbC_Rep_1 -style histone -region -size 150 -minDist 530 -o ${OUTDIR}/129-41_ChIP_cac-3_K27me3_AbC_Rep_1_peaks.txt -i ${TAGDIR}/129-46_ChIP_cac-3_input
findPeaks ${TAGDIR}/129-42_ChIP_set-7_K27me3_AbC_Rep_1 -style histone -region -size 150 -minDist 530 -o ${OUTDIR}/129-42_ChIP_set-7_K27me3_AbC_Rep_1_peaks.txt -i ${TAGDIR}/129-47_ChIP_set-7_input
