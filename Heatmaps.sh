#!/bin/bash
#SBATCH --job-name=HeatMapPeaks
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=../HeatMapPeaks.%j.out
#SBATCH --error=../HeatMapPeaks.%j.err

OUTDIR="/scratch/ry00555/OutputRun137/Neurospora_Output/BigWigs"
MATRICESDIR="/scratch/ry00555/OutputRun137/Neurospora_Output/Matrices"
HEATMAPDIR="/scratch/ry00555/OutputRun137/Neurospora_Output/Heatmaps"
ml deepTools

#conda install -c bioconda deeptools

computeMatrix reference-point --referencePoint center -b 1500 -a 1500 -S $OUTDIR/137-22_CUTANDRUN_WT_H3K27me3_Rep1_S22_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw $OUTDIR/137-22_CUTANDRUN_WT_H3K27me3_Rep1_S22_R1_001_val_1.fq.gz.bin_25.smooth_75_MNase.bw $OUTDIR/137-2_CUTANDRUN_WT_H3K27me3_Rep1_S2_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw $OUTDIR/137-2_CUTANDRUN_WT_H3K27me3_Rep1_S2_R1_001_val_1.fq.gz.bin_25.smooth_75_MNase.bw $OUTDIR/137-25_CUTANDRUN_set-7_H3K27me3_Rep1_S25_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw $OUTDIR/137-25_CUTANDRUN_set-7_H3K27me3_Rep1_S25_R1_001_val_1.fq.gz.bin_25.smooth_75_MNase.bw $OUTDIR/137-7_CUTANDRUN_set-7_H3K27me3_Rep1_S7_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw $OUTDIR/137-7_CUTANDRUN_set-7_H3K27me3_Rep1_S7_R1_001_val_1.fq.gz.bin_25.smooth_75_MNase.bw -R /scratch/ry00555/heatmapPRC2genes.bed --skipZeros -o $MATRICESDIR/matrix_PRC2Genes.gz

plotHeatmap -m $MATRICESDIR/matrix_PRC2Genes.gz -out $HEATMAPDIR/ZLRun137CutandRunWTset7Only_Prc2Genes.png --samplesLabel WT-22Bulk WT-22Mnase WT2-Bulk WT2-Mnase set7-25Bulk set7-25MNase set7-7Bulk set7-7MNase --hclust 1 --colorMap Reds

computeMatrix reference-point --referencePoint center -b 1500 -a 1500 -S $OUTDIR/137-22_CUTANDRUN_WT_H3K27me3_Rep1_S22_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw $OUTDIR/137-22_CUTANDRUN_WT_H3K27me3_Rep1_S22_R1_001_val_1.fq.gz.bin_25.smooth_75_MNase.bw $OUTDIR/137-2_CUTANDRUN_WT_H3K27me3_Rep1_S2_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw $OUTDIR/137-2_CUTANDRUN_WT_H3K27me3_Rep1_S2_R1_001_val_1.fq.gz.bin_25.smooth_75_MNase.bw $OUTDIR/137-25_CUTANDRUN_set-7_H3K27me3_Rep1_S25_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw $OUTDIR/137-25_CUTANDRUN_set-7_H3K27me3_Rep1_S25_R1_001_val_1.fq.gz.bin_25.smooth_75_MNase.bw $OUTDIR/137-7_CUTANDRUN_set-7_H3K27me3_Rep1_S7_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw $OUTDIR/137-7_CUTANDRUN_set-7_H3K27me3_Rep1_S7_R1_001_val_1.fq.gz.bin_25.smooth_75_MNase.bw -R /scratch/ry00555/WT_H3K27me3peaks.bed --skipZeros -o $MATRICESDIR/matrix_PRC2Domains.gz

plotHeatmap -m $MATRICESDIR/matrix_PRC2Domains.gz -out $HEATMAPDIR/ZLRun137CutandRunWTset7Only_Prc2Domains.png --samplesLabel WT-22Bulk WT-22Mnase WT2-Bulk WT2-Mnase set7-25Bulk set7-25MNase set7-7Bulk set7-7MNase --hclust 1 --colorMap Reds
