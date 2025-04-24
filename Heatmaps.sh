#!/bin/bash
#SBATCH --job-name=Heatmaps
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --time=08:00:00
#SBATCH --output=Heatmaps.%j.out
#SBATCH --error=Heatmaps.%j.err

cd $SLURM_SUBMIT_DIR

ml deepTools/3.5.2-foss-2022a

computeMatrix scale-regions -p 12 --startLabel "5'" --endLabel "3'" \
-S /scratch/ry00555/RTT109PaperFigures/BigWigs/131-22_ChIP_WT_H3K27me2me3_Rep1.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/131-52_ChIP_WT_H3K27me3_Rep1.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/131-67_ChIP_WT_H3K27me2me3_Rep2.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/131-87_ChIP_WT_H3K27me2me3_Rep3.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/135-18_ChIP_WT_H3K27me2_3.bin_25.smooth_75Bulk.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/147-3_ChIP_S2_WT_H3K27me3_Rep5_Nc_24hrVMMON_S3_L001_R1_001_val_1.fq.gz.bin_25.smooth_50_Q30.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/134-28_ChIP_WT_H3K27me3_Rep2.bin_25.smooth_75Bulk.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/134-2_ChIP_WT_H3K27me3_Rep2.bin_25.smooth_75Bulk.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/135-27_ChIP_WT_H3K27me3_Rep1.bin_25.smooth_75Bulk.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/137-27_ChIP_WT_H3K27me3_Rep3_S27_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/142-10_ChIP_WT_H3K27me3_Rep3.bin_25.smooth_50_Q30.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/142-76_ChIP_WT_H3K27me3__S76_L007_R1_001_val_1.fq.gz.bin_25.smooth_50_Q30.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/137-28_ChIP_WT_H3K36me3_Rep3_S28_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/142-11_ChIP_WT_H3K36me3_Rep3.bin_25.smooth_50_Q30.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/142-77_ChIP_WT_H3K36me3__S77_L007_R1_001_val_1.fq.gz.bin_25.smooth_50_Q30.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/131-27_ChIP_WT_H3K36me2_Rep1.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/131-54_ChIP_WT_H3K36me3_Rep1.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/131-72_ChIP_WT_H3K36me2_Rep2.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/131-37_ChIP_WT_input_Rep1.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/131-77_ChIP_WT_input_Rep2.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/131-92_ChIP_WT_input_Rep1.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/132-17_ChIP_WT_H4K16ac.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/135-16_ChIP_WT_H3K9me3.bin_25.smooth_75Bulk.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/135-23_ChIP_WT_H3K9me3.bin_25.smooth_75Bulk.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/135-24_ChIP_WT_H3K9me3.bin_25.smooth_75Bulk.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/135-25_ChIP_WT_H3K9me3.bin_25.smooth_75Bulk.bw -R "/scratch/ry00555/neurospora.bed" --skipZeros -o /scratch/ry00555/RTT109PaperFigures/Matrices/AllWTallgenes.gz --outFileNameMatrix /scratch/ry00555/RTT109PaperFigures/Matrices/AllWTallgenes.matrix.txt --sortRegions keep --missingDataAsZero -bs 10

plotHeatmap -m /scratch/ry00555/RTT109PaperFigures/Matrices/AllWTallgenes.gz -out /scratch/ry00555/RTT109PaperFigures/Heatmaps/AllWTAllModsAllGenes_V1_nosorting.png --samplesLabel 131-22-H3K27me2me3 131-52-H3K27me3 131-67-H3K27me2me3 131-87-H3K27me2me3 135-18-H3K27me2 147-3-H3K27me3 134-28-H3K27me3 134-2-H3K27me3 135-27-H3K27me3 137-27-H3K27me3 142-10-H3K27me3 142-76-H3K27me3 137-28-H3K36me3 142-11-H3K36me3 142-77-H3K36me3 131-27-H3K36me2 131-54-H3K36me3 131-72-H3K36me2 131-37-input 131-77-input 131-92-input 132-17-H4K16ac 135-16-H3K9me3 135-23-H3K9me3 135-24-H3K9me3 135-25-H3K9me3 --hclust 1
