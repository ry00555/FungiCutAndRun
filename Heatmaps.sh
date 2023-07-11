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

OUTDIR="/scratch/ry00555/HeatMapPeaks"

module load deepTools/3.5.1-intel-2020b-Python-3.8.6

#conda install -c bioconda deeptools

computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S $OUTDIR/124_1_ChIP_WT_K4me2_Rep1_4901.bin_25.smooth_50_BPM.bw $OUTDIR/124_2_ChIP_cac_1_K4me2_Rep1_4901.bin_25.smooth_50_BPM.bw $OUTDIR/124_3_ChIP_cac_2_K4me2_Rep1_4901.bin_25.smooth_50_BPM.bw $OUTDIR/124_4_ChIP_cac_3_K4me2_Rep1_4901.bin_25.smooth_50_BPM.bw $OUTDIR/124_5_ChIP_set_7_K4me2_Rep1_4901.bin_25.smooth_50_BPM.bw -R heatmapPRC2genes.bed --skipZeros -o matrix_PRC2.gz

plotHeatmap -m matrix_PRC2.gz -out cacheatmap_H3K4me2_Run124_hclust.png --samplesLabel WT cac-1 cac-2 cac-3 set-7 --hclust 1 --colorMap Reds


#133-23_ChIP_WT_H3K9me3_Rep2.bin_25.smooth_75.bw
#133-24_ChIP_cac-1_H3K9me3_Rep2.bin_25.smooth_75.bw
#133-25_ChIP_cac-2_H3K9me3_Rep2.bin_25.smooth_75.bw

#129-90_ChIP_WT_H3K36me3_Rep1_S71_L001_R1_001_val_1.fq.gz.bin_25.smooth_75.bw
#129-91_ChIP_cac-1_H3K36me3_Rep1_S72_L001_R1_001_val_1.fq.gz.bin_25.smooth_75.bw
#129-92_ChIP_cac-2_H3K36me3_Rep1_S73_L001_R1_001_val_1.fq.gz.bin_25.smooth_75.bw
#129-93_ChIP_cac-3_H3K36me3_Rep1_S74_L001_R1_001_val_1.fq.gz.bin_25.smooth_75.bw
#129-94_ChIP_set-7_H3K36me3_Rep1_S75_L001_R1_001_val_1.fq.gz.bin_25.smooth_75.bw

#124_1_ChIP_WT_K4me2_Rep1_4901.bin_25.smooth_50_BPM.bw
#124_2_ChIP_cac_1_K4me2_Rep1_4901.bin_25.smooth_50_BPM.bw
#124_3_ChIP_cac_2_K4me2_Rep1_4901.bin_25.smooth_50_BPM.bw
#124_4_ChIP_cac_3_K4me2_Rep1_4901.bin_25.smooth_50_BPM.bw
#124_5_ChIP_set_7_K4me2_Rep1_4901.bin_25.smooth_50_BPM.bw
