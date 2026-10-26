#!/bin/bash
#SBATCH --job-name=EPR1heatmap
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100gb
#SBATCH --time=36:00:00
#SBATCH --output=../EPR1heatmap.%j.out
#SBATCH --error=../EPR1heatmap.%j.err

ml deepTools

computeMatrix reference-point -p 12 -R "/home/ry00555/Research/Genomes/K27genes_NCU_n573_FINAL_fixed.bed" -S ../BigWigs/157-43_ChIP_tetRsuz12_36hr_H3K27me3_Rep2.bin_25.smooth_75Bulk.bw ../BigWigs/157-47_ChIP_tetRsuz12_36hr_H3K36me3_Rep2.bin_25.smooth_75Bulk.bw ../BigWigs/157-53_ChIP_epr1_tetRsuz12_36hr_H3K27me3_Rep3.bin_25.smooth_75Bulk.bw ../BigWigs/157-57_ChIP_epr1_tetRsuz12_36hr_H3K36me3_Rep3.bin_25.smooth_75Bulk.bw	-o EPR1_K27Genes_K27me3_K36me3_V1.mat 	--sortRegions keep --missingDataAsZero -bs 10 -a 2000 -b 1000

plotHeatmap -m EPR1_K27Genes_K27me3_K36me3_V1.mat  -o EPR1_K27Genes_K27me3_K36me3_V2.png --sortRegions descend  --sortUsingSamples 3 4  --outFileSortedRegions EPR1_K27Genes_K27me3_K36me3_V2.bed --startLabel "5'" -max 15 8 15 8 --endLabel "3'"  --boxAroundHeatmaps no  --samplesLabel "WT K27me3" "WT K36me3" "epr1 K27me3" "epr1 H3K36me3" --colorMap 'Greens' 'YlOrBr' 'Greens' 'YlOrBr' --sortUsing mean --heatmapHeight 8  --heatmapWidth 4


157-40_ChIP_tetRsuz12_0hr_H3K27me3_Rep2.bin_25.smooth_75Bulk.bw
157-41_ChIP_tetRsuz12_12hr_H3K27me3_Rep2.bin_25.smooth_75Bulk.bw
157-42_ChIP_tetRsuz12_24hr_H3K27me3_Rep2.bin_25.smooth_75Bulk.bw
157-43_ChIP_tetRsuz12_36hr_H3K27me3_Rep2.bin_25.smooth_75Bulk.bw
157-44_ChIP_tetRsuz12_0hr_H3K36me3_Rep2.bin_25.smooth_75Bulk.bw
157-45_ChIP_tetRsuz12_12hr_H3K36me3_Rep2.bin_25.smooth_75Bulk.bw
157-46_ChIP_tetRsuz12_24hr_H3K36me3_Rep2.bin_25.smooth_75Bulk.bw
157-47_ChIP_tetRsuz12_36hr_H3K36me3_Rep2.bin_25.smooth_75Bulk.bw
157-48_ChIP_epr1_tetRsuz12_0hr_Input_Rep3.bin_25.smooth_75Bulk.bw
157-49_ChIP_epr1_tetRsuz12_36hr_Input_Rep3.bin_25.smooth_75Bulk.bw
157-50_ChIP_epr1_tetRsuz12_0hr_H3K27me3_Rep3.bin_25.smooth_75Bulk.bw
157-51_ChIP_epr1_tetRsuz12_12hr_H3K27me3_Rep3.bin_25.smooth_75Bulk.bw
157-52_ChIP_epr1_tetRsuz12_24hr_H3K27me3_Rep3.bin_25.smooth_75Bulk.bw
#157-53_ChIP_epr1_tetRsuz12_36hr_H3K27me3_Rep3.bin_25.smooth_75Bulk.bw
157-54_ChIP_epr1_tetRsuz12_0hr_H3K36me3_Rep3.bin_25.smooth_75Bulk.bw
157-55_ChIP_epr1_tetRsuz12_12hr_H3K36me3_Rep3.bin_25.smooth_75Bulk.bw
157-56_ChIP_epr1_tetRsuz12_24hr_H3K36me3_Rep3.bin_25.smooth_75Bulk.bw
157-57_ChIP_epr1_tetRsuz12_36hr_H3K36me3_Rep3.bin_25.smooth_75Bulk.bw
