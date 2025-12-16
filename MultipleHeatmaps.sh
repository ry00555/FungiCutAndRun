#!/bin/bash
#SBATCH --job-name=Heatmaps_MultiBed
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --time=08:00:00
#SBATCH --output=../MultipleHeatmaps.%j.out
#SBATCH --error=../MultipleHeatmaps.%j.err
cd $SLURM_SUBMIT_DIR

OUTDIR="/scratch/ry00555/RTT109PaperFigures"
GENEDIR="/home/ry00555/Research/Genomes/HeatmapGeneFiles"
BWDIR="$OUTDIR/PreviouslyMappedBigWigs"

#mkdir -p "$OUTDIR/Heatmaps"

ml deepTools

 #computeMatrix reference-point -p 12 -R "$GENEDIR/K27genes.bed" -S $BWDIR/147-3_WT_H3K27me3_rep9.bin_25.smooth_50.bw $BWDIR/147-59_rtt109_H3K27me3_rep9.bin_25.smooth_50.bw $BWDIR/150-65_rtt109_H3K27me3_rep11.bin_25.smooth_50.bw $BWDIR/150-81_rtt109flag_H3K27me3_rep3.bin_25.smooth_50.bw	-o $OUTDIR/Heatmaps/RTT109_H3K27me3_K27genes_V2.mat 	--sortRegions keep --missingDataAsZero -bs 10 -a 2000 -b 1000

 #plotHeatmap -m $OUTDIR/Heatmaps/RTT109_H3K27me3_K27genes_V2.mat  -o $OUTDIR/Heatmaps/RTT109_H3K27me3_K27genes_V4.png --sortRegions descend  --sortUsingSamples 1 4  --outFileSortedRegions $OUTDIR/Heatmaps/RTT109_H3K27me3_K27genes_sortOrder_V2.bed --startLabel "5'" -max 50 40 40 40 --endLabel "3'"  --boxAroundHeatmaps no  --samplesLabel "WT K27me3" "∆rtt109 K27me3" "∆rtt109 K27me3" "rtt109-3xflag H3K27me3" --colorMap 'Greens' --sortUsing mean --heatmapHeight 8  --heatmapWidth 4

 #computeMatrix reference-point -p 12 -R "$GENEDIR/K27genes.bed" -S  "$BWDIR/147-3_WT_H3K27me3_rep9.bin_25.smooth_50.bw" "$BWDIR/147-47_rtt109_H3K27me3_rep7.bin_25.smooth_50.bw"  "$BWDIR/150-81_rtt109flag_H3K27me3_rep3.bin_25.smooth_50.bw" "$BWDIR/142-88_rtt109flag_H3K27me3_rep1.bin_25.smooth_50.bw" "$BWDIR/150-77_rtt109flag_H3K27me3_rep2.bin_25.smooth_50.bw" -o $OUTDIR/Heatmaps/RTT109_H3K27me3_K27genes_FilteredReps_V3.mat --sortRegions keep --missingDataAsZero -bs 10 -a 2000 -b 1000

#plotHeatmap -m $OUTDIR/Heatmaps/RTT109_H3K27me3_K27genes_FilteredReps_V3.mat -o $OUTDIR/Heatmaps/RTT109_H3K27me3_K27genes_FilteredReps_V3.png --sortRegions descend --sortUsingSamples 1 --heatmapHeight 8  --heatmapWidth 4  --outFileSortedRegions $OUTDIR/Heatmaps/RTT109_H3K27me3_K27genes_FilteredReps_V3.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --samplesLabel  "WT K27me3" "∆rtt109 K27me3" "rtt109-3xflag K27me3" "rtt109-3xflag K27me3" "rtt109-3xflag K27me3"  --colorMap 'Greens'



 #computeMatrix scale-regions -p 12 -R "$GENEDIR/Figure2_K9_Peaks.txt" -S $BWDIR/149-104_WT_H3K9me3_rep5.bin_25.smooth_50.bw  $BWDIR/149-43_rtt109_H3K9me3_rep3.bin_25.smooth_50.bw  $BWDIR/150-80_rtt109flag_H3K9me3_rep2.bin_25.smooth_50.bw --sortRegions keep --missingDataAsZero -bs 10	-b 500 -a 500  	-o $OUTDIR/Heatmaps/RTT109_H3K9me3_K9Peaks_V2.mat

 #plotHeatmap -m $OUTDIR/Heatmaps/RTT109_H3K9me3_K9Peaks_V2.mat -o $OUTDIR/Heatmaps/RTT109_H3K9me3_K9Peaks_V3.png --sortRegions descend --sortUsingSamples 1 --heatmapHeight 8  --heatmapWidth 4  --outFileSortedRegions $OUTDIR/Heatmaps/RTT109_H3K9me3_K9Peaks_sortOrder_V2.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --samplesLabel "WT K9me3" "∆rtt109 K9me3" "rtt109-3xflag K9me3"  --colorMap 'Blues'

 #computeMatrix scale-regions -p 12 -R "$GENEDIR/Figure2_K9_Peaks.txt" -S ""$BWDIR/135-16_WT_H3K9me3_rep1.bin_25.smooth_50.bw"  ""$BWDIR/148-117_WT_H3K9me3_rep4.bin_25.smooth_50.bw" ""$BWDIR/135-23_WT_H3K9me3_rep2.bin_25.smooth_50.bw" ""$BWDIR/149-104_WT_H3K9me3_rep5.bin_25.smooth_50.bw" ""$BWDIR/135-25_WT_H3K9me3_rep3.bin_25.smooth_50.bw"  ""$BWDIR/150-60_WT_H3K9me3_rep6.bin_25.smooth_50.bw" ""$BWDIR/111-5_rtt109_H3K9me3_rep1.bin_25.smooth_50.bw"   ""$BWDIR/150-64_rtt109_H3K9me3_rep4.bin_25.smooth_50.bw" ""$BWDIR/135-81_rtt109_H3K9me3_rep2.bin_25.smooth_50.bw"  ""$BWDIR/150-68_rtt109_H3K9me3_rep5.bin_25.smooth_50.bw" ""$BWDIR/149-43_rtt109_H3K9me3_rep3.bin_25.smooth_50.bw"  ""$BWDIR/150-72_rtt109_H3K9me3_rep6.bin_25.smooth_50.bw" --sortRegions keep --missingDataAsZero -bs 10	-b 500 -a 500  	-o $OUTDIR/Heatmaps/RTT109_H3K9me3_K9Peaks_AllReps.mat

 #plotHeatmap -m $OUTDIR/Heatmaps/RTT109_H3K9me3_K9Peaks_AllReps.mat -o $OUTDIR/Heatmaps/RTT109_H3K9me3_K9Peaks_AllReps.png --sortRegions descend --sortUsingSamples 1 --heatmapHeight 8  --heatmapWidth 4  --outFileSortedRegions $OUTDIR/Heatmaps/RTT109_H3K9me3_K9Peaks_AllReps.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --samplesLabel "WT K9me3" "WT K9me3" "WT K9me3" "WT K9me3" "WT K9me3" "WT K9me3" "∆rtt109 K9me3" "∆rtt109 K9me3" "∆rtt109 K9me3" "∆rtt109 K9me3" "∆rtt109 K9me3" "∆rtt109 K9me3"  --colorMap 'Blues'


#computeMatrix scale-regions -p 12 -R "$GENEDIR/Figure2_K9_Peaks.txt" -S  $BWDIR/150-60_WT_H3K9me3_rep6.bin_25.smooth_50.bw $BWDIR/150-64_rtt109_H3K9me3_rep4.bin_25.smooth_50.bw $BWDIR/150-80_rtt109flag_H3K9me3_rep2.bin_25.smooth_50.bw --sortRegions keep --missingDataAsZero -bs 10	-b 500 -a 500  	-o $OUTDIR/Heatmaps/RTT109_H3K9me3_K9Peaks_V3.mat

#plotHeatmap -m $OUTDIR/Heatmaps/RTT109_H3K9me3_K9Peaks_V3.mat -o $OUTDIR/Heatmaps/RTT109_H3K9me3_K9Peaks_V4.png --sortRegions descend --sortUsingSamples 3 --heatmapHeight 8  --heatmapWidth 4  --outFileSortedRegions $OUTDIR/Heatmaps/RTT109_H3K9me3_K9Peaks_sortOrder_V3.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --samplesLabel "WT K9me3" "∆rtt109 K9me3" "rtt109-3xflag K9me3"  --colorMap 'Blues'


#computeMatrix reference-point -p 12 -R "$GENEDIR/K27genes.bed" -S $BWDIR/153-124_ChIP_WT_H3K27me3_.bin_25.smooth_50Bulk.bw $BWDIR/153-126_ChIP_rtt109_H3K27me3_.bin_25.smooth_50Bulk.bw $BWDIR/150-62_WT_H3K36me3_rep8.bin_25.smooth_50.bw $BWDIR/150-70_rtt109_H3K36me3_rep7.bin_25.smooth_50.bw  -o $OUTDIR/Heatmaps/K27_K36onK27genes.gz --outFileNameMatrix $OUTDIR/Heatmaps/K27_K36onK27genes.tab 	--sortRegions keep  --missingDataAsZero -bs 10 -a 2000 -b 1000


#plotHeatmap -m $OUTDIR/Heatmaps/K27_K36onK27genes.gz -o $OUTDIR/Heatmaps/K27_K36onK27genes.png --sortRegions descend --sortUsingSamples 2 --outFileSortedRegions $OUTDIR/Heatmaps/K27_K36onK27genes_sorted.bed --startLabel "5'" --endLabel "3'" --heatmapHeight 8 --heatmapWidth 4  --boxAroundHeatmaps no --samplesLabel "WT 153-124" "∆rtt109 153-126"  "WT K36me3" "∆rtt109 K36me3"  --colorMap  'Greens' 'Greens' 'YlOrBr' 'YlOrBr' --hclust 2


#computeMatrix reference-point -p 12 -R "$OUTDIR/Heatmaps/K27_K36_onK27genes_cluster2_V1_sorted_filtered2_again2.bed" -S $BWDIR/153-124_ChIP_WT_H3K27me3_.bin_25.smooth_50Bulk.bw $BWDIR/153-126_ChIP_rtt109_H3K27me3_.bin_25.smooth_50Bulk.bw $BWDIR/150-62_WT_H3K36me3_rep8.bin_25.smooth_50.bw $BWDIR/150-70_rtt109_H3K36me3_rep7.bin_25.smooth_50.bw  -o $OUTDIR/Heatmaps/K27_K36_onK27genes_cluster2_V1_sorted_filtered2_again2.gz --outFileNameMatrix $OUTDIR/Heatmaps/K27_K36_onK27genes_cluster2_V1_sorted_filtered2_again2.tab 	--sortRegions keep  --missingDataAsZero -bs 10 -a 2000 -b 1000

#plotHeatmap -m $OUTDIR/Heatmaps/K27_K36_onK27genes_cluster2_V1_sorted_filtered2_again2.gz -o $OUTDIR/Heatmaps/K27_K36_onK27genes_cluster2_V1_sorted_filtered2_again2.png --sortRegions descend --sortUsingSamples 2 --outFileSortedRegions $OUTDIR/Heatmaps/K27_K36_onK27genes_cluster2_V1_sorted_filtered2_again3.bed --startLabel "5'" --endLabel "3'" --heatmapHeight 8 --heatmapWidth 4  --boxAroundHeatmaps no --samplesLabel "WT 153-124" "∆rtt109 153-126"  "WT K36me3" "∆rtt109 K36me3"  --colorMap  'Greens' 'Greens' 'YlOrBr' 'YlOrBr' --hclust 2

#computeMatrix reference-point -p 12 -R "$OUTDIR/Heatmaps/K27me3keptinRtt109.bed" "$OUTDIR/Heatmaps/K27me3lostinrtt109.bed" -S $BWDIR/153-124_ChIP_WT_H3K27me3_.bin_25.smooth_50Bulk.bw $BWDIR/153-126_ChIP_rtt109_H3K27me3_.bin_25.smooth_50Bulk.bw $BWDIR/150-62_WT_H3K36me3_rep8.bin_25.smooth_50.bw $BWDIR/150-70_rtt109_H3K36me3_rep7.bin_25.smooth_50.bw  -o $OUTDIR/Heatmaps/K27me3lostNgainrtt109.gz --outFileNameMatrix $OUTDIR/Heatmaps/K27me3lostNgainrtt109.tab 	--sortRegions keep  --missingDataAsZero -bs 10 -a 2000 -b 1000

#plotHeatmap -m $OUTDIR/Heatmaps/K27me3lostNgainrtt109.gz -o $OUTDIR/Heatmaps/K27me3lostNgainrtt109.png --perGroup  --sortRegions descend --sortUsingSamples 2 --startLabel "5'" --endLabel "3'" --heatmapHeight 8 --heatmapWidth 4  --boxAroundHeatmaps no --samplesLabel "WT 153-124" "∆rtt109 153-126"  "WT K36me3" "∆rtt109 K36me3"  --colorMap  'Greens' 'Greens' 'YlOrBr' 'YlOrBr'

#computeMatrix reference-point -p 12 -R "$GENEDIR/NonK27genes.bed" -S $BWDIR/150-62_WT_H3K36me3_rep8.bin_25.smooth_50.bw $BWDIR/150-70_rtt109_H3K36me3_rep7.bin_25.smooth_50.bw $BWDIR/150-82_rtt109flag_H3K36me3_rep3.bin_25.smooth_50.bw -o $OUTDIR/Heatmaps/RTT109_H3K36me3_nonK27genes_V2.mat 	--sortRegions keep --missingDataAsZero -bs 10 -a 2000 -b 1000

#plotHeatmap -m $OUTDIR/Heatmaps/RTT109_H3K36me3_nonK27genes_V2.mat -o $OUTDIR/Heatmaps/RTT109_H3K36me3_nonK27genes_V2.png --sortRegions descend --sortUsingSamples 1 --outFileSortedRegions $OUTDIR/Heatmaps/RTT109_H3K36me3_nonK27genes_nonK27genes_V2_sortOrder.bed --startLabel "5'" --endLabel "3'" --heatmapHeight 8 --heatmapWidth 4  --boxAroundHeatmaps no --samplesLabel "WT K36me3" "∆rtt109 K36me3" "rtt109-3xflag K36me3" --colorMap 'YlOrBr'


#computeMatrix scale-regions -p 12 -R "$GENEDIR/2024_04_23_WT_peaks.txt" -S $BWDIR/153-124_ChIP_WT_H3K27me3_.bin_25.smooth_50Bulk.bw $BWDIR/153-126_ChIP_rtt109_H3K27me3_.bin_25.smooth_50Bulk.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/150-81_ChIP_rtt109FLAG_H3K27me3__S81_L002_R1_001_val_1.fq.gz.bin_25.smooth_50_Q30.bw --skipZeros -b 250 -a 250 --sortUsingSamples 1 --sortRegions descend -o $OUTDIR/Heatmaps/RTT109_K27domains_2024_04_23_WT_peaks_V1.gz --outFileNameMatrix $OUTDIR/Heatmaps/RTT109_K27domains_2024_04_23_WT_peaks_V1.tab

#plotHeatmap -m $OUTDIR/Heatmaps/RTT109_K27domains_2024_04_23_WT_peaks_V1.gz -o $OUTDIR/Heatmaps/RTT109_K27domains_2024_04_23_WT_peaks_V3.png --sortRegions descend --sortUsingSamples 1 --heatmapHeight 8  --heatmapWidth 4  --outFileSortedRegions $OUTDIR/Heatmaps/RTT109_K27domains_2024_04_23_WT_peaks_V1_sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --samplesLabel  "WT 153-124" "∆rtt109 153-126"  "rtt109-3xflag 150-81" --colorMap 'Greens'  -max 15 15 15



#computeMatrix scale-regions -p 12 -R "$OUTDIR/Heatmaps/RTT109_K27domains_MyceliaK27me3_peaks_V1_sorted_filtered.bed" -S $BWDIR/153-124_ChIP_WT_H3K27me3_.bin_25.smooth_50Bulk.bw $BWDIR/153-126_ChIP_rtt109_H3K27me3_.bin_25.smooth_50Bulk.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/150-81_ChIP_rtt109FLAG_H3K27me3__S81_L002_R1_001_val_1.fq.gz.bin_25.smooth_50_Q30.bw --skipZeros -b 250 -a 250 --sortUsingSamples 1 --sortRegions descend -o $OUTDIR/Heatmaps/RTT109_K27domains_MyceliaK27me3_peaks_V1_sorted_filtered_V1.gz --outFileNameMatrix $OUTDIR/Heatmaps/RTT109_K27domains_MyceliaK27me3_peaks_V1_sorted_filtered_V1.tab

#plotHeatmap -m $OUTDIR/Heatmaps/RTT109_K27domains_MyceliaK27me3_peaks_V1_sorted_filtered_V1.gz -o $OUTDIR/Heatmaps/RTT109_K27domains_MyceliaK27me3_peaks_V1_sorted_filtered_V1.png --sortRegions descend --sortUsingSamples 1 --heatmapHeight 8  --heatmapWidth 4  --outFileSortedRegions $OUTDIR/Heatmaps/RTT109_K27domains_MyceliaK27me3_peaks_V1_sorted_filtered_V2.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --samplesLabel  "WT 153-124" "∆rtt109 153-126"  "rtt109-3xflag 150-81" --colorMap 'Greens'  -max 15 15 15


#computeMatrix scale-regions -p 12 -R "/$GENEDIR/MyceliaK27me3_peaks.bed" -S $BWDIR/153-124_ChIP_WT_H3K27me3_.bin_25.smooth_50Bulk.bw $BWDIR/153-126_ChIP_rtt109_H3K27me3_.bin_25.smooth_50Bulk.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/150-81_ChIP_rtt109FLAG_H3K27me3__S81_L002_R1_001_val_1.fq.gz.bin_25.smooth_50_Q30.bw --skipZeros -b 250 -a 250 --sortUsingSamples 1 --sortRegions descend -o $OUTDIR/Heatmaps/RTT109_K27domains_MyceliaK27me3_peaks_V1.gz --outFileNameMatrix $OUTDIR/Heatmaps/RTT109_K27domains_MyceliaK27me3_peaks_V1.tab


#plotHeatmap -m $OUTDIR/Heatmaps/RTT109_K27domains_MyceliaK27me3_peaks_V1.gz -o $OUTDIR/Heatmaps/RTT109_K27domains_MyceliaK27me3_peaks_V3.png --sortRegions descend --sortUsingSamples 1 --heatmapHeight 8  --heatmapWidth 4  --outFileSortedRegions $OUTDIR/Heatmaps/RTT109_K27domains_MyceliaK27me3_peaks_V1_sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --samplesLabel  "WT 153-124" "∆rtt109 153-126"  "rtt109-3xflag 150-81" --colorMap 'Greens'  -max 15 15 15

computeMatrix scale-regions -p 12 -R "$OUTDIR/Heatmaps/RTT109_K27domains_MyceliaK27me3_peaks_V1_sorted_filtered.bed" -S $BWDIR/138-62_ChIP_WT_H3K36me3_Rep2_6252_S61_L001_R1_001_val_1.fq.gz.bin_25.smooth_50Bulk.bw $BWDIR/137-28_ChIP_WT_H3K36me3_Rep3_S28_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw $BWDIR/137-84_ChIP_rtt109_H3K36me3_Rep3_S79_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw $BWDIR/147-48_ChIP_KOrtt109P13A6P4_rtt109hph_H3K36me3_Rep1_Nc_24hrVMMON_S48_L001_R1_001_val_1.fq.gz.bin_25.smooth_50_Q30.bw $BWDIR/147-56_ChIP_337rtt109P3_rtt109hph_H3K36me3_Rep1_Nc_24hrVMMON_S56_L001_R1_001_val_1.fq.gz.bin_25.smooth_50_Q30.bw $BWDIR/147-60_ChIP_534rtt109P7_rtt109hph_H3K36me3_Rep1_Nc_24hrVMMON_S60_L001_R1_001_val_1.fq.gz.bin_25.smooth_50_Q30.bw $BWDIR/137-86_ChIP_K56R13_H3K36me3_Rep1_S81_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw $BWDIR/137-88_ChIP_K56R40_H3K36me3_Rep1_S83_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw $BWDIR/138-49_ChIP_H3K56R13_H3K36me3_Rep2_6252_S48_L001_R1_001_val_1.fq.gz.bin_25.smooth_50Bulk.bw $BWDIR/138-52_ChIP_H3K56R40_H3K36me3_Rep2_6252_S51_L001_R1_001_val_1.fq.gz.bin_25.smooth_50Bulk.bw $BWDIR/147-44_ChIP_H3K56R13P4_H3K56R13_H3K36me3_Rep1_Nc_24hrVMMON_S44_L001_R1_001_val_1.fq.gz.bin_25.smooth_50_Q30.bw --skipZeros -b 250 -a 250 --sortUsingSamples 1 --sortRegions descend -o $OUTDIR/Heatmaps/H3K36me3_K27domains_MyceliaK27me3_peaks_V1_sorted_filtered_V1.gz --outFileNameMatrix $OUTDIR/Heatmaps/H3K36me3_K27domains_MyceliaK27me3_peaks_V1_sorted_filtered_V1.tab

plotHeatmap -m $OUTDIR/Heatmaps/H3K36me3_K27domains_MyceliaK27me3_peaks_V1_sorted_filtered_V1.gz -o $OUTDIR/Heatmaps/H3K36me3_K27domains_MyceliaK27me3_peaks_V1_sorted_filtered_V1.png --sortRegions descend --sortUsingSamples 1 --heatmapHeight 8  --heatmapWidth 4  --outFileSortedRegions $OUTDIR/Heatmaps/H3K36me3_K27domains_MyceliaK27me3_peaks_V1_sorted_filtered_V1.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --samplesLabel  "138-62_WT_H3K36me3" "137-28_WT_H3K36me3" "137-84_rtt109_H3K36me3" "147-48_ChIP_KOrtt109P13A6P4" "147-56_ChIP_337rtt109P3"  "147-60_ChIP_534rtt109P7" "137-86_ChIP_K56R13_H3K36me3" "137-88_ChIP_K56R40_H3K36me3" "138-49_ChIP_H3K56R13_H3K36me3" "138-52_ChIP_H3K56R40_H3K36me3" "147-44_ChIP_H3K56R13P4_H3K56R13" --colorMap 'YlOrBr'
