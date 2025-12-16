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

OUTDIR="/scratch/ry00555/RTT109PaperFigures/"
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


#computeMatrix reference-point -p 12 -R "$GENEDIR/K27genes.bed" -S 4 -o $OUTDIR/Heatmaps/RTT109_H3K36me3_K27genes_V3.mat	--sortRegions keep  --missingDataAsZero -bs 10 -a 1000 -b 1000



#plotHeatmap -m $OUTDIR/Heatmaps/RTT109_H3K36me3_K27genes_V3.mat -o $OUTDIR/Heatmaps/RTT109_H3K36me3_K27genes_V3.png --sortRegions descend --sortUsingSamples 1 --outFileSortedRegions $OUTDIR/Heatmaps/RTT109_H3K36me3_K27genes_K27genes_V3_sortOrder.bed --startLabel "5'" --endLabel "3'" --heatmapHeight 8 --heatmapWidth 4  --boxAroundHeatmaps no --samplesLabel "WT K36me3" "∆rtt109 K36me3" "rtt109-3xflag K36me3" --colorMap 'YlOrBr'

#computeMatrix reference-point -p 12 -R "$GENEDIR/NonK27genes.bed" -S $BWDIR/150-62_WT_H3K36me3_rep8.bin_25.smooth_50.bw $BWDIR/150-70_rtt109_H3K36me3_rep7.bin_25.smooth_50.bw $BWDIR/150-82_rtt109flag_H3K36me3_rep3.bin_25.smooth_50.bw -o $OUTDIR/Heatmaps/RTT109_H3K36me3_nonK27genes_V2.mat 	--sortRegions keep --missingDataAsZero -bs 10 -a 2000 -b 1000

#plotHeatmap -m $OUTDIR/Heatmaps/RTT109_H3K36me3_nonK27genes_V2.mat -o $OUTDIR/Heatmaps/RTT109_H3K36me3_nonK27genes_V2.png --sortRegions descend --sortUsingSamples 1 --outFileSortedRegions $OUTDIR/Heatmaps/RTT109_H3K36me3_nonK27genes_nonK27genes_V2_sortOrder.bed --startLabel "5'" --endLabel "3'" --heatmapHeight 8 --heatmapWidth 4  --boxAroundHeatmaps no --samplesLabel "WT K36me3" "∆rtt109 K36me3" "rtt109-3xflag K36me3" --colorMap 'YlOrBr'


computeMatrix reference-point -p 12 -R "$GENEDIR/2024_04_23_WT_peaks.txt" -S $BWDIR/153-124_ChIP_WT_H3K27me3_.bin_25.smooth_50Bulk.bw $BWDIR/153-126_ChIP_rtt109_H3K27me3_.bin_25.smooth_50Bulk.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/150-81_ChIP_rtt109FLAG_H3K27me3__S81_L002_R1_001_val_1.fq.gz.bin_25.smooth_50_Q30.bw --skipZeros -b 1000 -a 1000 --sortUsingSamples 1 --sortRegions descend -o $OUTDIR/Heatmaps/RTT109_K27domains_2024_04_23_WT_peaks_V1.gz --outFileNameMatrix $OUTDIR/Heatmaps/RTT109_K27domains_2024_04_23_WT_peaks_V1.tab

plotHeatmap -m $OUTDIR/Heatmaps/RTT109_K27domains_2024_04_23_WT_peaks_V1.gz -o $OUTDIR/Heatmaps/RTT109_K27domains_2024_04_23_WT_peaks_V1.png --sortRegions descend --sortUsingSamples 1 --heatmapHeight 8  --heatmapWidth 4  --outFileSortedRegions $OUTDIR/Heatmaps/RTT109_K27domains_2024_04_23_WT_peaks_V1_sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --samplesLabel  "WT 153-124" "∆rtt109 153-126"  "rtt109-3xflag 150-81" --colorMap 'Greens'



computeMatrix reference-point -p 12 -R "/$GENEDIR/MyceliaK27me3_peaks.bed" -S $BWDIR/153-124_ChIP_WT_H3K27me3_.bin_25.smooth_50Bulk.bw $BWDIR/153-126_ChIP_rtt109_H3K27me3_.bin_25.smooth_50Bulk.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/150-81_ChIP_rtt109FLAG_H3K27me3__S81_L002_R1_001_val_1.fq.gz.bin_25.smooth_50_Q30.bw --skipZeros -b 1000 -a 1000 --sortUsingSamples 1 --sortRegions descend -o $OUTDIR/Heatmaps/RTT109_K27domains_MyceliaK27me3_peaks_V1.gz --outFileNameMatrix $OUTDIR/Heatmaps/RTT109_K27domains_MyceliaK27me3_peaks_V1.tab


plotHeatmap -m $OUTDIR/Heatmaps/RTT109_K27domains_MyceliaK27me3_peaks_V1.gz -o $OUTDIR/Heatmaps/RTT109_K27domains_MyceliaK27me3_peaks_V1.png --sortRegions descend --sortUsingSamples 1 --heatmapHeight 8  --heatmapWidth 4  --outFileSortedRegions $OUTDIR/Heatmaps/RTT109_K27domains_MyceliaK27me3_peaks_V1_sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --samplesLabel  "WT 153-124" "∆rtt109 153-126"  "rtt109-3xflag 150-81" --colorMap 'Greens'


computeMatrix reference-point -p 12 -R "$GENEDIR/MyceliaK27me3_peaks.txt" -S $BWDIR/153-124_ChIP_WT_H3K27me3_.bin_25.smooth_50Bulk.bw $BWDIR/153-126_ChIP_rtt109_H3K27me3_.bin_25.smooth_50Bulk.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/150-81_ChIP_rtt109FLAG_H3K27me3__S81_L002_R1_001_val_1.fq.gz.bin_25.smooth_50_Q30.bw --skipZeros -b 1000 -a 1000 --sortUsingSamples 1 --sortRegions descend -o $OUTDIR/Heatmaps/RTT109_K27domains_MyceliaK27me3_peakstxt_V1.gz --outFileNameMatrix $OUTDIR/Heatmaps/RTT109_K27domains_MyceliaK27me3_peakstxt_V1.tab

plotHeatmap -m $OUTDIR/Heatmaps/RTT109_K27domains_MyceliaK27me3_peakstxt_V1.gz -o $OUTDIR/Heatmaps/RTT109_K27domains_MyceliaK27me3_peakstxt_V1.png --sortRegions descend --sortUsingSamples 1 --heatmapHeight 8  --heatmapWidth 4  --outFileSortedRegions $OUTDIR/Heatmaps/RTT109_K27domains_MyceliaK27me3_peakstxt_V1sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --samplesLabel  "WT 153-124" "∆rtt109 153-126"  "rtt109-3xflag 150-81" --colorMap 'Greens'
