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

OUTDIR="/scratch/ry00555/RNASeqPaper/Oct2025"
GENEDIR="/home/ry00555/Research/Genomes/HeatmapGeneFiles"
BWDIR="$OUTDIR/BigWigs"

mkdir -p "$OUTDIR/Heatmaps"

ml deepTools

 computeMatrix reference-point -p 12 -R "$GENEDIR/K27genes.bed" -S $BWDIR/147-3_WT_H3K27me3_rep9.bin_25.smooth_50.bw $BWDIR/147-59_rtt109_H3K27me3_rep9.bin_25.smooth_50.bw $BWDIR/150-65_rtt109_H3K27me3_rep11.bin_25.smooth_50.bw $BWDIR/150-81_rtt109flag_H3K27me3_rep3.bin_25.smooth_50.bw	-o $OUTDIR/Heatmaps/RTT109_H3K27me3_K27genes_V2.mat 	--sortRegions keep --missingDataAsZero -bs 10 -a 1000 -b 1000

 plotHeatmap -m $OUTDIR/Heatmaps/RTT109_H3K27me3_K27genes_V2.mat \
 -o $OUTDIR/Heatmaps/RTT109_H3K27me3_K27genes_V2.png \
 --sortRegions descend \
 --sortUsingSamples 1 4  \
 --outFileSortedRegions $OUTDIR/Heatmaps/RTT109_H3K27me3_K27genes_sortOrder_V2.bed \
 --startLabel "5'" \ --Zmax 20 --endLabel "3'"  --boxAroundHeatmaps no  --samplesLabel "WT K27me3" "∆rtt109 K27me3" "∆rtt109 K27me3" "rtt109-3xflag H3K27me3" --colorMap 'Greens'



 computeMatrix scale-regions -p 12 -R "$GENEDIR/Figure2_K9_Peaks.txt" -S $BWDIR/149-104_WT_H3K9me3_rep5.bin_25.smooth_50.bw  $BWDIR/149-43_rtt109_H3K9me3_rep3.bin_25.smooth_50.bw  $BWDIR/150-80_rtt109flag_H3K9me3_rep2.bin_25.smooth_50.bw --sortRegions keep --missingDataAsZero -bs 10	-b 500 -a 500  	-o $OUTDIR/Heatmaps/RTT109_H3K9me3_K9Peaks_V2.mat

 plotHeatmap -m $OUTDIR/Heatmaps/RTT109_H3K9me3_K9Peaks_V2.mat \
 -o $OUTDIR/Heatmaps/RTT109_H3K9me3_K9Peaks_V2.png \
 --sortRegions descend --sortUsingSamples 1 --heatmapHeight 8 \
 --heatmapWidth 1  \
 --outFileSortedRegions $OUTDIR/Heatmaps/RTT109_H3K9me3_K9Peaks_sortOrder_V2.bed \
 --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no \
 --samplesLabel "WT K9me3" "∆rtt109 K9me3" "rtt109-3xflag K9me3"  --colorMap 'Blues'


 computeMatrix reference-point -p 12 -R "$GENEDIR/K27genes.bed" -S $BWDIR/150-62_WT_H3K36me3_rep8.bin_25.smooth_50.bw $BWDIR/150-70_rtt109_H3K36me3_rep7.bin_25.smooth_50.bw $BWDIR/150-82_rtt109flag_H3K36me3_rep3.bin_25.smooth_50.bw -o $OUTDIR/Heatmaps/RTT109_H3K36me3_K27genes_V2.mat	--sortRegions keep \
 --missingDataAsZero -bs 10 -a 2000 -b 1000

 plotHeatmap -m $OUTDIR/Heatmaps/RTT109_H3K36me3_K27genes_V2.mat -o $OUTDIR/Heatmaps/RTT109_H3K36me3_K27genes_V2.png \
 --sortRegions descend --sortUsingSamples 1 --outFileSortedRegions $OUTDIR/Heatmaps/RTT109_H3K36me3_K27genes_K27genes_V2_sortOrder.bed --startLabel "5'" --endLabel "3'" --heatmapHeight 8 --heatmapWidth 4  --boxAroundHeatmaps no --samplesLabel "WT K36me3" "∆rtt109 K36me3" "rtt109-3xflag K36me3" --colorMap 'YlOrBr'

computeMatrix reference-point -p 12 \
-R "$GENEDIR/NonK27genes.bed" -S $BWDIR/150-62_WT_H3K36me3_rep8.bin_25.smooth_50.bw $BWDIR/150-70_rtt109_H3K36me3_rep7.bin_25.smooth_50.bw $BWDIR/150-82_rtt109flag_H3K36me3_rep3.bin_25.smooth_50.bw -o $OUTDIR/Heatmaps/RTT109_H3K36me3_nonK27genes_V2.mat 	--sortRegions keep --missingDataAsZero -bs 10 -a 2000 -b 1000

plotHeatmap -m $OUTDIR/Heatmaps/RTT109_H3K36me3_nonK27genes_V2.mat \
-o $OUTDIR/Heatmaps/RTT109_H3K36me3_nonK27genes_V2.png \
--sortRegions descend \
--sortUsingSamples 1 \
--outFileSortedRegions $OUTDIR/Heatmaps/RTT109_H3K36me3_nonK27genes_nonK27genes_V2_sortOrder.bed \
--startLabel "5'" \
--endLabel "3'" \
--heatmapHeight 8 \
--heatmapWidth 4  \
--boxAroundHeatmaps no --samplesLabel "WT K36me3" "∆rtt109 K36me3" "rtt109-3xflag K36me3" --colorMap 'YlOrBr'
