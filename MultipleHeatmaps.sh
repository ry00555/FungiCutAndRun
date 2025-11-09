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

 computeMatrix reference-point -p 12 -R "$GENEDIR/K27genes.bed" -S $BWDIR/148-116_WT_H3K27me3_rep11.bin_25.smooth_50.bw $BWDIR/138-57_WT_H3K27me3_rep6.bin_25.smooth_50.bw $BWDIR/147-3_WT_H3K27me3_rep9.bin_25.smooth_50.bw $BWDIR/111-3_rtt109_H3K27me3_rep1.bin_25.smooth_50.bw $BWDIR/135-83_rtt109_H3K27me3_rep2.bin_25.smooth_50.bw $BWDIR/147-59_rtt109_H3K27me3_rep9.bin_25.smooth_50.bw $BWDIR/150-65_rtt109_H3K27me3_rep11.bin_25.smooth_50.bw $BWDIR/142-88_rtt109flag_H3K27me3_rep1.bin_25.smooth_50.bw   $BWDIR/150-81_rtt109flag_H3K27me3_rep3.bin_25.smooth_50.bw $BWDIR/150-77_rtt109flag_H3K27me3_rep2.bin_25.smooth_50.bw	-o $OUTDIR/Heatmaps/RTT109_H3K27me3_K27genes.mat 	--sortRegions keep --missingDataAsZero -bs 10 -a 2000 -b 1000

 plotHeatmap -m $OUTDIR/Heatmaps/RTT109_H3K27me3_K27genes.mat \
 -o $OUTDIR/Heatmaps/RTT109_H3K27me3_K27genes_V1.png \
 --sortRegions descend \
 --sortUsingSamples 1 2 3  \
 --outFileSortedRegions $OUTDIR/Heatmaps/RTT109_H3K27me3_K27genes_sortOrder.bed \
 --startLabel "5'" \
 --endLabel "3'" \
 --heatmapHeight 5  --heatmapWidth 4   --boxAroundHeatmaps no  --samplesLabel "WT K27me3" "WT K27me3" "WT K27me3" "∆rtt109 K27me3" "∆rtt109 K27me3" "∆rtt109 H3K27me3" "∆rtt109 H3K27me3" "rtt109-3xflag H3K27me3" "rtt109-3xflag H3K27me3" "rtt109-3xflag H3K27me3" --colorMap 'Greens' --zMax 20 20 \



 computeMatrix reference-point -p 12 -R "$GENEDIR/Figure2_K9_Peaks.txt" -S  $BWDIR/148-117_WT_H3K9me3_rep4.bin_25.smooth_50.bw $BWDIR/149-104_WT_H3K9me3_rep5.bin_25.smooth_50.bw  $BWDIR/135-25_WT_H3K9me3_rep3.bin_25.smooth_50.bw $BWDIR/149-43_rtt109_H3K9me3_rep3.bin_25.smooth_50.bw  $BWDIR/150-64_rtt109_H3K9me3_rep4.bin_25.smooth_50.bw $BWDIR/150-72_rtt109_H3K9me3_rep6.bin_25.smooth_50.bw  $BWDIR/150-76_rtt109flag_H3K9me3_rep1.bin_25.smooth_50.bw   $BWDIR/150-80_rtt109flag_H3K9me3_rep2.bin_25.smooth_50.bw --sortRegions keep --missingDataAsZero -bs 10 -a 2000 -b 1000 	-o $OUTDIR/Heatmaps/RTT109_H3K9me3_K9Peaks.mat

 plotHeatmap -m $OUTDIR/Heatmaps/RTT109_H3K9me3_K9Peaks.mat \
 -o $OUTDIR/Heatmaps/RTT109_H3K9me3_K9Peaks_V1.png \
 --sortRegions descend \
 --sortUsingSamples 1 2 3 \
 --outFileSortedRegions $OUTDIR/Heatmaps/RTT109_H3K9me3_K9Peaks_sortOrder.bed \
 --startLabel "5'"  --endLabel "3'" --heatmapHeight 5 \
 --heatmapWidth 4   --boxAroundHeatmaps no \
 --samplesLabel "WT K9me3" "WT K9me3" "WT K9me3" "∆rtt109 K9me3" "∆rtt109 K9me3" "∆rtt109 K9me3" "rtt109-3xflag K9me3" "rtt109-3xflag K9me3"  --colorMap 'Blues'


 computeMatrix reference-point -p 12 -R "$GENEDIR/K27genes.bed" -S $BWDIR/138-62_WT_H3K36me3_rep4.bin_25.smooth_50.bw $BWDIR/149-99_WT_H3K36me3_rep7.bin_25.smooth_50.bw $BWDIR/150-62_WT_H3K36me3_rep8.bin_25.smooth_50.bw $BWDIR/150-70_rtt109_H3K36me3_rep7.bin_25.smooth_50.bw $BWDIR/150-74_rtt109_H3K36me3_rep8.bin_25.smooth_50.bw $BWDIR/147-56_rtt109_H3K36me3_rep4.bin_25.smooth_50.bw $BWDIR/137-84_rtt109_H3K36me3_rep2.bin_25.smooth_50.bw $BWDIR/147-48_rtt109_H3K36me3_rep3.bin_25.smooth_50.bw $BWDIR/147-60_rtt109_H3K36me3_rep5.bin_25.smooth_50.bw $BWDIR/150-66_rtt109_H3K36me3_rep6.bin_25.smooth_50.bw $BWDIR/142-89_rtt109flag_H3K36me3_rep1.bin_25.smooth_50.bw  $BWDIR/150-82_rtt109flag_H3K36me3_rep3.bin_25.smooth_50.bw $BWDIR/150-78_rtt109flag_H3K36me3_rep2.bin_25.smooth_50.bw -o $OUTDIR/Heatmaps/RTT109_H3K36me3_K27genes.mat 	--sortRegions keep \
 --missingDataAsZero -bs 10 -a 2000 -b 1000

 plotHeatmap -m $OUTDIR/Heatmaps/RTT109_H3K36me3_K27genes.mat -o $OUTDIR/Heatmaps/RTT109AllK36reps_K27genes_V1.png \
 --sortRegions descend --sortUsingSamples 1 2 3 --outFileSortedRegions $OUTDIR/Heatmaps/RTT109_H3K36me3_K27genes_K27genes_V1_sortOrder.bed --startLabel "5'" --endLabel "3'" --heatmapHeight 5 --heatmapWidth 4 --boxAroundHeatmaps no --samplesLabel "WT K36me3" "WT K36me3" "WT K36me3" "∆rtt109 K36me3" "∆rtt109 K36me3" "∆rtt109 K36me3" "∆rtt109 K36me3" "∆rtt109 K36me3" "∆rtt109 K36me3" "∆rtt109 K36me3" "rtt109-3xflag K36me3" "rtt109-3xflag K36me3" "rtt109-3xflag K36me3" --colorMap 'YlOrBr'

computeMatrix reference-point -p 12 \
-R "$GENEDIR/NonK27genes.bed" -S -S $BWDIR/138-62_WT_H3K36me3_rep4.bin_25.smooth_50.bw $BWDIR/149-99_WT_H3K36me3_rep7.bin_25.smooth_50.bw $BWDIR/150-62_WT_H3K36me3_rep8.bin_25.smooth_50.bw $BWDIR/150-70_rtt109_H3K36me3_rep7.bin_25.smooth_50.bw $BWDIR/150-74_rtt109_H3K36me3_rep8.bin_25.smooth_50.bw $BWDIR/147-56_rtt109_H3K36me3_rep4.bin_25.smooth_50.bw $BWDIR/137-84_rtt109_H3K36me3_rep2.bin_25.smooth_50.bw $BWDIR/147-48_rtt109_H3K36me3_rep3.bin_25.smooth_50.bw $BWDIR/147-60_rtt109_H3K36me3_rep5.bin_25.smooth_50.bw $BWDIR/150-66_rtt109_H3K36me3_rep6.bin_25.smooth_50.bw $BWDIR/142-89_rtt109flag_H3K36me3_rep1.bin_25.smooth_50.bw  $BWDIR/150-82_rtt109flag_H3K36me3_rep3.bin_25.smooth_50.bw $BWDIR/150-78_rtt109flag_H3K36me3_rep2.bin_25.smooth_50.bw -o $OUTDIR/Heatmaps/RTT109_H3K36me3_nonK27genes.mat 	--sortRegions keep --missingDataAsZero -bs 10 -a 2000 -b 1000

plotHeatmap -m $OUTDIR/Heatmaps/RTT109_H3K36me3_nonK27genes.mat \
-o $OUTDIR/Heatmaps/RTT109_H3K36me3_nonK27genes_nonK27genes_V1.png \
--sortRegions descend \
--sortUsingSamples 1 2 3 4 \
--outFileSortedRegions $OUTDIR/Heatmaps/RTT109_H3K36me3_nonK27genes_nonK27genes_V1_sortOrder.bed \
--startLabel "5'" \
--endLabel "3'" \
--heatmapHeight 5 \
--heatmapWidth 4  \
--boxAroundHeatmaps no --samplesLabel "WT K36me3" "WT K36me3" "WT K36me3" "∆rtt109 K36me3" "∆rtt109 K36me3" "∆rtt109 K36me3" "∆rtt109 K36me3" "∆rtt109 K36me3" "∆rtt109 K36me3" "∆rtt109 K36me3" "rtt109-3xflag K36me3" "rtt109-3xflag K36me3" "rtt109-3xflag K36me3" --colorMap 'YlOrBr'
