#!/bin/bash
#SBATCH --job-name=Heatmaps_MultiBed
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --mem=20gb
#SBATCH --time=08:00:00
#SBATCH --output=../MultipleHeatmaps.%j.out
#SBATCH --error=../MultipleHeatmaps.%j.err

cd $SLURM_SUBMIT_DIR

# Directories
OUTDIR="/scratch/ry00555/RTT109PaperFigures"
GENEDIR="/home/ry00555/Research/Genomes/HeatmapGeneFiles"
BWDIR="$OUTDIR/BigWigs/NormalizedBigWigs/L2FCBigWigs"

# Load deepTools
ml deepTools

#######################################
# K27me3 Heatmaps
#######################################
computeMatrix reference-point -p 12 \
-R "$GENEDIR/K27genes.bed" \
-S $BWDIR/WT_H3K27me3_R2_foldchange.bw \
   $BWDIR/WT_H3K27me3_R3_foldchange.bw \
   $BWDIR/WT_H3K27me3_R4_foldchange.bw \
   $BWDIR/WT_H3K27me3_R5_foldchange.bw \
   $BWDIR/rtt109_H3K27me3_R2_foldchange.bw \
   $BWDIR/rtt109_H3K27me3_R4_foldchange.bw \
   $BWDIR/rtt109_H3K27me3_R5_foldchange.bw \
   $BWDIR/rtt109_H3K27me3_R6_foldchange.bw \
   $BWDIR/rtt109_H3K27me3_R7_foldchange.bw \
   $BWDIR/rtt109_H3K27me3_R8_foldchange.bw \
   $BWDIR/rtt109_H3K27me3_R9_foldchange.bw \
   $BWDIR/rtt109flag_H3K27me3_R1_foldchange.bw \
   $BWDIR/rtt109flag_H3K27me3_R2_foldchange.bw \
   $BWDIR/rtt109flag_H3K27me3_R3_foldchange.bw \
-o $OUTDIR/Heatmaps/RTT109AllK27reps_K27genes.mat \
--sortRegions keep --missingDataAsZero -bs 10 -a 2000 -b 1000

plotHeatmap -m $OUTDIR/Heatmaps/RTT109AllK27reps_K27genes.mat \
-o $OUTDIR/Heatmaps/RTT109AllK27reps_K27genes_V1.png \
--sortRegions descend \
--sortUsingSamples 1 2 3 4 \
--outFileSortedRegions $OUTDIR/Heatmaps/RTT109AllK27reps_K27genes_sortOrder.bed \
--startLabel "5'" --endLabel "3'" \
--heatmapHeight 5 --heatmapWidth 4 \
--boxAroundHeatmaps no \
--samplesLabel "WT K27me3" "WT K27me3" "WT K27me3" "WT K27me3" \
"∆rtt109 K27me3" "∆rtt109 K27me3" "∆rtt109 K27me3" "∆rtt109 K27me3" "∆rtt109 K27me3" \
"rtt109-3xflag K27me3" "rtt109-3xflag K27me3" "rtt109-3xflag K27me3" \
--colorMap 'Greens'


#######################################
# K9me3 Heatmaps
#######################################
computeMatrix reference-point -p 12 \
-R "$GENEDIR/Figure2_K9_Peaks.txt" \
-S $BWDIR/WT_H3K9me3_R1_foldchange.bw \
   $BWDIR/WT_H3K9me3_R2_foldchange.bw \
   $BWDIR/WT_H3K9me3_R3_foldchange.bw \
   $BWDIR/WT_H3K9me3_R4_foldchange.bw \
   $BWDIR/rtt109_H3K9me3_R1_foldchange.bw \
   $BWDIR/rtt109_H3K9me3_R2_foldchange.bw \
   $BWDIR/rtt109_H3K9me3_R3_foldchange.bw \
   $BWDIR/rtt109_H3K9me3_R4_foldchange.bw \
   $BWDIR/rtt109flag_H3K9me3_R1_foldchange.bw \
   $BWDIR/rtt109flag_H3K9me3_R2_foldchange.bw \
-o $OUTDIR/Heatmaps/RTT109All9reps_K9Peaks.mat \
--sortRegions keep --missingDataAsZero -bs 10 -a 2000 -b 1000

plotHeatmap -m $OUTDIR/Heatmaps/RTT109All9reps_K9Peaks.mat \
-o $OUTDIR/Heatmaps/RTT109All9reps_K9Peaks_V1.png \
--sortRegions descend \
--sortUsingSamples 1 2 3 4 \
--outFileSortedRegions $OUTDIR/Heatmaps/RTT109All9reps_K9Peaks_sortOrder.bed \
--startLabel "5'" --endLabel "3'" \
--heatmapHeight 5 --heatmapWidth 4 \
--boxAroundHeatmaps no \
--samplesLabel "WT K9me3" "WT K9me3" "WT K9me3" "WT K9me3" \
"∆rtt109 K9me3" "∆rtt109 K9me3" "∆rtt109 K9me3" "∆rtt109 K9me3" \
"rtt109-3xflag K9me3" "rtt109-3xflag K9me3" \
--colorMap 'Blues'


#######################################
# K36me3 Heatmaps
#######################################
computeMatrix reference-point -p 12 \
-R "$GENEDIR/K27genes.bed" \
-S $BWDIR/WT_H3K36me3_R1_foldchange.bw \
   $BWDIR/WT_H3K36me3_R2_foldchange.bw \
   $BWDIR/WT_H3K36me3_R3_foldchange.bw \
   $BWDIR/WT_H3K36me3_R4_foldchange.bw \
   $BWDIR/rtt109_H3K36me3_R1_foldchange.bw \
   $BWDIR/rtt109_H3K36me3_R2_foldchange.bw \
   $BWDIR/rtt109_H3K36me3_R3_foldchange.bw \
   $BWDIR/rtt109flag_H3K36me3_R1_foldchange.bw \
   $BWDIR/rtt109flag_H3K36me3_R2_foldchange.bw \
   $BWDIR/rtt109flag_H3K36me3_R3_foldchange.bw \
-o $OUTDIR/Heatmaps/RTT109AllK36reps_K27genes.mat \
--sortRegions keep --missingDataAsZero -bs 10 -a 2000 -b 1000

plotHeatmap -m $OUTDIR/Heatmaps/RTT109AllK36reps_K27genes.mat \
-o $OUTDIR/Heatmaps/RTT109AllK36reps_K27genes_V1.png \
--sortRegions descend \
--sortUsingSamples 1 2 3 4 \
--outFileSortedRegions $OUTDIR/Heatmaps/RTT109AllK36reps_K27genes_sortOrder.bed \
--startLabel "5'" --endLabel "3'" \
--heatmapHeight 5 --heatmapWidth 4 \
--boxAroundHeatmaps no \
--samplesLabel "WT K36me3" "WT K36me3" "WT K36me3" "WT K36me3" \
"∆rtt109 K36me3" "∆rtt109 K36me3" "∆rtt109 K36me3" \
"rtt109-3xflag K36me3" "rtt109-3xflag K36me3" "rtt109-3xflag K36me3" \
--colorMap 'YlOrBr'
