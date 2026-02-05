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

OUTDIR="/lustre2/scratch/ry00555/RNASeqPaper2026"
GENEDIR="/home/ry00555/Research/Genomes/HeatmapGeneFiles"
BWDIR="$OUTDIR/BigWigs"

ml deepTools

computeMatrix scale-regions -p 12 -R $OUTDIR/MACSPeaksDONE/WT_consensus_merged.bed \
-S $BWDIR/133-78_ChIP_WT_H3K27me3_Rep1.bin_25.smooth_75Bulk.bw \
$BWDIR/134-28_ChIP_WT_H3K27me3_Rep2.bin_25.smooth_75Bulk.bw \
$BWDIR/134-2_ChIP_WT_H3K27me3_Rep2.bin_25.smooth_75Bulk.bw \
$BWDIR/135-27_WT_H3K27me3_rep3.bin_25.smooth_50.bw \
$BWDIR/135-9_ChIP_WT_H3K27me3_Rep3.bin_25.smooth_75Bulk.bw \
$BWDIR/137-27_ChIP_WT_H3K27me3_Rep3.bin_25.smooth_75Bulk.bw \
$BWDIR/138-57_WT_H3K27me3_rep6.bin_25.smooth_50.bw \
$BWDIR/139-3_ChIP_WT_H3K27me3_.bin_25.smooth_75Bulk.bw \
$BWDIR/142-10_ChIP_WT_H3K27me3_Rep3.bin_25.smooth_75Bulk.bw \
$BWDIR/145-30_WT_H3K27me3_rep8.bin_25.smooth_50.bw \
$BWDIR/147-3_WT_H3K27me3_rep9.bin_25.smooth_50.bw \
$BWDIR/148-102_WT_H3K27me3_rep10.bin_25.smooth_50.bw \
$BWDIR/148-116_WT_H3K27me3_rep11.bin_25.smooth_50.bw \
$BWDIR/149-55_WT_H3K27me3_rep12.bin_25.smooth_50.bw \
$BWDIR/149-63_WT_H3K27me3_rep13.bin_25.smooth_50.bw --skipZeros -b 500 -a 500 --sortRegions descend -o $OUTDIR/Heatmaps/WT_K27domains_Feb2026_V1.gz --outFileNameMatrix $OUTDIR/Heatmaps/WT_K27domains_Feb2026_V1.tab
#--sortUsingSamples 1

plotHeatmap -m $OUTDIR/Heatmaps/WT_K27domains_Feb2026_V1.gz -o $OUTDIR/Heatmaps/WT_K27domains_Feb2026_V1.png --sortRegions descend --sortUsingSamples 1 --heatmapHeight 5  --heatmapWidth 3  --outFileSortedRegions $OUTDIR/Heatmaps/WT_K27domains_Feb2026_V1_sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --colorMap 'Greens'
