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

computeMatrix scale-regions -p 12 -R $OUTDIR/MACSPeaksDONE/WT_consensus_merged_V3.bed \
-S $BWDIR/142-10_ChIP_WT_H3K27me3_Rep3.bin_25.smooth_75Bulk.bw \
$BWDIR/145-30_WT_H3K27me3_rep8.bin_25.smooth_50.bw \
$BWDIR/147-3_WT_H3K27me3_rep9.bin_25.smooth_50.bw \
$BWDIR/148-102_WT_H3K27me3_rep10.bin_25.smooth_50.bw \
$BWDIR/148-116_WT_H3K27me3_rep11.bin_25.smooth_50.bw \
$BWDIR/149-55_WT_H3K27me3_rep12.bin_25.smooth_50.bw \
$BWDIR/132-25_ChIP_ncu06787_H3K27me3_Rep_1.bin_25.smooth_75Bulk.bw \
$BWDIR/133-90_ChIP_ncu06787_hph_H3K27me3_Rep1.bin_25.smooth_75Bulk.bw \
$BWDIR/134-11_ChIP_ncu06787_H3K27me3_Rep2.bin_25.smooth_75Bulk.bw \
$BWDIR/135-75_ChIP_ncu06787_H3K27me3_Rep2.bin_25.smooth_75Bulk.bw \
$BWDIR/142-4_ChIP_NCU06787KO_H3K27me3_Rep1_bin_25.smooth_75Bulk.bw \
$BWDIR/147-12_mrg-15_H3K27me3_rep1.bin_25.smooth_50.bw \
$BWDIR/147-15_mrg-15_H3K27me3_rep2.bin_25.smooth_50.bw \
$BWDIR/153-120_ChIP_mrg15_H3K27me3_Rep3_S112_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw \
$BWDIR/134-14_ChIP_ncu06788_H3K27me3_Rep2.bin_25.smooth_75Bulk.bw \
$BWDIR/135-79_cdp-6_H3K27me3_rep4.bin_25.smooth_50.bw \
$BWDIR/142-20_ChIP_NCU06788KOA_H3K27me3_Rep3.bin_25.smooth_75Bulk.bw \
$BWDIR/142-8_ChIP_NCU06788KOa_H3K27me3_Rep3_bin_25.smooth_75Bulk.bw \
$BWDIR/147-18_cdp-6_H3K27me3_rep5.bin_25.smooth_50.bw \
$BWDIR/153-118_ChIP_cdp6_H3K27me3_Rep3_S110_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw --skipZeros -b 500 -a 500 --sortRegions descend -o $OUTDIR/Heatmaps/WT_K27domains_Feb2026_V3.gz --outFileNameMatrix $OUTDIR/Heatmaps/WT_K27domains_Feb2026_V3.tab
#--sortUsingSamples 1

plotHeatmap -m $OUTDIR/Heatmaps/WT_K27domains_Feb2026_V3.gz -o $OUTDIR/Heatmaps/WT_K27domains_Feb2026_V3.png --sortRegions descend --sortUsingSamples 5 --heatmapHeight 5  --heatmapWidth 3  --outFileSortedRegions $OUTDIR/Heatmaps/WT_K27domains_Feb2026_V3_sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --colorMap 'Greens'

computeMatrix reference-point -p 12 -R "$GENEDIR/NonK27genes.bed" \
-S $BWDIR/142-5_ChIP_NCU06787KO_H3K36me3_Rep1.bin_25.smooth_75Bulk.bw \
$BWDIR/153-123_ChIP_mrg15_H3K36me3_Rep2_S115_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw \
$BWDIR/134-15_ChIP_ncu06788_H3K36me3_Rep2.bin_25.smooth_75Bulk.bw \
$BWDIR/142-17_ChIP_NCU06788KOa_H3K36me3_Rep3.bin_25.smooth_75Bulk.bw \
$BWDIR/142-21_ChIP_NCU06788KOA_H3K36me3_Rep3.bin_25.smooth_75Bulk.bw \
$BWDIR/134-29_ChIP_WT_H3K36me3_Rep2.bin_25.smooth_75Bulk.bw \
$BWDIR/134-3_ChIP_WT_H3K36me3_Rep2.bin_25.smooth_75Bulk.bw \
$BWDIR/142-11_ChIP_WT_H3K36me3_Rep3.bin_25.smooth_75Bulk.bw \
$BWDIR/142-77_WT_H3K36me3_rep5.bin_25.smooth_50.bw \
--skipZeros -b 2000 -a 1000 --sortRegions descend -o $OUTDIR/Heatmaps/WT_nonK27genes_H3K36me3_Feb2026_V2.gz --outFileNameMatrix $OUTDIR/Heatmaps/WT_nonK27genes_H3K36me3_Feb2026_V2.tab

plotHeatmap -m $OUTDIR/Heatmaps/WT_nonK27genes_H3K36me3_Feb2026_V2.gz -o $OUTDIR/Heatmaps/WT_nonK27genes_H3K36me3_Feb2026_V3.png --sortRegions descend --sortUsingSamples 1 --heatmapHeight 5  --heatmapWidth 3  --outFileSortedRegions $OUTDIR/Heatmaps/WT_nonK27genes_H3K36me3_Feb2026_V2_sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --colorMap 'YlOrBr' --zMax 8

computeMatrix reference-point -p 12 -R "$GENEDIR/K27genes.bed" \
-S $BWDIR/142-5_ChIP_NCU06787KO_H3K36me3_Rep1.bin_25.smooth_75Bulk.bw \
$BWDIR/153-123_ChIP_mrg15_H3K36me3_Rep2_S115_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw \
$BWDIR/134-15_ChIP_ncu06788_H3K36me3_Rep2.bin_25.smooth_75Bulk.bw \
$BWDIR/142-17_ChIP_NCU06788KOa_H3K36me3_Rep3.bin_25.smooth_75Bulk.bw \
$BWDIR/142-21_ChIP_NCU06788KOA_H3K36me3_Rep3.bin_25.smooth_75Bulk.bw \
$BWDIR/134-29_ChIP_WT_H3K36me3_Rep2.bin_25.smooth_75Bulk.bw \
$BWDIR/134-3_ChIP_WT_H3K36me3_Rep2.bin_25.smooth_75Bulk.bw \
$BWDIR/142-11_ChIP_WT_H3K36me3_Rep3.bin_25.smooth_75Bulk.bw \
$BWDIR/142-77_WT_H3K36me3_rep5.bin_25.smooth_50.bw \
--skipZeros -b 2000 -a 1000 --sortRegions descend -o $OUTDIR/Heatmaps/WT_K27genes_H3K36me3_Feb2026_V2.gz --outFileNameMatrix $OUTDIR/Heatmaps/WT_K27genes_H3K36me3_Feb2026_V2.tab

plotHeatmap -m $OUTDIR/Heatmaps/WT_K27genes_H3K36me3_Feb2026_V2.gz -o $OUTDIR/Heatmaps/WT_K27genes_H3K36me3_Feb2026_V2.png --sortRegions descend --sortUsingSamples 1 --heatmapHeight 5  --heatmapWidth 3  --outFileSortedRegions $OUTDIR/Heatmaps/WT_K27genes_H3K36me3_Feb2026_V2_sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --colorMap 'YlOrBr'

computeMatrix scale-regions -p 12 -R ../MACSPeaksDONE/WT_consensus_merged_V5.bed \
-S ../BigWigs/148-102_WT_H3K27me3_rep10.bin_25.smooth_50.bw \
../BigWigs/147-12_mrg-15_H3K27me3_rep1.bin_25.smooth_50.bw \
../BigWigs/142-4_ChIP_NCU06787KO_H3K27me3_Rep1_bin_25.smooth_75Bulk.bw \
 --skipZeros -b 500 -a 500 --sortRegions descend -o WT_K27domains_Feb2026_V4.gz --outFileNameMatrix WT_K27domains_Feb2026_V4.tab


plotHeatmap -m WT_K27domains_Feb2026_V4.gz -o WT_K27domains_Feb2026_V4.png --sortRegions descend --sortUsingSamples 1 --heatmapHeight 5  --heatmapWidth 3  --outFileSortedRegions WT_K27domains_Feb2026_V4_sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --colorMap 'Greens' --zMax 30

computeMatrix reference-point -p 12 -R "/home/ry00555/Research/Genomes/HeatmapGeneFiles/NonK27genes.bed" \
-S ../BigWigs/142-5_ChIP_NCU06787KO_H3K36me3_Rep1.bin_25.smooth_75Bulk.bw \
../BigWigs/153-123_ChIP_mrg15_H3K36me3_Rep2_S115_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw \
../BigWigs/142-17_ChIP_NCU06788KOa_H3K36me3_Rep3.bin_25.smooth_75Bulk.bw \
../BigWigs/142-21_ChIP_NCU06788KOA_H3K36me3_Rep3.bin_25.smooth_75Bulk.bw \
../BigWigs/142-11_ChIP_WT_H3K36me3_Rep3.bin_25.smooth_75Bulk.bw \
../BigWigs/142-77_WT_H3K36me3_rep5.bin_25.smooth_50.bw \
--skipZeros -b 2000 -a 1000 --sortRegions descend -o WT_nonK27genes_H3K36me3_Feb2026_V2.gz --outFileNameMatrix WT_nonK27genes_H3K36me3_Feb2026_V2.tab

plotHeatmap -m WT_nonK27genes_H3K36me3_Feb2026_V2.gz -o WT_nonK27genes_H3K36me3_Feb2026_V3.png --sortRegions descend --sortUsingSamples 1 --heatmapHeight 5  --heatmapWidth 3  --zMax 8 --outFileSortedRegions WT_nonK27genes_H3K36me3_Feb2026_V2_sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --colorMap 'YlOrBr'

computeMatrix reference-point -p 12 -R "/home/ry00555/Research/Genomes/HeatmapGeneFiles/K27genes.bed" \
-S ../BigWigs/142-5_ChIP_NCU06787KO_H3K36me3_Rep1.bin_25.smooth_75Bulk.bw \
../BigWigs/153-123_ChIP_mrg15_H3K36me3_Rep2_S115_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw \
../BigWigs/142-17_ChIP_NCU06788KOa_H3K36me3_Rep3.bin_25.smooth_75Bulk.bw \
../BigWigs/142-21_ChIP_NCU06788KOA_H3K36me3_Rep3.bin_25.smooth_75Bulk.bw \
../BigWigs/142-11_ChIP_WT_H3K36me3_Rep3.bin_25.smooth_75Bulk.bw \
../BigWigs/142-77_WT_H3K36me3_rep5.bin_25.smooth_50.bw \
--skipZeros -b 2000 -a 1000 --sortRegions descend -o WT_K27genes_H3K36me3_Feb2026_V2.gz --outFileNameMatrix WT_K27genes_H3K36me3_Feb2026_V2.tab

plotHeatmap -m WT_K27genes_H3K36me3_Feb2026_V2.gz -o WT_K27genes_H3K36me3_Feb2026_V2.png --sortRegions descend --sortUsingSamples 1 --heatmapHeight 5  --heatmapWidth 3  --outFileSortedRegions WT_K27genes_H3K36me3_Feb2026_V2_sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --colorMap 'YlOrBr'
