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

computeMatrix scale-regions -p 12 -R $OUTDIR/MACSPeaksDONE/WT_consensus_merged_V2.bed \
-S $BWDIR/133-78_ChIP_WT_H3K27me3_Rep1.bin_25.smooth_75Bulk.bw \
$BWDIR/134-28_ChIP_WT_H3K27me3_Rep2.bin_25.smooth_75Bulk.bw \
$BWDIR/135-27_WT_H3K27me3_rep3.bin_25.smooth_50.bw \
$BWDIR/137-27_ChIP_WT_H3K27me3_Rep3.bin_25.smooth_75Bulk.bw \
$BWDIR/138-57_WT_H3K27me3_rep6.bin_25.smooth_50.bw \
$BWDIR/139-3_ChIP_WT_H3K27me3_.bin_25.smooth_75Bulk.bw \
$BWDIR/142-10_ChIP_WT_H3K27me3_Rep3.bin_25.smooth_75Bulk.bw \
$BWDIR/145-30_WT_H3K27me3_rep8.bin_25.smooth_50.bw \
$BWDIR/147-3_WT_H3K27me3_rep9.bin_25.smooth_50.bw \
$BWDIR/148-102_WT_H3K27me3_rep10.bin_25.smooth_50.bw \
$BWDIR/148-116_WT_H3K27me3_rep11.bin_25.smooth_50.bw \
$BWDIR/149-55_WT_H3K27me3_rep12.bin_25.smooth_50.bw \
$BWDIR/149-63_WT_H3K27me3_rep13.bin_25.smooth_50.bw --skipZeros -b 500 -a 500 --sortRegions descend -o $OUTDIR/Heatmaps/WT_K27domains_Feb2026_V2.gz --outFileNameMatrix $OUTDIR/Heatmaps/WT_K27domains_Feb2026_V2.tab
#--sortUsingSamples 1

plotHeatmap -m $OUTDIR/Heatmaps/WT_K27domains_Feb2026_V2.gz -o $OUTDIR/Heatmaps/WT_K27domains_Feb2026_V2.png --sortRegions descend --sortUsingSamples 1 --heatmapHeight 5  --heatmapWidth 3  --outFileSortedRegions $OUTDIR/Heatmaps/WT_K27domains_Feb2026_V2_sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --colorMap 'Greens'

computeMatrix reference-point -p 12 -R "$GENEDIR/NonK27genes.bed" \
-S $BWDIR/132-27_ChIP_ncu06787_H3K36me3_Rep_1.bin_25.smooth_75Bulk.bw \
$BWDIR/133-91_ChIP_ncu06787_hph_H3K36me3_Rep1.bin_25.smooth_75Bulk.bw \
$BWDIR/134-12_ChIP_ncu06787_H3K36me3_Rep2.bin_25.smooth_75Bulk.bw \
$BWDIR/142-5_ChIP_NCU06787KO_H3K36me3_Rep1.bin_25.smooth_75Bulk.bw \
$BWDIR/153-123_ChIP_mrg15_H3K36me3_Rep2_S115_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw \
$BWDIR/134-15_ChIP_ncu06788_H3K36me3_Rep2.bin_25.smooth_75Bulk.bw \
$BWDIR/142-17_ChIP_NCU06788KOa_H3K36me3_Rep3.bin_25.smooth_75Bulk.bw \
$BWDIR/142-21_ChIP_NCU06788KOA_H3K36me3_Rep3.bin_25.smooth_75Bulk.bw \
$BWDIR/147-19_cdp-6_H3K36me3_rep3.bin_25.smooth_50.bw \
$BWDIR/131-54_ChIP_WT_H3K36me3_Rep1.bin_25.smooth_75Bulk.bw \
$BWDIR/133-79_ChIP_WT_H3K36me3_Rep1.bin_25.smooth_75Bulk.bw \
$BWDIR/134-29_ChIP_WT_H3K36me3_Rep2.bin_25.smooth_75Bulk.bw \
$BWDIR/134-3_ChIP_WT_H3K36me3_Rep2.bin_25.smooth_75Bulk.bw \
$BWDIR/137-28_ChIP_WT_H3K36me3_Rep3.bin_25.smooth_75Bulk.bw \
$BWDIR/138-62_WT_H3K36me3_rep4.bin_25.smooth_50.bw \
$BWDIR/139-4_ChIP_WT_H3K36me3_.bin_25.smooth_75Bulk.bw \
$BWDIR/142-11_ChIP_WT_H3K36me3_Rep3.bin_25.smooth_75Bulk.bw \
$BWDIR/142-77_WT_H3K36me3_rep5.bin_25.smooth_50.bw \
$BWDIR/145-32_WT_H3K36me3_rep6.bin_25.smooth_50.bw \
$BWDIR/149-99_WT_H3K36me3_rep7.bin_25.smooth_50.bw \  --skipZeros -b 2000 -a 1000 --sortRegions descend -o $OUTDIR/Heatmaps/WT_nonK27genes_H3K36me3_Feb2026_V1.gz --outFileNameMatrix $OUTDIR/Heatmaps/WT_nonK27genes_H3K36me3_Feb2026_V1.tab

plotHeatmap -m $OUTDIR/Heatmaps/WT_nonK27genes_H3K36me3_Feb2026_V1.gz -o $OUTDIR/Heatmaps/WT_nonK27genes_H3K36me3_Feb2026_V1.png --sortRegions descend --sortUsingSamples 1 --heatmapHeight 5  --heatmapWidth 3  --outFileSortedRegions $OUTDIR/Heatmaps/WT_nonK27genes_H3K36me3_Feb2026_V1_sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --colorMap 'YlOrBr'

computeMatrix reference-point -p 12 -R "$GENEDIR/K27genes.bed" \
-S $BWDIR/132-27_ChIP_ncu06787_H3K36me3_Rep_1.bin_25.smooth_75Bulk.bw \
$BWDIR/133-91_ChIP_ncu06787_hph_H3K36me3_Rep1.bin_25.smooth_75Bulk.bw \
$BWDIR/134-12_ChIP_ncu06787_H3K36me3_Rep2.bin_25.smooth_75Bulk.bw \
$BWDIR/142-5_ChIP_NCU06787KO_H3K36me3_Rep1.bin_25.smooth_75Bulk.bw \
$BWDIR/153-123_ChIP_mrg15_H3K36me3_Rep2_S115_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw \
$BWDIR/134-15_ChIP_ncu06788_H3K36me3_Rep2.bin_25.smooth_75Bulk.bw \
$BWDIR/142-17_ChIP_NCU06788KOa_H3K36me3_Rep3.bin_25.smooth_75Bulk.bw \
$BWDIR/142-21_ChIP_NCU06788KOA_H3K36me3_Rep3.bin_25.smooth_75Bulk.bw \
$BWDIR/147-19_cdp-6_H3K36me3_rep3.bin_25.smooth_50.bw \
$BWDIR/131-54_ChIP_WT_H3K36me3_Rep1.bin_25.smooth_75Bulk.bw \
$BWDIR/133-79_ChIP_WT_H3K36me3_Rep1.bin_25.smooth_75Bulk.bw \
$BWDIR/134-29_ChIP_WT_H3K36me3_Rep2.bin_25.smooth_75Bulk.bw \
$BWDIR/134-3_ChIP_WT_H3K36me3_Rep2.bin_25.smooth_75Bulk.bw \
$BWDIR/137-28_ChIP_WT_H3K36me3_Rep3.bin_25.smooth_75Bulk.bw \
$BWDIR/138-62_WT_H3K36me3_rep4.bin_25.smooth_50.bw \
$BWDIR/139-4_ChIP_WT_H3K36me3_.bin_25.smooth_75Bulk.bw \
$BWDIR/142-11_ChIP_WT_H3K36me3_Rep3.bin_25.smooth_75Bulk.bw \
$BWDIR/142-77_WT_H3K36me3_rep5.bin_25.smooth_50.bw \
$BWDIR/145-32_WT_H3K36me3_rep6.bin_25.smooth_50.bw \
$BWDIR/149-99_WT_H3K36me3_rep7.bin_25.smooth_50.bw \  --skipZeros -b 2000 -a 1000 --sortRegions descend -o $OUTDIR/Heatmaps/WT_K27genes_H3K36me3_Feb2026_V1.gz --outFileNameMatrix $OUTDIR/Heatmaps/WT_K27genes_H3K36me3_Feb2026_V1.tab

plotHeatmap -m $OUTDIR/Heatmaps/WT_K27genes_H3K36me3_Feb2026_V1.gz -o $OUTDIR/Heatmaps/WT_K27genes_H3K36me3_Feb2026_V1.png --sortRegions descend --sortUsingSamples 1 --heatmapHeight 5  --heatmapWidth 3  --outFileSortedRegions $OUTDIR/Heatmaps/WT_K27genes_H3K36me3_Feb2026_V1_sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --colorMap 'YlOrBr'
