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

computeMatrix reference-point -p 12 -R "../Heatmaps/eaf3_EctopicGenes_V1.bed" -S "H3K36me3_WT_ChIP.bw" "H3K36me3_eaf3_ChIP.bw"  "H3K36me3_cdp6_ChIP.bw" "H3K36me3_mrg15_ChIP.bw"  --skipZeros -b 2000 -a 1000 -o  ../Heatmaps/eaf3_EctopicGenes_V1_H3K36me3.gz --outFileNameMatrix ../Heatmaps/eaf3_EctopicGenes_V1_H3K36me3.tab

plotHeatmap -m ../Heatmaps/eaf3_EctopicGenes_V1_H3K36me3.gz -o ../Heatmaps/eaf3_EctopicGenes_V3_H3K36me3.png \
--heatmapHeight 5 \
--heatmapWidth 3 \
    --startLabel "5'" \
    --endLabel "3'" \
    --boxAroundHeatmaps no \
    --colorMap 'YlOrBr' --zMax 30  \
      --samplesLabel "WT H3K36me3" "eaf3" "cdp6" "mrg15"




computeMatrix reference-point -p 12 -R "/home/ry00555/Research/Genomes/HeatmapGeneFiles/EctopicH3K27me3GenesGained_Final_wNCU.bed" -S "H3K27me3_WT_ChIP.bw" "H3K27me3_eaf3_ChIP.bw"  "H3K27me3_cdp6_ChIP.bw" "H3K27me3_mrg15_ChIP.bw"  --skipZeros -b 1000 -a 2000 --sortRegions descend -o ../Heatmaps/eaf3_EctopicGenes_V1.gz --outFileNameMatrix ../Heatmaps/eaf3_EctopicGenes_V1.tab

plotHeatmap -m ../Heatmaps/eaf3_EctopicGenes_V1.gz -o ../Heatmaps/eaf3_EctopicGenes_V2.png   \
    --heatmapHeight 5 \
    --heatmapWidth 3 \
    --outFileSortedRegions ../Heatmaps/eaf3_EctopicGenes_V1_sorted.bed \
    --startLabel "5'" \
    --endLabel "3'" \
    --boxAroundHeatmaps no \
    --colorMap 'Greens' 'Greens' 'Greens' 'Greens'  --zMax 40 40 40 40  \
      --samplesLabel "WT H3K27me3" "eaf3" "cdp6" "mrg15"

      computeMatrix scale-regions -p 12 -R "/home/ry00555/Research/Genomes/HeatmapGeneFiles/EctopicH3K27me3Domains.bed" -S "H3K27me3_WT_ChIP.bw" "H3K27me3_eaf3_ChIP.bw"  "H3K27me3_cdp6_ChIP.bw" "H3K27me3_mrg15_ChIP.bw"  --skipZeros -b 500 -a 500 --sortRegions descend -o ../Heatmaps/Ectopic_eaf3_K27domains_V1.gz --outFileNameMatrix ../Heatmaps/Ectopic_eaf3_K27domains_V1.tab


      plotHeatmap -m ../Heatmaps/Ectopic_eaf3_K27domains_V1.gz -o ../Heatmaps/Ectopic_eaf3_K27domains_V1.png --sortRegions descend --sortUsingSamples 4 --heatmapHeight 1  --heatmapWidth 3  --outFileSortedRegions ../Heatmaps/Ectopic_eaf3_K27domains_V1.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --colorMap 'Greens' --zMax 6 --samplesLabel "WT H3K27me3" "eaf3" "cdp6" "mrg15"



computeMatrix scale-regions -p 12 -R "../Heatmaps/WT_H3K27me3_domains_V3_eaf3.bed" -S "H3K27me3_WT_ChIP.bw" "H3K27me3_eaf3_ChIP.bw"  "H3K27me3_cdp6_ChIP.bw" "H3K27me3_mrg15_ChIP.bw"  --skipZeros -b 500 -a 500 --sortRegions descend -o ../Heatmaps/eaf3_K27domains_V3.gz --outFileNameMatrix ../Heatmaps/eaf3_K27domains_V3.tab


plotHeatmap -m ../Heatmaps/eaf3_K27domains_V3.gz -o ../Heatmaps/eaf3_K27domains_V4.png --sortRegions descend --sortUsingSamples 2 --heatmapHeight 5  --heatmapWidth 3  --outFileSortedRegions ../Heatmaps/WT_H3K27me3_domains_V3_eaf3_sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --colorMap 'Greens' --zMax 8 --samplesLabel "WT H3K27me3" "eaf3" "cdp6" "mrg15"


computeMatrix reference-point -p 12 -R "/home/ry00555/Research/Genomes/HeatmapGeneFiles/neurospora.bed" -S "H3K27me3_WT_ChIP.bw" "H3K27me3_eaf3_ChIP.bw"  "H3K27me3_cdp6_ChIP.bw" "H3K27me3_mrg15_ChIP.bw"  --skipZeros -b 2000 -a 1000 --sortRegions descend -o ../Heatmaps/eaf3_AllGenes_V1.gz --outFileNameMatrix ../Heatmaps/eaf3_AllGenes_V1.tab

plotHeatmap -m ../Heatmaps/eaf3_AllGenes_V1.gz -o ../Heatmaps/eaf3_AllGenes_V1.png   --sortRegions descend --sortUsingSamples 1 \
    --heatmapHeight 5 \
    --heatmapWidth 3 \
    --outFileSortedRegions ../Heatmaps/eaf3_AllGenes_V1.bed \
    --startLabel "5'" \
    --endLabel "3'" \
    --boxAroundHeatmaps no \
    --colorMap 'Greens' 'Greens' 'Greens' 'Greens'  --zMax 20 20 20 20  \
      --samplesLabel "WT H3K27me3" "eaf3" "cdp6" "mrg15"

      computeMatrix reference-point -p 12 -R "../Heatmaps/eaf3_H3K27me3_WTgenesonly_V4.bed" -S "H3K27me3_WT_ChIP.bw" "H3K27me3_eaf3_ChIP.bw"  "H3K27me3_cdp6_ChIP.bw" "H3K27me3_mrg15_ChIP.bw"  --skipZeros -b 1000 -a 2000 --sortRegions descend -o ../Heatmaps/eaf3_H3K27me3genes_V8.gz --outFileNameMatrix ../Heatmaps/eaf3_H3K27me3genes_V8.tab

      plotHeatmap -m ../Heatmaps/eaf3_H3K27me3genes_V8.gz \
          -o ../Heatmaps/eaf3_H3K27me3genes_V11.png \
          --sortRegions descend \
          --sortUsingSamples 3 \
          --heatmapHeight 5 \
          --heatmapWidth 3 \
          --outFileSortedRegions ../Heatmaps/eaf3_H3K27me3_WTgenesonly_V4_cdp6sorted.bed \
          --startLabel "5'" \
          --endLabel "3'" \
          --boxAroundHeatmaps no \
          --colorMap 'Greens' 'Greens' 'Greens' 'Greens'  --zMax 40 40 40 40  \
            --samplesLabel "WT H3K27me3" "eaf3" "cdp6" "mrg15"

            plotHeatmap -m ../Heatmaps/eaf3_H3K27me3genes_V4.gz \
                -o ../Heatmaps/eaf3_H3K27me3genes_V5.png \
                --sortRegions descend \
                --sortUsingSamples 1 \
                --heatmapHeight 5 \
                --heatmapWidth 3 \
                --outFileSortedRegions ../Heatmaps/eaf3_H3K27me3_top729list_WTsorted.bed \
                --startLabel "5'" \
                --endLabel "3'" \
                --boxAroundHeatmaps no \
                --colorMap 'Greens' 'Greens' 'Greens' 'Greens'  --zMax 40 40 40 40  \
                  --samplesLabel "WT H3K27me3" "eaf3" "cdp6" "mrg15"

REPLICATES=("134-28_ChIP_WT_H3K27me3_Rep2_peaks.xls"
"147-3_WT_H3K27me3_rep9_peaks.xls"
"145-30_ChIP_peaks.xls")
OUT_DIR="BedToolsMergedPeaks"
INTERMEDIATE_BEDS=()
for file in "${REPLICATES[@]}"; do
    base_name="${file%.xls}"
    output_path="$OUT_DIR/${base_name}_WT_.bed"
    awk 'BEGIN{OFS="\t"}
         /^[^#]/ && $1!="" && $1!="chr" {
             print $1, $2-1, $3, $4
         }' "$file" > "$output_path"
        INTERMEDIATE_BEDS+=("$output_path")
    echo "Converted $file -> $output_path"
done
MERGED_FILE="$OUT_DIR/merged_peaks_master.bed"

echo "Merging overlapping peaks across all replicates..."
cat "${INTERMEDIATE_BEDS[@]}" | \
    bedtools sort -i stdin | \
    bedtools merge -i stdin > "$MERGED_FILE"


computeMatrix reference-point -p 12 -R "/home/ry00555/Research/Genomes/HeatmapGeneFiles/K27genes_NCU_n573.bed" -S "/scratch/ry00555/Run153/BigWigs/153-26_ChIP_WT_H3K27me3_Rep1_S26_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw" "/scratch/ry00555/Run153/BigWigs/153-107_ChIP_WT_H3K36me3__S100_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw" "/scratch/ry00555/Run153/BigWigs/153-38_ChIP_set7_H3K27me3__S38_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw" "/scratch/ry00555/Run153/BigWigs/153-40_ChIP_set7_H3K36me3__S40_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw" --skipZeros -b 2000 -a 1000 --sortRegions descend -o ASH1_SET7_K27me3_H3K36me3_H3K27me3genes_V1.gz --outFileNameMatrix ASH1_SET7_K27me3_H3K36me3_H3K27me3genes_V1.tab

plotHeatmap -m ASH1_SET7_K27me3_H3K36me3_H3K27me3genes_V1.gz \
    -o ASH1_SET7_K27me3_H3K36me3_H3K27me3genes_V1.png \
    --sortRegions descend \
    --sortUsingSamples 4 \
    --heatmapHeight 5 \
    --heatmapWidth 3 \
    --outFileSortedRegions ASH1_SET7_K27me3_H3K36me3_H3K27me3genes_V1.bed \
    --startLabel "5'" \
    --endLabel "3'" \
    --boxAroundHeatmaps no \
    --colorMap 'Greens' 'YlOrBr' 'Greens' 'YlOrBr'  --zMax 20 10 20 10  \
      --samplesLabel "WT H3K27me3" "WT H3K36me3" "set7 H3K27me3" "set7 H3K36me3"

computeMatrix reference-point -p 12 -R "/home/ry00555/Research/Genomes/HeatmapGeneFiles/K27genes_NCU_n573.bed" -S "/scratch/ry00555/Run153/BigWigs/153-26_ChIP_WT_H3K27me3_Rep1_S26_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw" "/scratch/ry00555/Run153/BigWigs/153-107_ChIP_WT_H3K36me3__S100_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw" "/scratch/ry00555/Run153/BigWigs/153-38_ChIP_set7_H3K27me3__S38_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw" "/scratch/ry00555/Run153/BigWigs/153-40_ChIP_set7_H3K36me3__S40_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw" "/scratch/ry00555/Run153/BigWigs/153-30_ChIP_set2_H3K27me3_Rep1_S30_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw"  "/scratch/ry00555/Run153/BigWigs/153-34_ChIP_set2_H3K27me3_Rep1_S34_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw"  "/scratch/ry00555/Run153/BigWigs/153-35_ChIP_set2_H3K27me3_Rep2_S35_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw" "/scratch/ry00555/Run153/BigWigs/153-90_ChIP_set2set7_H3K36me3_Rep1_S84_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw" "/lustre2/scratch/ry00555/RNASeqPaper2026/BigWigs/142-119_set-2_H3K36me3_rep1.bin_25.smooth_50.bw" "ash1_Y888F_H3K27me2_3_renamed.bigwig" "ash1_Y888F_H3K27me2_3_scaled.bigwig" "cat_ash1_H3K36me3_renamed.bigwig" "cat_ash1_H3K36me3_scaled.bigwig" --skipZeros -b 2000 -a 1000 --sortRegions descend -o ASH1_SET2_K27me3_H3K36me3_H3K27me3genes_V3.gz --outFileNameMatrix ASH1_SET2_K27me3_H3K36me3_H3K27me3genes_V3.tab

plotHeatmap -m ASH1_SET2_K27me3_H3K36me3_H3K27me3genes_V3.gz \
    -o ASH1_SET2_K27me3_H3K36me3_H3K27me3genes_V4.png \
    --sortRegions descend \
    --sortUsingSamples 1 \
    --heatmapHeight 5 \
    --heatmapWidth 3 \
    --outFileSortedRegions ASH1_SET2_K27me3_H3K36me3_H3K27me3genes_V4.bed \
    --startLabel "5'" \
    --endLabel "3'" \
    --boxAroundHeatmaps no \
    --colorMap 'Greens' 'YlOrBr' 'Greens' 'YlOrBr' 'Greens' 'Greens' 'Greens' 'YlOrBr' 'YlOrBr' 'Greens' 'Greens' 'YlOrBr' 'YlOrBr' \
    --zMax 20 10 20 10 20 20 20 10 10 20 20 10 10 \
    --samplesLabel "WT H3K27me3" "WT H3K36me3" "set7 H3K27me3" "set7 H3K36me3" "set2 H3K27me3 Rep1" "set2 H3K27me3 Rep1b" "set2 H3K27me3 Rep2" "set2set7 H3K36me3" "set2 H3K36me3" "Y888F H3K27me2/3" "Y888F H3K27me2/3 scaled" "ASH1cat H3K36me3" "ASH1cat H3K36me3 scaled"

    computeMatrix reference-point -p 12 \
        -R "/home/ry00555/Research/Genomes/HeatmapGeneFiles/K27genes_NCU_n573.bed" \
        -S "/scratch/ry00555/Run153/BigWigs/153-26_ChIP_WT_H3K27me3_Rep1_S26_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw" \
        "/scratch/ry00555/Run153/BigWigs/153-107_ChIP_WT_H3K36me3__S100_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw" \
        "/scratch/ry00555/Run153/BigWigs/153-38_ChIP_set7_H3K27me3__S38_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw" \
        "/scratch/ry00555/Run153/BigWigs/153-40_ChIP_set7_H3K36me3__S40_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw" \
        "/scratch/ry00555/Run153/BigWigs/153-35_ChIP_set2_H3K27me3_Rep2_S35_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw" \
        "/scratch/ry00555/Run153/BigWigs/153-90_ChIP_set2set7_H3K36me3_Rep1_S84_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw" \
        "../BigWigs/ash1_Y888F_H3K27me2_3_scaled.bigwig" \
        "../BigWigs/cat_ash1_H3K36me3_scaled.bigwig" \
        --skipZeros -b 2000 -a 2000 \
        --sortRegions descend \
        -o ASH1_SET2_K27me3_H3K36me3_H3K27me3genes_V7.gz \
        --outFileNameMatrix ASH1_SET2_K27me3_H3K36me3_H3K27me3genes_V7.tab

    plotHeatmap -m ASH1_SET2_K27me3_H3K36me3_H3K27me3genes_V7.gz \
        -o ASH1_SET2_K27me3_H3K36me3_H3K27me3genes_V7.png \
        --sortRegions descend \
        --sortUsingSamples 1 2 \
        --heatmapHeight 8 \
        --heatmapWidth 3 \
        --outFileSortedRegions ASH1_SET2_K27me3_H3K36me3_H3K27me3genes_V7.bed \
        --startLabel "5'" \
        --endLabel "3'" \
        --boxAroundHeatmaps no \
        --colorMap 'Greens' 'YlOrBr' 'Greens' 'YlOrBr' 'Greens' 'YlOrBr' 'Greens' 'YlOrBr' \
        --zMax 20 10 20 8 20 10 20 10 \
        --samplesLabel "WT H3K27me3" "WT H3K36me3" "set7 H3K27me3" "set7 H3K36me3" "set2 H3K27me3" "set2set7 H3K36me3" "Y888F H3K27me2/3" "ASH1cat H3K36me3"

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

computeMatrix reference-point -p 12 -R "/home/ry00555/Research/Genomes/HeatmapGeneFiles/K27genes_NCU_n573.bed" -S ../BigWigs/153-15_ChIP_WT_H3K27me3_Rep1.bin_25.smooth_50Bulk.bw ../BigWigs/153-3_ChIP_WT_H3K27me3_Rep1_S3_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw ../BigWigs/155-76_ChIP_C8_H3K27me3__S76.bin_25.smooth_75Bulk.bw ../BigWigs/155-79_ChIP_C9_H3K27me3__S79.bin_25.smooth_75Bulk.bw ../BigWigs/155-82_ChIP_C10a_H3K27me3__S82.bin_25.smooth_75Bulk.bw ../BigWigs/155-88_ChIP_eaf3bar_H3K27me3__S88.bin_25.smooth_75Bulk.bw  ../BigWigs/155-91_ChIP_eaf3bar_H3K27me3__S91.bin_25.smooth_75Bulk.bw --skipZeros -b 2000 -a 1000 --sortRegions descend -o EAF3_K27genes_H3K27me3_May2026_V1.gz --outFileNameMatrix EAF3_K27genes_H3K27me3_May2026_V1.tab

plotHeatmap -m EAF3_K27genes_H3K27me3_May2026_V1.gz -o EAF3_K27genes_H3K27me3_May2026_V1.png --sortRegions descend --sortUsingSamples 7 --heatmapHeight 5  --heatmapWidth 3  --outFileSortedRegions EAF3_K27genes_H3K27me3_May2026_V1_sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --colorMap 'Greens' --zMax 30


computeMatrix reference-point -p 12 -R "/home/ry00555/Research/Genomes/HeatmapGeneFiles/K27genes_NCU_n573.bed" -S ../BigWigs/153-15_ChIP_WT_H3K27me3_Rep1.bin_25.smooth_50Bulk.bw  ../BigWigs/155-82_ChIP_C10a_H3K27me3__S82.bin_25.smooth_75Bulk.bw ../BigWigs/155-76_ChIP_C8_H3K27me3__S76.bin_25.smooth_75Bulk.bw ../BigWigs/155-91_ChIP_eaf3bar_H3K27me3__S91.bin_25.smooth_75Bulk.bw  --skipZeros -b 2000 -a 2000 --sortRegions descend -o EAF3_K27genes_H3K27me3_May2026_V3.gz --outFileNameMatrix EAF3_K27genes_H3K27me3_May2026_V3.tab

plotHeatmap -m EAF3_K27genes_H3K27me3_May2026_V3.gz -o EAF3_K27genes_H3K27me3_May2026_V6.png --sortRegions descend --sortUsingSamples 4 --heatmapHeight 5  --heatmapWidth 3  --outFileSortedRegions EAF3_K27genes_H3K27me3_May2026_V3.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no --samplesLabel "WT" "cdp6" "mrg15"  "eaf3" --colorMap 'Greens' --zMax 30

computeMatrix scale-regions -p 12 -R "/home/ry00555/Research/Genomes/HeatmapGeneFiles/WT_H3K27me3_consensus_MACSpeaks.bed" -S ../BigWigs/153-15_ChIP_WT_H3K27me3_Rep1.bin_25.smooth_50Bulk.bw  ../BigWigs/155-82_ChIP_C10a_H3K27me3__S82.bin_25.smooth_75Bulk.bw ../BigWigs/155-76_ChIP_C8_H3K27me3__S76.bin_25.smooth_75Bulk.bw ../BigWigs/155-91_ChIP_eaf3bar_H3K27me3__S91.bin_25.smooth_75Bulk.bw  --skipZeros -b 500 -a 500 --sortRegions descend -o EAF3_K27domains_H3K27me3_May2026_V2.gz --outFileNameMatrix EAF3_K27domains_H3K27me3_May2026_V2.tab

plotHeatmap -m EAF3_K27domains_H3K27me3_May2026_V2.gz -o EAF3_K27domains_H3K27me3_May2026_V2.png --sortRegions descend --sortUsingSamples 4 --heatmapHeight 5  --heatmapWidth 3  --outFileSortedRegions vEAF3_K27domains_H3K27me3_May2026_V2_sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no --samplesLabel "WT" "cdp6" "mrg15"  "eaf3" --colorMap 'Greens' --zMax 15

computeMatrix reference-point -p 12 -R "/home/ry00555/Research/Genomes/HeatmapGeneFiles/K27genes_NCU_n573.bed" -S ../BigWigs/153-107_ChIP_WT_H3K36me3_.bin_25.smooth_50Bulk.bw ../BigWigs/155-77_ChIP_C8_H3K36me3__S77.bin_25.smooth_75Bulk.bw ../BigWigs/155-80_ChIP_C9_H3K36me3__S80.bin_25.smooth_75Bulk.bw  ../BigWigs/155-83_ChIP_C10a_H3K36me3__S83.bin_25.smooth_75Bulk.bw ../BigWigs/155-86_ChIP_C10b_H3K36me3__S86.bin_25.smooth_75Bulk.bw ../BigWigs/155-89_ChIP_eaf3bar_H3K36me3__S89.bin_25.smooth_75Bulk.bw ../BigWigs/155-92_ChIP_eaf3bar_H3K36me3__S92.bin_25.smooth_75Bulk.bw --skipZeros -b 2000 -a 1000 --sortRegions descend -o EAF3_K27genes_H3K36me3_May2026_V1.gz --outFileNameMatrix EAF3_K27genes_H3K36me3_May2026_V1.tab

plotHeatmap -m EAF3_K27genes_H3K36me3_May2026_V1.gz -o EAF3_K27genes_H3K36me3_May2026_V1.png --sortRegions descend --sortUsingSamples 7 --heatmapHeight 5  --heatmapWidth 3  --outFileSortedRegions EAF3_K27genes_H3K36me3_May2026_V1_sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --colorMap 'YlOrBr' --zMax 10

computeMatrix reference-point -p 12 -R "/home/ry00555/Research/Genomes/HeatmapGeneFiles/NonK27genes.bed" -S ../BigWigs/153-107_ChIP_WT_H3K36me3_.bin_25.smooth_50Bulk.bw ../BigWigs/155-77_ChIP_C8_H3K36me3__S77.bin_25.smooth_75Bulk.bw ../BigWigs/155-80_ChIP_C9_H3K36me3__S80.bin_25.smooth_75Bulk.bw  ../BigWigs/155-83_ChIP_C10a_H3K36me3__S83.bin_25.smooth_75Bulk.bw ../BigWigs/155-86_ChIP_C10b_H3K36me3__S86.bin_25.smooth_75Bulk.bw ../BigWigs/155-89_ChIP_eaf3bar_H3K36me3__S89.bin_25.smooth_75Bulk.bw ../BigWigs/155-92_ChIP_eaf3bar_H3K36me3__S92.bin_25.smooth_75Bulk.bw --skipZeros -b 2000 -a 1000 --sortRegions descend -o EAF3_nonK27genes_H3K36me3_May2026_V1.gz --outFileNameMatrix EAF3_nonK27genes_H3K36me3_May2026_V1.tab

plotHeatmap -m EAF3_nonK27genes_H3K36me3_May2026_V1.gz -o EAF3_nonK27genes_H3K36me3_May2026_V1.png --sortRegions descend --sortUsingSamples 1 --heatmapHeight 5  --heatmapWidth 3  --outFileSortedRegions EAF3_nonK27genes_H3K36me3_May2026_V1_sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --colorMap 'YlOrBr' --zMax 10

computeMatrix reference-point -p 12 -R "/home/ry00555/Research/Genomes/HeatmapGeneFiles/K27genes_NCU_n573.bed" -S ../BigWigs/153-107_ChIP_WT_H3K36me3_.bin_25.smooth_50Bulk.bw ../BigWigs/155-80_ChIP_C9_H3K36me3__S80.bin_25.smooth_75Bulk.bw  ../BigWigs/155-77_ChIP_C8_H3K36me3__S77.bin_25.smooth_75Bulk.bw  ../BigWigs/155-92_ChIP_eaf3bar_H3K36me3__S92.bin_25.smooth_75Bulk.bw  --skipZeros -b 2000 -a 1000 --sortRegions descend -o EAF3_K27genes_H3K36me3_May2026_V2.gz --outFileNameMatrix EAF3_K27genes_H3K36me3_May2026_V2.tab

plotHeatmap -m EAF3_K27genes_H3K36me3_May2026_V2.gz -o EAF3_K27genes_H3K36me3_May2026_V2.png --sortRegions descend --sortUsingSamples 1 --heatmapHeight 5  --heatmapWidth 3  --outFileSortedRegions EAF3_K27genes_H3K36me3_May2026_V2_sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --colorMap 'YlOrBr' --zMax 8 --samplesLabel "WT" "cdp6" "mrg15" "eaf3"

computeMatrix reference-point -p 12 -R "/home/ry00555/Research/Genomes/HeatmapGeneFiles/nonASH1_noK27_wCM.bed" -S ../BigWigs/153-107_ChIP_WT_H3K36me3_.bin_25.smooth_50Bulk.bw  ../BigWigs/155-80_ChIP_C9_H3K36me3__S80.bin_25.smooth_75Bulk.bw ../BigWigs/155-77_ChIP_C8_H3K36me3__S77.bin_25.smooth_75Bulk.bw ../BigWigs/155-92_ChIP_eaf3bar_H3K36me3__S92.bin_25.smooth_75Bulk.bw --skipZeros -b 2000 -a 1000 --sortRegions descend -o EAF3_nonASH1_noK27_H3K36me3_May2026_V1.gz --outFileNameMatrix EAF3_nonASH1_noK27_H3K36me3_May2026_V1.tab

plotHeatmap -m EAF3_nonASH1_noK27_H3K36me3_May2026_V1.gz -o EAF3_nonASH1_noK27_H3K36me3_May2026_V2.png --sortRegions descend --sortUsingSamples 1 --heatmapHeight 8  --heatmapWidth 3  --outFileSortedRegions EAF3_nonASH1_noK27_H3K36me3_May2026_V1_sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --colorMap 'YlOrBr' --zMax 8 --samplesLabel "WT"  "cdp6" "mrg15" "eaf3"

computeMatrix reference-point -p 12 -R "/home/ry00555/Research/Genomes/HeatmapGeneFiles/Bicocca_ASH1_noK27_wCM.bed" -S ../BigWigs/153-107_ChIP_WT_H3K36me3_.bin_25.smooth_50Bulk.bw  ../BigWigs/155-80_ChIP_C9_H3K36me3__S80.bin_25.smooth_75Bulk.bw ../BigWigs/155-77_ChIP_C8_H3K36me3__S77.bin_25.smooth_75Bulk.bw ../BigWigs/155-92_ChIP_eaf3bar_H3K36me3__S92.bin_25.smooth_75Bulk.bw --skipZeros -b 2000 -a 1000 --sortRegions descend -o EAF3_ASH1_nonK27genes_H3K36me3_May2026_V1.gz --outFileNameMatrix EAF3_ASH1_nonK27genes_H3K36me3_May2026_V1.tab

plotHeatmap -m EAF3_ASH1_nonK27genes_H3K36me3_May2026_V1.gz -o EAF3_ASH1_nonK27genes_H3K36me3_June2026_V2.png --sortRegions descend --sortUsingSamples 2 --heatmapHeight 5  --heatmapWidth 3  --outFileSortedRegions EAF3_ASH1_nonK27genes_H3K36me3_May2026_V1_sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --colorMap 'YlOrBr' --zMax 8 --samplesLabel "WT"  "cdp6" "mrg15" "eaf3"

computeMatrix reference-point -p 12 -R "/home/ry00555/Research/Genomes/HeatmapGeneFiles/Bicocca_ASH1_noK27_wCM.bed" -S ../BigWigs/153-107_ChIP_WT_H3K36me3_.bin_25.smooth_50Bulk.bw  ../BigWigs/155-80_ChIP_C9_H3K36me3__S80.bin_25.smooth_75Bulk.bw ../BigWigs/155-77_ChIP_C8_H3K36me3__S77.bin_25.smooth_75Bulk.bw ../BigWigs/155-92_ChIP_eaf3bar_H3K36me3__S92.bin_25.smooth_75Bulk.bw --skipZeros -b 2000 -a 2000 --sortRegions descend -o EAF3_ASH1_nonK27genes_H3K36me3_June2026_V2.gz --outFileNameMatrix EAF3_ASH1_nonK27genes_H3K36me3_June2026_V2.tab

plotHeatmap -m EAF3_ASH1_nonK27genes_H3K36me3_June2026_V2.gz -o EAF3_ASH1_nonK27genes_H3K36me3_June2026_V2.png --sortRegions descend --sortUsingSamples 2 --heatmapHeight 5  --heatmapWidth 3  --outFileSortedRegions EAF3_ASH1_nonK27genes_H3K36me3_June2026_V2_sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --colorMap 'YlOrBr' --zMax 8 --samplesLabel "WT"  "cdp6" "mrg15" "eaf3"

computeMatrix reference-point -p 12 -R "/home/ry00555/Research/Genomes/HeatmapGeneFiles/Bicocca_ASH1_noK27_wCM.bed" -S ../BigWigs/153-107_ChIP_WT_H3K36me3_.bin_25.smooth_50Bulk.bw  ../BigWigs/155-80_ChIP_C9_H3K36me3__S80.bin_25.smooth_75Bulk.bw ../BigWigs/155-77_ChIP_C8_H3K36me3__S77.bin_25.smooth_75Bulk.bw ../BigWigs/155-92_ChIP_eaf3bar_H3K36me3__S92.bin_25.smooth_75Bulk.bw --skipZeros -b 2000 -a 2000 --sortRegions descend -o EAF3_ASH1_nonK27genes_H3K36me3_June2026_V2.gz --outFileNameMatrix EAF3_ASH1_nonK27genes_H3K36me3_June2026_V2.tab

plotHeatmap -m EAF3_ASH1_nonK27genes_H3K36me3_June2026_V2.gz -o EAF3_ASH1_nonK27genes_H3K36me3_June2026_V2.png --sortRegions descend --sortUsingSamples 2 --heatmapHeight 5  --heatmapWidth 3  --outFileSortedRegions EAF3_ASH1_nonK27genes_H3K36me3_June2026_V2_sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --colorMap 'YlOrBr' --zMax 8 --samplesLabel "WT"  "cdp6" "mrg15" "eaf3"

computeMatrix reference-point -p 12 -R "/home/ry00555/Research/Genomes/HeatmapGeneFiles/NonK27genes.bed" -S ../BigWigs/153-107_ChIP_WT_H3K36me3_.bin_25.smooth_50Bulk.bw  ../BigWigs/155-80_ChIP_C9_H3K36me3__S80.bin_25.smooth_75Bulk.bw ../BigWigs/155-77_ChIP_C8_H3K36me3__S77.bin_25.smooth_75Bulk.bw ../BigWigs/155-92_ChIP_eaf3bar_H3K36me3__S92.bin_25.smooth_75Bulk.bw --skipZeros -b 2000 -a 2000 --sortRegions descend -o EAF3_nonK27genes_H3K36me3_June2026_V2.gz --outFileNameMatrix EAF3_nonK27genes_H3K36me3_June2026_V2.tab

plotHeatmap -m EAF3_nonK27genes_H3K36me3_June2026_V2.gz -o EAF3_nonK27genes_H3K36me3_June2026_V2.png --sortRegions descend --sortUsingSamples 2 --heatmapHeight 5  --heatmapWidth 3  --outFileSortedRegions EAF3_nonK27genes_H3K36me3_June2026_V2.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --colorMap 'YlOrBr' --zMax 8 --samplesLabel "WT"  "cdp6" "mrg15" "eaf3"


computeMatrix reference-point -p 12 \
    -R "/home/ry00555/Research/Genomes/HeatmapGeneFiles/Bicocca_ASH1_noK27_wCM.bed" \
    -S "/scratch/ry00555/Run153/BigWigs/153-26_ChIP_WT_H3K27me3_Rep1_S26_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw" \
    "/scratch/ry00555/Run153/BigWigs/153-107_ChIP_WT_H3K36me3__S100_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw" \
    "/scratch/ry00555/Run153/BigWigs/153-38_ChIP_set7_H3K27me3__S38_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw" \
    "/scratch/ry00555/Run153/BigWigs/153-40_ChIP_set7_H3K36me3__S40_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw" \
    "/scratch/ry00555/Run153/BigWigs/153-35_ChIP_set2_H3K27me3_Rep2_S35_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw" \
    "/scratch/ry00555/Run153/BigWigs/153-90_ChIP_set2set7_H3K36me3_Rep1_S84_L002_R1_001_val_1.bin_25.smooth_50Bulk.bw" \
    "../BigWigs/ash1_Y888F_H3K27me2_3_scaled.bigwig" \
    "../BigWigs/cat_ash1_H3K36me3_scaled.bigwig" \
    --skipZeros -b 2000 -a 2000 \
    --sortRegions descend \
    -o ASH1_SET2_K27me3_H3K36me3_nonH3K27me3genes_V1.gz \
    --outFileNameMatrix ASH1_SET2_K27me3_H3K36me3_nonH3K27me3genes_V1.tab

plotHeatmap -m ASH1_SET2_K27me3_H3K36me3_nonH3K27me3genes_V1.gz \
    -o ASH1_SET2_K27me3_H3K36me3_nonH3K27me3genes_V1.png \
    --sortRegions descend \
    --sortUsingSamples 1 2 \
    --heatmapHeight 8 \
    --heatmapWidth 3 \
    --outFileSortedRegions ASH1_SET2_K27me3_H3K36me3_nonH3K27me3genes_V1.bed \
    --startLabel "5'" \
    --endLabel "3'" \
    --boxAroundHeatmaps no \
    --colorMap 'Greens' 'YlOrBr' 'Greens' 'YlOrBr' 'Greens' 'YlOrBr' 'Greens' 'YlOrBr' \
    --zMax 20 10 20 8 20 10 20 10 \
    --samplesLabel "WT H3K27me3" "WT H3K36me3" "set7 H3K27me3" "set7 H3K36me3" "set2 H3K27me3" "set2set7 H3K36me3" "Y888F H3K27me2/3" "ASH1cat H3K36me3"


# Set your file paths
BW1="sample1.bigwig"
BW2="sample2.bigwig"
OUTDIR="correlation_output"

mkdir -p $OUTDIR

# Step 1: Compute genome-wide bin scores
multiBigwigSummary bins \
  -b ASH1cat_H3K36me3_GSM3330994_Galaxy198-_bamCoverage_on_GTGCCA_N6875_ash1_H3K36me3_.bigwig 155-92_ChIP_eaf3bar_H3K36me3__S92.bin_25.smooth_75Bulk.bw \
  -o ash1cat_eaf3_May26_summary.npz \
  --binSize 10000 \
  -p 4 \
  --outRawCounts ash1cat_eaf3_May26_summary.tab

# Step 3: Plot scatter plot
plotCorrelation \
  -in ash1cat_eaf3_May26_summary.npz \
  --corMethod pearson \
  --skipZeros \
  --whatToPlot scatterplot \
  -o ash1cat_eaf3_May26_correlation_scatter.png
