# Step 1: Create a chromosome name mapping file (old_name new_name)
cat > /home/ry00555/Research/Genomes/Ncrassa_chrom_rename.txt << 'EOF'
Supercontig_12.1 CM002236.1
Supercontig_12.2 CM002237.1
Supercontig_12.3 CM002238.1
Supercontig_12.4 CM002239.1
Supercontig_12.5 CM002240.1
Supercontig_12.6 CM002241.1
Supercontig_12.7 CM002242.1
Supercontig_12.8 KI440765.1
Supercontig_12.9 KI440766.1
Supercontig_12.10 KI440767.1
Supercontig_12.11 KI440768.1
Supercontig_12.12 KI440769.1
Supercontig_12.13 KI440770.1
Supercontig_12.14 KI440771.1
Supercontig_12.15 KI440772.1
Supercontig_12.16 KI440773.1
Supercontig_12.17 KI440774.1
Supercontig_12.18 KI440775.1
Supercontig_12.19 KI440776.1
Supercontig_12.20 KI440777.1
EOF

# Step 2: Rename chromosomes in the ash1 bigwig using ucsc-bigwiginfo + bigWigLiftOver
# Easiest approach: use bigWigRenameChroms (if available) or convert via bedGraph

ASH1_BW="ASH1cat_H3K36me3_GSM3330994_Galaxy198-_bamCoverage_on_GTGCCA_N6875_ash1_H3K36me3_.bigwig"
EAF3_BW="155-92_ChIP_eaf3bar_H3K36me3__S92.bin_25.smooth_75Bulk.bw"

# Convert ash1 bigwig to bedGraph, rename chroms, convert back
bigWigToBedGraph ASH1cat_H3K36me3_GSM3330994_Galaxy198-_bamCoverage_on_GTGCCA_N6875_ash1_H3K36me3_.bigwig cat_ash1_H3K36me3_temp.bedGraph

# Rename using awk
awk 'NR==FNR{map[$1]=$2; next} ($1 in map){$1=map[$1]} 1' OFS="\t" \
    /home/ry00555/Research/Genomes/Ncrassa_chrom_rename.txt cat_ash1_H3K36me3_temp.bedGraph > cat_ash1_H3K36me3_renamed.bedGraph

    awk 'NR==FNR{map[$1]=$2; next} ($1 in map){$1=map[$1]} 1' OFS="\t" \
        /home/ry00555/Research/Genomes/Ncrassa_chrom_rename.txt \
        /home/ry00555/Research/Genomes/HeatmapGeneFiles/GSM3330995_ASH-1-marked_FullTranscript_1970_.bed \
        > /home/ry00555/Research/Genomes/HeatmapGeneFiles/ASH1markedgenes_renamed.bed

# Get chrom sizes from the eaf3 bw (already uses NCBI names)
bigWigInfo -chroms 155-92_ChIP_eaf3bar_H3K36me3__S92.bin_25.smooth_75Bulk.bw 2>&1 | grep -v "^#" | awk 'NF==3{print $1"\t"$3}' \
    > eaf3_chrom_sizes.txt
#sort bedgraph
    LC_COLLATE=C sort -k1,1 -k2,2n cat_ash1_H3K36me3_renamed.bedGraph > cat_ash1_H3K36me3_renamed_sorted.bedGraph


# Convert renamed bedGraph back to bigwig
bedGraphToBigWig cat_ash1_H3K36me3_renamed_sorted.bedGraph eaf3_chrom_sizes.txt cat_ash1_H3K36me3_renamed.bigwig

# Step 3: Run the correlation
multiBigwigSummary bins \
    -b cat_ash1_H3K36me3_renamed.bigwig 155-92_ChIP_eaf3bar_H3K36me3__S92.bin_25.smooth_75Bulk.bw \
    -o cat_ash1_eaf3_H3K36me3_summary.npz \
    --binSize 10000 \
    -p 4 \
    --outRawCounts cat_ash1_eaf3_H3K36me3_summary.tab

    #You will get this error: The following chromosome names did not match between the bigwig files
#chromosome	length
#     KC683708.1	     64840
#*Warning*
#The resulting bed file does not contain information for the chromosomes that were not common between the bigwig files
#Number of bins found: 4114


plotCorrelation \
    -in cat_ash1_eaf3_H3K36me3_summary.npz \
    --corMethod spearman \
    --skipZeros \
    --whatToPlot scatterplot \
    -o cat_ash1_eaf3_H3K36me3_correlation_scatter_May2026__spearman_V1.png

    # Get total signal and genome size for each bigwig
    bigWigInfo cat_ash1_H3K36me3_renamed.bigwig
    bigWigInfo 155-92_ChIP_eaf3bar_H3K36me3__S92.bin_25.smooth_75Bulk.bw
    bigWigInfo cat_ash1_H3K36me3_renamed.bigwig
bigWigInfo 155-92_ChIP_eaf3bar_H3K36me3__S92.bin_25.smooth_75Bulk.bw
version: 4
isCompressed: yes
isSwapped: 0
primaryDataSize: 16,530,054
primaryIndexSize: 119,620
zoomLevels: 9
chromCount: 20
basesCovered: 41,037,538
mean: 266.925392
min: 0.000000
max: 2759.459961
std: 228.904615
version: 4
isCompressed: yes
isSwapped: 0
primaryDataSize: 3,924,000
primaryIndexSize: 10,716
zoomLevels: 4
chromCount: 21
basesCovered: 41,110,515
mean: 4.179062
min: 0.000000
max: 1753.050049
std: 6.309173

bigwigCompare \
    -b1 cat_ash1_H3K36me3_renamed.bigwig \
    -b2 cat_ash1_H3K36me3_renamed.bigwig \
    --operation mean \
    --scaleFactors "0.01566:0.01566" \
    -o cat_ash1_H3K36me3_scaled.bigwig \
    -p 4
    multiBigwigSummary BED-file \
        -b cat_ash1_H3K36me3_scaled.bigwig 155-92_ChIP_eaf3bar_H3K36me3__S92.bin_25.smooth_75Bulk.bw \
        -o cat_ash1_eaf3_H3K36me3_ASH1genes_scaled_summary.npz \
        --BED /home/ry00555/Research/Genomes/HeatmapGeneFiles/ASH1markedgenes_renamed.bed \
        -p 4 \
        --outRawCounts cat_ash1_eaf3_H3K36me3_ASH1genes_scaled_summary.tab



        # Sort the .tab output by ash1 signal, take top quartile
        awk 'NR==1{header=$0; next} {print}' cat_ash1_eaf3_H3K36me3_ASH1genes_scaled_summary.tab | \
            sort -k4,4rn | \
            head -n 492 | \
            sort -k1,1 -k2,2n > top25_ash1_eaf3.tab
                awk '{print $1"\t"$2"\t"$3}' top25_ash1_eaf3.tab > top25_regions.bed

multiBigwigSummary BED-file \
    -b cat_ash1_H3K36me3_scaled.bigwig 155-92_ChIP_eaf3bar_H3K36me3__S92.bin_25.smooth_75Bulk.bw \
    -o top25_ash1_eaf3_summary.npz \
    --BED top25_regions.bed \
    -p 4 \
    --outRawCounts top25_ash1_eaf3_summary.tab

# Plot
plotCorrelation \
    -in top25_ash1_eaf3_summary.npz \
    --corMethod pearson \
    --skipZeros \
    --whatToPlot scatterplot \
    -o top25_ash1_eaf3_scatter.png

plotCorrelation \
    -in top25_ash1_eaf3_summary.npz \
    --corMethod spearman \
    --skipZeros \
    --whatToPlot scatterplot \
    -o top25_ash1_eaf3_spearman_scatter.png


    multiBigwigSummary BED-file \
        -b cat_ash1_H3K36me3_renamed.bigwig 155-92_ChIP_eaf3bar_H3K36me3__S92.bin_25.smooth_75Bulk.bw \
        -o cat_ash1_eaf3_H3K36me3_genes_summary.npz \
        --BED /home/ry00555/Research/Genomes/HeatmapGeneFiles/ASH1markedgenes_renamed.bed \
        -p 4 \
        --outRawCounts cat_ash1_eaf3_ASH1marked_genesONLY_summary.tab

        plotCorrelation \
            -in cat_ash1_eaf3_H3K36me3_genes_summary.npz \
            --corMethod pearson \
            --skipZeros \
            --whatToPlot scatterplot \
            -o cat_ash1_eaf3_ASH1marked_genesONLY_pearson_May2026_V1.png
