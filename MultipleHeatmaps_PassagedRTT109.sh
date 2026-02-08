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

 computeMatrix scale-regions -p 12 -R "$OUTDIR/Heatmaps/RTT109_K27domains_WT_H3K27me3_consensus_MACSpeaks_V1_RTT109sorted.bed" -S $BWDIR/153-27_ChIP_WT_H3K27me3_Rep2_.bin_25.smooth_75Bulk.bw $BWDIR/147-47_ChIP_KOrtt109P13A6P4_rtt109hph_H3K27me3_Rep1_Nc_24hrVMMON_S47.bin_25.smooth_50_Q30.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/149-42__rtt109hph_H3K27me3_Rep5_S42_L004.bin_25.smooth_50_Q30.bw $BWDIR/150-73_rtt109_H3K27me3_rep13.bin_25.smooth_50.bw $BWDIR/150-65_rtt109_H3K27me3_rep11.bin_25.smooth_50.bw $BWDIR/150-69_rtt109_H3K27me3_rep12.bin_25.smooth_50.bw --skipZeros -b 500 -a 500 --sortUsingSamples 2 --sortRegions descend -o $OUTDIR/Heatmaps/RTT109_Passged_K27domains_V2.gz --outFileNameMatrix $OUTDIR/Heatmaps/RTT109_Passged_K27domains_V2.tab


plotHeatmap -m $OUTDIR/Heatmaps/RTT109_Passged_K27domains_V2.gz -o $OUTDIR/Heatmaps/RTT109_Passged_K27domains_V2.png --sortRegions descend --sortUsingSamples 2 --heatmapHeight 5  --heatmapWidth 2  --outFileSortedRegions $OUTDIR/Heatmaps/RTT109_Passged_K27domains_V2_sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --samplesLabel  "WT 153-127" "∆rtt109 147-47"   "∆rtt109 149-42" "∆rtt109 150-73" "∆rtt109 150-65" "∆rtt109 150-69" --colorMap 'Greens'  --zMax 15

computeMatrix reference-point -p 12 -R "$OUTDIR/Heatmaps/K36ONLYonK27genes_StrainsTogether_n573_sortRTT109.bed"  -S $BWDIR/153-27_ChIP_WT_H3K27me3_Rep2_.bin_25.smooth_75Bulk.bw $BWDIR/147-47_ChIP_KOrtt109P13A6P4_rtt109hph_H3K27me3_Rep1_Nc_24hrVMMON_S47.bin_25.smooth_50_Q30.bw /scratch/ry00555/RTT109PaperFigures/BigWigs/149-42__rtt109hph_H3K27me3_Rep5_S42_L004.bin_25.smooth_50_Q30.bw $BWDIR/150-73_rtt109_H3K27me3_rep13.bin_25.smooth_50.bw $BWDIR/150-65_rtt109_H3K27me3_rep11.bin_25.smooth_50.bw $BWDIR/150-69_rtt109_H3K27me3_rep12.bin_25.smooth_50.bw  --skipZeros  -a 1000  -b 2000 --sortUsingSamples 2 --sortRegions descend -o $OUTDIR/Heatmaps/RTT109_Passged_K27genes_V2.gz --outFileNameMatrix $OUTDIR/Heatmaps/RTT109_Passged_K27genes_V2.tab

plotHeatmap -m $OUTDIR/Heatmaps/RTT109_Passged_K27genes_V2.gz -o $OUTDIR/Heatmaps/RTT109_Passged_K27genes_V2.png --sortRegions descend --sortUsingSamples 2 --heatmapHeight 5  --heatmapWidth 2  --outFileSortedRegions $OUTDIR/Heatmaps/RTT109_Passged_K27genes_V2_sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --samplesLabel  "WT 153-127" "∆rtt109 147-47"  "∆rtt109 149-42" "∆rtt109 150-73" "∆rtt109 150-65" "∆rtt109 150-69" --colorMap 'Greens'  --zMax 20

computeMatrix reference-point -p 12 -R "$OUTDIR/Heatmaps/K36ONLYonK27genes_StrainsTogether_n573_sortRTT109.bed" \
-S $BWDIR/150-62_WT_H3K36me3_rep8.bin_25.smooth_50.bw $BWDIR/150-70_rtt109_H3K36me3_rep7.bin_25.smooth_50.bw $BWDIR/147-48_ChIP_KOrtt109P13A6P4_rtt109hph_H3K36me3_Rep1_Nc_24hrVMMON_S48.bin_25.smooth_50_Q30.bw $BWDIR/147-56_ChIP_337rtt109P3_rtt109hph_H3K36me3_Rep1_Nc_24hrVMMON_S56.bin_25.smooth_50_Q30.bw $BWDIR/150-66_rtt109_H3K36me3_rep6.bin_25.smooth_50.bw $BWDIR/150-74_rtt109_H3K36me3_rep8.bin_25.smooth_50.bw  --skipZeros  -a 1000  -b 2000 --sortUsingSamples 3 --sortRegions descend -o $OUTDIR/Heatmaps/RTT109_Passged_K27genes_V2.gz --outFileNameMatrix $OUTDIR/Heatmaps/RTT109_Passged_K27genes_V2.tab

plotHeatmap -m $OUTDIR/Heatmaps/RTT109_Passged_K27genes_V2.gz -o $OUTDIR/Heatmaps/RTT109_Passged_K27genes_V2.png --sortRegions descend --sortUsingSamples 3 --heatmapHeight 5  --heatmapWidth 2  --outFileSortedRegions $OUTDIR/Heatmaps/RTT109_Passaged_H3K36_K27enes_V2_sorted.bed  --startLabel "5'"  --endLabel "3'" --boxAroundHeatmaps no  --samplesLabel  "WT 153-62" "∆rtt109 150-70"  "∆rtt109 147-48" "∆rtt109 147-56" "∆rtt109 150-66" "∆rtt109 150-74" --colorMap 'YlOrBr'  --zMax 8
