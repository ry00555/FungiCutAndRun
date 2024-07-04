#!/bin/bash
#SBATCH --job-name=zl_heatmap
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zlewis@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=20gb
#SBATCH --time=1:00:00
#SBATCH --output=./heatmap.%j.out
#SBATCH --error=./heatmap.%j.err

cd $SLURM_SUBMIT_DIR

ml deepTools/3.5.2-foss-2022a

#Peakfile for H3K9me3 is based on CutNRun Sample 114-19
####Plotting scaled enrichment across all het domains
computeMatrix scale-regions -p 12 -R Figure2_K27PeakFile.txt -S 138-57_ChIP_WT_H3K27me3_Rep3_6252_S56_L001_R1_001_val_1.fq.gz.bin_25.smooth_50Bulk.bw 137-75_ChIP_qa-suz12_H3K27me3_Rep1_S72_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw 138-61_ChIP_set-7_H3K27me3_Rep3_6252_S60_L001_R1_001_val_1.fq.gz.bin_25.smooth_50Bulk.bw 135-83_ChIP_rtt109_H3K27me3_Rep2_S79_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw 138-51_ChIP_H3K56R40_H3K27me3_Rep2_6252_S50_L001_R1_001_val_1.fq.gz.bin_25.smooth_50Bulk.bw	-b 1000 -a 1000 -o Rtt109vs24hrDenovo_qasuz12sort.matrix --missingDataAsZero -bs 10 --sortUsingSamples 2
     # --sortUsing sum --sortUsingSamples 5 2
    	#--sortRegions descend --sortUsing sum \

#https://deeptools.readthedocs.io/en/develop/content/tools/plotHeatmap.html
plotHeatmap -m Rtt109vs24hrDenovo_qasuz12sort.matrix \
	-o Rtt109vs24hrDenovo_qasuz12sort_recolor.pdf \
			--outFileSortedRegions Rtt109vs24hrDenovo_qasuz12sort.txt \
			--startLabel "5'" \
				--endLabel "3'" \
			--samplesLabel "WT" "24hr qa-suz-12" "∆set-7" "∆rtt109" "H3K56R" \
					--heatmapHeight 7 \
            --colorMap 'Blues' 'Blues' 'Greens' 'Reds' 'Purples' --zMax 10 10 10 10 10 10
