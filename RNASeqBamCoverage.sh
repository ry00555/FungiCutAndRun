#!/bin/bash
#SBATCH --job-name=RNASeqBamCoverage	                      # Job name
#SBATCH --partition=batch		                            # Partition (queue) name
#SBATCH --ntasks=1			                                # Single task job
#SBATCH --cpus-per-task=6		                            # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=24gb			                                # Total memory for job
#SBATCH --time=72:00:00  		                            # Time limit hrs:min:sec
#SBATCH --output=../RNASeqBamCoverage.%j.out
#SBATCH --error=../RNASeqBamCoverage.%j.err
#SBATCH --mail-user=ry00555@uga.edu                    # Where to send mail - # replace cbergman with your myid
#SBATCH --mail-type=ALL                            # Mail events (BEGIN, END, FAIL, ALL)


#set input and output directory variables
OUTDIR="/scratch/ry00555/RNASeqBamCoverage/Eaf3"
FASTQ="/scratch/ry00555/RNASeqBamCoverage/Eaf3/FASTQ"
GENOME="/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_00182925.2plusHphplusBarplusTetO_his3masked.fna"

cd $SLURM_SUBMIT_DIR

#if output directory doesn't exist, create it
 if [ ! -d $OUTDIR ]
 then
     mkdir -p $OUTDIR
 fi

#loading modules
ml SRA-Toolkit

#You can use a loop to prefetch each SRR ID and then subsequently run fastq-dump for each downloaded file. Here's how you can do it:


# List of SRR IDs
SRR_IDS=(
     "SRR9027634" "SRR9027635" "SRR9027636" "SRR9027653" "SRR9027655"
     "SRR9027701" "SRR9044213" "SRR9044244" "SRR9044324" "SRR10916182"
     "SRR10916183" "SRR10916184" "SRR10916163" "SRR10916164" "SRR10916165"
     "SRR8444005" "SRR8444042" "SRR8443998" "SRR12614222" "SRR12614223"
     "SRR12614224" "SRR12614225" "SRR12614226")
#
# # Prefetch SRA files
#for SRR_ID in "${SRR_IDS[@]}"; do
#    prefetch -O "${OUTDIR}" "${SRR_ID}"
#done

#downloading reference genome
#curl -s ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/Neurospora_crassa/latest_assembly_versions/GCA_000182925.2_NC12/GCA_000182925.2_NC12_genomic.fna.gz  | gunzip -c > ${OUTDIR}/NC12_genome.fna
#curl -s ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/182/925/GCA_000182925.2_NC12/GCA_000182925.2_NC12_genomic.gtf.gz | gunzip -c > ${OUTDIR}/NC12_annotation.gtf

#DOWNLOADING SRA
#prefetch -O SRR9027634 SRR9027635 SRR9027636 SRR9027653 SRR9027655 SRR9027701 SRR9044213 SRR9044244 SRR9044324 SRR10916182 SRR10916183 SRR10916184 SRR10916163 SRR10916164 SRR10916165 SRR8444005 SRR8444042 SRR8443998 SRR12614222 SRR12614223 SRR12614224 SRR12614225 SRR12614226

##convert to fastq
 #for SRR_ID in "${SRR_IDS[@]}"; do
 #fastq-dump --split-files --gzip ${OUTDIR}/${SRR_ID}/${SRR_ID}.sra -O ${FASTQ}
 #done
#single commmand
# fastq-dump --split-files --gzip	${OUTDIR}/SRR9027634
#

#Finish the reminaing SRR's as the script was interrupted halfway due to GACRC maintainence
#ml SRA-Toolkit
#fastq-dump --split-files --gzip "${OUTDIR}/SRR12614226/SRR12614226.sra" -O ${OUTDIR}/FASTQ
#fastq-dump --split-files --gzip "${OUTDIR}/SRR12614227/SRR12614227.sra" -O ${OUTDIR}/FASTQ

# #process reads using trimGalore

#
mkdir "${OUTDIR}/SortedBamFiles"
mkdir "${OUTDIR}/BigWigs"
mkdir "${OUTDIR}/Peaks"
mkdir "${OUTDIR}/TrimmedReads"
mkdir "$OUTDIR/Matrices"
mkdir "$OUTDIR/Heatmaps"

 ml BWA
ml SAMtools
ml Trim_Galore
 trim_galore --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${FASTQ}/*fastq\.gz
 #trim_galore --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${FASTQ}/SRR12614227*fastq\.gz
 #trim_galore --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${FASTQ}/SRR12614226*fastq\.gz


#
FILES="${OUTDIR}/TrimmedReads/*.fq.gz" #Don't forget the *
# #
#
# #
# #Iterate over the files
for f in $FILES#
do
#
# 	#Examples to Get Different parts of the file name
# 		#See here for details: http://tldp.org/LDP/abs/html/refcards.html#AEN22664
		#${string//substring/replacement}
# 		#dir=${f%/*}

	file=${f##*/}
# 	#remove ending from file name to create shorter names for bam files and other downstream output
 	name=${file/*_1.fq.gz/}
#
# #
# # 	# File Vars
# # 	#use sed to get the name of the second read matching the input file
 	read2=$(echo "$f" | sed 's/_1\.fq\.gz/_2\.fq\.gz/g')
# 	#variable for naming bam file
  	bam="${OUTDIR}/SortedBamFiles/${name}.bam"
# 	#variable name for bigwig output
 	bigwig="${OUTDIR}/BigWigs/${name}"
 QualityBam="${OUTDIR}/SortedBamFiles/${name}_Q30.bam"
# #
#
# #
 bwa mem -M -v 3 -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIR/SortedBamFiles/tempReps -o "$bam" -
 samtools index "$bam"
# #
samtools view -b -q 30 $bam > "$QualityBam"
samtools index "$QualityBam"
# #
# # ############################
# # # #deeptools
#
 ml deepTools
 bamCoverage -p $THREADS -bs $BIN --normalizeUsing CPM --ignoreDuplicates --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"
# #
  done
#

#Wildtype

#SRR8444005_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw
#SRR8444042_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw
#SRR8443998_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw
#NCU07496

#SRR10916163_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw
#SRR10916164_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw
#SRR10916165_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw
#NCU06787

#SRR9044213_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw
#SRR9044244_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw
#SRR9044324_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw
#NCU06788

#SRR9027653_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw
#SRR9027655_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw
#SRR9027701_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw


#"../BigWigs/SRR8444005_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw",
#"../BigWigs/SRR8444042_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw",
#"../BigWigs/SRR8443998_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw",
#"../BigWigs/SRR10916163_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw",
#"../BigWigs/SRR10916164_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw",
#"../BigWigs/SRR10916165_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw",
#"../BigWigs/SRR9044213_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw",
#"../BigWigs/SRR9044244_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw",
#"../BigWigs/SRR9044324_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw",
#"../BigWigs/SRR9027653_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw",
#"../BigWigs/SRR9027655_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw",
#"../BigWigs/SRR9027701_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw"

Wildtype Wildtype Wildtype ∆set-7 ∆set-7 ∆set-7 ∆mrg15 ∆mrg15 ∆mrg15 ∆cdp-6 ∆cdp-6 ∆cdp-6
for file_path in "${OUTDIR}/BigWigs"/*.bw; do
#   # Get the base name of the file
   BW_name=$(basename "${file_path}" .bw)
#
#   # Remove the file extension to create the sample ID
   BW_id="${BW_name%.*}"
#   # Replace special characters with underscores in the sample ID
  BW_id=${BW_id//[^a-zA-Z0-9]/_}
#
#  # Limit the length of the sample ID to avoid long filenames
  BW_id=${BW_id:0:50}
#  # Compute matrix for the reference-point TSS
  computeMatrix scale-regions --startLabel TSS --endLabel TES -a 500 -b 500 --unscaled5prime 500 --unscaled3prime 500 -S ../BigWigs/SRR8444005_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR8444042_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR10916163_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR10916164_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR10916165_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR9044213_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR9044244_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR9044324_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR9027653_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR9027655_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR9027701_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw -R /scratch/ry00555/neurospora.bed --skipZeros --missingDataAsZero -o EAF3_V3.gz
plotHeatmap --matrixFile ../Matrices/EAF3_V3.gz --outFileName EAF3_processed_V3_regions_hclust.png --samplesLabel Wildtype Wildtype Wildtype ∆set-7 ∆set-7 ∆set-7 ∆mrg15 ∆mrg15 ∆mrg15 ∆cdp-6 ∆cdp-6 ∆cdp-6 --hclust 1 --regionsLabel "CM002236.1" "CM002237.1" "CM002238.1" "CM002239.1" "CM002240.1" "CM002241.1" "CM002242.1" --colorMap spring summer autumn winter --missingDataColor white --sortUsingSamples 1 2 3

plotHeatmap --matrixFile ../Matrices/EAF3_V3.gz --outFileName EAF3_processed_V3_hclust.png --samplesLabel Wildtype Wildtype Wildtype ∆set-7 ∆set-7 ∆set-7 ∆mrg15 ∆mrg15 ∆mrg15 ∆cdp-6 ∆cdp-6 ∆cdp-6 --hclust 1 --colorMap spring summer autumn winter --sortRegions descend --missingDataColor white --sortUsingSamples 1 2 3

#Remove ../BigWigs/SRR8443998_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw bc is bad WT

computeMatrix scale-regions -p 12 --startLabel TSS --endLabel TES -a 1000 -b 1000 --unscaled5prime 500 --unscaled3prime 500 -S ../BigWigs/SRR8444005_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR8444042_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR10916163_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR10916164_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR10916165_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR9044213_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR9044244_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR9044324_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR9027653_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR9027655_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR9027701_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw -R /scratch/ry00555/neurospora.bed --skipZeros --missingDataAsZero -o EAF3_V4.gz

plotHeatmap --matrixFile ../Matrices/EAF3_V4.gz --outFileName EAF3_processed_V6_wholeGenome_hclust.png --samplesLabel Wildtype Wildtype ∆set-7 ∆set-7 ∆set-7 ∆mrg15 ∆mrg15 ∆mrg15 ∆cdp-6 ∆cdp-6 ∆cdp-6 --hclust 1 --sortRegions descend --colorMap coolwarm --missingDataColor white --sortUsingSamples 1 2 --zMax auto --zMin auto --heatmapHeight 12 --whatToShow 'heatmap and colorbar'

  #zcat EAF3_V2.gz | awk '{for (i=1; i<=NF; i++) if ($i == "nan") $i=0; print}' | gzip > EAF3_processed_V2.gz


#   # Compute matrix for the reference-point TSS with specific regions (e.g., PRC2 genes)
#   computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S "${file_path}" -R "/scratch/ry00555/heatmapPRC2genes.bed" --skipZeros -o "${OUTDIR}/Matrices/PRC2genes_matrix_${BW_id}.gz"
computeMatrix scale-regions -p 12  --startLabel TSS --endLabel TES -a 1000 -b 1000 -S ../BigWigs/SRR8444005_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR8444042_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR10916163_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR10916164_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR10916165_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR9044213_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR9044244_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR9044324_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR9027653_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR9027655_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw ../BigWigs/SRR9027701_2_val_2.fq.gz.bin_25.smooth_50Bulk.bw -R /scratch/ry00555/heatmapPRC2genes.bed --skipZeros -o EAF3_K27me3Domains2_V3.gz

plotHeatmap --matrixFile EAF3_K27me3Domains2_V3.gz --outFileName ../Heatmaps/EAF3_processed_V5_K27me3regions_hclust.png --samplesLabel Wildtype Wildtype ∆set-7 ∆set-7 ∆set-7 ∆mrg15 ∆mrg15 ∆mrg15 ∆cdp-6 ∆cdp-6 ∆cdp-6 --hclust 1 --sortRegions descend --colorMap coolwarm --missingDataColor white --sortUsingSamples 1 2 --zMax auto --zMin auto --heatmapHeight 12 --whatToShow 'heatmap and colorbar'

zcat EAF3_K27me3Domains2_V2.gz| awk '{for (i=1; i<=NF; i++) if ($i == "nan") $i=0; print}' | gzip > EAF3_K27me3Domains2_processed_V2.gz

#   # Preprocess the matrix files to replace nan values with zeros
   zcat "${OUTDIR}/Matrices/matrix_${BW_id}_V2.gz" | awk '{for (i=1; i<=NF; i++) if ($i == "nan") $i=0; print}' | gzip > "${OUTDIR}/Matrices/matrix_${BW_id}_processed_V2.gz"
   zcat "${OUTDIR}/Matrices/matrix_${BW_id}_K27me3Domains2_V2.gz" | awk '{for (i=1; i<=NF; i++) if ($i == "nan") $i=0; print}' | gzip > "${OUTDIR}/Matrices/PRC2genes_matrix_${BW_id}_processed_V2.gz"
#
#   # Plot heatmaps for the reference-point TSS
 plotHeatmap --matrixFile ../Matrices/EAF3_K27me3Domains2_processed_V2.gz --outFileName EAF3_K27me3Domains2_processed_V2_hclust.png --samplesLabel Wildtype Wildtype Wildtype ∆set-7 ∆set-7 ∆set-7 ∆mrg15 ∆mrg15 ∆mrg15 ∆cdp-6 ∆cdp-6 ∆cdp-6 --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1 2 3
#
#   # Plot heatmaps for the reference-point TSS with specific regions (e.g., PRC2 genes)
# plotHeatmap --matrixFile "${OUTDIR}/Matrices/PRC2genes_matrix_${BW_id}.gz" --outFileName "${OUTDIR}/Heatmaps/${BW_id}_hclust_PRC2genes.png" --samplesLabel "${BW_name}" --hclust 1 --colorMap Greens
#
#  plotHeatmap --matrixFile "${OUTDIR}/Matrices/matrix_${BW_id}_processed.gz" --outFileName "${OUTDIR}/Heatmaps/${BW_id}_hclust.png" --samplesLabel "${BW_name}" --hclust 1 --colorMap Greens --sortRegions descend --missingDataColor white --sortUsingSamples 1
#
#   # Plot heatmaps for the reference-point TSS with specific regions (e.g., PRC2 genes)
 plotHeatmap --matrixFile "${OUTDIR}/Matrices/PRC2genes_matrix_${BW_id}_processed_V2.gz" --outFileName "${OUTDIR}/Heatmaps/${BW_id}_hclust_PRC2genes_ZLMethod_V2.png" --samplesLabel Wildtype Wildtype Wildtype ∆set-7 ∆set-7 ∆set-7 ∆mrg15 ∆mrg15 ∆mrg15 ∆cdp-6 ∆cdp-6 ∆cdp-6 --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1 2 3
#
 done
