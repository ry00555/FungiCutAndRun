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

cd $SLURM_SUBMIT_DIR

source config.txt

#if output directory doesn't exist, create it
# if [ ! -d $OUTDIR ]
# then
#     mkdir -p $OUTDIR
# fi

#loading modules
#ml SRA-Toolkit

#You can use a loop to prefetch each SRR ID and then subsequently run fastq-dump for each downloaded file. Here's how you can do it:


# List of SRR IDs
# #SRR_IDS=(
#     "SRR9027634" "SRR9027635" "SRR9027636" "SRR9027653" "SRR9027655"
#     "SRR9027701" "SRR9044213" "SRR9044244" "SRR9044324" "SRR10916182"
#     "SRR10916183" "SRR10916184" "SRR10916163" "SRR10916164" "SRR10916165"
#     "SRR8444005" "SRR8444042" "SRR8443998" "SRR12614222" "SRR12614223"
#     "SRR12614224" "SRR12614225" "SRR12614226"
# )
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
# for SRR_ID in "${SRR_IDS[@]}"; do
# fastq-dump --split-files --gzip ${OUTDIR}/${SRR_ID}/${SRR_ID}.sra -O ${FASTQ}
# done
#single commmand
# fastq-dump --split-files --gzip	${OUTDIR}/SRR9027634
#

#Finish the reminaing SRR's as the script was interrupted halfway due to GACRC maintainence
ml SRA-Toolkit
fastq-dump --split-files --gzip "${OUTDIR}/SRR12614226/SRR12614226.sra" -O ${OUTDIR}/FASTQ
fastq-dump --split-files --gzip "${OUTDIR}/SRR12614227/SRR12614227.sra" -O ${OUTDIR}/FASTQ

# #process reads using trimGalore

#
#mkdir "${OUTDIR}/SortedBamFiles"
#mkdir "${OUTDIR}/BigWigs"
#mkdir "${OUTDIR}/Peaks"
#mkdir "${OUTDIR}/TrimmedReads"
#mkdir "$OUTDIR/Matrices"
#mkdir "$OUTDIR/Heatmaps"

 ml BWA
ml SAMtools
ml Trim_Galore
 trim_galore --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${FASTQ}/*fastq\.gz
 trim_galore --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${FASTQ}/SRR12614227*fastq\.gz
 trim_galore --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${FASTQ}/SRR12614226*fastq\.gz


#
FILES="${OUTDIR}/TrimmedReads/*.fq.gz" #Don't forget the *
# #
#
# #
# #Iterate over the files
for f in $FILES
do
#
# 	#Examples to Get Different parts of the file name
# 		#See here for details: http://tldp.org/LDP/abs/html/refcards.html#AEN22664
		#${string//substring/replacement}
# 		#dir=${f%/*}

	file=${f##*/}
	#remove ending from file name to create shorter names for bam files and other downstream output
	name=${file/*_1.fq.gz/}

#
# 	# File Vars
# 	#use sed to get the name of the second read matching the input file
	read2=$(echo "$f" | sed 's/_1\.fq\.gz/_2\.fq\.gz/g')
	#variable for naming bam file
 	bam="${OUTDIR}/SortedBamFiles/${name}.bam"
	#variable name for bigwig output
	bigwig="${OUTDIR}/BigWigs/${name}"
	#QualityBam="${OUTDIR}/SortedBamFiles/${name}_Q30.bam"
#

#
bwa mem -M -v 3 -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIR/SortedBamFiles/tempReps -o "$bam" -
samtools index "$bam"
#
# #samtools view -b -q 30 $bam > "$QualityBam"
# #samtools index "$QualityBam"
#
# ############################
# # #deeptools

ml deepTools
# #use these parameters for ChIP data
bamCoverage -p $THREADS -bs $BIN --normalizeUsing CPM --ignoreDuplicates --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"
#
#  done
#
# for file_path in "${OUTDIR}/BigWigs"/*.bw; do
#   # Get the base name of the file
#   BW_name=$(basename "${file_path}" .bw)
#
#   # Remove the file extension to create the sample ID
#   BW_id="${BW_name%.*}"
#   # Replace special characters with underscores in the sample ID
#  BW_id=${BW_id//[^a-zA-Z0-9]/_}
#
#  # Limit the length of the sample ID to avoid long filenames
#  BW_id=${BW_id:0:50}
#  # Compute matrix for the reference-point TSS
#  computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S "${file_path}" -R "/scratch/ry00555/neurospora.bed" --skipZeros -o "${OUTDIR}/Matrices/matrix_${BW_id}.gz"
#
#   # Compute matrix for the reference-point TSS with specific regions (e.g., PRC2 genes)
#   computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S "${file_path}" -R "/scratch/ry00555/heatmapPRC2genes.bed" --skipZeros -o "${OUTDIR}/Matrices/PRC2genes_matrix_${BW_id}.gz"
#
#   # Preprocess the matrix files to replace nan values with zeros
#   zcat "${OUTDIR}/Matrices/matrix_${BW_id}.gz" | awk '{for (i=1; i<=NF; i++) if ($i == "nan") $i=0; print}' | gzip > "${OUTDIR}/Matrices/matrix_${BW_id}_processed.gz"
#   zcat "${OUTDIR}/Matrices/PRC2genes_matrix_${BW_id}.gz" | awk '{for (i=1; i<=NF; i++) if ($i == "nan") $i=0; print}' | gzip > "${OUTDIR}/Matrices/PRC2genes_matrix_${BW_id}_processed.gz"
#
#   # Plot heatmaps for the reference-point TSS
#  plotHeatmap --matrixFile "${OUTDIR}/Matrices/matrix_${BW_id}.gz" --outFileName "${OUTDIR}/Heatmaps/${BW_id}_hclust.png" --samplesLabel "${BW_name}" --hclust 1 --colorMap Greens --sortRegions descend --missingDataColor white --sortUsingSamples 1
#
#   # Plot heatmaps for the reference-point TSS with specific regions (e.g., PRC2 genes)
#  plotHeatmap --matrixFile "${OUTDIR}/Matrices/PRC2genes_matrix_${BW_id}.gz" --outFileName "${OUTDIR}/Heatmaps/${BW_id}_hclust_PRC2genes.png" --samplesLabel "${BW_name}" --hclust 1 --colorMap Greens
#
#  plotHeatmap --matrixFile "${OUTDIR}/Matrices/matrix_${BW_id}_processed.gz" --outFileName "${OUTDIR}/Heatmaps/${BW_id}_hclust.png" --samplesLabel "${BW_name}" --hclust 1 --colorMap Greens --sortRegions descend --missingDataColor white --sortUsingSamples 1
#
#   # Plot heatmaps for the reference-point TSS with specific regions (e.g., PRC2 genes)
#  plotHeatmap --matrixFile "${OUTDIR}/Matrices/PRC2genes_matrix_${BW_id}_processed.gz" --outFileName "${OUTDIR}/Heatmaps/${BW_id}_hclust_PRC2genes.png" --samplesLabel "${BW_name}" --hclust 1 --colorMap Greens
#
#
# done
