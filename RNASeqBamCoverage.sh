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
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#set input and output directory variables
OUTDIR="/scratch/ry00555/RNASeqBamCoverage/Eaf3"

cd $SLURM_SUBMIT_DIR



#if output directory doesn't exist, create it
# if [ ! -d $OUTDIR ]
# then
#     mkdir -p $OUTDIR
# fi

#loading modules
ml SRA-Toolkit Subread



#!/bin/bash

# Assuming ${OUTDIR} is already defined

# Check if the input file exists
if [ ! -f "${OUTDIR}/SRAforEaf3Ash1Set7WTRco1withMeta.txt" ]; then
    echo "Error: Input file SRAforEaf3Ash1Set7WTRco1.txt not found!"
    exit 1
fi

# Read SRR IDs and corresponding descriptive titles from the input file
while IFS= read -r SRR_ID; do
    prefetch -O "${OUTDIR}" "${SRR_ID}"
done < "${OUTDIR}/SRAforEaf3Ash1Set7WTRco1withMeta.txt"

while IFS=$'\t' read -r SRR_ID Descriptive_Title; do
    # Perform fastq-dump for each SRR ID
    fastq-dump --split-files --gzip "${OUTDIR}/${SRR_ID}" -O "${OUTDIR}"

    # Rename the output files
    mv "${OUTDIR}/${SRR_ID}_1.fastq.gz" "${OUTDIR}/${Descriptive_Title}_1.fastq.gz"
done < "${OUTDIR}/SRAforEaf3Ash1Set7WTRco1withMeta.txt"

#downloading reference genome
#curl -s ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/Neurospora_crassa/latest_assembly_versions/GCA_000182925.2_NC12/GCA_000182925.2_NC12_genomic.fna.gz  | gunzip -c > ${OUTDIR}/NC12_genome.fna
#curl -s ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/182/925/GCA_000182925.2_NC12/GCA_000182925.2_NC12_genomic.gtf.gz | gunzip -c > ${OUTDIR}/NC12_annotation.gtf

#DOWNLOADING SRA IswiPaperFilz  aesFromMasayuki
#prefetch -O ${OUTDIR} SRR8730382 SRR8730383 SRR8730380 SRR8730381 SRR8730378 SRR8730379 SRR8730376 SRR8730377

# fastq-dump --split-files --gzip	${OUTDIR}/SRR8730382
# fastq-dump --split-files --gzip	${OUTDIR}/SRR8730383
# fastq-dump --split-files --gzip	${OUTDIR}/SRR8730380
# fastq-dump --split-files --gzip	${OUTDIR}/SRR8730381
# fastq-dump --split-files --gzip	${OUTDIR}/SRR8730378
# fastq-dump --split-files --gzip	${OUTDIR}/SRR8730379
# fastq-dump --split-files --gzip	${OUTDIR}/SRR8730376
# fastq-dump --split-files --gzip	${OUTDIR}/SRR8730377
#
#
# ##rename with descriptive titles
# mv	${OUTDIR}/SRR8730382_1.fastq.gzip	${OUTDIR}/EPR-1GFP_BAH_N7752_1.fastq.gzip
# mv	${OUTDIR}/SRR8730383_1.fastq.gzip	${OUTDIR}/EPR-1GFP_BAH_N7753_1.fastq.gzip
# mv	${OUTDIR}/SRR8730380_1.fastq.gzip	${OUTDIR}/EPR-1GFP_PHD_N7754_1.fastq.gzip
# mv	${OUTDIR}/SRR8730381_1.fastq.gzip	${OUTDIR}/EPR-1GFP_PHD_N7755_1.fastq.gzip
# mv	${OUTDIR}/SRR8730378_1.fastq.gzip	${OUTDIR}/EPR-1GFP_N7689_1.fastq.gzip
# mv	${OUTDIR}/SRR8730379_1.fastq.gzip	${OUTDIR}/EPR-1GFP_N7690_1.fastq.gzip
# mv	${OUTDIR}/SRR8730376_1.fastq.gzip	${OUTDIR}/EPR-1_H3K27me2_3_N7479_1.fastq.gzip
# mv	${OUTDIR}/SRR8730377_1.fastq.gzip	${OUTDIR}/EPR-1_H3K27me2_3_N7480_1.fastq.gzip


##convert to fastq
# fastq-dump --split-files --gzip ${OUTDIR}/EPR-1GFP_BAH_N7752 -O ${OUTDIR}
# fastq-dump --split-files --gzip ${OUTDIR}/EPR-1GFP_BAH_N7753 -O ${OUTDIR}
# fastq-dump --split-files --gzip ${OUTDIR}/EPR-1GFP_PHD_N7754 -O ${OUTDIR}
# fastq-dump --split-files --gzip ${OUTDIR}/EPR-1GFP_PHD_N7755.sra -O ${OUTDIR}
# fastq-dump --split-files --gzip ${OUTDIR}/EPR-1GFP_N7689.sra -O ${OUTDIR}
# fastq-dump --split-files --gzip ${OUTDIR}/EPR-1GFP_N7690.sra -O ${OUTDIR}
# fastq-dump --split-files --gzip ${OUTDIR}/EPR-1_H3K27me2_3_N7479.sra -O ${OUTDIR}
# fastq-dump --split-files --gzip ${OUTDIR}/EPR-1_H3K27me2_3_N7480.sra -O ${OUTDIR}

# #process reads using trimGalore

source config.txt

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
#
FILES="${OUTDIR}/TrimmedReads/*R1_001_val_1\.fq\.gz" #Don't forget the *
#

#
#Iterate over the files
for f in $FILES
do
#
# 	#Examples to Get Different parts of the file name
# 		#See here for details: http://tldp.org/LDP/abs/html/refcards.html#AEN22664
		#${string//substring/replacement}
# 		#dir=${f%/*}

	file=${f##*/}
	#remove ending from file name to create shorter names for bam files and other downstream output
	name=${file/%_S[1-12]*_L001_R1_001_val_1.fq.gz/}

#
# 	# File Vars
# 	#use sed to get the name of the second read matching the input file
	read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R2_001_val_2\.fq\.gz/g')
	#variable for naming bam file
 	bam="${OUTDIR}/SortedBamFiles/${name}.bam"
	#variable name for bigwig output
	bigwig="${OUTDIR}/BigWigs/${name}"
	#QualityBam="${OUTDIR}/SortedBamFiles/${name}_Q30.bam"
#

#
bwa mem -M -v 3 -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIR/SortedBamFiles/tempReps -o "$bam" -
samtools index "$bam"

#samtools view -b -q 30 $bam > "$QualityBam"
#samtools index "$QualityBam"

############################
# # #deeptools

ml deepTools
# #use these parameters for ChIP data
bamCoverage -p $THREADS $MNase -bs $BIN --normalizeUsing BPM --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}${MN}.bw"


done

for file_path in "${OUTDIR}/BigWigs"/*.bw; do
  # Get the base name of the file
  BW_name=$(basename "${file_path}" .bw)

  # Remove the file extension to create the sample ID
  BW_id="${BW_name%.*}"
  # Replace special characters with underscores in the sample ID
 BW_id=${BW_id//[^a-zA-Z0-9]/_}

 # Limit the length of the sample ID to avoid long filenames
 BW_id=${BW_id:0:50}
 # Compute matrix for the reference-point TSS
 computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S "${file_path}" -R "/scratch/ry00555/neurospora.bed" --skipZeros -o "${OUTDIR}/Matrices/matrix_${BW_id}.gz"

  # Compute matrix for the reference-point TSS with specific regions (e.g., PRC2 genes)
  computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S "${file_path}" -R "/scratch/ry00555/heatmapPRC2genes.bed" --skipZeros -o "${OUTDIR}/Matrices/PRC2genes_matrix_${BW_id}.gz"

  # Preprocess the matrix files to replace nan values with zeros
  zcat "${OUTDIR}/Matrices/matrix_${BW_id}.gz" | awk '{for (i=1; i<=NF; i++) if ($i == "nan") $i=0; print}' | gzip > "${OUTDIR}/Matrices/matrix_${BW_id}_processed.gz"
  zcat "${OUTDIR}/Matrices/PRC2genes_matrix_${BW_id}.gz" | awk '{for (i=1; i<=NF; i++) if ($i == "nan") $i=0; print}' | gzip > "${OUTDIR}/Matrices/PRC2genes_matrix_${BW_id}_processed.gz"

  # Plot heatmaps for the reference-point TSS
 plotHeatmap --matrixFile "${OUTDIR}/Matrices/matrix_${BW_id}.gz" --outFileName "${OUTDIR}/Heatmaps/${BW_id}_hclust.png" --samplesLabel "${BW_name}" --hclust 1 --colorMap Reds --sortRegions descend --missingDataColor white --sortUsingSamples 1

  # Plot heatmaps for the reference-point TSS with specific regions (e.g., PRC2 genes)
 plotHeatmap --matrixFile "${OUTDIR}/Matrices/PRC2genes_matrix_${BW_id}.gz" --outFileName "${OUTDIR}/Heatmaps/${BW_id}_hclust_PRC2genes.png" --samplesLabel "${BW_name}" --hclust 1 --colorMap Reds
done
