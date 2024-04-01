#!/bin/bash
#SBATCH --job-name=Run132ChIP
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=../MapCutAndRun132.%j.out
#SBATCH --error=../MapCutAndRun132.%j.err

cd $SLURM_SUBMIT_DIR

#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )

source config.txt

OUTDIR=/scratch/ry00555/OutputRun138
mkdir "${OUTDIR}/SortedBamFiles"
mkdir "${OUTDIR}/BigWigs"
mkdir "${OUTDIR}/Peaks"
mkdir "$OUTDIR/HomerTagDirectories"
mkdir "$OUTDIR/TdfFiles"

PEAKDIR="${OUTDIR}/Peaks"
TAGDIR="${OUTDIR}/HomerTagDirectories"
BAMDIR="${OUTDIR}/SortedBamFiles"
BEDDIR="${OUTDIR}/Beds"

# #process reads using trimGalore
#
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
	name=${file/%_S[1-12]*_R1_001_val_1.fq.gz/}

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

ml SAMtools
ml BWA
#
bwa mem -M -v 3 -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIR/SortedBamFiles/tempReps -o "$bam" -
samtools index "$bam"

#samtools view -b -q 30 $bam > "$QualityBam"
#samtools index "$QualityBam"

############################
# # #deeptools

ml deepTools
#Plot all reads
bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"

#plot mononucleosomes
bamCoverage -p $THREADS --MNase -bs 1 --normalizeUsing BPM --smoothLength 25 -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}_MNase.bw"
done

for bam_file in "${BAMDIR}"/*.bam; do
  # Get the sample ID from the BAM file name
  sample_id=$(basename "${bam_file}" .bam)
  # Remove everything after "Rep_1" in the sample ID
  sample_id="${sample_id%%_Rep_1*}"

ml Homer
ml Perl

ml BEDTools

makeTagDirectory "${TAGDIR}/${sample_id}" "${bam_file}"
#
#   # Call peaks
#
  findPeaks "${TAGDIR}/${sample_id}" -style histone -region -size 150 -minDist 530 -o "${TAGDIR}/${sample_id}_peaks.txt"
#
#  # Calculate coverage for 1kb windows
 bedtools coverage -a "/scratch/ry00555/${Genome}_1kbWindows.bed" -b "${bam_file}" > "${BEDDIR}/${sample_id}_coverage.bed"
done
