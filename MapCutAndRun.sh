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

OUTDIR="/scratch/ry00555/OutputRun135"


 # mkdir "${OUTDIR}/TrimmedReads"
 # mkdir "${OUTDIR}/BigWigs"
 # mkdir "${OUTDIR}/Peaks"
 # mkdir "$OUTDIR/HomerTagDirectories"
 # mkdir "$OUTDIR/TdfFiles"
 # mkdir "$OUTDIR/SortedBamFiles"
#
#
PEAKDIR="${OUTDIR}/Peaks"
TAGDIR="${OUTDIR}/HomerTagDirectories"
BAMDIR="${OUTDIR}/SortedBamFiles"
BEDDIR="${OUTDIR}/Beds"
#
# # #process reads using trimGalore
# #
  ml Trim_Galore
  trim_galore --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${FASTQ}/*fastq\.gz
# #
 FILES="${OUTDIR}/TrimmedReads/*R1_001_val_1\.fq\.gz" #Don't forget the *
# #
#
# #
# #Iterate over the files
for f in $FILES
 do
# #
# # 	#Examples to Get Different parts of the file name
# # 		#See here for details: http://tldp.org/LDP/abs/html/refcards.html#AEN22664
# 		#${string//substring/replacement}
# # 		#dir=${f%/*}
#
 	file=${f##*/}
# 	#remove ending from file name to create shorter names for bam files and other downstream output
 	name=${file/%_S[1-12]*_R1_001_val_1.fq.gz/}
#
# #
# # 	# File Vars
# # 	#use sed to get the name of the second read matching the input file
 	read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R2_001_val_2\.fq\.gz/g')
# 	#variable for naming bam file
bam="${OUTDIR}/SortedBamFiles/${name}.bam"
# 	#variable name for bigwig output
 	bigwig="${OUTDIR}/BigWigs/${name}"
# 	#QualityBam="${OUTDIR}/SortedBamFiles/${name}_Q30.bam"
# #
#
ml SAMtools
ml BWA
# #
 bwa mem -M -v 3 -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIR/SortedBamFiles/tempReps -o "$bam" -
 samtools index "$bam"
#
# #samtools view -b -q 30 $bam > "$QualityBam"
# #samtools index "$QualityBam"
#
# ############################
# # # #deeptools
#
 ml deepTools
# #Plot all reads
 bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"
#
 #plot mononucleosomes
 bamCoverage -p $THREADS --MNase -bs 1 --normalizeUsing BPM --smoothLength 25 -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}_MNase.bw"
 done

ml Homer
ml Perl
ml SAMtools
ml BEDTools
for bam_file in "${BAMDIR}"/*.bam; do
  # Get the sample ID from the BAM file name
  sample_id=$(basename "${bam_file}" .bam)
  # Remove everything after "Rep_1" in the sample ID
  sample_id="${sample_id%%_Rep_1*}"


#makeTagDirectory "${TAGDIR}/${sample_id}" "${bam_file}"
#
#   # Call peaks
#
findPeaks "${TAGDIR}/${sample_id}" -style histone -region -size 150 -minDist 530 -o "${PEAKDIR}/${sample_id}_peaks.txt"
#
done
#changing peak txt files to bed files to input into chipr
ml ChIP-R
 for infile in ${PEAKDIR}/${sample_id}_peaks.txt
do
  base=$(basename ${infile} .txt)
  sed '/^#/d' $infile | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $8 "\t" $5 "\t" $6 "\t" $12 "\t" "-1"}' | sed 's/\.000000//g' > ${PEAKDIR}/${base}.peaks.bed
done

ml Homer
ml Perl
##annotating peak files with masked reference (use HOMER module)
#curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/925/GCF_000182925.2_NC12/GCF_000182925.2_NC12_genomic.gtf.gz | gunzip -c > Ncrassa_refann.gtf
 annotatePeaks.pl ${PEAKDIR}/${base}.peaks.bed -gtf scratch/ry00555/Ncrassa_refann.gtf > ${PEAKDIR}/${base}_ann.txt

#now filtering for only peaks that are w/i 1000bps of their annotation:
 for infile in ${PEAKDIR}/${base}_ann.txt
 do
   base=$(basename ${infile} _masked_ann.txt)
   awk -F'\t' 'sqrt($10*$10) <=1000' $infile > ${PEAKDIR}/${base}.1000bp_ann.txt
 done
