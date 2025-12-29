#!/bin/bash
#SBATCH --job-name=Run153ChIP
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100gb
#SBATCH --time=20:00:00
#SBATCH --output=../MapCutAndRun153.%j.out
#SBATCH --error=../MapCutAndRun153.%j.err

cd $SLURM_SUBMIT_DIR

#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )

source config.txt
OUTDIR="/lustre2/scratch/ry00555/Run147/TrimmedReads/RemapTemp/Missed"

# if output directory doesn't exist, create it
  mkdir -p $OUTDIR
   mkdir -p "${OUTDIR}/TrimmedReads"
      mkdir -p "${OUTDIR}/BigWigs"
  mkdir -p "$OUTDIR/HomerTagDirectories"
#    mkdir -p "$OUTDIR/TdfFiles"
 mkdir -p "$OUTDIR/SortedBamFiles"

 TAGDIR="${OUTDIR}/HomerTagDirectories"
 BAMDIR="${OUTDIR}/SortedBamFiles"
 BEDDIR="${OUTDIR}/Beds"
#   process reads using trimGalore
module load Trim_Galore
trim_galore --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${FASTQ}/*fastq\.gz

 #FILES="${OUTDIR}/TrimmedReads/*_L006_R1_001_val_1\.fq\.gz"

# FILES="${OUTDIR}/TrimmedReads/*_R1_001_val_1\.fq\.gz"
#  Iterate over the files
# for f in $FILES
#do

#   	Examples to Get Different parts of the file name
#   		See here for details: http://tldp.org/LDP/abs/html/refcards.htmlAEN22664
#  		${string//substring/replacement}
#   		dir=${f%/*}
#
#file=${f##*/}
#  	remove ending from file name to create shorter names for bam files and other downstream output
# name=${file/%_S[1-150]*_L001_R1_001_val_1.fq.gz/}
# name=${file/%_S[1-990]*_L002_R1_001_val_1.fq.gz/}
 #name=${file/%_S[1-190]*R1_001_val_1.fq.gz/}
#
#
#  	 File Vars
#  	use sed to get the name of the second read matching the input file
#read2=$(echo "$f" | sed 's/_L001_R1_001_val_1\.fq\.gz/_L001_R2_001_val_2\.fq\.gz/g')
#
#   	 File Vars
#   	use sed to get the name of the second read matching the input file
#read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R2_001_val_2\.fq\.gz/g')
#  	variable for naming bam file
#bam="${OUTDIR}/SortedBamFiles/${name}.bam"
#  	variable name for bigwig output
#bigwig="${OUTDIR}/BigWigs/${name}"
#QualityBam="${OUTDIR}/SortedBamFiles/${name}_Q30.bam"


#ml SAMtools/1.16.1-GCC-11.3.0
#ml BWA/0.7.17-GCCcore-11.3.0

#  bwa mem -M -v 3 -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIR/SortedBamFiles/tempReps -o "$bam" -
#samtools index "$bam"
#
 #samtools view -b -q 30 $bam > "$QualityBam"
# samtools index "$QualityBam"
#
#
#    deeptools
#
#ml deepTools
#  Plot all reads
# bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --minMappingQuality 10 --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"
#
 #bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --minMappingQuality 10 --smoothLength $SMOOTH -of bigwig -b "$QualityBam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}_Q30.bw"
 #done
mkdir $OUTDIR/MACSPeaks
PEAKDIR="${OUTDIR}/MACSPeaks"

ml MACS3
# command line
#macs3 callpeak -t 137-11_CUTANDRUN_rtt109_H3K36me3_Rep1_S11_Ecoli.sorted.bam -f BAMPE -n 137-11_CUTANDRUN_rtt109_H3K36me3_Rep1_S11_Ecoli -c 137-9_CUTANDRUN_rtt109_IgG_Rep1_S9_Ecoli.sorted.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir /scratch/ry00555/OutputRun137/CutandRun/MACSPeaks --min-length 800 --max-gap 500

   for infile in $BAMDIR/*.bam
  do
     base=$(basename ${infile} .bam)
#    Input=$BAMDIR/ ${infile} Input_Q30.bam
  macs3 callpeak -t $infile -f BAMPE -n $base --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500 #-c $Input
 done
#
# HOMERPEAKSDIR="${OUTDIR}/HomerPeaks"
#    ml Homer
#    ml Perl
#   ml SAMtools
#     ml BEDTools
#     for bam_file in "${BAMDIR}"/*_Q30.bam; do
# #       Get the sample ID from the BAM file name
#   sample_id=$(basename "${bam_file}" _Q30.bam)
# #       Remove everything after "Rep_1" in the sample ID
#   HOMERINPUT="${TAGDIR}/${sample_id}_Input*"
# #
# #
#   makeTagDirectory "${TAGDIR}/${sample_id}" "${bam_file}"
#
# #       Call peaks
#
  # findPeaks "${TAGDIR}/${sample_id}" -style histone -region -size 150 -minDist 530 -o "${HOMERPEAKSDIR}/${sample_id}_Homerpeaks.txt" -i $HOMERINPUT
#
  #  done
#  changing peak txt files to bed files to input into chipr
#ml ChIP-R

#    for infile in ${HOMERPEAKSDIR}/*_Homerpeaks.txt
#   do
#     base=$(basename ${infile} _Homerpeaks.txt)
#     sed '/^/d' $infile | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $8 "\t" $5 "\t" $6 "\t" $12 "\t" "-1"}' | sed 's/\.000000//g' > ${HOMERPEAKSDIR}/${base}.peaks.bed
#
#   done
#
# for infile in ${HOMERPEAKSDIR}/*.peaks.bed
# do
# base=$(basename ${infile} .peaks.bed)
# annotatePeaks.pl ${HOMERPEAKSDIR}/${base}.peaks.bed "/home/ry00555/Research/Genomes/GCA_00182925.2plusHphplusBarplusTetO_his3masked.fna" -gff "/scratch/ry00555/GCA_000182925.2_NC12_genomic_WithExtras.gff" > ${HOMERPEAKSDIR}/${base}_ann.txt
# done
# #  annotating peak files with masked reference (use HOMER module)
# #  curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/925/GCF_000182925.2_NC12/GCF_000182925.2_NC12_genomic.gtf.gz | gunzip -c > Ncrassa_refann.gtf
#   annotatePeaks.pl ${HOMERPEAKSDIR}/${base}.peaks.bed -gtf /scratch/ry00555/Ncrassa_refann.gtf > ${HOMERPEAKSDIR}/${base}_ann.txt
#
# #  now filtering for only peaks that are w/i 1000bps of their annotation:
# for infile in ${HOMERPEAKSDIR}/*_ann.txt
#    do
#      base=$(basename ${infile} _ann.txt)
#      awk -F'\t' 'sqrt($10*$10) <=1000' $infile > ${HOMERPEAKSDIR}/${base}.1000bp_ann.txt
#    done

# ml deepTools
#  bamCompare -b1 treatment.bam -b2 control.bam -o log2ratio.bw
#  bamCompare -p max -b1 "${BAMDIR}"/147-3_ChIP_S2_WT_H3K27me3_Rep5_Nc_24hrVMMON_S3_L001_R1_001_val_1.fq.gz_Q30.bam -b2 "${BAMDIR}"/147-1_ChIP_S2_WT_Input_Rep5_Nc_24hrVMMON_Q30.bam -of bigwig -o "${BAMDIR}"/147-WT-H3K27me3-ReadCountInputNorm.bw
#
#   bamCompare -p max -b1 "${BAMDIR}"/135-83_ChIP_rtt109_H3K27me3_Rep2_S79_L001_R1_001_val_1.fq.gz_Q30.bam -b2 "${BAMDIR}"/135-80_ChIP_rtt109_Input_Rep2_S76_L001_R1_001_val_1.fq.gz_Q30.bam -of bigwig -o "${BAMDIR}"/135-RTT109-H3K27me3-ReadCountInputNorm.bw
#
#     bamCompare -p max -b1 "${BAMDIR}"/142-70_ChIP_naf2_H3K27me3__S70_L007_R1_001_val_1.fq.gz_Q30.bam -b2 "${BAMDIR}"/142-69_ChIP_naf2_Input__S69_L007_R1_001_val_1.fq.gz_Q30.bam -of bigwig -o "${BAMDIR}"/142-NAF2-H3K27me3-ReadCountInputNorm.bw
#
# bamCompare -p max -b1 "${BAMDIR}"/142-88_ChIP_rtt1093xFLAG_H3K27me3__S88_L007_R1_001_val_1.fq.gz_Q30.bam -b2  "${BAMDIR}"/142-87_ChIP_rtt1093xFLAG_Input__S87_L007_R1_001_val_1.fq.gz_Q30.bam -of bigwig -o "${BAMDIR}"/142-RTT1093xFLAG-H3K27me3-ReadCountInputNorm.bw
# bamCompare -p max -b1 "${BAMDIR}"/147-43_ChIP_H3K56R13P4_H3K56R13_H3K27me3_Rep1_Nc_24hrVMMON_S43_L001_R1_001_val_1.fq.gz_Q30.bam -b2  "${BAMDIR}"/147-41_ChIP_H3K56R13P4_H3K56R13_Input_Rep1_Nc_24hrVMMON_S41_L001_R1_001_val_1.fq.gz_Q30.bam -of bigwig -o "${BAMDIR}"/147-H3K56R13-H3K27me3-ReadCountInputNorm.bw
