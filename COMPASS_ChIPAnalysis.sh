#!/bin/bash
#SBATCH --job-name=COMPASSAnalysis
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=../MapCutAndRun145_142.%j.out
#SBATCH --error=../MapCutAndRun145_142.%j.err

cd $SLURM_SUBMIT_DIR
source config.txt

OUTDIR="/scratch/ry00555/RNASeqPaper/COMPASS/ChIP"

if [ ! -d $OUTDIR ]
then
mkdir -p $OUTDIR
mkdir -p "${OUTDIR}/TrimmedReads"
mkdir -p "${OUTDIR}/BigWigs"
mkdir -p "$OUTDIR/HomerTagDirectories"
mkdir -p "$OUTDIR/TdfFiles"
mkdir -p "$OUTDIR/SortedBamFiles"

fi

TAGDIR="${OUTDIR}/HomerTagDirectories"
BAMDIR="${OUTDIR}/SortedBamFiles"
BEDDIR="${OUTDIR}/Beds"
PEAKDIR="${OUTDIR}/MACSPeaks"

# FILES="${OUTDIR}/TrimmedReads/*_R1_001_val_1\.fq\.gz"
#
# for f in $FILES
# do
# file=${f##*/}
# name=${file/%_S[1-990]*_R1_001_val_1.fq.gz/}
# read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R2_001_val_2\.fq\.gz/g')
#
# bam="${OUTDIR}/SortedBamFiles/${name}.bam"
# variable name for bigwig output
#  	bigwig="${OUTDIR}/BigWigs/${name}"
# QualityBam="${OUTDIR}/SortedBamFiles/${name}_Q30.bam"
#
# ml SAMtools/1.16.1-GCC-11.3.0
# ml BWA/0.7.17-GCCcore-11.3.0
#
# bwa mem -M -v 3 -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIR/SortedBamFiles/tempReps -o "$bam" -
# samtools index "$bam"
#
# samtools view -b -q 30 $bam > "$QualityBam"
# samtools index "$QualityBam"
#
# ml deepTools/3.5.2-foss-2022a
# bamCoverage -p $THREADS -bs 10 --normalizeUsing BPM --minMappingQuality 10 --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_10.smooth_${SMOOTH}Bulk.bw"
#
# done


#142-94_ChIP_set7_H3K27me3__Q30.bam 142-93_ChIP_set7_Input__Q30.bam
#142-106_ChIP_swd1_H3K27me3__Q30.bam       142-105_ChIP_swd1_Input__Q30.bam
#142-10_ChIP_WT_H3K27me3_Rep3_Q30.bam 142-75_ChIP_WT_Input__Q30.bam

#142-115_ChIP_set2_H3K27me3__Q30.bam
#142-118_ChIP_set2_H3K27me3__Q30.bam
#142-121_ChIP_set1_H3K27me3__Q30.bam       142-123_ChIP_set1_Input__Q30.bam
#142-124_ChIP_set1_H3K27me3__Q30.bam         142-123_ChIP_set1_Input__Q30.bam
#142-127_ChIP_sgr9_H3K27me3__Q30.bam       142-126_ChIP_sgr9_Input__Q30.bam
#142-75_ChIP_WT_Input__Q30.bam


mkdir $OUTDIR/MACSPeaks


module load MACS3/3.0.0b1-foss-2022a-Python-3.10.4
 #command line
macs3 callpeak -t $BAMDIR/142-94_ChIP_set7_H3K27me3__Q30.bam -f BAMPE -n 142-94_ChIP_set7_H3K27me3 -c $BAMDIR/142-93_ChIP_set7_Input__Q30.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
macs3 callpeak -t $BAMDIR/142-106_ChIP_swd1_H3K27me3__Q30.bam -f BAMPE -n 142-106_ChIP_swd1_H3K27me3 -c $BAMDIR/142-105_ChIP_swd1_Input__Q30.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
macs3 callpeak -t $BAMDIR/142-121_ChIP_set1_H3K27me3__Q30.bam -f BAMPE -n 142-121_ChIP_set1_H3K27me3 -c $BAMDIR/142-123_ChIP_set1_Input__Q30.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
macs3 callpeak -t $BAMDIR/142-124_ChIP_set1_H3K27me3__Q30.bam -f BAMPE -n 142-124_ChIP_set1_H3K27me3 -c $BAMDIR/142-123_ChIP_set1_Input__Q30.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
macs3 callpeak -t $BAMDIR/142-127_ChIP_sgr9_H3K27me3__Q30.bam  -f BAMPE -n 142-127_ChIP_sgr9_H3K27me3 -c $BAMDIR/142-126_ChIP_sgr9_Input__Q30.bam  --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
macs3 callpeak -t $BAMDIR/142-10_ChIP_WT_H3K27me3_Rep3_Q30.bam  -f BAMPE -n 142-10_ChIP_WT_H3K27me3 -c $BAMDIR/142-75_ChIP_WT_Input__Q30.bam  --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500

macs3 callpeak -t $BAMDIR/142-115_ChIP_set2_H3K27me3__Q30.bam -f BAMPE -n 142-115_ChIP_set2_H3K27me3  --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
macs3 callpeak -t $BAMDIR/142-118_ChIP_set2_H3K27me3__Q30.bam -f BAMPE -n 142-118_ChIP_set2_H3K27me3  --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500

mkdir ${OUTDIR}/HomerPeaks
HOMERPEAKSDIR="${OUTDIR}/HomerPeaks"
ml Homer
ml Perl
ml SAMtools
ml BEDTools

for bam_file in "${BAMDIR}"/*__Q30.bam; do
  sample_id=$(basename "${bam_file}" __Q30.bam)
makeTagDirectory "${TAGDIR}/${sample_id}" "${bam_file}"
done
