#!/bin/bash
#SBATCH --job-name=Run148ChIP
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100gb
#SBATCH --time=72:00:00
#SBATCH --output=../MapCutAndRun148.%j.out
#SBATCH --error=../MapCutAndRun148.%j.err

cd $SLURM_SUBMIT_DIR

#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )

source config.txt
OUTDIR="/scratch/ry00555/Run148"

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
    mkdir -p "${OUTDIR}/TrimmedReads"
    mkdir -p "${OUTDIR}/BigWigs"
   mkdir -p "$OUTDIR/HomerTagDirectories"
#   mkdir -p "$OUTDIR/TdfFiles"
  mkdir -p "$OUTDIR/SortedBamFiles"

fi

TAGDIR="${OUTDIR}/HomerTagDirectories"
BAMDIR="${OUTDIR}/SortedBamFiles"
BEDDIR="${OUTDIR}/Beds"
#
# # #process reads using trimGalore
ml Trim_Galore
trim_galore --illumina --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${FASTQ}/*fastq\.gz
