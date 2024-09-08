#!/bin/bash
#SBATCH --job-name=Run142ChIP
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=../MapCutAndRun142.%j.out
#SBATCH --error=../MapCutAndRun142.%j.err

cd $SLURM_SUBMIT_DIR

#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )

source config.txt

OUTDIR="/scratch/ry00555/OutputRun142"


  # mkdir "${OUTDIR}/TrimmedReads"
  # mkdir "${OUTDIR}/BigWigs"
  # mkdir "$OUTDIR/HomerTagDirectories"
 #mkdir "$OUTDIR/TdfFiles"
 #mkdir "$OUTDIR/SortedBamFiles"
#
#
#TAGDIR="${OUTDIR}/HomerTagDirectories"
#BAMDIR="${OUTDIR}/SortedBamFiles"
#BEDDIR="${OUTDIR}/Beds"
#
# # #process reads using trimGalore
# #
 ml Trim_Galore/0.6.7-GCCcore-11.2.0
 trim_galore --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${FASTQ}/142-9*fastq\.gz
# #
