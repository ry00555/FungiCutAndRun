#!/bin/bash
#SBATCH --job-name=CUTandRun137
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=90gb
#SBATCH --time=48:00:00
#SBATCH --output=../MapQual137.%j.out
#SBATCH --error=../MapQual137.%j.err
cd $SLURM_SUBMIT_DIR
OUTDIR="/scratch/ry00555/OutputRun137/MapQual_CutNRun"
NCTRIMMED="/scratch/ry00555/OutputRun137/MapQual_CutNRun/TrimmedReads/"

FASTQ="/scratch/ry00555/OutputRun137/CutandRun"

ml STAR
 for file in $FASTQ/*fastq\.gz;
 do
   if [[ $prefix ]]; then
         base=$(basename ${first} _R1_001_val_1.fq.gz)
         sh /home/ry00555/Research/FungiCutAndRun/CUTandRUNAnalysis/PE_trim_and_star_RY.sh -o $OUTDIR -n $base -m ONE $first $file
         prefix=
     else
         first=$file
         prefix=${file%%_*}
     fi
 done
#
# #aligning to ecoli genome
# ml STAR
#  for file in $FASTQ/*fastq\.gz;
#  do
#    if [[ $prefix ]]; then
#          base=$(basename ${first} _R1_001_val_1.fq.gz)
#          sh /home/ry00555/Research/FungiCutAndRun/CUTandRUNAnalysis/PE_trim_and_star_e_coli_RY.sh -o $OUTDIR/Ecoli_Aligned -n $base -m one $first $file
#          prefix=
#      else
#          first=$file
#          prefix=${file%%_*}
#      fi
#  done
