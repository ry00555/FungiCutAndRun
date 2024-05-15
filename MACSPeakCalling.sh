#!/bin/bash
#SBATCH --job-name=CUTandRun137
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=90gb
#SBATCH --time=48:00:00
#SBATCH --output=../MacsPeakCalling.%j.out
#SBATCH --error=../MacsPeakCalling.%j.err
cd $SLURM_SUBMIT_DIR

# 137-10_CUTANDRUN_rtt109_H3K27me3_Rep1_S10.sorted.bam
# 137-19_CUTANDRUN_ncu00423_H3K27me3_Rep1_S19.sorted.bam
# 137-2_CUTANDRUN_WT_H3K27me3_Rep1_S2.sorted.bam
# 137-10_CUTANDRUN_rtt109_H3K27me3_Rep1_S10_Ecoli.sorted.bam
# 137-19_CUTANDRUN_ncu00423_H3K27me3_Rep1_S19_Ecoli.sorted.bam
# 137-2_CUTANDRUN_WT_H3K27me3_Rep1_S2_Ecoli.sorted.bam
# 137-11_CUTANDRUN_rtt109_H3K36me3_Rep1_S11.sorted.bam
# 137-1_CUTANDRUN_WT_IgG_Rep1_S1.sorted.bam
# 137-3_CUTANDRUN_WT_H3K36me3_Rep1_S3.sorted.bam
# 137-11_CUTANDRUN_rtt109_H3K36me3_Rep1_S11_Ecoli.sorted.bam
# 137-1_CUTANDRUN_WT_IgG_Rep1_S1_Ecoli.sorted.bam
# 137-3_CUTANDRUN_WT_H3K36me3_Rep1_S3_Ecoli.sorted.bam
# 137-12_CUTANDRUN_ncu06787_IgG_Rep1_S12.sorted.bam
# 137-20_CUTANDRUN_ncu00423_H3K36me3_Rep1_S20.sorted.bam
# 137-4_CUTANDRUN_WT_IgG_Rep1_S4.sorted.bam
# 137-12_CUTANDRUN_ncu06787_IgG_Rep1_S12_Ecoli.sorted.bam
# 137-20_CUTANDRUN_ncu00423_H3K36me3_Rep1_S20_Ecoli.sorted.bam
# 137-4_CUTANDRUN_WT_IgG_Rep1_S4_Ecoli.sorted.bam
# 137-13_CUTANDRUN_ncu06787_H3K27me3_Rep1_S13.sorted.bam
# 137-21_CUTANDRUN_WT_Input_Rep1_S21.sorted.bam
# 137-5_CUTANDRUN_WT_IgG_Rep1_S5.sorted.bam
# 137-13_CUTANDRUN_ncu06787_H3K27me3_Rep1_S13_Ecoli.sorted.bam
# 137-21_CUTANDRUN_WT_Input_Rep1_S21_Ecoli.sorted.bam
# 137-5_CUTANDRUN_WT_IgG_Rep1_S5_Ecoli.sorted.bam
# 137-14_CUTANDRUN_ncu06787_H3K36me3_Rep1_S14.sorted.bam
# 137-22_CUTANDRUN_WT_H3K27me3_Rep1_S22.sorted.bam
# 137-6_CUTANDRUN_set-7_IgG_Rep1_S6.sorted.bam
# 137-14_CUTANDRUN_ncu06787_H3K36me3_Rep1_S14_Ecoli.sorted.bam
# 137-22_CUTANDRUN_WT_H3K27me3_Rep1_S22_Ecoli.sorted.bam
# 137-6_CUTANDRUN_set-7_IgG_Rep1_S6_Ecoli.sorted.bam
# 137-15_CUTANDRUN_ncu06788_IgG_Rep1_S15.sorted.bam
# 137-23_CUTANDRUN_WT_H3K36me3_Rep1_S23.sorted.bam
# 137-7_CUTANDRUN_set-7_H3K27me3_Rep1_S7.sorted.bam
# 137-15_CUTANDRUN_ncu06788_IgG_Rep1_S15_Ecoli.sorted.bam
# 137-23_CUTANDRUN_WT_H3K36me3_Rep1_S23_Ecoli.sorted.bam
# 137-7_CUTANDRUN_set-7_H3K27me3_Rep1_S7_Ecoli.sorted.bam
# 137-16_CUTANDRUN_ncu06788_H3K27me3_Rep1_S16.sorted.bam
# 137-24_CUTANDRUN_set-7_Input_Rep1_S24.sorted.bam
# 137-8_CUTANDRUN_set-7_H3K36me3_Rep1_S8.sorted.bam
# 137-16_CUTANDRUN_ncu06788_H3K27me3_Rep1_S16_Ecoli.sorted.bam
# 137-24_CUTANDRUN_set-7_Input_Rep1_S24_Ecoli.sorted.bam
# 137-8_CUTANDRUN_set-7_H3K36me3_Rep1_S8_Ecoli.sorted.bam
# 137-17_CUTANDRUN_ncu06788_H3K36me3_Rep1_S17.sorted.bam
# 137-25_CUTANDRUN_set-7_H3K27me3_Rep1_S25.sorted.bam
# 137-9_CUTANDRUN_rtt109_IgG_Rep1_S9.sorted.bam
# 137-17_CUTANDRUN_ncu06788_H3K36me3_Rep1_S17_Ecoli.sorted.bam
# 137-25_CUTANDRUN_set-7_H3K27me3_Rep1_S25_Ecoli.sorted.bam
# 137-9_CUTANDRUN_rtt109_IgG_Rep1_S9_Ecoli.sorted.bam
# 137-18_CUTANDRUN_ncu00423_IgG_Rep1_S18.sorted.bam
# 137-26_CUTANDRUN_set-7_H3K36me3_Rep1_S26.sorted.bam           _Ecoli.sorted.bam
# 137-18_CUTANDRUN_ncu00423_IgG_Rep1_S18_Ecoli.sorted.bam
# 137-26_CUTANDRUN_set-7_H3K36me3_Rep1_S26_Ecoli.sorted.bam

module load MACS3
#command line
#macs3 callpeak -t 137-11_CUTANDRUN_rtt109_H3K36me3_Rep1_S11_Ecoli.sorted.bam -f BAMPE -n 137-11_CUTANDRUN_rtt109_H3K36me3_Rep1_S11_Ecoli -c 137-9_CUTANDRUN_rtt109_IgG_Rep1_S9_Ecoli.sorted.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir /scratch/ry00555/OutputRun137/CutandRun/MACSPeaks --min-length 800 --max-gap 500

for infile in $OUTDIR/SortedBamFiles/TagDirectories/*rtt109*.sorted.bam
 do
  base=$(basename ${infile} .sorted.bam)
  base2=$(basename ${infile} .EColi.sorted.bam)
macs3 callpeak -t $infile -f BAMPE -n $base -c 137-9_CUTANDRUN_rtt109_IgG_Rep1_S9.sorted.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir /scratch/ry00555/OutputRun137/CutandRun/MACSPeaks --min-length 800 --max-gap 500
macs3 callpeak -t $infile -f BAMPE -n $base2 -c 137-9_CUTANDRUN_rtt109_IgG_Rep1_S9_Ecoli.sorted.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir /scratch/ry00555/OutputRun137/CutandRun/MACSPeaks --min-length 800 --max-gap 500
done
