#!/bin/bash
#SBATCH --job-name=MacsPeakCalling
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


OUTDIR="/scratch/ry00555/ChIPQC/"
module load MACS3
#command line
macs3 callpeak -t 111_13_ChIP_K56R_40_K27me2_3_Rep1.bam -f BAMPE -n 111-13_H3K56R40_H3K27me3 -c 111_31_ChIP_K56R_40_input_Rep1_S31_R1_001_val_1.fq.gz.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir ../MACSPeaks --min-length 800 --max-gap 500
macs3 callpeak -t 111_2_ChIP_dRTT109_K27me2_3_Rep1.bam -f BAMPE -n 111-2_H3K56R40_H3K27me3 -c 111_1_ChIP_dRTT109_input_Rep1.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir ../MACSPeaks --min-length 800 --max-gap 500
macs3 callpeak -t 111_22_ChIP_wt_K27me3_Rep1.bam -f BAMPE -n 111-22_WT_H3K27me3 -c 111_20_ChIP_wt_input_Rep1.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir ../MACSPeaks --min-length 800 --max-gap 500

macs3 callpeak -t 111_9_ChIP_K56R_13_K27me3_Rep1_S9_R1_001_val_1.fq.gz.bam -f BAMPE -n 111-9_H3K56R13_H3K27me3 -c 111_28_ChIP_K56R_13_input_Rep1.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir ../MACSPeaks --min-length 800 --max-gap 500
macs3 callpeak -t 133-81_ChIP_rtt109_hph_H3K27me3_Rep1_S78_L001_R1_001_val_1.fq.gz.bam -f BAMPE -n 133-81_rtt109_H3K27me3 -c 133-80_ChIP_rtt109_hph_Input_Rep1_S77_L001_R1_001_val_1.fq.gz.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir ../MACSPeaks --min-length 800 --max-gap 500


macs3 callpeak -t 135-27_ChIP_WT_H3K27me3_Rep1.bam -f BAMPE -n 135-27_WT_H3K27me3 -c 135-26_ChIP_WT_Input_Rep1.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir ../MACSPeaks --min-length 800 --max-gap 500
macs3 callpeak -t 135-83_ChIP_rtt109_H3K27me3_Rep2_S79_L001_R1_001_val_1.fq.gz.bam -f BAMPE -n 135-27_rtt109_H3K27me3 -c 135-80_ChIP_rtt109_Input_Rep2_S76_L001_R1_001_val_1.fq.gz.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir ../MACSPeaks --min-length 800 --max-gap 500

macs3 callpeak -t 137-83_ChIP_rtt109_H3K27me3_Rep3_S78_R1_001_val_1.fq.gz.bam -f BAMPE -n 137-83_rtt109_H3K27me3 -c 135-80_ChIP_rtt109_Input_Rep2_S76_L001_R1_001_val_1.fq.gz.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir ../MACSPeaks --min-length 800 --max-gap 500

macs3 callpeak -t 137-85_ChIP_K56R13_H3K27me3_Rep1_S80_R1_001_val_1.fq.gz.bam -f BAMPE -n 137-85_H3K56R13_H3K27me3 -c 138-47_ChIP_H3K56R13_Input_Rep2_6252_S46_L001_R1_001_val_1.fq.gz.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir ../MACSPeaks --min-length 800 --max-gap 500

macs3 callpeak -t 138-57_ChIP_WT_H3K27me3_Rep3_6252_S56_L001_R1_001_val_1.fq.gz.bam -f BAMPE -n 138-57_WT_H3K27me3 -c 138-72_ChIP_WT_input__6252_S71_L001_R1_001_val_1.fq.gz.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir ../MACSPeaks --min-length 800 --max-gap 500







for infile in $OUTDIR/bams/*rtt109*.sorted.bam
 do
  base=$(basename ${infile} .sorted.bam)
  base2=$(basename ${infile} .EColi.sorted.bam)
macs3 callpeak -t $infile -f BAMPE -n $base -c $OUTDIR/SortedBamFiles/137-9_CUTANDRUN_rtt109_IgG_Rep1_S9.sorted.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $OUTDIR/MACSPeaks --min-length 800 --max-gap 500
macs3 callpeak -t $infile -f BAMPE -n $base2 -c $OUTDIR/SortedBamFiles/137-9_CUTANDRUN_rtt109_IgG_Rep1_S9_Ecoli.sorted.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $OUTDIR/MACSPeaks --min-length 800 --max-gap 500
done

for infile in $OUTDIR/SortedBamFiles/*WT*.sorted.bam
 do
  base=$(basename ${infile} .sorted.bam)
  base2=$(basename ${infile} .EColi.sorted.bam)
macs3 callpeak -t $infile -f BAMPE -n $base -c $OUTDIR/SortedBamFiles/137-1_CUTANDRUN_WT_IgG_Rep1_S1.sorted.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $OUTDIR/MACSPeaks --min-length 800 --max-gap 500
macs3 callpeak -t $infile -f BAMPE -n $base2.Ecoli -c $OUTDIR/SortedBamFiles/137-1_CUTANDRUN_WT_IgG_Rep1_S1_Ecoli.sorted.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $OUTDIR/MACSPeaks --min-length 800 --max-gap 500
done


for infile in $OUTDIR/SortedBamFiles/*set-7*.sorted.bam
 do
  base=$(basename ${infile} .sorted.bam)
  base2=$(basename ${infile} .EColi.sorted.bam)
macs3 callpeak -t $infile -f BAMPE -n $base -c $OUTDIR/SortedBamFiles/137-6_CUTANDRUN_set-7_IgG_Rep1_S6.sorted.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $OUTDIR/MACSPeaks --min-length 800 --max-gap 500
macs3 callpeak -t $infile -f BAMPE -n $base2.Ecoli -c $OUTDIR/SortedBamFiles/137-6_CUTANDRUN_set-7_IgG_Rep1_S6_Ecoli.sorted.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $OUTDIR/MACSPeaks --min-length 800 --max-gap 500
done

for infile in $OUTDIR/SortedBamFiles/*set-7*.sorted.bam
 do
  base=$(basename ${infile} .sorted.bam)
  base2=$(basename ${infile} .EColi.sorted.bam)
macs3 callpeak -t $infile -f BAMPE -n $base -c $OUTDIR/SortedBamFiles/137-6_CUTANDRUN_set-7_IgG_Rep1_S6.sorted.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $OUTDIR/MACSPeaks --min-length 800 --max-gap 500
macs3 callpeak -t $infile -f BAMPE -n $base2.Ecoli -c $OUTDIR/SortedBamFiles/137-6_CUTANDRUN_set-7_IgG_Rep1_S6_Ecoli.sorted.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $OUTDIR/MACSPeaks --min-length 800 --max-gap 500
done

for infile in $OUTDIR/SortedBamFiles/*ncu06787*.sorted.bam
 do
  base=$(basename ${infile} .sorted.bam)
  base2=$(basename ${infile} .EColi.sorted.bam)
macs3 callpeak -t $infile -f BAMPE -n $base -c $OUTDIR/SortedBamFiles/137-12_CUTANDRUN_ncu06787_IgG_Rep1_S12.sorted.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $OUTDIR/MACSPeaks --min-length 800 --max-gap 500
macs3 callpeak -t $infile -f BAMPE -n $base2.Ecoli -c $OUTDIR/SortedBamFiles/137-12_CUTANDRUN_ncu06787_IgG_Rep1_S12_Ecoli.sorted.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $OUTDIR/MACSPeaks --min-length 800 --max-gap 500
done

for infile in $OUTDIR/SortedBamFiles/*ncu06788*.sorted.bam
 do
  base=$(basename ${infile} .sorted.bam)
  base2=$(basename ${infile} .EColi.sorted.bam)
macs3 callpeak -t $infile -f BAMPE -n $base -c $OUTDIR/SortedBamFiles/137-15_CUTANDRUN_ncu06788_IgG_Rep1_S15.sorted.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $OUTDIR/MACSPeaks --min-length 800 --max-gap 500
macs3 callpeak -t $infile -f BAMPE -n $base2.Ecoli -c $OUTDIR/SortedBamFiles/137-15_CUTANDRUN_ncu06788_IgG_Rep1_S15_Ecoli.sorted.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $OUTDIR/MACSPeaks --min-length 800 --max-gap 500
done

for infile in $OUTDIR/SortedBamFiles/*ncu00423*.sorted.bam
 do
  base=$(basename ${infile} .sorted.bam)
  base2=$(basename ${infile} .EColi.sorted.bam)
macs3 callpeak -t $infile -f BAMPE -n $base -c $OUTDIR/SortedBamFiles/137-18_CUTANDRUN_ncu00423_IgG_Rep1_S18.sorted.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $OUTDIR/MACSPeaks --min-length 800 --max-gap 500
macs3 callpeak -t $infile -f BAMPE -n $base2.Ecoli -c $OUTDIR/SortedBamFiles/137-18_CUTANDRUN_ncu00423_IgG_Rep1_S18_Ecoli.sorted.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $OUTDIR/MACSPeaks --min-length 800 --max-gap 500
done
