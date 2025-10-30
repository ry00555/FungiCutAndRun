#!/bin/bash
#SBATCH --job-name=Macs3_Rerun
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=10:00:00
#SBATCH --output=../Macs3_Rerun.%j.out
#SBATCH --error=../Macs3_Rerun.%j.err

# --- Paths ---
BAMDIR=/scratch/ry00555/RNASeqPaper/Oct2025/SortedBamFiles
OUTDIR=/scratch/ry00555/RNASeqPaper/Oct2025/MACSPeaks
GENOME_SIZE=41037538
ml MACS3/3.0.1-gfbf-2023a

# Convert potential Windows line endings
dos2unix "$META" 2>/dev/null || true
echo "🚀 Recalling MACS3 for missing peak files individually..."

# 1. 131-72_WT_H3K36me2_rep1
macs3 callpeak \
    -t "${BAMDIR}/131-72_ChIP_WT_H3K36me2_Rep2_S60_L001_R1_001_val_1.fq.gz_Q30.bam" \
    -c "${BAMDIR}/131-37_ChIP_WT_input_Rep1_S27_L001_R1_001_val_1.fq.gz_Q30.bam" \
    -f BAMPE \
    -n "131-72_WT_H3K36me2_rep1" \
    --broad --broad-cutoff 0.1 \
    -g $GENOME_SIZE --outdir "$OUTDIR" --min-length 800 --max-gap 500

# 2. 142-105_swd-1_Input_rep1 (no control)
macs3 callpeak \
    -t "${BAMDIR}/142-105_ChIP_swd1_Input__Q30.bam" \
    -f BAMPE \
    -n "142-105_swd-1_Input_rep1" \
    --broad --broad-cutoff 0.1 \
    -g $GENOME_SIZE --outdir "$OUTDIR" --min-length 800 --max-gap 500

# 3. 131-54_WT_H3K36me3_rep2
macs3 callpeak \
    -t "${BAMDIR}/131-54_ChIP_WT_H3K36me3_Rep1_S42_L001_R1_001_val_1.fq.gz_Q30.bam" \
    -c "${BAMDIR}/131-37_ChIP_WT_input_Rep1_S27_L001_R1_001_val_1.fq.gz_Q30.bam" \
    -f BAMPE \
    -n "131-54_WT_H3K36me3_rep2" \
    --broad --broad-cutoff 0.1 \
    -g $GENOME_SIZE --outdir "$OUTDIR" --min-length 800 --max-gap 500

# 4. 134-12_cdp-6_H3K36me3_rep1
macs3 callpeak \
    -t "${BAMDIR}/134-12_ChIP_ncu06787_H3K36me3_Rep2_S9_L001_R1_001_val_1.fq.gz.bam" \
    -c "${BAMDIR}/134-13_ChIP_ncu06788_Input_Rep2.bam" \
    -f BAMPE \
    -n "134-12_cdp-6_H3K36me3_rep1" \
    --broad --broad-cutoff 0.1 \
    -g $GENOME_SIZE --outdir "$OUTDIR" --min-length 800 --max-gap 500

# 5. 134-15_cdp-6_H3K36me3_rep2
macs3 callpeak \
    -t "${BAMDIR}/134-15_ChIP_ncu06788_H3K36me3_Rep2.bam" \
    -c "${BAMDIR}/134-13_ChIP_ncu06788_Input_Rep2.bam" \
    -f BAMPE \
    -n "134-15_cdp-6_H3K36me3_rep2" \
    --broad --broad-cutoff 0.1 \
    -g $GENOME_SIZE --outdir "$OUTDIR" --min-length 800 --max-gap 500

# 6. 142-106_swd-1_H3K27me3_rep1
macs3 callpeak \
    -t "${BAMDIR}/142-106_ChIP_swd1_H3K27me3__Q30.bam" \
    -c "${BAMDIR}/142-105_ChIP_swd1_Input__Q30.bam" \
    -f BAMPE \
    -n "142-106_swd-1_H3K27me3_rep1" \
    --broad --broad-cutoff 0.1 \
    -g $GENOME_SIZE --outdir "$OUTDIR" --min-length 800 --max-gap 500

# 7. 132-29_rco-1_Input_rep1 (no control)
macs3 callpeak \
    -t "${BAMDIR}/132-29_ChIP_ncu00423_Input_Rep_1_Q30.bam" \
    -f BAMPE \
    -n "132-29_rco-1_Input_rep1" \
    --broad --broad-cutoff 0.1 \
    -g $GENOME_SIZE --outdir "$OUTDIR" --min-length 800 --max-gap 500

# 8. 148-130_nst4_H3K4ac_rep1
macs3 callpeak \
    -t "${BAMDIR}/148-130_ChIP_nst4_H3K4ac_Rep1_Q30.bam" \
    -c "${BAMDIR}/148-129_CHIP_nst4_Input_Rep1_Q30.bam" \
    -f BAMPE \
    -n "148-130_nst4_H3K4ac_rep1" \
    --broad --broad-cutoff 0.1 \
    -g $GENOME_SIZE --outdir "$OUTDIR" --min-length 800 --max-gap 500

    macs3 callpeak \
        -t "${BAMDIR}/148-133_ChIP_nst4_H3K4ac_Rep1_Q30.bam" \
        -c "${BAMDIR}/148-132_ChIP_nst4_Input_Rep1_Q30.bam" \
        -f BAMPE \
        -n "148-133_nst4_H3K4ac_rep2" \
        --broad --broad-cutoff 0.1 \
        -g $GENOME_SIZE --outdir "$OUTDIR" --min-length 800 --max-gap 500

# 9. 149-32_set-7_H3K36me3_rep8
macs3 callpeak \
    -t "${BAMDIR}/149-32__set7hph_H3K36me3_Rep1_S32_L004_R1_001_val_1.fq.gz_Q30.bam" \
    -c "${BAMDIR}/149-23__set1set7_Input_Rep1_S23_L004_R1_001_val_1.fq.gz_Q30.bam" \
    -f BAMPE \
    -n "149-32_set-7_H3K36me3_rep8" \
    --broad --broad-cutoff 0.1 \
    -g $GENOME_SIZE --outdir "$OUTDIR" --min-length 800 --max-gap 500

# 10. 137-27_WT_H3K27me3_rep5
macs3 callpeak \
    -t "${BAMDIR}/137-27_ChIP_WT_H3K27me3_Rep3_S27_R1_001_val_1.fq.gz.bam" \
    -c "${BAMDIR}/131-37_ChIP_WT_input_Rep1_S27_L001_R1_001_val_1.fq.gz_Q30.bam" \
    -f BAMPE \
    -n "137-27_WT_H3K27me3_rep5" \
    --broad --broad-cutoff 0.1 \
    -g $GENOME_SIZE --outdir "$OUTDIR" --min-length 800 --max-gap 500

# 11. 131-37_WT_Input_rep3 (no control)
macs3 callpeak \
    -t "${BAMDIR}/131-37_ChIP_WT_input_Rep1_S27_L001_R1_001_val_1.fq.gz_Q30.bam" \
    -f BAMPE \
    -n "131-37_WT_Input_rep3" \
    --broad --broad-cutoff 0.1 \
    -g $GENOME_SIZE --outdir "$OUTDIR" --min-length 800 --max-gap 500

# 12. 133-78_WT_H3K27me3_rep2
macs3 callpeak \
    -t "${BAMDIR}/133-78_ChIP_WT_H3K27me3_Rep1_Q30.bam" \
    -c "${BAMDIR}/131-37_ChIP_WT_input_Rep1_S27_L001_R1_001_val_1.fq.gz_Q30.bam" \
    -f BAMPE \
    -n "133-78_WT_H3K27me3_rep2" \
    --broad --broad-cutoff 0.1 \
    -g $GENOME_SIZE --outdir "$OUTDIR" --min-length 800 --max-gap 500

echo "🎯 All individual MACS3 peak calls complete!"
