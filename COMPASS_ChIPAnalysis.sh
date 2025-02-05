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
# mkdir -p "${OUTDIR}/TrimmedReads"
 mkdir -p "${OUTDIR}/BigWigs"
 mkdir -p "$OUTDIR/HomerTagDirectories"
 mkdir -p "$OUTDIR/TdfFiles"
 mkdir -p "$OUTDIR/SortedBamFiles"
#
 fi



#ml Trim_Galore
#trim_galore --illumina --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${OUTDIR}/FASTQ/132*fastq\.gz
#trim_galore --illumina --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${OUTDIR}/FASTQ/133*fastq\.gz
#trim_galore --illumina --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${OUTDIR}/FASTQ/106*fastq\.gz


TAGDIR="${OUTDIR}/HomerTagDirectories"
BAMDIR="${OUTDIR}/SortedBamFiles"
BEDDIR="${OUTDIR}/Beds"
PEAKDIR="${OUTDIR}/MACSPeaks"

 #FILES="${OUTDIR}/TrimmedReads/106*_R1_001_val_1\.fq\.gz"

# for f in $FILES
#  do
#  file=${f##*/}
# name=${file/%_S[1-990]*_R1_001_val_1.fq.gz/}
  #read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R2_001_val_2\.fq\.gz/g')

 #bam="${OUTDIR}/SortedBamFiles/${name}.bam"
# variable name for bigwig output
 	#bigwig="${OUTDIR}/BigWigs/${name}"
 #QualityBam="${OUTDIR}/SortedBamFiles/${name}_Q30.bam"
  #ml SAMtools/1.16.1-GCC-11.3.0
  #ml BWA/0.7.17-GCCcore-11.3.0
#  bwa mem -M -v 3 -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIR/SortedBamFiles/tempReps -o "$bam" -
  #bwa mem -M -v 3 -t $THREADS $GENOME $f | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIR/SortedBamFiles/tempReps -o "$bam" -

  #samtools index "$bam"
  #samtools view -b -q 30 $bam > "$QualityBam"
  #samtools index "$QualityBam"
 #ml deepTools/3.5.2-foss-2022a
 #bamCoverage -p $THREADS -bs 10 --normalizeUsing BPM --minMappingQuality 10 --smoothLength $SMOOTH -of bigwig -b "$QualityBam" -o "${bigwig}.bin_10.smooth_${SMOOTH}Bulk.bw"
#  done

#ml deepTools
#computeMatrix scale-regions -p 12 -R /scratch/ry00555/neurospora.bed -S ${OUTDIR}/BigWigs/


#mkdir $OUTDIR/MACSPeaks


#module load MACS3/3.0.0b1-foss-2022a-Python-3.10.4
 #command line
# macs3 callpeak -t $BAMDIR/142-14_ChIP_iswKO_H3K27me3_Rep2_Q30.bam -f BAMPE -n 142-14_ChIP_iswKO_H3K27me3_Rep2 -c $BAMDIR/135-66_ChIP_isw_Input_Rep1_Q30.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
 #macs3 callpeak -t $BAMDIR/133-78_ChIP_WT_H3K27me3_Rep1_Q30.bam -f BAMPE -n 133-78_ChIP_WT_H3K27me3 --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
 #macs3 callpeak -t $BAMDIR/142-45_ChIP_rco1KOA_H3K27me3__Q30.bam -f BAMPE -n 142-45_ChIP_rco1KOA_H3K27me3 -c $BAMDIR/142-44_ChIP_rco1KOA_Input__Q30.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
 #macs3 callpeak -t $BAMDIR/142-42_ChIP_rco1KOa_H3K27me3__Q30.bam -f BAMPE -n 142-42_ChIP_rco1KOa_H3K27me3 -c $BAMDIR/142-41_ChIP_rco1KOa_Input__Q30.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
 #macs3 callpeak -t $BAMDIR/132-35_ChIP_ncu00423_H3K27me3_Rep_1_Q30.bam -f BAMPE -n 132-35_ChIP_ncu00423_H3K27me3 -c $BAMDIR/132-34_ChIP_ncu00423_Input_Rep_1_Q30.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500




# macs3 callpeak -t $BAMDIR/135-30_ChIP_pot-1__nat_H3K27me3_Rep1_Q30.bam -f BAMPE -n 135-30_ChIP_pot-1__nat_H3K27me3_Rep1 -c $BAMDIR/135-29_ChIP_pot-1__nat_Input_Rep1_Q30.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
#  macs3 callpeak -t $BAMDIR/137-89_ChIP_pot-1_H3K27me3_Rep2_Q30.bam -f BAMPE -n 137-89_ChIP_pot-1_H3K27me3_Rep2 -c $BAMDIR/137-90_ChIP_pot-1_H3K36me3_Rep2_Q30.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
# macs3 callpeak -t $BAMDIR/138-55_ChIP_EPR1_H3K27me3_Rep1_6252_Q30.bam -f BAMPE -n 138-55_ChIP_EPR1_H3K27me3_Rep1 -c $BAMDIR/138-53_ChIP_EPR1_Input_Rep1_6252_Q30.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
# macs3 callpeak -t $BAMDIR/135-2_ChIP_ncu00423_H3K27me3_Rep3_Q30.bam  -f BAMPE -n 135-2_ChIP_ncu00423_H3K27me3_Rep3 -c $BAMDIR/135-1_ChIP_ncu00423_Input_Rep3_Q30.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
 #macs3 callpeak -t $BAMDIR/135-27_ChIP_WT_H3K27me3_Rep1_Q30.bam   -f BAMPE -n 135-27_ChIP_WT_H3K27me3_Rep1 -c $BAMDIR/135-26_ChIP_WT_Input_Rep1_Q30.bam  --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
# macs3 callpeak -t $BAMDIR/138-57_ChIP_WT_H3K27me3_Rep3_6252_Q30.bam  -f BAMPE -n 138-57_ChIP_WT_H3K27me3_Rep3 -c $BAMDIR/138-72_ChIP_WT_input__6252_Q30.bam  --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
#macs3 callpeak -t $BAMDIR/137-70_ChIP_WT_H3K27me3_Rep1_Q30.bam   -f BAMPE -n 137-70_ChIP_WT_H3K27me3_Rep1 -c $BAMDIR/137-63_ChIP_WT_input__Q30.bam  --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500





#  macs3 callpeak -t $BAMDIR/142-115_ChIP_set2_H3K27me3__Q30.bam -f BAMPE -n 142-115_ChIP_set2_H3K27me3  --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
# macs3 callpeak -t $BAMDIR/142-118_ChIP_set2_H3K27me3__Q30.bam -f BAMPE -n 142-118_ChIP_set2_H3K27me3  --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
#macs3 callpeak -t $BAMDIR/145-110_ChIP_set1E8_H3K27me3_Rep2_Q30.bam -f BAMPE -n 145-110_ChIP_set1E8_H3K27me3_Rep2 -c $BAMDIR/145-41_ChIP_set1E8_Input_Rep2_Q30.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
#macs3 callpeak -t $BAMDIR/145-118_ChIP_swd1_H3K27me3_Rep2_Q30.bam -f BAMPE -n 145-118_ChIP_swd1_H3K27me3_Rep2 -c $BAMDIR/145-116_ChIP_swd1_Input_Rep2_Q30.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
#macs3 callpeak -t $BAMDIR/145-35_ChIP_set7_H3K27me3_Rep2_Q30.bam -f BAMPE -n 145-35_ChIP_set7_H3K27me3_Rep2 -c $BAMDIR/145-33_ChIP_set7_Input_Rep2_Q30.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
#macs3 callpeak -t $BAMDIR/145-114_ChIP_sgr9_H3K27me3_Rep2_Q30.bam -f BAMPE -n 145-114_ChIP_sgr9_H3K27me3_Rep2 -c $BAMDIR/145-112_ChIP_sgr9_Input_Rep2_Q30.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
#macs3 callpeak -t $BAMDIR/145-39_ChIP_set1E7_H3K27me3_Rep2_Q30.bam -f BAMPE -n 145-39_ChIP_set1E7_H3K27me3_Rep2 -c $BAMDIR/145-37_ChIP_set1E7_Input_Rep2_Q30.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
#macs3 callpeak -t $BAMDIR/145-30_ChIP_WT_H3K27me3_Q30.bam -f BAMPE -n 145-30_ChIP_WT_H3K27me3 -c $BAMDIR/145-29_ChIP_WT_Input_Q30.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500



#ml ChIP-R
#set7
#chipr -i 142-94_ChIP_set7_H3K27me3_peaks.broadPeak 145-35_ChIP_set7_H3K27me3_Rep2_peaks.broadPeak -m 2 -o Intersected_set7_H3K27me3
#WT
#chipr -i ${PEAKDIR}/142-10_ChIP_WT_H3K27me3_peaks.broadPeak ${PEAKDIR}/145-30_ChIP_WT_H3K27me3_peaks.broadPeak -m 2 -o ${PEAKDIR}/Intersected_WT_H3K27me3
#set1
#chipr -i ${PEAKDIR}/142-121_ChIP_set1_H3K27me3_peaks.broadPeak ${PEAKDIR}/142-124_ChIP_set1_H3K27me3_peaks.broadPeak ${PEAKDIR}/145-110_ChIP_set1E8_H3K27me3_Rep2_peaks.broadPeak ${PEAKDIR}/145-39_ChIP_set1E7_H3K27me3_Rep2_peaks.broadPeak -m 2 -o ${PEAKDIR}/Intersected_set1_H3K27me3
#sgr9 1
#chipr -i ${PEAKDIR}/142-127_ChIP_sgr9_H3K27me3_peaks.broadPeak ${PEAKDIR}/145-114_ChIP_sgr9_H3K27me3_Rep2_peaks.broadPeak -m 2 -o ${PEAKDIR}/Intersected_sgr9_H3K27me3
#swd1
#chipr -i ${PEAKDIR}/145-118_ChIP_swd1_H3K27me3_Rep2_peaks.broadPeak ${PEAKDIR}/142-106_ChIP_swd1_H3K27me3_peaks.broadPeak -m 2 -o ${PEAKDIR}/Intersected_swd1_H3K27me3
#set2
#chipr -i ${PEAKDIR}/142-115_ChIP_set2_H3K27me3_peaks.broadPeak ${PEAKDIR}/142-118_ChIP_set2_H3K27me3_peaks.broadPeak -m 2 -o ${PEAKDIR}/Intersected_set2_H3K27me3

#awk '{print $1, $2, $3, $5}' Intersected_WT_H3K27me3_all.bed > BedGraphs/Intersected_WT_H3K27me3_all.bedgraph
#awk '{print $1, $2, $3, $5}' Intersected_set1_H3K27me3_all.bed > BedGraphs/Intersected_set1_H3K27me3_all.bedgraph
#awk '{print $1, $2, $3, $5}' Intersected_swd1_H3K27me3_all.bed > BedGraphs/Intersected_swd1_H3K27me3_all.bedgraph
#awk '{print $1, $2, $3, $5}' Intersected_set2_H3K27me3_all.bed > BedGraphs/Intersected_set2_H3K27me3_all.bedgraph
#awk '{print $1, $2, $3, $5}' Intersected_sgr9_H3K27me3_all.bed > BedGraphs/Intersected_sgr9_H3K27me3_all.bedraph
#awk '{print $1, $2, $3, $5}' Intersected_set7_H3K27me3_all.bed  > BedGraphs/Intersected_set7_H3K27me3_all.bedgraph


#ml BEDTools
#ml ucsc
#bedGraphToBigWig [options] in.bedGraph chrom.sizes out.bw
# for infile in $PEAKDIR/BedGraphs/*bedgraph
#do
  #base=$(basename ${infile} .bedgraph)
  #bedSort $infile $PEAKDIR/BedGraphs/${base}.sort.bedgraph
 #done

#mkdir $OUTDIR/bigwigs
#ml deepTools
#for infile in $PEAKDIR/BedGraphs/*.sort.bedgraph
# do
#  base=$(basename ${infile} .sort.bedgraph)
# bedGraphToBigWig $infile /scratch/ry00555/Run137CutandRun/ref/Ncrassa_ref/chrNameLength.txt $PEAKDIR/BedGraphs/${base}.bw
 #done
#bedGraphToBigWig Intersected_set7_H3K27me3_all.sort.bedgraph /scratch/ry00555/Run137CutandRun/ref/Ncrassa_ref/chrNameLength.txt Intersected_set7_H3K27me3_all.sort.bw

#bedtools intersect -wa -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GenesOnly.bed -b $PEAKDIR/Intersected_WT_H3K27me3_all.bed > $BEDDIR/MACS_WT_IntersectedH3K27me3_peaks.bed
#bedtools intersect -wa -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GenesOnly.bed -b $PEAKDIR/Intersected_set2_H3K27me3_all.bed > $BEDDIR/MACS_set2_IntersectedH3K27me3_peaks.bed
#bedtools intersect -wa -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GenesOnly.bed -b $PEAKDIR/Intersected_set1_H3K27me3_all.bed > $BEDDIR/MACS_set1_IntersectedH3K27me3_peaks.bed
#bedtools intersect -wa -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GenesOnly.bed -b $PEAKDIR/Intersected_swd1_H3K27me3_all.bed > $BEDDIR/MACS_swd1_IntersectedH3K27me3_peaks.bed
#bedtools intersect -wa -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GenesOnly.bed -b $PEAKDIR/Intersected_sgr9_H3K27me3_all.bed > $BEDDIR/MACS_sgr9_IntersectedH3K27me3_peaks.bed
#bedtools intersect -wa -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GenesOnly.bed -b $PEAKDIR/Intersected_set7_H3K27me3_all.bed > $BEDDIR/MACS_set7_IntersectedH3K27me3_peaks.bed


#mkdir ${OUTDIR}/HomerPeaks
#HOMERPEAKSDIR="${OUTDIR}/HomerPeaks"
#ml Homer
#ml Perl

#for bam_file in "${BAMDIR}"/*_Q30.bam; do
#  sample_id=$(basename "${bam_file}" _Q30.bam)
#makeTagDirectory "${TAGDIR}/${sample_id}" "${bam_file}"
#done
#makeTagDirectory 145-116_ChIP_swd1_Input_Rep2 ../SortedBamFiles/145-116_ChIP_swd1_Input_Rep2_Q30.bam

#findPeaks ${TAGDIR}/142-127_ChIP_sgr9_H3K27me3 -style histone -region -size 150 -minDist 530 -o ${HOMERPEAKSDIR}/142-127_ChIP_sgr9_H3K27me3_Homerpeaks.txt -i ${TAGDIR}/142-126_ChIP_sgr9_Input
#findPeaks ${TAGDIR}/145-35_ChIP_set7_H3K27me3_Rep2 -style histone -region -size 150 -minDist 530 -o ${HOMERPEAKSDIR}/145-35_ChIP_set7_H3K27me3_Rep2_Homerpeaks.txt -i ${TAGDIR}/145-33_ChIP_set7_Input_Rep2
#findPeaks ${TAGDIR}/142-106_ChIP_swd1_H3K27me3 -style histone -region -size 150 -minDist 530 -o ${HOMERPEAKSDIR}/142-106_ChIP_swd1_H3K27me3_Homerpeaks.txt
#findPeaks ${TAGDIR}/142-121_ChIP_set1_H3K27me3 -style histone -region -size 150 -minDist 530 -o ${HOMERPEAKSDIR}/142-121_ChIP_set1_H3K27me3_Homerpeaks.txt -i ${TAGDIR}/142-123_ChIP_set1_Input
#findPeaks ${TAGDIR}/142-124_ChIP_set1_H3K27me3 -style histone -region -size 150 -minDist 530 -o ${HOMERPEAKSDIR}/142-124_ChIP_set1_H3K27me3_Homerpeaks.txt -i ${TAGDIR}/142-123_ChIP_set1_Input
#findPeaks ${TAGDIR}/142-10_ChIP_WT_H3K27me3_Rep3 -style histone -region -size 150 -minDist 530 -o ${HOMERPEAKSDIR}/142-10_ChIP_WT_H3K27me3_Rep3_Homerpeaks.txt -i ${TAGDIR}/142-75_ChIP_WT_Input
#findPeaks ${TAGDIR}/142-115_ChIP_set2_H3K27me3 -style histone -region -size 150 -minDist 530 -o ${HOMERPEAKSDIR}/142-115_ChIP_set2_H3K27me3_Homerpeaks.txt
#findPeaks ${TAGDIR}/142-118_ChIP_set2_H3K27me3 -style histone -region -size 150 -minDist 530 -o ${HOMERPEAKSDIR}/142-118_ChIP_set2_H3K27me3_Homerpeaks.txt
#findPeaks ${TAGDIR}/145-110_ChIP_set1E8_H3K27me3_Rep2 -style histone -region -size 150 -minDist 530 -o ${HOMERPEAKSDIR}/145-110_ChIP_set1E8_H3K27me3_Rep2_Homerpeaks.txt -i ${TAGDIR}/145-41_ChIP_set1E8_Input_Rep2
#findPeaks ${TAGDIR}/145-118_ChIP_swd1_H3K27me3_Rep2 -style histone -region -size 150 -minDist 530 -o ${HOMERPEAKSDIR}/145-118_ChIP_swd1_H3K27me3_Rep2_Homerpeaks.txt -i ${TAGDIR}/145-116_ChIP_swd1_Input_Rep2
#findPeaks ${TAGDIR}/145-35_ChIP_set7_H3K27me3_Rep2 -style histone -region -size 150 -minDist 530 -o ${HOMERPEAKSDIR}/145-35_ChIP_set7_H3K27me3_Rep2_Homerpeaks.txt -i ${TAGDIR}/145-33_ChIP_set7_Input_Rep2
#findPeaks ${TAGDIR}/145-114_ChIP_sgr9_H3K27me3_Rep2 -style histone -region -size 150 -minDist 530 -o ${HOMERPEAKSDIR}/145-114_ChIP_sgr9_H3K27me3_Rep2_Homerpeaks.txt -i ${TAGDIR}/145-112_ChIP_sgr9_Input_Rep2
#findPeaks ${TAGDIR}/145-39_ChIP_set1E7_H3K27me3_Rep2 -style histone -region -size 150 -minDist 530 -o ${HOMERPEAKSDIR}/145-39_ChIP_set1E7_H3K27me3_Rep2_Homerpeaks.txt -i ${TAGDIR}/145-37_ChIP_set1E7_Input_Rep2
#findPeaks ${TAGDIR}/145-30_ChIP_WT_H3K27me3 -style histone -region -size 150 -minDist 530 -o ${HOMERPEAKSDIR}/145-30_ChIP_WT_H3K27me3_Homerpeaks.txt -i ${TAGDIR}/145-29_ChIP_WT_Input
#findPeaks ${TAGDIR}/142-94_ChIP_set7_H3K27me3  -style histone -region -size 150 -minDist 530 -o ${HOMERPEAKSDIR}/142-94_ChIP_set7_H3K27me3_Homerpeaks.txt -i ${TAGDIR}/142-93_ChIP_set7_Input
#Genome="/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_00182925.2plusHphplusBarplusTetO_his3masked.fna"
#GTF="/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_WithExtras_GFFtoGTFconversion.gtf"
#GFF="/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic.gff"
 #for infile in ${HOMERPEAKSDIR}/*_Homerpeaks.txt
#do
#  base=$(basename ${infile} _Homerpeaks.txt)
#pos2bed.pl ${HOMERPEAKSDIR}/${base}_Homerpeaks.txt > ${HOMERPEAKSDIR}/${base}.bed
 #annotatePeaks.pl ${HOMERPEAKSDIR}/${base}.peaks.bed $Genome -gff3 $GFF3 > ${HOMERPEAKSDIR}/${base}_ann.txt
#done

# 142-124_ChIP_set1_H3K27me3.bed
# 142-121_ChIP_set1_H3K27me3.bed
# 142-118_ChIP_set2_H3K27me3.bed
# 142-115_ChIP_set2_H3K27me3.bed
# 142-10_ChIP_WT_H3K27me3_Rep3.bed
# 142-106_ChIP_swd1_H3K27me3.bed
# 145-39_ChIP_set1E7_H3K27me3_Rep2.bed
# 145-35_ChIP_set7_H3K27me3_Rep2.bed
# 145-30_ChIP_WT_H3K27me3.bed
# 145-118_ChIP_swd1_H3K27me3_Rep2.bed
# 145-114_ChIP_sgr9_H3K27me3_Rep2.bed
# -rw-r--r-- 1 ry00555 zallab  32881 Jan 14 16:27 145-110_ChIP_set1E8_H3K27me3_Rep2.bed
# -rw-r--r-- 1 ry00555 zallab   8139 Jan 14 16:27 142-94_ChIP_set7_H3K27me3.bed
# -rw-r--r-- 1 ry00555 zallab  35419 Jan 14 16:27 142-127_ChIP_sgr9_H3K27me3.bed


#mergePeaks
#ml Homer
#mergePeaks 142-106_ChIP_swd1_H3K27me3.peaks.bed 145-118_ChIP_swd1_H3K27me3_Rep2.bed >  ../Beds/Homer_merged_swd1_h3K27me3.bed

#mergePeaks 142-10_ChIP_WT_H3K27me3_Rep3.peaks.bed 145-30_ChIP_WT_H3K27me3.peaks.bed > ../Beds/Homer_merged_WT_h3K27me3.bed

#mergePeaks 142-115_ChIP_set2_H3K27me3.peaks.bed 142-118_ChIP_set2_H3K27me3.peaks.bed > ../Beds/Homer_merged_set2_h3K27me3.bed

#mergePeaks 142-121_ChIP_set1_H3K27me3.peaks.bed 142-124_ChIP_set1_H3K27me3.peaks.bed 145-110_ChIP_set1E8_H3K27me3_Rep2.peaks.bed 145-39_ChIP_set1E7_H3K27me3_Rep2.peaks.bed > ../Beds/Homer_merged_set1_h3K27me3.bed

#mergePeaks 145-114_ChIP_sgr9_H3K27me3_Rep2.peaks.bed 142-127_ChIP_sgr9_H3K27me3.peaks.bed > ../Beds/Homer_merged_sgr9_h3K27me3.bed

#mergePeaks 142-94_ChIP_set7_H3K27me3.peaks.bed 145-35_ChIP_set7_H3K27me3_Rep2.peaks.bed > ../Beds/Homer_merged_set7_h3K27me3.bed


#ml BEDTools
 # bedtools multiinter -header -i HOMER_Annotated_142-10_ChIP_WT_H3K27me3_Rep3.peaks.bed  HOMER_Annotated_145-30_ChIP_WT_H3K27me3.peaks.bed  > WT_H3K27me3_Homer_Merged_peaks.bed
 #
 # bedtools multiinter -header -i HOMER_Annotated_142-115_ChIP_set2_H3K27me3.peaks.bed  HOMER_Annotated_142-118_ChIP_set2_H3K27me3.peaks.bed  > set2_H3K27me3_Homer_Merged_peaks.bed
 #
 # bedtools multiinter -header -i HOMER_Annotated_142-121_ChIP_set1_H3K27me3.peaks.bed  HOMER_Annotated_142-124_ChIP_set1_H3K27me3.peaks.bed HOMER_Annotated_145-110_ChIP_set1E8_H3K27me3_Rep2.peaks.bed HOMER_Annotated_145-39_ChIP_set1E7_H3K27me3_Rep2.peaks.bed > set1_H3K27me3_Homer_Merged_peaks.bed
 #
 # bedtools multiinter -header -i HOMER_Annotated_142-127_ChIP_sgr9_H3K27me3.peaks.bed  HOMER_Annotated_145-114_ChIP_sgr9_H3K27me3_Rep2.peaks.bed  > sgr9_H3K27me3_Homer_Merged_peaks.bed
 # bedtools multiinter -header -i HOMER_Annotated_142-94_ChIP_set7_H3K27me3.peaks.bed  HOMER_Annotated_145-35_ChIP_set7_H3K27me3_Rep2.peaks.bed > set7_H3K27me3_Homer_Merged_peaks.bed


#HOMER_Annotated_142-106_ChIP_swd1_H3K27me3.peaks.bed

# bedtools intersect -wa -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GenesOnly.bed -b WT_H3K27me3_Homer_Merged_peaks.bed > WT_H3K27me3_Homer_Merged_Annotated.peaks.bed
#
#  bedtools intersect -wa -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GenesOnly.bed -b set2_H3K27me3_Homer_Merged_peaks.bed > set2_H3K27me3_Homer_Merged_Annotated.peaks.bed
# #
#  bedtools intersect -wa -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GenesOnly.bed -b set1_H3K27me3_Homer_Merged_peaks.bed  > set1_H3K27me3_Homer_Merged_Annotated.peaks.bed
# #
# bedtools intersect -wa -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GenesOnly.bed -b sgr9_H3K27me3_Homer_Merged_peaks.bed  > sgr9_H3K27me3_Homer_Merged_Annotated.peaks.bed
# #
#  bedtools intersect -wa -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GenesOnly.bed -b set7_H3K27me3_Homer_Merged_peaks.bed  > set7_H3K27me3_Homer_Merged_Annotated.peaks.bed


 #grep -i "rRNA" /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_WithExtras.gff | awk '{print $1"\t"$4-1"\t"$5"\t"$9}' > rDNA_regions.bed
 #sort rDNA_regions.bed | uniq > rDNA_regions_no_duplicates.bed
 #sort -k1,1 -k2,2n rDNA_regions_no_duplicates.bed > rDNA_regions_sorted.bed
 #bedtools merge -i rDNA_regions_sorted.bed -c 4 -o distinct > rDNA_regions_merged.bed
# then manually narrow down to the 3 rDNA genes
#bigWigAverageOverBed ChIP_mutant.bw rDNA_regions.bed mutant_rDNA_enrichment.txt
#
# bigWigAverageOverBed Intersected_sgr9_H3K27me3_all.sort.bw /scratch/ry00555/rDNA_regions_merged.bed sgr9_rDNA_enrichment.txt
#
# bigWigAverageOverBed Intersected_WT_H3K27me3_all.bw /scratch/ry00555/neurospora.bed WT_WG_enrichment.txt
# bigWigAverageOverBed Intersected_swd1_H3K27me3_all.bw /scratch/ry00555/heatmapPRC2genes.bed swd1_K27Genes_enrichment.txt
#


ml deepTools
#bamCompare -b1 treatment.bam -b2 control.bam -o log2ratio.bw -of bigwig
#142-119_ChIP_set2_H3K36me3__Q30.bam
#bamCompare -p max -b1 ${BAMDIR}/142-125_ChIP_set1_H3K36me3__Q30.bam	-b2	${BAMDIR}/142-123_ChIP_set1_Input__Q30.bam --minMappingQuality 10  -of bigwig --effectiveGenomeSize 41037538 --skipZeroOverZero --binSize 30  --smoothLength 60 -o ${OUTDIR}/BigWigs/142-125_ChIP_set1_H3K36me3_log2ratio_testReadCountScaleFactor.bw

#bamCompare -p max -b1 ${BAMDIR}/142-125_ChIP_set1_H3K36me3__Q30.bam	-b2	${BAMDIR}/142-123_ChIP_set1_Input__Q30.bam --minMappingQuality 10  -of bigwig --scaleFactorsMethod None --effectiveGenomeSize 41037538 --skipZeroOverZero --binSize 30 --normalizeUsing BPM --smoothLength 60 -o ${OUTDIR}/BigWigs/142-125_ChIP_set1_H3K36me3_log2ratio_smooth60.bw

#bamCompare -p max -b1 ${BAMDIR}/145-40_ChIP_set1E7_H3K36me3_Rep2_Q30.bam	-b2	${BAMDIR}/145-37_ChIP_set1E7_Input_Rep2_Q30.bam  --minMappingQuality 10 -of bigwig  --scaleFactorsMethod None --effectiveGenomeSize 41037538 --skipZeroOverZero --binSize 30 --normalizeUsing BPM --smoothLength 60 -o ${OUTDIR}/BigWigs/145-40_ChIP_set1E7_H3K36me3_Rep2_log2ratio_smooth60.bw

#bamCompare -p max -b1 ${BAMDIR}/142-128_ChIP_sgr9_H3K36me3__Q30.bam	-b2	${BAMDIR}/142-126_ChIP_sgr9_Input__Q30.bam -of bigwig  --minMappingQuality 10 --scaleFactorsMethod None --effectiveGenomeSize 41037538 --skipZeroOverZero --binSize 30 --normalizeUsing BPM --smoothLength 60 -o ${OUTDIR}/BigWigs/142-128_ChIP_sgr9_H3K36me3_log2ratio_smooth60.bw

#bamCompare -p max -b1 ${BAMDIR}/142-77_ChIP_WT_H3K36me3__S77_L007_R1_001_val_1.fq.gz_Q30.bam	-b2	${BAMDIR}/142-75_ChIP_WT_Input__Q30.bam  --minMappingQuality 10  --scaleFactorsMethod None --effectiveGenomeSize 41037538 -of bigwig --skipZeroOverZero --binSize 30 --normalizeUsing BPM --smoothLength 60 -o ${OUTDIR}/BigWigs/142-77_ChIP_WT_H3K36me3_log2ratio_smooth60.bw

#bamCompare -p max -b1 ${BAMDIR}/145-115_ChIP_sgr9_H3K36me3_Rep2_Q30.bam -b2 ${BAMDIR}/145-112_ChIP_sgr9_Input_Rep2_Q30.bam -of bigwig  --minMappingQuality 10  --scaleFactorsMethod None --effectiveGenomeSize 41037538 --skipZeroOverZero --binSize 30 --normalizeUsing BPM --smoothLength 60 -o ${OUTDIR}/BigWigs/145-115_ChIP_sgr9_H3K36me3_log2ratio_smooth60.bw

#bamCompare -p max -b1 ${BAMDIR}/145-32_ChIP_S1_H3K36me3_Rep2_Q30.bam	-b2	${BAMDIR}/145-29_ChIP_WT_Input_Q30.bam -of bigwig --minMappingQuality 10  --scaleFactorsMethod None --effectiveGenomeSize 41037538 --skipZeroOverZero --binSize 30 --normalizeUsing BPM --smoothLength 60 -o ${OUTDIR}/BigWigs/145-32_ChIP_S1_H3K36me3_log2ratio_smooth60.bw

#bamCompare -p max -b1 ${BAMDIR}/145-36_ChIP_set7_H3K36me3_Rep2_Q30.bam	-b2	${BAMDIR}/145-33_ChIP_set7_Input_Rep2_Q30.bam -of bigwig --minMappingQuality 10  --scaleFactorsMethod None --effectiveGenomeSize 41037538 --skipZeroOverZero --binSize 30 --normalizeUsing BPM --smoothLength 60 -o ${OUTDIR}/BigWigs/145-36_ChIP_set7_H3K36me3_log2ratio_smooth60.bw

#computeMatrix scale-regions  -p max --startLabel "TSS'" --endLabel "TES'" -b 1000 -a 1000  --regionsFileName /scratch/ry00555/heatmapPRC2genes.bed /scratch/ry00555/TRUErDNA_genes_with_names_fixed.bed /scratch/ry00555/neighboringK27genes.bed /scratch/ry00555/rDNA_regions_no_duplicates_genes_format.bed /scratch/ry00555/ash1depgenes.bed --skipZeros  --scoreFileName ${OUTDIR}/BigWigs/142-77_ChIP_WT_H3K36me3_log2ratio.bw ${OUTDIR}/BigWigs/145-32_ChIP_S1_H3K36me3_log2ratio.bw ${OUTDIR}/BigWigs/145-36_ChIP_set7_H3K36me3_log2ratio.bw ${OUTDIR}/BigWigs/H3K36me3/142-116_ChIP_set2_H3K36me3_.bin_25.smooth_50_Q30.bw ${OUTDIR}/BigWigs/H3K36me3/142-119_ChIP_set2_H3K36me3_.bin_25.smooth_50_Q30.bw ${OUTDIR}/BigWigs/142-125_ChIP_set1_H3K36me3_log2ratio.bw ${OUTDIR}/BigWigs/145-40_ChIP_set1E7_H3K36me3_Rep2_log2ratio.bw ${OUTDIR}/BigWigs/142-128_ChIP_sgr9_H3K36me3_log2ratio.bw ${OUTDIR}/BigWigs/145-115_ChIP_sgr9_H3K36me3_log2ratio.bw  --outFileName ${OUTDIR}/Matrices/H3K36me3_log2ratiomatrix.gz

#computeMatrix scale-regions  -p max --startLabel "TSS'" --endLabel "TES'" -b 1000 -a 1000  --regionsFileName /scratch/ry00555/heatmapPRC2genes.bed /scratch/ry00555/TRUErDNA_genes_with_names_fixed.bed /scratch/ry00555/neighboringK27genes.bed /scratch/ry00555/rDNA_regions_no_duplicates_genes_format.bed /scratch/ry00555/ash1depgenes.bed --skipZeros  --scoreFileName ${OUTDIR}/BigWigs/142-77_ChIP_WT_H3K36me3_log2ratio.bw ${OUTDIR}/BigWigs/145-32_ChIP_S1_H3K36me3_log2ratio.bw ${OUTDIR}/BigWigs/145-36_ChIP_set7_H3K36me3_log2ratio.bw ${OUTDIR}/BigWigs/142-125_ChIP_set1_H3K36me3_log2ratio.bw ${OUTDIR}/BigWigs/145-40_ChIP_set1E7_H3K36me3_Rep2_log2ratio.bw ${OUTDIR}/BigWigs/142-128_ChIP_sgr9_H3K36me3_log2ratio.bw ${OUTDIR}/BigWigs/145-115_ChIP_sgr9_H3K36me3_log2ratio.bw  --outFileName ${OUTDIR}/Matrices/H3K36me3_log2ratiomatrix_NoSet2.gz

#plotProfile --startLabel "TSS'" --endLabel "TES'" --averageType mean  --matrixFile ${OUTDIR}/Matrices/H3K36me3_log2ratiomatrix_NoSet2.gz --outFileName ${OUTDIR}/Heatmaps/H3K36me3_log2ratiomatrix_mean_noSET2_V1.png  --outFileNameData ${OUTDIR}/Heatmaps/H3K36me3_log2ratiomatrix_mean_noSET2.tab --numPlotsPerRow 2 --plotType=fill  --perGroup --legendLocation lower-right --samplesLabel "142-77-WT" "145-32-WT" "145-36-set7" "142-125-set1" "145-40-set1E7" "142-128-sgr9" "145-115-sgr9" --regionsLabel "H3K27me3 marked genes" "True rDNA genes" "Neighboring H3K27me3 genes" "Pseudo rDNA genes" "Ash1 H3K36me genes" --plotTitle "L2FC Input Normalized H3K36me3_COMPASS Mean profile" --yMin 0

#plotProfile --startLabel "TSS'" --endLabel "TES'" --averageType median  --matrixFile ${OUTDIR}/Matrices/H3K36me3_log2ratiomatrix_NoSet2.gz --outFileName ${OUTDIR}/Heatmaps/H3K36me3_log2ratiomatrix_median_noSET2_V1.png  --outFileNameData ${OUTDIR}/Heatmaps/H3K36me3_log2ratiomatrix_median_noSET2.tab --numPlotsPerRow 2 --plotType=fill  --perGroup --legendLocation lower-right --samplesLabel "142-77-WT" "145-32-WT" "145-36-set7" "142-125-set1" "145-40-set1E7" "142-128-sgr9" "145-115-sgr9" --regionsLabel "H3K27me3 marked genes" "True rDNA genes" "Neighboring H3K27me3 genes" "Pseudo rDNA genes" "Ash1 H3K36me genes" --plotTitle "L2FC Input Normalized H3K36me3_COMPASS Median profile" --yMin 0

#computeMatrix scale-regions  -p max --startLabel "TSS'" --endLabel "TES'" -b 1000 -a 1000  --regionsFileName /scratch/ry00555/heatmapPRC2genes.bed /scratch/ry00555/TRUErDNA_genes_with_names_fixed.bed /scratch/ry00555/neighboringK27genes.bed /scratch/ry00555/rDNA_regions_no_duplicates_genes_format.bed /scratch/ry00555/ash1depgenes.bed --skipZeros  --scoreFileName ../BigWigs/142-77_ChIP_WT_H3K36me3_log2ratio.bw ../BigWigs/145-32_ChIP_S1_H3K36me3_log2ratio.bw ../BigWigs/145-36_ChIP_set7_H3K36me3_log2ratio.bw ../BigWigs/142-119_ChIP_set2_H3K36me3_.bin_25.smooth_50_Q30.bw ../BigWigs/142-125_ChIP_set1_H3K36me3_log2ratio.bw ../BigWigs/145-40_ChIP_set1E7_H3K36me3_Rep2_log2ratio.bw ../BigWigs/142-128_ChIP_sgr9_H3K36me3_log2ratio.bw ../BigWigs/145-115_ChIP_sgr9_H3K36me3_log2ratio.bw  --outFileName H3K36me3_log2ratiomatrix.gz

  #plotProfile --startLabel "TSS'" --endLabel "TES'" --averageType mean  --matrixFile ${OUTDIR}/Matrices/H3K36me3_log2ratiomatrix.gz --outFileName ${OUTDIR}/Heatmaps/H3K36me3_log2ratiomatrix_mean_V3.png  --outFileNameData ${OUTDIR}/Heatmaps/H3K36me3_log2ratiomatrix_mean.tab --numPlotsPerRow 2 --plotType=fill  --perGroup --legendLocation lower-right --samplesLabel "142-77-WT" "145-32-WT" "145-36-set7" "142-116-set2" "142-119-set2" "142-125-set1" "145-40-set1E7" "142-128-sgr9" "145-115-sgr9" --regionsLabel "H3K27me3 marked genes" "True rDNA genes" "Neighboring H3K27me3 genes" "Pseudo rDNA genes" "Ash1 H3K36me genes" --plotTitle "L2FC Input Normalized H3K36me3_COMPASS Mean profile" --yMin 0 --yMax 10 400 10 15 10

  #plotProfile --startLabel "TSS'" --endLabel "TES'" --averageType median  --matrixFile ${OUTDIR}/Matrices/H3K36me3_log2ratiomatrix.gz --outFileName ${OUTDIR}/Heatmaps/H3K36me3_log2ratiomatrix_median_V3.png  --outFileNameData ${OUTDIR}/Heatmaps/H3K36me3_log2ratiomatrix_median.tab --numPlotsPerRow 2 --plotType=fill  --perGroup --legendLocation lower-right --samplesLabel "142-77-WT" "145-32-WT" "145-36-set7" "142-116-set2" "142-119-set2" "142-125-set1" "145-40-set1E7" "142-128-sgr9" "145-115-sgr9" --regionsLabel "H3K27me3 marked genes" "True rDNA genes" "Neighboring H3K27me3 genes" "Pseudo rDNA genes" "Ash1 H3K36me genes" --plotTitle "L2FC Input Normalized H3K36me3_COMPASS Median profile" --yMin 0 --yMax 10 500 10 10 10

#  computeMatrix scale-regions  -p max --startLabel "TSS'" --endLabel "TES'" -b 1000 -a 1000  --regionsFileName /scratch/ry00555/heatmapPRC2genes.bed /scratch/ry00555/TRUErDNA_genes_with_names_fixed.bed /scratch/ry00555/neighboringK27genes.bed /scratch/ry00555/rDNA_regions_no_duplicates_genes_format.bed /scratch/ry00555/ash1depgenes.bed --skipZeros  --scoreFileName ${OUTDIR}/BigWigs/142-77_ChIP_WT_H3K36me3_log2ratio_smooth60.bw ${OUTDIR}/BigWigs/145-32_ChIP_S1_H3K36me3_log2ratio_smooth60.bw ${OUTDIR}/BigWigs/145-36_ChIP_set7_H3K36me3_log2ratio_smooth60.bw ${OUTDIR}/BigWigs/142-125_ChIP_set1_H3K36me3_log2ratio_smooth60.bw ${OUTDIR}/BigWigs/145-40_ChIP_set1E7_H3K36me3_Rep2_log2ratio_smooth60.bw ${OUTDIR}/BigWigs/142-128_ChIP_sgr9_H3K36me3_log2ratio_smooth60.bw ${OUTDIR}/BigWigs/145-115_ChIP_sgr9_H3K36me3_log2ratio_smooth60.bw  --outFileName ${OUTDIR}/Matrices/H3K36me3_log2ratiomatrix_NoSet2_smooth60.gz

#  plotProfile --startLabel "TSS'" --endLabel "TES'" --averageType mean  --matrixFile ${OUTDIR}/Matrices/H3K36me3_log2ratiomatrix_NoSet2_smooth60.gz --outFileName ${OUTDIR}/Heatmaps/H3K36me3_log2ratiomatrix_NoSet2_smooth60_V1_mean.png  --outFileNameData ${OUTDIR}/Heatmaps/H3K36me3_log2ratiomatrix_NoSet2_smooth60_V1_mean.tab --numPlotsPerRow 2 --plotType=fill  --perGroup --legendLocation lower-right --samplesLabel "142-77-WT" "145-32-WT" "145-36-set7" "142-125-set1" "145-40-set1E7" "142-128-sgr9" "145-115-sgr9" --regionsLabel "H3K27me3 marked genes" "True rDNA genes" "Neighboring H3K27me3 genes" "Pseudo rDNA genes" "Ash1 H3K36me genes" --plotTitle "L2FC Input Normalized H3K36me3_COMPASS Mean profile" --yMin 0

#  bamCompare -p max -b1 ${BAMDIR}/142-125_ChIP_set1_H3K36me3__Q30.bam	-b2	${BAMDIR}/142-123_ChIP_set1_Input__Q30.bam -of bigwig --scaleFactorsMethod None --effectiveGenomeSize 41037538 --skipZeroOverZero --binSize 30 --normalizeUsing CPM --smoothLength 60 -o ${OUTDIR}/BigWigs/142-125_ChIP_set1_H3K36me3_log2ratio_smooth60_CPM.bw

#bamCompare -p max -b1 ${BAMDIR}/145-40_ChIP_set1E7_H3K36me3_Rep2_Q30.bam -b2	${BAMDIR}/145-37_ChIP_set1E7_Input_Rep2_Q30.bam -of bigwig --scaleFactorsMethod None --effectiveGenomeSize 41037538 --skipZeroOverZero --binSize 30 --normalizeUsing CPM --smoothLength 60 -o ${OUTDIR}/BigWigs/145-40_ChIP_set1E7_H3K36me3_Rep2_log2ratio_smooth60_CPM.bw

#bamCompare -p max -b1 ${BAMDIR}/142-128_ChIP_sgr9_H3K36me3__Q30.bam	-b2	${BAMDIR}/142-126_ChIP_sgr9_Input__Q30.bam -of bigwig  --scaleFactorsMethod None --effectiveGenomeSize 41037538 --skipZeroOverZero --binSize 30 --normalizeUsing CPM --smoothLength 60 -o ${OUTDIR}/BigWigs/142-128_ChIP_sgr9_H3K36me3_log2ratio_smooth60_CPM.bw

#  bamCompare -p max -b1 ${BAMDIR}/142-77_ChIP_WT_H3K36me3__S77_L007_R1_001_val_1.fq.gz_Q30.bam	-b2	${BAMDIR}/142-75_ChIP_WT_Input__Q30.bam --scaleFactorsMethod None --effectiveGenomeSize 41037538 -of bigwig --skipZeroOverZero --binSize 30 --normalizeUsing CPM --smoothLength 60 -o ${OUTDIR}/BigWigs/142-77_ChIP_WT_H3K36me3_log2ratio_smooth60_CPM.bw

#  bamCompare -p max -b1 ${BAMDIR}/145-115_ChIP_sgr9_H3K36me3_Rep2_Q30.bam -b2 ${BAMDIR}/145-112_ChIP_sgr9_Input_Rep2_Q30.bam -of bigwig  --minMappingQuality 10  --scaleFactorsMethod None --effectiveGenomeSize 41037538 --skipZeroOverZero --binSize 30 --normalizeUsing CPM --smoothLength 60 -o ${OUTDIR}/BigWigs/145-115_ChIP_sgr9_H3K36me3_log2ratio_smooth60_CPM.bw

#bamCompare -p max -b1 ${BAMDIR}/145-32_ChIP_S1_H3K36me3_Rep2_Q30.bam	-b2	${BAMDIR}/145-29_ChIP_WT_Input_Q30.bam -of bigwig --minMappingQuality 10  --scaleFactorsMethod None --effectiveGenomeSize 41037538 --skipZeroOverZero --binSize 30 --normalizeUsing CPM --smoothLength 60 -o ${OUTDIR}/BigWigs/145-32_ChIP_S1_H3K36me3_log2ratio_smooth60_CPM.bw

#bamCompare -p max -b1 ${BAMDIR}/145-36_ChIP_set7_H3K36me3_Rep2_Q30.bam	-b2	${BAMDIR}/145-33_ChIP_set7_Input_Rep2_Q30.bam -of bigwig --minMappingQuality 10  --scaleFactorsMethod None --effectiveGenomeSize 41037538 --skipZeroOverZero --binSize 30 --normalizeUsing CPM --smoothLength 60 -o ${OUTDIR}/BigWigs/145-36_ChIP_set7_H3K36me3_log2ratio_smooth60_CPM.bw

#computeMatrix scale-regions  -p max --startLabel "TSS'" --endLabel "TES'" -b 1000 -a 1000  --regionsFileName /scratch/ry00555/heatmapPRC2genes.bed /scratch/ry00555/TRUErDNA_genes_with_names_fixed.bed /scratch/ry00555/neighboringK27genes.bed /scratch/ry00555/rDNA_regions_no_duplicates_genes_format.bed /scratch/ry00555/ash1depgenes.bed --skipZeros --scoreFileName ${OUTDIR}/BigWigs/142-77_ChIP_WT_H3K36me3_log2ratio_smooth60_CPM.bw ${OUTDIR}/BigWigs/145-32_ChIP_S1_H3K36me3_log2ratio_smooth60_CPM.bw ${OUTDIR}/BigWigs/145-36_ChIP_set7_H3K36me3_log2ratio_smooth60_CPM.bw ${OUTDIR}/BigWigs/142-125_ChIP_set1_H3K36me3_log2ratio_smooth60_CPM.bw ${OUTDIR}/BigWigs/145-40_ChIP_set1E7_H3K36me3_Rep2_log2ratio_smooth60_CPM.bw ${OUTDIR}/BigWigs/142-128_ChIP_sgr9_H3K36me3_log2ratio_smooth60_CPM.bw ${OUTDIR}/BigWigs/145-115_ChIP_sgr9_H3K36me3_log2ratio_smooth60_CPM.bw --outFileName ${OUTDIR}/Matrices/H3K36me3_log2ratiomatrix_NoSet2_smooth60_CPM.gz

#plotProfile --startLabel "TSS'" --endLabel "TES'" --averageType mean  --matrixFile ${OUTDIR}/Matrices/H3K36me3_log2ratiomatrix_NoSet2_smooth60_CPM.gz --outFileName ${OUTDIR}/Heatmaps/H3K36me3_log2ratiomatrix_NoSet2_smooth60_V1_mean_CPM.png  --outFileNameData ${OUTDIR}/Heatmaps/H3K36me3_log2ratiomatrix_NoSet2_smooth60_V1_mean_CPM.tab --numPlotsPerRow 2 --plotType=fill  --perGroup --legendLocation lower-right --samplesLabel "142-77-WT" "145-32-WT" "145-36-set7" "142-125-set1" "145-40-set1E7" "142-128-sgr9" "145-115-sgr9" --regionsLabel "H3K27me3 marked genes" "True rDNA genes" "Neighboring H3K27me3 genes" "Pseudo rDNA genes" "Ash1 H3K36me genes" --plotTitle "L2FC Input Normalized H3K36me3_COMPASS Mean profile" --yMin 0

#bamCompare -p max -b1 ${BAMDIR}/142-125_ChIP_set1_H3K36me3__Q30.bam	-b2	${BAMDIR}/142-123_ChIP_set1_Input__Q30.bam -of bigwig --skipZeroOverZero --binSize 30 --smoothLength 60 -o ${OUTDIR}/BigWigs/142-125_ChIP_set1_H3K36me3_log2ratio_smooth60_readcount.bw

#bamCompare -p max -b1 ${BAMDIR}/145-40_ChIP_set1E7_H3K36me3_Rep2_Q30.bam -b2	${BAMDIR}/145-37_ChIP_set1E7_Input_Rep2_Q30.bam -of bigwig --skipZeroOverZero --binSize 30 --smoothLength 60 -o ${OUTDIR}/BigWigs/145-40_ChIP_set1E7_H3K36me3_Rep2_log2ratio_smooth60_readcount.bw

#bamCompare -p max -b1 ${BAMDIR}/142-128_ChIP_sgr9_H3K36me3__Q30.bam	-b2	${BAMDIR}/142-126_ChIP_sgr9_Input__Q30.bam -of bigwig --effectiveGenomeSize 41037538 --skipZeroOverZero --binSize 30 --smoothLength 60 -o ${OUTDIR}/BigWigs/142-128_ChIP_sgr9_H3K36me3_log2ratio_smooth60_readcount.bw

#bamCompare -p max -b1 ${BAMDIR}/142-77_ChIP_WT_H3K36me3__S77_L007_R1_001_val_1.fq.gz_Q30.bam	-b2	${BAMDIR}/142-75_ChIP_WT_Input__Q30.bam -of bigwig --skipZeroOverZero --binSize 30 --smoothLength 60 -o ${OUTDIR}/BigWigs/142-77_ChIP_WT_H3K36me3_log2ratio_smooth60_readcount.bw

#bamCompare -p max -b1 ${BAMDIR}/145-115_ChIP_sgr9_H3K36me3_Rep2_Q30.bam -b2 ${BAMDIR}/145-112_ChIP_sgr9_Input_Rep2_Q30.bam -of bigwig  --minMappingQuality 10  --skipZeroOverZero --binSize 30 --smoothLength 60 -o ${OUTDIR}/BigWigs/145-115_ChIP_sgr9_H3K36me3_log2ratio_smooth60_readcount.bw

#bamCompare -p max -b1 ${BAMDIR}/145-32_ChIP_S1_H3K36me3_Rep2_Q30.bam	-b2	${BAMDIR}/145-29_ChIP_WT_Input_Q30.bam -of bigwig --minMappingQuality 10 --skipZeroOverZero --binSize 30 --smoothLength 60 -o ${OUTDIR}/BigWigs/145-32_ChIP_S1_H3K36me3_log2ratio_smooth60_readcount.bw

#bamCompare -p max -b1 ${BAMDIR}/145-36_ChIP_set7_H3K36me3_Rep2_Q30.bam	-b2	${BAMDIR}/145-33_ChIP_set7_Input_Rep2_Q30.bam -of bigwig --minMappingQuality 10 --skipZeroOverZero --binSize 30 --smoothLength 60 -o ${OUTDIR}/BigWigs/145-36_ChIP_set7_H3K36me3_log2ratio_smooth60_readcount.bw

#computeMatrix scale-regions  -p max --startLabel "TSS'" --endLabel "TES'" -b 1000 -a 1000  --regionsFileName /scratch/ry00555/heatmapPRC2genes.bed /scratch/ry00555/TRUErDNA_genes_with_names_fixed.bed /scratch/ry00555/neighboringK27genes.bed /scratch/ry00555/rDNA_regions_no_duplicates_genes_format.bed /scratch/ry00555/ash1depgenes.bed --skipZeros --scoreFileName ${OUTDIR}/BigWigs/142-77_ChIP_WT_H3K36me3_log2ratio_smooth60_readcount.bw ${OUTDIR}/BigWigs/145-32_ChIP_S1_H3K36me3_log2ratio_smooth60_readcount.bw ${OUTDIR}/BigWigs/145-36_ChIP_set7_H3K36me3_log2ratio_smooth60_readcount.bw ${OUTDIR}/BigWigs/142-125_ChIP_set1_H3K36me3_log2ratio_smooth60_readcount.bw ${OUTDIR}/BigWigs/145-40_ChIP_set1E7_H3K36me3_Rep2_log2ratio_smooth60_readcount.bw ${OUTDIR}/BigWigs/142-128_ChIP_sgr9_H3K36me3_log2ratio_smooth60_readcount.bw ${OUTDIR}/BigWigs/145-115_ChIP_sgr9_H3K36me3_log2ratio_smooth60_readcount.bw --outFileName ${OUTDIR}/Matrices/H3K36me3_log2ratiomatrix_NoSet2_smooth60_readcount.gz

#plotProfile --startLabel "TSS'" --endLabel "TES'" --averageType mean  --matrixFile ${OUTDIR}/Matrices/H3K36me3_log2ratiomatrix_NoSet2_smooth60_readcount.gz --outFileName ${OUTDIR}/Heatmaps/H3K36me3_log2ratiomatrix_NoSet2_smooth60_V1_mean_readcount.png  --outFileNameData ${OUTDIR}/Heatmaps/H3K36me3_log2ratiomatrix_NoSet2_smooth60_V1_mean_readcount.tab --numPlotsPerRow 2 --plotType=fill  --perGroup --legendLocation lower-right --samplesLabel "142-77-WT" "145-32-WT" "145-36-set7" "142-125-set1" "145-40-set1E7" "142-128-sgr9" "145-115-sgr9" --regionsLabel "H3K27me3 marked genes" "True rDNA genes" "Neighboring H3K27me3 genes" "Pseudo rDNA genes" "Ash1 H3K36me genes" --plotTitle "L2FC Input Normalized H3K36me3_COMPASS Mean profile" --yMin 0

plotProfile --startLabel "TSS'" --endLabel "TES'" --averageType sum --matrixFile ${OUTDIR}/Matrices/H3K36me3_log2ratiomatrix_NoSet2_smooth60_readcount.gz --outFileName ${OUTDIR}/Heatmaps/H3K36me3_log2ratiomatrix_NoSet2_smooth60_V1_sum_readcount_V2.png  --outFileNameData ${OUTDIR}/Heatmaps/H3K36me3_log2ratiomatrix_NoSet2_smooth60_V1_sum_readcount_V2.tab --numPlotsPerRow 2 --plotType=fill  --perGroup --legendLocation lower-right --samplesLabel "142-77-WT" "145-32-WT" "145-36-set7" "142-125-set1" "145-40-set1E7" "142-128-sgr9" "145-115-sgr9" --regionsLabel "H3K27me3 marked genes" "True rDNA genes" "Neighboring H3K27me3 genes" "Pseudo rDNA genes" "Ash1 H3K36me genes" --plotTitle "L2FC Input Normalized H3K36me3_COMPASS Mean profile" --yMin 0 --yMax 1000 5 300 200 300
