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

 #FILES="${OUTDIR}/TrimmedReads/*_R1_001_val_1\.fq\.gz"
#
# for f in $FILES
# do
# file=${f##*/}
#name=${file/%_S[1-990]*_R1_001_val_1.fq.gz/}
# read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R2_001_val_2\.fq\.gz/g')
#
# bam="${OUTDIR}/SortedBamFiles/${name}.bam"
# variable name for bigwig output
 #	bigwig="${OUTDIR}/BigWigs/${name}"
 #QualityBam="${OUTDIR}/SortedBamFiles/${name}_Q30.bam"
#
 #ml SAMtools/1.16.1-GCC-11.3.0
 #ml BWA/0.7.17-GCCcore-11.3.0
#
# bwa mem -M -v 3 -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIR/SortedBamFiles/tempReps -o "$bam" -
# samtools index "$bam"
#
# samtools view -b -q 30 $bam > "$QualityBam"
# samtools index "$QualityBam"
#
 #ml deepTools/3.5.2-foss-2022a
 #bamCoverage -p $THREADS -bs 10 --normalizeUsing BPM --minMappingQuality 10 --smoothLength $SMOOTH -of bigwig -b "$QualityBam" -o "${bigwig}.bin_10.smooth_${SMOOTH}Bulk.bw"
#
 #done

#ml deepTools
#computeMatrix scale-regions -p 12 -R /scratch/ry00555/neurospora.bed -S ${OUTDIR}/BigWigs/


#mkdir $OUTDIR/MACSPeaks


#module load MACS3/3.0.0b1-foss-2022a-Python-3.10.4
 #command line
# macs3 callpeak -t $BAMDIR/142-94_ChIP_set7_H3K27me3__Q30.bam -f BAMPE -n 142-94_ChIP_set7_H3K27me3 -c $BAMDIR/142-93_ChIP_set7_Input__Q30.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
# macs3 callpeak -t $BAMDIR/142-106_ChIP_swd1_H3K27me3__Q30.bam -f BAMPE -n 142-106_ChIP_swd1_H3K27me3 -c $BAMDIR/142-105_ChIP_swd1_Input__Q30.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
# macs3 callpeak -t $BAMDIR/142-121_ChIP_set1_H3K27me3__Q30.bam -f BAMPE -n 142-121_ChIP_set1_H3K27me3 -c $BAMDIR/142-123_ChIP_set1_Input__Q30.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
# macs3 callpeak -t $BAMDIR/142-124_ChIP_set1_H3K27me3__Q30.bam -f BAMPE -n 142-124_ChIP_set1_H3K27me3 -c $BAMDIR/142-123_ChIP_set1_Input__Q30.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
# macs3 callpeak -t $BAMDIR/142-127_ChIP_sgr9_H3K27me3__Q30.bam  -f BAMPE -n 142-127_ChIP_sgr9_H3K27me3 -c $BAMDIR/142-126_ChIP_sgr9_Input__Q30.bam  --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
# macs3 callpeak -t $BAMDIR/142-10_ChIP_WT_H3K27me3_Rep3_Q30.bam  -f BAMPE -n 142-10_ChIP_WT_H3K27me3 -c $BAMDIR/142-75_ChIP_WT_Input__Q30.bam  --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
# macs3 callpeak -t $BAMDIR/142-115_ChIP_set2_H3K27me3__Q30.bam -f BAMPE -n 142-115_ChIP_set2_H3K27me3  --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
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


ml BEDTools
ml UCSC
#bedGraphToBigWig [options] in.bedGraph chrom.sizes out.bw
module load ucsc/359
 for infile in $PEAKDIR/BedGraphs/*bedgraph
do
  base=$(basename ${infile} .bedgraph)
  bedSort $infile $PEAKDIR/BedGraphs/${base}.sort.bedgraph
 done

#mkdir $OUTDIR/bigwigs
ml deepTools
for infile in $PEAKDIR/BedGraphs/*.sort.bedgraph
 do
  base=$(basename ${infile} .kmet_sort.bga)
 bedGraphToBigWig $infile /scratch/ry00555/Run137CutandRun/ref/Ncrassa_ref/chrNameLength.txt $PEAKDIR/BedGraphs/${base}.bw
 done

#bedtools intersect -wa -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GenesOnly.bed -b $PEAKDIR/Intersected_WT_H3K27me3_all.bed > $BEDDIR/MACS_WT_IntersectedH3K27me3_peaks.bed
#bedtools intersect -wa -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GenesOnly.bed -b $PEAKDIR/Intersected_set2_H3K27me3_all.bed > $BEDDIR/MACS_set2_IntersectedH3K27me3_peaks.bed
#bedtools intersect -wa -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GenesOnly.bed -b $PEAKDIR/Intersected_set1_H3K27me3_all.bed > $BEDDIR/MACS_set1_IntersectedH3K27me3_peaks.bed
#bedtools intersect -wa -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GenesOnly.bed -b $PEAKDIR/Intersected_swd1_H3K27me3_all.bed > $BEDDIR/MACS_swd1_IntersectedH3K27me3_peaks.bed
#bedtools intersect -wa -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GenesOnly.bed -b $PEAKDIR/Intersected_sgr9_H3K27me3_all.bed > $BEDDIR/MACS_sgr9_IntersectedH3K27me3_peaks.bed
#bedtools intersect -wa -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GenesOnly.bed -b $PEAKDIR/Intersected_set7_H3K27me3_all.bed > $BEDDIR/MACS_set7_IntersectedH3K27me3_peaks.bed


#mkdir ${OUTDIR}/HomerPeaks
HOMERPEAKSDIR="${OUTDIR}/HomerPeaks"
ml Homer
ml Perl

#for bam_file in "${BAMDIR}"/145*_Q30.bam; do
#  sample_id=$(basename "${bam_file}" _Q30.bam)
#makeTagDirectory "${TAGDIR}/${sample_id}" "${bam_file}"
#done

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
Genome="/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_00182925.2plusHphplusBarplusTetO_his3masked.fna"
GTF="/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_WithExtras_GFFtoGTFconversion.gtf"
#GFF="/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic.gff"
 #for infile in ${HOMERPEAKSDIR}/*_Homerpeaks.txt
#do
  #base=$(basename ${infile} _Homerpeaks.txt)
#  sed '/^#/d' $infile | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $8 "\t" $5 "\t" $6 "\t" $12 "\t" "-1"}' | sed 's/\.000000//g' > ${HOMERPEAKSDIR}/${base}.peaks.bed
 #annotatePeaks.pl ${HOMERPEAKSDIR}/${base}.peaks.bed $Genome -gff3 $GFF3 > ${HOMERPEAKSDIR}/${base}_ann.txt

#done

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
