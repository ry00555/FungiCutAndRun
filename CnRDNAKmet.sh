#!/bin/bash
#SBATCH --job-name=CUTandRun137
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=../MapCutAndRun132.%j.out
#SBATCH --error=../MapCutAndRun132.%j.err
#cd $SLURM_SUBMIT_DIR
#$OUTDIR="/scratch/ry00555/OutputRun137/CutandRun"
#$FASTQ="/scratch/ry00555/OutputRun137/CutandRun"

# ml Trim_Galore
# #mkdir "$OUTDIR/TrimmedReads"
# trim_galore --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${OUTDIR}/*fastq\.gz
# #in line commands
# trim_galore --paired --length 20 --fastqc --gzip -o TrimmedReads *fastq\.gz

FILES="/scratch/ry00555/OutputRun137/CutandRun/TrimmedReads/*R1_001_val_1\.fq\.gz" #Don't forget the *

#mkdir "$OUTDIR/ref"
# curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz | gunzip -c > $OUTDIR/ref/ecoli_refseq.fa
# curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/925/GCF_000182925.2_NC12/GCF_000182925.2_NC12_genomic.fna.gz | gunzip -c > $OUTDIR/ref/Ncrassa_refseq.fa
# #in line commands
# curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz | gunzip -c > ref/ecoli_refseq.fa
# curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/925/GCF_000182925.2_NC12/GCF_000182925.2_NC12_genomic.fna.gz | gunzip -c > ref/Ncrassa_refseq.fa

#module load Bowtie2
# bowtie2-build -f $OUTDIR/ref/Ncrassa_refseq.fa $OUTDIR/ref/Ncrassa_ref
# #in line commands
# bowtie2-build -f ref/Ncrassa_refseq.fa ref/Ncrassa_ref
# bowtie2-build -f ref/ecoli_refseq.fa ref/Ecoli_ref

#Use STAR to create chrNameLength.txt files
# ml STAR
# STAR --runThreadN 20 --genomeSAindexNbases 8 --runMode genomeGenerate --genomeDir $OUTDIR/ref/Ncrassa_ref --genomeFastaFiles $OUTDIR/ref/Ncrassa_refseq.fa
# #command line
# STAR --runThreadN 20 --genomeSAindexNbases 8 --runMode genomeGenerate --genomeFastaFiles Ncrassa_refseq.fa

#for f in $FILES
  #do
 #file=${f##*/}
# # 	#remove ending from file name to create shorter names for bam files and other downstream output
 #name=${file/%_S[1-12]*_R1_001_val_1.fq.gz/}
# #
# # #
# # # 	# File Vars
# # # 	#use sed to get the name of the second read matching the input file
 #read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R2_001_val_2\.fq\.gz/g')
# #mkdir $OUTDIR/sam_files
# bowtie2 --local --very-sensitive-local --phred33 --no-unal -p 24 -x "/scratch/ry00555/OutputRun137/CutandRun/ref/Ncrassa_ref" -1 $f -2 $read2 -S "/scratch/ry00555/OutputRun137/CutandRun/sam_files/${name}.sam"
#
# #bowtie2-build -f ${OUTDIR}/ref/ecoli_refseq.fa $OUTDIR/ref/ecoli_ref
# bowtie2 --local --very-sensitive-local --phred33 --no-unal -p 24 -x "/scratch/ry00555/OutputRun137/CutandRun/ref/Ecoli_ref" -1 $f -2 $read2 -S "/scratch/ry00555/OutputRun137/CutandRun/sam_files/${name}_Ecoli.sam"
#
# module load SAMtools
# #mkdir "$OUTDIR/bam_files"
# samtools view -bS -h "/scratch/ry00555/OutputRun137/CutandRun/sam_files/${name}.sam"  > "/scratch/ry00555/OutputRun137/CutandRun/bam_files/${name}.bam"
# samtools view -bS -h "/scratch/ry00555/OutputRun137/CutandRun/sam_files/${name}_Ecoli.sam"  > "/scratch/ry00555/OutputRun137/CutandRun/bam_files/${name}_Ecoli.bam"
#
# #mkdir "$OUTDIR/SortedBamFiles"
# samtools sort "/scratch/ry00555/OutputRun137/CutandRun/bam_files/${name}.bam" -o "/scratch/ry00555/OutputRun137/CutandRun/SortedBamFiles/${name}.sorted.bam"
# samtools sort "/scratch/ry00555/OutputRun137/CutandRun/bam_files/${name}_Ecoli.bam" -o "/scratch/ry00555/OutputRun137/CutandRun/SortedBamFiles/${name}_Ecoli.sorted.bam"
#
# done
#no need to samtools merge at the moment because I only have one of each sample
##turning sorted bam files into bed graphs for DNA spike in
#mkdir $OUTDIR/bed_files
# for f in $FILES
#  do
# file=${f##*/}
# # 	#remove ending from file name to create shorter names for bam files and other downstream output
# name=${file/%_S[1-12]*_R1_001_val_1.fq.gz/}
#
# ml BEDTools
# bedtools bamtobed -i /scratch/ry00555/OutputRun137/CutandRun/SortedBamFiles/${name}.sorted.bam | awk -v OFS='\t' '{len = $3 - $2; print $0, len }' > /scratch/ry00555/OutputRun137/CutandRun/bed_files/${name}.btb.bed
# bedtools bamtobed -i /scratch/ry00555/OutputRun137/CutandRun/SortedBamFiles/${name}_Ecoli.sorted.bam | awk -v OFS='\t' '{len = $3 - $2; print $0, len }' > /scratch/ry00555/OutputRun137/CutandRun/bed_files/${name}_Ecoli.btb.bed
#
# #DNA-spike in normalization
# #mkdir $OUTDIR/bedgraphs
#sh /home/ry00555/Research/FungiCutAndRun/CUTandRUNAnalysis/DNAspike_in.kd.sh /scratch/ry00555/OutputRun137/CutandRun/bed_files/${name}.bed /scratch/ry00555/OutputRun137/CutandRun/bed_files/${name}_Ecoli.btb.bed 100000 bga "/scratch/ry00555/OutputRun137/CutandRun/ref/GenomeDir/chrNameLength.txt" 1 1000 /scratch/ry00555/OutputRun137/CutandRun/bedgraphs/${name}.norm.bga
# done
#sort bga files from  DNA spike in
 # ml ucsc
 # for infile in /scratch/ry00555/OutputRun137/CutandRun/bedgraphs/*norm.bga
 #  do
 #    base=$(basename ${infile} .norm.bga)
 #    bedSort $infile /scratch/ry00555/OutputRun137/CutandRun/bedgraphs/${base}.norm_sort.bga
 # done

#module load Homer
 #calling peaks
 # mkdir $OUTDIR/Peaks
#   for infile in /scratch/ry00555/OutputRun137/CutandRun/bedgraphs/*.norm_sort.bga
#     do base=$(basename ${infile} .norm_sort.bga)
# cat $infile | awk '{print $1 "\t" $2 "\t" $3 "\t" "+" "\t" "+" "\t" "+"}' > /scratch/ry00555/OutputRun137/CutandRun/Peaks/${base}.bgato.bed
#  done
#
#   for infile in /scratch/ry00555/OutputRun137/CutandRun/Peaks/*bgato.bed
#    do base=$(basename ${infile} .bgato.bed)
#     makeTagDirectory /scratch/ry00555/OutputRun137/CutandRun/TagDirectories/${base}.BtB.tagdir $infile -format bed
#  done

#!/bin/bash

# Set the directory containing the tag directories
OUTDIR="/scratch/ry00555/OutputRun137/CutandRun"
# # this code works where I have to type in strain  so do it for the rest change rtt109 to WT , set-7, ncu00423, ncu006787, ncu06788
#
# ##using IgG as input
# for infile in $OUTDIR/TagDirectories/*WT*.tagdir
# do
#   base=$(basename ${infile} .tagdir)
# findPeaks $infile -style histone -minDist 1000 -i $OUTDIR/TagDirectories/137-1_CUTANDRUN_WT_IgG_Rep1.BtB.tagdir -F 4 -gsize 4.5e7 -o $OUTDIR/Peaks/${base}_IgGNorm.txt
# done
#
# ##using IgG as input
# for infile in $OUTDIR/TagDirectories/*set-7*.tagdir
# do
#   base=$(basename ${infile} .tagdir)
# findPeaks $infile -style histone -minDist 1000 -i $OUTDIR/TagDirectories/137-6_CUTANDRUN_set-7_IgG_Rep1_S6_R1_001_val_1.fq.gz.BtB.tagdir -F 4 -gsize 4.5e7 -o $OUTDIR/Peaks/${base}_IgGNorm.txt
# done
# ##using IgG as input
# for infile in $OUTDIR/TagDirectories/*ncu00423*.tagdir
# do
#   base=$(basename ${infile} .tagdir)
# findPeaks $infile -style histone -minDist 1000 -i $OUTDIR/TagDirectories/137-18_CUTANDRUN_ncu00423_IgG_Rep1.BtB.tagdir -F 4 -gsize 4.5e7 -o $OUTDIR/Peaks/${base}_IgGNorm.txt
# done
# ##using IgG as input
# for infile in $OUTDIR/TagDirectories/*ncu06787*.tagdir
# do
#   base=$(basename ${infile} .tagdir)
# findPeaks $infile -style histone -minDist 1000 -i $OUTDIR/TagDirectories/137-12_CUTANDRUN_ncu06787_IgG_Rep1.BtB.tagdir -F 4 -gsize 4.5e7 -o $OUTDIR/Peaks/${base}_IgGNorm.txt
# done
# ##using IgG as input
# for infile in $OUTDIR/TagDirectories/*ncu06788*.tagdir
# do
#   base=$(basename ${infile} .tagdir)
# findPeaks $infile -style histone -minDist 1000 -i $OUTDIR/TagDirectories/137-15_CUTANDRUN_ncu06788_IgG_Rep1.BtB.tagdir -F 4 -gsize 4.5e7 -o $OUTDIR/Peaks/${base}_IgGNorm.txt
# done

#changing peak txt files to bed files to input into chipr
# for infile in $OUTDIR/Peaks/*_IgGNorm.txt
# do
# base=$(basename ${infile} _IgGNorm.txt)
# sed '/^#/d' $infile | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $8 "\t" $5 "\t" $6 "\t" $12 "\t" "-1"}' | sed 's/\.000000//g' > $OUTDIR/Peaks/${base}.peaks_IgGNorm.bed
#  done

 # don't need ChipR right now since I don't have replicates as of 4/10/24

 ##annotating peak files with masked reference (use HOMER module)
#  ml Perl
#  for infile in $OUTDIR/Peaks/*.peaks_IgGNorm.bed
#  do
#    base=$(basename ${infile} .peaks_IgGNorm.bed)
# annotatePeaks.pl $OUTDIR/Peaks/${base}.peaks_IgGNorm.bed /scratch/ry00555/OutputRun137/CutandRun/ref/Ncrassa_refseq.fa -gtf /scratch/ry00555/Ncrassa.gtf > $OUTDIR/Peaks/${base}.peaks_IgGNorm_ann.txt
# # you can analyze the peaks in excel now lets turn this into big wigs so we can make meta plots
# done

# ml ucsc
# for infile in $OUTDIR/bedgraphs/*.norm_sort.bga
# do
# base=$(basename ${infile} .norm_sort.bga)
# bedGraphToBigWig $infile $OUTDIR/ref/GenomeDir/chrNameLength.txt $OUTDIR/BigWigs/${base}_DNASpikeinNorm.bw
# done
#ml deepTools
#computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S ${OUTDIR}/BigWigs/137-2_CUTANDRUN_WT_H3K27me3_Rep1_DNASpikeinNorm.bw ${OUTDIR}/BigWigs/137-25_CUTANDRUN_set-7_H3K27me3_Rep1_DNASpikeinNorm.bw ${OUTDIR}/BigWigs/137-10_CUTANDRUN_rtt109_H3K27me3_Rep1_DNASpikeinNorm.bw ${OUTDIR}/BigWigs/137-19_CUTANDRUN_ncu00423_H3K27me3_Rep1_DNASpikeinNorm.bw ${OUTDIR}/BigWigs/137-13_CUTANDRUN_ncu06787_H3K27me3_Rep1_DNASpikeinNorm.bw ${OUTDIR}/BigWigs/137-16_CUTANDRUN_ncu06788_H3K27me3_Rep1_DNASpikeinNorm.bw -R "/scratch/ry00555/neurospora.bed" --skipZeros -o "${OUTDIR}/Matrices/matrix_CnR_H3K27me3.gz"
## command line

#computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S ../BigWigs/137-2_CUTANDRUN_WT_H3K27me3_Rep1_DNASpikeinNorm.bw ../BigWigs/137-25_CUTANDRUN_set-7_H3K27me3_Rep1_DNASpikeinNorm.bw ../BigWigs/137-10_CUTANDRUN_rtt109_H3K27me3_Rep1_DNASpikeinNorm.bw ../BigWigs/137-19_CUTANDRUN_ncu00423_H3K27me3_Rep1_DNASpikeinNorm.bw ../BigWigs/137-13_CUTANDRUN_ncu06787_H3K27me3_Rep1_DNASpikeinNorm.bw ../BigWigs/137-16_CUTANDRUN_ncu06788_H3K27me3_Rep1_DNASpikeinNorm.bw -R "/scratch/ry00555/Ncrassagenes.bed" --skipZeros -o matrix_CnR_H3K27me3.gz

#plotHeatmap -m "${OUTDIR}/Matrices/matrix_CnR_H3K27me3.gz" -out ${OUTDIR}/Heatmaps/CnR_H3K27me3_wholegenome_hclust.png --samplesLabel WT set-7 rtt109 ncu00423 ncu06787 ncu06788 --hclust 2 --colorMap Reds
#command line below too launch in Heatmaps directory
#plotHeatmap -m ../Matrices/matrix_CnR_H3K27me3.gz -out CnR_H3K27me3_wholegenome_hclust.png --samplesLabel WT set-7 rtt109 ncu00423 ncu06787 ncu06788 --hclust 2 --colorMap Reds
#
#scp -r ry00555@xfer.gacrc.uga.edu:CnR_H3K27me3_wholegenome_hclust.png /Users/ry00555/Desktop/

# computeMatrix reference-point --referencePoint TSS -b 1500 -a 1500 -S "${file_path}" -R "/scratch/ry00555/neurospora.bed" --skipZeros -o "${OUTDIR}/Matrices/matrix_CnR_H3K36me3.gz"
#

#kmet spike in for both reg sorted and ecoli sorted
for file in $OUTDIR/SortedBamFiles/*sorted.bam
 do
    base=$(basename "${file}" .sorted.bam)
 sh /home/ry00555/Research/FungiCutAndRun/CUTandRUNAnalysis/kmet_spike.sh $OUTDIR/KmetSpikeIn/bedgraphs $base $OUTDIR/TrimmedReads/${base}*R1_001_val_1.fq.gz \
$OUTDIR/TrimmedReads/${base}*R2_001_val_2.fq.gz $file bga $OUTDIR/ref/GenomeDir/chrNameLength.txt
done
#
# #sort bga files from spike in
#  ml ucsc
# for infile in $OUTDIR/KmetSpikeIn/bedgraphs/*kmet.bga
#  do
#  base=$(basename ${infile} _kmet.bga)
#      bedSort $infile $OUTDIR/KmetSpikeIn/bedgraphs/${base}.kmet_sort.bga
#    done

  #mkdir $OUTDIR/KmetSpikeIn/Peaks
 # for infile in $OUTDIR/KmetSpikeIn/bedgraphs/*kmet_sort.bga
 #   do base=$(basename ${infile} .kmet_sort.bga)
 #   cat $infile | awk '{print $1 "\t" $2 "\t" $3 "\t" "+" "\t" "+" "\t" "+"}' > $OUTDIR/KmetSpikeIn/Peaks/${base}.bgato.bed
 # done

 module load Homer
# mkdir $OUTDIR/KmetSpikeIn/TagDirectories
#  for infile in $OUTDIR/KmetSpikeIn/Peaks/*bgato.bed
#  do
# base=$(basename ${infile} .bgato.bed)
# makeTagDirectory $OUTDIR/KmetSpikeIn/TagDirectories/${base}.BtB.tagdir $infile -format bed
# done

# Do the same IgG analysis above for the KmetSpikeIn normalized for each strain at a time, then skip Chip-R and go straight to annotating peaks with input bed files use the same command as above and make sure to ml Perl then make bigwigs from the bga files with the end .kmet.sort.bga
##using IgG as input
# for infile in $OUTDIR/KmetSpikeIn/TagDirectories/*WT*.tagdir
# do
#  base=$(basename ${infile} .tagdir)
#  findPeaks $infile -style histone -minDist 1000 -i $OUTDIR/KmetSpikeIn/TagDirectories/137-1_CUTANDRUN_WT_IgG_Rep1.BtB.tagdir -F 4 -gsize 4.5e7 -o $OUTDIR/KmetSpikeIn/Peaks/${base}_IgGNorm.txt
#  done

 # for infile in $OUTDIR/KmetSpikeIn/TagDirectories/*rtt109*.tagdir
 # do
 #  base=$(basename ${infile} .tagdir)
 #  findPeaks $infile -style histone -minDist 1000 -i $OUTDIR/KmetSpikeIn/TagDirectories/137-1_CUTANDRUN_WT_IgG_Rep1.BtB.tagdir -F 4 -gsize 4.5e7 -o $OUTDIR/KmetSpikeIn/Peaks/${base}_IgGNorm.txt
 #  done
#set-7 is also missing ../
 # for infile in $OUTDIR/KmetSpikeIn/TagDirectories/*set-7*.tagdir
 # do
 #  base=$(basename ${infile} .tagdir)
 # findPeaks $infile -style histone -minDist 1000 -i $OUTDIR/KmetSpikeIn/TagDirectories/137-1_CUTANDRUN_WT_IgG_Rep1.BtB.tagdir -F 4 -gsize 4.5e7 -o $OUTDIR/KmetSpikeIn/Peaks/${base}_IgGNorm.txt
 # done
