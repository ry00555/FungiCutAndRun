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

# for f in $FILES
#  do
# file=${f##*/}
# # 	#remove ending from file name to create shorter names for bam files and other downstream output
# name=${file/%_S[1-12]*_R1_001_val_1.fq.gz/}
# #
# # #
# # # 	# File Vars
# # # 	#use sed to get the name of the second read matching the input file
# read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R2_001_val_2\.fq\.gz/g')
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
for f in $FILES
 do
file=${f##*/}
# 	#remove ending from file name to create shorter names for bam files and other downstream output
name=${file/%_S[1-12]*_R1_001_val_1.fq.gz/}

ml BEDTools
bedtools bamtobed -i /scratch/ry00555/OutputRun137/CutandRun/SortedBamFiles/${name}.sorted.bam | awk -v OFS='\t' '{len = $3 - $2; print $0, len }' > /scratch/ry00555/OutputRun137/CutandRun/bed_files/${name}.btb.bed
bedtools bamtobed -i /scratch/ry00555/OutputRun137/CutandRun/SortedBamFiles/${name}_Ecoli.sorted.bam | awk -v OFS='\t' '{len = $3 - $2; print $0, len }' > /scratch/ry00555/OutputRun137/CutandRun/bed_files/${name}_Ecoli.btb.bed

#DNA-spike in normalization
#mkdir $OUTDIR/bedgraphs
sh /home/ry00555/Research/FungiCutAndRun/CUTandRUNAnalysis/DNAspike_in.kd.sh /scratch/ry00555/OutputRun137/CutandRun/bed_files/${name}.bed /scratch/ry00555/OutputRun137/CutandRun/bed_files/${name}_Ecoli.btb.bed 100000 bga "/scratch/ry00555/OutputRun137/CutandRun/ref/GenomeDir/chrNameLength.txt" 1 1000 /scratch/ry00555/OutputRun137/CutandRun/bedgraphs/${name}.norm.bga
done
#sort bga files from  DNA spike in
 #ml ucsc
