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

OUTDIR= "/scratch/ry00555/OutputRun137/CutandRun"
FASTQ= "/scratch/ry00555/OutputRun137/CutandRun"

ml Trim_Galore
#starting with raw files in $OUTDIR/raw
mkdir "$OUTDIR/TrimmedReads"
trim_galore --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${FASTQ}/*fastq\.gz

FILES="${OUTDIR}/TrimmedReads/*R1_001_val_1\.fq\.gz" #Don't forget the *

mkdir "$OUTDIR/ref"
curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz | gunzip -c > $OUTDIR/ref/ecoli_refseq.fa
curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/925/GCF_000182925.2_NC12/GCF_000182925.2_NC12_genomic.fna.gz | gunzip -c > $OUTDIR/ref/Ncrassa_refseq.fa

module load Bowtie2
bowtie2-build -f $OUTDIR/ref/Ncrassa_refseq.fa $OUTDIR/ref/Ncrassa_ref
file=${f##*/}
# 	#remove ending from file name to create shorter names for bam files and other downstream output
name=${file/%_S[1-12]*_R1_001_val_1.fq.gz/}
#
# #
# # 	# File Vars
# # 	#use sed to get the name of the second read matching the input file
read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R2_001_val_2\.fq\.gz/g')
mkdir $OUTDIR/sam_files
bowtie2 --local --very-sensitive-local --phred33 --no-unal -p 24 -x $OUTDIR/ref/Ncrassa_ref -1 $f -2 $read2 -S $OUTDIR/sam_files/${name}.sam

bowtie2-build -f $OUTDIR/ref/ecoli_refseq.fa $OUTDIR/ref/ecoli_ref
bowtie2 --local --very-sensitive-local --phred33 --no-unal -p 24 -x $OUTDIR/ref/ecoli_ref -1 $f -2 $read2 -S $OUTDIR/sam_files/${name}_Ecoli.sam

module load SAMtools
mkdir "$OUTDIR/bam_files"
samtools view -bS -h $OUTDIR/sam_files/${name}.sam  > $OUTDIR/bam_files/${name}.bam
samtools view -bS -h $OUTDIR/sam_files/${name}_Ecoli.sam  > $OUTDIR/bam_files/${name}_Ecoli.bam

mkdir "$OUTDIR/SortedBamFiles"
samtools sort $OUTDIR/bam_files/${name}.bam -o $OUTDIR/SortedBamFiles/${name}.sorted.bam
samtools sort $OUTDIR/bam_files/${name}_Ecoli.bam -o $OUTDIR/SortedBamFiles/${name}_Ecoli.sorted.bam

#ml SAMtools
# samtools merge $OUTDIR/bam_files/timecourse_IgG_danio_merged.bam $OUTDIR/bam_files/2.5hpf_IgG.sorted.bam $OUTDIR/bam_files/4.5hpf_IgG.sorted.bam $OUTDIR/bam_files/24hpf_IgG.sorted.bam
#samtools merge
