#!/bin/bash
#SBATCH --job-name=j_GATK
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --time=08:00:00
#SBATCH --output=GATK.%j.out
#SBATCH --error=GATK.%j.err

cd $SLURM_SUBMIT_DIR

OUTDIR= "/scratch/ry00555/McEachern/"
FASTQ="/scratch/ry00555/McEachern/FastQ"
GENOME="/scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic"
#binsize for windows
BIN="25"

#smoothlength setting, for smoothing the enrichment curve
SMOOTH="50"
#number of CPUs
THREADS=12

# mkdir TrimmedReads
# mkdir BigWigs
# mkdir SortedBamFiles


#Load these modules that are compatible with GATK version 4.3
#module load BEDOPS/2.4.39-foss-2019b

#ml R/3.6.2-foss-2019b
# curl -s http://ftp.ensemblgenomes.org/pub/fungi/release-58/fasta/fungi_ascomycota1_collection/kluyveromyces_lactis_gca_000002515/dna/Kluyveromyces_lactis_gca_000002515.ASM251v1.dna.toplevel.fa.gz | gunzip -c > klactis.fasta
#
# curl -s http://ftp.ensemblgenomes.org/pub/fungi/release-58/fasta/fungi_ascomycota1_collection/kluyveromyces_lactis_gca_000002515/dna/Kluyveromyces_lactis_gca_000002515.ASM251v1.dna.toplevel.fa.gz | gunzip -c > $WorkDir/Genome/klactis.fasta
#
# module load SAMtools
# samtools faidx $WorkDir/Genome/klactis.fasta
# samtools faidx klactis.fasta
#
# module load Bowtie2
#
# faidx --transform bed klactis.fasta > klactis.bed
#
# ml picard/2.27.4-Java-13.0.2
#
# java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary \
#       -R GCF_000002515.2_ASM251v1_genomic.fna \
#       -O GCF_000002515.2_ASM251v1_genomic.dict
#
#       java -jar $EBROOTPICARD/picard.jar BedToIntervalList \
#       -I /scratch/ry00555/McEachern/Genome/Klactisgenes.bed \
#       -R /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic.fna \
#       -SD /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic.dict \
#       -O /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic.interval_list
#
# ml GATK/4.3.0.0-GCCcore-8.3.0-Java-1.8
#
# #WGS uses 1000 bp bins
#       gatk PreprocessIntervals \
#       -R /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic.fna \
#       -L /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic.interval_list \
#       --interval-merging-rule OVERLAPPING_ONLY \
#       --bin-length 1000 \
#       --padding 0 \
#       -O /scratch/ry00555/McEachern/Genome/klactis_preprocessed1000_intervals.interval_list
#
#       gatk AnnotateIntervals \
#       -R GCF_000002515.2_ASM251v1_genomic.fna \
#       -L klactis_preprocessed1000_intervals.interval_list \
#        --interval-merging-rule OVERLAPPING_ONLY \
#        -O GCF_000002515.2_ASM251v1_genomic_preprocessed10_annotated_intervals.tsv

    #   ml Trim_Galore
    #   trim_galore --paired --length 20 --fastqc --gzip -o /scratch/ry00555/McEachern/TrimmedReads /scratch/ry00555/McEachern/FastQ/*fastq\.gz
       # #
       FILES="/scratch/ry00555/McEachern/TrimmedReads/*R1_001_val_1\.fq\.gz" #Don't forget the *
       # #
       #
       # #
       # #Iterate over the files
       for f in $FILES
       do
         file=${f##*/}
       # 	#remove ending from file name to create shorter names for bam files and other downstream output
        	name=${file/%_S[1-12]*_R1_001_val_1.fq.gz/}
       #
       # #
       # # 	# File Vars
       # # 	#use sed to get the name of the second read matching the input file
        	read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R2_001_val_2\.fq\.gz/g')
       # 	#variable for naming bam file
       bam="/scratch/ry00555/McEachern/SortedBamFiles/${name}.bam"
       # 	#variable name for bigwig output
        	bigwig="/scratch/ry00555/McEachern/BigWigs/${name}"
       # 	#QualityBam="${OUTDIR}/SortedBamFiles/${name}_Q30.bam"
       # #
       #
       ml SAMtools
       # #
bwa mem -M -v 3 -t 12 "/scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic.fna" $f $read2 | samtools view -bhSu - | samtools sort -@ 12 -T /scratch/ry00555/McEachern/SortedBamFiles/tempReps -o "$bam" -
samtools index "$bam"
        ml deepTools
       # #Plot all reads
bamCoverage -p 12 -bs 25 --normalizeUsing BPM --smoothLength 50 -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"
       #
        #plot mononucleosomes
        #bamCoverage -p $THREADS --MNase -bs 1 --normalizeUsing BPM --smoothLength 25 -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}_MNase.bw"
        done
