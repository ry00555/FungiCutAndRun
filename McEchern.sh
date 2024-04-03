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

WorkDir = /scratch/ry00555/McEachern/

#Load these modules that are compatible with GATK version 4.3
ml GATK/4.3.0.0-GCCcore-8.3.0-Java-1.8
ml picard/2.27.4-Java-13.0.2
ml BWA/0.7.17-GCC-8.3.0
module load BEDTools/2.29.2-GCC-8.3.0
module load SAMtools/1.10-GCC-8.3.0
module load BEDOPS/2.4.39-foss-2019b

ml R/3.6.2-foss-2019b
curl -s http://ftp.ensemblgenomes.org/pub/fungi/release-58/fasta/fungi_ascomycota1_collection/kluyveromyces_lactis_gca_000002515/dna/Kluyveromyces_lactis_gca_000002515.ASM251v1.dna.toplevel.fa.gz | gunzip -c > klactis.fasta

curl -s http://ftp.ensemblgenomes.org/pub/fungi/release-58/fasta/fungi_ascomycota1_collection/kluyveromyces_lactis_gca_000002515/dna/Kluyveromyces_lactis_gca_000002515.ASM251v1.dna.toplevel.fa.gz | gunzip -c > $WorkDir/Genome/klactis.fasta

module load SAMtools
samtools faidx $WorkDir/Genome/klactis.fasta
samtools faidx klactis.fasta

module load Bowtie2
# bowtie2-build -f $WorkDir/Genome/klactis.fasta $WorkDir/Genome/klactis
bowtie2 -x  reference -U file.fasta -S file.sam


1) Index the reference genome and map your reads or FASTA sequences to it (for example with bowtie2)

# index reference genome (should be precomputed)
bowtie2-build klactis.fasta klactis
# map reads
bowtie2 -p 4 -x klactis -U klactis.fasta -S klactis.sam

faidx --transform bed klactis.fasta > klactis.bed


# compress SAM to a BAM (binary) file
samtools view -bS file.sam > file.bam
2) Create a BED files from the BAM file using bedtools

# extract BED file
bedtools bamtobed -i file.bam > file.bed
Alternatively, for sliding windows you can generate these from a reference sequence provided that you know the length of each chromosome (perhaps there is a way to extract these directly from reference.fasta):

# length per chromosome
samtools view -H file.bam | grep "SQ" | cut -d":" -f2-3 | sed 's/LN://' > file.chr.txt

# create sliding windows
bedtools makewindows -g file.chr.txt -w 1000 > windows.bed




java -jar picard.jar CreateSequenceDictionary \
      R=$WorkDir/Genome/klactis.fasta \
      O=$WorkDir/Genome/klactis.dict
