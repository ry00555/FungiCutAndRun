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

OUTDIR="/scratch/ry00555/ParpMus30_Run139"
SORTED_BAM_DIR="/scratch/ry00555/ParpMus30_Run139/OutputBams"
BAMDIR="/scratch/ry00555/OutputRun139/SortedBamFiles"
genome=/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic.fna

 curl -s
"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/925/GCF_000182925.2_NC12/GCF_000182925.2_NC12_genomic.fna.gz" | gunzip -c > $OUTDIR/Genome/Ncrassa.fasta

module load SAMtools
samtools faidx $OUTDIR/Genome/Ncrassa.fasta

module load Bowtie2
faidx --transform bed $OUTDIR/Genome/Ncrassa.fasta > $OUTDIR/Genome/Ncrassa.bed

curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/925/GCF_000182925.2_NC12/GCF_000182925.2_NC12_genomic.gtf.gz | gunzip -c > $OUTDIR/Genome/Ncrassa.gtf
module load BEDOPS/2.4.39-foss-2019b
 gtf2bed < $OUTDIR/Genome/Ncrassa.gtf > $OUTDIR/Genome/Ncrassa2.bed --max-mem 10G
awk '$8 == "gene"' $OUTDIR/Genome/Ncrassa2.bed > $OUTDIR/Genome/Ncrassa2_genes.bed


ml picard
#
 java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary \
       -R $genome \
       -O $OUTDIR/Genome/Ncrassa.dict

  java -jar $EBROOTPICARD/picard.jar BedToIntervalList \
 -I $OUTDIR/Genome/Ncrassa.bed \
-R $genome \
-SD $OUTDIR/Genome/Ncrassa.dict \
 -O $OUTDIR/Genome/Ncrassa.interval_list


 ml GATK
 #
 # #WGS uses 1000 bp bins
        gatk PreprocessIntervals \
       -R $genome \
              -L $OUTDIR/Genome/Ncrassa.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
      --bin-length 1000 \
      --padding 0 \
        -O $OUTDIR/Genome/Ncrassa_1000_intervals.interval_list

 gatk AnnotateIntervals \
 -R $genome  \
-L $OUTDIR/Genome/Ncrassa_1000_intervals.interval_list \
--interval-merging-rule OVERLAPPING_ONLY \
 -O $OUTDIR/Genome/Ncrassa_preprocessed10_annotated_intervals.tsv

 ml picard
 for bam_file in $BAMDIR/139*WGS*.bam
 do
 # # #     # Get the base name of the BAM file
   base_name=$(basename "$bam_file" .bam)
 # # #
 # # #     # Define the output file path
  output_file="${SORTED_BAM_DIR}/${base_name}_output.bam"
 # # #
 # # #     # Run Picard to add or replace read groups
 java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
 -I "$bam_file" \
  -O "$output_file" \
  -RGID 2 \
 -RGLB lib1 \
-RGPL illumina \
-RGPU S34 \
-RGSM "${base_name%.*}"
  done
