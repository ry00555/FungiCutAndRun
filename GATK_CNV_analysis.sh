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

 curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/925/GCF_000182925.2_NC12/GCF_000182925.2_NC12_genomic.fna.gz | gunzip -c > $OUTDIR/Genome/Ncrassa.fasta

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

#  ml picard
#  for bam_file in $BAMDIR/139*WGS*.bam
#  do
#  # # #     # Get the base name of the BAM file
#    base_name=$(basename "$bam_file" .bam)
#  # # #
#  # # #     # Define the output file path
#   output_file="${SORTED_BAM_DIR}/${base_name}_output.bam"
#  # # #
#  # # #     # Run Picard to add or replace read groups
#  java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
#  -I "$bam_file" \
#   -O "$output_file" \
#   -RGID 2 \
#  -RGLB lib1 \
# -RGPL illumina \
# -RGPU S34 \
# -RGSM "${base_name%.*}"
#   done


  ml GATK
for bam_file in $SORTED_BAM_DIR/139*WGS*.bam
 do
# Get the base name of the BAM file
      base_name=$(basename "$bam_file" _output.bam)
  # # # # #   # Define the output file path
  ml SAMtools
samtools index "$SORTED_BAM_DIR/139*WGS*_output.bam"
  # # # # #
gatk CollectReadCounts \
 -I "$bam_file" \
 -R $genome  \
-L /scratch/ry00555/McEachern/Genome/klactis_preprocessed1000_intervals.interval_list \
 --interval-merging-rule OVERLAPPING_ONLY \
 -O /scratch/ry00555/McEachern/CountTSVs/$base_name.counts.tsv
  done
CountTSVsDIR="$OUTDIR/CountTSVs"


 #gatk CreateReadCountPanelOfNormals \
# # # -I ${CountTSVsDIR}/113-1-gDNA-CBS2359_merged.counts.tsv \
# # # -I ${CountTSVsDIR}/113-12-gDNA-7B520_merged.counts.tsv  \
# # # --annotated-intervals /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic_preprocessed10_annotated_intervals.tsv \
# # # -O ${OUTDIR}/PanelofNormals/K_Samples.pon.hdf5
# # #
#  for count_files in $CountTSVsDIR/113*.counts.tsv
#   do
# # # # Get the base name of the counts file
# base_name=$(basename "$count_files" .counts.tsv)
# # #    # Define the output file path
# gatk DenoiseReadCounts \
#  -I "$count_files" \
# --annotated-intervals /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic_preprocessed10_annotated_intervals.tsv \
# --count-panel-of-normals ${OUTDIR}/PanelofNormals/K_Samples.pon.hdf5 \
# --standardized-copy-ratios ${OUTDIR}/CopyRatios/${base_name}.standardizedCR.tsv \
# --denoised-copy-ratios ${OUTDIR}/CopyRatios/${base_name}.denoisedCR.tsv
# done
