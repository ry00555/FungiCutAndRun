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
genome="/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic.fna"

#  curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/925/GCF_000182925.2_NC12/GCF_000182925.2_NC12_genomic.fna.gz | gunzip -c > $OUTDIR/Genome/Ncrassa.fasta
#
# module load SAMtools
# samtools faidx $OUTDIR/Genome/Ncrassa.fasta
# awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' $OUTDIR/Genome/Ncrassa.fasta.fai > $OUTDIR/GenomeNcrassa.bed
# module load Bowtie2
# faidx --transform bed $OUTDIR/Genome/Ncrassa.fasta > $OUTDIR/Genome/Ncrassa.bed
#wrong chromosomes
#curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/925/GCF_000182925.2_NC12/GCF_000182925.2_NC12_genomic.gtf.gz | gunzip -c > $OUTDIR/Genome/Ncrassa.gtf
# module load BEDOPS
#  gtf2bed < $OUTDIR/Genome/Ncrassa.gtf > $OUTDIR/Genome/Ncrassa2.bed --max-mem 10G
# #
#  awk '$8 == "gene"' $OUTDIR/Genome/Ncrassa2.bed > $OUTDIR/Genome/Ncrassa2_genes.bed
# #
#
# ml picard
#
 # java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary \
 #       -R $genome \
 #       -O $OUTDIR/Genome/Ncrassa.dict

# java -jar $EBROOTPICARD/picard.jar BedToIntervalList \
#   -I $OUTDIR/Genome/Ncrassa2_genes.bed \
#  -R $genome \
#  -SD $OUTDIR/Genome/Ncrassa.dict \
#   -O $OUTDIR/Genome/Ncrassa.interval_list
# #
# #
#   ml GATK
# #  #
# #  # #WGS uses 1000 bp bins
#          gatk PreprocessIntervals \
#         -R $genome \
# -L $OUTDIR/Genome/Ncrassa.interval_list \
#      --interval-merging-rule OVERLAPPING_ONLY \
#        --bin-length 1000 \
#        --padding 0 \
#          -O $OUTDIR/Genome/Ncrassa_1000_intervals.interval_list
# #
#  gatk AnnotateIntervals \
#   -R $genome  \
#  -L $OUTDIR/Genome/Ncrassa_1000_intervals.interval_list \
#  --interval-merging-rule OVERLAPPING_ONLY \
#   -O $OUTDIR/Genome/Ncrassa_preprocessed10_annotated_intervals.tsv

  #ml picard
#   for bam_file in $BAMDIR/139*WGS*.bam
#   do
# #  # # #     # Get the base name of the BAM file
#     base_name=$(basename "$bam_file" .bam)
# #  # # #
# #  # # #     # Define the output file path
#    output_file="${SORTED_BAM_DIR}/${base_name}_output.bam"
# #  # # #
# #  # # #     # Run Picard to add or replace read groups
#   java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
#   -I "$bam_file" \
#    -O "$output_file" \
#       -RGID 2 \
#  -RGLB lib1 \
# -RGPL illumina \
# -RGPU S34 \
# -RGSM "${base_name%.*}"
#   done


  ml GATK
  OUTPUTBAM="$SORTED_BAM_DIR/*Q30_output.bam"

# for bam_file in $OUTPUTBAM
#  do
# # Get the base name of the BAM file
#       base_name=$(basename "$bam_file" Q30_output.bam)
#   # # # # #   # Define the output file path
#  ml SAMtools
# #samtools index "$SORTED_BAM_DIR/*_output.bam"
#   # # # # #
# gatk CollectReadCounts \
#  -I "$bam_file" \
#  -R $genome  \
# -L $OUTDIR/Genome/Ncrassa_1000_intervals.interval_list \
#  --interval-merging-rule OVERLAPPING_ONLY \
#  -O" $OUTDIR/CountTSVs/$base_name.counts.tsv"
#   done


 gatk CreateReadCountPanelOfNormals \
-I $OUTDIR/CountTSVs/139-20_WGS_WT___.counts.tsv \
--annotated-intervals $OUTDIR/Genome/Ncrassa_preprocessed10_annotated_intervals.tsv \
-O ${OUTDIR}/PanelofNormals/139-20_WGS_WT___.pon.hdf5
# # #
CountTSVsDIR="$OUTDIR/CountTSVs"
for count_files in $CountTSVsDIR/*.counts.tsv
   do
# # # # Get the base name of the counts file
 base_name=$(basename "$count_files" .counts.tsv)
# # #    # Define the output file path
 gatk DenoiseReadCounts \
  -I "$count_files" \
  --annotated-intervals $OUTDIR/Genome/Ncrassa_preprocessed10_annotated_intervals.tsv \
 --count-panel-of-normals ${OUTDIR}/PanelofNormals/139-20_WGS_WT___.pon.hdf5 \
 --standardized-copy-ratios ${OUTDIR}/CopyRatios/${base_name}.standardizedCR.tsv \
 --denoised-copy-ratios ${OUTDIR}/CopyRatios/${base_name}.denoisedCR.tsv
 done

 for copy_ratios in ${OUTDIR}/CopyRatios/*.standardizedCR.tsv
 do
 # # #Get the base name of the counts file
   base_name=$(basename "$copy_ratios" .standardizedCR.tsv)
   Denoised=$(echo "$copy_ratios" | sed 's/\.standardizedCR\.tsv/\.denoisedCR\.tsv/g')

 # # Define the output file path
   gatk PlotDenoisedCopyRatios \
  --standardized-copy-ratios "$copy_ratios" \
  --denoised-copy-ratios "${OUTDIR}/CopyRatios/$Denoised" \
 --sequence-dictionary $OUTDIR/Genome/Ncrassa.dict \
  --point-size-copy-ratio 1 \
  --output-prefix ${base_name} \
  --output ${OUTDIR}/PlotDenoisedCopyRatios
  done

  for copy_ratios in ${OUTDIR}/CopyRatios/*.denoisedCR.tsv
    do
  # # #
  base_name=$(basename "$copy_ratios" .denoisedCR.tsv)
  # # #
 gatk ModelSegments \
    --denoised-copy-ratios ${OUTDIR}/CopyRatios/${base_name}.denoisedCR.tsv \
  --output-prefix ${base_name} \
    -O ${OUTDIR}/ModelSegments
  #
   gatk PlotModeledSegments \
  --denoised-copy-ratios ${OUTDIR}/CopyRatios/${base_name}.denoisedCR.tsv \
 --segments ${OUTDIR}/ModelSegments/${base_name}.modelFinal.seg \
 --sequence-dictionary $OUTDIR/Genome/Ncrassa.dict \
    --point-size-copy-ratio 1 \
      --output-prefix ${base_name} \
     -O ${OUTDIR}/PlotModelSegments
     done
