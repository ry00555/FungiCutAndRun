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

#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )

source config.txt

#Make output directory
OUTDIR="/scratch/ry00555/McEachern/"


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
 #ml picard
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
#ml GATK
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

#ml Trim_Galore
#trim_galore --paired --length 20 --fastqc --gzip -o /scratch/ry00555/McEachern/TrimmedReads /scratch/ry00555/McEachern/FastQ/113*fastq\.gz
       # #
     #
      # mkdir "${OUTDIR}/SortedBamFiles"
      # mkdir "${OUTDIR}/BigWigs"
      # mkdir "${OUTDIR}/Peaks"
     #mkdir "$OUTDIR/HomerTagDirectories"
     #mkdir "$OUTDIR/TdfFiles"


 # FILES="${OUTDIR}/TrimmedReads/*R1_001_val_1.fq.gz" # Don't forget the *
#
# # Iterate over the files
#  for f in $FILES
#  do
#      file=${f##*/}
#      name=${file/%_S[1-12]*R1_001_val_1.fq.gz/}
#      read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R1_001_val_2\.fq\.gz/g')
#      bam="/scratch/ry00555/McEachern/SortedBamFiles/${name}.bam"
# bigwig="/scratch/ry00555/McEachern/BigWigs/${name}"
# ml SAMtools
# ml BWA
#  bwa mem -M -v 3 -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T /scratch/ry00555/McEachern/SortedBamFiles/tempReps -o "$bam" -
# samtools index "$bam"
# ml deepTools
# bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"
#  done

# bam="/scratch/ry00555/McEachern/SortedBamFiles/113-1-gDNA-CBS2359_merged.bam"
# bigwig="/scratch/ry00555/McEachern/BigWigs/113-1-gDNA-CBS2359_merged"
# bam2="/scratch/ry00555/McEachern/SortedBamFiles/113-12-gDNA-7B520_merged.bam"
# bigwig2="/scratch/ry00555/McEachern/BigWigs/113-12-gDNA-7B520_merged"
#
#ml SAMtools
# ml BWA
#  bwa mem -M -v 3 -t $THREADS $GENOME ${OUTDIR}/TrimmedReads/113-1-gDNA-CBS2359_merged_trimmed.fq.gz | samtools view -bhSu - | samtools sort -@ $THREADS -T /scratch/ry00555/McEachern/SortedBamFiles/tempReps -o "$bam" -
#  bwa mem -M -v 3 -t $THREADS $GENOME ${OUTDIR}/TrimmedReads/113-12-gDNA-7B520_merged_trimmed.fq.gz | samtools view -bhSu - | samtools sort -@ $THREADS -T /scratch/ry00555/McEachern/SortedBamFiles/tempReps -o "$bam2" -
#
# samtools index "$bam"
# samtools index "$bam2"
#
# ml deepTools
# bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"
# bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --smoothLength $SMOOTH -of bigwig -b "$bam2" -o "${bigwig2}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"

# Set the directory containing the sorted BAM files
SORTED_BAM_DIR="/scratch/ry00555/McEachern/SortedBamFiles"

# Iterate over all BAM files in the directory
#  ml picard
# for bam_file in $SORTED_BAM_DIR/*.bam
# do
#     # Get the base name of the BAM file
#     base_name=$(basename "$bam_file" .bam)
#
#     # Define the output file path
#     output_file="${SORTED_BAM_DIR}/${base_name}_output.bam"
#
#     # Run Picard to add or replace read groups
#     java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
#     -I "$bam_file" \
#     -O "$output_file" \
#     -RGID 2 \
#     -RGLB lib1 \
#     -RGPL illumina \
#     -RGPU S34 \
#     -RGSM "${base_name%.*}"
# done
CountTSVsDIR="/scratch/ry00555/McEachern/CountTSVs"

#mkdir CountTSVs
#did this in command line
#   input_file="${SORTED_BAM_DIR}/*_output.bam"

#   ml SAMtools
#samtools index "$input_file"

#ml GATK
# for bam_file in ${SORTED_BAM_DIR}/*_output.bam; do
#     # Get the base name of the BAM file
#     base_name=$(basename "$bam_file" _output.bam)
#     gatk CollectReadCounts \
#     -I "$bam_file" \
#     -R /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic.fna \
#     -L /scratch/ry00555/McEachern/Genome/klactis_preprocessed1000_intervals.interval_list \
#     --interval-merging-rule OVERLAPPING_ONLY \
#     -O "${CountTSVsDIR}/${base_name}.counts.tsv"
# done




# gatk CreateReadCountPanelOfNormals \
# -I ${CountTSVsDIR}/113-1-gDNA-CBS2359_merged.counts.tsv \
# -I ${CountTSVsDIR}/113-12-gDNA-7B520_merged.counts.tsv  \
# --annotated-intervals /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic_preprocessed10_annotated_intervals.tsv \
# -O ${OUTDIR}/PanelofNormals/113_WT_Samples.pon.hdf5
#
# for count_files in $CountTSVsDIR/*.counts.tsv
# do
#
# #   # Get the base name of the counts file
# base_name=$(basename "$count_files" .counts.tsv)
# #  #   # Define the output file path
# gatk DenoiseReadCounts \
# -I "$count_files" \
# --annotated-intervals /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic_preprocessed10_annotated_intervals.tsv \
# --count-panel-of-normals ${OUTDIR}/PanelofNormals/113_WT_Samples.pon.hdf5 \
# --standardized-copy-ratios ${OUTDIR}/CopyRatios/${base_name}.standardizedCR.tsv \
# -denoised-copy-ratios ${OUTDIR}/CopyRatios/${base_name}.denoisedCR.tsv
#
# done
#for some reason even with the new WT panel of normals they don't work just proceed with K samples error is
#13:25:48.915 INFO  H5 - HDF5 library:
# 13:25:48.915 INFO  H5 -  successfully loaded.
# 13:25:48.919 INFO  DenoiseReadCounts - Reading read-counts file (/scratch/ry00555/McEachern/CountTSVs/138-9_Genomic_M6__Rep1_6252_S9_L001_R1_001_val_1.fq.gz.counts.tsv)...
# 13:25:49.068 WARN  DenoiseReadCounts - Panel of normals was provided; ignoring input GC-content annotations...
# 13:25:49.084 INFO  SVDDenoisingUtils - Validating sample intervals against original intervals used to build panel of normals...
# 13:25:49.092 INFO  SVDDenoisingUtils - Preprocessing and standardizing sample read counts...
# 13:25:49.096 INFO  SVDDenoisingUtils - Preprocessing read counts...
# 13:25:49.097 INFO  SVDDenoisingUtils - Transforming read counts to fractional coverage...
# 13:25:49.099 INFO  SVDDenoisingUtils - Performing GC-bias correction...
# 13:25:49.164 INFO  SVDDenoisingUtils - Subsetting sample intervals to post-filter panel intervals...
# 13:25:49.193 INFO  SVDDenoisingUtils - Dividing by interval medians from the panel of normals...
# 13:25:49.194 INFO  SVDDenoisingUtils - Sample read counts preprocessed.
# 13:25:49.194 INFO  SVDDenoisingUtils - Standardizing read counts...
# 13:25:49.194 INFO  SVDDenoisingUtils - Dividing by sample medians and transforming to log2 space...
# 13:25:49.199 INFO  DenoiseReadCounts - Shutting down engine
# [April 12, 2024 at 1:25:49 PM EDT] org.broadinstitute.hellbender.tools.copynumber.DenoiseReadCounts done. Elapsed time: 0.01 minutes.
# Runtime.totalMemory()=162267136
# java.lang.IllegalArgumentException: Sample does not have a positive sample median.
#         at org.broadinstitute.hellbender.utils.Utils.validateArg(Utils.java:798)
#         at org.broadinstitute.hellbender.utils.param.ParamUtils.isPositive(ParamUtils.java:165)
#         at org.broadinstitute.hellbender.tools.copynumber.denoising.SVDDenoisingUtils.lambda$divideBySampleMedianAndTransformToLog2$27(SVDDenoisingUtils.java:484)
#         at java.base/java.util.stream.Streams$RangeIntSpliterator.forEachRemaining(Streams.java:104)
#         at java.base/java.util.stream.IntPipeline$Head.forEach(IntPipeline.java:617)
#         at org.broadinstitute.hellbender.tools.copynumber.denoising.SVDDenoisingUtils.divideBySampleMedianAndTransformToLog2(SVDDenoisingUtils.java:483)
#         at org.broadinstitute.hellbender.tools.copynumber.denoising.SVDDenoisingUtils.preprocessAndStandardizeSample(SVDDenoisingUtils.java:406)
#         at org.broadinstitute.hellbender.tools.copynumber.denoising.SVDDenoisingUtils.denoise(SVDDenoisingUtils.java:123)
#         at org.broadinstitute.hellbender.tools.copynumber.denoising.SVDReadCountPanelOfNormals.denoise(SVDReadCountPanelOfNormals.java:88)
ml R/4.3.1-foss-2022a
ml GATK/4.3.0.0-GCCcore-11.3.0-Java-11
for copy_ratios in ${OUTDIR}/CopyRatios/*.standardizedCR.tsv
do

# # Get the base name of the counts file
   base_name=$(basename "$copy_ratios" .standardizedCR.tsv)
# #   # Define the output file path
#
 gatk PlotDenoisedCopyRatios \
 --standardized-copy-ratios ${OUTDIR}/CopyRatios/${base_name}.standardizedCR.tsv \
 --denoised-copy-ratios ${OUTDIR}/CopyRatios/${base_name}.denoisedCR.tsv  \
 --sequence-dictionary /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic.dict \
 --point-size-copy-ratio 1 \
 --output-prefix ${base_name} \
 --output ${OUTDIR}/PlotDenoisedCopyRatios
 done

 #R not working try command line
# files: 113-1-gDNA-CBS2359_merged.denoisedCR.tsv      138-1_Genomic_K1__Rep1_6252.denoisedCR.tsv      138-3_Genomic_K3__Rep1_6252_S3_L001_R1_001_val_1.fq.gz.denoisedCR.tsv
# 113-1-gDNA-CBS2359_merged.standardizedCR.tsv  138-1_Genomic_K1__Rep1_6252.standardizedCR.tsv  138-3_Genomic_K3__Rep1_6252_S3_L001_R1_001_val_1.fq.gz.standardizedCR.tsv
# 113-12-gDNA-7B520_merged.denoisedCR.tsv       138-2_Genomic_K2__Rep1_6252.denoisedCR.tsv
# 113-12-gDNA-7B520_merged.standardizedCR.tsv   138-2_Genomic_K2__Rep1_6252.standardizedCR.tsv

# gatk PlotDenoisedCopyRatios \
# --standardized-copy-ratios 138-1_Genomic_K1__Rep1_6252.standardizedCR.tsv \
# --denoised-copy-ratios 138-1_Genomic_K1__Rep1_6252.denoisedCR.tsv  \
# --sequence-dictionary /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic.dict \
# --point-size-copy-ratio 1 \
# --output-prefix 138-1_Genomic_K1__Rep1_6252 \
# --output /scratch/ry00555/McEachern/PlotDenoisedCopyRatios
