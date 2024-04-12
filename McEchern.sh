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


 FILES="${OUTDIR}/TrimmedReads/*R1_001_val_1.fq.gz" # Don't forget the *
#
# # Iterate over the files
 for f in $FILES
 do
     file=${f##*/}
     name=${file/%_S[1-12]*R1_001_val_1.fq.gz/}
     read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R1_001_val_2\.fq\.gz/g')
     bam="/scratch/ry00555/McEachern/SortedBamFiles/${name}.bam"
bigwig="/scratch/ry00555/McEachern/BigWigs/${name}"
ml SAMtools
ml BWA
 bwa mem -M -v 3 -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T /scratch/ry00555/McEachern/SortedBamFiles/tempReps -o "$bam" -
samtools index "$bam"
ml deepTools
bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"
 done

 FILES113="${OUTDIR}/TrimmedReads/113*_merged_trimmed.fq.gz" # Don't forget the *

 for f in $FILES113
 do
     file=${f##*/}
     name=${file/*_merged_trimmed.fq.gz/}
     #read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R1_001_val_2\.fq\.gz/g')
     bam="/scratch/ry00555/McEachern/SortedBamFiles/${name}.bam"
bigwig="/scratch/ry00555/McEachern/BigWigs/${name}"
ml SAMtools
ml BWA
 bwa mem -M -v 3 -t $THREADS $GENOME $f | samtools view -bhSu - | samtools sort -@ $THREADS -T /scratch/ry00555/McEachern/SortedBamFiles/tempReps -o "$bam" -
samtools index "$bam"
ml deepTools
bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"
 done
# Set the directory containing the sorted BAM files
SORTED_BAM_DIR="/scratch/ry00555/McEachern/SortedBamFiles/"

# Iterate over all BAM files in the directory
 ml picard
for bam_file in $SORTED_BAM_DIR/*.bam
do
    # Get the base name of the BAM file
    base_name=$(basename "$bam_file" .bam)

    # Define the output file path
    output_file="${SORTED_BAM_DIR}/${base_name}_output.bam"

    # Run Picard to add or replace read groups
    java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
    -I "$bam_file" \
    -O "$output_file" \
    -RGID 2 \
    -RGLB lib1 \
    -RGPL illumina \
    -RGPU S34 \
    -RGSM "${base_name%.*}"
done

#mkdir CountTSVs
ml GATK
 for bam_file in $SORTED_BAM_DIR/*_output.bam
 do
#   # Get the base name of the BAM file
   base_name=$(basename "$bam_file" _output.bam)
#   # Define the output file path
   input_file="${SORTED_BAM_DIR}/${base_name}"
samtools index "$input_file"
#
 gatk CollectReadCounts \
 -I "$input_file" \
 -R /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic.fna \
 -L /scratch/ry00555/McEachern/Genome/klactis_preprocessed1000_intervals.interval_list \
 --interval-merging-rule OVERLAPPING_ONLY \
 -O /scratch/ry00555/McEachern/CountTSVs/$base_name.counts.tsv

done

CountTSVsDIR="/scratch/ry00555/McEachern/CountTSVs/"

# gatk CreateReadCountPanelOfNormals \
# -I ${CountTSVsDIR}/138-1_Genomic_K1__Rep1_6252.bam_output.bam.counts.tsv \
# -I ${CountTSVsDIR}/138-2_Genomic_K2__Rep1_6252.bam_output.bam.counts.tsv  \
# -I ${CountTSVsDIR}/138-3_Genomic_K3__Rep1_6252_S3_L001_R1_001_val_1.fq.gz.bam_output.bam.counts.tsv \
# --annotated-intervals /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic_preprocessed10_annotated_intervals.tsv \
# -O ${OUTDIR}/PanelofNormals/K_Samples.pon.hdf5
#
# for count_files in $CountTSVsDIR/*M*tsv
# do
#
#   # Get the base name of the counts file
#      base_name=$(basename "$count_files" .counts.tsv)
#  #   # Define the output file path
#   input_file="${CountTSVsDIR}/${base_name}"
# gatk DenoiseReadCounts \
# -I "$input_file" \
# --annotated-intervals /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic_preprocessed10_annotated_intervals.tsv \
# --count-panel-of-normals ${OUTDIR}/PanelofNormals/K_Samples.pon.hdf5 \
# --standardized-copy-ratios ${OUTDIR}/CopyRatios/${base_name}.standardizedCR.tsv \
# --denoised-copy-ratios ${OUTDIR}/CopyRatios/${base_name}.denoisedCR.tsv
#
# done
#
# for copy_ratios in ${OUTDIR}/CopyRatios/
# do
# # Get the base name of the counts file
#    base_name=$(basename "$copy_ratios")
# #   # Define the output file path
# input_file="${OUTDIR}/CopyRatios/${base_name}"
#
# gatk PlotDenoisedCopyRatios \
# --standardized-copy-ratios ${OUTDIR}/CopyRatios/${base_name}.standardizedCR.tsv \
# --denoised-copy-ratios ${OUTDIR}/CopyRatios/${base_name}.denoisedCR.tsv  \
# --sequence-dictionary /scratch/ry00555/McEachern/Genome/GCF_000002515.2_ASM251v1_genomic.dict \
# --point-size-copy-ratio 1 \
# --output-prefix ${base_name} \
# --output ${OUTDIR}/PlotDenoisedCopyRatios
# done
