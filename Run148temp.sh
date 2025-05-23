#!/bin/bash
#SBATCH --job-name=Run148ChIP
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100gb
#SBATCH --time=72:00:00
#SBATCH --output=../MapCutAndRun148.%j.out
#SBATCH --error=../MapCutAndRun148.%j.err

cd $SLURM_SUBMIT_DIR

#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )

source config.txt
OUTDIR="/scratch/ry00555/Run148"

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
    mkdir -p "${OUTDIR}/TrimmedReads"
    mkdir -p "${OUTDIR}/BigWigs"
   mkdir -p "$OUTDIR/HomerTagDirectories"
#   mkdir -p "$OUTDIR/TdfFiles"
  mkdir -p "$OUTDIR/SortedBamFiles"

fi

TAGDIR="${OUTDIR}/HomerTagDirectories"
BAMDIR="${OUTDIR}/SortedBamFiles"
BEDDIR="${OUTDIR}/Beds"
#
FILES="${OUTDIR}/temp/*_L003_R1_001_val_1\.fq\.gz"
#FILES="${OUTDIR}/TrimmedReads/*R1_001_val_1\.fq\.gz"#



# #Iterate over the files
for f in $FILES
 do
# #
# # 	#Examples to Get Different parts of the file name
# # 		#See here for details: http://tldp.org/LDP/abs/html/refcards.html#AEN22664
# 		#${string//substring/replacement}
# # 		#dir=${f%/*}
#
file=${f##*/}
 	#remove ending from file name to create shorter names for bam files and other downstream output
#name=${file/%_S[1-150]*_L002_R1_001_val_1.fq.gz/}#
#name=${file/%_S[1-990]*_L002_R1_001_val_1.fq.gz/}
name=${file/%_S[1-190]*_L003_R1_001_val_1.fq.gz/}

#
# 	# File Vars
# 	#use sed to get the name of the second read matching the input file
read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R2_001_val_2\.fq\.gz/g')

# # 	# File Vars
# # 	#use sed to get the name of the second read matching the input file
#read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R2_001_val_2\.fq\.gz/g')
# 	#variable for naming bam file
bam="${OUTDIR}/SortedBamFiles/${name}.bam"
# 	#variable name for bigwig output
bigwig="${OUTDIR}/BigWigs/${name}"
QualityBam="${OUTDIR}/SortedBamFiles/${name}_Q30.bam"
# #
#
ml SAMtools/1.16.1-GCC-11.3.0
ml BWA/0.7.17-GCCcore-11.3.0
# #
 bwa mem -M -v 3 -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIR/SortedBamFiles/tempReps -o "$bam" -
 samtools index "$bam"
#
samtools view -b -q 30 $bam > "$QualityBam"
samtools index "$QualityBam"
#
# ############################
# # # #deeptools
#
ml deepTools/3.5.2-foss-2022a
 #Plot all reads
bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --minMappingQuality 10 --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"
#
bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --minMappingQuality 10 --smoothLength $SMOOTH -of bigwig -b "$QualityBam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}_Q30.bw"
done
mkdir $OUTDIR/MACSPeaks
PEAKDIR="${OUTDIR}/MACSPeaks"

module load MACS3/3.0.0b1-foss-2022a-Python-3.10.4
 #command line
#macs3 callpeak -t 137-11_CUTANDRUN_rtt109_H3K36me3_Rep1_S11_Ecoli.sorted.bam -f BAMPE -n 137-11_CUTANDRUN_rtt109_H3K36me3_Rep1_S11_Ecoli -c 137-9_CUTANDRUN_rtt109_IgG_Rep1_S9_Ecoli.sorted.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir /scratch/ry00555/OutputRun137/CutandRun/MACSPeaks --min-length 800 --max-gap 500

  for infile in $BAMDIR/*_Q30.bam
 do
    base=$(basename ${infile} _Q30.bam)
# #  Input=$BAMDIR/ ${infile} Input_Q30.bam
  macs3 callpeak -t $infile -f BAMPE -n $base --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500 #-c $Input
  done
#
  HOMERPEAKSDIR="${OUTDIR}/HomerPeaks"
   ml Homer
  ml Perl
  ml SAMtools
   ml BEDTools
     for bam_file in "${BAMDIR}"/*__Q30.bam; do
# # #   # Get the sample ID from the BAM file name
   sample_id=$(basename "${bam_file}" __Q30.bam)
# # #   # Remove everything after "Rep_1" in the sample ID
# #HOMERINPUT="${TAGDIR}/${sample_id}_Input*"
#
# # #
   makeTagDirectory "${TAGDIR}/${sample_id}" "${bam_file}"
# # # #
# # # #   # Call peaks
# # # #
  findPeaks "${TAGDIR}/${sample_id}" -style histone -region -size 150 -minDist 530 -o "${HOMERPEAKSDIR}/${sample_id}_Homerpeaks.txt" #-i $HOMERINPUT
# # # #
   done
