#!/bin/bash
#SBATCH --job-name=Trichoderma
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=../Trichoderma.%j.out
#SBATCH --error=../Trichoderma.%j.err

cd $SLURM_SUBMIT_DIR

#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )

source configforDM.txt

OUTDIR="/scratch/ry00555/Run144Trichoderma"

if [ ! -d $OUTDIR ]
then
mkdir -p $OUTDIR
mkdir -p "${OUTDIR}/TrimmedReads"
mkdir -p "${OUTDIR}/BigWigs"
mkdir -p "$OUTDIR/TdfFiles"
   mkdir -p "$OUTDIR/HomerTagDirectories"
 mkdir -p "$OUTDIR/SortedBamFiles"
 mkdir -p "$OUTDIR/MACSPeaks"
fi


#
#
TAGDIR="${OUTDIR}/HomerTagDirectories"
BAMDIR="${OUTDIR}/SortedBamFiles"
BEDDIR="${OUTDIR}/Beds"
#
# # #process reads using trimGalore
# #
 ml Trim_Galore/0.6.7-GCCcore-11.2.0
 trim_galore --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${FASTQ}/*fastq\.gz
# #
 FILES="${OUTDIR}/TrimmedReads/*_L002_R1_001\.fq\.gz" #Don't forget the *
# #
#
# #
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
# 	#remove ending from file name to create shorter names for bam files and other downstream output
name=${file/%_S[1-150]*_L002_R1_001.fq.gz/}
# #
# # 	# File Vars
# # 	#use sed to get the name of the second read matching the input file
read2=$(echo "$f" | sed 's/R1_001\.fq\.gz/R2_001\.fq\.gz/g')
 	#variable for naming bam file
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
 #plot mononucleosomes
bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --minMappingQuality 10 --smoothLength $SMOOTH -of bigwig -b "$QualityBam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}_Q30.bw"
done
mkdir -p $OUTDIR/MACSPeaks
PEAKDIR="${OUTDIR}/MACSPeaks"

module load MACS3/3.0.0b1-foss-2022a-Python-3.10.4
 #command line
 #macs3 callpeak -t 137-11_CUTANDRUN_rtt109_H3K36me3_Rep1_S11_Ecoli.sorted.bam -f BAMPE -n 137-11_CUTANDRUN_rtt109_H3K36me3_Rep1_S11_Ecoli -c 137-9_CUTANDRUN_rtt109_IgG_Rep1_S9_Ecoli.sorted.bam --broad -g 41037538 --broad-cutoff 0.1 --outdir /scratch/ry00555/OutputRun137/CutandRun/MACSPeaks --min-length 800 --max-gap 500

for infile in $BAMDIR/*__Q30.bam
do
   base=$(basename ${infile} __Q30.bam)
      Input=$BAMDIR/${infile} Input__Q30.bam
macs3 callpeak -t $infile -f BAMPE -n $base -c $Input --broad -g 41037538 --broad-cutoff 0.1 --outdir $PEAKDIR --min-length 800 --max-gap 500
done

HOMERPEAKSDIR="${OUTDIR}/HomerPeaks"
  ml Homer
 ml Perl
 ml SAMtools
  ml BEDTools

    for bam_file in ${BAMDIR}/*_L007_R1_001_val_1.fq.gz_Q30.bam; do
  sample_id=$(basename "${bam_file}" _L007_R1_001_val_1.fq.gz_Q30.bam)
input_id=$(basename  *Input*_L007_R1_001_val_1.fq.gz_Q30.bam)

  # Define HOMERINPUT by matching the control input based on the sample_id
  #makeTagDirectory "${TAGDIR}/${sample_id}" "${bam_file}"

  HOMERINPUT="${TAGDIR}/${input_id}"
 findPeaks "${TAGDIR}/${sample_id}" -style factor -o "${HOMERPEAKSDIR}/${sample_id}_Homerpeaks.txt" -i $HOMERINPUT

  done
# #changing peak txt files to bed files to input into chipr
 ml ChIP-R
  for infile in ${HOMERPEAKSDIR}/*_Homerpeaks.txt
 do
   sample_id=$(basename "${infile}" _Homerpeaks.txt)

   sed '/^#/d' $infile | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $8 "\t" $5 "\t" $6 "\t" $12 "\t" "-1"}' | sed 's/\.000000//g' > ${HOMERPEAKSDIR}/${sample_id}.peaks.bed
#   sed '/^#/d' 142-142_ChIP_WT_GFP_Rep1_S142_Homerpeaks.txt | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $8 "\t" $5 "\t" $6 "\t" $12 "\t" "-1"}' | sed 's/\.000000//g' > 142-142_ChIP_WT_GFP_Rep1_S142_Homerpeaks.bed
  # sed '/^#/d' 142-144_ChIP_cre1-GFP_GFP_Rep1_S144_Homerpeaks.txt | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $8 "\t" $5 "\t" $6 "\t" $12 "\t" "-1"}' | sed 's/\.000000//g' > 142-144_ChIP_cre1-GFP_GFP_Rep1_S144_Homerpeaks.bed


 done
#
 ml Homer
 ml Perl
# ##annotating peak files with masked reference (use HOMER module)
# #curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/925/GCF_000182925.2_NC12/GCF_000182925.2_NC12_genomic.gtf.gz | gunzip -c > Ncrassa_refann.gtf
 annotatePeaks.pl ${HOMERPEAKSDIR}/${base}.peaks.bed /home/ry00555/Research/Genomes/TrichodermaReesiQM6a/TrichodermaReesiQM6a_GCA_000167675.2_v2.0_genomic.fna -gtf ry00555/Research/Genomes/TrichodermaReesiQM6a_GCA_000167675.2_v2.0_genomic.gtf > ${HOMERPEAKSDIR}/${base}_ann.txt
 #annotatePeaks.pl 142-144_ChIP_cre1-GFP_GFP_Rep1_S144_Homerpeaks.bed /home/ry00555/Research/Genomes/TrichodermaReesiQM6a/TrichodermaReesiQM6a_GCA_000167675.2_v2.0_genomic.fna -gtf  /home/ry00555/Research/Genomes/TrichodermaReesiQM6a/TrichodermaReesiQM6a_GCA_000167675.2_v2.0_genomic.gtf > 142-144_ChIP_cre1-GFP_GFP_Rep1_S144_Homerpeaks_ann.txt
# annotatePeaks.pl 142-142_ChIP_WT_GFP_Rep1_S142_Homerpeaks.bed /home/ry00555/Research/Genomes/TrichodermaReesiQM6a/TrichodermaReesiQM6a_GCA_000167675.2_v2.0_genomic.fna -gtf /home/ry00555/Research/Genomes/TrichodermaReesiQM6a/TrichodermaReesiQM6a_GCA_000167675.2_v2.0_genomic.gtf > 142-142_ChIP_WT_GFP_Rep1_S142_Homerpeaks_ann.txt



# #now filtering for only peaks that are w/i 1000bps of their annotation:
  for infile in ${HOMERPEAKSDIR}/${base}_ann.txt
  do
    base=$(basename ${infile} _masked_ann.txt)
    awk -F'\t' 'sqrt($10*$10) <=1000' $infile > ${HOMERPEAKSDIR}/${base}.1000bp_ann.txt
  done
