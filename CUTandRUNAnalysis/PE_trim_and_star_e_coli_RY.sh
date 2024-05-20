#!/bin/bash
#######################
##Help
##########################################
Help()
{
  #Display Help
  echo "PE Trim and STAR is used to trim and align paired end reads"
  echo "This program uses STAR for alignment"
  echo "To use Bowtie2 instead, please use Align and Bowtie"
  echo
  echo "SYNTAX: $0 [options] read_file_1.fq.gz read_file_2.fq.gz"
  echo
  echo "OPTIONS:"
  echo "-o      output directory, required"
  echo "-g      path to reference genome file, optional, default is z11 release 98"
  echo "-n      basename for output files, optional"
  echo "-m      multimapping options, default is one
                  'none' - uniquely aligned reads only
                  'one' - allow multimappers to align to only their best alignment
                  'random' - allow multimappers to align to one randomly selected alignment
                  'all' - allow multimappers to align to all of their alignments
                regardless of choice, multimapping reads will be marked with the NH flag"
  echo
  echo "-h      print this help message"
  echo
}
#########################################
if [ "$#" == 0 ]; then
    echo "Error"
    Help
    exit
fi

while getopts ":ho:n:g:m:" option; do
  case "$option" in
    h) # display Help
      Help
      exit;;
    o) #setting output directory
      OUTDIR=$OPTARG;;
    :) echo "Option -$OPTARG requires an argument"
      echo
      Help
      exit
      ;;
    n) outname=$OPTARG;;
    g) genome=$OPTARG ;;
    m) mode=$OPTARG ;;
    \?) # invalid option
      echo "Invalid option -$OPTARG"
      echo
      echo "Help:"
      Help
      exit
      ;;
  esac
done

read_file1=${@:$OPTIND:1}
read_file2=${@:$OPTIND+1:1}

echo
echo "Working directory is $OUTDIR"

echo "Read files given:"
for infile in $read_file1 $read_file2
do
  if [ -f "$infile" ] ; then
  echo "$(basename $infile)"
  fi
done

if [ -z "$outname" ]; then
  outname=$(basename $read_file1)
  echo "output name is $outname"
else
  echo "output name is $outname"
fi

echo "Determining Multimapping Mode"
if [ $mode == 'none' ] ; then
  analysis_mode=NONE
  echo "  Multimapping Mode is $analysis_mode"
elif [ $mode == 'one' ] ; then
  analysis_mode=ONE
  echo "  Multimapping Mode is $analysis_mode"
elif [ $mode == 'random' ] ; then
  analysis_mode=RANDOM
  echo "  Multimapping Mode is $analysis_mode"
elif [ $mode == 'all' ] ; then
  analysis_mode=ALL
  echo "  Multimapping Mode is $analysis_mode"
else
  analysis_mode=ONE
  echo "  Multimapping Mode is $analysis_mode by default"
fi

echo
echo "All data input checks complete."
echo
echo "Starting Analysis"

if [ -d "$OUTDIR/TrimmedReads" ]
then
    echo "Directory $OUTDIR/TrimmedReads exists."
else
  mkdir $OUTDIR/TrimmedReads
fi

echo "Trimmed reads and fastQC files going into $OUTDIR/Ecoli_Aligned/TrimmedReads"
echo "... loading TrimGalore"
module load Trim_Galore
echo "...starting trimming"
#adapted between our old lab and Goll lab
trim_galore --paired --length 20 --fastqc --length 20 --gzip -j 24 --output_dir $OUTDIR/TrimmedReads --paired $read_file1 $read_file2

echo "... trimming complete"
echo "...loading MultiQC"
###need to edit this to only happen when the # of trimmed files = the number of input files?
module load MultiQC
echo "... running MultiQC"
multiqc -o $OUTDIR/TrimmedReads -f $OUTDIR/TrimmedReads
echo "... QC analysis complete"

echo
echo "... counting trimmed reads"
for infile in $OUTDIR/TrimmedReads/"$outname"*R1*val*.fq.gz
do
  base=$(basename ${infile} _val_1.fq.gz)
  echo $base >> $OUTDIR/TrimmedReads/trimmed_read_stats.txt
  echo $(zcat $infile |wc -l)/4|bc >> $OUTDIR/TrimmedReads/trimmed_read_stats_ecoli.txt
done
echo "... trimmed reads counted and reported in $OUTDIR/TrimmedReads/trimmed_read_stats"
echo

if [ -f "$genome" ]; then
  echo "genome file provided by user"
  echo "genome file is $genome"
elif [ -f "$OUTDIR/ref/ecoli_refseq.fa" ]; then
  echo "genome already downloaded"
  genome=$OUTDIR/ref/ecoli_refseq.fa
  echo "genome file is $genome"
else
  echo "no genome file supplied"
  echo "downloading ecoli genome"
  curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz | gunzip -c > $OUTDIR/ref/ecoli_refseq.fa
  genome=$OUTDIR/ref/ecoli_refseq.fa
  echo "genome file is $genome"
fi

if [ -d "$OUTDIR/SortedBamFiles" ]
then
    echo "Directory $OUTDIR/SortedBamFiles exists."
else
  mkdir $OUTDIR/SortedBamFiles
fi

echo "Alignment files going into $OUTDIR/SortedBamFiles"
echo "... loading STAR"
module load STAR

if [ -f "$OUTDIR/ref/Genome" ]
then
    echo "Reference genome index exists."
else
  echo "... building reference genome index"
  mkdir $OUTDIR/ref/ecoli_genome
  STAR --runThreadN 20 --genomeSAindexNbases 8 --runMode genomeGenerate \
  --genomeDir $OUTDIR/ref/ecoli_genome --genomeFastaFiles $OUTDIR/ref/ecoli_refseq.fa
  echo "... aligning reads to "$genome" genome"
fi


for file in $OUTDIR/TrimmedReads/"$outname"*_val_*.fq.gz;
do
  if [[ $prefix ]]; then
    if [ $analysis_mode == 'NONE' ] ; then
      STAR --runThreadN 24 --genomeDir $OUTDIR/ref/ecoli_genome --outFileNamePrefix $OUTDIR/SortedBamFiles/"$outname".Ecoli \
      --readFilesCommand zcat --readFilesIn "$first" "$file" --outSAMtype BAM SortedByCoordinate \
      --outFilterMultimapNmax 1 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000
    elif [ $analysis_mode == 'ONE' ] ; then

      STAR --runThreadN 24 --genomeDir $OUTDIR/ref/ecoli_genome --outFileNamePrefix $OUTDIR/SortedBamFiles/"$outname".Ecoli \
      --readFilesCommand zcat --readFilesIn "$first" "$file" --outSAMtype BAM SortedByCoordinate \
      --outSAMmultNmax 1 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000
    elif [ $analysis_mode == 'RANDOM'] ; then
    STAR --runThreadN 24 --genomeDir $OUTDIR/ref/ecoli_genome --outFileNamePrefix $OUTDIR/SortedBamFiles/"$outname".Ecoli \
    --readFilesCommand zcat --readFilesIn "$first" "$file" --outSAMtype BAM SortedByCoordinate \
    --outMultimapperOrder Random --outSAMmultNmax 1 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000
    elif [ $analysis_mode == 'ALL' ] ; then
      STAR --runThreadN 24 --genomeDir $OUTDIR/ref/ecoli_genome --outFileNamePrefix $OUTDIR/SortedBamFiles/"$outname".Ecoli \
      --readFilesCommand zcat --readFilesIn "$first" "$file" --outSAMtype BAM SortedByCoordinate \
      --outSAMprimaryFlag AllBestScore --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000
    fi
    echo "$outname is finished aligning"
    prefix=
  else
    first=$file
    prefix=${file%%_*}
  fi
done

rm $OUTDIR/SortedBamFiles/"$outname"*SJ.out.tab

if [ -d "$OUTDIR/SortedBamFiles/logs" ]
then
    mv $OUTDIR/SortedBamFiles/*Log* $OUTDIR/SortedBamFiles/logs
    echo "Log files moved to $OUTDIR/Ecoli_Aligned/SortedBamFiles/logs"
else
  mkdir $OUTDIR/SortedBamFiles/logs
  mv $OUTDIR/SortedBamFiles/*Log* $OUTDIR/SortedBamFiles/logs
  echo "Log files moved to $OUTDIR/SortedBamFiles/logs"
fi

echo "... loading SAMtools"
module load SAMtools
echo "... performing a MAPQ filter of 1"

for file in $OUTDIR/SortedBamFiles/"$outname"*sortedByCoord.out.bam
do
  base=$(basename ${file} Aligned.sortedByCoord.out.bam)
  #Our lab
  samtools view -bhSu -q 30 $file | samtools sort - > $OUTDIR/SortedBamFiles/"$base"ecoli_q30.bam

#Goll lab
  #samtools view -bq1 $file | samtools sort - > $OUTDIR/bams/"$base"_ecoli_q1.bam
done
echo "... quality filtering done"
echo

echo "... counting aligned reads"

for infile in $OUTDIR/SortedBamFiles/"$outname"*ecoli_q30.bam
do
  base=$(basename ${infile} ecoli_q30.bam)
  echo "$base total aligned reads -" >> $OUTDIR/SortedBamFiles/bam_stats.txt
  samtools view -@ 24 -F 0x4 $OUTDIR/SortedBamFiles/"$base"Aligned.sortedByCoord.out.bam | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/SortedBamFiles/bam_stats.txt
  echo "  $base total aligned reads (unique mappers) -" >> $OUTDIR/SortedBamFiles/bam_stats.txt
  samtools view -@ 24 -F 0x4 $OUTDIR/SortedBamFiles/"$base"Aligned.sortedByCoord.out.bam | grep "NH:i:1" | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/SortedBamFiles/bam_stats.txt
  echo "  $base total aligned reads (multi mappers) -" >> $OUTDIR/SortedBamFiles/bam_stats.txt
  samtools view -@ 24 -F 0x4 $OUTDIR/SortedBamFiles/"$base"Aligned.sortedByCoord.out.bam | grep -v "NH:i:1" | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/SortedBamFiles/bam_stats.txt
  echo "$base q1 aligned reads -" >> $OUTDIR/SortedBamFiles/bam_stats.txt
  samtools view -@ 24 -F 0x4 $infile | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/SortedBamFiles/bam_stats.txt
  echo "  $base q1 aligned reads (unique mappers) -" >> $OUTDIR/SortedBamFiles/bam_stats.txt
  samtools view -@ 24 -F 0x4 $infile | grep "NH:i:1" | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/SortedBamFiles/bam_stats.txt
  echo "  $base q1 aligned reads (multi mappers) -" >> $OUTDIR/SortedBamFiles/bam_stats.txt
  samtools view -@ 24 -F 0x4 $infile | grep -v "NH:i:1" | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/SortedBamFiles/bam_stats.txt
done

echo "... aligned reads counted and reported in $OUTDIR/SortedBamFiles/bam_stats"
echo

echo "$0 is complete"
