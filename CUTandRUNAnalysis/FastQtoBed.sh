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

echo "TrimmedReads reads and fastQC files going into $OUTDIR/TrimmedReads"
echo "... loading TrimGalore"
module load Trim_Galore
echo "...starting trimming"

trim_galore --fastqc -j 24 --length 20 --output_dir $OUTDIR/TrimmedReads --paired $read_file1 $read_file2

echo "... trimming complete"
echo "...loading MultiQC"
###need to edit this to only happen when the # of TrimmedReadsReads files = the number of input files?
module load MultiQC
echo "... running MultiQC"
multiqc -o $OUTDIR/TrimmedReads -f $OUTDIR/TrimmedReads
echo "... QC analysis complete"

echo
echo "... counting TrimmedReads reads"
for infile in $OUTDIR/TrimmedReads/"$outname"*R1*val*.fq.gz
do
  base=$(basename ${infile} _val_1.fq.gz)
  echo $base >> $OUTDIR/TrimmedReads/TrimmedReads_read_stats.txt
  echo $(zcat $infile |wc -l)/4|bc >> $OUTDIR/TrimmedReads/TrimmedReads_read_stats.txt
done
echo "... TrimmedReads reads counted and reported in $OUTDIR/TrimmedReads/TrimmedReads_read_stats"
echo

if [ -f "$genome" ]; then
  echo "genome file provided by user"
  echo "genome file is $genome"
elif [ -f "/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic.fna" ]; then
  echo "genome already downloaded"
  genome=/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic.fna
  echo "genome file is $genome"
else
  echo "no genome file supplied"
  echo "downloading z11 genome release 98"
#  curl -s ftp://ftp.ensembl.org/pub/release-98/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz | gunzip -c > $OUTDIR/danio_refseq.fa
  genome=/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic.fna
  echo "genome file is $genome"
fi

if [ -f "$ecoli" ]; then
  echo "genome file provided by user"
  echo "genome file is $ecoli"
elif [ -f "$OUTDIR/ecoli_refseq.fa" ]; then
  echo "genome already downloaded"
  ecoli=$OUTDIR/ecoli_refseq.fa
  echo "genome file is $ecoli"
else
  echo "no genome file supplied"
  echo "downloading ecoli genome"
  curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz | gunzip -c > $OUTDIR/ecoli_refseq.fa
  ecoli=$OUTDIR/ecoli_refseq.fa
  echo "genome file is $ecoli"

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
if [ -f "$OUTDIR/ref/Ncrassa_ref/Genome" ]
then
    echo "Reference genome index exists."
else
  echo "... building reference genome index"
  mkdir $OUTDIR/ref
  STAR --runThreadN 20 --runMode genomeGenerate \
  --genomeDir $OUTDIR/ref/Ncrassa_ref --genomeFastaFiles $genome
  STAR --runThreadN 20 --runMode genomeGenerate \
  --genomeDir $OUTDIR/ref/ecoli_ref --genomeFastaFiles $ecoli
  ml Bowtie2
  bowtie2-build -f ${ecoli} $OUTDIR/ref/ecoli_ref

  bowtie2-build -f $genome $OUTDIR/ref/Ncrassa_ref
echo "... aligning reads to "$genome" genome"
fi

FILES="$OUTDIR/TrimmedReads/*_R1_001_val_1.fq.gz" #Don't forget the *

for f in $FILES
 do
base=$(basename "${f}" _R1_001_val_1.fq.gz)
read2="$OUTDIR/TrimmedReads/${base}*_R2_001_val_2.fq.gz"
 #mkdir $OUTDIR/sam_files
 #mkdir "$OUTDIR/EColi_Aligned/sam_files"
 ml Bowtie2
bowtie2 --local --very-sensitive-local --phred33 --no-unal -p 24 -x "$OUTDIR/ref/Ncrassa_ref" -1 $f -2 $read2 -S "$OUTDIR/sam_files/${base}.sam"
bowtie2 --local --very-sensitive-local --phred33 --no-unal -p 24 -x "$OUTDIR/ref/ecoli_ref" -1 $f -2 $read2 -S "$OUTDIR/EColi_Aligned/sam_files/${base}_Ecoli.sam"
module load SAMtools
#  mkdir "$OUTDIR/bam_files"
#  mkdir "$OUTDIR/EColi_Aligned/bam_files"

 samtools view -bS -h -bq 30 "$OUTDIR/sam_files/${base}.sam" > "$OUTDIR/bam_files/${base}.bam"
 samtools view -bS -h -bq 30 "$OUTDIR/EColi_Aligned/sam_files/${base}_Ecoli.sam" > "$OUTDIR/EColi_Aligned/bam_files/${base}_Ecoli.bam"

#  mkdir "$OUTDIR/SortedBamFiles"
#  mkdir "$OUTDIR/EColi_Aligned/SortedBamFiles"

 samtools sort "$OUTDIR/bam_files/${base}.bam" -o "$OUTDIR/SortedBamFiles/${base}.sorted_q30.bam"
 samtools sort "$OUTDIR/EColi_Aligned/bam_files/${base}_Ecoli.bam" -o "$OUTDIR/EColi_Aligned/SortedBamFiles/${base}_Ecoli.sorted_q30.bam"

 samtools flagstat $OUTDIR/SortedBamFiles/${base}.sorted_q30.bam >> ${OUTDIR}/SortedBamFiles/${base}_flagstat.txt
 samtools flagstat $OUTDIR/EColi_Aligned/SortedBamFiles/${base}_Ecoli.sorted_q30.bam >> ${OUTDIR}/EColi_Aligned/SortedBamFiles/${base}_flagstat.txt

#bam="${OUTDIR}/bam_files/${base}.bam"
#mkdir ${OUTDIR}/BigWigs
#mkdir ${OUTDIR}/EColi_Aligned/BigWigs

 #	bigwig="${OUTDIR}/BigWigs/${base}"
  #ECbigwig="${OUTDIR}/EColi_Aligned/BigWigs/${base}"
  #samtools index "$bam"
#samtools index $OUTDIR/EColi_Aligned/bam_files/${base}_Ecoli.bam
#ml deepTools
#bamCoverage -p 12 -bs 25 --normalizeUsing BPM --smoothLength 50 -of bigwig -b "$bam" -o "${bigwig}.bin12.smooth25.Bulk.bw"
#bamCoverage -p 12 -bs 25 --normalizeUsing BPM --smoothLength 50 -of bigwig -b $OUTDIR/EColi_Aligned/bam_files/${base}_Ecoli.bam -o "${ECbigwig}.bin12.smooth25.Bulk.bw"
done

for infile in $OUTDIR/SortedBamFiles/*.sorted_q30.bam
do
  base=$(basename ${infile} .sorted_q30.bam)
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
ml picard
module load SAMtools
  for infile in $OUTDIR/SortedBamFiles/*.sorted_q30.bam
  do
    base=$(basename ${infile} .sorted_q30.bam)
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I $infile -M $OUTDIR/bams/"$base"_dupmetrics.txt -O $OUTDIR/SortedBamFiles/"$base"_nodups.bam --REMOVE_DUPLICATES true
    ml BEDTools
#    mkdir $OUTDIR/bed_files
    bedtools bamtobed -i $infile | awk -v OFS='\t' '{len = $3 - $2; print $0, len }' > $OUTDIR/bed_files/$base.btb.bed
done
#
# #Ncrassa nodups.bam conversion to bed files to call peaks and then pass through KMetSpikeIn Script
ml picard
module load SAMtools
  for infile in $OUTDIR/EColi_Aligned/SortedBamFiles/*.sorted_q30.bam
  do
    base=$(basename ${infile} .sorted_q30.bam)
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I $infile -M $OUTDIR/EColi_Aligned/SortedBamFiles/"$base"_dupmetrics.txt -O $OUTDIR/EColi_Aligned/SortedBamFiles/"$base"_nodups.bam --REMOVE_DUPLICATES true
    ml BEDTools
#    mkdir $OUTDIR/EColi_Aligned/bed_files
    bedtools bamtobed -i $infile | awk -v OFS='\t' '{len = $3 - $2; print $0, len }' > $OUTDIR/EColi_Aligned/bed_files/$base.btb.bed
done

 #ml BEDTools
 # #Ecoli nodups.bam
#Feeding output of no duplicates to DNA Spike In normalization by taking the bam to convert to the bed and get a bedgraph
#for infile in $OUTDIR/Ecoli_Aligned/*_ecoli_nodups.bam
#do
# base=$(basename ${infile} _ecoli_nodups.bam)
 #bedtools bamtobed -i $infile | awk -v OFS='\t' '{len = $3 - $2; print $0, len }' > $OUTDIR/bed_files/${base}_ecoli.btb.bed
#  done


# for file in $OUTDIR/TrimmedReads/"$outname"*_val_*.fq.gz;
# do
#   if [[ $prefix ]]; then
#     if [ $analysis_mode == 'NONE' ] ; then
#       STAR --runThreadN 24 --genomeDir $OUTDIR/genome --outFileNamePrefix $OUTDIR/SortedBamFiles/"$outname" \
#       --readFilesCommand zcat --readFilesIn "$first" "$file" --outSAMtype BAM SortedByCoordinate \
#       --outFilterMultimapNmax 1 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000
#     elif [ $analysis_mode == 'ONE' ] ; then
#       STAR --runThreadN 24 --genomeDir $OUTDIR/genome --outFileNamePrefix $OUTDIR/SortedBamFiles/"$outname" \
#       --readFilesCommand zcat --readFilesIn "$first" "$file" --outSAMtype BAM SortedByCoordinate \
#       --outSAMmultNmax 1 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000
#     elif [ $analysis_mode == 'RANDOM'] ; then
#     STAR --runThreadN 24 --genomeDir $OUTDIR/genome --outFileNamePrefix $OUTDIR/SortedBamFiles/"$outname" \
#     --readFilesCommand zcat --readFilesIn "$first" "$file" --outSAMtype BAM SortedByCoordinate \
#     --outMultimapperOrder Random --outSAMmultNmax 1 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000
#     elif [ $analysis_mode == 'ALL' ] ; then
#       STAR --runThreadN 24 --genomeDir $OUTDIR/genome --outFileNamePrefix $OUTDIR/SortedBamFiles/"$outname" \
#       --readFilesCommand zcat --readFilesIn "$first" "$file" --outSAMtype BAM SortedByCoordinate \
#       --outSAMprimaryFlag AllBestScore --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000
#     fi
#     echo "$outname is finished aligning"
#     prefix=
#   else
#     first=$file
#     prefix=${file%%_*}
#   fi
# done
#
# rm $OUTDIR/SortedBamFiles/"$outname"*SJ.out.tab
#
# if [ -d "$OUTDIR/SortedBamFiles/logs" ]
# then
#     mv $OUTDIR/SortedBamFiles/*Log* $OUTDIR/SortedBamFiles/logs
#     echo "Log files moved to $OUTDIR/SortedBamFiles/logs"
# else
#   mkdir $OUTDIR/SortedBamFiles/logs
#   mv $OUTDIR/SortedBamFiles/*Log* $OUTDIR/SortedBamFiles/logs
#   echo "Log files moved to $OUTDIR/SortedBamFiles/logs"
# fi
#
# echo "... loading SAMtools"
# module load SAMtools
# echo "... performing a MAPQ filter of 1"
#
# for file in $OUTDIR/SortedBamFiles/"$outname"*sortedByCoord.out.bam
# do
#   base=$(basename ${file} Aligned.sortedByCoord.out.bam)
#   samtools view -bS -h -bq 30 $file | samtools sort - > $OUTDIR/SortedBamFiles/"$base"_q30.bam
# done
echo "... quality filtering done"
echo

echo "... counting aligned reads"

# for infile in $OUTDIR/SortedBamFiles/"$outname"*_q1.bam
# do
#   base=$(basename ${infile} _q1.bam)
#   echo "$base total aligned reads -" >> $OUTDIR/SortedBamFiles/bam_stats.txt
#   samtools view -@ 24 -F 0x4 $OUTDIR/SortedBamFiles/"$base"Aligned.sortedByCoord.out.bam | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/SortedBamFiles/bam_stats.txt
#   echo "  $base total aligned reads (unique mappers) -" >> $OUTDIR/SortedBamFiles/bam_stats.txt
#   samtools view -@ 24 -F 0x4 $OUTDIR/SortedBamFiles/"$base"Aligned.sortedByCoord.out.bam | grep "NH:i:1" | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/SortedBamFiles/bam_stats.txt
#   echo "  $base total aligned reads (multi mappers) -" >> $OUTDIR/SortedBamFiles/bam_stats.txt
#   samtools view -@ 24 -F 0x4 $OUTDIR/SortedBamFiles/"$base"Aligned.sortedByCoord.out.bam | grep -v "NH:i:1" | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/SortedBamFiles/bam_stats.txt
#   echo "$base q1 aligned reads -" >> $OUTDIR/SortedBamFiles/bam_stats.txt
#   samtools view -@ 24 -F 0x4 $infile | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/SortedBamFiles/bam_stats.txt
#   echo "  $base q1 aligned reads (unique mappers) -" >> $OUTDIR/SortedBamFiles/bam_stats.txt
#   samtools view -@ 24 -F 0x4 $infile | grep "NH:i:1" | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/SortedBamFiles/bam_stats.txt
#   echo "  $base q1 aligned reads (multi mappers) -" >> $OUTDIR/SortedBamFiles/bam_stats.txt
#   samtools view -@ 24 -F 0x4 $infile | grep -v "NH:i:1" | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/SortedBamFiles/bam_stats.txt
# done

echo "... aligned reads counted and reported in $OUTDIR/SortedBamFiles/bam_stats"
echo

echo "$0 is complete"
