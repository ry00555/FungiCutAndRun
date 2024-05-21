# !/bin/bash
# SBATCH --job-name=CnR_setab_MO		                         Job name
# SBATCH --partition=batch		                             Partition (queue) name
# SBATCH --ntasks=1	                                 Single task job
# SBATCH --cpus-per-task=24		                             Number of cores per task - match this to the num_threads used by BLAST
# SBATCH --mem=120gb			                                 Total memory for job
# SBATCH --time=48:00:00  		                             Time limit hrs:min:sec
# SBATCH --output=/scratch/ara67776/CnR_setAB			     Location of standard output and error log files (replace cbergman with your myid)
# SBATCH --mail-user=ara67776@uga.edu                     Where to send mail (replace cbergman with your myid)
# SBATCH --mail-type=ALL                             Mail events (BEGIN, END, FAIL, ALL)
# SBATCH --output=/home/ara67776/work/error/log.%j			     Location of standard output and error log files (replace cbergman with your myid)
#
# set output directory
OUTDIR="/scratch/ara67776/CnR_timecourse"

 ml STAR
 for file in $OUTDIR/raw/CnR_timecourse_K9_2_to_4h/*_R*.fastq.gz;
 do
   if [[ $prefix ]]; then
         base=$(basename ${first} _R1.fastq.gz)
         sh /home/ara67776/scripts/PE_trim_and_star.sh -o $OUTDIR -n $base -m one $first $file
         prefix=
     else
         first=$file
         prefix=${file%%_*}
     fi
 done

 #aligning to ecoli genome
 ml STAR
 for file in $OUTDIR/raw/CnR_timecourse_K9_2_to_4h/*_R*.fastq.gz;
 do
   if [[ $prefix ]]; then
         base=$(basename ${first} _R1.fastq.gz)
         sh /home/ara67776/scripts/PE_trim_and_star_2.sh -o $OUTDIR -n $base -m one $first $file
         prefix=
     else
         first=$file
         prefix=${file%%_*}
     fi
 done

# aligning to ecoli genome
 curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz | gunzip -c > $OUTDIR/ecoli_refseq.fa
#  note here that STAR suggests SAindex = 10 but that makes the alignment FAIL, do 8 instead
 STAR --runThreadN 20 --genomeSAindexNbases 8 --runMode genomeGenerate --genomeDir $OUTDIR/ecoli_genome --genomeFastaFiles $OUTDIR/ecoli_refseq.fa

 for file in $OUTDIR/trimmed/*_val_*.fq.gz;
 do
   if [[ $prefix ]]; then
         base=$(basename ${first} _R1_val_1.fq.gz)
         STAR --runThreadN 20 --genomeDir $OUTDIR/ecoli_genome --outFileNamePrefix $OUTDIR/bams/"$base"_ecoli \
         --readFilesCommand zcat --readFilesIn "$first" "$file" --outSAMtype BAM SortedByCoordinate \
         --outSAMmultNmax 1 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000

         STAR --runThreadN 20 --genomeDir $OUTDIR/genome --outFileNamePrefix $OUTDIR/bams/"$base"_ecoli \
       --readFilesCommand zcat --readFilesIn "$first" "$file" --outSAMtype BAM SortedByCoordinate \
       --outSAMmultNmax 1 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000

       STAR --runThreadN 20 --genomeDir $OUTDIR/genome --outFileNamePrefix $OUTDIR/bams/"$base"_ecoli \
     --readFilesCommand zcat --readFilesIn "$first" "$file" --outSAMtype BAM SortedByCoordinate \
     --outMultimapperOrder Random --outSAMmultNmax 1 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000

     STAR --runThreadN 20 --genomeDir $OUTDIR/genome --outFileNamePrefix $OUTDIR/bams/"$base"_ecoli \
       --readFilesCommand zcat --readFilesIn "$first" "$file" --outSAMtype BAM SortedByCoordinate \
       --outSAMprimaryFlag AllBestScore --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000
         prefix=
     else
         first=$file
         prefix=${file%%_*}
     fi
 done

 rm $OUTDIR/bams/"$base"*SJ.out.tab

 if [ -d "$OUTDIR/bams/logs" ]
 then
     mv $OUTDIR/bams/*Log* $OUTDIR/bams/logs
 else
   mkdir $OUTDIR/bams/logs
   mv $OUTDIR/bams/*Log* $OUTDIR/bams/logs
 fi

 module load SAMtools

 for file in $OUTDIR/bams/"$base"*ecoliAligned.sortedByCoord.out.bam
 do
   base=$(basename ${file} ecoliAligned.sortedByCoord.out.bam)
   samtools view -bq1 $file | samtools sort - > $OUTDIR/bams/"$base"_ecoli_q1.bam
 done

 for infile in $OUTDIR/bams/"$base"*_ecoli_q1.bam
 do
   base=$(basename ${infile} _ecoli_q1.bam)
   echo "$base total aligned reads -" >> $OUTDIR/bams/bam_stats.txt
   samtools view -@ 24 -F 0x4 $OUTDIR/bams/"$base"ecoliAligned.sortedByCoord.out.bam | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/bams/bam_stats.txt
   echo "  $base total aligned reads (unique mappers) -" >> $OUTDIR/bams/bam_stats.txt
   samtools view -@ 24 -F 0x4 $OUTDIR/bams/"$base"ecoliAligned.sortedByCoord.out.bam | grep "NH:i:1" | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/bams/bam_stats.txt
   echo "  $base total aligned reads (multi mappers) -" >> $OUTDIR/bams/bam_stats.txt
   samtools view -@ 24 -F 0x4 $OUTDIR/bams/"$base"ecoliAligned.sortedByCoord.out.bam | grep -v "NH:i:1" | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/bams/bam_stats.txt
   echo "$base q1 aligned reads -" >> $OUTDIR/bams/bam_stats.txt
   samtools view -@ 24 -F 0x4 $infile | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/bams/bam_stats.txt
   echo "  $base q1 aligned reads (unique mappers) -" >> $OUTDIR/bams/bam_stats.txt
   samtools view -@ 24 -F 0x4 $infile | grep "NH:i:1" | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/bams/bam_stats.txt
   echo "  $base q1 aligned reads (multi mappers) -" >> $OUTDIR/bams/bam_stats.txt
   samtools view -@ 24 -F 0x4 $infile | grep -v "NH:i:1" | cut -f 1 | sort | uniq | wc -l >> $OUTDIR/bams/bam_stats.txt
 done


 #Remove PCR duplicates
 ml picard
 module load SAMtools
 for infile in $OUTDIR/bams/*_q1.bam
 do
   base=$(basename ${infile} _q1.bam)
   java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I $infile -M $OUTDIR/bams/"$base"_dupmetrics.txt -O $OUTDIR/bams/"$base"_nodups.bam --REMOVE_DUPLICATES true
 done


 for infile in $OUTDIR/bams/*_ecoli_q1.bam
 do
   base=$(basename ${infile} _ecoli_q1.bam)
   java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I $infile -M $OUTDIR/bams/"$base"_dupmetrics.txt -O $OUTDIR/bams/"$base"_ecoli_nodups.bam --REMOVE_DUPLICATES true
 done

 #merging IgG samples from all time points to create uniformity in peak calling later
 ml SAMtools
 samtools merge $OUTDIR/bams/2hpf_IgG_merge_nodup.bam $OUTDIR/bams/2hpf_IgG_1_nodups.bam $OUTDIR/bams/3.5hpf_IgG_1_nodups.bam
 samtools merge $OUTDIR/bams/2.5hpf_IgG_merge_nodup.bam $OUTDIR/bams/2.5hpf_IgG_1_nodups.bam $OUTDIR/bams/2.5hpf_IgG_2_nodups.bam
 samtools merge $OUTDIR/bams/3hpf_IgG_merge_nodup.bam $OUTDIR/bams/3hpf_IgG_1_nodups.bam $OUTDIR/bams/3hpf_IgG_2_nodups.bam
 samtools merge $OUTDIR/bams/4hpf_IgG_merge_nodup.bam $OUTDIR/bams/4hpf_IgG_1_nodups.bam $OUTDIR/bams/4hpf_IgG_2_nodups.bam
 samtools merge $OUTDIR/bams/4.5hpf_IgG_merge_nodup.bam $OUTDIR/bams/4.5hpf_IgG_1_nodups.bam $OUTDIR/bams/4.5hpf_IgG_2_nodups.bam

 samtools merge $OUTDIR/bams/2hpf_IgG_ecoli_merge_nodup.bam $OUTDIR/bams/2hpf_IgG_1_ecoli_nodups.bam $OUTDIR/bams/3.5hpf_IgG_1_ecoli_nodups.bam
 samtools merge $OUTDIR/bams/2.5hpf_IgG_ecoli_merge_nodup.bam $OUTDIR/bams/2.5hpf_IgG_1_ecoli_nodups.bam $OUTDIR/bams/2.5hpf_IgG_2_ecoli_nodups.bam
 samtools merge $OUTDIR/bams/3hpf_IgG_ecoli_merge_nodup.bam $OUTDIR/bams/3hpf_IgG_1_ecoli_nodups.bam $OUTDIR/bams/3hpf_IgG_2_ecoli_nodups.bam
 samtools merge $OUTDIR/bams/4hpf_IgG_ecoli_merge_nodup.bam $OUTDIR/bams/4hpf_IgG_1_ecoli_nodups.bam $OUTDIR/bams/4hpf_IgG_2_ecoli_nodups.bam
 samtools merge $OUTDIR/bams/4.5hpf_IgG_ecoli_merge_nodup.bam $OUTDIR/bams/4.5hpf_IgG_1_ecoli_nodups.bam $OUTDIR/bams/4.5hpf_IgG_2_ecoli_nodups.bam

 samtools merge $OUTDIR/bams/IgG_danio_nodups.bam $OUTDIR/bams/2hpf_IgG_merge_nodup.bam $OUTDIR/bams/2.5hpf_IgG_merge_nodup.bam $OUTDIR/bams/3hpf_IgG_merge_nodup.bam $OUTDIR/bams/4hpf_IgG_merge_nodup.bam $OUTDIR/bams/4.5hpf_IgG_merge_nodup.bam
 samtools merge $OUTDIR/bams/IgG_ecoli_nodups.bam $OUTDIR/bams/2hpf_IgG_ecoli_merge_nodup.bam $OUTDIR/bams/2.5hpf_IgG_ecoli_merge_nodup.bam $OUTDIR/bams/3hpf_IgG_ecoli_merge_nodup.bam $OUTDIR/bams/4hpf_IgG_ecoli_merge_nodup.bam $OUTDIR/bams/4.5hpf_IgG_ecoli_merge_nodup.bam

 #Now we need to extract all the aligned reads in preperation for spike in normalization
 ml BEDTools
 #DANIO nodups.bam
 for infile in $OUTDIR/bams/*nodups.bam
 do
   base=$(basename ${infile} .bam)
   bedtools bamtobed -i $infile | awk -v OFS='\t' '{len = $3 - $2; print $0, len }' > $OUTDIR/bams/$base.btb.bed
 done
 #Ecoli nodups.bam
 for infile in $OUTDIR/bams/*_ecoli_nodups.bam
 do
   base=$(basename ${infile} _ecoli_nodups.bam)
   bedtools bamtobed -i $infile | awk -v OFS='\t' '{len = $3 - $2; print $0, len }' > $OUTDIR/bams/${base}_ecoli.btb.bed
 done

 #spike in normalization
 ml BEDTools
 mkdir $OUTDIR/bdgrphs
 sh /home/ara67776/scripts/DNAspike_in.kd.sh $OUTDIR/bams/IgG_danio_nodups.btb.bed $OUTDIR/bams/IgG_ecoli_nodups.btb.bed 100000 bga $OUTDIR/genome/chrNameLength.txt 1 1000 $OUTDIR/bdgrphs/timecourse_IgG_merged.norm.bga
 sh /home/ara67776/scripts/DNAspike_in.kd.sh $OUTDIR/bams/IgG_danio_nodups.btb.bed $OUTDIR/bams/IgG_ecoli_nodups.btb.bed 100000 bga $OUTDIR/genome/chrNameLength.txt 1 1000 $OUTDIR/bdgrphs/timecourse_IgG_merged.norm.bga
 sh /home/ara67776/scripts/DNAspike_in.kd.sh $OUTDIR/bams/2hpf_K9_1_nodups.btb.bed $OUTDIR/bams/2hpf_K9_1_ecoli_nodups.btb.bed 100000 bga $OUTDIR/genome/chrNameLength.txt 1 1000 $OUTDIR/bdgrphs/2hpf_K9_1.norm.bga
 sh /home/ara67776/scripts/DNAspike_in.kd.sh $OUTDIR/bams/2hpf_K9_2_nodups.btb.bed $OUTDIR/bams/2hpf_K9_2_ecoli_nodups.btb.bed 100000 bga $OUTDIR/genome/chrNameLength.txt 1 1000 $OUTDIR/bdgrphs/2hpf_K9_2.norm.bga
 sh /home/ara67776/scripts/DNAspike_in.kd.sh $OUTDIR/bams/2.5hpf_K9_1_nodups.btb.bed $OUTDIR/bams/2.5hpf_K9_1_ecoli.btb.bed 100000 bga $OUTDIR/genome/chrNameLength.txt 1 1000 $OUTDIR/bdgrphs/2.5hpf_K9_1.norm.bga
 sh /home/ara67776/scripts/DNAspike_in.kd.sh $OUTDIR/bams/2.5hpf_K9_2_nodups.btb.bed $OUTDIR/bams/2.5hpf_K9_2_ecoli.btb.bed 100000 bga $OUTDIR/genome/chrNameLength.txt 1 1000 $OUTDIR/bdgrphs/2.5hpf_K9_2.norm.bga
 sh /home/ara67776/scripts/DNAspike_in.kd.sh $OUTDIR/bams/2.5hpf_K9_3_nodups.btb.bed $OUTDIR/bams/2.5hpf_K9_3_ecoli.btb.bed 100000 bga $OUTDIR/genome/chrNameLength.txt 1 1000 $OUTDIR/bdgrphs/2.5hpf_K9_3.norm.bga
 sh /home/ara67776/scripts/DNAspike_in.kd.sh $OUTDIR/bams/3hpf_K9_1_nodups.btb.bed $OUTDIR/bams/3hpf_K9_1_ecoli_nodups.btb.bed 100000 bga $OUTDIR/genome/chrNameLength.txt 1 1000 $OUTDIR/bdgrphs/3hpf_K9_1.norm.bga
 sh /home/ara67776/scripts/DNAspike_in.kd.sh $OUTDIR/bams/3hpf_K9_2_nodups.btb.bed $OUTDIR/bams/3hpf_K9_2_ecoli_nodups.btb.bed 100000 bga $OUTDIR/genome/chrNameLength.txt 1 1000 $OUTDIR/bdgrphs/3hpf_K9_2.norm.bga
 sh /home/ara67776/scripts/DNAspike_in.kd.sh $OUTDIR/bams/3hpf_K9_3_nodups.btb.bed $OUTDIR/bams/3hpf_K9_3_ecoli_nodups.btb.bed 100000 bga $OUTDIR/genome/chrNameLength.txt 1 1000 $OUTDIR/bdgrphs/3hpf_K9_3.norm.bga
 sh /home/ara67776/scripts/DNAspike_in.kd.sh $OUTDIR/bams/3.5hpf_K9_1_nodups.btb.bed $OUTDIR/bams/3.5hpf_K9_1_ecoli_nodups.btb.bed 100000 bga $OUTDIR/genome/chrNameLength.txt 1 1000 $OUTDIR/bdgrphs/3.5hpf_K9_1.norm.bga
 sh /home/ara67776/scripts/DNAspike_in.kd.sh $OUTDIR/bams/3.5hpf_K9_2_nodups.btb.bed $OUTDIR/bams/3.5hpf_K9_2_ecoli_nodups.btb.bed 100000 bga $OUTDIR/genome/chrNameLength.txt 1 1000 $OUTDIR/bdgrphs/3.5hpf_K9_2.norm.bga
 sh /home/ara67776/scripts/DNAspike_in.kd.sh $OUTDIR/bams/3.5hpf_K9_3_nodups.btb.bed $OUTDIR/bams/3.5hpf_K9_3_ecoli_nodups.btb.bed 100000 bga $OUTDIR/genome/chrNameLength.txt 1 1000 $OUTDIR/bdgrphs/3.5hpf_K9_3.norm.bga
 sh /home/ara67776/scripts/DNAspike_in.kd.sh $OUTDIR/bams/4hpf_K9_1_nodups.btb.bed $OUTDIR/bams/4hpf_K9_1_ecoli_nodups.btb.bed 100000 bga $OUTDIR/genome/chrNameLength.txt 1 1000 $OUTDIR/bdgrphs/4hpf_K9_1.norm.bga
 sh /home/ara67776/scripts/DNAspike_in.kd.sh $OUTDIR/bams/4hpf_K9_2_nodups.btb.bed $OUTDIR/bams/4hpf_K9_2_ecoli_nodups.btb.bed 100000 bga $OUTDIR/genome/chrNameLength.txt 1 1000 $OUTDIR/bdgrphs/4hpf_K9_2.norm.bga
 sh /home/ara67776/scripts/DNAspike_in.kd.sh $OUTDIR/bams/4hpf_K9_3_nodups.btb.bed $OUTDIR/bams/4hpf_K9_3_ecoli_nodups.btb.bed 100000 bga $OUTDIR/genome/chrNameLength.txt 1 1000 $OUTDIR/bdgrphs/4hpf_K9_3.norm.bga
 sh /home/ara67776/scripts/DNAspike_in.kd.sh $OUTDIR/bams/4.5hpf_K9_1_nodups.btb.bed $OUTDIR/bams/4.5hpf_K9_1_ecoli.btb.bed 100000 bga $OUTDIR/genome/chrNameLength.txt 1 1000 $OUTDIR/bdgrphs/4.5hpf_K9_1.norm.bga
 sh /home/ara67776/scripts/DNAspike_in.kd.sh $OUTDIR/bams/4.5hpf_K9_2_nodups.btb.bed $OUTDIR/bams/4.5hpf_K9_2_ecoli.btb.bed 100000 bga $OUTDIR/genome/chrNameLength.txt 1 1000 $OUTDIR/bdgrphs/4.5hpf_K9_2.norm.bga
 sh /home/ara67776/scripts/DNAspike_in.kd.sh $OUTDIR/bams/4.5hpf_K9_3_nodups.btb.bed $OUTDIR/bams/4.5hpf_K9_3_ecoli.btb.bed 100000 bga $OUTDIR/genome/chrNameLength.txt 1 1000 $OUTDIR/bdgrphs/4.5hpf_K9_3.norm.bga

 #peak calling
module load Homer/
 mkdir $OUTDIR/peaks_5

 for infile in $OUTDIR/bdgrphs/*.norm.bga
   do base=$(basename ${infile} .norm.bga)
   cat $infile | awk '{print $1 "\t" $2 "\t" $3 "\t" "+" "\t" "+" "\t" "+"}' > $OUTDIR/peaks_5/${base}.bgato.bed
 done

 for infile in $OUTDIR/peaks_5/*bgato.bed
   do base=$(basename ${infile} .bgato.bed)
   makeTagDirectory $OUTDIR/peaks_5/${base}.BtB.tagdir $infile -format bed
 done

 for infile in $OUTDIR/peaks_5/*K9*.tagdir w/ IgG as an input
   do base=$(basename ${infile} .BtB.tagdir)
   findPeaks $infile -style histone -minDist 1000 -gsize 1.5e9 -F 4 -i $OUTDIR/peaks_5/timecourse_IgG_merged.BtB.tagdir -o $OUTDIR/peaks_5/${base}.txt
 done

 for infile in $OUTDIR/peaks_5/*.txt
 do
   base=$(basename ${infile} .txt)
   sed '/^/d' $infile | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $8 "\t" $5 "\t" $6 "\t" $12 "\t" "-1"}' | sed 's/\.000000//g' > $OUTDIR/peaks_5/${base}.peaks.bed
 done

 #inter secting peaks across replicates using Chipr
 ml ChIP-R
 chipr -i $OUTDIR/peaks_5/2hpf_K9_1.peaks.bed $OUTDIR/peaks_5/2hpf_K9_2.peaks.bed $OUTDIR/peaks_5/2hpf_K9_3.peaks.bed -m 2 -o $OUTDIR/peaks_5/2hpf_K9_repPeaks
 chipr -i $OUTDIR/peaks_5/2.5hpf_K9_1.peaks.bed $OUTDIR/peaks_5/2.5hpf_K9_2.peaks.bed $OUTDIR/peaks_5/2.5hpf_K9_3.peaks.bed -m 2 -o $OUTDIR/peaks_5/2.5hpf_K9_repPeaks
 chipr -i $OUTDIR/peaks_5/3hpf_K9_1.peaks.bed $OUTDIR/peaks_5/3hpf_K9_2.peaks.bed $OUTDIR/peaks_5/3hpf_K9_3.peaks.bed -m 2 -o $OUTDIR/peaks_5/3hpf_K9_repPeaks
 chipr -i $OUTDIR/peaks_5/3.5hpf_K9_1.peaks.bed $OUTDIR/peaks_5/3.5hpf_K9_2.peaks.bed $OUTDIR/peaks_5/3.5hpf_K9_3.peaks.bed -m 2 -o $OUTDIR/peaks_5/3.5hpf_K9_repPeaks
 chipr -i $OUTDIR/peaks_5/4hpf_K9_1.peaks.bed $OUTDIR/peaks_5/4hpf_K9_2.peaks.bed $OUTDIR/peaks_5/4hpf_K9_3.peaks.bed -m 2 -o $OUTDIR/peaks_5/4hpf_K9_repPeaks
 chipr -i $OUTDIR/peaks_5/4.5hpf_K9_1.peaks.bed $OUTDIR/peaks_5/4.5hpf_K9_2.peaks.bed $OUTDIR/peaks_5/4.5hpf_K9_3.peaks.bed -m 2 -o $OUTDIR/peaks_5/4.5hpf_K9_repPeaks

 #make a blacklist file
 findPeaks $OUTDIR/peaks_5/timecourse_IgG_merged.BtB.tagdir -style factor -o $OUTDIR/peaks_5/IgG.txt
 sed '/^/d' $OUTDIR/peaks_5/IgG.txt | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" "1" "\t" $5 "\t" $6 "\t" $12 "\t" "-1"}' > $OUTDIR/peaks_5/blacklist.bed

 ml BEDTools
  #intersect the peaks with the blacklist file to make sure we aren't looking at sticky regions before this step
 for infile in $OUTDIR/peaks_5/*all.bed
 do
   base=$( basename ${infile} _repPeaks_all.bed)
   bedtools intersect -a $infile -b $OUTDIR/peaks_5/blacklist.bed -v > $OUTDIR/peaks_5/"${base}"_final.bed
 done

 #peak annotation
 curl -s ftp://ftp.ensembl.org/pub/release-98/gtf/danio_rerio/Danio_rerio.GRCz11.98.gtf.gz | gunzip -c > $OUTDIR/refann.gtf
 mkdir $OUTDIR/peaks_5/ann

 for infile in $OUTDIR/peaks_5/*final.bed
 do
   base=$( basename ${infile} final.bed)
   annotatePeaks.pl $infile danRer11 -gtf $OUTDIR/refann.gtf > $OUTDIR/peaks_5/ann/${base}.maskann.txt
 done

 for infile in $OUTDIR/peaks_5/ann/*maskann.txt
 do
   base=$(basename ${infile} .maskann.txt)
   awk -F'\t' 'sqrt($10*$10) <=1000' $infile > $OUTDIR/peaks_5/ann/${base}.1000bp_ann.txt
 done

 for infile in $OUTDIR/peaks_4/ann/*maskann.txt
 do
   base=$(basename ${infile} .maskann.txt)
   awk -F'\t' 'sqrt($10*$10) >=1000' $infile | awk '{print $2 "\t" $3 "\t" $4 }' > $OUTDIR/peaks_4/ann/${base}.MOREthan1000bp.bed
 done


# making bigwigs of the gene list curated
 #sort bga files from  DNA spike in
 ml ucsc
 for infile in $OUTDIR/bdgrphs/*norm.bga
 do
   base=$(basename ${infile} .norm.bga)
   bedSort $infile $OUTDIR/bdgrphs/${base}.norm_sort.bga
 done

 mkdir $OUTDIR/bigwigs_4
 ml ucsc
 for infile in $OUTDIR/bdgrphs/*.norm_sort.bga
 do
  base=$(basename ${infile} .norm_sort.bga)
 bedGraphToBigWig $infile $OUTDIR/genome/chrNameLength.txt $OUTDIR/bigwigs_4/${base}.bw
 done

#making bws of the mean of the replicates for profile plots
 ml deepTools
 bigwigCompare -b1 $OUTDIR/bigwigs_4/2hpf_K9_1.bw -b2 $OUTDIR/bigwigs_4/2hpf_K9_2.bw --operation add -bs 10 -p 20 -o $OUTDIR/bigwigs_4/2hpf_K9_AVG.bw

 bigwigCompare -b1 $OUTDIR/bigwigs_4/2.5hpf_K9_1.bw -b2 $OUTDIR/bigwigs_4/2.5hpf_K9_2.bw --operation add --scaleFactors 0.333:0.333 -bs 10 -p 20 -o $OUTDIR/bigwigs_4/2.5hpf_K9_rep1rep2.bw
 bigwigCompare -b1 $OUTDIR/bigwigs_4/2.5hpf_K9_rep1rep2.bw -b2 $OUTDIR/bigwigs_4/2.5hpf_K9_3.bw --operation add --scaleFactors 1:0.333 -bs 10 -p 20 -o $OUTDIR/bigwigs_4/2.5hpf_K9_AVG.bw

 bigwigCompare -b1 $OUTDIR/bigwigs_4/3hpf_K9_1.bw -b2 $OUTDIR/bigwigs_4/3hpf_K9_2.bw --operation add --scaleFactors 0.333:0.333 -bs 10 -p 20 -o $OUTDIR/bigwigs_4/3hpf_K9_rep1rep2.bw
 bigwigCompare -b1 $OUTDIR/bigwigs_4/3hpf_K9_rep1rep2.bw -b2 $OUTDIR/bigwigs_4/3hpf_K9_3.bw --operation add --scaleFactors 1:0.333 -bs 10 -p 20 -o $OUTDIR/bigwigs_4/3hpf_K9_AVG.bw

 bigwigCompare -b1 $OUTDIR/bigwigs_4/3.5hpf_K9_1.bw -b2 $OUTDIR/bigwigs_4/3.5hpf_K9_2.bw --operation add --scaleFactors 0.333:0.333 -bs 10 -p 20 -o $OUTDIR/bigwigs_4/3.5hpf_K9_rep1rep2.bw
 bigwigCompare -b1 $OUTDIR/bigwigs_4/3.5hpf_K9_rep1rep2.bw -b2 $OUTDIR/bigwigs_4/3.5hpf_K9_3.bw --operation add --scaleFactors 1:0.333 -bs 10 -p 20 -o $OUTDIR/bigwigs_4/3.5hpf_K9_AVG.bw

 bigwigCompare -b1 $OUTDIR/bigwigs_4/4hpf_K9_1.bw -b2 $OUTDIR/bigwigs_4/4hpf_K9_2.bw --operation add --scaleFactors 0.333:0.333 -bs 10 -p 20 -o $OUTDIR/bigwigs_4/4hpf_K9_rep1rep2.bw
 bigwigCompare -b1 $OUTDIR/bigwigs_4/4hpf_K9_rep1rep2.bw -b2 $OUTDIR/bigwigs_4/4hpf_K9_3.bw --operation add --scaleFactors 1:0.333 -bs 10 -p 20 -o $OUTDIR/bigwigs_4/4hpf_K9_AVG.bw

 bigwigCompare -b1 $OUTDIR/bigwigs_4/4.5hpf_K9_1.bw -b2 $OUTDIR/bigwigs_4/4.5hpf_K9_2.bw --operation add --scaleFactors 0.333:0.333 -bs 10 -p 20 -o $OUTDIR/bigwigs_4/4.5hpf_K9_rep1rep2.bw
 bigwigCompare -b1 $OUTDIR/bigwigs_4/4.5hpf_K9_rep1rep2.bw -b2 $OUTDIR/bigwigs_4/4.5hpf_K9_3.bw --operation add --scaleFactors 1:0.333 -bs 10 -p 20 -o $OUTDIR/bigwigs_4/4.5hpf_K9_AVG.bw

 #intersect TE annotation file with K9 "gene" peak files to identify genic H3K9 enrichment
 ml BEDTools
 mkdir $OUTDIR/peaks_5/TE_ann
 for infile in $OUTDIR/peaks_5/*final.bed
 do
   base=$( basename ${infile} _final.bed)
   bedtools intersect -a $infile -b $OUTDIR/final_beds/TEann_35_0.1filt.bed -f 0.50 -v > $OUTDIR/peaks_5/TE_ann/"${base}"_TEann_final.bed
 done
 mkdir $OUTDIR/peaks_4/TE_ann

 for infile in $OUTDIR/peaks_4/*TEann_final.bed
 do
   base=$( basename ${infile} TEann_final.bed)
   annotatePeaks.pl $infile danRer11 -gtf $OUTDIR/refann.gtf > $OUTDIR/peaks_4/TE_ann/${base}.TE_maskann.txt
 done

 for infile in $OUTDIR/peaks_4/TE_ann/*TE_maskann.txt
 do
   base=$(basename ${infile} .TE_maskann.txt)
   awk -F'\t' 'sqrt($10*$10) <=1000' $infile > $OUTDIR/peaks_4/TE_ann/${base}.1000bp_TEann.txt
 done


 #intersecting TE file with 4.5hpf annotated txt file
 for infile in $OUTDIR/peaks_4/ann/*ann.txt
 do
   base=$(basename ${infile} _ann.txt)
   cat $infile | awk '{print $2 "\t" $3 "\t" $4 "\t" "+" "\t" "+" "\t" "+"}' > $OUTDIR/peaks_4/ann/"${base}"_TEann.txt
 done

 #converting txt file to bed file for bedtools intersect input
 cut -f 1,2,3,4 $OUTDIR/peaks_4/ann/4.5hpf_K9_.1000bp_TEann.txt > $OUTDIR/peaks_4/ann/4.5hpf_K9_.1000bp_TEann.bed

 #intersecting TEann files with TE annotated file
 ml BEDTools
 bedtools intersect -a $OUTDIR/peaks_4/ann/4.5hpf_K9_genic_peaks.bed -b $OUTDIR/peaks_4/TEann_35_0.1filt.bed -v > $OUTDIR/peaks_4/ann/4.5hpf_K9_.1000bp_TEann_final.bed







#converting refann file into a bed file for hox genes
 ml BEDTools
 grep "hox" $OUTDIR/refann.gtf > $OUTDIR/refann_hox.bed

#creating heatmap
ml deepTools
computeMatrix scale-regions -S /scratch/ara67776/H3K9me3_rep1.bw /scratch/ara67776/H3K9me3_rep2.bw /scratch/ara67776/H3K27me3_rep1.bw /scratch/ara67776/H3K27me3_rep2.bw -R $OUTDIR/refann_hox.bed -b 10000 -out $OUTDIR/figures/K9_K27_hox_genes.gz
