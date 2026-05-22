#!/bin/bash
#SBATCH --job-name=EAF3_fulltranscriptome
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --time=02:00:00
#SBATCH --output=EAF3_fulltranscriptome.%j.out
#SBATCH --error=EAF3_fulltranscriptome.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu

ml GCC/12.3.0 SAMtools/1.23.1-GCC-13.3.0 StringTie/3.0.0-GCC-13.3.0 gffread/0.12.7-GCCcore-12.3.0 TransDecoder/5.7.1-GCC-12.3.0

GENOME="/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic.fna"
GTF="/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic.gff"
OUTDIR="/scratch/ry00555/EpigeneticMemoryPaper2026/RNASeq/EAF3_assembly"
WTBAMs="/lustre2/scratch/ry00555/Run155/SortedBamFiles/"

WT1="$WTBAMs/155-93_RNA_WT__Rep1_S93_L002/155-93_RNA_WT__Rep1_S93_L002_Aligned.sortedByCoord.out.bam"
WT2="$WTBAMs/155-94_RNA_WT__Rep2_S94_L002/155-94_RNA_WT__Rep2_S94_L002_Aligned.sortedByCoord.out.bam"
WT3="$WTBAMs/155-95_RNA_WT__Rep3_S95_L002/155-95_RNA_WT__Rep3_S95_L002_Aligned.sortedByCoord.out.bam"

# Step 1: merge ALL three WT BAMs (full genome, no region extraction)
samtools merge -f -@ 8 \
    $OUTDIR/WT_fulltranscriptome_merged.bam \
    $WT1 $WT2 $WT3
samtools index $OUTDIR/WT_fulltranscriptome_merged.bam

# Step 2: assemble full transcriptome — StringTie uses all reads
# -m 50 minimum transcript length
# -c 2 minimum coverage
stringtie $OUTDIR/WT_fulltranscriptome_merged.bam \
    -G $GTF \
    -o $OUTDIR/WT_fulltranscriptome_assembly.gtf \
    -p 8 -v \
    -m 50 \
    -c 2

# Step 3: extract ALL transcript sequences
gffread $OUTDIR/WT_fulltranscriptome_assembly.gtf \
    -g $GENOME \
    -w $OUTDIR/WT_fulltranscriptome_transcripts.fa

echo "Done — full transcriptome assembled"

# Step 4: extract JUST the EAF3 locus transcripts from the full assembly
# Filter the GTF for just the EAF3 region
awk '$1 == "CM002237.1" && $4 >= 1930000 && $5 <= 1940000' \
    $OUTDIR/WT_fulltranscriptome_assembly.gtf > $OUTDIR/EAF3_locus_fromfull.gtf

gffread $OUTDIR/EAF3_locus_fromfull.gtf \
    -g $GENOME \
    -w $OUTDIR/EAF3_locus_fromfull_transcripts.fa

echo "EAF3 locus transcripts extracted"
cat $OUTDIR/EAF3_locus_fromfull_transcripts.fa

# Step 5: TransDecoder on just the EAF3 transcripts
cd $OUTDIR

TransDecoder.LongOrfs -t EAF3_locus_fromfull_transcripts.fa -m 50
TransDecoder.Predict -t EAF3_locus_fromfull_transcripts.fa --no_refine_starts

cp EAF3_locus_fromfull_transcripts.fa.transdecoder.pep EAF3_locus_fromfull_proteins.fa

echo ""
echo "Predicted proteins:"
grep ">" $OUTDIR/EAF3_locus_fromfull_proteins.fa

echo ""
echo "Partial ORFs:"
grep "5prime_partial" $OUTDIR/EAF3_locus_fromfull_proteins.fa | wc -l

echo ""
echo "Protein lengths:"
grep ">" $OUTDIR/EAF3_locus_fromfull_proteins.fa | grep -oP "len:\d+" | sort -t: -k2 -rn
