#!/bin/bash
#SBATCH --job-name=EAF3_assembly
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=5gb
#SBATCH --time=00:30:00
#SBATCH --output=EAF3_assembly.%j.out
#SBATCH --error=EAF3_assembly.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu

ml GCC/12.3.0 SAMtools/1.23.1-GCC-13.3.0 StringTie/3.0.0-GCC-13.3.0 gffread/0.12.7-GCCcore-12.3.0 TransDecoder/5.7.1-GCC-12.3.0

GENOME="/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic.fna"
GTF="/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic.gff"
OUTDIR="/scratch/ry00555/EpigeneticMemoryPaper2026/RNASeq/EAF3_assembly"
WTBAMs="/lustre2/scratch/ry00555/Run155/SortedBamFiles/"
mkdir -p $OUTDIR

WT1="$WTBAMs/155-93_RNA_WT__Rep1_S93_L002/155-93_RNA_WT__Rep1_S93_L002_Aligned.sortedByCoord.out.bam"
WT2="$WTBAMs/155-94_RNA_WT__Rep2_S94_L002/155-94_RNA_WT__Rep2_S94_L002_Aligned.sortedByCoord.out.bam"
WT3="$WTBAMs/155-95_RNA_WT__Rep3_S95_L002/155-95_RNA_WT__Rep3_S95_L002_Aligned.sortedByCoord.out.bam"

# Step 1: extract EAF3 locus with 3kb upstream extension
for bam in $WT1 $WT2 $WT3; do
    name=$(basename $bam .bam)
    samtools view -b $bam "CM002237.1:1931000-1937500" > $OUTDIR/${name}_EAF3_5prime_3kbextension.bam
    samtools index $OUTDIR/${name}_EAF3_5prime_3kbextension.bam
done

# Step 2: merge only the 3kb extension BAMs
samtools merge -f -@ 8 \
    $OUTDIR/EAF3_WT_merged_5prime_3kbextension.bam \
    $OUTDIR/*_EAF3_5prime_3kbextension.bam
samtools index $OUTDIR/EAF3_WT_merged_5prime_3kbextension.bam

# Step 3: assemble transcripts
stringtie $OUTDIR/EAF3_WT_merged_5prime_3kbextension.bam \
    -G $GTF \
    -o $OUTDIR/EAF3_assembly_5prime_3kbextension.gtf \
    -p 8 -v

# Step 4: extract transcript sequences
gffread $OUTDIR/EAF3_assembly_5prime_3kbextension.gtf \
    -g $GENOME \
    -w $OUTDIR/EAF3_transcripts_5prime_3kbextension.fa

echo "Done — transcripts written to $OUTDIR/EAF3_transcripts_5prime_3kbextension.fa"

cd $OUTDIR

# Step 5: TransDecoder
TransDecoder.LongOrfs -t EAF3_transcripts_5prime_3kbextension.fa -m 50
TransDecoder.Predict -t EAF3_transcripts_5prime_3kbextension.fa --no_refine_starts

cp EAF3_transcripts_5prime_3kbextension.fa.transdecoder.pep EAF3_proteins_5prime_3kbextension.fa

echo "Done — protein sequences written to $OUTDIR/EAF3_proteins_5prime_3kbextension.fa"
echo ""
echo "Predicted proteins:"
grep ">" $OUTDIR/EAF3_proteins_5prime_3kbextension.fa

echo ""
echo "Partial ORFs remaining:"
grep "5prime_partial" $OUTDIR/EAF3_proteins_5prime_3kbextension.fa | wc -l

echo ""
echo "Protein lengths:"
grep ">" $OUTDIR/EAF3_proteins_5prime_3kbextension.fa | grep -oP "len:\d+" | sort -t: -k2 -rn
