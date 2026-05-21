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

ml SAMtools/1.23.1-GCC-13.3.0 StringTie/3.0.0-GCC-13.3.0 gffread/0.12.7-GCCcore-12.3.0 TransDecoder/5.7.1-GCC-12.3.0

GENOME="/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic.fna"
GTF="/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic.gff"
OUTDIR="/scratch/ry00555/EpigeneticMemoryPaper2026/RNASeq/EAF3_assembly"
WTBAMs="/lustre2/scratch/ry00555/Run155/SortedBamFiles/"
mkdir -p $OUTDIR

WT1="$WTBAMs/155-93_RNA_WT__Rep1_S93_L002/155-93_RNA_WT__Rep1_S93_L002_Aligned.sortedByCoord.out.bam"
WT2="$WTBAMs/155-94_RNA_WT__Rep2_S94_L002/155-94_RNA_WT__Rep2_S94_L002_Aligned.sortedByCoord.out.bam"
WT3="$WTBAMs/155-95_RNA_WT__Rep3_S95_L002/155-95_RNA_WT__Rep3_S95_L002_Aligned.sortedByCoord.out.bam"

# Step 1: extract just the EAF3 locus from each BAM (with padding)
for bam in $WT1 $WT2 $WT3; do
    name=$(basename $bam .bam)
    samtools view -b $bam "CM002237.1:1931000-1937500" > $OUTDIR/${name}_EAF3_5prime_3kbextension.bam # CM002237.1:1934000-1937500 is based on genome coordinates from fungidb, but the long isoform is missing a start codon, so run again but with 3kb upstream 5'end
    samtools index $OUTDIR/${name}_EAF3_5prime_3kbextension.bam
done

# Step 2: merge the 3 extracted BAMs
samtools merge -f -@ 8 \
    $OUTDIR/EAF3_WT_merged_5prime_3kbextension.bam  \
    $OUTDIR/*_EAF3_WT_merged_5prime_3kbextension.bam
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

#TransDecoder writes its output files to whatever directory you run it from (your submit directory), not to $OUTDIR
cd $OUTDIR
#figure out amino acid sequence using transdecoder identify long ORFs (minimum 50 aa to catch short isoform)
TransDecoder.LongOrfs -t  EAF3_transcripts_5prime_3kbextension.fa -m 50

# predict coding regions, keep all ORFs per transcript (not just best)
TransDecoder.Predict -t  EAF3_transcripts_5prime_3kbextension.fa --no_refine_starts # which skips the  Position Weight Matrix model training step, this is okay for this situation because I only care about one gene probably not okay in de novo genome assembly

# copy final protein sequences to a clean output file
p EAF3_transcripts_5prime_3kbextension.fa.transdecoder.pep EAF3_proteins_5prime_3kbextension.fa

echo "Done — protein sequences written to $OUTDIR/EAF3_proteins_5prime_3kbextension.fa"
echo ""
echo "Predicted proteins:"
grep ">" $OUTDIR/EAF3_proteins_5prime_3kbextension.fa

echo ""
echo "Checking for remaining partial ORFs:"
grep "5prime_partial" $OUTDIR/EAF3_proteins_5prime_3kbextension.fa | wc -l
echo "partial ORFs remaining (should be 0 if upstream extension worked)"

echo ""
echo "Protein lengths:"
grep ">" $OUTDIR/EAF3_proteins_5prime_3kbextension.fa | grep -oP "len:\d+" | sort -t: -k2 -rn
