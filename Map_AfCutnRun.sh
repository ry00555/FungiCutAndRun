#!/bin/bash
#SBATCH --partition=batch
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40gb
#SBATCH --time=12:00:00

## CUT&RUN mapping for one sample, E. coli DNA spike-in (no K-MET).
## Submitted by:  bash submitCutRunJob.sh sampleList_CutRun.txt
##
## Steps: trim -> Bowtie2 to A. fumigatus + E. coli -> MAPQ filter -> remove duplicates
##        -> bigWigs (BPM and E. coli spike-in scaled) -> MACS3 peaks (no control)
##
## Output (under $outdir/CutRun/):
##   TrimmedFastQs/<name>/        trimmed reads + FastQC
##   bamFiles/<name>.bam          A. fumigatus, MAPQ-filtered, deduplicated
##   bamFiles/<name>_ecoli.bam    E. coli spike-in reads
##   bigWig/<name>.BPM.bw         normalized to sequencing depth
##   bigWig/<name>.EcoliSpikeIn.bw  scaled by E. coli spike-in
##   Peaks/<name>/                MACS3 broad peaks
##   stats/<name>.stats.txt       read counts and spike-in scale factor

cd $SLURM_SUBMIT_DIR
source Af_config.txt
set -euo pipefail

################### input files ###################
if [[ -s "${fastqPath}/${accession}_R1.fastq.gz" ]]; then sfx=".fastq.gz"; else sfx="_001.fastq.gz"; fi
read1="${fastqPath}/${accession}_R1${sfx}"
read2="${fastqPath}/${accession}_R2${sfx}"
for f in "$read1" "$read2"; do
  [[ -s "$f" ]] || { echo "ERROR: missing $f"; exit 1; }
done

################### output folders ###################
base="${outdir}/CutNRun"
trimmed="${base}/TrimmedFastQs/${accession}"
BAMDIR="${base}/SortedBamFiles"
BWDIR="${base}/BigWigs"
PEAKDIR="${base}/MACS3Peaks/${accession}"
statsDir="${base}/stats"
tmp="${base}/tempFile/${accession}"
mkdir -p "$trimmed" "$BAMDIR" "$BWDIR" "$PEAKDIR" "$statsDir" "$tmp"

bam="${BAMDIR}/${accession}.bam"
ecbam="${BAMDIR}/${accession}_ecoli.bam"
stats="${statsDir}/${accession}.stats.txt"

################### trimming ###################
module load Trim_Galore/0.6.10-GCCcore-12.3.0
trim_galore --illumina --fastqc --paired --length 20 --quality 20 --cores 2 \
  --basename "${accession}" --gzip -o "$trimmed" "$read1" "$read2"
t1="${trimmed}/${accession}_val_1.fq.gz"
t2="${trimmed}/${accession}_val_2.fq.gz"

################### alignment ###################
module load Bowtie2
module load SAMtools/1.21-GCC-13.3.0

# Standard CUT&RUN settings: local alignment, proper pairs only, fragments 10-700 bp
BT2_OPTS="--local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -p $THREADS"

# Keep only properly paired reads at MAPQ >= cutoff, then sort.
# fixmate -m adds the tags samtools markdup needs.
align () {   # $1 = index, $2 = output bam
  bowtie2 $BT2_OPTS -x "$1" -1 "$t1" -2 "$t2" \
    | samtools view -@ $THREADS -bh -f 2 -q $MAPQ - \
    | samtools sort -@ $THREADS -n -T "${tmp}/n" - \
    | samtools fixmate -@ $THREADS -m - - \
    | samtools sort -@ $THREADS -T "${tmp}/c" -o "${2%.bam}.sorted.bam" -
  if [[ "$REMOVE_DUPS" == "true" ]]; then
    samtools markdup -@ $THREADS -r -s "${2%.bam}.sorted.bam" "$2" 2> "${2%.bam}.markdup.txt"
    rm "${2%.bam}.sorted.bam"
  else
    mv "${2%.bam}.sorted.bam" "$2"
  fi
  samtools index "$2"
}

[[ -s "${BT2_GENOME}.1.bt2" ]] || { echo "ERROR: no Bowtie2 index at ${BT2_GENOME} (check config.txt)"; exit 1; }
align "$BT2_GENOME" "$bam"

# E. coli spike-in: only if its index has been built
if [[ -s "${BT2_ECOLI}.1.bt2" ]]; then
  align "$BT2_ECOLI" "$ecbam"
  ec_frags=$(samtools view -c -f 0x42 "$ecbam")
else
  echo "WARNING: no E. coli index at ${BT2_ECOLI}; skipping spike-in"
  ec_frags=0
fi

################### spike-in scale factor ###################
# Count fragments (= read 1 of each proper pair)
af_frags=$(samtools view -c -f 0x42 "$bam")

{
  echo -e "sample\t${accession}"
  echo -e "Afum_fragments\t${af_frags}"
  echo -e "Ecoli_fragments\t${ec_frags}"
} > "$stats"

################### bigWigs ###################
module load deepTools/3.5.5-gfbf-2023a

# Depth-normalized (BPM), comparable to the ChIP bigWigs
bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --smoothLength $SMOOTH \
  --extendReads -of bigwig -b "$bam" -o "${BWDIR}/${accession}.BPM.bw"

# E. coli spike-in scaled
if (( ec_frags > 0 )); then
  scale=$(awk -v c=$SPIKE_CONSTANT -v e=$ec_frags 'BEGIN {printf "%.6f", c/e}')
  echo -e "spikein_scale_factor\t${scale}" >> "$stats"
  bamCoverage -p $THREADS -bs $BIN --scaleFactor $scale --smoothLength $SMOOTH \
    --extendReads -of bigwig -b "$bam" -o "${BWDIR}/${accession}.EcoliSpikeIn.bw"
else
  echo -e "spikein_scale_factor\tNA (no E. coli reads; spike-in bigWig not made)" >> "$stats"
  echo "WARNING: no E. coli reads for ${accession}; skipped spike-in bigWig"
fi

################### peaks ###################
module load MACS3/3.0.1-gfbf-2023a
# No IgG control, so --nolambda (same settings as the ChIP pipeline)
macs3 callpeak -t "$bam" -f BAMPE -n "${accession}" --broad -g $GSIZE \
  --broad-cutoff 0.1 --outdir "$PEAKDIR" --min-length 800 --max-gap 500 --nolambda

rm -r "$tmp"
echo "${accession} done"
