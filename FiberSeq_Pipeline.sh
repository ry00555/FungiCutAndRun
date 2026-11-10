#!/bin/bash
#SBATCH --job-name=FiberSeq
#SBATCH --partition=batch
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=24:00:00

# =============================================================================
#  FiberSeq_pipeline.sh — complete Fiber-seq pipeline in one script
# =============================================================================
#  Run with bash (NOT sbatch) from the folder holding this script and
#  fiber_nuc_spacing.py:
#
#     bash FiberSeq_pipeline.sh                          # submit everything
#     RUN_LIST="ONTRun14" bash FiberSeq_pipeline.sh      # a different run, no editing
#     POOL_RUNS=true bash FiberSeq_pipeline.sh           # pool the same DNA sequenced in several runs
#     MODE=resummary bash FiberSeq_pipeline.sh           # re-submit ONLY the summary job (TSS plots etc.)
#
#  It reads each run's metadata CSV, then submits
#    1. one job PER SAMPLE (all in parallel), each doing:
#         merge split BAMs (+ re-align to the chosen genome if needed, keeping m6A/5mC)
#         -> ft add-nucleosomes -> ft fire (FIRE elements) -> ft pileup -> BigWigs
#         -> nucleosome + accessible-patch (MSP) positions, fiber/footprint lengths
#         -> mononucleosome dyad track, dinucleosome footprints, NDR track + NDR calls
#         -> FIRE peaks -> IGV decorators -> ft qc -> spacing / NRL
#         -> single-fiber tables for PLOT_REGIONS (for FiberSeq_figures.Rmd)
#    2. one SUMMARY job after all samples finish:
#         combined QC + spacing tables, comparison plots, strain-averaged BigWigs,
#         TSS profiles + heatmaps (nucleosome occupancy and 6mA)
#
#  Every step skips work that already exists, so re-running after a failure
#  only redoes what is missing.
# =============================================================================


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  SETTINGS — the only part you should need to edit                         ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

# Runs to process together (space-separated). Override without editing:
#   RUN_LIST="ONTRun14" bash FiberSeq_pipeline.sh
RUN_LIST="${RUN_LIST:-ONT13_CEA10align_megabam}"

# Where each run's MinKNOW BAMs and metadata live ({RUN} is replaced by the run name)
BAM_PASS_PATTERN="/scratch/ry00555/Afum/{RUN}/bam_pass"          # contains barcode01/, barcode02/ ...
META_PATTERN="/scratch/ry00555/Afum/{RUN}/{RUN}_meta.csv"        # columns: run,barcode,sample_name,genome,genotype,condition,bio_rep,tech_rep,include,...

# ── GENOME: un-comment ONE block, comment out the others ─────────────────────
# GENOME_SIZES and TSS_BED are optional ("" = make sizes from the FASTA / skip TSS plots).
# All samples are aligned to this genome (re-aligned if MinKNOW used a different one).

#GENOME_NAME="Ncrassa"
#GENOME_FASTA="/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic.fna"
#GENOME_SIZES="/home/ry00555/Research/Genomes/GenBankNcrassachromsizes.txt"
#TSS_BED="/home/ry00555/Research/Genomes/neurospora_TSS_stranded.bed"

GENOME_NAME="Afum_CEA10"
GENOME_FASTA="/home/ry00555/Research/Genomes/Afum_CEA10/GCA_051225625.1_ASM5122562v1_genomic.fna"
GENOME_SIZES="/home/ry00555/Research/Genomes/Afum_CEA10/GCA_051225625.1.chrom.sizes.txt"
TSS_BED="/home/ry00555/Research/Genomes/Afum_CEA10/GCA_051225625.1_TSS_CEA10.bed"


#make a new .bed file from a list of genes in .txt
#grep -Ff my_gene_list.txt reference.gtf > extracted_genes.gtf
#awk '$3=="gene" {print $1 "\t" $4-1 "\t" $5 "\t" $9}' extracted_genes.gtf > extracted_genes.bed


#module load SAMtools

#cd /scratch/ry00555/ONTRun13/bam_pass          # change to your run

# list every barcode's split BAMs (unclassified is left out)
# ls barcode0[7-9]/*.bam > all_bams.txt
# wc -l < all_bams.txt                            # how many files will be merged
#
# # concatenate, sort by position, index
# samtools cat -b all_bams.txt -o - \
#   | samtools sort -@ 8 -m 2G -o ONTRun13_all_barcodes.bam -
# samtools index ONTRun13_all_barcodes.bam
# samtools flagstat ONTRun13_all_barcodes.bam

# Pool re-sequenced DNA: true = rows with the same sample_name in different runs are merged
# and analysed as ONE sample (more fibers). false = every run/barcode analysed separately.
# Pooled and unpooled results go to different folders, so you can run both and compare.
POOL_RUNS="${POOL_RUNS:-false}"

# All output goes here: one folder per set of runs + genome (+ "_pooled")
PROJECT="${PROJECT:-$(echo "$RUN_LIST" | tr ' ' '_')}"
OUT="/scratch/ry00555/FiberSeq/${PROJECT}_${GENOME_NAME}$([[ "$POOL_RUNS" == true ]] && echo _pooled)"

# Analysis settings
MIN_MAPQ=10        # alignment quality filter when merging
FIRE_OPTS="--ont"  # Nanopore reads
MIN_COV=5          # min fibers per bp for occupancy / FIRE / methylation tracks
MIN_NUC_LEN=50     # ignore nucleosome calls shorter than this (bp)
MAX_LAG=1000       # max dyad distance for the spacing ladder (bp)
KEEP_MERGED_BAM=false   # true = keep the pre-nucleosome merged BAM too (fire.bam has everything)

# Nucleosome-depleted regions (NDRs): accessible patches (MSPs) at least this long
NDR_MIN_MSP=150    # bp; an MSP this long has no nucleosome over it
NDR_MIN_FRAC=0.3   # call an NDR where >= this fraction of fibers have such an MSP
NDR_MIN_LEN=100    # bp; minimum called NDR width
# Dinucleosomes: protected footprints in this size range (two nucleosomes, no resolvable linker)
DINUC_MIN=230
DINUC_MAX=400
# FIRE peaks: regions where >= this fraction of fibers carry a FIRE element
FIRE_PEAK_FRAC=0.2

# Optional BED of regions for single-fiber plots (chrom start end name), "" = none.
# Keep regions small (≤ ~20 kb). Used by FiberSeq_figures.Rmd.
PLOT_REGIONS=""

# Software
CONDA_ENV="/home/ry00555/fibertools"
load_fibertools () { module load Miniforge3/24.11.3-0; source activate "$CONDA_ENV"; }
MINIMAP2_MODULE="minimap2"
UCSC_MODULE="ucsc/443"
DEEPTOOLS_MODULE="deepTools/3.5.5-gfbf-2023a"
BEDTOOLS_MODULE="BEDTools"

# ═════════════════════════════════════════════════════════════════════════════
#  Nothing below needs editing
# ═════════════════════════════════════════════════════════════════════════════

set -euo pipefail
SAMPLE_TABLE="$OUT/samples.tsv"
D_BAM="$OUT/bams"; D_BW="$OUT/bigwigs"; D_BED="$OUT/beds"; D_STAT="$OUT/stats"
D_SUM="$OUT/summary"; D_LOG="$OUT/logs"; D_GEN="$OUT/genome_sizes"

# chrom sizes for the chosen genome (given file, or made once from the FASTA)
genome_sizes () {
  if [[ -n "$GENOME_SIZES" && -s "$GENOME_SIZES" ]]; then echo "$GENOME_SIZES"; return; fi
  local out="$D_GEN/$GENOME_NAME.sizes"
  if [[ ! -s "$out" ]]; then
    mkdir -p "$D_GEN"
    ( [[ "$GENOME_FASTA" == *.gz ]] && zcat "$GENOME_FASTA" || cat "$GENOME_FASTA" ) | awk '
      /^>/ { if (n) print n "\t" l; n = substr($1, 2); l = 0; next } { l += length($0) }
      END  { if (n) print n "\t" l }' > "$out"
  fi
  echo "$out"
}

# bedgraph -> bigWig (clip to chrom ends, drop unknown contigs, sort)
to_bw () {
  local bg="$1" bw="$2" sizes="$3"
  [[ -s "$bg" ]] || { rm -f "$bg"; echo "     ⚠️  no data for $(basename "$bw")"; return 0; }
  awk -v OFS='\t' 'NR==FNR{len[$1]=$2; next}
       ($1 in len){ if($2<0)$2=0; if($3>len[$1])$3=len[$1]; if($2<$3) print $1,$2,$3,$4 }' \
       "$sizes" "$bg" | LC_ALL=C sort -k1,1 -k2,2n > "$bg.sorted"
  bedGraphToBigWig "$bg.sorted" "$sizes" "$bw"
  rm -f "$bg" "$bg.sorted"
  echo "     ✅ $(basename "$bw")"
}

# =============================================================================
#  SETUP (bash on the login node): build sample table + submit jobs
# =============================================================================
setup () {
  [[ -f fiber_nuc_spacing.py ]] || { echo "ERROR: run from the folder containing fiber_nuc_spacing.py"; exit 1; }
  [[ -s "$GENOME_FASTA" ]] || { echo "ERROR: genome FASTA not found: $GENOME_FASTA"; exit 1; }
  mkdir -p "$OUT" "$D_LOG"
  local SZ; SZ=$(genome_sizes)

  printf "run\tbarcode\tsample\tstrain\trep\tgenome\tfasta\tsizes\tbam_dir\n" > "$SAMPLE_TABLE"
  local RUN META BP
  for RUN in $RUN_LIST; do
    META="${META_PATTERN//\{RUN\}/$RUN}"; BP="${BAM_PASS_PATTERN//\{RUN\}/$RUN}"
    [[ -s "$META" ]] || { echo "ERROR: metadata not found: $META"; exit 1; }
    [[ -d "$BP" ]]   || { echo "ERROR: bam_pass not found: $BP"; exit 1; }
    # Metadata columns used: barcode, sample_name, genome, genotype, condition, bio_rep, tech_rep, include
    # Only rows with include=yes AND genome = the active GENOME_NAME are processed.
    tr -d '\r' < "$META" | awk -F',' -v G="$GENOME_NAME" '
      NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i
                if (!("sample_name" in c) || !("genome" in c)) {
                  print "ERROR: metadata needs columns barcode,sample_name,genome,genotype,condition,bio_rep,tech_rep,include" > "/dev/stderr"; exit 1 }
                next }
      { inc = tolower($c["include"]); if (inc != "yes") next
        if ($c["genome"] != G) { skipped++; next }
        grp = $c["genotype"]; if ($c["condition"] != "") grp = grp "_" $c["condition"]
        rep = "rep" $c["bio_rep"] "_tech" $c["tech_rep"]
        gsub(/[^A-Za-z0-9_.-]/, "_", grp)
        print $c["barcode"] "\t" $c["sample_name"] "\t" grp "\t" rep }
      END { if (skipped) printf "  (%d sample(s) in this run use another genome — skipped)\n", skipped > "/dev/stderr" }' \
    | while IFS=$'\t' read -r BC NAME GROUP REP; do
        local S="${RUN}_${NAME}"
        if cut -f3 "$SAMPLE_TABLE" | grep -qx "$S"; then
          echo "ERROR: sample_name '$NAME' appears twice in $META — names must be unique"; exit 1
        fi
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$RUN" "$BC" "$S" "$GROUP" "$REP" "$GENOME_NAME" \
          "$GENOME_FASTA" "$SZ" "$BP/$BC" >> "$SAMPLE_TABLE"
      done
  done

  # pool rows with the same sample_name across runs into one sample
  if [[ "$POOL_RUNS" == true ]]; then
    awk -F'\t' -v OFS='\t' '
      NR == 1 { print; next }
      { name = substr($3, length($1) + 2)                       # sample_name without the run prefix
        if (!(name in seen)) { seen[name] = 1; order[++n] = name; grp[name] = $4; rep[name] = $5
                               gen[name] = $6; fa[name] = $7; sz[name] = $8 }
        runs[name] = runs[name] (runs[name] ? "+" : "") $1
        bcs[name]  = bcs[name]  (bcs[name]  ? "+" : "") $1 ":" $2
        dirs[name] = dirs[name] (dirs[name] ? "," : "") $9 }
      END { for (i = 1; i <= n; i++) { k = order[i]
              print runs[k], bcs[k], k, grp[k], rep[k], gen[k], fa[k], sz[k], dirs[k] } }' \
      "$SAMPLE_TABLE" > "$SAMPLE_TABLE.pooled" && mv "$SAMPLE_TABLE.pooled" "$SAMPLE_TABLE"
  fi

  local N; N=$(($(wc -l < "$SAMPLE_TABLE") - 1))
  echo "============================================================"
  echo " FiberSeq pipeline — $PROJECT ($N samples, genome: $GENOME_NAME, pooled runs: $POOL_RUNS)"
  echo " Output: $OUT"
  echo "============================================================"
  awk -F'\t' 'NR>1{ t = 0; k = split($9, d, ",")
                    for (j = 1; j <= k; j++) { n = 0; cmd = "ls " d[j] "/*.bam 2>/dev/null | wc -l"
                                               cmd | getline n; close(cmd); t += n }
                    printf "  %-40s %-30s %s split BAMs\n", $3, $2, t }' "$SAMPLE_TABLE"
  (( N > 0 )) || { echo "ERROR: no samples with include=yes and genome=$GENOME_NAME in the metadata"; exit 1; }

  export RUN_LIST PROJECT POOL_RUNS
  local J1 J2
  J1=$(sbatch --parsable --array=1-"$N" --job-name="FS_${PROJECT}" \
         --output="$D_LOG/sample_%a_%A.out" --error="$D_LOG/sample_%a_%A.err" \
         --export=ALL,MODE=sample "$0")
  J2=$(sbatch --parsable --dependency=afterany:"$J1" --job-name="FS_${PROJECT}_summary" \
         --cpus-per-task=8 --mem=32G --time=8:00:00 --mail-type=END,FAIL \
         --output="$D_LOG/summary_%j.out" --error="$D_LOG/summary_%j.err" \
         --export=ALL,MODE=summary "$0")
  echo ""
  echo "Submitted: sample jobs $J1 (array of $N, in parallel); summary job $J2 (after all samples)"
  echo "Logs: $D_LOG   Progress: squeue -u \$USER"
}

# =============================================================================
#  SAMPLE (one Slurm array task)
# =============================================================================
run_sample () {
  local row; row=$(awk -F'\t' -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR==n' "$SAMPLE_TABLE")
  [[ -n "$row" ]] || { echo "No row $SLURM_ARRAY_TASK_ID in $SAMPLE_TABLE"; exit 1; }
  local RUN BC S STRAIN REP G FASTA SIZES SRC
  IFS=$'\t' read -r RUN BC S STRAIN REP G FASTA SIZES SRC <<< "$row"
  local T="${SLURM_CPUS_PER_TASK:-8}" TMP="$OUT/tmp/$S"
  mkdir -p "$D_BAM" "$D_BW" "$D_BED" "$D_STAT" "$TMP"

  local MERGED="$D_BAM/$S.merged.bam" FIRE="$D_BAM/$S.fire.bam"
  echo "============================================================"
  echo " $S   ($RUN $BC, strain $STRAIN, genome $G)"
  echo "============================================================"

  load_fibertools
  module load "$UCSC_MODULE" "$BEDTOOLS_MODULE"

  # ── 1. merge split BAMs (+ re-align if needed) ──────────────────────────────
  echo ""; echo "[1/11] merge split BAMs"
  if [[ ! -s "$FIRE" && ! -s "$MERGED" ]]; then
    # SRC is one barcode folder, or several (comma-separated) when POOL_RUNS=true
    local DIRS PARTS=() i=0 dir part
    IFS=',' read -ra DIRS <<< "$SRC"
    for dir in "${DIRS[@]}"; do
      local BAMS=( "$dir"/*.bam )
      [[ -f "${BAMS[0]}" ]] || { echo "❌ no BAMs in $dir"; exit 1; }
      part="$TMP/part$i.bam"; i=$((i + 1))
      echo "     $dir: ${#BAMS[@]} split BAMs"
      if diff -q <(samtools view -H "${BAMS[0]}" | awk '$1=="@SQ"{sub(/^SN:/,"",$2); print $2}' | sort) \
                 <(cut -f1 "$SIZES" | sort) > /dev/null; then
        echo "       already aligned to $G → MAPQ≥$MIN_MAPQ filter + sort"
        samtools cat "${BAMS[@]}" | samtools view -u -q "$MIN_MAPQ" - \
          | samtools sort -@ "$T" -T "$TMP/s$i" -o "$part" -
      else
        echo "       not aligned to $G → re-aligning with minimap2 (keeping MM/ML/MN tags)"
        module load "$MINIMAP2_MODULE"
        [[ -s "$FASTA" ]] || { echo "❌ FASTA not found: $FASTA"; exit 1; }
        samtools cat "${BAMS[@]}" | samtools fastq -F 0x900 -T MM,ML,MN - \
          | minimap2 -ax map-ont -y -t "$T" "$FASTA" - \
          | samtools view -u -q "$MIN_MAPQ" - \
          | samtools sort -@ "$T" -T "$TMP/s$i" -o "$part" -
      fi
      PARTS+=("$part")
    done
    if (( ${#PARTS[@]} == 1 )); then
      mv "${PARTS[0]}" "$MERGED"
    else
      echo "     pooling ${#PARTS[@]} runs"
      samtools merge -@ "$T" -o "$MERGED" "${PARTS[@]}"
      rm -f "${PARTS[@]}"
    fi
    samtools index "$MERGED"
    echo "     ✅ $(samtools view -c "$MERGED") reads"
  else
    echo "     ⏭️  done"
  fi

  # ── 2. nucleosomes + FIRE ───────────────────────────────────────────────────
  echo ""; echo "[2/11] ft add-nucleosomes + ft fire (FIRE elements)"
  if [[ ! -s "$FIRE" ]]; then
    ft add-nucleosomes -t "$T" "$MERGED" "$TMP/nucs.bam"
    ft fire -t "$T" $FIRE_OPTS "$TMP/nucs.bam" "$FIRE"
    samtools index "$FIRE"
    rm -f "$TMP/nucs.bam"
    [[ "$KEEP_MERGED_BAM" == true ]] || rm -f "$MERGED" "$MERGED.bai"
  else
    echo "     ⏭️  done"
  fi
  [[ -s "$D_BED/$S.FIRE.bed.gz" ]] || ft fire --extract "$FIRE" "$D_BED/$S.FIRE.bed.gz"

  # ── 3. pileup -> tracks ─────────────────────────────────────────────────────
  echo ""; echo "[3/11] ft pileup -> BigWigs"
  local PU="$D_STAT/$S.pileup.bed.gz"
  [[ -s "$PU" ]] || ft pileup -t "$T" --m6a --cpg --fiber-coverage "$FIRE" | gzip > "$PU"
  if [[ ! -s "$D_BW/$S.nuc_occupancy.bw" ]]; then
    zcat "$PU" | awk -v OFS='\t' -v mc="$MIN_COV" -v p="$TMP/$S" '
      /^#/ { for (i = 1; i <= NF; i++) { h = tolower($i); sub(/^#/, "", h); c[h] = i
               if (!m6 && h ~ /m6a/) m6 = i; if (!cg && h ~ /cpg|5mc/) cg = i }; next }
      { cv = $c["coverage"]; if (cv < mc) next
        print $1,$2,$3,cv > p ".coverage.bg"
        if ("nuc_coverage"  in c) print $1,$2,$3,$c["nuc_coverage"]/cv  > p ".nuc_occupancy.bg"
        if ("msp_coverage"  in c) print $1,$2,$3,$c["msp_coverage"]/cv  > p ".accessibility.bg"
        if ("fire_coverage" in c) print $1,$2,$3,$c["fire_coverage"]/cv > p ".FIRE_fraction.bg"
        if (m6 && $m6 > 0)        print $1,$2,$3,$m6/cv                 > p ".m6A_fraction.bg"
        if (cg && $cg > 0)        print $1,$2,$3,$cg/cv                 > p ".5mC_fraction.bg" }'
    local t
    for t in coverage nuc_occupancy accessibility FIRE_fraction m6A_fraction 5mC_fraction; do
      [[ -f "$TMP/$S.$t.bg" ]] && to_bw "$TMP/$S.$t.bg" "$D_BW/$S.$t.bw" "$SIZES"
    done
  else
    echo "     ⏭️  done"
  fi

  # ── 4. nucleosomes, MSPs (accessible patches), fiber lengths ────────────────
  echo ""; echo "[4/11] nucleosome + accessible-patch (MSP) positions, fiber lengths"
  local NUCS="$D_BED/$S.nucleosomes.bed.gz" MSPS="$D_BED/$S.msps.bed.gz"
  if [[ ! -s "$NUCS" || ! -s "$MSPS" ]]; then
    ft extract -t "$T" --reference --nuc "$TMP/nuc.bed12" --msp "$TMP/msp.bed12" "$FIRE"
    # bed12 blocks -> one line per feature: chrom start end fiber length strand fiber_start fiber_end
    # (ft adds 1-bp placeholder blocks at fiber ends; those are dropped)
    expand_bed12 () { awk -v OFS='\t' -v ml="$1" '/^#/ || NF < 12 { next }
      { n = split($11, sz, ","); split($12, st, ",")
        for (i = 1; i <= n; i++) { if (sz[i] == "" || sz[i] < ml) continue
          s = $2 + st[i]; if (s < 0) continue
          print $1, s, s + sz[i], $4, sz[i], $6, $2, $3 } }' "$2"; }
    expand_bed12 "$MIN_NUC_LEN" "$TMP/nuc.bed12" | gzip > "$NUCS"
    expand_bed12 2              "$TMP/msp.bed12" | gzip > "$MSPS"
    rm -f "$TMP/nuc.bed12" "$TMP/msp.bed12"
  fi
  # length histograms (small tables for the Rmd figures)
  hist_col () { awk -v c="$1" '!/^#/ { h[$c]++ } END { for (k in h) print k "\t" h[k] }' | sort -k1,1n; }
  [[ -s "$D_STAT/$S.nuc_length_hist.tsv" ]]  || { echo -e "length\tcount"; zcat "$NUCS" | hist_col 5; } > "$D_STAT/$S.nuc_length_hist.tsv"
  [[ -s "$D_STAT/$S.msp_length_hist.tsv" ]]  || { echo -e "length\tcount"; zcat "$MSPS" | hist_col 5; } > "$D_STAT/$S.msp_length_hist.tsv"
  [[ -s "$D_STAT/$S.fire_length_hist.tsv" ]] || { echo -e "length\tcount"; zcat -f "$D_BED/$S.FIRE.bed.gz" \
      | awk '!/^#/ { print $3 - $2 }' | hist_col 1; } > "$D_STAT/$S.fire_length_hist.tsv"
  [[ -s "$D_STAT/$S.fiber_lengths.tsv.gz" ]] || samtools view -F 0x904 "$FIRE" \
      | awk 'BEGIN { print "fiber\tlength" } { print $1 "\t" length($10) }' | gzip > "$D_STAT/$S.fiber_lengths.tsv.gz"
  echo "     $(zcat "$NUCS" | wc -l) nucleosomes, $(zcat "$MSPS" | wc -l) MSPs, $(($(zcat "$D_STAT/$S.fiber_lengths.tsv.gz" | wc -l) - 1)) fibers"

  # coverage bedgraph (fibers per bp), needed for the fraction tracks below
  local COVBG="$TMP/$S.cov.bg"
  zcat "$PU" | awk -v OFS='\t' '/^#/ { for (i = 1; i <= NF; i++) { h = tolower($i); sub(/^#/, "", h); c[h] = i }; next }
                                { print $1, $2, $3, $c["coverage"] }' | LC_ALL=C sort -k1,1 -k2,2n > "$COVBG"

  # interval_cov <bed.gz> <filter-awk> <out.bg> : per-bp count of features passing the filter
  interval_cov () {
    zcat "$1" | awk -v OFS='\t' "$2" | awk -v OFS='\t' 'NR==FNR { ok[$1] = 1; next } ($1 in ok)' "$SIZES" - \
      | LC_ALL=C sort -k1,1 -k2,2n | bedtools genomecov -bg -i - -g "$SIZES" > "$3"
  }
  # fraction_bg <count.bg> <out.bg> : count / fiber coverage, where coverage >= MIN_COV
  fraction_bg () {
    bedtools unionbedg -i "$1" "$COVBG" \
      | awk -v OFS='\t' -v mc="$MIN_COV" '$5 >= mc { print $1, $2, $3, $4 / $5 }' > "$2"
  }
  # call_regions <fraction.bg> <min fraction> <min length> <out.bed> : merge high-signal bp into regions
  call_regions () {
    awk -v OFS='\t' -v f="$2" '$4 >= f' "$1" | bedtools merge -d 50 -c 4 -o max -i - \
      | awk -v OFS='\t' -v ml="$3" '$3 - $2 >= ml { print $1, $2, $3, "region" NR, $4 }' > "$4"
  }

  # ── 5. mononucleosome footprint (dyad) track ────────────────────────────────
  echo ""; echo "[5/11] mononucleosome dyad / footprint track"
  if [[ ! -s "$D_BW/$S.nuc_dyads.bw" ]]; then
    local tot; tot=$(zcat "$NUCS" | awk -v d="$DINUC_MIN" '$5 < d' | wc -l)
    zcat "$NUCS" | awk -v OFS='\t' -v d="$DINUC_MIN" '$5 < d { print $1, int(($2 + $3) / 2) }' \
      | LC_ALL=C sort -k1,1 -k2,2n | uniq -c \
      | awk -v OFS='\t' -v t="$tot" '{ print $2, $3, $3 + 1, $1 * 1e6 / t }' > "$TMP/$S.dyads.bg"
    to_bw "$TMP/$S.dyads.bg" "$D_BW/$S.nuc_dyads.bw" "$SIZES"
  fi

  # ── 6. dinucleosome footprints ──────────────────────────────────────────────
  echo ""; echo "[6/11] dinucleosome footprints (${DINUC_MIN}-${DINUC_MAX} bp protected)"
  if [[ ! -s "$D_BED/$S.dinucleosomes.bed.gz" ]]; then
    zcat "$NUCS" | awk -v a="$DINUC_MIN" -v b="$DINUC_MAX" '$5 >= a && $5 <= b' | gzip > "$D_BED/$S.dinucleosomes.bed.gz"
  fi
  echo "     $(zcat "$D_BED/$S.dinucleosomes.bed.gz" | wc -l) dinucleosome footprints"
  if [[ ! -s "$D_BW/$S.dinuc_fraction.bw" ]]; then
    interval_cov "$D_BED/$S.dinucleosomes.bed.gz" '{ print $1, $2, $3 }' "$TMP/$S.dinuc.count.bg"
    fraction_bg "$TMP/$S.dinuc.count.bg" "$TMP/$S.dinuc_fraction.bg"
    to_bw "$TMP/$S.dinuc_fraction.bg" "$D_BW/$S.dinuc_fraction.bw" "$SIZES"
  fi

  # ── 7. nucleosome-depleted regions (NDRs) ───────────────────────────────────
  echo ""; echo "[7/11] nucleosome-depleted regions (MSPs >= ${NDR_MIN_MSP} bp)"
  if [[ ! -s "$D_BED/$S.NDRs.bed" ]]; then
    interval_cov "$MSPS" "\$5 >= $NDR_MIN_MSP { print \$1, \$2, \$3 }" "$TMP/$S.ndr.count.bg"
    fraction_bg "$TMP/$S.ndr.count.bg" "$TMP/$S.NDR_fraction.bg"
    call_regions "$TMP/$S.NDR_fraction.bg" "$NDR_MIN_FRAC" "$NDR_MIN_LEN" "$D_BED/$S.NDRs.bed"
    to_bw "$TMP/$S.NDR_fraction.bg" "$D_BW/$S.NDR_fraction.bw" "$SIZES"
  fi
  echo "     $(wc -l < "$D_BED/$S.NDRs.bed") NDRs called"

  # ── 8. FIRE peaks ───────────────────────────────────────────────────────────
  echo ""; echo "[8/11] FIRE peaks (FIRE in >= ${FIRE_PEAK_FRAC} of fibers)"
  if [[ ! -s "$D_BED/$S.FIRE_peaks.bed" ]]; then
    zcat "$PU" | awk -v OFS='\t' -v mc="$MIN_COV" '
      /^#/ { for (i = 1; i <= NF; i++) { h = tolower($i); sub(/^#/, "", h); c[h] = i }; next }
      ("fire_coverage" in c) && $c["coverage"] >= mc { print $1, $2, $3, $c["fire_coverage"] / $c["coverage"] }' \
      | LC_ALL=C sort -k1,1 -k2,2n > "$TMP/$S.fire_frac.bg"
    call_regions "$TMP/$S.fire_frac.bg" "$FIRE_PEAK_FRAC" 50 "$D_BED/$S.FIRE_peaks.bed"
  fi
  echo "     $(wc -l < "$D_BED/$S.FIRE_peaks.bed") FIRE peaks"

  # ── 9. IGV per-fiber decorators ─────────────────────────────────────────────
  echo ""; echo "[9/11] IGV decorators"
  if [[ ! -s "$D_BED/$S.fiber_tracks.bed" ]]; then
    ft track-decorators --bed12 "$D_BED/$S.fiber_tracks.bed" \
      --decorator "$D_BED/$S.fiber_decorations.bed" "$FIRE" || echo "     ⚠️  track-decorators failed"
  fi

  # ── 10. QC + spacing ────────────────────────────────────────────────────────
  echo ""; echo "[10/11] ft qc (+ m6A autocorrelation) and nucleosome spacing"
  [[ -s "$D_STAT/$S.qc.txt" ]] || ft qc -t "$T" --acf --m6a-per-msp "$FIRE" "$D_STAT/$S.qc.txt"
  [[ -s "$D_STAT/$S.spacing_summary.tsv" ]] || \
    python "$SLURM_SUBMIT_DIR/fiber_nuc_spacing.py" sample \
      --nucs "$NUCS" --sample "$S" --out "$D_STAT" --max-lag "$MAX_LAG" --max-nuc-len "$((DINUC_MIN - 1))"

  # ── 11. single-fiber data for plotting regions ──────────────────────────────
  echo ""; echo "[11/11] single-fiber data for plot regions"
  if [[ -n "$PLOT_REGIONS" && -s "$PLOT_REGIONS" ]]; then
    mkdir -p "$OUT/regions"
    tr -d '\r' < "$PLOT_REGIONS" | awk -v OFS='\t' '!/^#/ && NF >= 3 { print $1, $2, $3, ($4 != "" ? $4 : $1 "_" $2 "_" $3) }' > "$OUT/regions/regions.bed"
    while IFS=$'\t' read -r c s e name _; do
      [[ -z "$c" || "$c" == \#* ]] && continue
      name="${name:-${c}_${s}_${e}}"
      local RF="$OUT/regions/$S.$name.fibers.tsv.gz"
      [[ -s "$RF" ]] && continue
      samtools view -b "$FIRE" "$c:$((s + 1))-$e" > "$TMP/region.bam"
      samtools index "$TMP/region.bam"
      ft extract -t "$T" --simplify --all "$TMP/region.tsv" "$TMP/region.bam"
      gzip -c "$TMP/region.tsv" > "$RF"
      echo "     $name: $(($(zcat "$RF" | wc -l) - 1)) fibers"
    done < <(tr -d '\r' < "$PLOT_REGIONS")
  else
    echo "     no PLOT_REGIONS set — skipped"
  fi

  rm -rf "$TMP"
  echo ""; echo "✅ $S complete"
}

# =============================================================================
#  SUMMARY (after all samples)
# =============================================================================
run_summary () {
  mkdir -p "$D_SUM"
  load_fibertools

  echo "== Spacing comparison"
  python "$SLURM_SUBMIT_DIR/fiber_nuc_spacing.py" compare --in "$D_STAT" --out "$D_SUM" \
    || echo "⚠️  spacing comparison failed"

  echo "== Combined QC table"
  { echo -e "sample\tstatistic\tvalue\tcount"
    for q in "$D_STAT"/*.qc.txt; do
      [[ -f "$q" ]] || continue
      grep -v "fiber_length\|msp_length\|nuc_length\|m6a_per_msp_size\|acf\|^statistic" "$q" \
        | awk -v s="$(basename "$q" .qc.txt)" '{ print s "\t" $0 }'
    done; } > "$D_SUM/all_samples.qc.tsv"
  echo "   $D_SUM/all_samples.qc.tsv"
  conda deactivate 2>/dev/null || source deactivate 2>/dev/null || true

  module load "$DEEPTOOLS_MODULE"
  local T="${SLURM_CPUS_PER_TASK:-8}" G STRAIN TR

  echo "== Group-averaged BigWigs (genotype + condition; replicates and runs combined)"
  mkdir -p "$D_SUM/strain_bigwigs"
  for STRAIN in $(tail -n +2 "$SAMPLE_TABLE" | cut -f4 | sort -u); do
    for TR in nuc_occupancy nuc_dyads dinuc_fraction NDR_fraction accessibility FIRE_fraction m6A_fraction; do
      local BWS=()
      while IFS=$'\t' read -r s; do
        [[ -s "$D_BW/$s.$TR.bw" ]] && BWS+=("$D_BW/$s.$TR.bw")
      done < <(awk -F'\t' -v st="$STRAIN" 'NR>1 && $4==st {print $3}' "$SAMPLE_TABLE")
      (( ${#BWS[@]} )) || continue
      local O="$D_SUM/strain_bigwigs/$STRAIN.$TR.bw"
      [[ -s "$O" ]] || bigwigAverage -b "${BWS[@]}" -bs 10 -p "$T" -o "$O" \
        && echo "   $STRAIN.$TR.bw (${#BWS[@]} samples)"
    done
  done

  echo "== TSS profiles"
  G="$GENOME_NAME"; local BED="$TSS_BED"
  if [[ -n "$BED" && -s "$BED" ]]; then
    local STRAINS; STRAINS=$(tail -n +2 "$SAMPLE_TABLE" | cut -f4 | sort -u)
    for TR in nuc_occupancy m6A_fraction nuc_dyads NDR_fraction accessibility FIRE_fraction; do
      local BWS=() LAB=()
      for STRAIN in $STRAINS; do
        [[ -s "$D_SUM/strain_bigwigs/$STRAIN.$TR.bw" ]] && { BWS+=("$D_SUM/strain_bigwigs/$STRAIN.$TR.bw"); LAB+=("$STRAIN"); }
      done
      (( ${#BWS[@]} )) || continue
      local M="$D_SUM/TSS_${G}.${TR}.matrix.gz"
      computeMatrix reference-point --referencePoint TSS -R "$BED" -S "${BWS[@]}" \
        --samplesLabel "${LAB[@]}" -b 1000 -a 1500 -bs 10 --missingDataAsZero --skipZeros \
        -p "$T" -o "$M" || continue
      plotProfile -m "$M" --perGroup --plotTitle "$G: $TR around TSS" -out "$D_SUM/TSS_${G}.${TR}.png"
      plotProfile -m "$M" --perGroup --plotTitle "$G: $TR around TSS" -out "$D_SUM/TSS_${G}.${TR}.pdf"
      # heatmaps: nucleosome occupancy (blue) and m6A / 6mA accessibility (purple)
      case "$TR" in
        nuc_occupancy) plotHeatmap -m "$M" --colorMap Blues   --heatmapHeight 15 \
                         --zMin 0 --yAxisLabel "nucleosome occupancy" \
                         -out "$D_SUM/TSS_${G}.${TR}.heatmap.png" ;;
        m6A_fraction)  plotHeatmap -m "$M" --colorMap Purples --heatmapHeight 15 \
                         --zMin 0 --yAxisLabel "6mA (m6A) fraction" \
                         -out "$D_SUM/TSS_${G}.${TR}.heatmap.png" ;;
      esac
      echo "   TSS_${G}.${TR}"
    done
  else
    echo "   no TSS_BED set — skipped"
  fi
  echo ""; echo "✅ Summary complete → $D_SUM"
}

# =============================================================================
case "${MODE:-setup}" in
  setup)   setup ;;
  sample)  run_sample ;;
  summary) set +e; run_summary ;;
  resummary)                                  # submit just the summary job (samples already done)
    mkdir -p "$D_LOG"; export RUN_LIST PROJECT POOL_RUNS
    sbatch --job-name="FS_${PROJECT}_summary" --cpus-per-task=8 --mem=32G --time=8:00:00 \
      --mail-type=END,FAIL --output="$D_LOG/summary_%j.out" --error="$D_LOG/summary_%j.err" \
      --export=ALL,MODE=summary "$0" ;;
  *) echo "Unknown MODE=$MODE"; exit 1 ;;
esac
