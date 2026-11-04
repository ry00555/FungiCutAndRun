#!/bin/bash
# =============================================================================
# FiberTools_config.sh — the ONLY file you need to edit between runs
# =============================================================================
# Every FiberTools script sources this file. To process a different run:
#   - change RUN below, or
#   - leave it and pass the run at submit time:
#       sbatch --export=ALL,RUN=ONTRun14 FiberTools_Pt1_MergeBams.sh
#
# Expected layout for each run:
#   /scratch/ry00555/<RUN>/bam_pass/barcode01/*.bam ...   (from MinKNOW)
#   /scratch/ry00555/<RUN>/<RUN>_meta.csv                 (sample sheet)
# Everything is written to /scratch/ry00555/<RUN>/FiberTools/
# =============================================================================

# ── Run ──────────────────────────────────────────────────────────────────────
RUN="${RUN:-ONTRun13}"

BASE="/scratch/ry00555/${RUN}"
BAM_PASS="${BASE}/bam_pass"
META="${BASE}/${RUN}_meta.csv"

# ── Output folders (no need to edit) ─────────────────────────────────────────
OUTROOT="${BASE}/FiberTools"
MERGED_DIR="${OUTROOT}/MergedBams"            # Pt1 output
RESULTS_DIR="${OUTROOT}/fibertools_results"   # Pt2 + Pt3 output
HEATMAP_DIR="${OUTROOT}/MergedForHeatmaps"    # MergeForHeatmaps + MergedQC
SAMPLE_TABLE="${OUTROOT}/${RUN}_samples.tsv"  # made by Pt1, read by Pt2/Pt3

# ── Filtering ────────────────────────────────────────────────────────────────
MIN_MAPQ=10

# ── Genomes ──────────────────────────────────────────────────────────────────
# Short name -> genome FASTA. Add a line here to use a new genome.
declare -A GENOME_FASTA=(
  [Ncrassa]="/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic.fna"
  [Afum_CEA10]="/home/ry00555/Research/Genomes/Afum_CEA10/Afum_A1163_edit6.fasta"
)

# Optional: an existing chrom-sizes file (chrom<TAB>length) per genome.
# Genomes not listed get one made automatically from their FASTA.
declare -A GENOME_SIZES=(
  [Ncrassa]="/home/ry00555/Research/Genomes/GenBankNcrassachromsizes.txt"
)

# Genome used when no rule below matches a sample's Strain
DEFAULT_GENOME="Ncrassa"

# Strain -> genome rules: "text in Strain column : genome short name".
# Checked in order; the first match wins.
GENOME_RULES=(
  "CEA10:Afum_CEA10"
  "_Af:Afum_CEA10"
)

# =============================================================================
# Helper functions (no need to edit)
# =============================================================================

# genome_for_strain <Strain>  -> prints the genome short name
genome_for_strain () {
  local strain="$1" rule
  for rule in "${GENOME_RULES[@]}"; do
    if [[ "$strain" == *"${rule%%:*}"* ]]; then echo "${rule#*:}"; return; fi
  done
  echo "$DEFAULT_GENOME"
}

# sizes_for_genome <genome>  -> prints a chrom-sizes file (made from the FASTA if needed)
sizes_for_genome () {
  local g="$1"
  if [[ -n "${GENOME_SIZES[$g]:-}" && -s "${GENOME_SIZES[$g]}" ]]; then
    echo "${GENOME_SIZES[$g]}"; return
  fi
  local fasta="${GENOME_FASTA[$g]:?no FASTA set for genome $g in FiberTools_config.sh}"
  local out="${OUTROOT}/genome_sizes/${g}.sizes"
  if [[ ! -s "$out" ]]; then
    mkdir -p "$(dirname "$out")"
    ( [[ "$fasta" == *.gz ]] && zcat "$fasta" || cat "$fasta" ) | awk '
      /^>/ { if (n) print n "\t" l; n = substr($1, 2); l = 0; next }
           { l += length($0) }
      END  { if (n) print n "\t" l }' > "$out"
  fi
  echo "$out"
}

# make_sample_table  -> writes $SAMPLE_TABLE from $META
# Columns: barcode, sample, strain, genome, fasta, sizes
# Rows with an empty Strain (unused barcodes) are skipped.
make_sample_table () {
  [[ -s "$META" ]] || { echo "ERROR: metadata not found: $META"; return 1; }
  mkdir -p "$OUTROOT"
  printf "barcode\tsample\tstrain\tgenome\tfasta\tsizes\n" > "$SAMPLE_TABLE"
  local barcode owner strain sample genome
  while IFS=$'\t' read -r barcode owner strain; do
    [[ -z "$strain" || -z "$barcode" ]] && continue
    sample="${RUN}_${strain}${owner:+_$owner}_${barcode}"
    genome=$(genome_for_strain "$strain")
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$barcode" "$sample" "$strain" "$genome" \
      "${GENOME_FASTA[$genome]}" "$(sizes_for_genome "$genome")" >> "$SAMPLE_TABLE"
  done < <(tr -d '\r' < "$META" | awk -F',' '
      NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
      { gsub(/ /, "_", $col["Strain"]); gsub(/ /, "_", $col["Owner"])
        print $col["Index"] "\t" $col["Owner"] "\t" $col["Strain"] }')
}

# table_value <sample> <column name>  -> one value from $SAMPLE_TABLE
table_value () {
  awk -F'\t' -v s="$1" -v c="$2" '
    NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
    $2 == s { print $col[c]; exit }' "$SAMPLE_TABLE"
}
