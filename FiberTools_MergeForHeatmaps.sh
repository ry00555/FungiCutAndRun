#!/usr/bin/env bash
#SBATCH --job-name=FiberTools_MergeHeatmaps
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=85G
#SBATCH --time=24:00:00
#SBATCH --output=FiberTools_MergeHeatmaps.%j.out
#SBATCH --error=FiberTools_MergeHeatmaps.%j.err

# =============================================================================
# Merge nucs.bam files across replicates + both runs, then generate
# normalized TSS heatmaps and multi-sample comparison plots
#
# Key fixes vs previous version:
#   - ft pileup WITHOUT --fiber-coverage → outputs fraction of reads
#     methylated/nucleosome-occupied at each base (0-1), not raw depth
#   - Stranded TSS BED (neurospora_TSS_stranded.bed) for correct anchoring
#   - --missingDataAsZero and symmetric 2kb window
#   - Sorted by m6A signal, plotProfile for average line plots
#
# Groups defined in fiber_groups.txt (must be in same dir as this script)
# Format: GROUP_NAME barcode01 barcode02 ...
# =============================================================================

set -euo pipefail

module load Miniforge3/24.11.3-0 ucsc/443
source activate /home/ry00555/fibertools

FT_RESULTS="/lustre2/scratch/ry00555/ONTRun9_10Combined/fibertools_results"
OUT_DIR="/lustre2/scratch/ry00555/ONTRun9_10Combined/MergedForHeatmaps"
GENOME="/home/ry00555/Research/Genomes/GenBankNcrassachromsizes.txt"
TSS_BED="/home/ry00555/Research/Genomes/neurospora_TSS_stranded.bed"
GROUPS_FILE="/home/ry00555/Research/FungiCutAndRun/fiber_groups.txt"

mkdir -p "$OUT_DIR"

if [ ! -f "$GROUPS_FILE" ]; then
    echo "❌  Groups file not found: $GROUPS_FILE"
    exit 1
fi

if [ ! -f "$TSS_BED" ]; then
    echo "❌  TSS BED not found: $TSS_BED"
    echo "    Generate it with:"
    echo "    awk '!/^#/ && \$3==\"gene\" && /gene_biotype=protein_coding/ && \$1~/^CM|^KC|^KI/' \\"
    echo "        GCA_000182925.2_NC12_genomic_WithExtras.gff | \\"
    echo "    awk 'BEGIN{OFS=\"\t\"} { match(\$9,/Name=([^;]+)/,n); name=(n[1]!=\"\")?n[1]:\".\"; if (\$7==\"+\") print \$1,\$4-1,\$4,name,\".\",\$7,\$4-1,\$5,\$4-1,\$5; else print \$1,\$5-1,\$5,name,\".\",\$7,\$4-1,\$5,\$4-1,\$5 }' \\"
    echo "        | sort -k1,1 -k2,2n > neurospora_TSS_stranded.bed"
    exit 1
fi

# =============================================================================
# Helper: bedgraph → BigWig
# =============================================================================
make_bigwig() {
    local BG="$1"
    local BW="$2"

    [ -f "$BW" ] && { echo "      ⏭️  $(basename $BW) exists"; return 0; }
    [ -f "$BG" ] || { echo "      ⚠️  $(basename $BG) not found — skipping"; return 0; }

    local BG4="${BG%.bedgraph}.4col.bedgraph"
    sort -k1,1 -k2,2n "$BG" \
        | cut -f1-4 \
        | awk -v genome="$GENOME" '
            BEGIN { while ((getline line < genome) > 0) {
                split(line,a,"\t"); chromlen[a[1]]=a[2] } }
            { if ($2<0) $2=0
              if (chromlen[$1] && $3>chromlen[$1]) $3=chromlen[$1]
              if ($2<$3) print $1"\t"$2"\t"$3"\t"$4 }
        ' > "$BG4"

    bedGraphToBigWig "$BG4" "$GENOME" "$BW" && rm -f "$BG4" \
        || { echo "      ❌ bedGraphToBigWig failed"; rm -f "$BG4"; return 1; }
    echo "      ✅  $(basename $BW)"
}

# =============================================================================
# Main: merge BAMs + pileups (fractional) + BigWigs per group
# =============================================================================
echo "============================================================"
echo " Merging replicates + runs for heatmap generation"
echo " Groups file : $GROUPS_FILE"
echo " TSS BED     : $TSS_BED"
echo " Output      : $OUT_DIR"
echo "============================================================"

PROCESSED_GROUPS=""

while read -r LINE; do
    [[ -z "$LINE" || "$LINE" == \#* ]] && continue

    GROUP=$(echo "$LINE" | awk '{print $1}')
    BARCODES=$(echo "$LINE" | awk '{$1=""; print $0}' | xargs)

    echo ""
    echo "╔══════════════════════════════════════════════════════╗"
    echo "  Group: $GROUP  (barcodes: $BARCODES)"
    echo "╚══════════════════════════════════════════════════════╝"

    GROUP_DIR="$OUT_DIR/${GROUP}"
    MERGED_BAM="$OUT_DIR/${GROUP}_merged.nucs.bam"
    mkdir -p "$GROUP_DIR"

    # ── Collect matching nucs.bam files from both runs ───────────
    BAMS_TO_MERGE=""
    for RUN in ONTRun9 ONTRun10; do
        for BARCODE in $BARCODES; do
            for SAMPLE_DIR in "$FT_RESULTS/$RUN"/*${BARCODE}*/; do
                [ -d "$SAMPLE_DIR" ] || continue
                NUCS_BAM=$(find "$SAMPLE_DIR" -maxdepth 1 -name "*.nucs.bam" | head -1)
                if [ -n "$NUCS_BAM" ] && [ -f "$NUCS_BAM" ]; then
                    echo "  Found: $(basename $NUCS_BAM)"
                    BAMS_TO_MERGE="$BAMS_TO_MERGE $NUCS_BAM"
                else
                    echo "  ⚠️  No nucs.bam in $SAMPLE_DIR"
                fi
            done
        done
    done

    BAMS_TO_MERGE="${BAMS_TO_MERGE# }"

    if [ -z "$BAMS_TO_MERGE" ]; then
        echo "  ❌  No nucs.bam files found for $GROUP — skipping"
        continue
    fi

    BAM_COUNT=$(echo "$BAMS_TO_MERGE" | wc -w)
    echo "  Total: $BAM_COUNT nucs.bam(s) to merge"

    # ── Merge ────────────────────────────────────────────────────
    if [ ! -f "$MERGED_BAM" ]; then
        if [ "$BAM_COUNT" -eq 1 ]; then
            cp $BAMS_TO_MERGE "$MERGED_BAM"
        else
            samtools merge -f -@ 8 "$MERGED_BAM" $BAMS_TO_MERGE
        fi
        samtools index "$MERGED_BAM"
        echo "  ✅  Merged: $(basename $MERGED_BAM)"
    else
        echo "  ⏭️  Merged BAM already exists"
    fi

    # ── Pileups — FRACTIONAL (no --fiber-coverage) ───────────────
    # Without --fiber-coverage, ft pileup outputs the FRACTION of
    # reads with m6A / nucleosome / CpG at each base (0-1 scale).
    # This is what you want for biological interpretation — not raw depth.
    echo "  Running fractional pileups..."

    M6A_BG="$GROUP_DIR/${GROUP}.m6Apileup.bedgraph"
    CPG_BG="$GROUP_DIR/${GROUP}.5mcpileup.bedgraph"
    NUC_BG="$GROUP_DIR/${GROUP}.nucspileup.bedgraph"

    if [ ! -f "$M6A_BG" ]; then
        echo "    ft pileup — m6A fraction..."
        ft pileup --m6a --per-base \
            --out "$M6A_BG" "$MERGED_BAM" \
            || echo "  ❌ m6A pileup failed"
    else
        echo "    ⏭️  m6A pileup exists"
    fi

    if [ ! -f "$CPG_BG" ]; then
        echo "    ft pileup — 5mC fraction..."
        ft pileup --cpg --per-base \
            --out "$CPG_BG" "$MERGED_BAM" \
            || echo "  ❌ 5mC pileup failed"
    else
        echo "    ⏭️  5mC pileup exists"
    fi

    if [ ! -f "$NUC_BG" ]; then
        echo "    ft pileup — nucleosome fraction..."
        ft pileup --nuc --per-base \
            --out "$NUC_BG" "$MERGED_BAM" \
            || echo "  ❌ nuc pileup failed"
    else
        echo "    ⏭️  nuc pileup exists"
    fi

    # ── BigWigs ──────────────────────────────────────────────────
    echo "  Converting to BigWig..."
    NUC_BW="$GROUP_DIR/${GROUP}.nuc.bw"
    M6A_BW="$GROUP_DIR/${GROUP}.m6A.bw"
    CPG_BW="$GROUP_DIR/${GROUP}.5mC.bw"

    make_bigwig "$NUC_BG" "$NUC_BW"
    make_bigwig "$M6A_BG" "$M6A_BW"
    make_bigwig "$CPG_BG" "$CPG_BW"

    PROCESSED_GROUPS="$PROCESSED_GROUPS $GROUP"

done < "$GROUPS_FILE"

# =============================================================================
# Per-group heatmaps — sorted by m6A signal, symmetric 2kb window
# =============================================================================
echo ""
echo "============================================================"
echo " Generating per-group heatmaps"
echo "============================================================"

for GROUP in $PROCESSED_GROUPS; do
    GROUP_DIR="$OUT_DIR/${GROUP}"
    NUC_BW="$GROUP_DIR/${GROUP}.nuc.bw"
    M6A_BW="$GROUP_DIR/${GROUP}.m6A.bw"
    CPG_BW="$GROUP_DIR/${GROUP}.5mC.bw"

    [ -f "$M6A_BW" ] || { echo "  ⚠️  $GROUP BigWigs missing — skipping"; continue; }

    echo ""
    echo "  ── $GROUP ──"

    MATRIX="$GROUP_DIR/${GROUP}.TSS.gz"
    if [ ! -f "$MATRIX" ]; then
        echo "    computeMatrix..."
        computeMatrix reference-point \
            --referencePoint TSS \
            -b 2000 -a 2000 \
            -R "$TSS_BED" \
            -S "$M6A_BW" "$NUC_BW" "$CPG_BW" \
            -o "$MATRIX" \
            --outFileNameMatrix "$GROUP_DIR/${GROUP}.TSS.tab" \
            --missingDataAsZero \
            --samplesLabel "m6A" "nucleosome" "5mC" \
            -p 12 \
            && echo "    ✅  matrix done" \
            || { echo "    ❌  computeMatrix failed"; continue; }
    else
        echo "    ⏭️  matrix exists"
    fi

    # Heatmap — sorted by m6A accessibility (most accessible genes on top)
    plotHeatmap \
        -m "$MATRIX" \
        -out "$GROUP_DIR/${GROUP}.TSS_heatmap.png" \
        --sortUsingSamples 1 \
        --sortRegions descend \
        --colorMap Greens Reds Blues \
        --whatToShow 'heatmap and colorbar' \
        --plotTitle "$GROUP" \
        --heatmapHeight 12 \
        --heatmapWidth 3 \
        --zMin 0 0 0 \
        && echo "    ✅  heatmap done" \
        || echo "    ❌  plotHeatmap failed"

    # Profile — average line plot showing NDR dip + nucleosome phasing
    plotProfile \
        -m "$MATRIX" \
        -out "$GROUP_DIR/${GROUP}.TSS_profile.png" \
        --plotTitle "$GROUP — signal around TSS" \
        --samplesLabel "m6A (accessibility)" "nucleosome" "5mC" \
        --colors green red blue \
        --plotHeight 6 \
        --plotWidth 8 \
        && echo "    ✅  profile done" \
        || echo "    ❌  plotProfile failed"

done

# =============================================================================
# Multi-sample comparison — all strains on one plot per signal type
# Excludes controls (gDNA, HMW)
# =============================================================================
echo ""
echo "============================================================"
echo " Generating multi-sample comparison plots"
echo "============================================================"

COMPARE_GROUPS="WT_Eddie WT_Rochelle cac-1 cac-2 rtt109 rtt109FLAG"

for SIGNAL in m6A nuc 5mC; do

    BW_LIST=""
    LABELS=""
    for GROUP in $COMPARE_GROUPS; do
        case $SIGNAL in
            m6A) BW="$OUT_DIR/${GROUP}/${GROUP}.m6A.bw" ;;
            nuc) BW="$OUT_DIR/${GROUP}/${GROUP}.nuc.bw" ;;
            5mC) BW="$OUT_DIR/${GROUP}/${GROUP}.5mC.bw" ;;
        esac
        [ -f "$BW" ] && BW_LIST="$BW_LIST $BW" && LABELS="$LABELS $GROUP"
    done

    [ -z "$BW_LIST" ] && { echo "  ⚠️  No BigWigs for $SIGNAL — skipping"; continue; }

    echo ""
    echo "  ── All strains: $SIGNAL ──"

    MATRIX="$OUT_DIR/allgroups.${SIGNAL}.TSS.gz"
    if [ ! -f "$MATRIX" ]; then
        computeMatrix reference-point \
            --referencePoint TSS \
            -b 2000 -a 2000 \
            -R "$TSS_BED" \
            -S $BW_LIST \
            -o "$MATRIX" \
            --missingDataAsZero \
            --samplesLabel $LABELS \
            -p 12 \
            && echo "    ✅  matrix done" \
            || { echo "    ❌  computeMatrix failed"; continue; }
    else
        echo "    ⏭️  matrix exists"
    fi

    case $SIGNAL in
        m6A) CMAP="Greens" ;;
        nuc) CMAP="Reds"   ;;
        5mC) CMAP="Blues"  ;;
    esac

    # All-strain profile (one line per strain — key comparison figure)
    plotProfile \
        -m "$MATRIX" \
        -out "$OUT_DIR/allgroups.${SIGNAL}.TSS_profile.png" \
        --plotTitle "${SIGNAL} signal around TSS — all strains" \
        --plotHeight 6 \
        --plotWidth 10 \
        --perGroup \
        && echo "    ✅  comparison profile done" \
        || echo "    ❌  plotProfile failed"

    # All-strain heatmap (one column per strain, same color scale)
    plotHeatmap \
        -m "$MATRIX" \
        -out "$OUT_DIR/allgroups.${SIGNAL}.TSS_heatmap.png" \
        --sortUsingSamples 1 \
        --sortRegions descend \
        --colorMap $CMAP \
        --whatToShow 'heatmap and colorbar' \
        --plotTitle "${SIGNAL} — all strains" \
        --heatmapHeight 12 \
        --heatmapWidth 2 \
        && echo "    ✅  comparison heatmap done" \
        || echo "    ❌  plotHeatmap failed"

done

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "============================================================"
echo " All done. Outputs:"
echo "============================================================"
find "$OUT_DIR" -name "*.png" | sort | while read -r PNG; do
    echo "  $(basename $PNG)"
done
