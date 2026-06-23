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
# normalized TSS heatmaps from both BigWigs AND bedgraphs directly
#
# Bedgraph column layout from ft pileup --m6a --cpg:
#   col1=chrom col2=start col3=end col4=coverage col5=fire_coverage
#   col6=score col7=nuc_coverage col8=msp_coverage col9=m6a_coverage col10=cpg_coverage
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
    echo "❌  Groups file not found: $GROUPS_FILE"; exit 1
fi
if [ ! -f "$TSS_BED" ]; then
    echo "❌  TSS BED not found: $TSS_BED"; exit 1
fi

# =============================================================================
# Helper: extract one signal column from multi-column bedgraph → BigWig
# Col 7=nuc, 8=msp, 9=m6a, 10=cpg
# =============================================================================
make_bigwig() {
    local BG="$1"
    local BW="$2"
    local COL="$3"

    [ -f "$BW" ] && { echo "      ⏭️  $(basename $BW) exists"; return 0; }
    [ -f "$BG" ] || { echo "      ⚠️  $(basename $BG) not found — skipping"; return 0; }

    echo "      Extracting col $COL → $(basename $BW)"

    local BG4="${BG%.bedgraph}.col${COL}.bedgraph"

    awk -v col="$COL" 'NR>1 && $col >= 0 {print $1"\t"$2"\t"$3"\t"$col}' "$BG" \
        | sort -k1,1 -k2,2n \
        | awk -v genome="$GENOME" '
            BEGIN { while ((getline line < genome) > 0) {
                split(line,a,"\t"); chromlen[a[1]]=a[2] } }
            { if ($2<0) $2=0
              if (chromlen[$1] && $3>chromlen[$1]) $3=chromlen[$1]
              if ($2<$3) print $1"\t"$2"\t"$3"\t"$4 }
        ' > "$BG4"

    if [ ! -s "$BG4" ]; then
        echo "      ❌  Empty after column extraction"
        rm -f "$BG4"; return 1
    fi

    bedGraphToBigWig "$BG4" "$GENOME" "$BW" && rm -f "$BG4" \
        || { echo "      ❌ bedGraphToBigWig failed"; rm -f "$BG4"; return 1; }
    echo "      ✅  $(basename $BW)"
}

# =============================================================================
# Helper: extract one signal column from bedgraph → 4-col bedgraph for deepTools
# =============================================================================
make_4col_bedgraph() {
    local BG="$1"
    local OUT="$2"
    local COL="$3"

    [ -f "$OUT" ] && { echo "      ⏭️  $(basename $OUT) exists"; return 0; }
    [ -f "$BG" ] || { echo "      ⚠️  $(basename $BG) not found"; return 0; }

    echo "      Extracting col $COL → $(basename $OUT)"

    awk -v col="$COL" 'NR>1 && $col >= 0 {print $1"\t"$2"\t"$3"\t"$col}' "$BG" \
        | sort -k1,1 -k2,2n \
        | awk -v genome="$GENOME" '
            BEGIN { while ((getline line < genome) > 0) {
                split(line,a,"\t"); chromlen[a[1]]=a[2] } }
            { if ($2<0) $2=0
              if (chromlen[$1] && $3>chromlen[$1]) $3=chromlen[$1]
              if ($2<$3) print $1"\t"$2"\t"$3"\t"$4 }
        ' > "$OUT"

    [ -s "$OUT" ] && echo "      ✅  $(basename $OUT)" || { echo "      ❌  Empty output"; return 1; }
}

# =============================================================================
# Helper: computeMatrix + plotHeatmap + plotProfile
# Works with both BigWig and bedgraph inputs (deepTools accepts both)
# =============================================================================
make_heatmap_and_profile() {
    local MATRIX="$1"
    local LABEL="$2"
    local OUTDIR="$3"
    local SORT_SAMPLE="$4"   # which sample index to sort by (1-based)
    shift 4
    local BW_OR_BG_LIST=("$@")

    if [ ! -f "$MATRIX" ]; then
        echo "    computeMatrix from ${#BW_OR_BG_LIST[@]} tracks..."
        computeMatrix reference-point \
            --referencePoint TSS \
            -b 2000 -a 2000 \
            -R "$TSS_BED" \
            -S "${BW_OR_BG_LIST[@]}" \
            -o "$MATRIX" \
            --outFileNameMatrix "${MATRIX%.gz}.tab" \
            --missingDataAsZero \
            -p 12 \
            && echo "    ✅  matrix done" \
            || { echo "    ❌  computeMatrix failed"; return 1; }
    else
        echo "    ⏭️  matrix exists"
    fi

    plotHeatmap \
        -m "$MATRIX" \
        -out "${OUTDIR}/${LABEL}.TSS_heatmap.png" \
        --sortUsingSamples "$SORT_SAMPLE" \
        --sortRegions descend \
        --colorMap Greens Greens Reds Blues \
        --whatToShow 'heatmap and colorbar' \
        --plotTitle "$LABEL" \
        --heatmapHeight 12 \
        --heatmapWidth 3 \
        --zMin 0 0 0 0 \
        && echo "    ✅  heatmap done" \
        || echo "    ❌  plotHeatmap failed"

    plotProfile \
        -m "$MATRIX" \
        -out "${OUTDIR}/${LABEL}.TSS_profile.png" \
        --plotTitle "${LABEL} — signal around TSS" \
        --colors darkgreen green red blue \
        --plotHeight 6 \
        --plotWidth 8 \
        && echo "    ✅  profile done" \
        || echo "    ❌  plotProfile failed"
}

# =============================================================================
# Main: merge BAMs + pileup + BigWigs + 4-col bedgraphs per group
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

    # ── Collect nucs.bam files ───────────────────────────────────
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
        echo "  ❌  No nucs.bam files found — skipping"
        continue
    fi

    echo "  Total: $(echo $BAMS_TO_MERGE | wc -w) nucs.bam(s)"

    # ── Merge ────────────────────────────────────────────────────
    if [ ! -f "$MERGED_BAM" ]; then
        BAM_COUNT=$(echo "$BAMS_TO_MERGE" | wc -w)
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

    # ── Single pileup with all signal flags ──────────────────────
    PILEUP_BG="$GROUP_DIR/${GROUP}.allpileup.bedgraph"

    if [ ! -f "$PILEUP_BG" ]; then
        echo "  Running pileup (all signals)..."
        ft pileup --m6a --cpg --per-base \
            --out "$PILEUP_BG" "$MERGED_BAM" \
            || { echo "  ❌ pileup failed"; continue; }
    else
        echo "  ⏭️  Pileup exists"
    fi

    # Print columns so we can verify
    echo "  Pileup columns: $(head -1 $PILEUP_BG)"

    # ── Extract per-signal 4-col bedgraphs ───────────────────────
    echo "  Extracting signal bedgraphs..."
    NUC_BG4="$GROUP_DIR/${GROUP}.nuc.4col.bedgraph"
    MSP_BG4="$GROUP_DIR/${GROUP}.msp.4col.bedgraph"
    M6A_BG4="$GROUP_DIR/${GROUP}.m6a.4col.bedgraph"
    CPG_BG4="$GROUP_DIR/${GROUP}.5mC.4col.bedgraph"

    make_4col_bedgraph "$PILEUP_BG" "$NUC_BG4" 7
    make_4col_bedgraph "$PILEUP_BG" "$MSP_BG4" 8
    make_4col_bedgraph "$PILEUP_BG" "$M6A_BG4" 9
    make_4col_bedgraph "$PILEUP_BG" "$CPG_BG4" 10

    # ── BigWigs ──────────────────────────────────────────────────
    echo "  Converting to BigWig..."
    NUC_BW="$GROUP_DIR/${GROUP}.nuc.bw"
    MSP_BW="$GROUP_DIR/${GROUP}.msp.bw"
    M6A_BW="$GROUP_DIR/${GROUP}.m6A.bw"
    CPG_BW="$GROUP_DIR/${GROUP}.5mC.bw"

    make_bigwig "$PILEUP_BG" "$NUC_BW" 7
    make_bigwig "$PILEUP_BG" "$MSP_BW" 8
    make_bigwig "$PILEUP_BG" "$M6A_BW" 9
    make_bigwig "$PILEUP_BG" "$CPG_BW" 10

    PROCESSED_GROUPS="$PROCESSED_GROUPS $GROUP"

done < "$GROUPS_FILE"

# =============================================================================
# Per-group heatmaps — from BigWig AND from bedgraph directly
# =============================================================================
echo ""
echo "============================================================"
echo " Generating per-group heatmaps (BigWig + bedgraph)"
echo "============================================================"

for GROUP in $PROCESSED_GROUPS; do
    GROUP_DIR="$OUT_DIR/${GROUP}"

    NUC_BW="$GROUP_DIR/${GROUP}.nuc.bw"
    MSP_BW="$GROUP_DIR/${GROUP}.msp.bw"
    M6A_BW="$GROUP_DIR/${GROUP}.m6A.bw"
    CPG_BW="$GROUP_DIR/${GROUP}.5mC.bw"

    NUC_BG4="$GROUP_DIR/${GROUP}.nuc.4col.bedgraph"
    MSP_BG4="$GROUP_DIR/${GROUP}.msp.4col.bedgraph"
    M6A_BG4="$GROUP_DIR/${GROUP}.m6a.4col.bedgraph"
    CPG_BG4="$GROUP_DIR/${GROUP}.5mC.4col.bedgraph"

    echo ""
    echo "  ── $GROUP ──"

    # --- Heatmaps from BigWig ---
    if [ -f "$M6A_BW" ]; then
        echo "    [BigWig heatmaps]"
        MATRIX_BW="$GROUP_DIR/${GROUP}.bw.TSS.gz"
        make_heatmap_and_profile \
            "$MATRIX_BW" "${GROUP}_BigWig" "$GROUP_DIR" 2 \
            "$M6A_BW" "$MSP_BW" "$NUC_BW" "$CPG_BW"
    else
        echo "    ⚠️  BigWigs missing — skipping BigWig heatmaps"
    fi

    # --- Heatmaps from bedgraph directly ---
    if [ -f "$M6A_BG4" ]; then
        echo "    [bedgraph heatmaps]"
        MATRIX_BG="$GROUP_DIR/${GROUP}.bg.TSS.gz"
        make_heatmap_and_profile \
            "$MATRIX_BG" "${GROUP}_bedgraph" "$GROUP_DIR" 2 \
            "$M6A_BG4" "$MSP_BG4" "$NUC_BG4" "$CPG_BG4"
    else
        echo "    ⚠️  4-col bedgraphs missing — skipping bedgraph heatmaps"
    fi

done

# =============================================================================
# Multi-sample comparison — all strains on one plot per signal type
# Both BigWig and bedgraph versions
# =============================================================================
echo ""
echo "============================================================"
echo " Multi-sample comparison plots"
echo "============================================================"

COMPARE_GROUPS="WT_Eddie WT_Rochelle cac-1 cac-2 rtt109 rtt109FLAG"

for SIGNAL in m6A msp nuc 5mC; do
    for MODE in bw bg; do

        BW_LIST=""
        LABELS=""
        for GROUP in $COMPARE_GROUPS; do
            if [ "$MODE" = "bw" ]; then
                case $SIGNAL in
                    m6A) F="$OUT_DIR/${GROUP}/${GROUP}.m6A.bw" ;;
                    msp) F="$OUT_DIR/${GROUP}/${GROUP}.msp.bw" ;;
                    nuc) F="$OUT_DIR/${GROUP}/${GROUP}.nuc.bw" ;;
                    5mC) F="$OUT_DIR/${GROUP}/${GROUP}.5mC.bw" ;;
                esac
            else
                case $SIGNAL in
                    m6A) F="$OUT_DIR/${GROUP}/${GROUP}.m6a.4col.bedgraph" ;;
                    msp) F="$OUT_DIR/${GROUP}/${GROUP}.msp.4col.bedgraph" ;;
                    nuc) F="$OUT_DIR/${GROUP}/${GROUP}.nuc.4col.bedgraph" ;;
                    5mC) F="$OUT_DIR/${GROUP}/${GROUP}.5mC.4col.bedgraph" ;;
                esac
            fi
            [ -f "$F" ] && BW_LIST="$BW_LIST $F" && LABELS="$LABELS $GROUP"
        done

        [ -z "$BW_LIST" ] && continue

        echo ""
        echo "  ── All strains: $SIGNAL ($MODE) ──"

        MATRIX="$OUT_DIR/allgroups.${SIGNAL}.${MODE}.TSS.gz"
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
            m6A|msp) CMAP="Greens" ;;
            nuc)     CMAP="Reds"   ;;
            5mC)     CMAP="Blues"  ;;
        esac

        plotProfile \
            -m "$MATRIX" \
            -out "$OUT_DIR/allgroups.${SIGNAL}.${MODE}.TSS_profile.png" \
            --plotTitle "${SIGNAL} (${MODE}) — all strains" \
            --plotHeight 6 --plotWidth 10 \
            --perGroup \
            && echo "    ✅  profile done" \
            || echo "    ❌  plotProfile failed"

        plotHeatmap \
            -m "$MATRIX" \
            -out "$OUT_DIR/allgroups.${SIGNAL}.${MODE}.TSS_heatmap.png" \
            --sortUsingSamples 1 \
            --sortRegions descend \
            --colorMap $CMAP \
            --whatToShow 'heatmap and colorbar' \
            --plotTitle "${SIGNAL} (${MODE}) — all strains" \
            --heatmapHeight 12 --heatmapWidth 2 \
            && echo "    ✅  heatmap done" \
            || echo "    ❌  plotHeatmap failed"

    done
done

echo ""
echo "============================================================"
echo " All done. Outputs:"
echo "============================================================"
find "$OUT_DIR" -name "*.png" | sort | while read -r PNG; do
    echo "  $(basename $PNG)"
done
