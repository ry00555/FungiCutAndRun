#!/usr/bin/env bash
#SBATCH --job-name=FiberTools_MergedQC
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=FiberTools_MergedQC.%j.out
#SBATCH --error=FiberTools_MergedQC.%j.err

# =============================================================================
# Run ft qc on strain-merged nucs.bam files + collect into merged_fiberqc.txt
# =============================================================================
# Depends on FiberTools_MergeForHeatmaps.sh having run first.
# Outputs:
#   - Per-strain QC files in MergedForHeatmaps/<strain>/<strain>_merged.qc.txt
#   - Combined summary: ~/merged_fiberqc.txt
# =============================================================================

set -euo pipefail

module load Miniforge3/24.11.3-0
source activate /home/ry00555/fibertools

MERGED_DIR="/lustre2/scratch/ry00555/ONTRun9_10Combined/MergedForHeatmaps"
COMBINED_QC="$HOME/merged_fiberqc.txt"

echo "============================================================"
echo " FiberTools merged QC"
echo " Input : $MERGED_DIR"
echo " Output: $COMBINED_QC"
echo "============================================================"

for BAM in "$MERGED_DIR"/*_merged.nucs.bam; do
    [ -f "$BAM" ] || continue

    GROUP=$(basename "$BAM" _merged.nucs.bam)
    GROUP_DIR="$MERGED_DIR/${GROUP}"
    QC_OUT="$GROUP_DIR/${GROUP}_merged.qc.txt"

    mkdir -p "$GROUP_DIR"

    echo ""
    echo "  ── $GROUP ──"

    if [ ! -f "${BAM}.bai" ]; then
        echo "    Indexing BAM..."
        samtools index "$BAM"
    fi

    if [ ! -f "$QC_OUT" ]; then
        echo "    ft qc..."
        ft qc --m6a-per-msp "$BAM" "$QC_OUT" \
            && echo "    ✅  $(basename $QC_OUT)" \
            || echo "    ❌  qc failed for $GROUP"
    else
        echo "    ⏭️  QC already exists — skipping"
    fi

done

# =============================================================================
# Collect all merged QC files into one combined table
# Columns: sample | statistic | value | count
# Excludes per-base histograms (fiber_length, msp_length, etc.)
# =============================================================================
echo ""
echo "============================================================"
echo " Collecting QC into $COMBINED_QC"
echo "============================================================"

# Write header
echo -e "sample\tstatistic\tvalue\tcount" > "$COMBINED_QC"

for QC in "$MERGED_DIR"/*/*_merged.qc.txt; do
    [ -f "$QC" ] || continue
    SAMPLE=$(basename "$QC" .qc.txt)
    grep -v "fiber_length\|msp_length\|nuc_length\|m6a_per_msp_size\|^statistic" "$QC" \
        | awk -v s="$SAMPLE" '{print s"\t"$0}'
done >> "$COMBINED_QC"

NSAMPLES=$(tail -n +2 "$COMBINED_QC" | awk '{print $1}' | sort -u | wc -l)
NROWS=$(wc -l < "$COMBINED_QC")
echo "  ✅  Written: $COMBINED_QC"
echo "      $NSAMPLES samples, $NROWS rows total"
echo ""
echo "  Key metrics to check:"
grep "m6a_ratio\|phased_reads" "$COMBINED_QC" | awk '{
    if ($2=="phased_reads") printf "  %-40s  total_fibers=%s\n", $1, $4
    if ($2=="m6a_ratio" && NR<20) printf "  %-40s  m6a_ratio=%s%%\n", $1, $3
}' | head -20
