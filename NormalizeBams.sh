#!/bin/bash
#SBATCH --job-name=NormalizeBams
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=8:00:00
#SBATCH --output=../NormalizeBams.%j.out
#SBATCH --error=../NormalizeBams.%j.err

ml deepTools

BAMDIR="/scratch/ry00555/RNASeqPaper/Oct2025/SortedBamFiles"
META="/scratch/ry00555/RNASeqPaper/Oct2025/BAM_File_Metadata_with_index_merged_V2.csv"
OUTDIR="/scratch/ry00555/RNASeqPaper/Oct2025/SortedBamFiles/DeepToolsNormalized"

mkdir -p "$OUTDIR"

echo "🚀 Starting deepTools normalization..."

tail -n +2 "$META" | while IFS=, read -r RunID bamReads BamIndex SampleID Factor Tissue Condition Replicate bamControl bamInputIndex ControlID Peaks PeakCaller DesiredPeakName MACS3minlength MACS3maxgap; do
    [[ -z "$RunID" ]] && continue

    chip_path="${BAMDIR}/${bamReads}"
    input_path="${BAMDIR}/${bamControl}"

    if [[ ! -f "$chip_path" ]]; then
        echo "⚠ Missing: $chip_path"
        continue
    fi

    # Output basename
    if [[ -f "$input_path" ]]; then
        outname="${DesiredPeakName}_log2FC.bw"
        outfile="${OUTDIR}/${outname}"

        # ✅ Skip if already exists
        if [[ -f "$outfile" ]]; then
            echo "✅ Skipping (already exists): $outfile"
            continue
        fi
        echo "📌 Normalizing: $chip_path  vs  $input_path"

        bamCompare \
            -b1 "$chip_path" \
            -b2 "$input_path" \
            --operation log2 \
            --pseudocount 1 \
            --binSize 25 \
            --skipZeroOverZero \
            --outFileFormat bigwig \
            -o "${OUTDIR}/${outname}"

    else
        outname="${DesiredPeakName}_norm.bw"
        echo "⚠ No input for $chip_path → coverage only"

        bamCoverage \
            -b "$chip_path" \
            --binSize 25 \
            --normalizeUsing BPM \
            -o "${OUTDIR}/${outname}"
    fi
done

multiBigwigSummary BED-file \
    --bwfiles ${OUTDIR}/*.bw \
    --BED "/scratch/ry00555/GeneList_BedFiles/K27genes.bed" \
    -out ${OUTDIR}/K27genes_signal_matrix.npz \
    --outRawCounts ${OUTDIR}/K27genes_signal_matrix.tab

    multiBigwigSummary BED-file \
        --bwfiles ${OUTDIR}/*.bw \
        --BED "/scratch/ry00555/GeneList_BedFiles/NonK27genes.bed" \
        -out ${OUTDIR}/NonK27genes_signal_matrix.npz \
        --outRawCounts ${OUTDIR}/NonK27genes_signal_matrix.tab

        multiBigwigSummary BED-file \
            --bwfiles ${OUTDIR}/*.bw \
            --BED "/scratch/ry00555/Figure2G_K27regions_Scaledcenter_FileToCheckOrderFINAL.txt" \
            -out ${OUTDIR}/K27regions_signal_matrix.npz \
            --outRawCounts ${OUTDIR}/K27regions_signal_matrix.tab
