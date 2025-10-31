#!/usr/bin/env python
#SBATCH --job-name=FixMetas
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=400gb
#SBATCH --time=4:00:00
#SBATCH --output=../FixMetas.%j.out
#SBATCH --error=../FixMetas.%j.err
module load Python/3.11.3-GCCcore-12.3.0
module load SciPy-bundle/2024.05-gfbf-2024a

META_FILE="/scratch/ry00555/RNASeqPaper/Oct2025/BAM_File_Metadata_with_index_merged_V2.csv"
BAM_DIR="/scratch/ry00555/RNASeqPaper/Oct2025/SortedBamFiles"
OUTPUT_REPORT="bam_meta_check.tsv"

# Get a list of BAM files in BAMDIR (basename only)
ls "$BAMDIR"/*.bam | xargs -n1 basename > /tmp/bamdir_files.txt

# Read meta file
meta_df = pd.read_csv(META_FILE)

# Get list of BAM files in directory
bam_files = os.listdir(BAM_DIR)

report_rows = []

for idx, row in meta_df.iterrows():
    meta_bam = str(row['bamReads']).strip()

    # Check if exact file exists
    if meta_bam in bam_files:
        continue  # all good

    # Try to find a close match in directory
    suggested = [f for f in bam_files if meta_bam.replace(".fq.gz","") in f]
    suggested_name = suggested[0] if suggested else ""

    report_rows.append({
        "Row": idx+2,  # +2 to match Excel row numbering (header + 1-index)
        "SampleID": row['SampleID'],
        "bamReads_meta": meta_bam,
        "Suggested_BAM": suggested_name
    })

# Create report
if report_rows:
    report_df = pd.DataFrame(report_rows)
    report_df.to_csv(OUTPUT_REPORT, index=False)
    print(f"✅ BAM check report saved to {OUTPUT_REPORT}")
else:
    print("✅ All BAMs in meta exist in BAM directory. No issues found.")
