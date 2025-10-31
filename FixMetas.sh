#!/bin/bash
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

META="/scratch/ry00555/RNASeqPaper/Oct2025/BAM_File_Metadata_with_index_merged_V2.csv"
BAMDIR="/scratch/ry00555/RNASeqPaper/Oct2025/SortedBamFiles"
REPORT="bam_meta_check.tsv"

# Get a list of BAM files in BAMDIR (basename only)
ls "$BAMDIR"/*.bam | xargs -n1 basename > /tmp/bamdir_files.txt

# Run Python to check
python3 <<'EOF'
import pandas as pd

# Load metadata
meta_file = "$META"
meta = pd.read_csv(meta_file)

# Load BAM filenames from BAMDIR
with open("/tmp/bamdir_files.txt") as f:
    bam_files = [x.strip() for x in f.readlines()]

# Clean BAM names in the directory for matching: remove .bam, .fq.gz, _Q30
bam_clean = {}
for bf in bam_files:
    clean = bf.replace('.bam','').replace('.fq.gz','').replace('_Q30','')
    bam_clean[clean] = bf

# Check each row
rows_to_fix = []
for idx, row in meta.iterrows():
    bam_name = row['bamReads']
    # Clean metadata bamReads similarly
    clean_meta = bam_name.replace('.bam','').replace('.fq.gz','').replace('_Q30','')
    if clean_meta not in bam_clean:
        rows_to_fix.append({'Row': idx+2, 'SampleID': row['SampleID'], 'bamReads_meta': bam_name})

# Save report
import csv
with open("$REPORT", 'w') as out:
    writer = csv.DictWriter(out, fieldnames=['Row','SampleID','bamReads_meta'])
    writer.writeheader()
    for r in rows_to_fix:
        writer.writerow(r)

print(f"✅ Report generated: {REPORT}")
EOF
