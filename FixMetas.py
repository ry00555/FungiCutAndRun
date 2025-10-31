#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pandas as pd

# ==== USER SETTINGS ====
META_FILE = "/scratch/ry00555/RNASeqPaper/Oct2025/BAM_File_Metadata_with_index_merged_V2.csv"
BAM_DIR = "/scratch/ry00555/RNASeqPaper/Oct2025/SortedBamFiles"
OUTPUT_REPORT = "bam_meta_check.tsv"
# ========================

# Get list of BAM files in directory (basename only)
bam_files = [f for f in os.listdir(BAM_DIR) if f.endswith(".bam")]

# Read meta CSV
meta_df = pd.read_csv(META_FILE)

report_rows = []

for idx, row in meta_df.iterrows():
    meta_bam = str(row['bamReads']).strip()

    # Check if exact file exists
    if meta_bam in bam_files:
        continue  # all good

    # Try to find a close match in directory
    suggested = [f for f in bam_files if meta_bam.replace(".fq.gz", "") in f]
    suggested_name = suggested[0] if suggested else ""

    report_rows.append({
        "Row": idx + 2,  # Excel row numbering (header + 1-index)
        "SampleID": row['SampleID'],
        "bamReads_meta": meta_bam,
        "Suggested_BAM": suggested_name
    })

# Create report
if report_rows:
    report_df = pd.DataFrame(report_rows)
    report_df.to_csv(OUTPUT_REPORT, sep="\t", index=False)
    print(f"✅ BAM check report saved to {OUTPUT_REPORT}")
else:
    print("✅ All BAMs in meta exist in BAM directory. No issues found.")
