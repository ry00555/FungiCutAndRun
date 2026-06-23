#!/usr/bin/env python
# -*- coding: utf-8 -*-

cd /scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek

python3 - << 'EOF'
import requests, time, re, tarfile, os
from pathlib import Path

# ── Config ────────────────────────────────────────────────────────────
ORGANISMS = {
    "ncr": ("ncr_proteome.fasta", "neurospora"),
    "fgr": ("fgr_proteome.fasta", "fusarium"),
    "mgr": ("mgr_proteome.fasta", "magnaporthe"),
    "zt":  ("zt_proteome.fasta",  "zymoseptoria"),
    "cne": ("cne_proteome.fasta", "cryptococcus"),
}
AFDB_API  = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
CIF_DIR   = Path("cif_fungi")
CIF_DIR.mkdir(exist_ok=True)

session = requests.Session()
session.headers["User-Agent"] = "FungalProteomeFetcher/1.0 (ry00555@uga.edu)"

def get_accs_from_fasta(fasta_path):
    """Extract UniProt accessions from FASTA headers."""
    accs = []
    for line in Path(fasta_path).read_text().splitlines():
        if not line.startswith(">"): continue
        # UniProt format: >sp|ACC|NAME or >tr|ACC|NAME
        m = re.search(r"[st][rp]\|([A-Z0-9]+)\|", line)
        if m:
            accs.append(m.group(1))
    return accs

for prefix, (fasta_file, label) in ORGANISMS.items():
    print(f"\n{'='*60}")
    print(f"Processing {label} ({fasta_file})")
    accs = get_accs_from_fasta(fasta_file)
    print(f"  {len(accs)} accessions found in FASTA")

    org_dir = CIF_DIR / label
    org_dir.mkdir(exist_ok=True)

    fetched = skipped = failed = 0
    for i, acc in enumerate(accs, 1):
        out_cif = org_dir / f"AF-{acc}-F1-model_v4.cif"
        if out_cif.exists():
            skipped += 1
            continue

        # Query AFDB API for real URL (handles merged/secondary accessions)
        try:
            r = session.get(AFDB_API.format(acc=acc), timeout=15)
            if r.status_code == 200:
                entry = r.json()
                if entry:
                    cif_url = entry[0].get("cifUrl", "")
                    if cif_url:
                        cif_r = session.get(cif_url, timeout=30)
                        if cif_r.status_code == 200 and cif_r.content[:5] == b"data_":
                            out_cif.write_bytes(cif_r.content)
                            fetched += 1
                        else:
                            failed += 1
                    else:
                        failed += 1
                else:
                    failed += 1
            else:
                failed += 1
        except Exception as e:
            failed += 1

        time.sleep(0.3)

        if i % 500 == 0:
            print(f"  [{i}/{len(accs)}] fetched={fetched} skipped={skipped} failed={failed}")

    print(f"  Done: fetched={fetched} skipped={skipped} failed={failed}")
    total_cif = len(list(org_dir.glob("*.cif")))
    print(f"  CIF files on disk: {total_cif}")

    # Pack into tar
    tar_path = Path(f"{label}_cif_v4.tar")
    print(f"  Packing -> {tar_path}...")
    with tarfile.open(tar_path, "w") as tf:
        for cif in sorted(org_dir.glob("*.cif")):
            tf.add(cif, arcname=cif.name)
    print(f"  Tar size: {tar_path.stat().st_size / 1e9:.2f} GB")

print("\nAll done.")
EOF
