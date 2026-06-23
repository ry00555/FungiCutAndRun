#!/bin/bash
#SBATCH --job-name=foldseek_readers
#SBATCH --partition=gpu_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --gres=gpu:A100:1
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/logs/foldseek_%j.out
#SBATCH --error=/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/logs/foldseek_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ry00555@uga.edu

# ══════════════════════════════════════════════════════════════════════
# 04_run_foldseek.sh
#
# STEP 0 — HMMER domain boundary scan (Pfam-A)
#   Runs hmmscan on all sequences to get per-protein domain coordinates.
#   Outputs:
#     pfam/hmmscan_domtbl.txt        raw hmmscan domain table
#     pfam/domain_boundaries.csv     cleaned per-protein domain coords
#                                    (gene, domain, env_start, env_end)
#
# STEP 1 — Organise structures into per-domain symlink subdirectories
#
# STEP 2 — Per-domain FoldSeek: TSV + HTML + superposed PDBs
#
# STEP 3 — Full all-vs-all (TSV + clustering, no HTML — too large)
#
# STEP 4 — Build HTML index page linking all domain viewers
# ══════════════════════════════════════════════════════════════════════

set -euo pipefail

module load Foldseek/10-941cd33-GPU
module load HMMER/3.4-gompi-2024a

SCRATCH="/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity"
PDB_DIR="${SCRATCH}/pdb_lists/pdbs"
MANIFEST="${SCRATCH}/pdb_lists/structure_manifest.tsv"
FASTA="${SCRATCH}/fastas/all_readers_v4.fasta"
PFAM_DIR="${SCRATCH}/pfam"
RESULTS="${SCRATCH}/results"
BY_DOMAIN="${RESULTS}/by_domain"
TMP="/tmp/foldseek_$$"
NCPU="${SLURM_CPUS_PER_TASK}"

mkdir -p "${BY_DOMAIN}" "${RESULTS}" "${TMP}" "${PFAM_DIR}"

echo "════════════════════════════════════════════════════════════"
echo " FoldSeek — Per-domain + Full all-vs-all"
echo " Scratch  : ${SCRATCH}"
echo " Job ID   : ${SLURM_JOB_ID}"
echo " Node     : ${SLURM_NODELIST}"
echo " CPUs     : ${NCPU}"
echo " Started  : $(date)"
echo "════════════════════════════════════════════════════════════"

N_STRUCT=$(find "${PDB_DIR}" \( -name "*.pdb" -o -name "*.cif" \) 2>/dev/null | wc -l)
echo "Structure files: ${N_STRUCT}"
[[ $N_STRUCT -eq 0 ]] && { echo "ERROR: No structures in ${PDB_DIR}"; exit 1; }

# ══════════════════════════════════════════════════════════════════════
# STEP 0 — HMMER domain boundary scan
# ══════════════════════════════════════════════════════════════════════
echo ""
echo "══ STEP 0: HMMER domain boundary scan ═════════════════════"

DOMTBL="${PFAM_DIR}/hmmscan_domtbl.txt"
BOUNDARIES="${PFAM_DIR}/domain_boundaries.csv"

# Download and press Pfam-A HMM database if not already there
if [[ ! -f "${PFAM_DIR}/Pfam-A.hmm.h3i" ]]; then
    echo "  Downloading Pfam-A HMM database..."
    if [[ ! -f "${PFAM_DIR}/Pfam-A.hmm" ]]; then
        wget -q -O "${PFAM_DIR}/Pfam-A.hmm.gz" \
            https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
        gunzip "${PFAM_DIR}/Pfam-A.hmm.gz"
        echo "  Downloaded Pfam-A.hmm"
    fi
    echo "  Pressing Pfam-A HMM database..."
    hmmpress "${PFAM_DIR}/Pfam-A.hmm"
    echo "  Pfam-A database ready"
else
    echo "  Pfam-A database already pressed"
fi

# Run hmmscan
if [[ ! -f "${DOMTBL}" ]]; then
    echo "  Running hmmscan on $(grep -c '^>' "${FASTA}") sequences..."
    hmmscan \
        --domtblout "${DOMTBL}" \
        --cpu "${NCPU}" \
        -E 0.001 \
        --domE 0.001 \
        --noali \
        "${PFAM_DIR}/Pfam-A.hmm" \
        "${FASTA}" \
        > "${PFAM_DIR}/hmmscan.log" 2>&1
    echo "  hmmscan done: $(grep -v '^#' "${DOMTBL}" | wc -l) domain hits"
else
    echo "  hmmscan output already exists: ${DOMTBL}"
fi

# Parse domain table into clean CSV with residue boundaries
python3 - << 'PYEOF'
import re, csv
from pathlib import Path

SCRATCH  = Path("/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity")
DOMTBL   = SCRATCH / "pfam/hmmscan_domtbl.txt"
MASTER   = SCRATCH / "master_reader_proteins_v4.csv"
OUT_CSV  = SCRATCH / "pfam/domain_boundaries.csv"

# Map Pfam accession/name -> our domain label
# Covers the 15 structural domain classes used in this study
PFAM_TO_DOMAIN = {
    # Chromodomain
    "Chromo":          "Chromo", "Chromo_shadow": "Chromo",
    # PWWP
    "PWWP":            "PWWP",
    # Tudor
    "Tudor":           "Tudor", "Tudor_2": "Tudor",
    # MRG
    "MRG":             "MRG",
    # SET
    "SET":             "SET",
    # WD40
    "WD40":            "WD40", "WD_repeats_1": "WD40",
    # JmjC
    "JmjC":            "JmjC",
    # CXXC
    "zf-CXXC":         "CXXC",
    # PHD
    "PHD":             "PHD", "PHD_2": "PHD", "PHD_3": "PHD",
    # BAH
    "BAH":             "BAH",
    # Bromodomain
    "Bromodomain":     "Bromodomain",
    # Zinc_finger (C2H2-type and related)
    "zf-C2H2":        "Zinc_finger", "zf-RING_2": "Zinc_finger",
    # LSD1
    "FAD_binding_2":   "LSD1", "Amino_oxidase": "LSD1",
    # MBT
    "MBT":             "MBT",
    # SWIRM
    "SWIRM":           "SWIRM",
}

# Load gene name lookup from FASTA header -> gene name in master CSV
import re as _re

def header_to_gene(header):
    """Extract gene name from FASTA header (first field before |)."""
    return header.lstrip(">").split("|")[0].strip()

# Parse hmmscan domtbl
records = []
with open(DOMTBL) as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        if len(parts) < 23:
            continue
        pfam_name  = parts[0]          # Pfam HMM name
        gene_id    = parts[3]          # query sequence name (from FASTA header)
        seq_eval   = float(parts[6])   # sequence e-value
        dom_eval   = float(parts[11])  # domain e-value
        env_start  = int(parts[19])    # envelope start (residue coords)
        env_end    = int(parts[20])    # envelope end
        dom_score  = float(parts[13])  # domain score

        if dom_eval > 1e-3:
            continue
        domain = PFAM_TO_DOMAIN.get(pfam_name)
        if domain is None:
            continue

        gene = header_to_gene(gene_id)
        records.append({
            "gene":       gene,
            "pfam_name":  pfam_name,
            "domain":     domain,
            "env_start":  env_start,
            "env_end":    env_end,
            "dom_score":  dom_score,
            "dom_evalue": dom_eval,
        })

# Write output
with open(OUT_CSV, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["gene","pfam_name","domain",
                                        "env_start","env_end",
                                        "dom_score","dom_evalue"])
    w.writeheader()
    w.writerows(records)

print(f"  Domain boundaries: {len(records)} domain hits across "
      f"{len(set(r['gene'] for r in records))} proteins")
print(f"  Saved -> {OUT_CSV}")

# Summary by domain
from collections import Counter
dom_counts = Counter(r["domain"] for r in records)
for dom, n in sorted(dom_counts.items(), key=lambda x: -x[1]):
    genes = len(set(r["gene"] for r in records if r["domain"] == dom))
    print(f"    {dom:<20} {n:>4} hits  ({genes} proteins)")
PYEOF

echo ""
echo "  Domain boundary CSV: ${BOUNDARIES}"

# ══════════════════════════════════════════════════════════════════════
# STEP 1 — Organise structures into per-domain symlink subdirectories
#          Uses predicted_functional_domains (multi-domain aware)
#          so EPR-1 (BAH;PHD) appears in both BAH/ and PHD/ dirs
# ══════════════════════════════════════════════════════════════════════
echo ""
echo "══ STEP 1: Organise structures by domain ══════════════════"

python3 - << 'PYEOF'
import csv, re
from pathlib import Path

SCRATCH  = Path("/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity")
MANIFEST = SCRATCH / "pdb_lists/structure_manifest.tsv"
MASTER   = SCRATCH / "master_reader_proteins_v4.csv"
BY_DOM   = SCRATCH / "pdb_lists/by_domain"

KNOWN_DOMAINS = {
    'BAH','Bromodomain','CXXC','Chromo','JmjC','LSD1','MBT',
    'MRG','PHD','PWWP','SET','SWIRM','Tudor','WD40','Zinc_finger'
}

def san(s):
    return re.sub(r"[^A-Za-z0-9_\-\.]", "_", s)

# Load multi-domain annotations from master CSV
# predicted_functional_domains = "BAH;PHD" etc.
gene_to_domains = {}
with open(MASTER, newline="") as fh:
    for row in csv.DictReader(fh):
        gene = row["protein_gene_name"].strip()
        raw  = row.get("predicted_functional_domains", "")
        domains = [d.strip() for d in raw.split(";") if d.strip() in KNOWN_DOMAINS]
        # Fall back to reader_domain if no known domains found
        if not domains:
            rd = row.get("reader_domain","").strip()
            if rd in KNOWN_DOMAINS:
                domains = [rd]
            else:
                domains = ["unknown"]
        gene_to_domains[gene] = domains

# Load structure file paths from manifest
gene_to_struct = {}
with open(MANIFEST, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        src  = row.get("structure_file", "MISSING").strip()
        gene = row.get("protein_gene_name", "").strip()
        if src in ("MISSING", "") or not gene:
            continue
        p = Path(src)
        if p.exists():
            gene_to_struct[gene] = p

# Create symlinks: each protein symlinked into all its domain dirs
domain_counts = {}
for gene, struct_path in gene_to_struct.items():
    domains = gene_to_domains.get(gene, ["unknown"])
    for dom in domains:
        dom_dir = BY_DOM / san(dom)
        dom_dir.mkdir(parents=True, exist_ok=True)
        link = dom_dir / struct_path.name
        if not link.exists():
            link.symlink_to(struct_path.resolve())
        domain_counts[dom] = domain_counts.get(dom, 0) + 1

for dom, n in sorted(domain_counts.items()):
    print(f"  {dom:<22} {n:>4} structures")

print(f"\nDomain subdirs -> {BY_DOM}")
PYEOF

# ══════════════════════════════════════════════════════════════════════
# STEP 2 — Per-domain FoldSeek: TSV + HTML + superposed PDBs
# ══════════════════════════════════════════════════════════════════════
echo ""
echo "══ STEP 2: Per-domain searches ════════════════════════════"

BY_DOMAIN_PDB="${SCRATCH}/pdb_lists/by_domain"
FMT="query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,alntmscore,qtmscore,ttmscore,rmsd,prob"

for DOM_DIR in "${BY_DOMAIN_PDB}"/*/; do
    DOM=$(basename "${DOM_DIR}")
    N=$(find "${DOM_DIR}" \( -name "*.pdb" -o -name "*.cif" \) 2>/dev/null | wc -l)

    if [[ $N -lt 2 ]]; then
        echo "  Skipping ${DOM}: only ${N} structure(s)"
        continue
    fi

    OUT="${BY_DOMAIN}/${DOM}"
    T="${TMP}/${DOM}"
    mkdir -p "${OUT}" "${T}"

    echo ""
    echo "  ── ${DOM}  (${N} proteins) ─────────────────────────"

    # TSV ─────────────────────────────────────────────────────────
    TSV="${OUT}/${DOM}_allvall.tsv"
    if [[ ! -f "${TSV}" ]]; then
        foldseek easy-search \
            "${DOM_DIR}" "${DOM_DIR}" "${TSV}" "${T}/tsv" \
            --format-output "${FMT}" \
            --alignment-type 1 \
            --tmscore-threshold 0.3 \
            --exhaustive-search 0 \
            --threads "${NCPU}" -v 2
        echo "    TSV    $(wc -l < "${TSV}") hits -> ${TSV}"
    else
        echo "    TSV already exists"
    fi

    # HTML interactive 3D viewer ──────────────────────────────────
    HTML="${OUT}/${DOM}_allvall.html"
    if [[ ! -f "${HTML}" ]]; then
        foldseek easy-search \
            "${DOM_DIR}" "${DOM_DIR}" "${HTML}" "${T}/html" \
            --format-mode 3 \
            --alignment-type 1 \
            --tmscore-threshold 0.3 \
            --exhaustive-search 0 \
            --threads "${NCPU}" -v 2
        echo "    HTML -> ${HTML}"
    else
        echo "    HTML already exists"
    fi

    # Superposed Cα PDB files (load into PyMOL / ChimeraX) ────────
    SUPER="${OUT}/${DOM}_superposed"
    mkdir -p "${SUPER}"
    if [[ $(find "${SUPER}" -name "*.pdb" 2>/dev/null | wc -l) -eq 0 ]]; then
        foldseek easy-search \
            "${DOM_DIR}" "${DOM_DIR}" \
            "${SUPER}/${DOM}_superposed.tsv" "${T}/super" \
            --format-mode 5 \
            --alignment-type 1 \
            --tmscore-threshold 0.3 \
            --exhaustive-search 0 \
            --threads "${NCPU}" -v 2
        NS=$(find "${SUPER}" -name "*.pdb" 2>/dev/null | wc -l)
        echo "    Superposed PDBs (${NS}) -> ${SUPER}"
    else
        echo "    Superposed PDBs already exist"
    fi

done

# ══════════════════════════════════════════════════════════════════════
# STEP 3 — Full all-vs-all (TSV + clustering, no HTML — too large)
# ══════════════════════════════════════════════════════════════════════
echo ""
echo "══ STEP 3: Full all-vs-all (${N_STRUCT} structures) ═══════════"

FULL="${RESULTS}/allvall_results.tsv"
if [[ ! -f "${FULL}" ]]; then
    foldseek easy-search \
        "${PDB_DIR}" "${PDB_DIR}" "${FULL}" "${TMP}/full" \
        --format-output "${FMT}" \
        --alignment-type 1 \
        --tmscore-threshold 0.3 \
        --exhaustive-search 0 \
        --threads "${NCPU}" -v 3
    echo "  Full TSV -> ${FULL}  ($(wc -l < "${FULL}") hits)"
else
    echo "  Full TSV already exists"
fi

CLUST="${RESULTS}/clusters"
if [[ ! -f "${CLUST}_cluster.tsv" ]]; then
    foldseek easy-cluster \
        "${PDB_DIR}" "${CLUST}" "${TMP}/cluster" \
        --min-seq-id 0.0 \
        --tmscore-threshold 0.4 \
        --alignment-type 1 \
        --cluster-mode 1 \
        --threads "${NCPU}" -v 3
    echo "  Clusters -> ${CLUST}_cluster.tsv"
    echo "  N clusters: $(cut -f1 "${CLUST}_cluster.tsv" | sort -u | wc -l)"
else
    echo "  Cluster TSV already exists"
fi

# ══════════════════════════════════════════════════════════════════════
# STEP 4 — Build HTML index page linking all domain viewers
# ══════════════════════════════════════════════════════════════════════
echo ""
echo "══ STEP 4: Building domain_index.html ═════════════════════"

python3 - << 'PYEOF'
from pathlib import Path

RESULTS = Path("/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/results")
BY_DOM  = RESULTS / "by_domain"

rows = []
for html_file in sorted(BY_DOM.rglob("*_allvall.html")):
    domain  = html_file.parent.name
    tsv     = html_file.parent / html_file.name.replace(".html", ".tsv")
    n_hits  = sum(1 for _ in open(tsv)) if tsv.exists() else 0
    super_d = html_file.parent / f"{domain}_superposed"
    n_super = len(list(super_d.glob("*.pdb"))) if super_d.exists() else 0
    rows.append((domain, html_file, tsv, n_hits, n_super))

table_rows = ""
for domain, hf, tsv, nh, ns in rows:
    rel_h = hf.relative_to(RESULTS)
    rel_t = tsv.relative_to(RESULTS) if tsv.exists() else "#"
    table_rows += (
        f"<tr>"
        f"<td><strong>{domain}</strong></td>"
        f"<td style='text-align:center'>{nh}</td>"
        f"<td style='text-align:center'>{ns}</td>"
        f'<td><a href="{rel_h}" target="_blank">🔬 Open 3D Viewer</a></td>'
        f'<td><a href="{rel_t}">📄 TSV</a></td>'
        f"</tr>\n"
    )

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>FoldSeek Domain Comparisons — Chromatin Reader/Writer Proteins</title>
<style>
  body  {{ font-family: 'Segoe UI', Arial, sans-serif; max-width: 920px;
           margin: 40px auto; color: #1a1a2e; background: #f8f9fa; }}
  .card {{ background: #fff; border-radius: 8px; padding: 28px 32px;
           box-shadow: 0 2px 8px rgba(0,0,0,.1); margin-bottom: 24px; }}
  h1    {{ color: #023E8A; margin-top: 0; }}
  h2    {{ color: #457B9D; border-bottom: 2px solid #e0e8f0;
           padding-bottom: 6px; margin-top: 0; }}
  table {{ border-collapse: collapse; width: 100%; }}
  th    {{ background: #023E8A; color: #fff; padding: 9px 14px;
           text-align: left; font-weight: 600; }}
  td    {{ padding: 8px 14px; border-bottom: 1px solid #e8edf2; }}
  tr:hover td {{ background: #eef4fb; }}
  a     {{ color: #023E8A; text-decoration: none; font-weight: 500; }}
  a:hover {{ text-decoration: underline; }}
  .tip  {{ background: #fff8e6; border-left: 4px solid #f4a261;
           padding: 11px 16px; border-radius: 4px; font-size: .92em; }}
  .badge {{ display: inline-block; background: #e0e8f0; color: #023E8A;
            border-radius: 12px; padding: 2px 10px; font-size: .82em;
            margin-left: 6px; font-weight: 600; }}
</style>
</head>
<body>

<div class="card">
  <h1>Chromatin Reader/Writer — Structural Similarity by Domain</h1>
  <p>FoldSeek all-vs-all comparisons per domain class using TM-align
  (alignment-type&nbsp;1, TM-score&nbsp;≥&nbsp;0.3, e-value&nbsp;≤&nbsp;0.01).
  Structures: AlphaFold&nbsp;3 (preferred) or AlphaFold&nbsp;2.
  Multi-domain proteins appear in all relevant domain groups.
  Domain boundaries determined by HMMER/Pfam-A scan.
  <span class="badge">{len(rows)} domain classes</span>
  </p>

  <div class="tip">
    <strong>How to use:</strong>
    Click <em>Open 3D Viewer</em> to launch the interactive FoldSeek HTML report
    for that domain class. The NGL viewer shows superimposed structures, TM-scores,
    RMSD, and alignment details. Works in any modern browser — no server needed.<br><br>
    Use the TSV link for downstream analysis in Python / R (05_parse_foldseek_results.py).
    Superposed Cα PDB files in each domain subfolder can be loaded directly into
    PyMOL or ChimeraX for figure generation.<br><br>
    Domain boundaries from HMMER/Pfam-A scan:
    <code>pfam/domain_boundaries.csv</code>
  </div>
</div>

<div class="card">
  <h2>Per-domain comparisons</h2>
  <table>
    <thead>
      <tr>
        <th>Domain class</th>
        <th style="text-align:center">Alignment hits</th>
        <th style="text-align:center">Superposed PDBs</th>
        <th>Interactive 3D HTML</th>
        <th>Raw TSV</th>
      </tr>
    </thead>
    <tbody>
{table_rows}
    </tbody>
  </table>
</div>

<div class="card">
  <h2>Full dataset (all proteins)</h2>
  <p>
    <a href="allvall_results.tsv">allvall_results.tsv</a>
    — all-vs-all alignment table (467 proteins)<br>
    <a href="clusters_cluster.tsv">clusters_cluster.tsv</a>
    — structural clusters (TM-score&nbsp;≥&nbsp;0.4)<br>
    <a href="figures/clustermap_tmscore.pdf">clustermap_tmscore.pdf</a>
    — hierarchical clustermap with organism / domain colour bars<br>
    <a href="figures/organism_tmscore_heatmap.pdf">organism_tmscore_heatmap.pdf</a>
    — cross-kingdom heatmap<br>
    <a href="figures/domain_tmscore_barplot.pdf">domain_tmscore_barplot.pdf</a>
    — cross-domain TM-score ranking
  </p>
</div>

<p style="font-size:.83em; color:#888; text-align:center">
  Generated by 04_run_foldseek.sh · RNASeqPaper2026 · ry00555@uga.edu
</p>
</body>
</html>
"""

idx = RESULTS / "domain_index.html"
idx.write_text(html)
print(f"Index page: {idx}")
for domain, hf, tsv, nh, ns in rows:
    print(f"  {domain:<22} hits={nh:<7} superposed={ns}")
PYEOF

rm -rf "${TMP}" 2>/dev/null || true

# ══════════════════════════════════════════════════════════════════════
# STEP 5 — Extract domain coordinates + domain-level FoldSeek
#
# Uses HMMER boundaries (pfam/domain_boundaries.csv) to extract
# per-domain coordinate windows from each structure as CA-only PDB files,
# then runs FoldSeek all-vs-all within each domain class.
# This gives true domain-level TM-scores (not whole-protein).
#
# Multi-domain proteins (e.g. EPR-1 BAH+PHD, NSD1 PWWP+PHD+SET)
# each contribute one extracted PDB per domain to the relevant group.
#
# Manually-fetched FASTAs: genes that 01_fetch_sequences.py could not
# reach automatically have per-gene .fasta files already on disk —
# these are included in HMMER scan and domain extraction automatically.
# ══════════════════════════════════════════════════════════════════════
echo ""
echo "══ STEP 5: Domain coordinate extraction + domain FoldSeek ═"

BOUNDARIES="${PFAM_DIR}/domain_boundaries.csv"
BY_DOM_EXT="${SCRATCH}/pdb_lists/by_domain_extracted"
DOM_RESULTS="${RESULTS}/by_domain_extracted"
mkdir -p "${BY_DOM_EXT}" "${DOM_RESULTS}"

if [[ ! -f "${BOUNDARIES}" ]]; then
    echo "  WARNING: domain_boundaries.csv not found — skipping domain FoldSeek"
    echo "  (HMMER step may have failed — check ${PFAM_DIR}/hmmscan.log)"
else
    # Extract domain coordinate windows from each structure
    echo "  Extracting domain PDBs from HMMER boundaries..."
    python3 - << 'PYEOF'
import csv, re
from pathlib import Path

SCRATCH   = Path("/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity")
BOUNDS    = SCRATCH / "pfam/domain_boundaries.csv"
MANIFEST  = SCRATCH / "pdb_lists/structure_manifest.tsv"
OUT_BASE  = SCRATCH / "pdb_lists/by_domain_extracted"

gene_to_struct = {}
with open(MANIFEST) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        src  = row.get("structure_file","").strip()
        gene = row.get("protein_gene_name","").strip()
        if src and src != "MISSING" and gene:
            p = Path(src)
            if p.exists():
                gene_to_struct[gene] = p

def parse_structure_ca(path):
    atoms = []
    if path.suffix == ".cif":
        col_map, cols = {}, []
        for line in path.read_text().splitlines():
            if line.startswith("_atom_site."):
                col = line.strip().lstrip("_atom_site.")
                col_map[col] = len(cols); cols.append(col)
            elif cols and line.startswith("ATOM"):
                parts = line.split()
                if len(parts) < len(cols): continue
                try:
                    if parts[col_map.get("label_atom_id",3)].strip('"') != "CA": continue
                    atoms.append((int(parts[col_map.get("label_seq_id",8)]),
                                  float(parts[col_map.get("Cartn_x",10)]),
                                  float(parts[col_map.get("Cartn_y",11)]),
                                  float(parts[col_map.get("Cartn_z",12)]),
                                  parts[col_map.get("label_comp_id",5)]))
                except (ValueError, IndexError, KeyError): continue
    else:
        for line in path.read_text().splitlines():
            if not line.startswith("ATOM"): continue
            if line[12:16].strip() != "CA": continue
            try:
                atoms.append((int(line[22:26]), float(line[30:38]),
                               float(line[38:46]), float(line[46:54]),
                               line[17:20].strip()))
            except ValueError: continue
    return atoms

def write_domain_pdb(atoms, out_path):
    lines = [f"ATOM  {i:5d}  CA  {rn:3s} A{rn_:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
             for i, (rn_, x, y, z, rn) in enumerate(atoms, 1)]
    lines.append("END")
    out_path.write_text("\n".join(lines))

written = skipped = 0
domain_counts = {}
with open(BOUNDS) as f:
    for row in csv.DictReader(f):
        gene, domain = row["gene"], row["domain"]
        env_start, env_end = int(row["env_start"]), int(row["env_end"])
        struct_path = gene_to_struct.get(gene)
        if not struct_path:
            skipped += 1; continue
        try:
            all_atoms = parse_structure_ca(struct_path)
        except Exception as e:
            print(f"  ERROR parsing {gene}: {e}"); skipped += 1; continue
        dom_atoms = [(r,x,y,z,rn) for r,x,y,z,rn in all_atoms
                     if env_start-5 <= r <= env_end+5]
        if len(dom_atoms) < 10:
            skipped += 1; continue
        out_dir = OUT_BASE / domain
        out_dir.mkdir(parents=True, exist_ok=True)
        safe = re.sub(r"[^A-Za-z0-9_-]", "_", gene)
        out_path = out_dir / f"{safe}_{domain}_{env_start}-{env_end}.pdb"
        write_domain_pdb(dom_atoms, out_path)
        written += 1
        domain_counts[domain] = domain_counts.get(domain, 0) + 1

print(f"  Written: {written}  Skipped: {skipped}")
for dom, n in sorted(domain_counts.items(), key=lambda x: -x[1]):
    print(f"    {dom:<20} {n:>4} domain PDBs")
PYEOF

    # Run FoldSeek on each domain-extracted subdirectory
    echo "  Running domain-level FoldSeek..."
    FMT="query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,alntmscore,qtmscore,ttmscore,rmsd,prob"

    for DOM_DIR in "${BY_DOM_EXT}"/*/; do
        DOM=$(basename "${DOM_DIR}")
        N=$(find "${DOM_DIR}" -name "*.pdb" 2>/dev/null | wc -l)
        if [[ $N -lt 2 ]]; then
            echo "    Skipping ${DOM}: only ${N} domain PDB(s)"
            continue
        fi
        TSV="${DOM_RESULTS}/${DOM}_domain_allvall.tsv"
        if [[ -f "${TSV}" ]]; then
            echo "    ${DOM}: TSV already exists (${N} PDBs)"
            continue
        fi
        echo "    ${DOM}: FoldSeek on ${N} domain PDBs..."
        mkdir -p "${TMP}/domext_${DOM}"
        foldseek easy-search \
            "${DOM_DIR}" "${DOM_DIR}" "${TSV}" "${TMP}/domext_${DOM}" \
            --format-output "${FMT}" \
            --alignment-type 1 \
            --tmscore-threshold 0.3 \
            --exhaustive-search 0 \
            --threads "${NCPU}" -v 1
        echo "    ${DOM}: $(wc -l < "${TSV}") hits -> ${TSV}"
    done

    echo "  Domain-level FoldSeek complete"
    echo "  Results: ${DOM_RESULTS}/"
fi

echo ""
echo "════════════════════════════════════════════════════════════"
echo " Done: $(date)"
echo ""
echo " Outputs:"
echo "   HMMER boundaries : ${PFAM_DIR}/domain_boundaries.csv"
echo "   Whole-protein TSV: ${RESULTS}/allvall_results.tsv"
echo "   Domain-level TSVs: ${DOM_RESULTS}/*_domain_allvall.tsv"
echo "   Domain 3D viewer : ${RESULTS}/domain_index.html"
echo ""
echo " Next (submitted automatically by run_reader_pipeline.sh):"
echo "   python3 05_parse_foldseek_results.py"
echo "════════════════════════════════════════════════════════════"
