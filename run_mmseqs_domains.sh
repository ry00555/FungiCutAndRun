#!/bin/bash
#SBATCH --job-name=seq_identity
#SBATCH --partition=inter_p
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=20:00:00
#SBATCH --output=/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek/chromatin_domain_results/logs/seq_identity_%j.out
#SBATCH --error=/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek/chromatin_domain_results/logs/seq_identity_%j.err

ml MMseqs2/18-8cc5c-gompi-2025a

# ── Paths ──────────────────────────────────────────────────────────────
BASE="/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek/chromatin_domain_results"
PULL_OUT="/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/results/proteome_domain_pull"
MASTER_CSV="${PULL_OUT}/domain_proteins_master.csv"   # from 07a
OUT="${BASE}/mmseqs_identity"
DOM_FASTA_DIR="${BASE}/domain_fastas"
PROT_FASTA_DIR="${BASE}/protein_fastas"
TMP="${BASE}/tmp_$$"

mkdir -p "${OUT}" "${DOM_FASTA_DIR}" "${PROT_FASTA_DIR}" "${TMP}" \
         "${OUT}/domain_results" "${OUT}/protein_results" \
         "${BASE}/logs"

echo "════════════════════════════════════════════════════════════"
echo " MMseqs2 Domain + Whole-Protein Sequence Identity"
echo " Job ID  : ${SLURM_JOB_ID}"
echo " Started : $(date)"
echo "════════════════════════════════════════════════════════════"

# ── Step 1: Build per-domain FASTAs from UniProt pull ─────────────────
# Uses domain_proteins_master.csv (from 07a_uniprot_domain_pull.py)
# Each protein gets:
#   - one entry in its primary domain FASTA (full protein sequence)
#   - also added to any co-occurring domain FASTA (multi-domain aware)
echo ""
echo "[Step 1] Building per-domain protein FASTAs from UniProt..."

python3 << 'PYEOF'
import csv, re, os, time
from pathlib import Path
from collections import defaultdict
import requests

PULL_OUT      = Path("/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/results/proteome_domain_pull")
MASTER_CSV    = PULL_OUT / "domain_proteins_master.csv"
PROT_FASTA    = Path("/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek/chromatin_domain_results/protein_fastas")
PROT_FASTA.mkdir(parents=True, exist_ok=True)

# Load master CSV
with open(MASTER_CSV, newline="") as fh:
    rows = list(csv.DictReader(fh))
print(f"Loaded {len(rows)} rows from master CSV")

# Group accessions by domain — multi-domain proteins go into all their domains
domain_to_accs = defaultdict(dict)   # domain -> {acc: row}
for row in rows:
    acc  = row["uniprot_acc"]
    pdom = row["query_domain"]
    domain_to_accs[pdom][acc] = row
    # Also add to co-occurring domains
    for codom in row.get("cooccurring_chromatin_domains","").split(";"):
        codom = codom.strip()
        if codom and codom != pdom:
            domain_to_accs[codom][acc] = row

# Fetch sequences from UniProt for all unique accessions
unique_accs = {row["uniprot_acc"] for row in rows}
print(f"Fetching sequences for {len(unique_accs)} unique accessions...")

session = requests.Session()
session.headers["User-Agent"] = "MMseqsDomainFetcher/1.0 (ry00555@uga.edu)"
FASTA_URL = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"

sequences = {}
failed = []
for i, acc in enumerate(sorted(unique_accs), 1):
    try:
        r = session.get(FASTA_URL.format(acc=acc), timeout=20)
        if r.status_code == 200 and r.text.startswith(">"):
            # Extract sequence only (skip header, will rewrite)
            lines = r.text.strip().splitlines()
            seq   = "".join(l for l in lines[1:] if not l.startswith(">"))
            if seq:
                sequences[acc] = seq
        else:
            failed.append(acc)
        time.sleep(0.15)
        if i % 200 == 0:
            print(f"  [{i}/{len(unique_accs)}] fetched={len(sequences)} failed={len(failed)}")
    except Exception as e:
        failed.append(acc)

print(f"Sequences fetched: {len(sequences)}  Failed: {len(failed)}")

# Write per-domain protein FASTA files
# Header format: >GENE|ACC|ORGANISM|DOMAIN|KINGDOM
for domain, acc_dict in sorted(domain_to_accs.items()):
    out_path = PROT_FASTA / f"{domain}_proteins.fasta"
    written  = 0
    with open(out_path, "w") as fh:
        for acc, row in sorted(acc_dict.items()):
            seq = sequences.get(acc)
            if not seq:
                continue
            gene    = row.get("gene_name", acc) or acc
            org     = row.get("organism", "").replace(" ", "_")
            kingdom = row.get("kingdom", "")
            all_dom = row.get("all_pfam_domains", "")
            header  = f">{gene}|{acc}|{org}|{domain}|{kingdom}|{all_dom[:60]}"
            fh.write(f"{header}\n{seq}\n")
            written += 1
    print(f"  {domain}: {written} sequences -> {out_path.name}")

# Write all-domains combined FASTA
all_path = PROT_FASTA / "all_chromatin_proteins.fasta"
seen_accs = set()
with open(all_path, "w") as fh:
    for row in rows:
        acc = row["uniprot_acc"]
        if acc in seen_accs:
            continue
        seen_accs.add(acc)
        seq = sequences.get(acc)
        if not seq:
            continue
        gene    = row.get("gene_name", acc) or acc
        org     = row.get("organism","").replace(" ","_")
        dom     = row.get("query_domain","")
        kingdom = row.get("kingdom","")
        codoms  = row.get("cooccurring_chromatin_domains","")
        all_dom = row.get("all_pfam_domains","")
        header  = f">{gene}|{acc}|{org}|{dom}|{kingdom}|codoms={codoms[:40]}|pfam={all_dom[:60]}"
        fh.write(f"{header}\n{seq}\n")
print(f"\nAll-domains combined: {all_path.name} ({len(seen_accs)} proteins)")
print("Done.")
PYEOF

echo ""
echo "[Step 2] Running MMseqs2 domain-level all-vs-all..."

for fasta in "${PROT_FASTA_DIR}"/*_proteins.fasta; do
    [[ ! -f "$fasta" ]] && continue
    domain=$(basename "$fasta" _proteins.fasta)
    N=$(grep -c "^>" "$fasta" 2>/dev/null || echo 0)
    [[ $N -lt 2 ]] && echo "  [SKIP] ${domain}: only ${N} sequence(s)" && continue

    DB="${OUT}/protein_results/${domain}_db"
    RESULT="${OUT}/protein_results/${domain}_identity.tsv"

    echo "  ${domain}: ${N} proteins"

    # Build MMseqs database — FASTA input, protein mode
    mmseqs createdb "${fasta}" "${DB}" \
        --dbtype 1 \
        --threads ${SLURM_CPUS_PER_TASK} \
        2>/dev/null

    # All-vs-all protein identity
    # --search-type 1 = amino acid (prevents nucleotide mode error)
    mmseqs easy-search \
        "${fasta}" \
        "${fasta}" \
        "${RESULT}" \
        "${TMP}/tmp_${domain}" \
        --search-type 1 \
        -e 1e-5 \
        -c 0.3 \
        --cov-mode 0 \
        --format-output "query,target,pident,alnlen,qcov,tcov,evalue,bits" \
        --threads ${SLURM_CPUS_PER_TASK} \
        2>/dev/null

    N_HITS=$(wc -l < "${RESULT}" 2>/dev/null || echo 0)
    echo "    -> ${N_HITS} hits"
done

echo ""
echo "[Step 3] All-chromatin-proteins combined identity..."

ALL_FASTA="${PROT_FASTA_DIR}/all_chromatin_proteins.fasta"
ALL_RESULT="${OUT}/all_protein_sequence_identity.tsv"
N_ALL=$(grep -c "^>" "${ALL_FASTA}" 2>/dev/null || echo 0)
echo "  Total proteins: ${N_ALL}"

mmseqs easy-search \
    "${ALL_FASTA}" \
    "${ALL_FASTA}" \
    "${ALL_RESULT}" \
    "${TMP}/tmp_all" \
    --search-type 1 \
    -e 1e-5 \
    -c 0.3 \
    --cov-mode 0 \
    --format-output "query,target,pident,alnlen,qcov,tcov,evalue,bits" \
    --threads ${SLURM_CPUS_PER_TASK} \
    2>/dev/null

echo "  All-vs-all hits: $(wc -l < ${ALL_RESULT})"

echo ""
echo "[Step 4] Combining results and merging with FoldSeek master..."

python3 << 'PYEOF'
import csv, os
from pathlib import Path
from collections import defaultdict

BASE        = Path("/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/FoldSeek/chromatin_domain_results")
PULL_OUT    = Path("/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/results/proteome_domain_pull")
OUT         = BASE / "mmseqs_identity"
PROT_RESULT = OUT / "protein_results"
MASTER_CSV  = PULL_OUT / "domain_proteins_master.csv"
FOLDSEEK_CSV= PULL_OUT / "foldseek_results/annotated_hits.csv"

# ── Parse gene|acc|org headers back to accession ─────────────────────
def header_to_acc(h: str) -> str:
    """Extract UniProt acc from >GENE|ACC|ORG|... header"""
    parts = h.lstrip(">").split("|")
    return parts[1] if len(parts) > 1 else parts[0]

# ── Combine per-domain identity TSVs ─────────────────────────────────
COMBINED = OUT / "all_protein_sequence_identity.tsv"
COLS = ["query_acc","target_acc","pident","alnlen","qcov","tcov","evalue","bits","domain"]

rows_out = []
for tsv in sorted(PROT_RESULT.glob("*_identity.tsv")):
    domain = tsv.stem.replace("_identity","")
    with open(tsv) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 8: continue
            q_acc = header_to_acc(parts[0])
            t_acc = header_to_acc(parts[1])
            rows_out.append({
                "query_acc":  q_acc,
                "target_acc": t_acc,
                "pident":     parts[2],
                "alnlen":     parts[3],
                "qcov":       parts[4],
                "tcov":       parts[5],
                "evalue":     parts[6],
                "bits":       parts[7],
                "domain":     domain,
            })

with open(COMBINED, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t")
    w.writeheader()
    w.writerows(rows_out)
print(f"Combined identity: {COMBINED}  ({len(rows_out)} rows)")

# ── Load master CSV for metadata ──────────────────────────────────────
meta = {}
with open(MASTER_CSV, newline="") as fh:
    for row in csv.DictReader(fh):
        meta[row["uniprot_acc"]] = row

# ── Merge MMseqs with FoldSeek results if available ───────────────────
if FOLDSEEK_CSV.exists():
    # Build mmseqs lookup: (q_acc, t_acc) -> best pident
    mmseqs_idx = {}
    for row in rows_out:
        key = (row["query_acc"], row["target_acc"])
        try:
            pi = float(row["pident"])
            if key not in mmseqs_idx or pi > float(mmseqs_idx[key]["pident"]):
                mmseqs_idx[key] = row
        except ValueError:
            pass

    # Load FoldSeek hits and join
    merged = []
    with open(FOLDSEEK_CSV, newline="") as fh:
        fs_rows = list(csv.DictReader(fh))

    for fs in fs_rows:
        q = fs.get("query_acc","")
        t = fs.get("target_acc","")
        mm = mmseqs_idx.get((q,t), mmseqs_idx.get((t,q), {}))
        merged.append({
            **fs,
            "mmseqs_pident":  mm.get("pident",""),
            "mmseqs_alnlen":  mm.get("alnlen",""),
            "mmseqs_qcov":    mm.get("qcov",""),
            "mmseqs_evalue":  mm.get("evalue",""),
            "seq_vs_struct":  "both" if mm else "struct_only",
        })

    merged_path = PULL_OUT / "foldseek_results/annotated_hits_with_seqid.csv"
    if merged:
        with open(merged_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(merged[0].keys()))
            w.writeheader()
            w.writerows(merged)
        print(f"Merged FoldSeek+MMseqs: {merged_path}  ({len(merged)} rows)")

        # Summary: structural similarity vs sequence identity
        paired = [r for r in merged if r["mmseqs_pident"]]
        print(f"\nPaired (both structural + sequence data): {len(paired)}")
        print("Cross-kingdom hits where struct TM>0.5 but seq_id<30%:")
        divergent = [r for r in paired
                     if float(r.get("mean_tmscore",0) or 0) > 0.5
                     and float(r.get("mmseqs_pident","100") or 100) < 30
                     and r.get("cross_kingdom","") == "yes"]
        print(f"  {len(divergent)} hits — structurally conserved but sequence-divergent")
        print("  (These are your most interesting cross-kingdom comparisons)")
        for r in sorted(divergent, key=lambda x: -float(x.get("mean_tmscore",0) or 0))[:15]:
            print(f"  {r.get('query_gene','?'):<12} vs {r.get('target_gene','?'):<12}"
                  f"  TM={r.get('mean_tmscore','?')}  seqid={r.get('mmseqs_pident','?')}%"
                  f"  {r.get('query_organism','').split()[0]} vs {r.get('target_organism','').split()[0]}")
else:
    print(f"FoldSeek annotated_hits.csv not found yet — run 07c first")
    print(f"MMseqs results saved to: {COMBINED}")

print("\nDone.")
PYEOF

rm -rf "${TMP}"

echo ""
echo "════════════════════════════════════════════════════════════"
echo " Complete: $(date)"
echo " Outputs:"
echo "   Per-domain identity : ${OUT}/protein_results/*_identity.tsv"
echo "   Combined identity   : ${OUT}/all_protein_sequence_identity.tsv"
echo "   Merged (if 07c done): foldseek_results/annotated_hits_with_seqid.csv"
echo "════════════════════════════════════════════════════════════"
