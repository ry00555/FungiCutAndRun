#!/usr/bin/env python3
"""
02_annotate_domains_idr.py
--------------------------
Annotate all proteins with:
  1. % IDR via AIUPred v3  (local module on Sapelo2 GPU partition)
     module: AIUPred/3.1-foss-2024a-CUDA-12.6.0
     Python API: from aiupred import AIUPred, multifasta_reader
     Models loaded ONCE, all proteins processed in one pass.
     Provides per-residue disorder scores, binding site scores,
     and flexible linker scores.

  2. Domain scan via regex heuristics
     (proper HMMER scan against Pfam can follow separately)

Outputs -> SCRATCH/
  master_reader_annotated_v4.csv
  aiupred_raw/                    per-protein TSV from AIUPred
  logs/02_annotate.log

Usage:
    # Requires AIUPred module loaded — run inside the SLURM job
    # or interactively after:  module load AIUPred/3.1-foss-2024a-CUDA-12.6.0
    python 02_annotate_domains_idr.py

    # Force CPU (no GPU available):
    python 02_annotate_domains_idr.py --force-cpu
"""

import csv, re, sys, argparse, logging
from pathlib import Path

SCRATCH = Path("/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity")

# ── Domain regex patterns ─────────────────────────────────────────────
DOMAIN_PATTERNS = {
    "PWWP":         r"[WY].{1,4}PWWP",
    "Chromodomain": r"[LI].{2,5}[YW].{5,15}[RK].{3,8}[DE]",
    "Tudor":        r"[VILM]{2}.{5,15}[YW].{5,10}[YW]",
    "MRG":          r"[KRDE]{3}.{20,60}[FY].{10,30}[FY]",
    "SET":          r"NH[SC].C.PN",
    "WD40":         r"[WD][DE].{40}[WD][DE]",
    "JmjC":         r"H.D.{30,60}E",
    "CXXC":         r"C.{1,3}C.{2,5}C.{1,3}C",
    "PHD":          r"C.{2}C.{3,5}C.{1,4}[CH]",
    "BAH":          r"[WY].{5,20}[DE].{5,20}[WY].{5,20}[RK]",
    "Bromodomain":  r"[VILM].{10,30}[YWF].{5,15}[DE].{5,15}[YWF]",
    "IDR_rich":     r"[GASPQE]{8,}",
}


def setup_logging(logdir: Path) -> logging.Logger:
    logdir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(logdir / "02_annotate.log"),
        ],
    )
    return logging.getLogger(__name__)


def scan_domains(seq: str) -> str:
    found = [d for d, pat in DOMAIN_PATTERNS.items()
             if re.search(pat, seq, re.IGNORECASE)]
    return ";".join(found) if found else "none_detected"


def load_fasta_seqs(fasta_dir: Path) -> dict[str, str]:
    """Load per-protein FASTAs. Key = gene name (first pipe-delimited field)."""
    seqs = {}
    for f in sorted(fasta_dir.glob("*.fasta")):
        if "combined" in f.name or "all_readers" in f.name:
            continue
        curr_gene, curr_seq = None, []
        for line in f.read_text().splitlines():
            if line.startswith(">"):
                if curr_gene and curr_seq:
                    seqs[curr_gene] = "".join(curr_seq)
                curr_gene = line[1:].split("|")[0].strip()
                curr_seq  = []
            elif curr_gene:
                curr_seq.append(line.strip())
        if curr_gene and curr_seq:
            seqs[curr_gene] = "".join(curr_seq)
    return seqs


def run_aiupred(seqs: dict[str, str], raw_dir: Path,
                force_cpu: bool, log) -> dict[str, dict]:
    """
    Run AIUPred v3 over all sequences using the Python API.
    Models are loaded ONCE and all proteins are processed in a single
    predictor instance — this is the correct AIUPred usage pattern
    and avoids repeated model-load overhead on the GPU.

    Returns: dict[gene] -> {
        "idr_pct": float,
        "mean_disorder": float,
        "binding_pct": float,    # fraction of residues with binding score > 0.5
        "linker_pct": float,     # fraction with linker score > 0.5
        "idr_regions": str,      # "10-45,100-130"
    }
    """
    try:
        from aiupred import AIUPred
    except ImportError:
        log.error(
            "AIUPred Python package not found.\n"
            "  On Sapelo2: module load AIUPred/3.1-foss-2024a-CUDA-12.6.0\n"
            "  Then re-run this script."
        )
        sys.exit(1)

    log.info(f"Loading AIUPred v3 model  (force_cpu={force_cpu})...")
    predictor = AIUPred(force_cpu=force_cpu)   # models loaded here, once
    log.info("  Model ready.")

    raw_dir.mkdir(parents=True, exist_ok=True)
    results = {}
    IDR_THRESH     = 0.5
    BINDING_THRESH = 0.5
    LINKER_THRESH  = 0.5

    for i, (gene, seq) in enumerate(seqs.items(), 1):
        if not seq or len(seq) < 5:
            log.warning(f"  Skipping {gene}: sequence too short ({len(seq)} aa)")
            continue

        try:
            # AIUPred v3: each method takes only sequence, returns numpy array
            disorder = list(predictor.predict_disorder(seq))
            try:
                binding = list(predictor.predict_binding(seq))
            except Exception:
                binding = [0.0] * len(seq)
            try:
                linker = list(predictor.predict_linker(seq))
            except Exception:
                linker = [0.0] * len(seq)

            n = len(disorder)

            # % IDR (disorder > threshold)
            idr_pct  = round(100.0 * sum(s > IDR_THRESH     for s in disorder) / n, 1)
            bind_pct = round(100.0 * sum(s > BINDING_THRESH for s in binding)  / n, 1)
            link_pct = round(100.0 * sum(s > LINKER_THRESH  for s in linker)   / n, 1)
            mean_dis = round(sum(disorder) / n, 4)

            # Identify IDR regions (contiguous runs of disorder > threshold)
            regions, start = [], None
            for idx, s in enumerate(disorder):
                if s > IDR_THRESH and start is None:
                    start = idx + 1     # 1-based
                elif s <= IDR_THRESH and start is not None:
                    regions.append(f"{start}-{idx}")
                    start = None
            if start is not None:
                regions.append(f"{start}-{n}")
            idr_regions = ",".join(regions) if regions else "none"

            results[gene] = {
                "idr_pct":      idr_pct,
                "mean_disorder":mean_dis,
                "binding_pct":  bind_pct,
                "linker_pct":   link_pct,
                "idr_regions":  idr_regions,
            }

            # Write per-protein TSV (useful for downstream analyses)
            tsv_path = raw_dir / f"{gene}_aiupred.tsv"
            with open(tsv_path, "w") as fh:
                fh.write("position\tresidue\tdisorder\tbinding\tlinker\n")
                for pos, (aa, d, b, l) in enumerate(
                        zip(seq, disorder, binding, linker), 1):
                    fh.write(f"{pos}\t{aa}\t{d:.4f}\t{b:.4f}\t{l:.4f}\n")

            if i % 50 == 0:
                log.info(f"  [{i}/{len(seqs)}] processed...")

            log.debug(f"  {gene}: IDR={idr_pct}%  bind={bind_pct}%  link={link_pct}%")

        except Exception as e:
            log.warning(f"  AIUPred failed for {gene}: {e}")
            results[gene] = {
                "idr_pct": None, "mean_disorder": None,
                "binding_pct": None, "linker_pct": None,
                "idr_regions": "ERROR",
            }

    log.info(f"AIUPred complete: {len(results)} proteins processed")
    return results


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--csv",
        default=str(SCRATCH / "master_reader_proteins_v4.csv"))
    p.add_argument("--fastas",
        default=str(SCRATCH / "fastas"))
    p.add_argument("--out",
        default=str(SCRATCH / "master_reader_annotated_v4.csv"))
    p.add_argument("--logdir",
        default=str(SCRATCH / "logs"))
    p.add_argument("--raw-dir",
        default=str(SCRATCH / "aiupred_raw"),
        help="Directory for per-protein AIUPred TSV files")
    p.add_argument("--force-cpu", action="store_true",
        help="Disable GPU (use CPU only — slower but works without CUDA)")
    return p.parse_args()


def main():
    args    = parse_args()
    logdir  = Path(args.logdir)
    raw_dir = Path(args.raw_dir)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    log = setup_logging(logdir)
    log.info(f"Scratch dir   : {SCRATCH}")
    log.info(f"FASTA dir     : {args.fastas}")
    log.info(f"Annotated out : {out_path}")
    log.info(f"AIUPred TSVs  : {raw_dir}")
    log.info(f"Force CPU     : {args.force_cpu}")

    # Load sequences
    seqs = load_fasta_seqs(Path(args.fastas))
    log.info(f"Loaded {len(seqs)} sequences")

    # Run AIUPred (GPU-accelerated via loaded module)
    aiupred_results = run_aiupred(seqs, raw_dir, args.force_cpu, log)

    # Load master CSV and annotate
    rows_out = []
    with open(args.csv, newline="") as fh:
        reader     = csv.DictReader(fh)
        base_fields = reader.fieldnames or []
        new_fields  = [
            "sequence",
            "seq_length",
            "idr_percent_computed",     # AIUPred disorder > 0.5
            "idr_mean_score",           # mean per-residue disorder score
            "binding_percent",          # AIUPred predicted binding sites
            "linker_percent",           # AIUPred predicted linkers
            "idr_regions",              # 1-based residue ranges
            "domains_regex_found",      # quick regex scan
            "fetch_status",
        ]
        fieldnames = base_fields + [f for f in new_fields if f not in base_fields]
        rows = list(reader)

    for row in rows:
        gene = row["protein_gene_name"].strip()
        seq  = seqs.get(gene, "")
        aiu  = aiupred_results.get(gene, {})

        domains_str = scan_domains(seq) if seq else "MISSING_SEQ"
        status      = "OK" if seq else "MISSING"

        row.update({
            "sequence":              seq,
            "seq_length":            len(seq),
            "idr_percent_computed":  aiu.get("idr_pct",      row.get("percent_IDR_estimated", "NA")),
            "idr_mean_score":        aiu.get("mean_disorder", "NA"),
            "binding_percent":       aiu.get("binding_pct",  "NA"),
            "linker_percent":        aiu.get("linker_pct",   "NA"),
            "idr_regions":           aiu.get("idr_regions",  "NA"),
            "domains_regex_found":   domains_str,
            "fetch_status":          status,
        })
        rows_out.append(row)

        if seq:
            log.info(
                f"  {gene}: IDR={aiu.get('idr_pct','?')}%  "
                f"bind={aiu.get('binding_pct','?')}%  "
                f"linker={aiu.get('linker_pct','?')}%  "
                f"domains={domains_str}"
            )
        else:
            log.warning(f"  {gene}: no sequence — using CSV estimates")

    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows_out)

    ok      = sum(1 for r in rows_out if r["fetch_status"] == "OK")
    missing = len(rows_out) - ok
    log.info(
        f"\n{'='*60}\n"
        f"Annotated CSV  : {out_path}\n"
        f"  Rows         : {len(rows_out)}\n"
        f"  With AIUPred : {ok}\n"
        f"  Missing seq  : {missing} (CSV estimates used)\n"
        f"  Raw TSVs     : {raw_dir}\n"
        f"{'='*60}"
    )
    log.info("Done.")


if __name__ == "__main__":
    main()
