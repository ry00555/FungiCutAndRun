#!/bin/bash
#SBATCH --job-name=readers_AF3
#SBATCH --partition=inter_p
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/logs/af3_%j.out
#SBATCH --error=/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity/logs/af3_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ry00555@uga.edu

# ══════════════════════════════════════════════════════════════════════
# 03b_run_alphafold3.sh
#
# Predicts structures for proteins MISSING from the AlphaFold DB.
# MSA (batch) runs as a throttled array (max MSA_THROTTLE concurrent).
# Each MSA task submits its own single INF job on completion, so INF
# jobs trickle in one at a time — never exceeding QOS limits.
#
# Usage:
#   sbatch 03b_run_alphafold3.sh          # normal run
#   sbatch 03b_run_alphafold3.sh recheck  # retry still-missing
# ══════════════════════════════════════════════════════════════════════

set -euo pipefail

SCRATCH="/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity"
MANIFEST="${SCRATCH}/pdb_lists/structure_manifest.tsv"
JSON_DIR="${SCRATCH}/pdb_lists/af3_json_inputs"
AF3_OUT="${SCRATCH}/pdb_lists/af3_outputs"
PDB_DIR="${SCRATCH}/pdb_lists/pdbs"
FASTA_DIR="${SCRATCH}/fastas"
SCRIPTS_DIR="/home/ry00555/Research/FungiCutAndRun"

MODEL_DIR="/home/ry00555/Research/AlphaFold3ModelParameters"
PUBLIC_DB="/db/AlphaFold3/20241114"
SIF="/apps/singularity-images/alphafold-3.0.1.sif"

# Max concurrent MSA jobs — keep well below QOS limit
# INF jobs submitted one-per-protein as MSA completes (never floods queue)
MSA_THROTTLE=8

mkdir -p "${JSON_DIR}" "${AF3_OUT}" "${PDB_DIR}" "${SCRATCH}/logs"

START=${1:-run}

echo "════════════════════════════════════════════════════════════"
echo " Reader Protein AF3 — ${START}"
echo " Job ID  : ${SLURM_JOB_ID}"
echo " Started : $(date)"
echo "════════════════════════════════════════════════════════════"

# ── Step 1: Generate JSON inputs for MISSING proteins ─────────────────
echo ""
echo "[Step 1] Generating JSON inputs for MISSING proteins..."

python3 - << 'PYEOF'
import csv, json, re
from pathlib import Path

SCRATCH   = Path("/scratch/ry00555/RNASeqPaper2026/Proteome/StructuralSimilarity")
MANIFEST  = SCRATCH / "pdb_lists/structure_manifest.tsv"
FASTA_DIR = SCRATCH / "fastas"
JSON_DIR  = SCRATCH / "pdb_lists/af3_json_inputs"
JSON_DIR.mkdir(parents=True, exist_ok=True)

missing = set()
with open(MANIFEST) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row.get("structure_source","").strip() == "MISSING":
            missing.add(row["protein_gene_name"].strip())
print(f"MISSING proteins: {len(missing)}")

def parse_fasta(path):
    seqs, curr = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            curr = line[1:].split("|")[0].strip()
            seqs[curr] = ""
        elif curr:
            seqs[curr] += line.strip()
    return seqs

generated = skipped = 0
for fasta_file in sorted(FASTA_DIR.glob("*.fasta")):
    if "all_readers" in fasta_file.name:
        continue
    for gene, seq in parse_fasta(fasta_file).items():
        if gene not in missing:
            continue
        seq_clean = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq.upper())
        if len(seq_clean) < 20:
            skipped += 1
            continue
        safe = re.sub(r"[^A-Za-z0-9_-]", "_", gene)
        out  = JSON_DIR / f"{safe}.json"
        if out.exists():
            skipped += 1
            continue
        out.write_text(json.dumps({
            "name": safe,
            "dialect": "alphafold3",
            "version": 1,
            "modelSeeds": [42, 137],
            "sequences": [{"protein": {"id": "A", "sequence": seq_clean}}],
        }, indent=2))
        generated += 1

print(f"Generated: {generated}  Skipped (exist/short): {skipped}")
PYEOF

N_JSON=$(ls -1 "${JSON_DIR}"/*.json 2>/dev/null | wc -l || echo 0)
echo "  ${N_JSON} JSON inputs ready"

if [[ $N_JSON -eq 0 ]]; then
    echo "No MISSING proteins to predict — exiting."
    exit 0
fi

# ── Step 2: Write the INF helper script ───────────────────────────────
# Each MSA task calls this script to submit its own INF job.
# Writing it here ensures all variables are correctly expanded.
cat > "${SCRIPTS_DIR}/03b_inf_single.sh" << INFSCRIPT
#!/bin/bash
# Called by each MSA array task to submit one INF job for a single protein.
# Usage: bash 03b_inf_single.sh <job_name>
set -euo pipefail

job_name="\$1"
SCRATCH="${SCRATCH}"
JSON_DIR="${JSON_DIR}"
AF3_OUT="${AF3_OUT}"
PDB_DIR="${PDB_DIR}"
MODEL_DIR="${MODEL_DIR}"
PUBLIC_DB="${PUBLIC_DB}"
SIF="${SIF}"

# Skip if already done
if [[ -f "\${PDB_DIR}/\${job_name}_AF3local.cif" ]]; then
    echo "INF already done, skipping: \${job_name}"
    exit 0
fi

sbatch \\
    --job-name=readers_INF \\
    --partition=gpu_p \\
    --ntasks=1 \\
    --cpus-per-task=16 \\
    --gres=gpu:A100:1 \\
    --mem=60gb \\
    --time=4:00:00 \\
    --output="\${SCRATCH}/logs/af3_inf_%j_\${job_name}.out" \\
    --error="\${SCRATCH}/logs/af3_inf_%j_\${job_name}.err" \\
    --mail-type=FAIL \\
    --mail-user=ry00555@uga.edu \\
    --wrap="
set -euo pipefail
job_name=${job_name}
MSA_JSON=${AF3_OUT}/\${job_name}/\${job_name}_data.json
echo INF: \${job_name}
[[ ! -f \${MSA_JSON} ]] && echo ERROR: MSA JSON not found: \${MSA_JSON} && exit 1
export XLA_FLAGS='--xla_disable_hlo_passes=custom-kernel-fusion-rewriter'
singularity exec \\
    --nv \\
    --bind ${JSON_DIR}:/root/af_input \\
    --bind ${AF3_OUT}:/root/af_output \\
    --bind ${MODEL_DIR}:/root/models \\
    --bind ${PUBLIC_DB}:/root/public_databases \\
    ${SIF} \\
    python /app/alphafold/run_alphafold.py \\
        --json_path=${AF3_OUT}/\${job_name}/\${job_name}_data.json \\
        --model_dir=/root/models \\
        --db_dir=/root/public_databases \\
        --output_dir=${AF3_OUT} \\
        --run_data_pipeline=false \\
        --run_inference=true
timestamped=\$(ls -dt ${AF3_OUT}/\${job_name}_*/ 2>/dev/null | head -1)
if [[ -n \"\$timestamped\" ]]; then
    mv \"\$timestamped\"/* ${AF3_OUT}/\${job_name}/
    rmdir \"\$timestamped\"
fi
for cif in \$(find ${AF3_OUT}/\${job_name} -name '*_model.cif' 2>/dev/null); do
    dest=${PDB_DIR}/\${job_name}_AF3local.cif
    [[ ! -f \"\$dest\" ]] && cp \"\$cif\" \"\$dest\" && echo Copied: \$(basename \$dest)
done
echo INF done: \${job_name}
"
echo "INF submitted for \${job_name}"
INFSCRIPT

chmod +x "${SCRIPTS_DIR}/03b_inf_single.sh"
echo "  INF helper script written: ${SCRIPTS_DIR}/03b_inf_single.sh"

# ── Step 3: Submit MSA array ───────────────────────────────────────────
# Each task runs MSA then calls 03b_inf_single.sh to submit its INF job.
echo ""
echo "[Step 3] Submitting MSA array (batch, throttle=${MSA_THROTTLE})..."

MSA_JOBID=$(sbatch --parsable \
    --job-name=readers_MSA \
    --partition=batch \
    --ntasks=1 \
    --cpus-per-task=32 \
    --mem=150gb \
    --time=8:00:00 \
    --output="${SCRATCH}/logs/af3_msa_%j_%a.out" \
    --error="${SCRATCH}/logs/af3_msa_%j_%a.err" \
    --mail-type=FAIL \
    --mail-user=ry00555@uga.edu \
    --array=1-${N_JSON}%${MSA_THROTTLE} \
    --wrap="
set -euo pipefail
file=\$(ls ${JSON_DIR}/*.json | awk \"NR==\${SLURM_ARRAY_TASK_ID}\")
job_name=\$(basename \"\$file\" .json | tr '[:upper:]' '[:lower:]')
echo \"MSA: \${job_name}  [task \${SLURM_ARRAY_TASK_ID}]\"

if [[ -f \"${AF3_OUT}/\${job_name}/\${job_name}_data.json\" ]]; then
    echo \"MSA already done, skipping: \${job_name}\"
else
    singularity exec \\
        --bind ${JSON_DIR}:/root/af_input \\
        --bind ${AF3_OUT}:/root/af_output \\
        --bind ${MODEL_DIR}:/root/models \\
        --bind ${PUBLIC_DB}:/root/public_databases \\
        ${SIF} \\
        python /app/alphafold/run_alphafold.py \\
            --json_path=/root/af_input/\$(basename \$file) \\
            --model_dir=/root/models \\
            --db_dir=/root/public_databases \\
            --output_dir=/root/af_output \\
            --run_data_pipeline=true \\
            --run_inference=false
    echo \"MSA done: \${job_name}\"
fi

# Submit INF for this protein via helper script
bash ${SCRIPTS_DIR}/03b_inf_single.sh \"\${job_name}\"
")

[[ -z "$MSA_JOBID" ]] && echo "MSA submission failed." && exit 1
echo "  MSA array: ${MSA_JOBID}  (1-${N_JSON}, max ${MSA_THROTTLE} concurrent)"
echo "  INF jobs will be submitted one-per-protein as MSA tasks complete"

echo ""
echo "════════════════════════════════════════════════════════════"
echo " Submitted:"
echo "   MSA array : ${MSA_JOBID}  (batch, 1-${N_JSON}, max ${MSA_THROTTLE} concurrent)"
echo "   INF jobs  : 1 per protein, submitted as each MSA task completes"
echo ""
echo " Monitor:  squeue -u ry00555"
echo " After all INF complete:"
echo "   python3 03_fetch_structures.py  # update manifest"
echo "   sbatch  04_run_foldseek.sh      # re-run to include new structs"
echo "════════════════════════════════════════════════════════════"
