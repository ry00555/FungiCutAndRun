#!/bin/bash
#SBATCH --job-name=PRC2_batcher
#SBATCH --partition=inter_p
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/ry00555/Boltz2/logs/boltz2_batch.%j.out
#SBATCH --error=/scratch/ry00555/Boltz2/logs/boltz2_batch.%j.err

cd $SLURM_SUBMIT_DIR

# ── Only change these lines when swapping proteins ────────────────────────────
PROJECT="PRC2"
TOTAL=3261       # total YAMLs minus 1 (0-based) — update after running generator
BATCH_SIZE=15
MAX_RUNNING=8
# ─────────────────────────────────────────────────────────────────────────────

SCRIPT="/home/ry00555/Research/FungiCutAndRun/Boltz2_Step1.sh"

START=${1:-0}
END=$((START + BATCH_SIZE - 1))
[ $END -gt $TOTAL ] && END=$TOTAL

echo "Submitting Boltz2 $PROJECT array ${START}-${END}%${MAX_RUNNING}"

JOBID=$(sbatch --parsable --array=${START}-${END}%${MAX_RUNNING} $SCRIPT)
[ -z "$JOBID" ] && echo "Submission failed." && exit 1
echo "Submitted job $JOBID"

NEXT_START=$((END + 1))
if [ $NEXT_START -le $TOTAL ]; then
    sbatch --dependency=afterany:${JOBID} Batches_Boltz2.sh $NEXT_START
else
    echo "All $PROJECT Boltz2 batches submitted."
fi
