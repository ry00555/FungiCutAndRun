#!/bin/bash
#SBATCH --job-name=EAF3
#SBATCH --partition=inter_p
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --time=01:00:00
#SBATCH --output=AF3_Pool.%j.out
#SBATCH --error=AF3_Pool.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu

# ── Settings ──────────────────────────────────────────────────────────────────
TOTAL=173         # update to your total number of JSONs
BATCH_SIZE=5      # number of pools per batch
MSA_SCRIPT="AF3_MSA_Pool.sh"
INF_SCRIPT="AF3_INF_Pool.sh"

START=${1:-1}
END=$((START + BATCH_SIZE - 1))
[ $END -gt $TOTAL ] && END=$TOTAL

echo "Submitting batch: pools ${START}-${END}"

# ── Submit MSA jobs for this batch ────────────────────────────────────────────
MSA_JOBID=$(sbatch --parsable --array=${START}-${END} $MSA_SCRIPT)

if [ -z "$MSA_JOBID" ]; then
    echo "MSA submission failed. Exiting."
    exit 1
fi
echo "MSA job submitted: $MSA_JOBID (pools ${START}-${END})"

# ── Submit INF jobs dependent on MSA finishing ────────────────────────────────
INF_JOBID=$(sbatch --parsable \
    --dependency=afterok:${MSA_JOBID} \
    --array=${START}-${END} \
    $INF_SCRIPT)

if [ -z "$INF_JOBID" ]; then
    echo "INF submission failed. Exiting."
    exit 1
fi
echo "INF job submitted: $INF_JOBID (depends on MSA $MSA_JOBID)"

# ── Submit next batch dependent on MSA finishing ──────────────────────────────
NEXT_START=$((END + 1))
if [ $NEXT_START -le $TOTAL ]; then
    echo "Queuing next batch starting at pool $NEXT_START (depends on MSA $MSA_JOBID)" #this was changed from INF_JOBID so if you want MSA 1-5 → INF 1-5 → THEN MSA 6-10 then switch to INF_JOBID
    sbatch --dependency=afterok:${MSA_JOBID} AF3_Pools.sh $NEXT_START
else
    echo "All $TOTAL pools submitted."
fi
