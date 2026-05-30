#!/bin/bash
#SBATCH --job-name=PRC2
#SBATCH --partition=inter_p
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --time=01:00:00
#SBATCH --output=AF3_Pool.%j.out
#SBATCH --error=AF3_Pool.%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ry00555@uga.edu

# ── Settings ──────────────────────────────────────────────────────────────────
TOTAL=3366
BATCH_SIZE=9
MSA_SCRIPT="PRC2_AF3_MSA_Pool.sh"
INF_SCRIPT="PRC2_AF3_INF_Pool.sh"
INPUT_DIR="/scratch/ry00555/RNASeqPaper2026/Proteome/PRC2_Proteome_pools/PooledPPI/PRC2_StringDB_pools/AF3_JSONs"
OUTPUT_DIR="/scratch/ry00555/RNASeqPaper2026/Proteome/PRC2_Proteome_pools/PooledPPI/PRC2_AF3_PooledJSON_output"
# ─────────────────────────────────────────────────────────────────────────────

START=${1:-1}

# ── If START is "recheck" run the missing pools check ─────────────────────────
if [ "$START" == "recheck" ]; then
    echo "Running missing pools check..."
    missing=()
    for i in $(seq 1 $TOTAL); do
        pool=$(printf "prc2_pool%03d" $i)
        dir="${OUTPUT_DIR}/${pool}"
        if ! ls $dir/seed-1_sample-0/summary_confidences.json &>/dev/null 2>&1; then
            missing+=($i)
        fi
    done

    if [ ${#missing[@]} -eq 0 ]; then
        echo "All pools complete!"
        exit 0
    fi

    echo "Found ${#missing[@]} missing pools: ${missing[@]}"

    # Submit missing in batches of 10
    CHUNK=200
   PREV_JOBID=""
   for ((start=0; start<${#missing[@]}; start+=CHUNK)); do
       chunk=("${missing[@]:$start:$CHUNK}")
       indices=$(IFS=,; echo "${chunk[*]}")

       if [ -z "$PREV_JOBID" ]; then
           MSA_JOBID=$(sbatch --parsable --array=${indices}%5 $MSA_SCRIPT)
       else
           MSA_JOBID=$(sbatch --parsable \
               --dependency=afterany:${PREV_JOBID} \
               --array=${indices}%5 $MSA_SCRIPT)
       fi
       [ -z "$MSA_JOBID" ] && echo "MSA submission failed at chunk $start." && exit 1
       echo "MSA chunk submitted: $MSA_JOBID (${#chunk[@]} pools)"

       INF_JOBID=$(sbatch --parsable \
           --dependency=afterok:${MSA_JOBID} \
           --array=${indices}%5 \
           $INF_SCRIPT)
       echo "INF chunk submitted: $INF_JOBID"

       PREV_JOBID=$INF_JOBID
   done

   sbatch --dependency=afterany:${PREV_JOBID} PRC2_AF3_Pools.sh recheck
    exit 0
fi

# ── Normal batch submission ───────────────────────────────────────────────────
END=$((START + BATCH_SIZE - 1))
[ $END -gt $TOTAL ] && END=$TOTAL

echo "Submitting batch: pools ${START}-${END}"

MSA_JOBID=$(sbatch --parsable --array=${START}-${END} $MSA_SCRIPT)
[ -z "$MSA_JOBID" ] && echo "MSA submission failed." && exit 1
echo "MSA job submitted: $MSA_JOBID (pools ${START}-${END})"

INF_JOBID=$(sbatch --parsable \
    --dependency=afterok:${MSA_JOBID} \
    --array=${START}-${END} \
    $INF_SCRIPT)
[ -z "$INF_JOBID" ] && echo "INF submission failed." && exit 1
echo "INF job submitted: $INF_JOBID"

NEXT_START=$((END + 1))
if [ $NEXT_START -le $TOTAL ]; then
    sbatch --dependency=afterany:${MSA_JOBID} PRC2_AF3_Pools.sh $NEXT_START
else
    # All batches submitted — schedule automatic recheck
    echo "All batches submitted — scheduling recheck for missing pools"
    sbatch --dependency=afterany:${INF_JOBID} PRC2_AF3_Pools.sh recheck
fi
