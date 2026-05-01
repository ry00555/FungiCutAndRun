#!/bin/bash

TOTAL=584
BATCH_SIZE=5
MAX_RUNNING=1
SCRIPT="EAF3_AF3_Pool.sh"

START=${1:-1}
END=$((START + BATCH_SIZE - 1))

if [ $END -gt $TOTAL ]; then
    END=$TOTAL
fi

echo "Submitting EAF3 AF3 array ${START}-${END}%${MAX_RUNNING}"

JOBID=$(sbatch --parsable --array=${START}-${END}%${MAX_RUNNING} $SCRIPT)

echo "Submitted job $JOBID"

NEXT_START=$((END + 1))

if [ $NEXT_START -le $TOTAL ]; then
    echo "Submitting dependency for next batch starting at $NEXT_START"

    sbatch --dependency=afterok:${JOBID} Batches_AF3.sh $NEXT_START
else
    echo "All batches submitted."
fi