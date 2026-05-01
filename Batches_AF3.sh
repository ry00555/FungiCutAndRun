#!/bin/bash

TOTAL=584
BATCH_SIZE=5
MAX_RUNNING=1
SCRIPT="EAF3_AF3_Pool.sh"

for START in $(seq 1 $BATCH_SIZE $TOTAL); do
    END=$((START + BATCH_SIZE - 1))

    if [ $END -gt $TOTAL ]; then
        END=$TOTAL
    fi

    echo "Submitting array ${START}-${END}%${MAX_RUNNING}"

    sbatch --array=${START}-${END}%${MAX_RUNNING} $SCRIPT

    sleep 2
done