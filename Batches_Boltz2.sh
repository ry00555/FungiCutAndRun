#!/bin/bash
#SBATCH --job-name=EAF3
#SBATCH --partition=inter_p
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/ry00555/Boltz2/logs/boltz2_batch.%j.out
#SBATCH --error=/scratch/ry00555/Boltz2/logs/boltz2_batch.%j.err

cd $SLURM_SUBMIT_DIR
PROJECT="EAF3" #change this and job name to dock protein
TOTAL=583 # total number of files minus one as Boltz is on a 0 scale not +1
BATCH_SIZE=12    # how many jobs to submit at once (stay under QOS limit of 20)
MAX_RUNNING=7    # how many run simultaneously (QOS limit is 8)
SCRIPT="/home/ry00555/Research/FungiCutAndRun/Boltz2_Step1.sh"

START=${1:-0}
END=$((START + BATCH_SIZE - 1))

if [ $END -gt $TOTAL ]; then
    END=$TOTAL
fi

echo "Submitting Boltz2 $PROJECT array ${START}-${END}%${MAX_RUNNING}"

JOBID=$(sbatch --parsable --array=${START}-${END}%${MAX_RUNNING} $SCRIPT)

if [ -z "$JOBID" ]; then
    echo "Submission failed. Exiting."
    exit 1
fi

echo "Submitted job $JOBID"

NEXT_START=$((END + 1))

if [ $NEXT_START -le $TOTAL ]; then
    echo "Submitting dependency for next batch starting at $NEXT_START"
    sbatch --dependency=afterany:${JOBID} Batches_Boltz2.sh $NEXT_START
    #sbatch --dependency=afterok:${JOBID} Batches_Boltz2.sh $NEXT_START use this one it you want to see where the chain fails (reasoning varies
else
    echo "All batches submitted."
fi
