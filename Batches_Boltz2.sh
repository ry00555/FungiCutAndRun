#!/bin/bash
#SBATCH --job-name=ASH1
#SBATCH --partition=inter_p
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/ry00555/Boltz2/ASH1/logs/boltz2_batch.%j.out
#SBATCH --error=/scratch/ry00555/Boltz2/ASH1/logs/boltz2_batch.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --output=/scratch/ry00555/Boltz2/ASH1/logs/boltz2_%A_%a.out
#SBATCH --error=/scratch/ry00555/Boltz2/ASH1/logs/boltz2_%A_%a.err

TOTAL=294       
BATCH_SIZE==10    # how many jobs to submit at once (stay under QOS limit of 20)
MAX_RUNNING=8    # how many run simultaneously (QOS limit is 8)
SCRIPT="/home/ry00555/Research/FungiCutAndRun/Boltz2_Step1.sh"

START=${1:-0}
END=$((START + BATCH_SIZE - 1))

if [ $END -gt $TOTAL ]; then
    END=$TOTAL
fi

echo "Submitting Boltz2 ASH1 array ${START}-${END}%${MAX_RUNNING}"

JOBID=$(sbatch --parsable --array=${START}-${END}%${MAX_RUNNING} $SCRIPT)

if [ -z "$JOBID" ]; then
    echo "Submission failed. Exiting."
    exit 1
fi

echo "Submitted job $JOBID"

NEXT_START=$((END + 1))

if [ $NEXT_START -le $TOTAL ]; then
    echo "Submitting dependency for next batch starting at $NEXT_START"
    sbatch --dependency=afterok:${JOBID} Boltz2_Step1.sh $NEXT_START
else
    echo "All batches submitted."
fi