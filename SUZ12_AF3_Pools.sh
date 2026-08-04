#!/bin/bash
#SBATCH --job-name=SUZ12_ctrl
#SBATCH --partition=inter_p
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/ry00555/RNASeqPaper2026/Proteome/SingleBaitPools/logs/SUZ12_ctrl.%j.out
#SBATCH --error=/scratch/ry00555/RNASeqPaper2026/Proteome/SingleBaitPools/logs/SUZ12_ctrl.%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ry00555@uga.edu

cd $SLURM_SUBMIT_DIR

# ══════════════════════════════════════════════════════════════════════
# SUZ12_AF3_Pools.sh
#
# Orchestrates MSA -> INF for the SUZ12 single-bait screen (1318 pools).
#
# ── Why the numbers below are what they are ────────────────────────────
# MSA runs on the CPU-only `batch` partition -- it is NOT subject to the
# gpu_p job cap, so we submit the bait's ENTIRE MSA array in one sbatch
# call, throttled to MSA_THROTTLE concurrent tasks. Per the workflow doc,
# GACRC is comfortable with up to ~200 simultaneous batch-partition MSA
# jobs; MSA_THROTTLE=50 here is this bait's 1/4 share so
# that running all 4 bait pipelines together still totals ~200.
# Bump this up if you're only running THIS bait's pipeline at a time.
#
# INF requires an A100 GPU and is subject to GACRC's ACCOUNT-WIDE cap of
# 8 running / 20 submitted (running+pending) gpu_p jobs -- shared across
# ALL 4 bait pipelines, not per bait. INF_THROTTLE=2 x 4 baits = 8
# (the hard running-job cap), and INF_BATCH=4 x 4 baits = 16
# pending array tasks (under the 20-job cap, leaving headroom). If you are
# only running one or two baits' pipelines at a time, raise INF_THROTTLE
# / INF_BATCH for those baits accordingly (e.g. INF_THROTTLE=8 if this is
# the only bait pipeline active).
# ══════════════════════════════════════════════════════════════════════

TOTAL=1318
MSA_THROTTLE=50
INF_THROTTLE=2
INF_BATCH=8
MSA_SCRIPT="SUZ12_AF3_MSA_Pool.sh"
INF_SCRIPT="SUZ12_AF3_INF_Pool.sh"
INPUT_DIR="/scratch/ry00555/RNASeqPaper2026/Proteome/SingleBaitPools/SUZ12/SUZ12_PooledPPI_InputJSONs"
OUTPUT_DIR="/scratch/ry00555/RNASeqPaper2026/Proteome/SingleBaitPools/SUZ12/SUZ12_PooledPPI_output"

MODE=${1:-start}

# ── start: submit the full MSA array once, then hand off to inf_scan ──────────
if [ "$MODE" == "start" ]; then
    mkdir -p "$OUTPUT_DIR" "/scratch/ry00555/RNASeqPaper2026/Proteome/SingleBaitPools/logs"

    n_json=$(ls "$INPUT_DIR"/*.json 2>/dev/null | wc -l)
    if [ "$n_json" -ne "$TOTAL" ]; then
        echo "WARNING: TOTAL=$TOTAL but found $n_json JSONs in $INPUT_DIR"
        echo "         Update TOTAL at the top of this script if this is not a re-run."
    fi

    echo "Submitting full MSA array for SUZ12: 1-${TOTAL}%${MSA_THROTTLE}"
    MSA_JOBID=$(sbatch --parsable --array=1-${TOTAL}%${MSA_THROTTLE} $MSA_SCRIPT)
    [ -z "$MSA_JOBID" ] && echo "MSA submission failed." && exit 1
    echo "MSA array submitted: $MSA_JOBID"

    sbatch --dependency=afterany:${MSA_JOBID} $0 inf_scan
    echo "inf_scan queued to start after MSA array $MSA_JOBID"
    exit 0
fi

# ── inf_scan: find MSA-complete / INF-incomplete pools, submit in small
#              GPU-budget-safe chunks, self-reschedule ────────────────────────
if [ "$MODE" == "inf_scan" ]; then
    echo "Scanning SUZ12 for pools ready for INF..."
    line_num=0
    ready=()
    while IFS= read -r json_file; do
        line_num=$((line_num + 1))
        pool=$(basename "$json_file" .json | tr '[:upper:]' '[:lower:]')
        msa_json="${OUTPUT_DIR}/${pool}/${pool}_data.json"
        inf_done="${OUTPUT_DIR}/${pool}/seed-1_sample-0/summary_confidences.json"

        [ ! -f "$msa_json" ] && continue
        [ -f "$inf_done" ] && continue

        ready+=($line_num)
    done < <(ls $INPUT_DIR/*.json)

    echo "Found ${#ready[@]} SUZ12 pools ready for INF"

    if [ ${#ready[@]} -eq 0 ]; then
        msa_running=$(squeue --me --partition=batch --name=SUZ12_MSA --noheader | wc -l)
        if [ "$msa_running" -gt 0 ]; then
            echo "No INF ready yet but $msa_running MSA jobs still running — reschedule scan in 1hr"
            sbatch --begin=now+1hour $0 inf_scan
        else
            echo "No MSAs running and no INF ready — triggering recheck"
            sbatch $0 recheck
        fi
        exit 0
    fi

    prev_job=""
    for ((start=0; start<${#ready[@]}; start+=INF_BATCH)); do
        chunk=("${ready[@]:$start:$INF_BATCH}")
        indices="${chunk[0]}"
        for idx in "${chunk[@]:1}"; do indices="${indices},${idx}"; done

        if [ -z "$prev_job" ]; then
            job=$(sbatch --parsable --array=${indices}%${INF_THROTTLE} $INF_SCRIPT)
        else
            job=$(sbatch --parsable --dependency=afterany:${prev_job} --array=${indices}%${INF_THROTTLE} $INF_SCRIPT)
        fi

        if [ -n "$job" ]; then
            echo "INF submitted: $job ($indices)"
            prev_job=$job
        else
            echo "FAILED (QOS limit): $indices — stopping here, will retry on next scan"
            break
        fi
    done

    if [ -n "$prev_job" ]; then
        sbatch --dependency=afterany:${prev_job} $0 inf_scan
        echo "Next inf_scan scheduled after job $prev_job"
    fi
    exit 0
fi

# ── recheck: catch anything that failed (missing MSA or missing INF) ──────────
if [ "$MODE" == "recheck" ]; then
    echo "Running missing-pools check for SUZ12..."
    missing_msa=()
    missing_inf_only=()

    for i in $(seq 1 $TOTAL); do
        pool_glob="${OUTPUT_DIR}/suz12_pool$(printf '%03d' $i)"
        # pool numbers >=1000 aren't zero-padded past 3 digits by the generator
        [ ! -d "$pool_glob" ] && pool_glob="${OUTPUT_DIR}/suz12_pool${i}"
        dir="$pool_glob"
        has_msa=0
        has_inf=0
        ls $dir/*_data.json &>/dev/null 2>&1 && has_msa=1
        ls $dir/seed-1_sample-0/summary_confidences.json &>/dev/null 2>&1 && has_inf=1
        [ $has_msa -eq 0 ] && missing_msa+=($i)
        [ $has_msa -eq 1 ] && [ $has_inf -eq 0 ] && missing_inf_only+=($i)
    done

    echo "Missing MSA: ${#missing_msa[@]}"
    echo "Missing INF only (MSA done): ${#missing_inf_only[@]}"

    if [ ${#missing_msa[@]} -eq 0 ] && [ ${#missing_inf_only[@]} -eq 0 ]; then
        echo "All SUZ12 pools complete!"
        exit 0
    fi

    if [ ${#missing_msa[@]} -gt 0 ]; then
        indices="${missing_msa[0]}"
        for idx in "${missing_msa[@]:1}"; do indices="${indices},${idx}"; done
        MSA_JOBID=$(sbatch --parsable --array=${indices}%${MSA_THROTTLE} $MSA_SCRIPT)
        echo "Resubmitted MSA: $MSA_JOBID (${#missing_msa[@]} pools)"
        sbatch --dependency=afterany:${MSA_JOBID} $0 inf_scan
    else
        sbatch $0 inf_scan
    fi
    exit 0
fi

echo "Unknown mode: $MODE (use start | inf_scan | recheck)"
exit 1
