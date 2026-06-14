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
MSA_THROTTLE=20   # MSA on batch partition — can run many at once
INF_THROTTLE=9    # INF on gpu_p — keep low for fairshare
MSA_SCRIPT="PRC2_AF3_MSA_Pool.sh"
INF_SCRIPT="PRC2_AF3_INF_Pool.sh"
INPUT_DIR="/scratch/ry00555/RNASeqPaper2026/Proteome/PRC2_Proteome_pools/PooledPPI/PRC2_AF3_PooledJSONs"
OUTPUT_DIR="/scratch/ry00555/RNASeqPaper2026/Proteome/PRC2_Proteome_pools/PooledPPI/PRC2_AF3_PooledJSON_output"
# ─────────────────────────────────────────────────────────────────────────────

START=${1:-1}

# ── If START is "recheck" run the missing pools check ─────────────────────────
if [ "$START" == "recheck" ]; then
    echo "Running missing pools check..."
    missing_msa=()
    missing_inf_only=()

    for i in $(seq 1 $TOTAL); do
        pool=$(printf "prc2_pool%03d" $i)
        dir="${OUTPUT_DIR}/${pool}"
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
        echo "All pools complete!"
        exit 0
    fi

    PREV_JOBID=""

    # Submit MSA only for truly missing pools
    if [ ${#missing_msa[@]} -gt 0 ]; then
        CHUNK=200
        for ((start=0; start<${#missing_msa[@]}; start+=CHUNK)); do
            chunk=("${missing_msa[@]:$start:$CHUNK}")
            indices="${chunk[0]}"
            for idx in "${chunk[@]:1}"; do indices="${indices},${idx}"; done
            if [ -z "$PREV_JOBID" ]; then
                MSA_JOBID=$(sbatch --parsable --array=${indices}%${MSA_THROTTLE} $MSA_SCRIPT)
            else
                MSA_JOBID=$(sbatch --parsable --dependency=afterany:${PREV_JOBID} --array=${indices}%${MSA_THROTTLE} $MSA_SCRIPT)
            fi
            [ -z "$MSA_JOBID" ] && echo "MSA submission failed." && exit 1
            echo "MSA submitted: $MSA_JOBID (${#chunk[@]} pools)"
            INF_JOBID=$(sbatch --parsable --dependency=afterok:${MSA_JOBID} --array=${indices}%${INF_THROTTLE} $INF_SCRIPT)
            echo "INF submitted: $INF_JOBID"
            PREV_JOBID=$INF_JOBID
        done
    fi

    # Trigger inf_scan for pools that have MSA but no INF
    # inf_scan uses line numbers so it handles these correctly
    if [ -n "$PREV_JOBID" ]; then
        sbatch --dependency=afterany:${PREV_JOBID} PRC2_AF3_Pools.sh inf_scan
    else
        sbatch PRC2_AF3_Pools.sh inf_scan
    fi

    exit 0
fi

# ── If START is "inf_scan" — find MSA-complete pools and submit INF ───────────
if [ "$START" == "inf_scan" ]; then
    echo "Scanning for pools ready for INF..."

    line_num=0
    ready=()
    while IFS= read -r json_file; do
        line_num=$((line_num + 1))
        pool=$(basename "$json_file" .json | tr '[:upper:]' '[:lower:]')
        msa_json="${OUTPUT_DIR}/${pool}/${pool}_data.json"
        inf_done="${OUTPUT_DIR}/${pool}/seed-1_sample-0/summary_confidences.json"
        if [ -f "$msa_json" ] && [ ! -f "$inf_done" ]; then
            ready+=($line_num)
        fi
    done < <(ls $INPUT_DIR/*.json)

    echo "Found ${#ready[@]} pools ready for INF"

    if [ ${#ready[@]} -eq 0 ]; then
        # Check if MSAs are still running
        msa_running=$(squeue --me --partition=batch --name=PRC2_MSA --noheader | wc -l)
        if [ "$msa_running" -gt 0 ]; then
      echo "No INF ready yet but $msa_running MSA jobs still running — reschedule scan in 2hrs"
      sbatch --begin=now+2hour PRC2_AF3_Pools.sh inf_scan
  else
            echo "No MSAs running and no INF ready — triggering recheck"
            sbatch PRC2_AF3_Pools.sh recheck
        fi
        exit 0
    fi

    # Submit in batches of 8 with dependencies
    prev_job=""
    for ((start=0; start<${#ready[@]}; start+=8)); do
        chunk=("${ready[@]:$start:8}")
        indices="${chunk[0]}"
        for idx in "${chunk[@]:1}"; do
            indices="${indices},${idx}"
        done

        if [ -z "$prev_job" ]; then
            job=$(sbatch --parsable --array=$indices $INF_SCRIPT)
        else
            job=$(sbatch --parsable --dependency=afterany:${prev_job} --array=$indices $INF_SCRIPT)
        fi

        if [ -n "$job" ]; then
            echo "INF submitted: $job ($indices)"
            prev_job=$job
        else
            echo "FAILED (QOS limit): $indices — stopping here"
            break
        fi
    done

    # Schedule next scan after last INF batch finishes to catch new MSAs
    if [ -n "$prev_job" ]; then
        sbatch --dependency=afterany:${prev_job} PRC2_AF3_Pools.sh inf_scan
        echo "Next INF scan scheduled after job $prev_job"
    fi

    exit 0
fi


# ── If START is "large_scan" — find oversized MSA-complete pools and submit with unifiedmem ──
if [ "$START" == "large_scan" ]; then
    echo "Scanning for oversized pools (>5120 tokens) ready for INF..."

    INF_LARGE_SCRIPT="PRC2_AF3_INF_Pool_large.sh"
    line_num=0
    ready=()
    while IFS= read -r json_file; do
        line_num=$((line_num + 1))
        pool=$(basename "$json_file" .json | tr '[:upper:]' '[:lower:]')
        msa_json="${OUTPUT_DIR}/${pool}/${pool}_data.json"
        inf_done="${OUTPUT_DIR}/${pool}/seed-1_sample-0/summary_confidences.json"

        # Skip if INF already done
        [ -f "$inf_done" ] && continue

        # Skip if no MSA
        [ ! -f "$msa_json" ] && continue

        # Check token count
        total=$(python3 -c "
import json
d = json.load(open('$json_file'))
print(sum(len(s['protein']['sequence']) for s in d['sequences']))
" 2>/dev/null)

        if [ -n "$total" ] && [ "$total" -gt 5120 ]; then
            ready+=($line_num)
        fi
    done < <(ls $INPUT_DIR/*.json)

    echo "Found ${#ready[@]} oversized pools ready for INF"

    if [ ${#ready[@]} -eq 0 ]; then
        echo "All oversized pools complete!"
        exit 0
    fi

    # Submit in batches of 8 with dependencies
    prev_job=""
    for ((start=0; start<${#ready[@]}; start+=8)); do
        chunk=("${ready[@]:$start:8}")
        indices="${chunk[0]}"
        for idx in "${chunk[@]:1}"; do
            indices="${indices},${idx}"
        done

        if [ -z "$prev_job" ]; then
            job=$(sbatch --parsable --array=$indices $INF_LARGE_SCRIPT)
        else
            job=$(sbatch --parsable --dependency=afterany:${prev_job} --array=$indices $INF_LARGE_SCRIPT)
        fi

        if [ -n "$job" ]; then
            echo "INF_large submitted: $job ($indices)"
            prev_job=$job
        else
            echo "FAILED (QOS limit): $indices — stopping here"
            break
        fi
    done

    # Schedule next scan after last batch finishes
    if [ -n "$prev_job" ]; then
        sbatch --dependency=afterany:${prev_job} PRC2_AF3_Pools.sh large_scan
        echo "Next large_scan scheduled after job $prev_job"
    fi

    exit 0
fi

# ── Normal batch submission ───────────────────────────────────────────────────
END=$((START + BATCH_SIZE - 1))
[ $END -gt $TOTAL ] && END=$TOTAL

echo "Submitting batch: pools ${START}-${END}"

MSA_JOBID=$(sbatch --parsable --array=${START}-${END}%${MSA_THROTTLE} $MSA_SCRIPT)
[ -z "$MSA_JOBID" ] && echo "MSA submission failed." && exit 1
echo "MSA job submitted: $MSA_JOBID (pools ${START}-${END})"

INF_JOBID=$(sbatch --parsable \
    --dependency=afterok:${MSA_JOBID} \
    --array=${START}-${END}%${INF_THROTTLE} \
    $INF_SCRIPT)
[ -z "$INF_JOBID" ] && echo "INF submission failed." && exit 1
echo "INF job submitted: $INF_JOBID"

NEXT_START=$((END + 1))
if [ $NEXT_START -le $TOTAL ]; then
    sbatch --dependency=afterany:${MSA_JOBID} PRC2_AF3_Pools.sh $NEXT_START
else
    echo "All batches submitted — scheduling recheck for missing pools"
    sbatch --dependency=afterany:${INF_JOBID} PRC2_AF3_Pools.sh recheck
fi
