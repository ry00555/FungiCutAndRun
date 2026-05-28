#!/bin/bash
#SBATCH --job-name=EAF3 #change this line for record keeping
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8gb
#SBATCH --time=1:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --output=AF3_Pool_Scraper.%j.out
#SBATCH --error=AF3_Pool_Scraper.%j.err

KEY_CSV="/scratch/ry00555/EpigeneticMemoryPaper2026/AlphaFold3/PooledPPI/EAF3_AF3_PooledJSONs/EAF3_pool_key_scraper_lower.csv"
BASE_DIR="/scratch/ry00555/EpigeneticMemoryPaper2026/AlphaFold3/PooledPPI/EAF3_AF3_PooledJSON_output"
OUTPUT_CSV="/scratch/ry00555/EpigeneticMemoryPaper2026/AlphaFold3/PooledPPI/EAF3_AF3_ipTM_scores.csv"
SCRIPTS_DIR="/home/ry00555/Research/FungiCutAndRun"

python3 ${SCRIPTS_DIR}/AF3pool_scraper_v1.py \
    $KEY_CSV \
    $BASE_DIR \
    $OUTPUT_CSV \
    --row-index 1

echo "Done! Scores written to $OUTPUT_CSV"
