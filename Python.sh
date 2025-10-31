#!/bin/bash
#SBATCH --job-name=FixMetas
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=1:00:00
#SBATCH --output=FixMetas.%j.out
#SBATCH --error=FixMetas.%j.err

# Load Python and required modules
module load Python/3.11.3-GCCcore-12.3.0
module load SciPy-bundle/2024.05-gfbf-2024a

# Run Python script
python /home/ry00555/Research/FungiCutAndRun/FixMetas.py
