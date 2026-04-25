#!/bin/bash
#SBATCH --job-name=alphaMSAs_step1
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --time=12:00:00
#SBATCH --output=%AlphaFold2_MSABuilder.%j.out
#SBATCH --error=%AlphaFold2_MSABuilder.%j.err
#SBATCH --array=1-200
#Substep1 before running check and stop elements script
cd $SLURM_SUBMIT_DIR

ml purge
ml AlphaFold/2.3.2-foss-2023a-CUDA-12.1.1
export ALPHAFOLD_DATA_DIR=/db/AlphaFold/2.3.1
WORKDIR="/scratch/ry00555/AlphaFold"

file=$(awk "NR==${SLURM_ARRAY_TASK_ID}" $WORKDIR/FastaforAlphaFold/input.lst)

alphafold \
--models_to_relax=none \
--data_dir=$ALPHAFOLD_DATA_DIR \
--uniref90_database_path=$ALPHAFOLD_DATA_DIR/uniref90/uniref90.fasta \
--mgnify_database_path=$ALPHAFOLD_DATA_DIR/mgnify/mgy_clusters.fa \
--bfd_database_path=$ALPHAFOLD_DATA_DIR/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
--uniref30_database_path=$ALPHAFOLD_DATA_DIR/uniref30/UniRef30_2021_03 \
--pdb_seqres_database_path=$ALPHAFOLD_DATA_DIR/pdb_seqres/pdb_seqres.txt \
--template_mmcif_dir=$ALPHAFOLD_DATA_DIR/pdb_mmcif/mmcif_files \
--obsolete_pdbs_path=$ALPHAFOLD_DATA_DIR/pdb_mmcif/obsolete.dat \
--uniprot_database_path=$ALPHAFOLD_DATA_DIR/uniprot/uniprot.fasta \
--model_preset=multimer \
--num_multimer_predictions_per_model=1 \
--max_template_date=2023-10-01 \
--db_preset=full_dbs \
--output_dir=$WORKDIR/output/AF2/$(basename $file .fa) \
--fasta_paths=$WORKDIR/FastaforAlphaFold/$file
