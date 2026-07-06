#!/bin/bash
#SBATCH --job-name=protein_prep
#SBATCH --nodes=1                
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=1000
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err

# Purge loaded modules
module purge 

# Activate conda environment
module load Miniconda3
source activate /path/to/env/biobb_md # e.g. /shared/work/BiobbWorkflows/envs/biobb_md

# Input files
DATA_FOLDER=../../data
INPUT_PDB=$DATA_FOLDER/1r9o.pdb
OUTPUT_PATH=output

# Launch workflow
protein_preparation --config input.yml \
                    --input_pdb $INPUT_PDB \
                    --output $OUTPUT_PATH \
                    --ph 7 \
                    --cap_ter \
                    --output_format gromacs \
                    --modeller_key <MODELLER-KEY>

# Remove report files
rm report_*
rm reduce_info.log
