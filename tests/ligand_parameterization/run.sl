#!/bin/bash
#SBATCH --job-name=ligand_param   # Job name
#SBATCH --nodes=1                 # node count
#SBATCH --ntasks=8                # total number of tasks across all nodes
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=1000
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err

# Purge loaded modules
module purge 

# Activate conda environment, see environment.yml
module load Miniconda3
source activate /path/to/env/biobb_md  # e.g. /shared/work/BiobbWorkflows/envs/biobb_md

# Input files
DATA_FOLDER=../../data
INPUT_PDB=$DATA_FOLDER/1r9o.pdb
CUSTOM_PARAMETERS=$DATA_FOLDER/custom_parameters
OUTPUT_PATH=output

# Launch workflow
ligand_parameterization --config input.yml \
                        --input_pdb $INPUT_PDB \
                        --forcefields protein.ff14SB \
                        --ligands HEM FLP \
                        --format gromacs \
                        --ligand_parameters $CUSTOM_PARAMETERS \
                        --output $OUTPUT_PATH

# Remove report files
rm report_*

# Remove tmp folders
rm -r sandbox_*
rm -r HEM.*