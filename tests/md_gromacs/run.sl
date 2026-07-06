#!/bin/bash
#SBATCH --job-name=md_gromacs
#SBATCH --nodes=1                
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=1000
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err
#SBATCH --gres=gpu:1

# Purge loaded modules
module purge 

# Load GROMACS module
ml GROMACS
GMX_BIN=$(which gmx)

# Activate conda environment
module load Miniconda3
source activate /path/to/env/biobb_md # e.g. /shared/work/BiobbWorkflows/envs/biobb_md

# Unset PYTHONPATH to avoid conflicts with GROMACS module
unset PYTHONPATH

# Input files
INPUT_PDB=../protein_preparation/output/1r9o.pdb
LIGANDS_FOLDER=../ligand_parameterization/output/topologies
OUTPUT_PATH=output

# Launch workflow
md_gromacs --input_pdb $INPUT_PDB \
           --ligands_folder $LIGANDS_FOLDER \
           --output $OUTPUT_PATH \
           --restart \
           --prod_time 0.1 \
           --equil_time 0.1 \
           --use_gpu \
           --gmx_bin $GMX_BIN