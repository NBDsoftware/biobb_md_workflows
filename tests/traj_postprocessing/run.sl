#!/bin/bash
#SBATCH --job-name=traj_postpro
#SBATCH --nodes=1                
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=500
#SBATCH --output=report_%j.out
#SBATCH --error=report_%j.err

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
# WARNING: This test requires an external raw trajectory not contained in the repository
INPUT_FOLDER=/path/to/MD/simulation/
INPUT_TPR=$INPUT_FOLDER/topology.tpr
INPUT_XTC=$INPUT_FOLDER/raw_md_1_part0.xtc
INPUT_STR=/path/to/topology/topology.pdb

OUTPUT_PATH=output
OUTPUT_TRAJ=trajectory.xtc
OUTPUT_STRUCTURE=structure.pdb

# Launch workflow
traj_postprocessing --input_top $INPUT_TPR \
                    --input_traj $INPUT_XTC \
                    --input_structure $INPUT_STR \
                    --output $OUTPUT_PATH \
                    --gmx_bin $GMX_BIN \
                    --output_traj $OUTPUT_TRAJ \
                    --output_str $OUTPUT_STRUCTURE