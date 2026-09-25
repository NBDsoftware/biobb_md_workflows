#!/bin/bash
#SBATCH --job-name=godmd_sim
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

# Adenylate Kinase (ADK): closed (1AKE) -> open (4AKE) conformational transition
OUTPUT_PATH=output

# Launch workflow
godmd_simulation --origin_pdb_code 1ake \
                  --origin_chain A \
                  --target_pdb_code 4ake \
                  --target_chain A \
                  --temp 400 \
                  --output $OUTPUT_PATH

# Remove report files
rm report_*
