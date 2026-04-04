#!/bin/bash

## run as
## sbatch -p cpu-g2 -A macc ./test_s5cmd.sh

## Nodes
#SBATCH --nodes=1

## Tasks per node
#CBATCH --ntasks=3
#SBATCH --ntasks-per-node=3

## Walltime 
#SBATCH --time=01:00:00

## Use all memory on the node [0 to use all, or 128G]
#SBATCH --mem=64G

# Do not return until the job is finished
#SBATCH --wait

#SBATCH --cpus-per-task=10

source /gscratch/macc/parker/miniconda3/etc/profile.d/conda.sh

conda activate loenv

# These exports, along with the cpus-per-task set above
# (and workers=10 in the nearest neighbor code, e.g. in tracker2)
# allow fast multi-threading in python.

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK      # For OpenMP
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK # For OpenBLAS
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK      # For Intel MKL
export VECLIB_MAXIMUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

# These lines were used in testing to see if the conda environment loaded correctly.
# echo -e "Pre: $(date)\n" > /gscratch/macc/parker/LO/driver/sbatch_test.txt
# conda list >> /gscratch/macc/parker/LO/driver/sbatch_test.txt

echo -e "Pre: $(date)\n" > /gscratch/macc/parker/LO/driver/sbatch_test.txt
which s3cmd >> /gscratch/macc/parker/LO/driver/sbatch_test.txt
which s5cmd >> /gscratch/macc/parker/LO/driver/sbatch_test.txt
