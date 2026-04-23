#!/bin/bash

## Nodes
#SBATCH --nodes=1

## Tasks per node
#SBATCH --ntasks-per-node=1

## Walltime 
#SBATCH --time=06:00:00

## Set memory use (what is best for compute nodes?)
#SBATCH --mem=0

# Do not return until the job is finished
#SBATCH --wait

#SBATCH --exclusive

#SBATCH --cpus-per-task=40

source /gscratch/macc/parker/miniconda3/etc/profile.d/conda.sh

conda activate loenv

# These exports, along with the cpus-per-task set above
# (and workers=-1 in the nearest neighbor code, e.g. in tracker2)
# allow fast multi-threading in python.

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK      # For OpenMP
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK # For OpenBLAS
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK      # For Intel MKL
export VECLIB_MAXIMUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

LOd=/gscratch/macc/parker/LO/driver

# For testing
python3 $LOd/driver_post1.py -gtx cas7_t2_x11b -ro 0 -r forecast -test True -override_cmd_list_test True < /dev/null > $LOd/post_K1_test.log
