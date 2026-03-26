#!/bin/bash

## Nodes
#SBATCH --nodes=1

## Tasks per node
#CBATCH --ntasks=3
#SBATCH --ntasks-per-node=3

## Walltime 
#SBATCH --time=12:00:00

## Use all memory on the node [0 to use all, or 128G]
#SBATCH --mem=64G

# Do not return until the job is finished
#SBATCH --wait

#SBATCH --cpus-per-task=10

source /gscratch/macc/parker/miniconda3/etc/profile.d/conda.sh

conda activate loenv

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK      # For OpenMP
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK # For OpenBLAS
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK      # For Intel MKL
export VECLIB_MAXIMUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

LOd=/gscratch/macc/parker/LO/driver

python3 $LOd/driver_post3.py -gtx oly2_t2_xn11b -ro 0 -r forecast < /dev/null > $LOd/post_K3.log
