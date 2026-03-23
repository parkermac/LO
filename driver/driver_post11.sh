#!/bin/bash

## Nodes
#SBATCH --nodes=1

## Tasks per node
#CBATCH --ntasks=10
#SBATCH --ntasks-per-node=10

## Walltime 
#SBATCH --time=01:00:00

## Use all memory on the node [0 to use all, or 128G]
#SBATCH --mem=64G

# Do not return until the job is finished
#SBATCH --wait

#SBATCH --cpus-per-task=4

source /gscratch/macc/parker/miniconda3/etc/profile.d/conda.sh

conda activate loenv

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK      # For OpenMP
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK # For OpenBLAS
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK      # For Intel MKL
export VECLIB_MAXIMUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

# echo -e "Pre: $(date)\n" > /gscratch/macc/parker/LO/driver/sbatch_test.txt
# conda list >> /gscratch/macc/parker/LO/driver/sbatch_test.txt

LOd=/gscratch/macc/parker/LO/driver

python3 $LOd/driver_post11.py -gtx cas7_t2_x11b -ro 0 -r forecast -test True -override_cmd_list_test True < /dev/null > $LOd/post11.log

# python3 $LOd/driver_post11.py -gtx cas7_t2_x11b -ro 0 -r forecast < /dev/null > $LOd/post11.log
