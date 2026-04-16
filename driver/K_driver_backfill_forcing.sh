#!/bin/bash

# This is a command line tool for forcing generation on klone.
# It is meant to be run as a batch job using a command like:
# sbatch -p cpu-g2 -A macc ./K_driver_backfill_forcing.sh wgh2 2019.07.04 2019.07.04 tide01
# Thus there are 4 required arguments that would usually go with the flags
# -g -0 -1 -f when using driver_forcing00.py.

## Nodes
#SBATCH --nodes=1

## Tasks per node
#SBATCH --ntasks-per-node=3

## Walltime 
#SBATCH --time=01:00:00

## Set memory use. Each slice (32 cores) has 256G
#SBATCH --mem=128G

# Do not return until the job is finished
## #SBATCH --wait

#SBATCH --cpus-per-task=10

# source /gscratch/macc/parker/miniconda3/etc/profile.d/conda.sh
source $conda_source

conda activate loenv

# These exports, along with the cpus-per-task set above
# (and workers=10 in the nearest neighbor code, e.g. in tracker2)
# allow fast multi-threading in python.

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK      # For OpenMP
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK # For OpenBLAS
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK      # For Intel MKL
export VECLIB_MAXIMUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

# How to make this automatic for all users?
# LOd=/gscratch/macc/parker/LO/driver

# Assumes you are running it in driver
# LOd=$PWD

gridname=$1
ds0=$2
ds1=$3
frc=$4

python3 $LOd/driver_forcing00.py -g $gridname -r backfill -0 $ds0 -1 $ds1 -f $frc -k True > $LOd"/"$gridname"_"$frc"_"$ds0"_"$ds1".log"

# python3 $LOd/driver_forcing00.py -g cas7 -r forecast -f atm02 -k True > $LOd/atm02_cas7.log

# python3 $LOd/driver_forcing00.py -g cas7 -r forecast -f ocnG01 -do_bio True -k True > $LOd/ocnG01_cas7.log

# python3 $LOd/driver_forcing00.py -g cas7 -r forecast -tP trapsP01 -f trapsN00 -k True > $LOd/trapsN00_cas7.log
