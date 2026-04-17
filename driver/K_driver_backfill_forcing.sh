#!/bin/bash

# This is a command line tool for forcing generation on klone.
# It is meant to be run as a batch job using a command like:
# sbatch -p cpu-g2 -A macc ./K_driver_backfill_forcing.sh wgh2 2019.07.04 2019.07.04 tide01 blank
#
# Thus there are 5 required arguments that would usually go with the flags
# -g -0 -1 -f -gtx when using driver_forcing00.py.

## Nodes
#SBATCH --nodes=1

## Tasks per node
#SBATCH --ntasks-per-node=3

## Walltime, long for multi-year forcing generation
#SBATCH --time=24:00:00

## Set memory use. Each slice (32 cores) has 256G
#SBATCH --mem=128G

# Do not return until the job is finished. Since we comment this out it will return to shell
# immediately, sending the sbatch job to another node. If we had kept it we could use "&" at
# the end of our command line, and then wwait for some screen output to tell us the job
# is done (useful for testing)
## #SBATCH --wait

#SBATCH --cpus-per-task=10

#SBATCH --job-name='frc_backfill'

# We defined and exported conda_source and LOd in .bashrc
source ${conda_source}
conda activate loenv

# These exports, along with the cpus-per-task set above
# (and workers=10 in the nearest neighbor code, e.g. in tracker2)
# allow fast multi-threading in python.
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK      # For OpenMP
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK # For OpenBLAS
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK      # For Intel MKL
export VECLIB_MAXIMUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

# rename command line arguments
gridname=$1
ds0=$2
ds1=$3
frc=$4
gtx=$5

python3 $LOd/driver_forcing00.py -g $gridname -r backfill -0 $ds0 -1 $ds1 -f $frc -k True -gtx $gtx > $LOd"/backfill_forcing_"$gridname"_"$frc"_"$ds0"_"$ds1".log"

# python3 $LOd/driver_forcing00.py -g cas7 -r forecast -f atm02 -k True > $LOd/atm02_cas7.log

# python3 $LOd/driver_forcing00.py -g cas7 -r forecast -f ocnG01 -do_bio True -k True > $LOd/ocnG01_cas7.log

# python3 $LOd/driver_forcing00.py -g cas7 -r forecast -tP trapsP01 -f trapsN00 -k True > $LOd/trapsN00_cas7.log
