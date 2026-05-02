#!/bin/bash

## run using
# sbatch -p cpu-g2 -A macc ./speed_test.sh &

## Nodes
#SBATCH --nodes=1

## Tasks per node
#SBATCH --ntasks-per-node=1

## Walltime 
#SBATCH --time=01:00:00

## Set memory use. Each slice (32 cores) has 256G
#SBATCH --mem=128G

# Do not return until the job is finished
#SBATCH --wait

## Assume xarray to_zarr() uses multi-threading.
#SBATCH --cpus-per-task=10

source /gscratch/macc/parker/miniconda3/etc/profile.d/conda.sh
conda activate loenv
# These exports, along with the cpus-per-task set above
# allow fast multi-threading in python.
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK      # For OpenMP
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK # For OpenBLAS
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK      # For Intel MKL
export VECLIB_MAXIMUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

dir0='/gscratch/macc/parker'
dir1=${dir0}'/LO/extract/to_zarr'

s5cmd cp -acl public-read 's3://liveocean-pmacc/LO_roms/cas7_t2_x11b/f2026.05.01/ocean_his_0001.nc' 's3://liveocean-pmacc/tmp/'
s5cmd sync 's3://liveocean-pmacc/LO_roms/cas7_t2_x11b_zarr/f2026.05.01/h_01.zarr/*' '/gscratch/macc/parker/tmp/h_01.zarr/'

python3 ${dir1}/speed_test.py > ${dir1}/speed_test.log

