#!/bin/bash

## run using
# sbatch --array1-25 ./driver_to_zarr.sh  cas7_t2_x11b 2026.05.01 &

## Group
#SBATCH -A macc

## Node type
#SBATCH -p cpu-g2

## Number of nodes
#SBATCH --nodes=1

## Tasks per node (number of history files in a day)
#SBATCH --ntasks-per-node=25

## Walltime 
#SBATCH --time=01:00:00

## Set memory use. Each slice (32 cores) has 256G. Use 0 to use all.
#SBATCH --mem=0

# Do not return until the job is finished
#SBATCH --wait

## Assume some python code uses multi-threading.
#SBATCH --cpus-per-task=6
## Note: 25*6 = 150, which is less than the 192 coes we have on a node.

## Reserve all the cores on this node.
#SBATCH --exclusive

source /gscratch/macc/parker/miniconda3/etc/profile.d/conda.sh
conda activate loenv
# These exports, along with the cpus-per-task set above
# allow fast multi-threading in python.
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK      # For OpenMP
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK # For OpenBLAS
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK      # For Intel MKL
export VECLIB_MAXIMUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Path to this directory
dir0='/gscratch/macc/parker'
dir1=${dir0}'/LO/extract/to_zarr'

gtx=$1
ds0=$2

# put the files for this day in /var/tmp
bucket='liveocean-pmacc'
s5cmd sync "s3://"${bucket}"/LO_roms/"${gtx}"/f"${ds0}"/*.nc" "/var/tmp/"

python3 ${dir1}/worker_to_zarr.py -his_str $SLURM_ARRAY_TASK_ID > ${dir1}"/test_"$SLURM_ARRAY_TASK_ID".log"

s5cmd sync --acl public-read "/var/tmp/*.zarr" "s3://"${bucket}"/LO_roms/"${gtx}"_zarr/f"${ds0}"/"
