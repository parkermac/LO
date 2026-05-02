#!/bin/bash

## run using
# pmsrun2
# sbatch -p cpu-g2 -A macc $LOd/driver_to_zarr.sh &

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
dir1=${dir0}'/extract/to_zarr'

gtx='cas7_t2_x11b'
ds0='2026.05.01'
his_srt='01'

# put the file in /var/tmp
bucket='liveocean-pmacc'
s5cmd cp "s3://"${bucket}"/LO_roms/"${gtx}"/f"${ds0}"/ocean_his_00"${his_str}".nc" "/var/tmp/h_"${his_str}".nc"

python3 ${dir1}/worker_to_zarr.py -his_str ${his_str} > ${dir1}/test.log

s5cmd cp "/var/tmp/h_"${his_num}".zarr" "s3://"${bucket}"/LO_roms/"${gtx}"_zarr/f"${ds0}"/ocean_his_00"${his_str}".zarr"
