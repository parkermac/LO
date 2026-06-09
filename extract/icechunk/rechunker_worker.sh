#!/bin/bash

## run using a command of the form
## sbatch --array=1-?? ./rechunker_worker.sh cas7_t2_x11b 2026.05.01 2026.05.31 average

## Group
#SBATCH -A macc

## Node type
#SBATCH -p cpu-g2

## Number of nodes
#SBATCH --nodes=1

## Tasks per node
#SBATCH --ntasks-per-node=1

## Walltime 
#SBATCH --time=00:10:00

## Set memory use. Each slice (32 cores) has 256G, so 8GB per core.
#SBATCH --mem=8G

# Do not return until the job is finished (use with & on the command line)
#SBATCH --wait

# activate our conda environment
source /gscratch/macc/parker/miniconda3/etc/profile.d/conda.sh
conda activate loenv

# make a unique directory for each task
export TMPDIR=/var/tmp/$USER/$SLURM_ARRAY_JOB_ID_$SLURM_ARRAY_TASK_ID
mkdir -p $TMPDIR
# NOTE: exported variables are available to all jobs launched from this shell

# paths
dir0='/gscratch/macc/parker'
dir1=${dir0}'/LO/extract/icechunk'

# rename command line arguments for clarity
indir=$1
outdir=$2
gtagex=$3
this_ds0=$4
this_ds1=$5
list_type=$6

# Run the worker
python3 ${dir1}/rechunker_worker_one_time.py -tid ${SLURM_ARRAY_TASK_ID} -gtagex ${gtagex} \
    -indir ${indir} -outdir ${outdir} \
    -this_ds0 ${this_ds0} -this_ds1 ${this_ds1} -list_type ${list_type}
