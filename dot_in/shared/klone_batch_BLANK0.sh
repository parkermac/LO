#!/bin/bash

## Job Name
#SBATCH --job-name=LiveOcean

## Resources
## Nodes
#SBATCH --nodes=$node_num$
## Tasks per node (Slurm assumes you want to run 28 tasks per node unless explicitly told otherwise)
#SBATCH --ntasks-per-node=$cores_per_node$

## Walltime 
#SBATCH --time=05:00:00

## Memory per node
#SBATCH --mem=100G

module purge
module load intel/oneAPI
NFDIR=/gscratch/macc/local/netcdf-ifort/
export LD_LIBRARY_PATH=${NFDIR}/lib:${LD_LIBRARY_PATH}

mpirun $roms_ex_dir$/oceanM $roms_out_dir$/liveocean.in > $roms_out_dir$/log.txt

