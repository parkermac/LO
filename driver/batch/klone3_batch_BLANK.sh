#!/bin/bash

## Job Name
#SBATCH --job-name=$jobname$

## Resources
## Nodes
#SBATCH --nodes=$node_num$
## Tasks per node
#SBATCH --ntasks-per-node=$cores_per_node$

## Walltime 
#SBATCH --time=02:00:00

## Memory per node
#SBATCH --mem=128G

## bug fix 2022.09.14
##. /opt/ohpc/admin/lmod/lmod/init/profile

module purge
module load intel/oneAPI
NFDIR=/gscratch/macc/local/netcdf-ifort/
export LD_LIBRARY_PATH=${NFDIR}/lib:${LD_LIBRARY_PATH}

mpirun -np $np_num$ $roms_ex_dir$/$roms_ex_name$ $roms_out_dir$/liveocean.in > $roms_out_dir$/log.txt

