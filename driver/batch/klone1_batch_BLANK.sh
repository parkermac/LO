#!/bin/bash

## Job Name
#SBATCH --job-name=LiveOcean

## Resources
## Nodes
#SBATCH --nodes=$node_num$
## Tasks per node
#SBATCH --ntasks-per-node=$cores_per_node$

## Walltime 
#SBATCH --time=02:00:00

## Memory per node
#SBATCH --mem=100G

module purge
## module load intel/oneAPI
## NFDIR=/gscratch/macc/local/netcdf-ifort/
## export LD_LIBRARY_PATH=${NFDIR}/lib:${LD_LIBRARY_PATH}

module load intel/oneAPI
LODIR=/gscratch/macc/local

NFDIR=${LODIR}/netcdf-ifort
NCDIR=${LODIR}//netcdf-icc
PIODIR=${LODIR}/pio2
PNDIR=${PNDIR}/pnetcdf

export LD_LIBRARY_PATH=${PIODIR}/lib:${PNDIR}:${NFDIR}/lib:${NCDIR}/lib:${LD_LIBRARY_PATH}

export PATH=/gscratch/macc/local/netcdf-ifort/bin:$PATH
export PATH=/gscratch/macc/local/netcdf-icc/bin:$PATH

mpirun -np $np_num$ $roms_ex_dir$/$roms_ex_name$ $roms_out_dir$/liveocean.in > $roms_out_dir$/log.txt

