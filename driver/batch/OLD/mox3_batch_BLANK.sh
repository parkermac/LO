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

module load icc_17-impi_2017
module load netcdf_fortran+c_4.4.1.1-icc_17

mpirun -np $np_num$ $roms_ex_dir$/$roms_ex_name$ $roms_out_dir$/liveocean.in > $roms_out_dir$/log.txt

