#!/bin/bash

## Job Name
#SBATCH --job-name=$jobname$

## Resources

## Nodes
#SBATCH --nodes=$number_of_nodes$

## Tasks per node
#SBATCH --ntasks-per-node=$ntasks_per_node$

## Walltime 
#SBATCH --time=02:00:00

## Use all memory on the node [0 to use all, or 128G]
#SBATCH --mem=$sbatch_mem$

## Request the entire node by using #SBATCH --exclusive
$sbatch_exclusive_line$

module purge
module load intel/oneAPI
NFDIR=/gscratch/macc/local/netcdf-ifort
NCDIR=/gscratch/macc/local/netcdf-icc
export LD_LIBRARY_PATH=${NFDIR}/lib:${NCDIR}/lib:${LD_LIBRARY_PATH}

echo -e "Pre: $(date)\n" # Timestamp to .out file before anything else
mpirun -np $np_num$ $roms_ex_dir$/$roms_ex_name$ $roms_out_dir$/liveocean.in > $roms_out_dir$/log.txt
MPIRUN_EXIT_CODE=$?  # Capture exit code of mpirun
echo -e "\nPost: $(date) with mpirun exit code: $MPIRUN_EXIT_CODE" # Timestamp to .out just before exit
exit $MPIRUN_EXIT_CODE  # Exit script with exit code from mpirun
