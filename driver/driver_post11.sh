#!/bin/bash

## Nodes
#SBATCH --nodes=1

## Tasks per node
#CBATCH --ntasks=64
#SBATCH --ntasks-per-node=64

## Walltime 
#SBATCH --time=01:00:00

## Use all memory on the node [0 to use all, or 128G]
#SBATCH --mem=64G

# Do not return until the job is finished
#SBATCH --wait

source /gscratch/macc/parker/miniconda3/etc/profile.d/conda.sh

conda activate loenv

# echo -e "Pre: $(date)\n" > /gscratch/macc/parker/LO/driver/sbatch_test.txt
# conda list >> /gscratch/macc/parker/LO/driver/sbatch_test.txt

LOd=/gscratch/macc/parker/LO/driver

python3 $LOd/driver_post11.py -gtx cas7_t2_x11b -ro 0 -r forecast -test True -override_cmd_list_test True < /dev/null > $LOd/post11.log
