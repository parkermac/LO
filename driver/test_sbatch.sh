#!/bin/bash

## Nodes
#SBATCH --nodes=1

## Tasks per node
#CBATCH --ntasks=10
#SBATCH --ntasks-per-node=10

## Walltime 
#SBATCH --time=01:00:00

## Use all memory on the node [0 to use all, or 128G]
#SBATCH --mem=32G

# Do not return until the job is finished
#SBATCH --wait

source /gscratch/macc/parker/miniconda3/etc/profile.d/conda.sh

conda activate loenv

# echo -e "Pre: $(date)\n" > /gscratch/macc/parker/LO/driver/sbatch_test.txt
# conda list >> /gscratch/macc/parker/LO/driver/sbatch_test.txt

LOd=/gscratch/macc/parker/LO/driver

python3 $LOd/driver_forcing00.py -g cas7 -r forecast -f tide01 > $LOd/tide01_cas7.log

python3 $LOd/driver_forcing00.py -g cas7 -r forecast -f atm02 > $LOd/atm02_cas7.log