#!/bin/bash

## Nodes
#SBATCH --nodes=1

## Tasks per node
#CBATCH --ntasks=1
#SBATCH --ntasks-per-node=1

## Walltime 
#SBATCH --time=00:01:00

## Use all memory on the node [0 to use all, or 128G]
#SBATCH --mem=1G

# Do not return until the job is finished
#SBATCH --wait

conda activate loenv

echo -e "Pre: $(date)\n" > /gscratch/macc/parker/LO/driver/sbatch_test.txt