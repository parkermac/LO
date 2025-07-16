#!/bin/bash

# Test of using klone compute nodes.

source ~/.bashrc
srun -p cpu-g2 -A macc -n1 -c1
conda activate loenv
echo "Date = $(date)" > job_test.txt
conda deactivate