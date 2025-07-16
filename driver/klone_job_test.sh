#!/bin/bash

# Test of using klone compute nodes.

source ~/.bashrc
srun -p compute -A macc --pty bash -l
conda activate loenv
echo "Date = $(date)" > job_test.txt
conda deactivate
exit