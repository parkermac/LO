#!/bin/bash

# Test of using klone compute nodes.
# run as:
# ./klone_job_test.sh < /dev/null > job.log &

source ~/.bashrc
srun -p cpu-g2 -A macc -n1 -c1
conda activate loenv
echo "Date = $(date)"
conda deactivate