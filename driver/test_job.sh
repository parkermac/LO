#!/bin/bash

#SBATCH --job-name=test
#SBATCH --account=macc
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=02:00:00

#SBATCH --export=ALL
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err

source ~/.bashrc
conda activate loenv
echo "Date = $(date)"
python test_job_worker.py &> job.log
conda deactivate