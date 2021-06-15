#!/bin/bash
#
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00
#SBATCH --mem=1G
#SBATCH --partition=slow

srun /home/jsa306/assignment2/heat_transfer_parallel --nThreads 10 --tSteps 100 --gSize 5000 --mTemp 800 --iCX 0.26 --iCY 0.25
