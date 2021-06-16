#!/bin/bash
#
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00
#SBATCH --mem=1G
#SBATCH --partition=slow

srun /home/jsa306/assignment3/page_rank_pull_parallel --nThreads 5 --nIterations 10 --inputFile /home/jsa306/assignment3/input_graphs/roadNet-CA