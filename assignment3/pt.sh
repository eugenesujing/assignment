#!/bin/bash
#
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00
#SBATCH --mem=1G
#SBATCH --partition=slow

./page_rank_push_parallel_atomic --nThreads 8 --nIterations 20 --inputFile input_graphs/roadNet-CA --strategy 2
