#!/bin/bash
#
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00
#SBATCH --mem=1G
#SBATCH --partition=slow

srun python /home/jsa306/assignment3/scripts/page_rank_tester.pyc --execPath=/home/jsa306/assignment3/ --scriptPath=/home/jsa306/assignment3/scripts/page_rank_evaluator.pyc --inputPath=/home/jsa306/assignment3/input_graphs/