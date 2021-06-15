#!/bin/bash
#
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00
#SBATCH --mem=1G
#SBATCH --partition=slow

srun /home/jsa306/assignment2/test_scripts/ellipse_area_tester.py --execPath=/home/jsa306/assignment2/ellipse_area_parallel --scriptPath=/home/jsa306/assignment2/test_scripts/ellipse_area_evaluator.py
