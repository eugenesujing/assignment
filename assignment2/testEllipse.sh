#!/bin/bash
#
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00
#SBATCH --mem=1G
#SBATCH --partition=slow

srun /home/jsa306/assignment2/ellipse_area_parallel --nThreads 37 --nPoints 2000000  --majorRad 2.0 --minorRad 1.5 --rSeed 0
