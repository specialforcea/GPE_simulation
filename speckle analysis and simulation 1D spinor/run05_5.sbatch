#!/bin/bash


#SBATCH -n 1
#SBATCH -p pml
#SBATCH -t 55:00:00
#SBATCH -J 05_5
#SBATCH -o 05_5-%j.out
#SBATCH -e 05_5-%j.err
#SBATCH --export=ALL


/nist/apps/matlab-r2018b/bin/matlab  -nodisplay -nojvm < kick_evolve_with_soc_cluster05_5.m



