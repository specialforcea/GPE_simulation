#!/bin/bash

#SBATCH -n 1
#SBATCH -p pml
#SBATCH -t 20:00:00
#SBATCH -J 7
#SBATCH -o 7-%j.out
#SBATCH -e 7-%j.err
#SBATCH --export=ALL


/nist/apps/matlab-r2018b/bin/matlab  -nodisplay -nojvm < kick_and_evolve_cluster7.m