#!/bin/bash

#SBATCH -n 1
#SBATCH -p pml
#SBATCH -t 20:00:00
#SBATCH -J 4
#SBATCH -o 4-%j.out
#SBATCH -e 4-%j.err
#SBATCH --export=ALL


/nist/apps/matlab-r2018b/bin/matlab  -nodisplay -nojvm < kick_and_evolve_cluster4.m