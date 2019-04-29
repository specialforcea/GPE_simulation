#!/bin/bash

#SBATCH -n 1
#SBATCH -p pml
#SBATCH -t 20:00:00
#SBATCH -J 1
#SBATCH -o 1-%j.out
#SBATCH -e 1-%j.err
#SBATCH --export=ALL


/nist/apps/matlab-r2018b/bin/matlab  -nodisplay -nojvm < kick_and_evolve_cluster1.m