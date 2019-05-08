#!/bin/bash

#SBATCH -n 1
#SBATCH -p pml
#SBATCH -t 55:00:00
#SBATCH -J 07_2
#SBATCH -o 07_2-%j.out
#SBATCH -e 07_2-%j.err
#SBATCH --export=ALL


/nist/apps/matlab-r2018b/bin/matlab  -nodisplay -nojvm < kick_evolve_with_soc_cluster07_2.m



