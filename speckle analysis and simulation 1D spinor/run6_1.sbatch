#!/bin/bash

#SBATCH -n 1
#SBATCH -p pml
#SBATCH -t 55:00:00
#SBATCH -J 6_1
#SBATCH -o 6_1-%j.out
#SBATCH -e 6_1-%j.err
#SBATCH --export=ALL


/nist/apps/matlab-r2018b/bin/matlab  -nodisplay -nojvm < kick_evolve_with_soc_cluster6_1.m



