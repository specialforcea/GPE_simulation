#!/bin/bash


#SBATCH -n 1
#SBATCH -p pml
#SBATCH -t 55:00:00
#SBATCH -J 1_3
#SBATCH -o 1_3-%j.out
#SBATCH -e 1_3-%j.err
#SBATCH --export=ALL


/nist/apps/matlab-r2018b/bin/matlab  -nodisplay -nojvm < kick_evolve_with_soc_cluster1_3.m



