#!/bin/bash


#SBATCH -n 1
#SBATCH -p pml
#SBATCH -t 55:00:00
#SBATCH -J 1_4
#SBATCH -o 1_4-%j.out
#SBATCH -e 1_4-%j.err
#SBATCH --export=ALL


/nist/apps/matlab-r2018b/bin/matlab  -nodisplay -nojvm < kick_evolve_with_soc_cluster1_4.m



