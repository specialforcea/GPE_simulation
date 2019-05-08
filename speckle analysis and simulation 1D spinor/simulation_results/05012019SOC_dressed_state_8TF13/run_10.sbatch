#!/bin/bash

#SBATCH -n 1
#SBATCH -p pml
#SBATCH -t 48:00:00
#SBATCH -J 10
#SBATCH -o 10-%j.out
#SBATCH -e 10-%j.err
#SBATCH --export=ALL


/nist/apps/matlab-r2018b/bin/matlab  -nodisplay -nojvm < SOC_dressed_state_cluster_10.m