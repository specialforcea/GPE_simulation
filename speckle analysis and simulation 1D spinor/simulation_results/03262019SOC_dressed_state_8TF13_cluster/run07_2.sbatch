#!/bin/bash

#SBATCH -n 1
#SBATCH -p pml
#SBATCH -t 10:00:00
#SBATCH -J 07_2
#SBATCH -o 07_2-%j.out
#SBATCH -e 07_2-%j.err
#SBATCH --export=ALL


/nist/apps/matlab-r2018b/bin/matlab  -nodisplay -nojvm < SOC_dressed_state_cluster07_2.m