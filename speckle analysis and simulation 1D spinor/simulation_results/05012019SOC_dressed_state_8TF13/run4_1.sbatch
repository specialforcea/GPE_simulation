#!/bin/bash

#SBATCH -n 1
#SBATCH -p pml
#SBATCH -t 48:00:00
#SBATCH -J 4_1
#SBATCH -o 4_1-%j.out
#SBATCH -e 4_1-%j.err
#SBATCH --export=ALL


/nist/apps/matlab-r2018b/bin/matlab  -nodisplay -nojvm < SOC_dressed_state_cluster4_1.m