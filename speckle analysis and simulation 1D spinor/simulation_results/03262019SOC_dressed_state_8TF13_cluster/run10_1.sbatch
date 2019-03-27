#!/bin/bash

#SBATCH -n 1
#SBATCH -p pml
#SBATCH -t 10:00:00
#SBATCH -J 10_1
#SBATCH -o 10_1-%j.out
#SBATCH -e 10_1-%j.err
#SBATCH --export=ALL


/nist/apps/matlab-r2018b/bin/matlab  -nodisplay -nojvm < SOC_dressed_state_cluster10_1.m