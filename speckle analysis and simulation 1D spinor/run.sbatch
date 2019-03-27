#!/bin/bash

#SBATCH -n 1
#SBATCH -p pml
#SBATCH -t 01:00:00
#SBATCH -J 0510
#SBATCH -o 0510-%j.out
#SBATCH -e 0510-%j.err
#SBATCH --export=ALL


/nist/apps/matlab-r2018b/bin/matlab  -nodisplay -nojvm < SOC_dressed_state_cluster10.m