#!/bin/bash

#SBATCH -n 1
#SBATCH -p pml
#SBATCH -t 10:00:00
#SBATCH -J 10_2
#SBATCH -o 10_2-%j.out
#SBATCH -e 10_2-%j.err
#SBATCH --export=ALL


/nist/apps/matlab-r2018b/bin/matlab  -nodisplay -nojvm < SOC_dressed_state_cluster10_2.m