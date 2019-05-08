#!/bin/bash

#SBATCH -n 1
#SBATCH -p pml
#SBATCH -t 48:00:00
#SBATCH -J 6_6
#SBATCH -o 6_6-%j.out
#SBATCH -e 6_6-%j.err
#SBATCH --export=ALL


/nist/apps/matlab-r2018b/bin/matlab  -nodisplay -nojvm < SOC_dressed_state_cluster6_6.m