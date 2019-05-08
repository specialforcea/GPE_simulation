#!/bin/bash

#SBATCH -n 1
#SBATCH -p pml
#SBATCH -t 55:00:00
#SBATCH -J 05_1
#SBATCH -o 05_1-%j.out
#SBATCH -e 05_1-%j.err
#SBATCH --export=ALL


/nist/apps/matlab-r2018b/bin/matlab  -nodisplay -nojvm < kick_evolve_with_soc_cluster05_1.m

#!/bin/bash


#SBATCH -n 1
#SBATCH -p pml
#SBATCH -t 55:00:00
#SBATCH -J 05_2
#SBATCH -o 05_2-%j.out
#SBATCH -e 05_2-%j.err
#SBATCH --export=ALL


/nist/apps/matlab-r2018b/bin/matlab  -nodisplay -nojvm < kick_evolve_with_soc_cluster05_2.m

#!/bin/bash


#SBATCH -n 1
#SBATCH -p pml
#SBATCH -t 55:00:00
#SBATCH -J 05_3
#SBATCH -o 05_3-%j.out
#SBATCH -e 05_3-%j.err
#SBATCH --export=ALL


/nist/apps/matlab-r2018b/bin/matlab  -nodisplay -nojvm < kick_evolve_with_soc_cluster05_3.m

#!/bin/bash


#SBATCH -n 1
#SBATCH -p pml
#SBATCH -t 55:00:00
#SBATCH -J 05_4
#SBATCH -o 05_4-%j.out
#SBATCH -e 05_4-%j.err
#SBATCH --export=ALL


/nist/apps/matlab-r2018b/bin/matlab  -nodisplay -nojvm < kick_evolve_with_soc_cluster05_4.m

#!/bin/bash


#SBATCH -n 1
#SBATCH -p pml
#SBATCH -t 55:00:00
#SBATCH -J 05_5
#SBATCH -o 05_5-%j.out
#SBATCH -e 05_5-%j.err
#SBATCH --export=ALL


/nist/apps/matlab-r2018b/bin/matlab  -nodisplay -nojvm < kick_evolve_with_soc_cluster05_5.m