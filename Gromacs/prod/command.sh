#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --account=def-bingalls
#SBATCH --mail-user=mrastwoo@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=16

module load  gcc/7.3.0 openmpi/3.1.2 gromacs/2019.3
gmx mdrun -deffnm nvt