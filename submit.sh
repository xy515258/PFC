#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH -t 1:00:00 
#SBATCH --mail-user=yangx@princeton.edu 

# Load openmpi environment 
srun ./pfcpure3d


