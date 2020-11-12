#!/bin/bash
#SBATCH --job-name=test       # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=20        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G per cpu-core is default)
#SBATCH --time=01:00:00          # total run time limit (HH:MM:SS)
# SBATCH --mail-type=all          # send email on job start, end and fail
# SBATCH --mail-user=yangx@princeton.edu

cd .

module purge
module load /projects/SOFTWARE/Modules/modulefiles/tigressdata/paraview-headless/5.8.1-py3.7

mpiexec -n 20 pvbatch animation.py
# pvpython plot.py
