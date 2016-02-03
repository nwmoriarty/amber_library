#!/bin/sh
#SBATCH -J amberphenix
#SBATCH -o log/amberphenix.%J.stdout
#SBATCH -e log/amberphenix.%J.stderr
#SBATCH -p long
#SBATCH -N 4
#SBATCH -t 96:00:00

# require: mpi4py
# (conda install mpi4py)

cd /project1/dacase-001/haichit/phenix_amber/run_amber_library

# Note: You need to submit this file to cluster in the folder having
# - ./amber_library
# - ./output/

# RUN:
# sbatch amber_library/source/casegroup/submit_job.sh

# this program split ~15K ligands to n_cores
# each core run a chunk of ligands
# 48 cores * 3 days should be sufficient to finish.

# if you want to update ligand code, edit: amber_library/source/casegroup/ligand_codes.dat

mpirun -n 96 python ./amber_library/source/casegroup/phenix_mpi.py
