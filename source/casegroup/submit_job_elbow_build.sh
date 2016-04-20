#!/bin/sh
#SBATCH -J amberphenix
#SBATCH -o log/amberphenix.%J.stdout
#SBATCH -e log/amberphenix.%J.stderr
#SBATCH -p development
#SBATCH -N 1
#SBATCH -t 2:00:00

# require: mpi4py
# (conda install mpi4py)

cd /project1/dacase-001/haichit/phenix_amber/run_amber_library

# if there is no LIGAND_CODES given, phenix_amber.py will
# use ligand_codes.dat in ./amber_library/source/casegroup/

# export LIGAND_CODES=./missing_ligands.dat
export LIGAND_CODES=./amber_library/missing_mol2.dat

# Note: You need to submit this file to cluster in the folder having
# - ./amber_library
# - ./output/

# RUN:
# sbatch amber_library/source/casegroup/submit_job.sh

# this program split ~15K ligands to n_cores
# each core run a chunk of ligands
# 48 cores * 3 days should be sufficient to finish.

# mpirun -n 24 python ./amber_library/source/casegroup/phenix_mpi.py
mpirun -n 24 python  ./amber_library/source/casegroup/phenix_mpi.py --elbow
