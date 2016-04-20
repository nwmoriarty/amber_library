#!/bin/sh
#SBATCH -J tleap_sander 
#SBATCH -o log/tleap_sander%J.stdout
#SBATCH -e log/tleap_sander.%J.stderr
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 4:00:00

export LIGAND_CODES=./amber_library/source/casegroup/ligand_codes.dat
mpirun -n 24 python  ./amber_library/source/casegroup/phenix_mpi.py --tleap
