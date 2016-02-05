#!/usr/bin/env python
import os
import subprocess
from mpi4py import MPI
comm = MPI.COMM_WORLD

rank, n_cores = comm.rank, comm.size

def split_range(n_chunks, start, stop):
    '''split a given range to n_chunks

    Examples
    --------
    >>> split_range(3, 0, 10)
    [(0, 3), (3, 6), (6, 10)]
    '''
    list_of_tuple = []
    chunksize = (stop - start) // n_chunks
    for i in range(n_chunks):
        if i < n_chunks - 1:
            _stop = start + (i + 1) * chunksize
        else:
            _stop = stop
        list_of_tuple.append((start + i * chunksize, _stop))
    return list_of_tuple

def run_each_core(ligand_codes):
    '''run a chunk of ligand codes in each core
    '''
    for code in partial_codes:
        run_elbow = (run_elbow_template
                     .strip()
                     .format(amber_library=amber_library, code=code) 
                    )
        
        subprocess.call(' '.join(run_elbow.split('\n')), shell=True)

def get_codes_for_my_rank(ligand_codes, rank):
    n_ligands = len(ligand_codes)
    start, stop = split_range(n_cores, 0, n_ligands)[rank]
    return ligand_codes[start: stop]


if __name__ == '__main__':
    '''How? mpirun -n 24 python this_script.py
    '''

    run_elbow_template = '''
    elbow.python
        {amber_library}/source/generate_all_chemical_component_restraint_files.py
        amber=True
        ignore_output_files=True
        skip_ligands_in_library=False
        only_type=NON-POLYMER
        pH=8
        only_code={code}
        only_i=None >& output/amber.{code}.output
    '''
    
    cwd = os.getcwd()
    amber_library = os.path.join(cwd, 'amber_library')
    ligand_list = os.environ.get('LIGAND_CODES', amber_library + '/source/casegroup/ligand_codes.dat')
    
    with open(ligand_list, 'r') as fh:
        ligand_codes = [line.split()[-1] for line in fh.readlines()]

    partial_codes = get_codes_for_my_rank(ligand_codes, rank)
    run_each_core(partial_codes)

