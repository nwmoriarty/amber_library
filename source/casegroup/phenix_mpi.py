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

def run_each_core(ligand_codes, run_elbow_template, amber_library_path):
    '''run a chunk of ligand codes in each core

    Parameters
    ----------
    ligand_codes : list of single code (e.g. ['OTT', ])
    run_elbow_template : str, bash template to run elbow
    amber_library_path : str, absolute path of amber_library
    '''
    for code in partial_codes:
        run_elbow = (run_elbow_template
                     .strip()
                     .format(amber_library=amber_library_path, code=code) 
                    )
        
        subprocess.call(' '.join(run_elbow.split('\n')), shell=True)

def get_codes_for_my_rank(ligand_codes, rank):
    n_ligands = len(ligand_codes)
    start, stop = split_range(n_cores, 0, n_ligands)[rank]
    return ligand_codes[start: stop]

def get_ligand_codes():
    '''use env LIGAND_CODES. If not, use ./amber_library/source/casegroup/ligand_codes.dat

    Returns
    -------
    ligands : list of str
    '''
    # ligand_fn = os.environ.get('LIGAND_CODES', amber_library + '/source/casegroup/ligand_codes.dat')
    ligand_fn = os.environ.get('LIGAND_CODES', '')
    if not ligand_fn:
      raise RuntimeError('you must export LIGAND_CODES')
      sys.exit(1)
     
    if rank == 0:
        print("using ligand codes from ", ligand_fn)
    
    with open(ligand_fn, 'r') as fh:
        ligand_codes = []
        for line in fh.readlines():
            if '/' not in line:
                ligand_codes.append(line.strip().split()[-1])
            else:
                ligand_codes.append(line.strip().split('/')[-1])
    return ligand_codes


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
    ligand_codes = get_ligand_codes()
    partial_codes = get_codes_for_my_rank(ligand_codes, rank)
    run_each_core(partial_codes, run_elbow_template, amber_library)

