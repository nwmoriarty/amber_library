#!/usr/bin/env python
import os
import subprocess
from mpi4py import MPI
comm = MPI.COMM_WORLD

rank, n_cores = comm.rank, comm.size

def run_each_core(ligand_codes, amber_library_path, template, tleap=False):
    '''run a chunk of ligand codes in each core

    Parameters
    ----------
    ligand_codes : list of single code (e.g. ['OTT', ])
    run_elbow_template : str, bash template to run elbow
    amber_library_path : str, absolute path of amber_library
    '''
    for code in ligand_codes:
        print('running {}'.format(code))
        if not tleap:
            run_elbow = (template
                         .strip()
                         .format(amber_library=amber_library_path, code=code) 
                        )
            subprocess.call(' '.join(run_elbow.split('\n')), shell=True)
        else:
            run_elbow = (template
                         .format(amber_library=amber_library_path, code=code) 
                        )
            subprocess.call(' '.join(run_elbow.split('\n')), shell=True)
        

def get_codes_for_my_rank(ligand_codes, rank):
    import numpy as np
    return np.array_split(ligand_codes, n_cores)[rank]

def get_ligand_codes():
    '''use env LIGAND_CODES. If not, use ./amber_library/source/casegroup/ligand_codes.dat

    Returns
    -------
    ligands : list of str
    '''
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

def main_elbow():
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

    n_codes = len(partial_codes)
    x = comm.gather(n_codes, root=0)
    if rank == 0:
        print('max n_codes per node = {}'.format(max(x)))
        print('min n_codes per node = {}'.format(min(x)))
    run_each_core(partial_codes, amber_library, template=run_elbow_template)

def main_tleap():
    run_elbow_template = '''
export source=`pwd`/amber_library/source/;
export opwd=`pwd`;

echo {code};
mycode=`elbow.python $source/casegroup/get_code.py output {code}`;
echo $mycode;
if [ $mycode ]; then
    cd {amber_library};
    elbow.python $source/casegroup/run_tleap_sander.py $mycode --force >& ../output/tleap.$mycode.output;
    cd $opwd;
fi
    '''
    cwd = os.getcwd()
    amber_library = os.path.join(cwd, 'amber_library')
    ligand_codes = get_ligand_codes()
    partial_codes = get_codes_for_my_rank(ligand_codes, rank)

    n_codes = len(partial_codes)
    x = comm.gather(n_codes, root=0)
    if rank == 0:
        print('max n_codes per node = {}'.format(max(x)))
        print('min n_codes per node = {}'.format(min(x)))
    run_each_core(partial_codes, amber_library, template=run_elbow_template, tleap=True)


if __name__ == '__main__':
    '''How? 
        mpirun -n 24 python this_script.py --elbow
        mpirun -n 24 python this_script.py --tleap
    '''
    import sys

    if len(sys.argv) != 2:
        if rank == 0:
            raise ValueError("Must follow: python {} --elbow (or --tleap)".format(sys.argv[0]))
    if sys.argv[-1] == '--elbow':
        main_elbow()
    elif sys.argv[-1] == '--tleap':
        main_tleap()
    else:
        if rank == 0:
            raise ValueError('must using --elbow or --tleap options')
