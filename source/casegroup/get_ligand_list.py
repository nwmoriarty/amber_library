#!/usr/bin/env phenix.python

from __future__ import print_function
import os
from glob import glob
from itertools import chain
import pytraj as pt
import parmed as pmd


cwd = os.getcwd()

root_dir = cwd.split('amber_library')[0] + '/amber_library/'
print('root_dir = {}'.format(root_dir))

dirs = [root_dir + letter
        for letter in '0 1 2 3 4 5 6 7 8 9 a b c d e f g h i j k l m n o p q r s t u v w x y z'.split()]

# get all ligand codes thave have mol2 file
def get_codes(dirs):
    for letter in dirs:
        for item in glob(letter + '/*.mol2'):
            try:
                parm = pmd.load_file(item)
            except pmd.exceptions.Mol2Error:
                print(item)
            # parm = pt.load_topology(item)
            # if 'B' not in set(atom.name for atom in parm.atoms):
            #     code = item.split('/')[-1].split('.')[0]
            #     # yield code
            #     print(code)

def write(codes, filename='ligand_codes.dat'):
    with open(filename, 'w') as fh:
        for index, code in enumerate(codes):
            fh.write('{code} \n'.format(code=code))

if __name__ == '__main__':
    get_codes(dirs)
