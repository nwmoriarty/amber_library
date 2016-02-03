# use phenix.python

from __future__ import print_function
import os
from glob import glob
from itertools import chain

# from mmtbx.chemical_components import generate_chemical_components_codes
# codes = list(generate_chemical_components_codes())

dirs = '0 1 2 3 4 5 6 7 8 9 a b c d e f g h i j k l m n o p q r s t u v w x y z'.split()

# get all ligand codes thave have mole file
codes = [item.split('/')[-1].split('.')[0] 
         for letter in dirs 
         for item in glob(letter + '/*.mol2')]

with open('ligand_codes.dat', 'w') as fh:
    for index, code in enumerate(codes):
        fh.write('{index} {code} \n'.format(index=index, code=code))
