#!/usr/bin/env python

'''compare missing mol2 file between old run (gaff) and new run (gaff2, Hai run)
'''

dirs = '0 1 2 3 4 5 6 7 8 9 a b c d e f g h i j k l m n o p q r s t u v w x y z'.split()

missing_mol2_file = 'missing_mol2.dat' 
case_dir = '/home/case/amber_library/'
hai_dir = '/home/haichit/research/phenix_amber/amber_library/'  

def get_missing_mol2s(fn):
    with open(fn) as fh:
        return {line.strip('\n').split('/')[-1] for line in fh.readlines()}

def get_union_of_mol2(case_dir, hai_dir, missing_mol2):
    mol2_missing_gaff_run = get_missing_mol2s(case_dir + missing_mol2_file)
    mol2_missing_gaff2_run = get_missing_mol2s(hai_dir + missing_mol2_file)
    return mol2_missing_gaff_run |  mol2_missing_gaff2_run

if __name__ == '__main__':
    hop_gaff_gaff2 = get_union_of_mol2(case_dir, hai_dir, missing_mol2_file)

    for i, code in enumerate(hop_gaff_gaff2):
        print(i ,code)
