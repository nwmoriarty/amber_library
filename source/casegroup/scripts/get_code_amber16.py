#!/usr/bin/env python
from glob import iglob
import subprocess
import linecache

"""print all unminized ligand code with amber16 (using gaff2)
"""

with open('../amber_library/source/casegroup/ligand_codes.dat') as fh:
    all_codes = {line.split()[-1] for line in fh.readlines()}

outfiles = iglob('../amber_library/*/*.min_igb.out')

def has_amber16(fn):
    # dummy search
    value = 'haichit' in linecache.getline(fn, 9)
    return value
 
finished_codes = {fn.split('.')[-3].split('/')[-1] for fn in outfiles if has_amber16(fn)}

for code in (all_codes ^ finished_codes):
    print(code)
