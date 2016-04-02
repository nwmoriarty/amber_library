#!/usr/bin/env python
from glob import iglob
import subprocess
import linecache

"""print all unminized ligand code with amber16 (using gaff2)
"""

def has_amber16(fn):
    # dummy search
    value = 'haichit' in linecache.getline(fn, 10)
    return value
 
outfiles = iglob('ante*.txt')
amber16_codes = {fn.split('.')[0].split('_')[-1] for fn in outfiles if has_amber16(fn)}

outfiles = iglob('ante*.txt')
all_codes = {fn.split('.')[0].split('_')[-1] for fn in outfiles}

print(amber16_codes)
print(len(amber16_codes), len(all_codes))
