#!/usr/bin/env python
from glob import glob
import subprocess

done_codes_sqm = {fn.split('.')[1] for fn in glob('amber.*.output')}
done_codes_tleap = {fn.split('.')[1] for fn in glob('tleap.*.output')}

with open('../amber_library/source/casegroup/ligand_codes.dat') as fh:
    all_codes = {line.split()[-1] for line in fh.readlines()}

def get_codes(all_codes, done_codes):
    for idx, code in enumerate(all_codes - done_codes):
        # take all not-run-yet codes
        # print(idx, code)
        print(code)

# get_codes(all_codes, done_codes_sqm)
get_codes(all_codes, done_codes_tleap)

unminized_codes = set(all_codes) - set(done_codes_tleap)
print('unminized_codes', unminized_codes)
