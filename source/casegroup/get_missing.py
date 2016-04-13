#!/usr/bin/env python

fn = 'missing_min_pdb.dat'
fn2 = '/home/case/amber_library/missing_min_pdb.dat'

def get_missing(fn):
    with open(fn) as fh:
        return set(line.split()[0] for line in fh.readlines())

print(get_missing(fn))
print(get_missing(fn2))
print(get_missing(fn) - get_missing(fn2))
