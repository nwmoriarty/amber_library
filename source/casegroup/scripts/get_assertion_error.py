#!/usr/bin/env python
from glob import glob
import subprocess

def get():
    for fn in glob('amber.*.output'):
        code = fn.split('.')[1]
        with open(fn) as fh:
            for line in fh.readlines():
                if 'Assertion' in line:
                    yield code

print(list(get()))
