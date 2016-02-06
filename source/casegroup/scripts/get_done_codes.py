#!/usr/bin/env python
from glob import glob

done_codes = [fn.split('.')[-2] for fn in glob('./amber.*.output')]

for code in done_codes:
    print(code)
