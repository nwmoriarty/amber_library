import os, sys
from mmtbx.chemical_components import get_smiles

def main(code):
  p = os.path.join('amber_library',
                   code[0].lower(),
                   '%s.mol2' % code)
  if os.path.exists(p):
    print 'False'
  else:
    print 'True'

if __name__ == '__main__':
  main(sys.argv[1])
