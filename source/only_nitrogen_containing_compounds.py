import os, sys
from mmtbx.chemical_components import get_smiles

def main(code):
  smiles = get_smiles(code)
  if smiles.find('N')>-1:
    print 'True'
  else:
    print 'False'

if __name__ == '__main__':
  main(sys.argv[1])
