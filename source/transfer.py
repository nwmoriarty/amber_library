import os, sys

def main(code):
  print 'main',code
  run_filename = os.path.realpath(__file__)
  working = os.path.dirname(os.path.dirname(run_filename))
  print working
  mol2_relative_path = os.path.join(code[0].lower(),'%s.mol2' % code.upper())
  mol2 = os.path.join(working, mol2_relative_path)
  print mol2
  if not os.path.exists(mol2): assert 0
  frcmod_relative_path = os.path.join(code[0].lower(),'%s.frcmod' % code.upper())
  frcmod = os.path.join(working, frcmod_relative_path)
  print frcmod
  if not os.path.exists(frcmod): assert 0
  #
  #
  #
  if not os.path.exists(code[0].lower()): os.mkdir(code[0].lower())
  os.rename(mol2, mol2_relative_path)
  os.rename(frcmod, frcmod_relative_path)

if __name__ == '__main__':
  main(*tuple(sys.argv[1:]))
