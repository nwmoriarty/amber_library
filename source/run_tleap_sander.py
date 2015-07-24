import os, sys

def run(only_code):
  print 'only_code',only_code
  os.chdir(only_code[0])
  assert os.path.exists("%s.mol2" % only_code)
  cmd = "../tleap_sander.sh %s" % only_code
  print cmd
  os.system(cmd)
  cmd = "../tleap_sander_igb.sh %s" % only_code
  print cmd
  os.system(cmd)
  
if __name__=="__main__":
  run(sys.argv[1])
  
