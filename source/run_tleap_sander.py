import os, sys
import StringIO
from libtbx import easy_run

tleap_in = """source leaprc.gaff
x = loadMol2 {0}.mol2
loadAmberParams {0}.frcmod
saveAmberParm x {0}.prmtop {0}.rst7
quit
"""

sander_in = """  short minimization
 &cntrl
   imin=1, maxcyc=100, ncyc=50, ntpr=1, cut=99., ntb=0,
 /
"""

def run_tleap(code, igb=False):
  inl = tleap_in.format(code)
  if igb:
    inl = inl.replace("loadAmber", "set default PBradii mbondi2\nloadAmber")
  f=file("%s_tleap.in" % code, "wb")
  f.write(inl)
  f.close()
  cmd='tleap -f {0}_tleap.in'.format(code)
  print cmd
  ero = easy_run.fully_buffered(command=cmd)
  std = StringIO.StringIO()
  ero.show_stdout(out=std)
  outl = ""
  for line in std.getvalue().split("\n"):
    outl += "%s\n" % line
  print outl
  #os.remove("%s_tleap.in" % code)
  os.remove("leap.log")
  return 0

def run_sander(code, igb=False):
  sin = sander_in
  if igb:
    sin = sin.replace("ntb=0,", "ntb=0, igb=2")
  f=file("%s_sander.in" % code, "wb")
  f.write(sin)
  f.close()
  cmd='sander -O -i {0}_sander.in'.format(code)
  cmd+=' -p {0}.prmtop'.format(code)
  cmd+=' -c {0}.rst7'.format(code)
  if igb:
    cmd+=' -o {0}.min_igb.out'.format(code)
    cmd+=' -r {0}.min_igb.rst7'.format(code)
  else:
    cmd+=' -o {0}.min.out'.format(code)
    cmd+=' -r {0}.min.rst7'.format(code)
  print cmd
  ero = easy_run.fully_buffered(command=cmd)
  std = StringIO.StringIO()
  ero.show_stdout(out=std)
  outl = ""
  for line in std.getvalue().split("\n"):
    outl += "%s\n" % line
  print outl
  #os.remove("%s_sander.in" % code)
  #os.remove("mdinfo")
  return 0

def run_ambpdb(code, igb=False):
  cmd='ambpdb'
  cmd+=' -p {0}.prmtop'.format(code)
  if igb:
    cmd+=' < {0}.min_igb.rst7'.format(code)
    cmd+=' > {0}.min_igb.pdb'.format(code)
  else:
    cmd+=' < {0}.min.rst7'.format(code)
    cmd+=' > {0}.min.pdb'.format(code)
  print cmd
  ero = easy_run.fully_buffered(command=cmd)
  std = StringIO.StringIO()
  ero.show_stdout(out=std)
  outl = ""
  for line in std.getvalue().split("\n"):
    outl += "%s\n" % line
  print outl
  
def run(only_code):
  print 'only_code',only_code
  os.chdir(only_code[0])
  assert os.path.exists("%s.mol2" % only_code)
  for igb in range(2):
    run_tleap(only_code, igb=igb)
    run_sander(only_code, igb=igb)
    run_ambpdb(only_code, igb=igb)
  
if __name__=="__main__":
  run(sys.argv[1])
  
