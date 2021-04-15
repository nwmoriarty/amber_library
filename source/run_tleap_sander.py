import os, sys
import StringIO
from libtbx import easy_run

tleap_in = """source leaprc.gaff2
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

def get_binary(s):
  return os.path.abspath(os.path.join(os.environ["LIBTBX_BUILD"],'..','conda_base','bin',s))

def files_exist_tleap(code, igb=False):
  for filename in ["%s.prmtop" % code,
                   "%s.rst7" % code,
                   ]:
    if not os.path.exists(filename):
      return False
  return True

def run_tleap(code, igb=False):
  inl = tleap_in.format(code)
  if igb:
    inl = inl.replace("loadAmber", "set default PBradii mbondi2\nloadAmber")
  f=file("%s_tleap.in" % code, "wb")
  f.write(inl)
  f.close()
  cmd = '{0} -f {1}_tleap.in'.format(get_binary('tleap'), code)
  print cmd
  ero = easy_run.fully_buffered(command=cmd)
  std = StringIO.StringIO()
  ero.show_stdout(out=std)
  outl = ""
  for line in std.getvalue().split("\n"):
    outl += "%s\n" % line
  print outl
  for filename in ["%s_tleap.in" % code, "leap.log"]:
    if os.path.exists(filename):
      os.remove(filename)
  return 0

def files_exist_sander(code, igb=False):
  for filename in ["%s.min.out" % code,
                   "%s.min.rst7" % code,
                   ]:
    if igb: filename = filename.replace("min", "min_igb")
    if not os.path.exists(filename):
      return False
  return True

def run_sander(code, igb=False):
  sin = sander_in
  if igb:
    sin = sin.replace("ntb=0,", "ntb=0, igb=2")
  f=file("%s_sander.in" % code, "wb")
  f.write(sin)
  f.close()
  cmd='{0} -O -i {1}_sander.in'.format(get_binary('sander'), code)
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
  os.remove("%s_sander.in" % code)
  if os.path.exists("mdinfo"): os.remove("mdinfo")
  return 0

def files_exist_ambpdb(code, igb=False):
  for filename in ["%s.min.pdb" % code,
                   ]:
    if igb: filename = filename.replace("min", "min_igb")
    if not os.path.exists(filename):
      return False
  return True

def run_ambpdb(code, igb=False):
  cmd=get_binary('ambpdb')
  cmd+=' -p {0}.prmtop'.format(code)
  if igb:
    cmd+=' -c {0}.min_igb.rst7'.format(code)
    cmd+=' > {0}.min_igb.pdb'.format(code)
  else:
    cmd+=' -c {0}.min.rst7'.format(code)
    cmd+=' > {0}.min.pdb'.format(code)
  print cmd
  ero = easy_run.fully_buffered(command=cmd)
  std = StringIO.StringIO()
  ero.show_stdout(out=std)
  outl = ""
  for line in std.getvalue().split("\n"):
    outl += "%s\n" % line
  print outl

def run(only_code, force=False):
  only_code=only_code.upper()
  print 'only_code',only_code
  print 'force',force
  os.chdir(only_code.lower()[0])
  assert os.path.exists("%s.mol2" % only_code), '%s %s' % (os.getcwd(),
                                                           only_code)
  for igb in range(2):
    if not files_exist_tleap(only_code, igb=igb) or force:
      print 'run tleap',only_code,igb
      run_tleap(only_code, igb=igb)

    if not files_exist_sander(only_code, igb=igb) or force:
      print 'run sander',only_code,igb
      run_sander(only_code, igb=igb)

    if not files_exist_ambpdb(only_code, igb=igb) or force:
      print 'run ambpdb',only_code,igb
      run_ambpdb(only_code, igb=igb)

if __name__=="__main__":
  run(*tuple(sys.argv[1:]))
