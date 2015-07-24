import os, sys
import pickle
import time
from math import sqrt

from multiprocessing import Pool

import iotbx.pdb

from elbow.utilities import mogul_utils, rmsd_utils
from elbow.chemistry.any_chemical_format_reader import \
  any_chemical_format_reader
from elbow.scripts.load_and_test_gaussian_files import validate_elbow_object

results = []

skip = {
  "amber_library" : [
    #"0EY",
    "0OC",
    "08S",
    #"0UC",
    "1BY",
    "LCO",
    "511",
    "9D9",
    "AXZ",
    "A7D",
    "DWN",
    "T07",
    "TSD",
    ],
  }
skip = {}

def results_callback(args):
#  print '\n\n\targs',args
  sys.stdout.flush()
  results.append(args)

def _is_newer(f1, f2):
  if not os.path.exists(f1): return False
  if not os.path.exists(f2): return False
  if os.stat(f1).st_mtime>os.stat(f2).st_mtime: return True
  return False

def mogul_process_filename(filename):
  print 'mogul_process_filename',filename
  if filename is None: return None
  preamble = filename.replace(".pdb", "")
  pf="%s_mogul.pickle" % preamble
  if _is_newer(pf,filename):
    if os.path.exists(pf): return ["mogul", filename]
  #
  if not os.path.exists(filename): return ["failed", filename]
  molecule = any_chemical_format_reader(filename,
                                        simple=False,
                                      )
  veo = validate_elbow_object(molecule)
  if 0 and filter(None, veo.values()):
    if os.path.exists(pf): os.remove(pf)
    return ["invalid", filename]
  #
  if _is_newer(pf,filename):
    if os.path.exists(pf): return ["pickle",filename]
  # optimise H
  format_str = "%s"
  if 0:
    format_str = "%s_minimized"
    code = os.path.basename(preamble)
    code = code.split(".")[0]
    print 'code',code
    mol = any_chemical_format_reader("%s.pdb" % preamble, simple=False)
    mol.UpdateBonded()
    for bond in mol.bonds:
      if bond[0].isH() or bond[1].isH():
        print bond
        bond.mogul_value=0.98
    print mol
    f=file("%s.cif" % code, "wb")
    f.write(mol.WriteCIF2String(ligand_code=code))
    f.close()
    print code
    cmd = 'phenix.geometry_minimization'
    cmd += ' selection="element H or element D"'
    cmd += ' %s.pdb' % preamble
    cmd += ' %s.cif' % code
    cmd += ' output=%s_minimized' % preamble
    print cmd
    os.system(cmd)
  #
  t0=time.time()
  rc = mogul_utils.run_mogul(format_str % preamble,
                             format="pdb",
                             defaults=["ON", 25, 25, 60, 25],
                             strict=False,
                             generate_hydrogens=False,
                             verbose=True,
                             )
  print 'Mogul took %0.1fs' % (time.time()-t0)
  if os.path.exists("%s_minimized.pdb" % preamble):
    os.remove("%s_minimized.pdb" % preamble)
  if type(rc)==type([]):
    print "\n\tMogul doesn't appear to be installed"
    print "\nOutput from Mogul command\n"
    for line in mogul_object:
      print line
    return ["failed", filename]

  f=file(pf, "wb")
  pickle.dump(rc, f)
  f.close()
  return ["mogul", filename]
    
def run(only_code=None,
        cpus=1,
        ):
  if only_code.lower()=="none":
    only_code=None
  try: cpus=int(cpus)
  except: cpus=1

  print 'only_code',only_code
  print 'cpus',cpus

  pool = None
  if cpus>1:
    pool = Pool(processes=cpus)

  if 1:
    for d in sorted(os.listdir(os.getcwd())):
      if len(d)!=1: continue
      if d!="0": continue
      for i, filename in enumerate(sorted(
          os.listdir(os.path.join(os.getcwd(), d)))):
        #if filename.find("min")==-1: continue
        if filename.find(".pdb")==-1: continue
        print i, filename
        if only_code is not None and filename.find(only_code.upper())!=0:
          continue
        filename = os.path.join(os.getcwd(), d, filename)
        print i,filename,len(results)
        #if len(results)>20: break
        if pool:
          rc = pool.apply_async(
            mogul_process_filename,
            [filename],
            callback=results_callback,
            )
        else:
          rc = mogul_process_filename(filename)
          print rc
          #assert rc[0]!="mogul"

        #if only_code is not None: break

    if pool:
      pool.close()
      pool.join()
      print '\nProcesses have joined\n'


  print results
  if only_code is not None: return

  df = "data.pickle"
  try:
    f=file(df, "rb")
    data = pickle.load(f)
    f.close()
  except:
    data = {}
  data.setdefault("zfiles", {})
  for d in sorted(os.listdir(os.getcwd())):
    if len(d)!=1: continue
    os.chdir(d)
    for filename in os.listdir(os.getcwd()):
      if filename.find("mogul")==-1: continue
      if filename.find(".pickle")==-1: continue
      print filename
      if filename in data["zfiles"]: continue
      key = "vacuum"
      if filename.find("igb")>-1:
        key = "solvent"
      data.setdefault(key, {})
      f=file(filename, "rb")
      mo = pickle.load(f)
      f.close()
      if type(mo)==type([]): continue
      rc = mo.squared_deviations_and_numbers()
      
      for frag in rc:
        if rc[frag][1]==0: continue
        if rc[frag][3]==0: continue
        if rc[frag][4]==0: continue
        msd = rc[frag][0]
        try:
          msd = rc[frag][0]/rc[frag][1]
        except:
          pass
        tmp = []
        tmp.append(sqrt(msd))
        tmp.append(rc[frag][1]/rc[frag][4])
        msd = rc[frag][2]
        try:
          msd = rc[frag][2]/rc[frag][3]
        except:
          pass
        tmp.append(sqrt(msd))
        tmp.append(rc[frag][3]/rc[frag][4])
        print 'tmp',tmp
        data[key].setdefault(frag, [])
        data[key][frag].append(tmp)
        data["zfiles"][filename] = None
      
    os.chdir("..")

  f=file(df, "wb")
  pickle.dump(data, f)
  f.close()

  outl = {}
  histo = {}
  for i, enviro in enumerate(sorted(data)):
    if enviro in ["zfiles"]: continue
    histo.setdefault(enviro, {})
    for frag in data[enviro]:
      outl.setdefault(frag, "")
      histo[enviro].setdefault(frag, [])
      outl[frag] += "@target G0.S%d\n@type xy\n" % i
      for pair in data[enviro][frag]:
        outl[frag] += " %s %s\n" % (pair[0], pair[1])
        histo[enviro][frag].append(pair[0])
      outl[frag] += "&\n"
      print outl[frag]

  for frag in outl:
    to = ""
    for i, enviro in enumerate(sorted(histo)):
      if enviro in ["zfiles"]: continue
      to += "@target G1.S%d\n&type xy\n" % i
      if 1:
        for pair in data[enviro][frag]:
          to += " %s %s\n" % (pair[1], pair[0])
        to += "&\n"
      else:
        th = rmsd_utils.histogram(histo[enviro][frag],
                                  bins=25,
          )
        for t in sorted(th):
          to += " %s %s\n" % (t, th[t])
        to += "&\n"
      m = rmsd_utils.mean(histo[enviro][frag])
    #print '-'*80
    #print outl[frag]
    print '-'*80
    print to
    gf = "%s.dat" % (frag)
    f=file(gf, "wb")
    f.write(outl[frag])
    f.write(to)
    f.close()
    print 'mean',m

    cmd = "xmgrace -param ../test.par -nxy %s" % gf
    print cmd
    os.system(cmd)


if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
