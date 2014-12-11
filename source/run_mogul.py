import os, sys
import pickle
from math import sqrt

from multiprocessing import Pool

from elbow.utilities import mogul_utils
from elbow.chemistry.any_chemical_format_reader import \
  any_chemical_format_reader
from elbow.scripts.load_and_test_gaussian_files import validate_elbow_object

results = []

def results_callback(args):
  print '\n\n\targs',args
  sys.stdout.flush()
  results.append(args)

def mogul_process_filename(filename):
  print '-'*80
  print 'mogul_process_filename',filename
  preamble = filename.replace(".pdb", "")
  pf="%s_mogul.pickle" % preamble
  #
  molecule = any_chemical_format_reader(filename,
                                        simple=False,
                                      )
  print molecule.DisplayBrief()  
  veo = validate_elbow_object(molecule)
  if filter(None, veo.values()):
    if os.path.exists(pf):
      os.remove(pf)
    return ["invalid", filename]
  #
  if os.path.exists(pf): return ["pickle",filename]
  #
  rc = mogul_utils.run_mogul(preamble,
                             format="pdb",
                             defaults=["ON", 25, 25, 60, 25],
                             strict=False,
                             generate_hydrogens=False,
                             verbose=True,
                             )
  if type(rc)==type([]):
    print "\n\tMogul doesn't appear to be installed"
    print "\nOutput from Mogul command\n"
    for line in mogul_object:
      print line
    return ["failed", filename]
  #print rc
  f=file(pf, "wb")
  pickle.dump(rc, f)
  f.close()
  return ["mogul", filename]

def run(only_code=None):

  cpus = 1
  pool = None
  if cpus>1:
    pool = Pool(processes=cpus)

  if 1:
    for d in sorted(os.listdir(os.getcwd())):
      if len(d)!=1: continue
      for i, filename in enumerate(sorted(
          os.listdir(os.path.join(os.getcwd(), d)))):
        if filename.find("min")==-1: continue
        if filename.find(".pdb")==-1: continue
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
          assert rc[0]!="mogul"

    if pool:
      pool.close()
      pool.join()
      print '\nProcesses have joined\n'

  print results
  if only_code is not None: return

  data = {}
  for d in sorted(os.listdir(os.getcwd())):
    if len(d)!=1: continue
    print 'D',d,os.getcwd()
    os.chdir(d)
    for filename in os.listdir(os.getcwd()):
      if filename.find("mogul")==-1: continue
      if filename.find(".pickle")==-1: continue
      print filename
      key = "vacuum"
      if filename.find("igb")>-1:
        key = "solvent"
      f=file(filename, "rb")
      mo = pickle.load(f)
      f.close()
      rc = mo.squared_deviations_and_numbers()
      print rc
      
      data.setdefault(key, {})
      for frag in rc:
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
      
    print data

    os.chdir("..")

  outl = {}
  for enviro in data:
    for frag in data[enviro]:
      outl.setdefault(frag, "")
      outl[frag] += "@type xy\n"
      for pair in data[enviro][frag]:
        outl[frag] += " %s %s\n" % (pair[1], pair[0])
      outl[frag] += "&\n"
      print outl[frag]

  for frag in outl:
    print '-'*80
    print outl[frag]
    print '-'*80
    gf = "%s.dat" % (frag)
    f=file(gf, "wb")
    f.write(outl[frag])
    f.close()

    cmd = "xmgrace -param plot.par -nxy %s" % gf
    print cmd
    os.system(cmd)


if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
