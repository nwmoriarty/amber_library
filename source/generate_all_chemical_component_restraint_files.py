#######################################################################
#                                                                     #
# This is a copy of the file found in                                 #
# $PHENIX/modules/elbow/elbow/command_line/                           #
# as of August 1, 2018. It is placed here for reference to have all   #
# but is now the real version.                                        #
#                                                                     #
#######################################################################

import os, sys
import time
import pickle
import time
import random

import libtbx.load_env
from libtbx import phil, Auto
import libtbx.phil.command_line
from libtbx.option_parser import OptionParser
from mmtbx.chemical_components import generate_chemical_components_codes
from mmtbx.chemical_components import get_smiles, get_type
from elbow.utilities.chemical_components_utils import \
  get_elbow_molecule_from_chemical_components

from elbow.command_line import builder
from elbow.chemistry.Chemistry import elements
from elbow.utilities import geostd_utils

#from amber_adaptbx.command_line import AmberPrepClass as AmberPrep

from libtbx import easy_run
import StringIO

skip_codes = [
  "H", # single hydrogen
  "UNL",
  "UNK",
  "DUM",
  "UNX",
  "1P1", # huge cyclic
  "KEG", # large W cluster causes CCTBX error
  "WO2", # large W cluster causes CCTBX error
  "GLX", # ambiguous [C@H](CC[C@H]([N])C=O)([F,Cl,Br,I])[F,Cl,Br,I]
  "ASX", # "[C@@H](C[C@H]([N])C=O)([F,Cl,Br,I])[F,Cl,Br,I]"
  ]

master_params = """
restraints
{
  input
  {
    only_i = None
      .type = int
    only_code = None
      .type = str
    only_start = None
      .type = str
    qm_package = None
      .type = str
    qm_method = None
      .type = str
    qm_basis = None
      .type = str
    qm_solvent_model = None
      .type = str
  }
  control
  {
    dry_run = False
      .type = bool
    amber = False
      .type = bool
    use_temp_dir = True
      .type = bool
    ignore_md5 = True
      .type = bool
    ignore_output_files = True
      .type = bool
    ignore_queue_jobs = False
      .type = bool
    min_length_smiles = 0
      .type = int
    max_length_smiles = 100
      .type = int
    ligand_list = None
      .type = path
    skip_ligands_in_library = True
      .type = bool
    chunks_n = None
      .type = int
    chunks_i = None
      .type = int
    run_if_pH_same = Auto
      .type = bool
    only_type = None
      .type = str
    only_non_polymers = False
      .type = bool
    only_external_program = None
      .type = path
    exclude_external_program = None
      .type = path
  }
  properties {
    chromophore = False
      .type = bool
  }
  elbow
  {
    validate = False
      .type = bool
    no_output = False
      .type = bool
    pH = neutral
      .type = str
    opt_tol = default
      .type = str
  }
  output
  {
    method_basis_dir = None
      .type = path
    list_skipped_ligands = False
      .type = bool
  }
}
"""
master_phil = phil.parse(master_params,
                         process_includes=True,
                         )

metals = []
for i, e in enumerate(elements):
  if i<21: continue
  if 31<i<=38: continue
  if 49<i<=56: continue
  if 81<i<=88: continue
  metals.append("[%s" % e)

geostd_codes = []
for code in geostd_utils.generate_geostd_codes():
  geostd_codes.append(code)

top_dir = os.path.dirname(libtbx.env.dist_path("elbow"))
mon_lib_path = libtbx.env.find_in_repositories("chem_data")
mon_lib_path = os.path.join(mon_lib_path, "mon_lib")
mon_lib_codes = []
for d in os.listdir(mon_lib_path):
  if not os.path.isdir(os.path.join(mon_lib_path, d)): continue
  if len(d)>1: continue
  for filename in os.listdir(os.path.join(mon_lib_path, d)):
    if filename.find("_")!=-1: continue
    mon_lib_codes.append(filename.split(".")[0])

def setup_parser():
  usage="""
  elbow.generate_all_chemical_component_restraint_files only_code=NAG

  """
  parser = OptionParser(
    prog="elbow.generate_all_chemical_component_restraint_files",
    version="""
  up-to-date version
  """,
    usage=usage,
    )
  # Input options
  parser.add_option("",
                    "--show_defaults",
                    dest="show_defaults",
                    default=False,
                    action="store_true",
                    help="Display defaults",
                    )
  parser.add_option("",
                    "--dry-run",
                    dest="dry_run",
                    default=False,
                    action="store_true",
                    help="Display which residues will be processed",
                   )
  if 0:
    parser.add_option("",
                      "--verbose",
                      dest="verbose",
                      default=False,
                      action="store_true",
                      help="Verbose output",
                      )
    parser.add_option("",
                      "--silent",
                      dest="silent",
                      default=False,
                      action="store_true",
                      help="No output to screen",
                      )
  return parser

def setup_options_args(rargs):
  rargs = list(rargs)
  parser = setup_parser()
  (options, args) = parser.parse_args(args=rargs)
  if options.show_defaults:
    tmp = phil.parse(master_params,
                     process_includes=True,
                     )
    tmp.show()
    sys.exit()
  if len(args)==0:
    parser.print_help()
    sys.exit()

  from phenix.utilities import citations
  citations.add_citation('elbow')
  #
  argument_interpreter = libtbx.phil.command_line.argument_interpreter(
    master_phil=master_phil,
    home_scope="restraints")
  #
  phils = []
  phil_args = []
  for arg in args:
    if os.path.isfile(arg) :
      try :
        file_phil = phil.parse(file_name=arg)
      except RuntimeError :
        pass
      else :
        phils.append(file_phil)
    else :
      phil_args.append(arg)
      phils.append(argument_interpreter.process(arg))
  working_phil = master_phil.fetch(sources=phils)
  working_phil.show()
  working_params = working_phil.extract()
  working_params.restraints.control.dry_run = options.dry_run
  #check_working_params(working_params)
  #preamble = get_output_preamble(working_params)
  preamble = "generate_all"
  #preamble = preamble.replace(".pdb","")
  #in_scope = working_params.guided_ligand_replacement.input
  #for attr in in_scope.__dict__:
  #  tmp = getattr(in_scope, attr)
  #  if tmp is None: continue
  #  if type(tmp)!=type(""): continue
  #  if os.path.exists(tmp):
  #    setattr(in_scope, attr, os.path.abspath(tmp))
  print "  Writing effective parameters to %s.eff\n" % preamble
  #working_phil.format(python_object=working_params).show()
  print "#phil __ON__"
  master_phil.fetch_diff(source=master_phil.format(
      python_object=working_params)).show()
  if working_params.restraints.control.ligand_list:
    working_params.restraints.control.ligand_list = os.path.abspath(working_params.restraints.control.ligand_list)
  print "#phil __OFF__\n\n"
  f=file("%s.eff" % preamble, "wb")
  f.write(working_phil.format(python_object=working_params).as_str())
  f.close()
  return working_params

def load_pickle_molecule(d, code):
  pf = os.path.join(d, code.lower()[0], "%s.pickle" % code.upper())
  print pf, os.path.exists(pf)
  if os.path.exists(pf):
    f=file(pf, "rb")
    mol = pickle.load(f)
    f.close()
  else:
    mol = None
  return mol

def last_qm_lines(kwds, n=100):
  if qm_job_not_run(kwds): return ''
  if kwds.get("gamess", None):
    filename = "%s.gamess.gam" % kwds["chemical_component"]
    f=file(filename, "rb")
    lines = f.read()
    f.close()
    outl = ""
    for line in lines.splitlines()[n*-1:]:
      outl += "%s\n" % line
    return outl
  else:
    assert 0

def qm_job_not_run(kwds):
  outl = 'qm_job_not_run ? '
  if kwds.get("gamess", None):
    filename = "%s.gamess.gam" % kwds["chemical_component"]
    if not os.path.exists(filename):
      outl += " filename %s does not exist" % filename
      print outl
      return True
  else:
    return True
  outl += ' QM job file exists'
  print outl
  return False

def qm_job_run_and_finished(kwds):
  outl = 'qm_job_run_and_finished ? '
  if qm_job_not_run(kwds):
    outl += ' QM job file does not exist'
    print outl
    return False
  if kwds.get("gamess", None):
    filename = "%s.gamess.gam" % kwds["chemical_component"]
    f=file(filename, "rb")
    lines = f.read()
    f.close()
    for e in ["ddikick.x: application process 0 quit unexpectedly.",
              "semget errno=ENOSPC -- check system limit for sysv",
              ]:
      if lines.find(e)>-1:
        outl += ' found error : %s' % e
        print outl
        return False
    if lines.find("----- accounting info -----")!=-1:
      outl += ' found end line'
      print outl
      return True
  else:
    assert 0
  outl += ' ???'
  print outl
  return False

def delete_dat_files(lines):
  reading = False
  for line in lines.splitlines():
    if reading:
      tmp = line.split(",")
      if os.path.exists(tmp[0].strip()):
        print 'deleting',tmp[0]
        os.remove(tmp[0].strip())
    if line.find("Please save, rename, or erase these")>-1: reading=True

def qm_job_run_and_failed(kwds):
  outl = 'qm_job_run_and_failed ?'
  if qm_job_not_run(kwds):
    outl += ' QM job not run'
    print outl
    return False
  if kwds.get("gamess", None):
    filename = "%s.gamess.gam" % kwds["chemical_component"]
    f=file(filename, "rb")
    lines = f.read()
    f.close()
    if lines.find("THE GEOMETRY SEARCH IS NOT CONVERGED!")>-1:
      outl += ' QM job not converged'
      print outl
      return True
    if lines.find("......END OF GEOMETRY SEARCH......")>-1:
      outl += ' QM job end of geometry search'
      print outl
      return False
    if lines.find("Please save, rename, or erase these files from a previous run")>-1:
      outl += ' need to delete files, deleting file'
      print outl
      delete_dat_files(lines)
      return True
  else:
    assert 0
  return True

def is_same_molecule_regardless_of_pH(code, verbose=True):
  if verbose:
    print 'is_same_molecule_regardless_of_pH'
    print os.getcwd(), code
  rc = []
  default = os.path.join("..", "elbow.default")
  if not os.path.exists(default): return None
  default = load_pickle_molecule(default, code)
  try:
    if verbose:
      print default.DisplayBrief()
  except: return None

  dirs = ["elbow.low",
          "elbow.8",
          ]
  for other_directory in dirs:
    low = os.path.join("..", other_directory)
    if verbose:
      print 'default',default
      print 'low',low
    if not os.path.exists(low):
      rc.append(None)
      continue
    low = load_pickle_molecule(low, code)
    try:
      if verbose:
        print low.DisplayBrief()
    except:
      rc.append(None)
      continue
    print rc
    if verbose:  print 'len(default), len(low)',len(default), len(low)
    if len(default)==len(low):
      rc.append(True)
      #return True
    else:
      rc.append(False)
  print rc
  print filter(None, rc)
  if None in rc: return None
  if len(filter(None, rc))!=len(dirs):
    return False
  return True

def get_cc_md5(code):
  import hashlib
  cc_filename = os.path.join(top_dir, #os.environ["PHENIX"],
                             "chem_data",
                             "chemical_components",
                             "%s" % code[0].lower(),
                             "data_%s.cif" % code.upper(),
                             )
  print cc_filename
  f=file(cc_filename, "rb")
  lines=f.read()
  f.close()
  #print lines
  m=hashlib.md5()
  m.update(lines)
  print "MD5",m.hexdigest()
  return m.hexdigest()

def parse_ligand_list(filename):
  f=file(filename, "rb")
  lines = f.readlines()
  f.close()
  tmp = []
  for line in lines:
    tmp.append(line.split()[1])
  return tmp

def get_directory_name_from_method_basis(qm_method, qm_basis, qm_solvent_model):
  lookup = {
    "(" : "_lp_",
    ")" : "_rp_",
    "*" : "_star_",
    "," : "_comma_",
    }
  if qm_method is None and qm_basis is None:
    s = "elbow.default"
  else:
    s = "%s.%s" % (qm_method, qm_basis)
    if qm_solvent_model:
      s += ".%s" % qm_solvent_model
    for l in lookup:
      s = s.replace(l, lookup[l])
  return s

def display_kwds(kwds):
  print '\n  Running eLBOW with the following options'
  for key in sorted(kwds):
    print "    %s : %s" % (key, kwds[key])
  print

def get_elbow_molecule_cif_filename(code, kwds={}, return_molecule=False):
  def _get_qm_package_kwds(kwds, qm_package):
    if qm_package is not None:
      if qm_package.lower()=="gamess":
        kwds["gamess"]=True
      elif qm_package.lower()=="gaussian":
        kwds["gaussian"]=True
      #else:
      #  assert 0
    del kwds["qm_package"]
    return kwds

  if code in skip_codes:
    print 'code %s in skip list' % code
    return None

  code=code.upper()
  filename = "%s" % code
  default_kwds = {
    "qm_package" : None,
    "qm_method"  : None,
    "qm_basis"   : None,
    "qm_solvent_model" : None,
    "opt_nprocs" : 1,
    "opt"        : False,
    "mogul"      : False,
    "geostd"     : False,
    "qsub"       : "",
    "automatic"  : False,
    "qm_check"   : False,
    "dry_run"    : False,
    "quiet"      : True,
    #"validate"   : True,
    "pH"         : "neutral",
    #"xml"        : True,
    "write_redundant_dihedrals_boolean" : False, # special case!!!
    #"initial_geometry" : code,
    #"smiles" : get_smiles(code),
    #"file" : sdf_filename,
    #"no_output" : True,
    #"output": "ANP_neutral",
    "overwrite" : True,
    }
  kwds["chemical_component"] = code
  for dk in default_kwds:
    if dk not in kwds:
      kwds[dk]=default_kwds[dk]

  count = 0
  for kwd in [
    "opt",
    #"qm_method",
    "mogul",
    "automatic",
    ]:
    if kwds[kwd]:
      count+=1
  assert count<=1

  if 0:
    for key in sorted(kwds):
      print " %-20s : %s" % (key, kwds[key])

  kwds = _get_qm_package_kwds(kwds, kwds["qm_package"])

  if "smiles" in kwds and not kwds["smiles"]: return None

  if not kwds["overwrite"] and os.path.exists("%s.cif" % filename):
    print "CIF found",os.path.abspath("%s.cif" % filename)
    if return_molecule: assert 0
    return filename
  if not kwds["overwrite"] and os.path.exists("%s.pickle" % filename):
    print "pickle found",os.path.abspath("%s.pickle" % filename)
    if return_molecule: assert 0
    f=file(os.path.abspath("%s.pickle" % filename), "rb")
    mol = pickle.load(f)
    del f
    mol.WriteCIF(filename)
    return filename
  else:
    display_kwds(kwds)
    if qm_job_not_run(kwds):
      print '\n\n\tNo previous job\n\n'
# WHY WAS THIS INDENTED
#    if kwds["qm_method"]:
#      if qm_job_not_run(kwds):
#        print '\n\n\tNo previous job\n\n'
#        mol = builder.run(**kwds)
#      elif qm_job_run_and_failed(kwds):
#        print last_qm_lines(kwds, 100)
#        print '\n\n\tPrevious job is failed\n\n'
#        mol = builder.run(**kwds)
#      elif not qm_job_run_and_finished(kwds):
#        print '\n\n\tJob is still running\n\n'
#        mol = None
#      elif qm_job_run_and_finished(kwds):
#        print '\n\n\tJob is finished running\n\n'
#        assert 0
#      else:
#        print 'else'
#        mol = None
#    else:
      mol = builder.run(**kwds)
    elif not qm_job_run_and_finished(kwds):
      print '\n\n\tJob is still running\n\n'
      mol = None
    elif qm_job_run_and_failed(kwds):
      print last_qm_lines(kwds, 100)
      print '\n\n\tPrevious job is failed\n\n'
      mol = builder.run(**kwds)
    else:
      print 'else'
      mol = None
    if mol:
      if not os.path.exists(os.path.dirname(filename)):
        try: os.mkdir(os.path.dirname(filename))
        except: pass
      mol.WriteCIF(filename)

  if return_molecule:
    return mol
  else:
    return filename

def get_elbow_molecule_cif_filename_from_directory_tree(code,
                                                        params,
                                                        ):
  pr = params.restraints
  code=code.upper()
  cwd = os.getcwd()
  filename = os.path.join(os.getcwd(),
                          code[0].lower(),
                          "data_%s" % code,
                          )
  try: os.mkdir(code[0].lower())
  except: pass
  os.chdir(code[0].lower())
  kwds = {}
  kwds['qm_package'] = pr.input.qm_package
  kwds['qm_method'] = pr.input.qm_method
  kwds['qm_basis'] = pr.input.qm_basis
  kwds['qm_solvent_model'] = pr.input.qm_solvent_model
  kwds['opt_tol'] = pr.elbow.opt_tol
  #assert kwds['qm_solvent_model']
  #kwds['mogul'] = mogul
  kwds["overwrite"] = pr.control.ignore_output_files
  for attr in pr.elbow.__dict__:
    if attr.find("__")==0: continue
    kwds[attr] = getattr(pr.elbow, attr)
  filename = get_elbow_molecule_cif_filename(code, kwds)
  os.chdir(cwd)
  return filename

def calculate_amber_files(code,
                          pH=8,
                          use_temp_file=True,
                          ):
  import tempfile
  cwd = os.getcwd()
  if use_temp_file:
    td = tempfile.mkdtemp()
    os.chdir(td)
  mol = builder.run(chemical_component=code,
                    final_geometry=code,
                    no_output=True,
                    quiet=True,
                    pH=pH,
                    )
  # a very dull problem
  if hasattr(mol, "restraint_class"): del mol.restraint_class
  mol.OptimiseHydrogens()
  mol.WritePDB('4antechamber_%s.pdb' % code,
                pymol_pdb_bond_order=False)
  mol.WriteTriposMol2('4antechamber_%s.mol2' % code)
  mol.Multiplicitise()
  try: os.remove('sqm.pdb')
  except OSError: pass

  cmd='antechamber -i 4antechamber_%s.mol2 -fi mol2 -o %s.mol2 -fo mol2 \
    -nc %d -s 2 -pf y -c bcc -at gaff2 \
    -ek "qm_theory='AM1', scfconv=1.d-10, ndiis_attempts=700,"maxcyc=0"' \
    %(code, code, mol.charge)
  print cmd
  ero=easy_run.fully_buffered(cmd)
  f=open('antelog_%s.txt' % code,'wb')
  ero.show_stdout(out=f)
  ero.show_stderr(out=f)

  cmd='parmchk2 -i %s.mol2 -f mol2 -o %s.frcmod -s 2' %(code,code)
  print cmd
  ero=easy_run.fully_buffered(cmd)
  ero.show_stdout(out=f)
  ero.show_stderr(out=f)
  f.close()

  print '\nMoving files from',os.getcwd()
  cmd = "mv %s %s " % ("antelog_%s.txt" % (code),
                       os.path.join(cwd,"..","output", "antelog_%s.txt" % (code))
            )
  print cmd
  easy_run.call(cmd)
  cmd = "mv %s %s " % ("sqm.out",
                       os.path.join(cwd,"..","output", "sqm_%s.out" % (code))
            )
  print cmd
  easy_run.call(cmd)

  # copy the geometry used to create the inputs
  cmd = "mv %s %s " % ("4antechamber_%s.pdb" % (code),
                       os.path.join(cwd, "%s.final.pdb" % (code))
            )
  print cmd
  easy_run.call(cmd)

  for ext in ["frcmod", "mol2"]:
    if not os.path.exists("%s.%s" % (code, ext)): return None
    cmd = "mv %s %s" % ("%s.%s" % (code, ext),
                       os.path.join(cwd, "%s.%s" % (code,ext))
                       )
    print cmd
    easy_run.call(cmd)

  os.chdir(cwd)

  import shutil
  if use_temp_file:
    print "Deleting temporary directory : %s" % td
    shutil.rmtree(td)
  return True

def get_amber_filenames_from_directory_tree(code,
                                            ignore_output_files=True,
                                            pH="8",
                                            use_temp_file=True,
                                            ):
  code=code.upper()
  cwd = os.getcwd()
  try: os.mkdir(code[0].lower())
  except: pass
  os.chdir(code[0].lower())
  run_amber = True
  if not ignore_output_files:
    for ext in ["frcmod", "mol2", "final.pdb"]:
      print "%s.%s" % (code, ext),
      print os.path.exists("%s.%s" % (code, ext))
      if not os.path.exists("%s.%s" % (code, ext)):
        break
    else:
      run_amber = False
  rc = True
  if run_amber:
    rc = calculate_amber_files(code, pH=pH, use_temp_file=use_temp_file)
  else:
    print '\n\tAll output files found'
  os.chdir(cwd)
  return rc

def ligand_is_polymer(code):
  print 'code',code
  lt = get_type(code)
  print lt
  assert 0

def generate_ligand_codes(ligand_list_filename=None):
  cwd = os.getcwd()
  #if True:
  #  ligand_list = os.path.join(cwd, "all_chemical_components.dat")
  if ligand_list_filename:
    ligand_list = ligand_list_filename
  else:
    ligand_list = os.path.join(cwd, "ligands.dat")
  if os.path.exists(ligand_list):
    ligand_codes = parse_ligand_list(ligand_list)
  else:
    ligand_codes = []
  for code in generate_chemical_components_codes():
    if code not in ligand_codes:
      ligand_codes.append(code)
  for code in ligand_codes:
    #if code in ligand_codes: continue
    yield code

def run_ligand_loop(params):
  pr = params.restraints
  only_i    = pr.input.only_i
  only_code = pr.input.only_code
  only_start= pr.input.only_start
  chunks_n  = pr.control.chunks_n
  chunks_i  = pr.control.chunks_i
  list_skip = pr.output.list_skipped_ligands
  amber     = pr.control.amber
  min_length_smiles = pr.control.min_length_smiles
  max_length_smiles = pr.control.max_length_smiles
  only_type = pr.control.only_type
  #
  chromophore = pr.properties.chromophore
  #
  only_external_program = pr.control.only_external_program
  exclude_external_program = pr.control.exclude_external_program
  #
  try: only_i = int(only_i)
  except: pass
  try: chunks_n = int(chunks_n)
  except: pass
  try: chunks_i = int(chunks_i)
  except: pass
  if only_i==-1: only_i=None
  try:
    if only_code.lower()=="none":
      only_code=None
  except:
    pass
  #
  print 'hostname',os.environ.get("HOSTNAME", None)
  #
  if list_skip:
    assert only_i is None, "Can only list skipped ligands if only_i is None"
    assert only_code is None, "Can only list skipped ligands if only_code is None"
    all_codes = set()
    run_codes = set()
  #
  cwd = os.getcwd()
  if pr.output.method_basis_dir:
    method_basis_dir = pr.output.method_basis_dir
  elif amber:
    method_basis_dir = "amber_library"
  else:
    method_basis_dir = get_directory_name_from_method_basis(
      pr.input.qm_method,
      pr.input.qm_basis,
      pr.input.qm_solvent_model,
      )
  print '\n  Directory for output',method_basis_dir
  try:
    os.mkdir(method_basis_dir)
  except: pass
  os.chdir(method_basis_dir)
  try: os.mkdir("output")
  except: pass
  t0=time.time()
  non_polymers=0

  stats = {}

  for i, ligand_code in enumerate(generate_ligand_codes(
    ligand_list_filename=pr.control.ligand_list
    )):
    stats.setdefault("all", set())
    stats["all"].add(ligand_code)
    if list_skip:
      all_codes.add(ligand_code)
      assert 0
    if only_i is not None and only_i!=i+1: continue
    if chunks_i is not None and chunks_n is not None:
      if chunks_i!=i%chunks_n+1: continue
    if only_code is not None and ligand_code.upper()!=only_code.upper(): continue
    if only_code is not None: print i+1, ligand_code, only_i, chunks_i, chunks_n
    if only_start is not None:
      if ligand_code.find(only_start)!=0: continue

    print '  %5d %3s %s' % (i, ligand_code, get_type(ligand_code))
    if only_type is not None and only_type.lower()!=get_type(ligand_code).lower():
      stats.setdefault("polymer", set())
      stats["polymer"].add(ligand_code)
      continue

    # never run
    if ligand_code in skip_codes:
      print "Ligand skipped for any of a number of reasons"
      stats.setdefault("skip", set())
      stats["skip"].add(ligand_code)
      continue

    # don't reproduce the geostd and monomer lib
    if not amber and pr.control.skip_ligands_in_library:
      if ligand_code in geostd_codes:
        print '\n\tCode %s in GeoStd' % ligand_code
        stats.setdefault("geostd", set())
        stats["geostd"].add(ligand_code)
        continue
      if ligand_code in mon_lib_codes:
        print '\n\tCode %s in Monomer Library' % ligand_code
        stats.setdefault("monlib", set())
        stats["monlib"].add(ligand_code)
        continue

    if only_external_program:
      lines = os.popen('iotbx.python %s %s' % (only_external_program, ligand_code)).read()
      lines = lines.splitlines()
      assert len(lines)==1
      if lines[0] not in ['True']:
        print 'external program skipping',ligand_code
        continue

    if chunks_n is None and chunks_i is not None and str(chunks_i)!=ligand_code:
      continue

    if pr.control.only_non_polymers:
      print pr.control.only_non_polymers
      print type(pr.control.only_non_polymers)
      if ligand_is_polymer(ligand_code):
        stats.setdefault("polymer", set())
        stats["polymer"].add(ligand_code)
        continue
      else:
        assert 0

    smiles = get_smiles(ligand_code)
    if not smiles:
      print 'no SMILES found'
      stats.setdefault("no smiles", set())
      stats["no smiles"].add(ligand_code)
      continue
    if len(smiles)>max_length_smiles:
      print 'SMILES too long',max_length_smiles,len(smiles),smiles
      stats.setdefault("smiles too long", set())
      stats["smiles too long"].add(ligand_code)
      continue
    if len(smiles)<=min_length_smiles:
      print 'SMILES too short',min_length_smiles,len(smiles),smiles
      stats.setdefault("smiles too short", set())
      stats["smiles too short"].add(ligand_code)
      continue
    for e in metals:
      if smiles.find(e)>-1:
        print 'SMILES has metal',e,smiles
        stats.setdefault("smiles has metal", set())
        stats["smiles has metal"].add(ligand_code)
        break
    else:
      print '%5d %3s "%s"' % (i+1, ligand_code, smiles)
      #
      if pr.control.run_if_pH_same!=Auto:
        same = is_same_molecule_regardless_of_pH(ligand_code,
                                                 pr.control.run_if_pH_same,
                                                 )
        if same is not None:
          if pr.control.run_if_pH_same:
            if not same:
              print '\n\tpH has an effect'
              continue
          else:
            if same:
              print '\n\tpH has no effect'
              continue
        else:
          assert 0
      #
      if list_skip:
        run_codes.add(ligand_code)
        continue
      #
      if not os.path.exists("md5s"): os.mkdir("md5s")
      md5_filename = os.path.join("md5s", "%s.md5" % ligand_code)
      cc_md5 = get_cc_md5(ligand_code)
      if not pr.control.ignore_md5:
        if os.path.exists(md5_filename):
          f=file(md5_filename, "rb")
          old_md5 = pickle.load(f)
          f.close()
        else:
          old_md5 = None
        if old_md5==cc_md5:
          print "\n\tChemical Component file unchanged\n"
          continue
      #
      stats.setdefault("run", set())
      stats["run"].add(ligand_code)
      stats.setdefault("final", set())
      ff = os.path.join(
        ligand_code[0].lower(),
        "%s.final.pdb" % ligand_code,
        )
      if os.path.exists(ff):
        stats["final"].add(ligand_code)
      if pr.control.dry_run:
        print "\n\tRunning %s\n" % ligand_code
        if i>99 and 0:
          print '\n\tLeaving loop'
          break
        continue
      #
      if chromophore: continue
      #
      if amber:
        rc = get_amber_filenames_from_directory_tree(
          ligand_code,
          ignore_output_files=pr.control.ignore_output_files,
          pH=pr.elbow.pH,
          )
        if rc is None:
          print 'Calculation of Amber files failed'
          assert 0
      else:
        pickle_filename = get_elbow_molecule_cif_filename_from_directory_tree(
          ligand_code,
          params,
          )

      # should only do this on success
      f=file(md5_filename, "wb")
      pickle.dump(cc_md5, f)
      f.close()

    if only_i is not None: break
    if only_code is not None: break
  os.chdir(cwd)
  #
  double_single = {}
  for ligand_type, ligand_set in stats.items():
    if ligand_type in ['all']: continue
    for i, ligand_code in enumerate(ligand_set):
      if chromophore:
        double_single.setdefault(ligand_code, [])
        mol = get_elbow_molecule_from_chemical_components(ligand_code)
        for bond1 in mol.bonds:
          if bond1.order!=2: continue
          for bond2 in mol.bonds:
            if bond2.order!=2: continue
            if bond1==bond2: continue
            for bond3 in mol.bonds:
              if bond3.order!=1: continue
              if((bond3[0] in bond1 or bond3[0] in bond2) and
                 (bond3[1] in bond1 or bond3[1] in bond2)):
                if bond3 not in double_single[ligand_code]:
                  double_single[ligand_code].append(bond3)
  if double_single:
    for ligand_code, singles in sorted(double_single.items()):
      if len(singles)==0: continue
      print 'CHROMOPHORE',ligand_code, len(singles)
  #
  if list_skip:
    print '\n\nSkipped ligands'
    for i, ligand_code in enumerate(sorted(all_codes.difference(run_codes))):
      outl = ""
      smiles = get_smiles(ligand_code)
      for e in metals:
        if smiles.find(e)>-1:
          #print 'SMILES has metal',e,smiles
          outl += "%s " % e.replace("[", "")
      print " %-3d %-3s %-10s %-3d %s" % (i+1, ligand_code, outl, len(smiles), smiles)


  print 'Statistics'
  total = 0
  for act in stats:
    #print dir(stats[act])
    print "  %-20s : %5d %s" % (act, len(stats[act]), list(stats[act])[:10])
    if act not in ["all", "final"]: total+=len(stats[act])
  print '  Total : %d' % total
  try: missing = list(stats["run"].difference(stats["final"]))
  except: missing = []
  print '        : %s %d' % (missing, len(missing))


def run(rargs):
  print rargs
  working_params = setup_options_args(rargs)
  run_ligand_loop(working_params)

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(args)
