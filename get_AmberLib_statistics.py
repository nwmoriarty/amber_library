#! /usr/bin/python
import os, sys
import glob
import subprocess

dirs = '0 1 2 3 4 5 6 7 8 9 a b c d e f g h i j k l m n o p q r s t u v w x y z'.split()
dirs = ['0']

def check_log_for_errors(filename, igb=False):
  if not os.path.exists(filename):
    min_pdb_err = "does not exist %s" % filename
  else:
    if os.stat(filename).st_size ==0:
      print filename
      assert 0
    else:
      with open(filename) as f:
        i=-1
        min_pdb_err = f.readlines()[i].strip()
  return min_pdb_err

def run(only_code=None):
  subprocess.call('rm tmp*', shell=True)
  missing_mol2=open('missing_mol2.dat','w')
  missing_pdb = open('missing_min_pdb.dat', 'w')

  n_pdbs = 0
  n_mol2s = 0
  n_frcmods = 0
  n_min_pdbs = 0
  n_min_igb_pdbs = 0
  for d in dirs:
    pdbs = glob.glob('%s/*final.pdb' %d)
    mol2s = glob.glob('%s/*mol2' %d)
    frcmods = glob.glob('%s/*frcmod' %d)
    min_pdbs = glob.glob('%s/*min.pdb' %d)
    min_igb_pdbs = glob.glob('%s/*min_igb.pdb' %d)

    n_pdbs += len(pdbs)
    n_mol2s += len(mol2s)
    n_frcmods += len(frcmods)
    n_min_pdbs += len(min_pdbs)
    n_min_igb_pdbs += len(min_igb_pdbs)

    for pdb in pdbs:
      code=pdb.split('.')[0]
      mol2=code+'.mol2'
      if mol2 not in mol2s:
        missing_mol2.write('%s\n' %code)

    for mol2 in mol2s:
      code = mol2.split('.')[0]
      is_missing = False
      min_pdb_err = ''
      min_igb_pdb_err = ''

      min_pdb = code+'.min.pdb'
      if min_pdb not in min_pdbs:
        is_missing = True
        min_pdb_err = check_log_for_errors("%s.min.out" % code)

      min_igb_pdb = code+'.min_igb.pdb'
      if min_igb_pdb not in min_igb_pdbs:
        is_missing = True
        min_igb_pdb_err = check_log_for_errors("%s.min_igb.out" % code)

        #with open('%s.min_igb.out' %code) as f:
        #  i=-1
        #  lines=f.readlines()
        #  while min_igb_pdb_err == '':
        #    min_igb_pdb_err = lines[i].strip()
        #    i-=1

      if is_missing:
        outl = '%s | %-80s | %-80s\n' %(code,
                                        min_pdb_err,
                                        min_igb_pdb_err)
        print outl
        missing_pdb.write(outl)

    print 'getting "%s" min energies' %d
    cmd = '''
      grep -H -A5 \"FINAL RESULTS\" %s/*.min.out | grep -v \"FINAL RESULTS\" | grep -v \"NSTEP\" | grep -v \'^.*.min.out-$\'|grep -v \'^..$\' |\
       awk \'{printf(\"%%s      %%7.2f     %%5.2f\\n\",  $1,$3,$4)}\' | \
       sort -k1 -u | sort -k2 -n >> tmp1
       ''' %(d)
    x=subprocess.call(cmd, shell=True)
    print 'getting "%s" min_igb energies' %d
    cmd = '''
      grep -H -A5 \"FINAL RESULTS\" %s/*.min_igb.out | grep -v \"FINAL RESULTS\" | grep -v \"NSTEP\" | grep -v \'^.*.min_igb.out-$\'|grep -v \'^..$\' |\
       awk \'{printf(\"%%s      %%7.2f     %%5.2f\\n\",  $1,$3,$4)}\' | \
       sort -k1 -u | sort -k2 -n >> tmp2
       ''' %(d)
    x=subprocess.call(cmd, shell=True)
    print 'getting "%s" initial min energies' %d
    cmd = '''
      grep -H \"^      1\" %s/*.min.out | \
       awk \'{printf(\"%%s      %%7.2f     %%5.2f\\n\",  $1,$3,$4)}\' | \
       sort -k1 -u | sort -k2 -n >> tmp3
       ''' %(d)
    x=subprocess.call(cmd, shell=True)
    print 'getting "%s" initial min_igb energies' %d
    cmd = '''
      grep -H \"^      1\" %s/*.min_igb.out | \
       awk \'{printf(\"%%s      %%7.2f     %%5.2f\\n\",  $1,$3,$4)}\' | \
       sort -k1 -u | sort -k2 -n >> tmp4
       ''' %(d)
    x=subprocess.call(cmd, shell=True)

  missing_mol2.close()
  missing_pdb.close()

  subprocess.call('sort -k2 -n tmp1 > ene_final_min.dat', shell=True)
  subprocess.call('sort -k2 -n tmp2 > ene_final_min_igb.dat', shell=True)
  subprocess.call('sort -k2 -n tmp3 > ene_initial_min.dat', shell=True)
  subprocess.call('sort -k2 -n tmp4 > ene_initial_min_igb.dat', shell=True)
  subprocess.call('rm tmp?', shell=True)

  print 'Total no of *.final.pdb files (made by eLBOW): %d' %n_pdbs
  print 'Total no of *.mol2 files (made by antechamber): %d' %n_mol2s
  print 'Total no of *.mol2 files (made by prmchk2): %d' %n_frcmods
  print 'Total no of *.min.pdb files (made by sander): %d' %n_min_pdbs
  print 'Total no of *.min_igb.pdb files (made by sander): %d' %n_min_igb_pdbs

  from numpy import genfromtxt
  import itertools
  outfile = open('outliers_min_energy.dat','w')
  with open('ene_final_min.dat','r') as f:
    lines=f.readlines()
  data = genfromtxt('ene_final_min.dat')
  ene = data[:,1]
  rmsd = data[:,2]
  mask = (rmsd>1000)|(ene<-1000)|(ene>1000)
  lines = list(itertools.compress(lines,mask))
  outfile.write(''.join(lines))

  outfile = open('outliers_min_igb_energy.dat','w')
  with open('ene_final_min_igb.dat','r') as f:
    lines=f.readlines()
  data = genfromtxt('ene_final_min_igb.dat')
  ene = data[:,1]
  rmsd = data[:,2]
  mask = (rmsd>1000)|(ene<-1000)|(ene>1000)
  lines = list(itertools.compress(lines,mask))
  outfile.write(''.join(lines))

  mol2s=[]
  for d in dirs:
    mol2s_t = glob.glob('%s/*mol2' %d)
    mol2s = mol2s+mol2s_t
  mol2s=[mol2.split('.')[0] for mol2 in mol2s]
  with open('missing_min_pdb.dat') as f:
    lines=f.readlines()
  missingpdb = [line.split()[0] for line in lines]
  with open('ene_final_min.dat') as f:
    lines=f.readlines()
  enepdb = [line.split()[0].split('.')[0] for line in lines]
  #import code; code.interact(local=dict(globals(), **locals()))
  print '_'*80
  for mol2 in mol2s:
    if mol2 not in missingpdb and mol2 not in enepdb:
      print mol2

if __name__=="__main__":
  run(*tuple(sys.argv[1:]))
