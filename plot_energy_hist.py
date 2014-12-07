#! /usr/bin/python
from numpy import *
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import sys



files = ['ene_final_min',
         'ene_final_min_igb',
         'ene_initial_min',
         'ene_initial_min_igb']
# files = ['ene_final_min']


plt.rc(('xtick.major','xtick.minor','ytick.major','ytick.minor'), pad=10)
plt.rc('axes',linewidth=3)
plt.rc('legend', fontsize=20) 
plt.rc('lines', markeredgewidth=2)
plt.rc('xtick.minor',size=5)
plt.rc('xtick.major',size=10)
plt.rc('lines', linewidth=3) 


for file in files:
  data = genfromtxt('%s.dat' %file)
  ene = data[:,1]
  rmsd = data[:,2]
  if 'final' in file:
    mask = (rmsd<1000) & (ene>-1000) & (ene<1000)
    print "%s.dat outliers (gradient RMSD or |energy|>1000): %d" %(file, len(rmsd)-sum(mask))
    ene = ene[mask]
  else:
    mask = (ene>-1000) & (ene<1000)
    print "%s.dat outliers (|energy|>1000): %d" %(file, len(rmsd)-sum(mask))
    ene = ene[mask]
  # import code; code.interact(local=dict(globals(), **locals()))
  fig=plt.figure(figsize=(16, 12))
  ax = fig.add_subplot(111)
  plt.hist(ene,100)
  plt.savefig('%s.pdf' %file)
	


