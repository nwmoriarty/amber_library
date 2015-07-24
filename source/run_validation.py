import os, sys
import pickle

def run(only_code):
  d={}
  if only_code.lower()=="none": return
  print only_code
  os.chdir(only_code[0])
  for ext in ["mol2",
              "prmtop",
              "frcmod",
              ]:
    filename = os.path.join(os.getcwd(), "%s.%s" % (only_code, ext))
    if os.path.exists(filename):
      d[ext]=filename
    else:
      d[ext]=None
  print d

if __name__=="__main__":
  run(sys.argv[1])
  
