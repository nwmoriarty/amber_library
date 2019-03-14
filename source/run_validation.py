import os, sys
import pickle

def move_to_errors(code):
  print 'moving files to errors directory'
  print code
  err = os.path.join('..', 'errors_directory')
  if not os.path.exists(err):
    os.mkdir(err)
  for filename in os.listdir(os.getcwd()):
    if filename.startswith(code):
      os.rename(filename, os.path.join(err, filename))

def run(only_code):
  d={}
  if only_code.lower()=="none": return
  print only_code
  os.chdir(only_code[0])
  for ext in ["mol2",
              #"prmtop",
              "frcmod",
              ]:
    filename = os.path.join(os.getcwd(), "%s.%s" % (only_code, ext))
    if os.path.exists(filename):
      d[ext]=filename
    #else:
    #  d[ext]=None
  #if len(d)!=2: move_to_errors(only_code.upper())
  for filename in os.listdir(os.getcwd()):
    if not filename.startswith(only_code): continue
    print filename
    if filename.endswith('.pickle'):
      f=file(filename, 'rb')
      mo = pickle.load(f)
      f.close()
      print mo
      print dir(mo)

if __name__=="__main__":
  run(sys.argv[1])
