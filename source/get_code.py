import os, sys

def run(d,i):
  filename = os.path.join(d, "amber.%s.out" % i)
  if os.path.exists(filename):
    f=file(filename, "rb")
    lines = f.readlines()
    f.close()
    for line in lines:
      line = line[:-1]
      #print line
      if line.find("chemical_components")>-1:
        return line.split("/")[-1].replace("data_", "").replace(".cif", "").lower()
  else:
    assert 0, 'no %s' % filename


if __name__=="__main__":
  print run(*tuple(sys.argv[1:]))
