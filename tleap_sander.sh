#!/bin/bash

i=$1

cat<<EOF > leap.in
source leaprc.gaff
x = loadMol2 $i.mol2
loadAmberParams $i.frcmod
saveAmberParm x $i.prmtop $i.rst7
quit
EOF

tleap -f leap.in > tleap.out 2>&1
err=$?
if [ $err -gt 0 ]; then
	echo "error in tleap for $i"
	exit 1
fi
/bin/rm leap.in leap.log tleap.out

cat<<EOF > sander.in
  short minimization
 &cntrl
   imin=1, maxcyc=100, ncyc=50, ntpr=1, cut=99., ntb=0,
 /
EOF

sander -O -i sander.in -p $i.prmtop -c $i.rst7 -o $i.min.out -r $i.min.rst7
err=$?
/bin/rm sander.in mdinfo
if [ $err -gt 0 ]; then
	echo "error in sander for $i"
	exit 1
fi

ambpdb -p $i.prmtop <$i.min.rst7 >$i.min.pdb 2>/dev/null