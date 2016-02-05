#!/bin/sh

#  Hai Nguyen adapted from Nigel Moriarty's script

# require
#  - $AMBERHOME
#  - $PHENIX

# python source code dir
export source=/project1/dacase-001/haichit/phenix_amber/run_amber_library/amber_library/source/
export opwd=$PWD

# code is in 2nd colum
# seq=`cat $source/casegroup/ligand_codes.dat | awk '{print $2}'`

# one column
seq=`cat output/done_codes.dat`

for code in $seq; do
    echo $code

    code=`elbow.python $source/casegroup/get_code.py output $code`
    if [ $code ]; then
        cd amber_library
        elbow.python $source/casegroup/run_tleap_sander.py $code --force >& ../output/tleap.$code.output
        cd $opwd
     fi
done
