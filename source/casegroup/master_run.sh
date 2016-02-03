#!/bin/sh

#  Hai Nguyen adapted from Nigel Moriarty's script

# require
#  - $AMBERHOME
#  - $PHENIX

# python source code dir
export source=/project1/dacase-001/haichit/phenix_amber/run_amber_library/amber_library/source/
export opwd=$PWD
# seq="0QE 0GY"
seq=`cat $source/casegroup/ligand_codes.dat | head -100`

for code in $seq; do
    echo $code

    elbow.python $source/generate_all_chemical_component_restraint_files.py \
        amber=True \
        ignore_output_files=True \
        skip_ligands_in_library=False \
        only_type=NON-POLYMER \
        pH=8 \
        only_code=$code \
        only_i=None >& output/amber.$code.output
    
    # code=`elbow.python $source/get_code.py output $code`
    # cd amber_library
    # 
    # elbow.python $source/run_tleap_sander.py $code --force >& ../output/tleap.$code.output
    # elbow.python $source/run_mogul.py $code >& ../output/mogul.$code.output
    # elbow.python $source/run_validation.py $code
    # cd $opwd
done
