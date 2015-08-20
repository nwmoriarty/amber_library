#######################################################################
#                                                                     #
# This is the bash script used to run AmberLibrary on the cci machine #
#                                                                     #
#  author: Nigel Moriarty                                             #
#                                                                     #
#######################################################################

setup gcc-4.3.4_fc8_x86_64
source /net/chevy/raid1/nigel/phenix_amber_build/setpaths.csh

elbow.generate_all_chemical_component_restraint_files \
  amber=True \
  ignore_output_files=False \
  skip_ligands_in_library=False \
  only_type=NON-POLYMER \
  pH=8 \
  only_i=$SGE_TASK_ID >& output/amber.$SGE_TASK_ID.output

setenv code `elbow.python amber_library/source/get_code.py output $SGE_TASK_ID`
setenv opwd $PWD
cd amber_library
source /net/chevy/raid1/nigel/bootstrap/external/amber14/amber.csh

elbow.python $opwd/amber_library/source/run_tleap_sander.py $code >& ../output/tleap.$SGE_TASK_ID.output
elbow.python $opwd/amber_library/source/run_mogul.py $code >& ../output/mogul.$SGE_TASK_ID.output
elbow.python $opwd/amber_library/source/run_validation.py $code

