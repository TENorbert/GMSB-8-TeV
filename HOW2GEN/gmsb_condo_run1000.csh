#! /bin/tcsh
source /uscmst1/prod/sw/cms/setup/cshrc prod
setenv SCRAM_ARCH slc5_amd64_gcc472
cd /uscms_data/d3/tnorbert/MC_GEN/WORKING_RELEASE/CMSSW_6_1_2/src/
eval `scramv1 runtime -csh`
cd ${_CONDOR_SCRATCH_DIR}
cmsRun /uscms_data/d3/tnorbert/MC_GEN/WORKING_RELEASE/CMSSW_6_1_2/src/GMSB_Lambda180_CTau1000_8TeV_pythia6_cff_py_GEN_SIM.py
 

