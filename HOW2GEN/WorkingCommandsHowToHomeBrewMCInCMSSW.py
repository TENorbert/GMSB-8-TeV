### Simple How to to Produce  MC Samples Privately in using CMSSW
### Assuming you do not have your own Private Configuration files
### If you have yours which you have produced from say ISAJET or other Model
## Implementation packages likemyself for SUSY GMSB which is here: https://github.com/TENorbert/GMSB-8-TeV/tree/master/mGMSB_CONFIG_FILES

## You might want to put those in the Configuration/Generator/data directory
### make the relevant configuration files (.py) files
## and put them in  Configuration/GenProduction/python/EightTeV

## Check whcih architecture has been configured

echo $SCRAM_ARCH
#### Choice of particular architecture:

## In bash shell
export SCRAM_ARCH=slc5_amd64_gcc462
## In tcsh/csh Shell
setenv SCRAM_ARCH slc5_amd64_gcc462

### Check available CMSSW Releases on a give acrhitecture say SCRAM_ARCH =slc5_amd64_gcc462
scram -arch slc5_amd64_gcc462 list CMSSW
### Get the Release & set configuration:
cmsrel CMSSW_5_3_17
cd CMSSW_5_3_17/src
cmsenv

## Get Relevant MC configurations and specific MC production interface eg Pythia6Interface
### other MC Generation interfaces are found in here: https://github.com/cms-sw/cmssw/tree/CMSSW_5_3_X/GeneratorInterface
git cms-addpkg Configuration/Generator
git cms-addpkg GeneratorInterface/Pythia6Interface

scram b -j9;

### Now Begin the MC Generation and RECO Process!!
##1)
##********************
##GEN-SIM:

	 cmsDriver.py Configuration/GenProduction/python/EightTeV/GMSB_Lambda180_CTau50_8TeV_pythia6_cff.py --fileout file:EXO-Summer12-02641.root --mc --eventcontent RAWSIM --datatier GEN-SIM --conditions START53_V7C::All --beamspot Realistic8TeVCollision --step GEN,SIM 

##2)
##********************
##HLT
    cmsDriver.py step1 --filein "dbs:/GMSB_Lambda-180_CTau-50_TuneZ2star_8TeV-pythia6/Summer12-START53_V7C_ext1-v1/GEN-SIM" --fileout file:EXO-Summer12DR53X-02697_step1.root --pileup_input "dbs:/MinBias_TuneZ2star_8TeV-pythia6/Summer12-START50_V13-v3/GEN-SIM" --mc --eventcontent RAWSIM --pileup 2012_Summer_50ns_PoissonOOTPU --datatier GEN-SIM-RAW --conditions START53_V19::All --step DIGI,L1,DIGI2RAW,HLT:7E33v2 

##3)
##********************
##HLT-AODSIM
   cmsDriver.py step2 --filein file:EXO-Summer12DR53X-02697_step1.root --fileout file:EXO-Summer12DR53X-02697.root --mc --eventcontent AODSIM,DQM --datatier AODSIM,DQM --conditions START53_V19::All --step RAW2DIGI,L1Reco,RECO,VALIDATION:validation_prod,DQM:DQMOfflinePOGMC 


##HLT-RECO (I think)
	 cmsDriver.py step2 --filein file:EXO-Summer12DR53X-02697_step1.root --fileout file:EXO-Summer12DR53X-02697.root --mc --eventcontent RECO,DQM --datatier RECO,DQM --conditions START53_V19::All --step RAW2DIGI,L1Reco,RECO,VALIDATION:validation_prod,DQM:DQMOfflinePOGMC 

###(or try replacing RECO with AODSIM, RECO to produce both)

###********************



Here’s where the I got the info from:
GENSIM: 
https://cms-pdmv.cern.ch/mcm/requests?prepid=EXO-Summer12-02641&page=0&shown=17179869311 

DIGI-RECO:
https://cms-pdmv.cern.ch/mcm/requests?prepid=EXO-Summer12DR53X-02697&page=0&shown=17179869311





### working Steps for GEN SIM

## Step 1 GEN-SIM
cmsDriver.py Configuration/GenProduction/python/EightTeV/GMSB_Lambda180_CTau50_8TeV_pythia6_cff.py --fileout file:EXO-Summer12-02641.root --mc --eventcontent RAWSIM --datatier GEN-SIM --conditions START53_V7C::All --beamspot Realistic8TeVCollision --step GEN,SIM -n 5 --no_exec

 cmsRun GMSB_Lambda180_CTau50_8TeV_pythia6_cff_py_GEN_SIM.py

## Step 2 DIGI to HLT

cmsDriver.py step1 --filein "file:./EXO-Summer12-02641.root" --fileout file:EXO-Summer12DR53X-02697_step1.root --mc --eventcontent RAWSIM --pileup 2012_Summer_50ns_PoissonOOTPU --datatier GEN-SIM-RAW --conditions START53_V19::All --step DIGI,L1,DIGI2RAW,HLT:7E33v2 -n 5 --no_exec

cmsRun step1_DIGI_L1_DIGI2RAW_HLT_PU.py

## step 3 HLT to RECO

cmsDriver.py step2 --filein file:EXO-Summer12DR53X-02697_step1.root --fileout file:EXO-Summer12DR53X-02697.root --mc --eventcontent RECOSIM,DQM --datatier RECOSIM,DQM --conditions START53_V19::All -s RAW2DIGI,L1Reco,RECO,DQM:DQMOfflinePOGMC -n 5 --no_exec

cmsRun step2_RAW2DIGI_L1Reco_RECO_DQM.py


### And Thats it!! Now you can go ahead and  Analyse your events(RECO) using
### you custom designed and Build Analyzer

