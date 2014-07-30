import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

# run the input file through the end;
# for a limited number of events, replace -1 with the desired number 
#
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load( "SimGeneral.HepPDTESSource.pythiapdt_cfi" )

process.source = cms.Source( "PoolSource",
                             fileNames = cms.untracked.vstring(
			     #'file:../../GMSB_Lambda140_CTau6000_8TeV_pythia6_cff_py_GEN_SIM.root'
			     #'file:../../GMSB_Lambda100_CTau7000_8TeV_pythia6_cff_py_GEN.root'
			     #'file:../../GMSB_Lambda180_CTau6000_8TeV_pythia6_cff_py_GEN.root'
			     #'file:../../GMSB_Lambda180_CTau4000_8TeV_pythia6_cff_py_GEN.root'
			     #'file:../../GMSB_Lambda180_CTau3000_8TeV_pythia6_cff_py_GEN.root'
			     #'file:../../GMSB_Lambda160_CTau250_8TeV_pythia6_cff_py_GEN.root'
			     #'file:../../GMSB_Lambda140_CTau500_8TeV_pythia6_cff_py_GEN.root'
			     #'file:../../GMSB_Lambda120_CTau1000_8TeV_pythia6_cff_py_GEN.root'
			     #'file:../../GMSB_Lambda120_CTau2000_8TeV_pythia6_cff_py_GEN.root'
			    # 'file:../../GMSB_Lambda160_CTau3000_8TeV_pythia6_cff_py_GEN.root'
			    #'file:../../GMSB_Lambda160_CTau4000_8TeV_pythia6_cff_py_GEN.root'
			     'file:../../GMSB_Lambda160_CTau6000_8TeV_pythia6_cff_py_GEN.root'
			     #'file:../../GMSB_Lambda180_CTau2000_8TeV_pythia6_cff_py_GEN.root'
			     #'file:../../GMSB_Lambda180_CTau1000_8TeV_pythia6_cff_py_GEN.root'
			     #'file:../../GMSB_Lambda180_CTau500_8TeV_pythia6_cff_py_GEN.root'
			    # 'file:../../GMSB_Lambda180_CTau250_8TeV_pythia6_cff_py_GEN.root'
			     #'file:GMSB_Lambda180_CTau6000_8TeV_pythia6_cff_py_GEN_SIM.root'
			     #'file:GMSB_Lambda100_CTau6000_8TeV_Cgrav329_pythia6_cff_py_GEN_SIM.root'
			     #'file:GMSB_Lambda100_CTau20000_8TeV_pythia6_cff_py_GEN_SIM.root'
			     # 'file:../GGM_M3_1420_Msq_1400_M1_375_slha_8TeV_pythia6_cff_py_GEN.root'
			     #'file:/uscms_data/d3/tnorbert/MC_GEN/WORKING_RELEASE/CMSSW_6_1_2/src/GMSB_Lambda140_CTau6000_8TeV_pythia6_cff_py_GEN.root'
			     )
                           )
	      
# FileService is mandatory, as the following analyzer module 
# will want it, to create output histogram file
# 
process.TFileService = cms.Service("TFileService",
        fileName = cms.string("Outputfile.root")
)

# the analyzer itself - empty parameter set 
#
process.TestHepMCEvt = cms.EDAnalyzer("MyGenAnalyzer")

process.p1 = cms.Path( process.TestHepMCEvt )

