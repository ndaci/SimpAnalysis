import FWCore.ParameterSet.Config as cms

process = cms.Process("SIMPTREE")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
			#'/store/mc/RunIISpring16MiniAODv2/QCD_Pt-15to3000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/FlatPU8to37HcalNZSRAW_withHLT_80X_mcRun2_asymptotic_v14-v1/00000/00E73272-0C34-E611-9DD4-0025905B85CA.root')
								#'/store/data/Run2016C/JetHT/MINIAOD/23Sep2016-v1/50000/7A50D439-BA89-E611-9BA5-0025905B8586.root')
								'/store/mc/RunIISummer16MiniAODv2/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/110000/0055B499-54B6-E611-9F86-FA163E1F94C5.root')
                 #'/store/mc/RunIISpring16MiniAODv2/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v2/00000/00223E2E-F827-E611-8974-0CC47A4D764A.root')
                 #'/store/mc/RunIISpring16MiniAODv2/XXTo4J_M-3000_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/00000/3A221244-A439-E611-A405-0025905B85DA.root', 
                 #'/store/mc/RunIISpring16MiniAODv2/XXTo4J_M-3000_CTau-1000mm_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/40000/D66F8F6C-053A-E611-A189-0025905A6134.root')
)

#process.triggerSelection = cms.EDFilter( "TriggerResultsFilter",
    #triggerConditions = cms.vstring(
      #'HLT_DiCentralPFJet170_CFMax0p1_v*',
      #'HLT_DiCentralPFJet220_CFMax0p3_v*',
      #'HLT_DiCentralPFJet330_CFMax0p5_v*',
      #'HLT_DiCentralPFJet430_v*',
      #'HLT_DiCentralPFJet170_v*',
      #'HLT_SingleCentralPFJet170_CFMax0p1_v*'),
    #hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
    #l1tResults = cms.InputTag( "" ),
    #l1tIgnoreMask = cms.bool( False ),
    #l1techIgnorePrescales = cms.bool( False ),
    #daqPartitions = cms.uint32( 1 ),
    #throw = cms.bool( True )
#)

#import HLTrigger.HLTfilters.hltHighLevel_cfi
#process.triggerSelection =  HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    #HLTPaths = ['HLT_DiCentralPFJet170_CFMax0p1_v*',
      #'HLT_DiCentralPFJet220_CFMax0p3_v*',
      #'HLT_DiCentralPFJet330_CFMax0p5_v*',
      #'HLT_DiCentralPFJet430_v*',
      #'HLT_DiCentralPFJet170_v*',
      #'HLT_SingleCentralPFJet170_CFMax0p1_v*'],
    #throw = False
#)

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDPhotonIdProducer(process, dataFormat)
my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring16_V2p2_cff']
for idmod in my_id_modules:
	setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)
    
# Tree producer
process.load("SimpAnalysis.TreeProducer_miniAOD.Treeproducer_miniAOD_cfi") 
process.tree.METfilter = cms.InputTag("TriggerResults", "", "PAT")
#process.tree.triggerResults = cms.InputTag("TriggerResults", "", "HLT2") #for XXTo4J
process.p = cms.Path(process.egmPhotonIDSequence * process.tree)

# Output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('QCD_PUMoriond17_test.root')
)
