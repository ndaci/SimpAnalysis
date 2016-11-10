import FWCore.ParameterSet.Config as cms

process = cms.Process("SIMPTREE")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
                 '/store/mc/RunIISpring16MiniAODv2/XXTo4J_M-3000_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/00000/3A221244-A439-E611-A405-0025905B85DA.root') 
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
    
# Tree producer
process.load("SimpAnalysis.TreeProducer_miniAOD.Treeproducer_miniAOD_cfi") 
process.p = cms.Path(process.tree)

# Output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('simptree.root')
)
