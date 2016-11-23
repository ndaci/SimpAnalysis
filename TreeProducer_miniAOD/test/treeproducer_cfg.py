import FWCore.ParameterSet.Config as cms

process = cms.Process("SIMPTREE")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
                 '/store/data/Run2016D/JetHT/MINIAOD/PromptReco-v2/000/276/315/00000/283C90B1-0445-E611-B0AF-02163E014393.root',
                 '/store/data/Run2016D/JetHT/MINIAOD/PromptReco-v2/000/276/315/00000/48DA42E9-0145-E611-B53F-02163E011A4F.root')
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

import HLTrigger.HLTfilters.hltHighLevel_cfi
process.triggerSelection =  HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    HLTPaths = ['HLT_DiCentralPFJet170_CFMax0p1_v*',
      'HLT_DiCentralPFJet220_CFMax0p3_v*',
      'HLT_DiCentralPFJet330_CFMax0p5_v*',
      'HLT_DiCentralPFJet430_v*',
      'HLT_DiCentralPFJet170_v*',
      'HLT_SingleCentralPFJet170_CFMax0p1_v*'],
    throw = False
)
    
# Tree producer
process.load("SimpAnalysis.TreeProducer_miniAOD.Treeproducer_miniAOD_cfi") 
process.p = cms.Path(process.triggerSelection * process.tree)

# Output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('simptree.root')
)
