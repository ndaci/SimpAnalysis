import FWCore.ParameterSet.Config as cms

process = cms.Process("SIMPTREE")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring(
        #'/store/user/ndaci/JetHT/Simps_JetHT_2012D_v4_13Jan2015/150113_110917/0000/skimjets_52.root'
        ##'/store/user/ndaci/JetHT/Simps_JetHT_2012D_08Jan2015_v1/150108_171434/0000/skimjets_122.root'
        ##'file:skimjets.root'
    #)
    fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/data/Run2016D/JetHT/MINIAOD/PromptReco-v2/000/276/315/00000/283C90B1-0445-E611-B0AF-02163E014393.root')
)

# Tree producer
process.load("SimpAnalysis.TreeProducer_miniAOD.Treeproducer_miniAOD_cfi") 
process.p = cms.Path(process.tree)

# Output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('simptree.root')
)
