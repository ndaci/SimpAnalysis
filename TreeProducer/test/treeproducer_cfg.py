import FWCore.ParameterSet.Config as cms

process = cms.Process("SIMPTREE")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/user/ndaci/JetHT/Simps_JetHT_2012D_08Jan2015_v1/150108_171434/0000/skimjets_122.root'
    )
)

# Tree producer
process.load("SimpAnalysis.TreeProducer.treeproducer_cfi") 
process.p = cms.Path(process.treep)

# Output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('file:simptree.root')
)
