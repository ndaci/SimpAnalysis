import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:myfile.root'
    )
)

process.load("SimpAnalysis.TreeProducer.treeproducer_cfi") 
process.p = cms.Path(process.treep)
