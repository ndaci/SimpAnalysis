import FWCore.ParameterSet.Config as cms

tree = cms.EDAnalyzer(
    'TreeProducer_miniAOD',
    hltProcessName   = cms.string("HLT"),
    pfjetCollection  = cms.InputTag("ak5PFJets"),
    vertexCollection = cms.InputTag("offlinePrimaryVertices")
)
