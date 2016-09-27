import FWCore.ParameterSet.Config as cms

tree = cms.EDAnalyzer(
    'TreeProducer_miniAOD',
    hltProcessName   = cms.string("HLT"),
    pfjetCollection  = cms.InputTag("slimmedJets"),
    vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices")
)
