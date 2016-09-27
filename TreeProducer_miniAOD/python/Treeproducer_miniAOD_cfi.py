import FWCore.ParameterSet.Config as cms

tree = cms.EDAnalyzer(
    'TreeProducer_miniAOD',
    triggerResults   = cms.InputTag("TriggerResults", "", "HLT"),
    pfjetCollection  = cms.InputTag("slimmedJets"),
    vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices")
)
