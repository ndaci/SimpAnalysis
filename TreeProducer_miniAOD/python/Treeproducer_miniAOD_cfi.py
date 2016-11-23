import FWCore.ParameterSet.Config as cms

tree = cms.EDAnalyzer(
    'TreeProducer_miniAOD',
    triggerResults   = cms.InputTag("TriggerResults", "", "HLT"),
    prescales        = cms.InputTag("patTrigger", "", "RECO"),
    METfilter        = cms.InputTag("TriggerResults", "", "RECO"),
    pfjetCollection  = cms.InputTag("slimmedJets"),
    vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
    METCollection    = cms.InputTag("slimmedMETs")
)
