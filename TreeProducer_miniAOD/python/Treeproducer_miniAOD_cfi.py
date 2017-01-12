import FWCore.ParameterSet.Config as cms

tree = cms.EDAnalyzer(
    'TreeProducer_miniAOD',
    triggerResults   = cms.InputTag("TriggerResults", "", "HLT"),
    #triggerResults2  = cms.InputTag("TriggerResults", "", "HLT2"),
    prescales        = cms.InputTag("patTrigger", "", "RECO"),
    prescales2        = cms.InputTag("patTrigger", "", "PAT"),
    METfilter        = cms.InputTag("TriggerResults", "", "RECO"),
    pfjetCollection  = cms.InputTag("slimmedJets"),
    genjetCollection  = cms.InputTag("slimmedGenJets"),
    vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
    METCollection    = cms.InputTag("slimmedMETs"),
    photonCollection    = cms.InputTag("slimmedPhotons"),
    packedPFCollection    = cms.InputTag("packedPFCandidates"),
    isData = cms.untracked.bool(True)
)
