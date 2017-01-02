import FWCore.ParameterSet.Config as cms

tree = cms.EDAnalyzer(
    'TreeProducer_AOD',
    triggerResults   = cms.InputTag("TriggerResults", "", "HLT"),
    #triggerResults2  = cms.InputTag("TriggerResults", "", "HLT2"),
    prescales        = cms.InputTag("patTrigger", "", "RECO"), #not found
    prescales2        = cms.InputTag("patTrigger", "", "PAT"), #not found
    METfilter        = cms.InputTag("TriggerResults", "", "RECO"),
    pfjetCollection  = cms.InputTag("ak4PFJetsCHS"),
    genjetCollection  = cms.InputTag("ak4GenJetsNoNu", "", "HLT"),
    vertexCollection = cms.InputTag("offlinePrimaryVertices"),
    METCollection    = cms.InputTag("pfMet"),
    photonCollection    = cms.InputTag("gedPhotons"),
    trackCollection    = cms.InputTag("generalTracks", "", "RECO")    
)
