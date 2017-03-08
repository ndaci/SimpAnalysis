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
    trackCollection    = cms.InputTag("generalTracks", "", "RECO"),
    isData = cms.untracked.bool(True),
    phoLooseIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-loose"),
    phoMediumIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-medium"),
    phoTightIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-tight"), 
    pfRho = cms.InputTag("fixedGridRhoFastjetAll")
)
