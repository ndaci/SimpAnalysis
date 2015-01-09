import FWCore.ParameterSet.Config as cms

treep = cms.EDAnalyzer(
    'TreeProducer',
    hltProcessName   = cms.string("HLT"),
    pfjetCollection    = cms.InputTag("ak5PFJets")
)
