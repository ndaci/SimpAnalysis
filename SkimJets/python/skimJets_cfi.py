import FWCore.ParameterSet.Config as cms

skimJets = cms.EDFilter(
    "SkimJets",
    jetCollection = cms.InputTag("ak5PFJets"),
)

