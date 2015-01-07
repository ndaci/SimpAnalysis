import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIMJETS")

# Max events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5)
)

# Source
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/user/ndaci/Data/JetHT/JetHT_2012D_AOD_22Jan2013-v1/02BA259B-2B91-E211-9EDC-E0CB4E55363D.root'
        )
)

# Standard setup
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = cms.string(autoCond['startup'])

# MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.MessageLogger.cerr.threshold = 'INFO'

# Our Skim
from SimpAnalysis.SkimJets.skimJets_cfi import *
process.mySkimJets = skimJets.clone()

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('skimjets.root'),
    outputCommands = cms.untracked.vstring(
        #'keep *_*_*_*'
        'drop *',
        "keep *_TriggerResults_*_*",
        "keep *_hltTriggerSummaryAOD_*_*",
        'keep *_hltL1GtObjectMap_*_*',
        #'keep FEDRawDataCollection_rawDataCollector_*_*',
        #'keep FEDRawDataCollection_source_*_*',
        'keep edmTriggerResults_*_*_*',
        'keep triggerTriggerEvent_*_*_*',
        #"keep *_*_*_SIM",
        "keep *_genParticles_*_*",
        "keep *_l1extraParticles_*_*",
        "keep *_*ak5*_*_*",
        "keep *_*ak4*_*_*",
        "keep *_*Met*_*_*",
        "keep *_*MET*_*_*",
        #"keep *_*muon*_*_*",
        #"keep *_*Muon*_*_*",
        #"keep *_*electron*_*_*",
        #"keep *_*Electron*_*_*",
        "keep *_*track*_*_*",
        "keep *_*Track*_*_*",
        #"keep *_*photon*_*_*",
        #"keep *_*Photon*_*_*",
        #"keep *_*Tau*_*_*",
        #"keep *_*tau*_*_*"
        ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p')
        )
    )

# Paths
process.p = cms.Path(
    process.mySkimJets
    )

process.o = cms.EndPath(process.out)