import FWCore.ParameterSet.Config as cms

process = cms.Process("SIMPTREE")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016SeptRepro_v7', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
										'/store/data/Run2016G/JetHT/AOD/23Sep2016-v1/100000/0006CE1E-9986-E611-8DFB-6C3BE5B5C0B0.root')
)

#process.triggerSelection = cms.EDFilter( "TriggerResultsFilter",
    #triggerConditions = cms.vstring(
      #'HLT_DiCentralPFJet170_CFMax0p1_v*',
      #'HLT_DiCentralPFJet220_CFMax0p3_v*',
      #'HLT_DiCentralPFJet330_CFMax0p5_v*',
      #'HLT_DiCentralPFJet430_v*',
      #'HLT_DiCentralPFJet170_v*',
      #'HLT_SingleCentralPFJet170_CFMax0p1_v*'),
    #hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
    #l1tResults = cms.InputTag( "" ),
    #l1tIgnoreMask = cms.bool( False ),
    #l1techIgnorePrescales = cms.bool( False ),
    #daqPartitions = cms.uint32( 1 ),
    #throw = cms.bool( True )
#)

import HLTrigger.HLTfilters.hltHighLevel_cfi
process.triggerSelection =  HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    HLTPaths = ['HLT_DiCentralPFJet170_CFMax0p1_v*',
      'HLT_DiCentralPFJet220_CFMax0p3_v*',
      'HLT_DiCentralPFJet330_CFMax0p5_v*',
      'HLT_DiCentralPFJet430_v*',
      'HLT_DiCentralPFJet170_v*',
      'HLT_SingleCentralPFJet170_CFMax0p1_v*'],
    throw = False
)

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.AOD
switchOnVIDPhotonIdProducer(process, dataFormat)
my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring16_V2p2_cff']
for idmod in my_id_modules:
	setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)
	
#from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup
#process.jec = cms.ESSource('PoolDBESSource',
    #CondDBSetup,
    #connect = cms.string('sqlite:Summer16_23Sep2016V4_MC.db'),
    #toGet = cms.VPSet(
        #cms.PSet(
            #record = cms.string('JetCorrectionsRecord'),
            #tag    = cms.string('JetCorrectorParametersCollection_Summer16_23Sep2016V4_MC_AK4PFchs'),
            #label  = cms.untracked.string('AK4PFchs')
        #),
    #)
#)
	
# Add an ESPrefer to override JEC that might be available from the global tag
#process.es_prefer_jec = cms.ESPrefer('PoolDBESSource', 'jec')

####Correct ak4PFCHS Jets
#process.load("JetMETCorrections.Configuration.CorrectedJetProducers_cff")
#process.load("JetMETCorrections.Configuration.JetCorrectors_cff")
#from JetMETCorrections.Configuration.JetCorrectors_cff import ak4PFCHSL1FastL2L3CorrectorChain
##process.ak4PFJetsL1L2L3Residual.src = cms.InputTag("ak4PFJetsCHS")
#process.ak4PFCHSJetsCorr = process.ak4PFCHSJetsL1.clone()
#process.ak4PFCHSJetsCorr.correctors = cms.VInputTag('ak4PFCHSL1FastL2L3Corrector')  

# Remove first PV
# Create goodOfflinePrimaryVertrices
# Run CHS
# Correct jets
# Tree producer
######################
from CommonTools.ParticleFlow.pfPileUp_cfi  import pfPileUp as _pfPileUp
from CommonTools.ParticleFlow.TopProjectors.pfNoPileUp_cfi import pfNoPileUp as _pfNoPileUp
from CommonTools.ParticleFlow.goodOfflinePrimaryVertices_cfi import *


process.pfPileUpSPV = _pfPileUp.clone(PFCandidates='particleFlowPtrs',
                                      Vertices = 'SPVgoodOfflinePrimaryVertices',
                                      checkClosestZVertex = False )
process.pfNoPileUpSPV = _pfNoPileUp.clone(topCollection = 'pfPileUpSPV',
                                          bottomCollection = 'particleFlowPtrs' )

process.SPVgoodOfflinePrimaryVertices = goodOfflinePrimaryVertices.clone()
process.SPVgoodOfflinePrimaryVertices.src = src = cms.InputTag("pfchsSecondPV:SecondOfflinePrimaryVerices")

process.load('SimpAnalysis.PFCHS_SecondPV.pfchs_SecondPV_cfi')

process.load('RecoJets.Configuration.RecoPFJets_cff')
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets

process.ak4PFCHSJetsSPV = ak4PFJets.clone(
    rParam = 0.4,
    jetPtMin = 5.0,
    src = "pfNoPileUpSPV",
)

#process.ak4SPVPFCHSJetsCorr = process.ak4PFCHSJetsL1.clone()
#process.ak4SPVPFCHSJetsCorr.correctors = cms.VInputTag('ak4PFCHSL1FastL2L3Corrector')
#process.ak4SPVPFCHSJetsCorr.src = cms.InputTag('ak4PFCHSJetsSPV')
    
# Tree producer
process.load("SimpAnalysis.TreeProducer_AOD.Treeproducer_AOD_cfi") 
#process.tree.triggerResults = cms.InputTag("TriggerResults", "", "HLT2") #for XXTo4J
process.tree.triggerName = cms.vstring(
"HLT_DiCentralPFJet170_v1","HLT_DiCentralPFJet170_v2","HLT_DiCentralPFJet170_v3", "HLT_DiCentralPFJet170_v4", "HLT_DiCentralPFJet170_v5",
"HLT_SingleCentralPFJet170_CFMax0p1_v1", "HLT_SingleCentralPFJet170_CFMax0p1_v1", "HLT_SingleCentralPFJet170_CFMax0p1_v3", "HLT_SingleCentralPFJet170_CFMax0p1_v4", "HLT_SingleCentralPFJet170_CFMax0p1_v5")

process.treeSPV = process.tree.clone()
process.treeSPV.pfjetCollection = cms.InputTag("ak4PFCHSJetsSPV")
process.treeSPV.vertexCollection = cms.InputTag("SPVgoodOfflinePrimaryVertices")
process.treeSPV.isData = cms.untracked.bool(True)
#############################
process.treeCorr = process.tree.clone()
process.treeCorr.pfjetCollection = cms.InputTag("ak4PFJetsCHS")
process.treeCorr.pfRho = cms.InputTag("fixedGridRhoFastjetAll")
process.treeCorr.isData = cms.untracked.bool(True)

process.p = cms.Path(process.triggerSelection
                    +process.egmPhotonIDSequence                          
                    +process.treeCorr
                    +process.pfchsSecondPV
                    +process.SPVgoodOfflinePrimaryVertices
                    +process.pfPileUpSPV
                    +process.pfNoPileUpSPV
                    +process.ak4PFCHSJetsSPV
                    +process.treeSPV)

# Output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('Data_AOD_test.root')
)
