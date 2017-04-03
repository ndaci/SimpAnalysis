import FWCore.ParameterSet.Config as cms

process = cms.Process("SIMPTREE")


process.load('Configuration.StandardSequences.GeometryRecoDB_cff')

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v8', '')
process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
 duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),

    fileNames = cms.untracked.vstring(
#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M1/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M1_1.root",
#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M1/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M1_2.root",
#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M1/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M1_3.root",
#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M1/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M1_4.root",
#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M1/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M1_6.root",
#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M1/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M1_7.root",
#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M1/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M1_8.root",
#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M1/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M1_9.root",
#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M1/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M1_11.root"

#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M10/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M10_1.root",
#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M10/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M10_2.root",
#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M10/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M10_3.root",
#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M10/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M10_4.root",
#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M10/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M10_5.root",
#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M10/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M10_6.root",
#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M10/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M10_7.root",
#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M10/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M10_8.root",
"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M10/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M10_9.root",
#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M10/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M10_10.root",
"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M10/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M10_11.root",
"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M10/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M10_12.root"

#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M100/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M100_2.root",
#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M100/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M100_3.root",
#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M100/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M100_4.root",
#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M100/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M100_5.root",
#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M100/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M100_6.root",
#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M100/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M100_8.root"


#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M200/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M200.root"

#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M400/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M400.root"

#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M700/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M700.root"

#"file:///pnfs/iihe/cms/store/user/gflouris/SIMPS/SignalProduction/SIMP_M1000/AOD_SIMPs/SUS-RunIISummer16DR80Premix-00068_SIMP_M1000.root"
    )
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
    #connect = cms.string('sqlite:JECs/Summer16_23Sep2016V4_MC.db'),
    #toGet = cms.VPSet(
        #cms.PSet(
            #record = cms.string('JetCorrectionsRecord'),
            #tag    = cms.string('JetCorrectorParametersCollection_Summer16_23Sep2016V4_MC_AK4PFchs'),
            #label  = cms.untracked.string('AK4PFchs')
        #),
    #)
#)

## Add an ESPrefer to override JEC that might be available from the global tag
#process.es_prefer_jec = cms.ESPrefer('PoolDBESSource', 'jec')

#process.load("JetMETCorrections.Configuration.CorrectedJetProducers_cff")
#process.load("JetMETCorrections.Configuration.JetCorrectors_cff")
#from JetMETCorrections.Configuration.JetCorrectors_cff import ak4PFCHSL1FastL2L3CorrectorChain
#process.ak4PFCHSJetsCorr = process.ak4PFCHSJetsL1.clone()
#process.ak4PFCHSJetsCorr.correctors = cms.VInputTag('ak4PFCHSL1FastL2L3Corrector')
#process.ak4PFCHSJetsCorr.src = cms.InputTag('ak4PFCHSJetsSIMPs')

## Tree producer
#process.load("SimpAnalysis.TreeProducer_AOD.Treeproducer_AOD_cfi")
#process.tree.pfjetCollection  = cms.InputTag("ak4PFCHSJetsCorr")
##process.tree.pfjetCollection  = cms.InputTag("ak4PFCHSJetsSIMPs") #uncorrected


#process.p = cms.Path(process.egmPhotonIDSequence *process.ak4PFCHSL1FastL2L3CorrectorChain*process.ak4PFCHSJetsCorr*process.tree)

# Remove first PV
# Create goodOfflinePrimaryVertrices
# Run CHS
# Correct jets
# Tree producer
######################
from CommonTools.ParticleFlow.pfPileUp_cfi  import pfPileUp as _pfPileUp
from CommonTools.ParticleFlow.TopProjectors.pfNoPileUp_cfi import pfNoPileUp as _pfNoPileUp
from CommonTools.ParticleFlow.goodOfflinePrimaryVertices_cfi import *


process.pfPileUpSPV = _pfPileUp.clone(PFCandidates='SIMPsPFPtrs',
                                      Vertices = 'SPVgoodOfflinePrimaryVertices',
                                      checkClosestZVertex = False )
process.pfNoPileUpSPV = _pfNoPileUp.clone(topCollection = 'pfPileUpSPV',
                                          bottomCollection = 'SIMPsPFPtrs' )

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

process.treeSPV = process.tree.clone()
process.treeSPV.pfjetCollection = cms.InputTag("ak4PFCHSJetsSPV")
process.treeSPV.vertexCollection = cms.InputTag("SPVgoodOfflinePrimaryVertices")
process.treeSPV.isData = cms.untracked.bool(False)
#############################
process.treeCorr = process.tree.clone()
process.treeCorr.pfjetCollection = cms.InputTag("ak4PFCHSJetsSIMPs")
process.treeCorr.pfRho = cms.InputTag("fixedGridRhoFastjetAll")
process.treeCorr.isData = cms.untracked.bool(False)

process.p = cms.Path(process.egmPhotonIDSequence                          
                    +process.treeCorr
                    +process.pfchsSecondPV
                    +process.SPVgoodOfflinePrimaryVertices
                    +process.pfPileUpSPV
                    +process.pfNoPileUpSPV
                    +process.ak4PFCHSJetsSPV
                    +process.treeSPV)

# Output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('SIMPs_PUMoriond17_AOD_M10_3.root')
)
