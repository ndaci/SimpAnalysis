from CRABClient.UserUtilities import config
config = config()

name = 'SIMPs_QCD_HT2000ToInf_ext_PUMoriond17_AOD'

# GENERAL
config.section_("General")
config.General.requestName = name 
config.General.workArea    = '/user/isdebruy/SIMPs/CMSSW_8_0_20/src/SimpAnalysis/TreeProducer_AOD/test/crab_ntuples/'
config.General.transferLogs = True
#config.General.transferOutput = True

# JOB TYPE
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../treeproducer_AOD_MC_cfg.py'
#config.JobType.pyCfgParams = []
config.JobType.inputFiles = ['../../../../../external/slc6_amd64_gcc530/data/RecoEgamma/PhotonIdentification/data/Spring16/photon_general_MVA_Spring16_EB_V3.weights.xml', '../../../../../external/slc6_amd64_gcc530/data/RecoEgamma/PhotonIdentification/data/Spring16/photon_general_MVA_Spring16_EE_V3.weights.xml']
#config.JobType.allowNonProductionCMSSW = True

# INPUT DATA
config.section_("Data")
config.Data.inputDataset = '/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/AODSIM'
config.Data.inputDBS  = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 4
config.Data.publication = False
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
config.Data.outputDatasetTag = name
config.Data.ignoreLocality = False # allows to process inputs on CE != site hosting inputs
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/Cert_271036-282037_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt'
#config.Data.runRange = 

#A custom string to insert in the output file name inside the CRAB-created directory path to allow organizing groups of tasks.
#config.Data.prefix =  

# USER
#config.section_("User")
#config.User.email = 'nadir.daci@cern.ch'
#config.User.voRole = 
#config.User.voGroup = 'becms'

# GRID
config.section_("Site")
config.Site.storageSite = 'T2_BE_IIHE'
config.Site.ignoreGlobalBlacklist = True
config.Site.whitelist = ['T2_ES_IFCA']
config.Site.blacklist = ['T3_US_Rice', 'T2_BR_UERJ', 'T2_TR_METU', 'T2_EE_Estonia', 'T2_PK_NCP', 'T2_RU_PNPI', 'T3_RU_FIAN', 'T3_US_MIT', 'T3_UK_London_UCL', 'T3_US_UCD', 'T3_CO_Uniandes', 'T3_TW_NTU_HEP', 'T3_ES_Oviedo', 'T3_US_NU', 'T3_US_NotreDame', 'T2_RU_SINP', 'T2_CH_CERN_AI', 'T3_IN_PUHEP', 'T2_CH_CSCS_HPC', 'T2_GR_Ioannina', 'T2_CH_CERN_HLT', 'T2_MY_UPM_BIRUNI', 'T3_UK_London_RHUL', 'T2_TH_CUNSTDA', 'T3_US_Kansas', 'T3_US_Princeton_ICSE', 'T3_IN_TIFRCloud', 'T0_CH_CERN', 'T3_GR_IASA', 'T3_CN_PKU', 'T3_US_Baylor', 'T2_PL_Warsaw', 'T2_RU_INR', 'T3_US_JHU', 'T3_BY_NCPHEP', 'T3_US_FSU', 'T3_KR_UOS', 'T3_CH_PSI']
