from CRABClient.UserUtilities import config
config = config()

name = 'SIMPs_XXTo4J_ctau2000_M1500'

# GENERAL
config.section_("General")
config.General.requestName = name 
config.General.workArea    = '/user/isdebruy/SIMPs/SignalMC/'
config.General.transferLogs = True
#config.General.transferOutput = True

# JOB TYPE
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../treeproducer_MC_cfg.py'
#config.JobType.pyCfgParams = []
config.JobType.inputFiles = ['../../../../../external/slc6_amd64_gcc530/data/RecoEgamma/PhotonIdentification/data/Spring16/photon_general_MVA_Spring16_EB_V3.weights.xml', '../../../../../external/slc6_amd64_gcc530/data/RecoEgamma/PhotonIdentification/data/Spring16/photon_general_MVA_Spring16_EE_V3.weights.xml']
#config.JobType.allowNonProductionCMSSW = True

# INPUT DATA
config.section_("Data")
config.Data.inputDataset = '/XXTo4J_M-1500_CTau-2000mm_TuneCUETP8M1_13TeV_pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM'
config.Data.inputDBS  = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
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
#config.Site.whitelist = 
config.Site.blacklist = ['T1_US_FNAL','T2_UA_KIPT','T2_UK_London_Brunel','T2_CH_CSCS','T2_US_*','T1_RU_JINR']
