from CRABClient.UserUtilities import config
config = config()
#from WMCore.Configuration import Configuration
#config = Configuration()

name = 'SIMPs_JetHT_2016'

# GENERAL
config.section_("General")
config.General.requestName = name 
config.General.workArea    = '/user/isdebruy/SIMPs/Data2016/'
config.General.transferLogs = True
#config.General.transferOutput = True

# JOB TYPE
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'treeproducer_cfg.py'
#config.JobType.pyCfgParams = []
#config.JobType.inputFiles = ''
#config.JobType.allowNonProductionCMSSW = True

# INPUT DATA
config.section_("Data")
config.Data.inputDataset = '/JetHT/Run2016B-PromptReco-v1/MINIAOD'
config.Data.inputDBS  = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 250
config.Data.publication = False
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
config.Data.outputDatasetTag = name
#config.Data.ignoreLocality = False # allows to process inputs on CE != site hosting inputs
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/Cert_271036-282037_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt'
#config.Data.runRange = 

#A custom string to insert in the output file name inside the CRAB-created directory path to allow organizing groups of tasks.
#config.Data.prefix =  

# USER
config.section_("User")
#config.User.email = 'nadir.daci@cern.ch'
#config.User.voRole = 
config.User.voGroup = 'becms'

# GRID
config.section_("Site")
config.Site.storageSite = 'T2_BE_IIHE'
#config.Site.whitelist = 
config.Site.blacklist = ['T1_US_FNAL']
