from WMCore.Configuration import Configuration
config = Configuration()

name = 'Simps_JetHT_2012D_treep_v3_13Jan2015'

# GENERAL
config.section_("General")
config.General.requestName = name 
config.General.workArea    = '/user/ndaci/CRABBY/Simps/Data/JetHT/'
config.General.transferLogs = True
#config.General.transferOutput = True

# JOB TYPE
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'treeproducer_cfg.py'
#config.JobType.pyCfgParams = []
#config.JobType.inputFiles = ''
config.JobType.allowNonProductionCMSSW = True

# INPUT DATA
config.section_("Data")
config.Data.inputDataset = '/JetHT/ndaci-Simps_JetHT_2012D_v4_13Jan2015-c6764b838ac38d8cf7e5bf224b39d63a/USER'
config.Data.inputDBS  = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSReader/'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 250
config.Data.publication = False
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
config.Data.publishDataName = name
config.Data.ignoreLocality = False # allows to process inputs on CE != site hosting inputs
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions12/8TeV/Reprocessing/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt'
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
