from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'tnp_study_data_SingleElectron_2'
config.General.workArea = 'crab_projects'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runElectrons_TNPStudy.py'

config.section_("Data")
# config.Data.inputDataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
config.Data.inputDataset = '/SingleElectron/Run2015D-16Dec2015-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.publication = False
config.Data.ignoreLocality = True

config.section_("Site")
# Limit to US Tier2 sites. It appears to give more reliable performance.
config.Site.whitelist = ['T2_US_MIT','T2_US_UCSD','T2_US_Florida','T2_US_Wisconsin','T2_US_Caltech','T2_US_Purdue']
config.Site.storageSite = 'T2_US_Nebraska'
