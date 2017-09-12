from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.workArea = 'crab_phase2electrons'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runElectrons_phase2.py'

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/%s/phase2electrons' % (getUsernameFromSiteDB())
config.Data.publication = False

#config.Site.storageSite = 'T2_US_Wisconsin'
config.Site.storageSite = 'T2_CH_CERN'

pds = [
    #'/DoubleElectron_FlatPt-1To300/PhaseIITDRSpring17MiniAOD-PU200CalAging300NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    #'/DoubleElectron_FlatPt-1To300/PhaseIITDRSpring17MiniAOD-PU200CalAging1000NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    #'/DoubleElectron_FlatPt-1To300/PhaseIITDRSpring17MiniAOD-PU200CalAging3000NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    #'/DoubleElectron_FlatPt-1To300/PhaseIITDRSpring17MiniAOD-PU200CalAging4500UltimateNewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/DoubleElectron_FlatPt-1To300/PhaseIITDRSpring17MiniAOD-noPUCaloAging300NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/DoubleElectron_FlatPt-1To300/PhaseIITDRSpring17MiniAOD-noPUCaloAging1000NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/DoubleElectron_FlatPt-1To300/PhaseIITDRSpring17MiniAOD-noPUCaloAging3000NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/DoubleElectron_FlatPt-1To300/PhaseIITDRSpring17MiniAOD-noPUCaloAging4500UltimateNewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
]

if __name__ == '__main__':
    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    for i, pd in enumerate(pds):
        (_, primaryDS, conditions, dataTier) = pd.split('/')
        config.General.requestName = 'p2electronsPU0_%d_%s' % (i, primaryDS)
        config.Data.outputDatasetTag = conditions
        config.Data.inputDataset = pd
        submit(config)
