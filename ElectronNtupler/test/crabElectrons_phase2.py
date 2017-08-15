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

config.Site.storageSite = 'T2_US_Wisconsin'

pds = [
    '/TTTo2L2Nu_TuneCUETP8M1_14TeV-powheg-pythia8/PhaseIITDRSpring17MiniAOD-PU200_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/DYJetsToLL_M-10to50_TuneCUETP8M1_14TeV-madgraphMLM-pythia8/PhaseIITDRSpring17MiniAOD-PU200_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
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
        config.General.requestName = 'phase2electrons_%d_%s' % (i, primaryDS)
        config.Data.outputDatasetTag = conditions
        config.Data.inputDataset = pd
        #submit(config)
