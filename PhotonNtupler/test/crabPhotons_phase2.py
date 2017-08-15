from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.workArea = 'crab_phase2photons'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runPhotons_phase2.py'

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2
config.Data.useParent = True
config.Data.outLFNDirBase = '/store/user/%s/phase2photons' % (getUsernameFromSiteDB())
config.Data.publication = False

#config.Site.storageSite = 'T2_US_Wisconsin'
config.Site.storageSite = 'T2_CH_CERN'

pds_old = [
    '/GluGluHToGG_M125_14TeV_amcatnloFXFX_pythia8/PhaseIITDRSpring17MiniAOD-PU200_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/QCD_Flat_Pt-15to7000_TuneCUETP8M1_14TeV_pythia8/PhaseIITDRSpring17MiniAOD-PU200_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/DiPhotonJetsBox_MGG-80toInf_14TeV-Sherpa/PhaseIITDRSpring17MiniAOD-PU200_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
]

pds = [
    '/SinglePhoton_Pt-5/PhaseIITDRSpring17MiniAOD-noPUCaloAging4500Ultimate_91X_upgrade2023_realistic_v3-v4/MINIAODSIM',
    '/SinglePhoton_Pt-5/PhaseIITDRSpring17MiniAOD-noPUCaloAging300_91X_upgrade2023_realistic_v3-v4/MINIAODSIM',
    '/SinglePhoton_Pt-5/PhaseIITDRSpring17MiniAOD-noPUCaloAging3000_91X_upgrade2023_realistic_v3-v4/MINIAODSIM',
    '/SinglePhoton_Pt-5/PhaseIITDRSpring17MiniAOD-noPUCaloAging1000_91X_upgrade2023_realistic_v3-v4/MINIAODSIM',
    '/SinglePhoton_Pt-5/PhaseIITDRSpring17MiniAOD-PU200CalAging4500Ultimate_91X_upgrade2023_realistic_v3-v4/MINIAODSIM',
    '/SinglePhoton_Pt-5/PhaseIITDRSpring17MiniAOD-PU200CalAging300_91X_upgrade2023_realistic_v3-v4/MINIAODSIM',
    '/SinglePhoton_Pt-5/PhaseIITDRSpring17MiniAOD-PU200CalAging3000_91X_upgrade2023_realistic_v3-v4/MINIAODSIM',
    '/SinglePhoton_Pt-5/PhaseIITDRSpring17MiniAOD-PU200CalAging1000_91X_upgrade2023_realistic_v3-v4/MINIAODSIM',
    '/SinglePhoton_Pt-10/PhaseIITDRSpring17MiniAOD-PU200CalAging1000_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    '/SinglePhoton_Pt-10/PhaseIITDRSpring17MiniAOD-PU200CalAging3000_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-10/PhaseIITDRSpring17MiniAOD-PU200CalAging300_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-10/PhaseIITDRSpring17MiniAOD-PU200CalAging4500Ultimate_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-10/PhaseIITDRSpring17MiniAOD-noPUCaloAging1000_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    '/SinglePhoton_Pt-10/PhaseIITDRSpring17MiniAOD-noPUCaloAging3000_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-10/PhaseIITDRSpring17MiniAOD-noPUCaloAging300_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-20/PhaseIITDRSpring17MiniAOD-PU200CalAging1000_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    '/SinglePhoton_Pt-20/PhaseIITDRSpring17MiniAOD-PU200CalAging3000_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    '/SinglePhoton_Pt-20/PhaseIITDRSpring17MiniAOD-PU200CalAging300_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    '/SinglePhoton_Pt-20/PhaseIITDRSpring17MiniAOD-PU200CalAging4500Ultimate_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    '/SinglePhoton_Pt-20/PhaseIITDRSpring17MiniAOD-noPUCaloAging1000_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    '/SinglePhoton_Pt-20/PhaseIITDRSpring17MiniAOD-noPUCaloAging3000_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    '/SinglePhoton_Pt-20/PhaseIITDRSpring17MiniAOD-noPUCaloAging300_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    '/SinglePhoton_Pt-20/PhaseIITDRSpring17MiniAOD-noPUCaloAging4500Ultimate_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    '/SinglePhoton_Pt-50/PhaseIITDRSpring17MiniAOD-noPUCaloAging4500Ultimate_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-50/PhaseIITDRSpring17MiniAOD-noPUCaloAging300_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-50/PhaseIITDRSpring17MiniAOD-noPUCaloAging3000_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-50/PhaseIITDRSpring17MiniAOD-noPUCaloAging1000_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    '/SinglePhoton_Pt-50/PhaseIITDRSpring17MiniAOD-PU200CalAging4500Ultimate_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-50/PhaseIITDRSpring17MiniAOD-PU200CalAging300_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-50/PhaseIITDRSpring17MiniAOD-PU200CalAging3000_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-50/PhaseIITDRSpring17MiniAOD-PU200CalAging1000_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    '/SinglePhoton_Pt-100/PhaseIITDRSpring17MiniAOD-noPUCaloAging4500Ultimate_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-100/PhaseIITDRSpring17MiniAOD-noPUCaloAging300_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-100/PhaseIITDRSpring17MiniAOD-noPUCaloAging3000_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-100/PhaseIITDRSpring17MiniAOD-noPUCaloAging1000_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    '/SinglePhoton_Pt-100/PhaseIITDRSpring17MiniAOD-PU200CalAging4500Ultimate_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-100/PhaseIITDRSpring17MiniAOD-PU200CalAging300_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-100/PhaseIITDRSpring17MiniAOD-PU200CalAging3000_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-100/PhaseIITDRSpring17MiniAOD-PU200CalAging1000_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    '/SinglePhoton_Pt-200/PhaseIITDRSpring17MiniAOD-noPUCaloAging4500Ultimate_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-200/PhaseIITDRSpring17MiniAOD-noPUCaloAging300_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-200/PhaseIITDRSpring17MiniAOD-noPUCaloAging3000_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-200/PhaseIITDRSpring17MiniAOD-noPUCaloAging1000_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    '/SinglePhoton_Pt-200/PhaseIITDRSpring17MiniAOD-PU200CalAging4500Ultimate_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-200/PhaseIITDRSpring17MiniAOD-PU200CalAging300_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-200/PhaseIITDRSpring17MiniAOD-PU200CalAging3000_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-200/PhaseIITDRSpring17MiniAOD-PU200CalAging1000_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    '/SinglePhoton_Pt-500/PhaseIITDRSpring17MiniAOD-noPUCaloAging4500Ultimate_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-500/PhaseIITDRSpring17MiniAOD-noPUCaloAging300_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-500/PhaseIITDRSpring17MiniAOD-noPUCaloAging3000_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-500/PhaseIITDRSpring17MiniAOD-noPUCaloAging1000_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    '/SinglePhoton_Pt-500/PhaseIITDRSpring17MiniAOD-PU200CalAging4500Ultimate_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-500/PhaseIITDRSpring17MiniAOD-PU200CalAging300_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-500/PhaseIITDRSpring17MiniAOD-PU200CalAging3000_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-500/PhaseIITDRSpring17MiniAOD-PU200CalAging1000_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    '/SinglePhoton_Pt-1000/PhaseIITDRSpring17MiniAOD-noPUCaloAging4500Ultimate_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-1000/PhaseIITDRSpring17MiniAOD-noPUCaloAging300_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-1000/PhaseIITDRSpring17MiniAOD-noPUCaloAging3000_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-1000/PhaseIITDRSpring17MiniAOD-noPUCaloAging1000_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    '/SinglePhoton_Pt-1000/PhaseIITDRSpring17MiniAOD-PU200CalAging4500Ultimate_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-1000/PhaseIITDRSpring17MiniAOD-PU200CalAging300_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-1000/PhaseIITDRSpring17MiniAOD-PU200CalAging3000_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-1000/PhaseIITDRSpring17MiniAOD-PU200CalAging1000_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    '/SinglePhoton_Pt-2000/PhaseIITDRSpring17MiniAOD-noPUCaloAging4500Ultimate_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-2000/PhaseIITDRSpring17MiniAOD-noPUCaloAging300_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-2000/PhaseIITDRSpring17MiniAOD-noPUCaloAging3000_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-2000/PhaseIITDRSpring17MiniAOD-noPUCaloAging1000_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    '/SinglePhoton_Pt-2000/PhaseIITDRSpring17MiniAOD-PU200CalAging4500Ultimate_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-2000/PhaseIITDRSpring17MiniAOD-PU200CalAging300_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-2000/PhaseIITDRSpring17MiniAOD-PU200CalAging3000_91X_upgrade2023_realistic_v3-v3/MINIAODSIM',
    '/SinglePhoton_Pt-2000/PhaseIITDRSpring17MiniAOD-PU200CalAging1000_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
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
        config.General.requestName = 'phase2photons_%d_%s' % (i, primaryDS)
        config.Data.outputDatasetTag = conditions
        config.Data.inputDataset = pd
        submit(config)
