from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.workArea = 'crab_phase2photons'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runPhotons_phase2.py'

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 4
config.Data.outLFNDirBase = '/store/user/%s/phase2photons' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.allowNonValidInputDataset = True
config.Data.totalUnits = 20

#config.Site.storageSite = 'T2_US_Wisconsin'
config.Site.storageSite = 'T2_CH_CERN'

# 7 pt * 3 pu * 4 aging = 84
pds_singlepho = [
    '/SinglePhoton_Pt-20/PhaseIITDRSpring17MiniAOD-PU140CaloAging1000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-20/PhaseIITDRSpring17MiniAOD-PU140CaloAging3000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-20/PhaseIITDRSpring17MiniAOD-PU140CaloAging300NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-20/PhaseIITDRSpring17MiniAOD-PU140CaloAging4500UltimateNewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-20/PhaseIITDRSpring17MiniAOD-PU200CalAging1000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-20/PhaseIITDRSpring17MiniAOD-PU200CalAging3000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-20/PhaseIITDRSpring17MiniAOD-PU200CalAging300NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-20/PhaseIITDRSpring17MiniAOD-PU200CalAging4500UltimateNewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-20/PhaseIITDRSpring17MiniAOD-noPUCaloAging1000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-20/PhaseIITDRSpring17MiniAOD-noPUCaloAging3000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-20/PhaseIITDRSpring17MiniAOD-noPUCaloAging300NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-20/PhaseIITDRSpring17MiniAOD-noPUCaloAging4500UltimateNewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-50/PhaseIITDRSpring17MiniAOD-PU140CaloAging1000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-50/PhaseIITDRSpring17MiniAOD-PU140CaloAging3000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-50/PhaseIITDRSpring17MiniAOD-PU140CaloAging300NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-50/PhaseIITDRSpring17MiniAOD-PU140CaloAging4500UltimateNewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-50/PhaseIITDRSpring17MiniAOD-PU200CalAging1000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-50/PhaseIITDRSpring17MiniAOD-PU200CalAging3000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-50/PhaseIITDRSpring17MiniAOD-PU200CalAging300NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-50/PhaseIITDRSpring17MiniAOD-PU200CalAging4500UltimateNewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-50/PhaseIITDRSpring17MiniAOD-noPUCaloAging1000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-50/PhaseIITDRSpring17MiniAOD-noPUCaloAging3000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-50/PhaseIITDRSpring17MiniAOD-noPUCaloAging300NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-50/PhaseIITDRSpring17MiniAOD-noPUCaloAging4500UltimateNewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-100/PhaseIITDRSpring17MiniAOD-PU140CaloAging1000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-100/PhaseIITDRSpring17MiniAOD-PU140CaloAging3000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-100/PhaseIITDRSpring17MiniAOD-PU140CaloAging300NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-100/PhaseIITDRSpring17MiniAOD-PU140CaloAging4500UltimateNewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-100/PhaseIITDRSpring17MiniAOD-PU200CalAging1000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-100/PhaseIITDRSpring17MiniAOD-PU200CalAging3000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-100/PhaseIITDRSpring17MiniAOD-PU200CalAging300NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-100/PhaseIITDRSpring17MiniAOD-PU200CalAging4500UltimateNewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-100/PhaseIITDRSpring17MiniAOD-noPUCaloAging1000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-100/PhaseIITDRSpring17MiniAOD-noPUCaloAging3000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-100/PhaseIITDRSpring17MiniAOD-noPUCaloAging300NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-100/PhaseIITDRSpring17MiniAOD-noPUCaloAging4500UltimateNewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-200/PhaseIITDRSpring17MiniAOD-PU140CaloAging1000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-200/PhaseIITDRSpring17MiniAOD-PU140CaloAging3000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-200/PhaseIITDRSpring17MiniAOD-PU140CaloAging300NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-200/PhaseIITDRSpring17MiniAOD-PU140CaloAging4500UltimateNewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-200/PhaseIITDRSpring17MiniAOD-PU200CalAging1000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-200/PhaseIITDRSpring17MiniAOD-PU200CalAging3000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-200/PhaseIITDRSpring17MiniAOD-PU200CalAging300NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-200/PhaseIITDRSpring17MiniAOD-PU200CalAging4500UltimateNewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-200/PhaseIITDRSpring17MiniAOD-noPUCaloAging1000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-200/PhaseIITDRSpring17MiniAOD-noPUCaloAging3000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-200/PhaseIITDRSpring17MiniAOD-noPUCaloAging300NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-200/PhaseIITDRSpring17MiniAOD-noPUCaloAging4500UltimateNewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-500/PhaseIITDRSpring17MiniAOD-PU140CaloAging1000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-500/PhaseIITDRSpring17MiniAOD-PU140CaloAging3000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-500/PhaseIITDRSpring17MiniAOD-PU140CaloAging300NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-500/PhaseIITDRSpring17MiniAOD-PU140CaloAging4500UltimateNewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-500/PhaseIITDRSpring17MiniAOD-PU200CalAging1000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-500/PhaseIITDRSpring17MiniAOD-PU200CalAging3000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-500/PhaseIITDRSpring17MiniAOD-PU200CalAging300NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-500/PhaseIITDRSpring17MiniAOD-PU200CalAging4500UltimateNewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-500/PhaseIITDRSpring17MiniAOD-noPUCaloAging1000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-500/PhaseIITDRSpring17MiniAOD-noPUCaloAging3000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-500/PhaseIITDRSpring17MiniAOD-noPUCaloAging300NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-500/PhaseIITDRSpring17MiniAOD-noPUCaloAging4500UltimateNewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-1000/PhaseIITDRSpring17MiniAOD-PU140CaloAging1000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-1000/PhaseIITDRSpring17MiniAOD-PU140CaloAging3000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-1000/PhaseIITDRSpring17MiniAOD-PU140CaloAging300NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-1000/PhaseIITDRSpring17MiniAOD-PU140CaloAging4500UltimateNewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-1000/PhaseIITDRSpring17MiniAOD-PU200CalAging1000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-1000/PhaseIITDRSpring17MiniAOD-PU200CalAging3000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-1000/PhaseIITDRSpring17MiniAOD-PU200CalAging300NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-1000/PhaseIITDRSpring17MiniAOD-PU200CalAging4500UltimateNewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-1000/PhaseIITDRSpring17MiniAOD-noPUCaloAging1000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-1000/PhaseIITDRSpring17MiniAOD-noPUCaloAging3000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-1000/PhaseIITDRSpring17MiniAOD-noPUCaloAging300NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-1000/PhaseIITDRSpring17MiniAOD-noPUCaloAging4500UltimateNewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-2000/PhaseIITDRSpring17MiniAOD-PU140CaloAging1000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-2000/PhaseIITDRSpring17MiniAOD-PU140CaloAging3000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-2000/PhaseIITDRSpring17MiniAOD-PU140CaloAging300NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-2000/PhaseIITDRSpring17MiniAOD-PU140CaloAging4500UltimateNewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-2000/PhaseIITDRSpring17MiniAOD-PU200CalAging1000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-2000/PhaseIITDRSpring17MiniAOD-PU200CalAging3000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-2000/PhaseIITDRSpring17MiniAOD-PU200CalAging300NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-2000/PhaseIITDRSpring17MiniAOD-PU200CalAging4500UltimateNewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-2000/PhaseIITDRSpring17MiniAOD-noPUCaloAging1000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-2000/PhaseIITDRSpring17MiniAOD-noPUCaloAging3000NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-2000/PhaseIITDRSpring17MiniAOD-noPUCaloAging300NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/SinglePhoton_Pt-2000/PhaseIITDRSpring17MiniAOD-noPUCaloAging4500UltimateNewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
]

pds_gjets = [
    # 4(+1 ext) aging * 2 pu * 3 phase-space
    # '/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-PU200CaloAging3000NewCT_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    # '/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-PU200CaloAging300NewCT_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    # '/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-PU200CaloAging4500UltimateNewCT_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    # '/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-PU200NewCT_91X_upgrade2023_realistic_v3_ext1-v1/MINIAODSIM',
    # '/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-PU200NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    # '/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-noPUCaloAging3000NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    # '/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-noPUCaloAging300NewCT_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    # '/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-noPUCaloAging4500UltimateNewCT_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    # '/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-noPUNewCT_91X_upgrade2023_realistic_v3_ext1-v1/MINIAODSIM',
    # '/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-noPUNewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    # '/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-PU200CaloAging3000NewCT_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    # '/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-PU200CaloAging300NewCT_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    # '/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-PU200CaloAging4500UltimateNewCT_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    # '/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-PU200NewCT_91X_upgrade2023_realistic_v3_ext1-v1/MINIAODSIM',
    # '/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-PU200NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    # '/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-noPUCaloAging3000NewCT_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    # '/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-noPUCaloAging300NewCT_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    # '/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-noPUCaloAging4500UltimateNewCT_91X_upgrade2023_realistic_v3-v2/MINIAODSIM',
    # '/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-noPUNewCT_91X_upgrade2023_realistic_v3_ext1-v1/MINIAODSIM',
    # '/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-noPUNewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-PU200CaloAging3000NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-PU200CaloAging300NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-PU200CaloAging4500UltimateNewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-PU200NewCT_91X_upgrade2023_realistic_v3_ext1-v1/MINIAODSIM',
    '/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-PU200NewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-noPUCaloAging3000NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-noPUCaloAging300NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-noPUCaloAging4500UltimateNewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    '/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-noPUNewCT_91X_upgrade2023_realistic_v3_ext1-v1/MINIAODSIM',
    '/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_14TeV_Pythia8/PhaseIITDRSpring17MiniAOD-noPUNewCT_NewCT_91X_upgrade2023_realistic_v3-v1/MINIAODSIM',
    # run 2
    # '/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',
    # '/GluGluHToGG_M-125_13TeV_powheg_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',
]

# 5 aging * 3 pu = 15
pds_relvalHgg = [
    ('/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU200ExtEA4500Ult-v1/MINIAODSIM', '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU200ExtEA4500Ult-v1/GEN-SIM-RECO'),
    ('/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU200ExtEA3000Ult-v1/MINIAODSIM', '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU200ExtEA3000Ult-v1/GEN-SIM-RECO'),
    ('/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU200ExtEA3000-v1/MINIAODSIM', '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU200ExtEA3000-v1/GEN-SIM-RECO'),
    ('/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU200ExtEA300-v1/MINIAODSIM', '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU200ExtEA300-v1/GEN-SIM-RECO'),
    ('/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU200ExtEA1000-v1/MINIAODSIM', '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU200ExtEA1000-v1/GEN-SIM-RECO'),
    # ('/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU140ExtEA4500Ult-v1/MINIAODSIM', '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU140ExtEA4500Ult-v1/GEN-SIM-RECO'),
    # ('/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU140ExtEA3000Ult-v1/MINIAODSIM', '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU140ExtEA3000Ult-v1/GEN-SIM-RECO'),
    # ('/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU140ExtEA3000-v1/MINIAODSIM', '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU140ExtEA3000-v1/GEN-SIM-RECO'),
    # ('/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU140ExtEA300-v1/MINIAODSIM', '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU140ExtEA300-v1/GEN-SIM-RECO'),
    # ('/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU140ExtEA1000-v1/MINIAODSIM', '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-PU25ns_91X_upgrade2023_realistic_v3_D17PU140ExtEA1000-v1/GEN-SIM-RECO'),
    ('/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-91X_upgrade2023_realistic_v3_D17noPUExtEA4500Ult-v1/MINIAODSIM', '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-91X_upgrade2023_realistic_v3_D17noPUExtEA4500Ult-v1/GEN-SIM-RECO'),
    ('/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-91X_upgrade2023_realistic_v3_D17noPUExtEA3000Ult-v1/MINIAODSIM', '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-91X_upgrade2023_realistic_v3_D17noPUExtEA3000Ult-v1/GEN-SIM-RECO'),
    ('/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-91X_upgrade2023_realistic_v3_D17noPUExtEA3000-v1/MINIAODSIM', '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-91X_upgrade2023_realistic_v3_D17noPUExtEA3000-v1/GEN-SIM-RECO'),
    ('/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-91X_upgrade2023_realistic_v3_D17noPUExtEA300-v1/MINIAODSIM', '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-91X_upgrade2023_realistic_v3_D17noPUExtEA300-v1/GEN-SIM-RECO'),
    ('/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-91X_upgrade2023_realistic_v3_D17noPUExtEA1000-v1/MINIAODSIM', '/RelValH125GGgluonfusion_14/CMSSW_9_1_2_patch1-91X_upgrade2023_realistic_v3_D17noPUExtEA1000-v1/GEN-SIM-RECO'),
]

pds = pds_gjets + pds_relvalHgg

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
        if type(pd) is tuple:
            pd, second = pd
        else:
            second = None
        (_, primaryDS, conditions, dataTier) = pd.split('/')
        config.General.requestName = 'p2phoID_%d_%s' % (i, primaryDS)
        config.Data.outputDatasetTag = conditions
        config.Data.inputDataset = pd
        if second:
            config.Data.secondaryInputDataset = second
            config.Data.useParent = False
        else:
            config.Data.secondaryInputDataset = None
            config.Data.useParent = True
        submit(config)
