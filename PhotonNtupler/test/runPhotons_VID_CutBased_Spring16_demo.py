import FWCore.ParameterSet.Config as cms

process = cms.Process("TestPhotons")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# NOTE: the pick the right global tag!
#    for Spring15 50ns MC: global tag is 'auto:run2_mc_50'
#    for Spring15 25ns MC: global tag is 'auto:run2_mc'
#    for Spring16 MC     : global tag is 'auto:run2_mc'
#    for Run 2 data: global tag is 'auto:run2_data'
#  as a rule, find the "auto" global tag in $CMSSW_RELEASE_BASE/src/Configuration/AlCa/python/autoCond.py
#  This auto global tag will look up the "proper" global tag
#  that is typically found in the DAS under the Configs for given dataset
#  (although it can be "overridden" by requirements of a given release)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

#
# Define input data to read
#
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

inputFilesAOD = cms.untracked.vstring(
    # AOD test files from a GJet PT40 dataset
    # /GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM
       '/store/mc/RunIISpring16DR80/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/001A4AF9-990C-E611-8024-008CFA5D275C.root',
       '/store/mc/RunIISpring16DR80/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/004F2F9D-E40D-E611-B3CA-02163E01646D.root',
       '/store/mc/RunIISpring16DR80/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/022D45D6-EA0A-E611-97B0-002590491B22.root',
       '/store/mc/RunIISpring16DR80/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/02708CED-A10C-E611-826B-008CFA197A60.root',
       '/store/mc/RunIISpring16DR80/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/06B34196-E609-E611-B608-549F35AD8B95.root',
    )    

inputFilesMiniAOD = cms.untracked.vstring(
    # MiniAOD test files from a GJet PT40 dataset
    # /GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM
       '/store/mc/RunIISpring16MiniAODv2/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/001E983A-AC2F-E611-8BFE-A0369F5BD91C.root',
       '/store/mc/RunIISpring16MiniAODv2/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/004C3C6A-962F-E611-AE7F-44A842B4D8BE.root',
       '/store/mc/RunIISpring16MiniAODv2/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/009B95ED-922F-E611-A2D3-44A842B4D8CB.root',
       '/store/mc/RunIISpring16MiniAODv2/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/02239BCA-B72F-E611-9668-0026B951BAB0.root',
       '/store/mc/RunIISpring16MiniAODv2/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/02EB1A11-932F-E611-B031-44A842B2D63E.root',
    )

# Set up input/output depending on the format
# You can list here either AOD or miniAOD files, but not both types mixed
#
useAOD = False

if useAOD == True :
    inputFiles = inputFilesAOD
    outputFile = "photon_ntuple.root"
    print("AOD input files are used")
else :
    inputFiles = inputFilesMiniAOD
    outputFile = "photon_ntuple_mini.root"
    print("MiniAOD input files are used")
process.source = cms.Source ("PoolSource", fileNames = inputFiles ) 

#
# Set up photon ID (VID framework)
#

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate
if useAOD == True :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

switchOnVIDPhotonIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring16_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

#
# Configure an example module for user analysis of photons
#

process.ntupler = cms.EDAnalyzer(
    'PhotonNtuplerVIDDemo',
    # The module automatically detects AOD vs miniAOD, so we configure both
    #
    # Common to all formats objects
    #                                    
    # ... 
    #
    # Objects specific to AOD format
    #
    photons = cms.InputTag("gedPhotons"),
    genParticles = cms.InputTag("genParticles"),
    #
    # Objects specific to MiniAOD format
    #
    photonsMiniAOD = cms.InputTag("slimmedPhotons"),
    genParticlesMiniAOD = cms.InputTag("prunedGenParticles"),
    #
    # ID decisions (common to all formats)
    #
    phoLooseIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V1-loose"),
    phoMediumIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V1-medium"),
    phoTightIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V1-tight"),
    #
    phoMediumIdFullInfoMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V1-medium"),
    # This is a fairly verbose mode if switched on, with full cut flow 
    # diagnostics for each candidate. Use it in a low event count test job.
    phoIdVerbose = cms.bool(False)
    )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( outputFile )
                                   )
process.p = cms.Path(process.egmPhotonIDSequence * process.ntupler)
