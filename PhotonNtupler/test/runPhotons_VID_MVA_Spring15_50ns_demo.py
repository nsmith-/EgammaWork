import FWCore.ParameterSet.Config as cms

process = cms.Process("TestPhotons")

process.load("FWCore.MessageService.MessageLogger_cfi")

#process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# NOTE: the pick the right global tag!
#    for Spring15 50ns MC: global tag is 'auto:run2_mc_50'
#    for Spring15 25ns MC: global tag is 'auto:run2_mc'
#    for Run 2 data: global tag is 'auto:run2_data'
#  as a rule, find the "auto" global tag in $CMSSW_RELEASE_BASE/src/Configuration/AlCa/python/autoCond.py
#  This auto global tag will look up the "proper" global tag
#  that is typically found in the DAS under the Configs for given dataset
#  (although it can be "overridden" by requirements of a given release)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc_50ns', '')

#
# Define input data to read
#
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

inputFilesAOD = cms.untracked.vstring(
    # AOD test files from 
    # /GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/AODSIM
       'root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/60000/0006B4F9-AB04-E511-B402-6CC2173C9150.root',
       'root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/60000/0034FBCC-7C04-E511-AD7C-00266CFCC214.root',
       'root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/60000/0082354B-8F04-E511-9822-B083FED76E05.root',
    )    

inputFilesMiniAOD = cms.untracked.vstring(
    # MiniAOD test files from 
    # /GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1/MINIAODSIM
       'root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/60000/02FC90C1-4A05-E511-9282-0025905521B2.root',
       'root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/60000/0421F25A-4705-E511-A0AC-E0CB4E19F981.root',
       'root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/60000/04D633A2-CC05-E511-9663-001E68862A32.root',
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
my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_PHYS14_PU20bx25_nonTrig_V1_cff',
                 'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring15_50ns_nonTrig_V2_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

#
# Configure an example module for user analysis with photons
#

process.ntupler = cms.EDAnalyzer(
    'PhotonNtuplerVIDwithMVADemo',
    # The module automatically detects AOD vs miniAOD, so we configure both
    #
    # Common to all formats objects
    #
    # ... none ...
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
    # (the names of the ValueMaps for just decision and full info are the same)
    phoMediumIdBoolMap = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring15-50ns-nonTrig-V2-wp90"),
    phoMediumIdFullInfoMap = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring15-50ns-nonTrig-V2-wp90"),
    # This is a fairly verbose mode if switched on, with full cut flow 
    # diagnostics for each candidate. Use it in a low event count test job.
    phoIdVerbose = cms.bool(False),
    #
    # ValueMaps with MVA results
    #
    mvaValuesMap     = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring15NonTrig50nsV2Values"),
    mvaCategoriesMap = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring15NonTrig50nsV2Categories")
    )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( outputFile )
                                   )

# Make sure to add the ID sequence upstream from the user analysis module
process.p = cms.Path(process.egmPhotonIDSequence * process.ntupler)

