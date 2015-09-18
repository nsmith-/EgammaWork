import FWCore.ParameterSet.Config as cms

process = cms.Process("TestElectrons")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.Geometry_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# NOTE: the pick the right global tag!
#    for Spring15 50ns MC: global tag is 'auto:run2_mc_50ns'
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
    # AOD test files from /DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/AODSIM
    'root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/0033A97B-8707-E511-9D3B-008CFA1980B8.root',
    'root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/003826EE-8807-E511-9628-008CFA06470C.root',
    'root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/025EA7DE-9907-E511-B054-002590593876.root'
    )    

inputFilesMiniAOD = cms.untracked.vstring(
    # MiniAOD test files from /DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM
    'root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/00C4781D-6B08-E511-8A0A-0025905A6084.root',
    'root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/02DE3B74-6C08-E511-ABE3-0025905A60D0.root',
    'root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/08CEA172-6C08-E511-8097-0025905A6090.root'
    )

# Set up input/output depending on the format
# You can list here either AOD or miniAOD files, but not both types mixed
#
useAOD = False

if useAOD == True :
    inputFiles = inputFilesAOD
    outputFile = "electron_ntuple.root"
    print("AOD input files are used")
else :
    inputFiles = inputFilesMiniAOD
    outputFile = "electron_ntuple_mini.root"
    print("MiniAOD input files are used")
process.source = cms.Source ("PoolSource", fileNames = inputFiles )                             

#
# Set up electron ID (VID framework)
#

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
if useAOD == True :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_50ns_Trig_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

#
# Configure an example module for user analysis with electrons
#

process.ntupler = cms.EDAnalyzer('ElectronNtuplerVIDwithMVADemo',
                                 # The module automatically detects AOD vs miniAOD, so we configure both
                                 #
                                 # Common to all formats objects
                                 #
                                 # ... none ...
                                 #
                                 # Objects specific to AOD format
                                 #
                                 electrons    = cms.InputTag("gedGsfElectrons"),
                                 genParticles = cms.InputTag("genParticles"),
                                 #
                                 # Objects specific to MiniAOD format
                                 #
                                 electronsMiniAOD    = cms.InputTag("slimmedElectrons"),
                                 genParticlesMiniAOD = cms.InputTag("prunedGenParticles"),
                                 #
                                 # ID decisions (common to all formats)
                                 eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-50ns-Trig-V1-wp90"),
                                 eleTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-50ns-Trig-V1-wp80"),
                                 #
                                 # ValueMaps with MVA results
                                 #
                                 mvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig50nsV1Values"),
                                 mvaCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig50nsV1Categories")
                                )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( outputFile )
                                   )

# Make sure to add the ID sequence upstream from the user analysis module
process.p = cms.Path(process.egmGsfElectronIDSequence * process.ntupler)
