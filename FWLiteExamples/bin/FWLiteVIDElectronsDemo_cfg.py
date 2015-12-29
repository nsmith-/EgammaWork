import FWCore.ParameterSet.Config as cms

process = cms.PSet()

# Define data format
useAOD = False

#define inputs for either AOD or miniAOD
inputFilesAOD = cms.vstring([
    # AOD test files from /DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3/AODSIM
       'root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/002F7FDD-BA13-E511-AA63-0026189437F5.root',
       'root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/00610EE7-C213-E511-842C-00304833529A.root',
       'root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/00623DCC-A813-E511-A302-0025905B85A2.root',
    ])    

inputFilesMiniAOD = cms.vstring([
    # MiniAOD test files from /DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3/MINIAODSIM
    'root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/009D49A5-7314-E511-84EF-0025905A605E.root',
    'root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/00C0BECF-6F14-E511-96F8-0025904B739A.root',
    'root://cms-xrd-global.cern.ch//store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/0260F225-7614-E511-A79F-00A0D1EE8EB4.root',
    ])

# Set up input/output names depending on the format
#
if useAOD == True :
    inputFiles = inputFilesAOD
    outputFile = "electron_ntuple.root"
    print("AOD input files are used")
else :
    inputFiles = inputFilesMiniAOD
    outputFile = "electron_ntuple_mini.root"
    print("MiniAOD input files are used")


process.fwliteInput = cms.PSet(
    fileNames = inputFiles,
    maxEvents   = cms.int32(100),
    outputEvery = cms.uint32(10)
    )

# Define the output ROOT file. It is required, but it won't be 
# created if the user analysis code doesn't create any output objects to save.
process.fwliteOutput = cms.PSet(
    fileName  = cms.string(outputFile) ## mandatory
    )

#
# Create a PSet for the user analysis module
#
process.electronAnalyzer = cms.PSet(
    # define the name of the input electron collection
    electronsAOD     = cms.InputTag("gedGsfElectrons"),
    electronsMiniAOD = cms.InputTag("slimmedElectrons")
    )

#
# START VID ID CONFIGURATION
#

# 
# Load the ID definition python code.
# All certified IDs for electrons are normally found in /RecoEgamma/ElectronIdentification/python/Identification/ area
#
# This example is a cut-based ID tuned on Spring15, 25ns conditions
from RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff import *

# Define which is the specific working point(s) of interest. For the list
# of all working points one can just look at the end of the _cff.py file 
# named above. Note that one can choose to work with just one of the working points,
# configuring the others is not required.

wpLoose  = cutBasedElectronID_Spring15_25ns_V1_standalone_loose
wpMedium = cutBasedElectronID_Spring15_25ns_V1_standalone_medium
wpTight  = cutBasedElectronID_Spring15_25ns_V1_standalone_tight

# The passage below is needed to remove untracked parameters which VID doesn't like.
# This is just technical details, no deeper meaning related to being approved or not.
if hasattr(wpLoose,'isPOGApproved'):
    del wpLoose.isPOGApproved
if hasattr(wpMedium,'isPOGApproved'):
    del wpMedium.isPOGApproved
if hasattr(wpTight,'isPOGApproved'):
    del wpTight.isPOGApproved

#
# The final top-level configuration PSet that defines the ID, passed to the analysis module
#
process.my_vid_configuration = cms.PSet(
    loose = wpLoose,
    medium = wpMedium,
    tight = wpTight
)

#
# END VID ID CONFIGURATION
#
