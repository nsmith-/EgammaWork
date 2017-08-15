import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("Ntupler")
options = VarParsing('analysis')
options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.load('Configuration.Geometry.GeometryExtended2023D4Reco_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.source = cms.Source ("PoolSource", fileNames = cms.untracked.vstring(options.inputFiles) )
#
# Configure the ntupler module
#

process.ntupler = cms.EDAnalyzer('SimpleElectronNtupler',
                                 #
                                 # Common to all formats objects
                                 #
                                 pileup   = cms.InputTag("slimmedAddPileupInfo"),
                                 rho      = cms.InputTag("fixedGridRhoFastjetAll"),
                                 beamSpot = cms.InputTag('offlineBeamSpot'),
                                 genEventInfoProduct = cms.InputTag('generator'),
                                 #
                                 # Objects specific to MiniAOD format
                                 #
                                 electronsMiniAOD    = cms.InputTag("slimmedElectrons"),
                                 genParticlesMiniAOD = cms.InputTag("prunedGenParticles"),
                                 verticesMiniAOD     = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                 conversionsMiniAOD  = cms.InputTag('reducedEgamma:reducedConversions'),
                                 # Effective areas for computing PU correction for isolations
                                 effAreasConfigFile = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Spring15/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt"),
                                 )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.outputFile)
                                   )


#process.p = cms.Path(process.electronTrackIsolationLcone+process.particleFlowRecHitHGCSeq+process.ntupler)
process.p = cms.Path(process.ntupler)
