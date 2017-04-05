// -*- C++ -*-
//
// Package:    EgammaWork/ElectronNtupler
// Class:      SimpleElectronNtupler
// 
/**\class SimpleElectronNtupler SimpleElectronNtupler.cc EgammaWork/ElectronNtupler/plugins/SimpleElectronNtupler.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ilya Kravchenko
//         Created:  Thu, 10 Jul 2014 09:54:13 GMT
//
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "Math/VectorUtil.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"

#include "RecoEgamma/Phase2InterimID/interface/HGCalIDTool.h"
//
// class declaration
//

class SimpleElectronNtupler : public edm::EDAnalyzer {
   public:
      explicit SimpleElectronNtupler(const edm::ParameterSet&);
      ~SimpleElectronNtupler();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  enum ElectronMatchType {UNMATCHED = 0, 
			  TRUE_PROMPT_ELECTRON, 
			  TRUE_ELECTRON_FROM_TAU,
			  TRUE_NON_PROMPT_ELECTRON}; // The last does not include tau parents

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      int matchToTruth(const edm::Ptr<reco::GsfElectron> el,  
		       const edm::Handle<edm::View<reco::GenParticle>>  &prunedGenParticles);
      void findFirstNonElectronMother(const reco::Candidate *particle, int &ancestorPID, int &ancestorStatus);
  
      // ----------member data ---------------------------
      // Data members that are the same for AOD and miniAOD
      edm::EDGetTokenT<edm::View<PileupSummaryInfo> > pileupToken_;
      edm::EDGetTokenT<double> rhoToken_;
      edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
      edm::EDGetTokenT<GenEventInfoProduct> genEventInfoProduct_;

      // AOD case data members
      edm::EDGetToken electronsToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;
      edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;

      // MiniAOD case data members
      edm::EDGetToken electronsMiniAODToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxMiniAODToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesMiniAODToken_;
      edm::EDGetTokenT<reco::ConversionCollection> conversionsMiniAODToken_;

  TTree *electronTree_;

  // Vars for PVs
  Int_t pvNTracks_;

  // Vars for weight (can be negative)
  Float_t genWeight_;

  // Vars for pile-up
  Int_t nPUTrue_;    // true pile-up
  Int_t nPU_;        // generated pile-up
  Int_t nPV_;        // number of reconsrtucted primary vertices
  Float_t rho_;      // the rho variable

  // all electron variables
  Int_t nElectrons_;

  std::vector<Float_t> pt_;
  std::vector<Float_t> etaSC_;
  std::vector<Float_t> phiSC_;
  std::vector<Float_t> dEtaIn_;
  std::vector<Float_t> dEtaSeed_;
  std::vector<Float_t> dPhiIn_;
  std::vector<Float_t> hOverE_;
  std::vector<Float_t> full5x5_sigmaIetaIeta_;
  std::vector<Float_t> isoChargedHadrons_;
  std::vector<Float_t> isoNeutralHadrons_;
  std::vector<Float_t> isoPhotons_;
  std::vector<Float_t> isoChargedFromPU_;
  std::vector<Float_t> relCombIsoWithEA_;
  std::vector<Float_t> ooEmooP_;
  std::vector<Float_t> d0_;
  std::vector<Float_t> dz_;
  std::vector<Int_t>   expectedMissingInnerHits_;
  std::vector<Int_t>   passConversionVeto_;     
  std::vector<Int_t>   isTrue_;
  // addl phase2
  std::vector<Int_t>   isTrackerDriven_;
  std::vector<Int_t>   isEcalDriven_;
  std::vector<Float_t> scR9_;
  std::vector<float>   eSuperClusterOverP ;        // the supercluster energy / track momentum at the PCA to the beam spot
  std::vector<float>   eSeedClusterOverP ;         // the seed cluster energy / track momentum at the PCA to the beam spot
  std::vector<float>   eSeedClusterOverPout ;      // the seed cluster energy / track momentum at calo extrapolated from the outermost track state
  std::vector<float>   eEleClusterOverPout ;       // the electron cluster energy / track momentum at calo extrapolated from the outermost track state
  std::vector<float>   deltaEtaSuperClusterAtVtx ; // the supercluster eta - track eta position at calo extrapolated from innermost track state
  std::vector<float>   deltaEtaSeedClusterAtCalo ; // the seed cluster eta - track eta position at calo extrapolated from the outermost track state
  std::vector<float>   deltaEtaEleClusterAtCalo ;  // the electron cluster eta - track eta position at calo extrapolated from the outermost state
  std::vector<float>   deltaPhiEleClusterAtCalo ;  // the electron cluster phi - track phi position at calo extrapolated from the outermost track state
  std::vector<float>   deltaPhiSuperClusterAtVtx ; // the supercluster phi - track phi position at calo extrapolated from the innermost track state
  std::vector<float>   deltaPhiSeedClusterAtCalo ; // the seed cluster phi - track phi position at calo extrapolated from the outermost track state
  std::vector<float>   hOverE_hgcalSafe_;
  std::vector<float>   hgcId_startPosition_;
  std::vector<float>   hgcId_sigmaietaieta_;
  std::vector<float>   hgcId_lengthCompatibility_;

  EffectiveAreas   effectiveAreas_;

  std::unique_ptr<HGCalIDTool> hgcEmId_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SimpleElectronNtupler::SimpleElectronNtupler(const edm::ParameterSet& iConfig):
  effectiveAreas_( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath() )
{

  //
  // Prepare tokens for all input collections and objects
  //

  // Universal tokens for AOD and miniAOD  
  genEventInfoProduct_ = consumes<GenEventInfoProduct> 
    (iConfig.getParameter <edm::InputTag>
     ("genEventInfoProduct"));

  pileupToken_ = consumes<edm::View<PileupSummaryInfo> >
    (iConfig.getParameter <edm::InputTag>
     ("pileup"));

  rhoToken_    = consumes<double> 
    (iConfig.getParameter <edm::InputTag>
     ("rho"));

  beamSpotToken_    = consumes<reco::BeamSpot> 
    (iConfig.getParameter <edm::InputTag>
     ("beamSpot"));


  // AOD tokens
  electronsToken_    = mayConsume<edm::View<reco::GsfElectron> >
    (iConfig.getParameter<edm::InputTag>
     ("electrons"));

  vtxToken_          = mayConsume<reco::VertexCollection>
    (iConfig.getParameter<edm::InputTag>
     ("vertices"));

  genParticlesToken_ = mayConsume<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticles"));

  conversionsToken_ = mayConsume< reco::ConversionCollection >
    (iConfig.getParameter<edm::InputTag>
     ("conversions"));

  // MiniAOD tokens
  // For electrons, use the fact that pat::Electron can be cast into 
  // GsfElectron
  electronsMiniAODToken_    = mayConsume<edm::View<reco::GsfElectron> >
    (iConfig.getParameter<edm::InputTag>
     ("electronsMiniAOD"));

  vtxMiniAODToken_          = mayConsume<reco::VertexCollection>
    (iConfig.getParameter<edm::InputTag>
     ("verticesMiniAOD"));

  genParticlesMiniAODToken_ = mayConsume<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticlesMiniAOD"));

  conversionsMiniAODToken_ = mayConsume< reco::ConversionCollection >
    (iConfig.getParameter<edm::InputTag>
     ("conversionsMiniAOD"));

  const edm::ParameterSet& hgcIdCfg = iConfig.getParameterSet("HGCalIDToolConfig");
  auto cc = consumesCollector();
  hgcEmId_.reset( new HGCalIDTool(hgcIdCfg, cc) );

  //
  // Set up the ntuple structure
  //
  edm::Service<TFileService> fs;
  electronTree_ = fs->make<TTree> ("ElectronTree", "Electron data");
  
  electronTree_->Branch("pvNTracks"    ,  &pvNTracks_ , "pvNTracks/I");

  electronTree_->Branch("nPV"        ,  &nPV_     , "nPV/I");
  electronTree_->Branch("nPU"        ,  &nPU_     , "nPU/I");
  electronTree_->Branch("nPUTrue"    ,  &nPUTrue_ , "nPUTrue/I");
  electronTree_->Branch("rho"        ,  &rho_ , "rho/F");

  electronTree_->Branch("genWeight"    ,  &genWeight_ , "genWeight/F");
  electronTree_->Branch("nEle"    ,  &nElectrons_ , "nEle/I");
  electronTree_->Branch("pt"    ,  &pt_    );
  electronTree_->Branch("etaSC" ,  &etaSC_ );
  electronTree_->Branch("phiSC" ,  &phiSC_ );
  electronTree_->Branch("dEtaIn",  &dEtaIn_);
  electronTree_->Branch("dEtaSeed",  &dEtaSeed_);
  electronTree_->Branch("dPhiIn",  &dPhiIn_);
  electronTree_->Branch("hOverE",  &hOverE_);
  electronTree_->Branch("full5x5_sigmaIetaIeta", &full5x5_sigmaIetaIeta_);
  electronTree_->Branch("isoChargedHadrons"      , &isoChargedHadrons_);
  electronTree_->Branch("isoNeutralHadrons"      , &isoNeutralHadrons_);
  electronTree_->Branch("isoPhotons"             , &isoPhotons_);
  electronTree_->Branch("relCombIsoWithEA"       , &relCombIsoWithEA_);
  electronTree_->Branch("isoChargedFromPU"       , &isoChargedFromPU_);
  electronTree_->Branch("ooEmooP", &ooEmooP_);
  electronTree_->Branch("d0"     , &d0_);
  electronTree_->Branch("dz"     , &dz_);
  electronTree_->Branch("expectedMissingInnerHits", &expectedMissingInnerHits_);
  electronTree_->Branch("passConversionVeto", &passConversionVeto_);
  electronTree_->Branch("isTrue"    , &isTrue_);
  electronTree_->Branch("isTrackerDriven"    , &isTrackerDriven_);
  electronTree_->Branch("isEcalDriven"    , &isEcalDriven_);
  electronTree_->Branch("scR9"    , &scR9_);
  electronTree_->Branch("eSuperClusterOverP", &eSuperClusterOverP );
  electronTree_->Branch("eSeedClusterOverP", &eSeedClusterOverP );
  electronTree_->Branch("eSeedClusterOverPout", &eSeedClusterOverPout );
  electronTree_->Branch("eEleClusterOverPout", &eEleClusterOverPout );
  electronTree_->Branch("deltaEtaSuperClusterAtVtx", &deltaEtaSuperClusterAtVtx );
  electronTree_->Branch("deltaEtaSeedClusterAtCalo", &deltaEtaSeedClusterAtCalo );
  electronTree_->Branch("deltaEtaEleClusterAtCalo", &deltaEtaEleClusterAtCalo );
  electronTree_->Branch("deltaPhiEleClusterAtCalo", &deltaPhiEleClusterAtCalo );
  electronTree_->Branch("deltaPhiSuperClusterAtVtx", &deltaPhiSuperClusterAtVtx );
  electronTree_->Branch("deltaPhiSeedClusterAtCalo", &deltaPhiSeedClusterAtCalo );
  electronTree_->Branch("hOverE_hgcalSafe", &hOverE_hgcalSafe_);
  electronTree_->Branch("hgcId_startPosition", &hgcId_startPosition_);
  electronTree_->Branch("hgcId_sigmaietaieta", &hgcId_sigmaietaieta_);
  electronTree_->Branch("hgcId_lengthCompatibility", &hgcId_lengthCompatibility_);
 
}


SimpleElectronNtupler::~SimpleElectronNtupler()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SimpleElectronNtupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;
  
  // Get gen weight info
  edm::Handle< GenEventInfoProduct > genWeightH;
  iEvent.getByToken(genEventInfoProduct_,genWeightH);
  genWeight_ = genWeightH->GenEventInfoProduct::weight();

  // Get Pileup info
  Handle<edm::View<PileupSummaryInfo> > pileupHandle;
  iEvent.getByToken(pileupToken_, pileupHandle);
  for( auto & puInfoElement : *pileupHandle){
    if( puInfoElement.getBunchCrossing() == 0 ){
      nPU_    = puInfoElement.getPU_NumInteractions();
      nPUTrue_= puInfoElement.getTrueNumInteractions();
    }
  }
  
  // Get rho value
  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  rho_ = *rhoH;

  // Get the beam spot
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(beamSpotToken_,theBeamSpot);  
  
  // Retrieve the collection of electrons from the event.
  // If we fail to retrieve the collection with the standard AOD
  // name, we next look for the one with the stndard miniAOD name.
  //   We use exactly the same handle for AOD and miniAOD formats
  // since pat::Electron objects can be recast as reco::GsfElectron objects.
  edm::Handle<edm::View<reco::GsfElectron> > electrons;
  bool isAOD = true;
  iEvent.getByToken(electronsToken_, electrons);
  if( !electrons.isValid() ){
    isAOD = false;
    iEvent.getByToken(electronsMiniAODToken_,electrons);
  }
  
  // Get the MC collection
  Handle<edm::View<reco::GenParticle> > genParticles;
  if( isAOD )
    iEvent.getByToken(genParticlesToken_,genParticles);
  else
    iEvent.getByToken(genParticlesMiniAODToken_,genParticles);
  
  // Get PV
  edm::Handle<reco::VertexCollection> vertices;
  if( isAOD )
    iEvent.getByToken(vtxToken_, vertices);
  else
    iEvent.getByToken(vtxMiniAODToken_, vertices);
  
  if (vertices->empty()) return; // skip the event if no PV found
  const reco::Vertex &pv = vertices->front();
  nPV_    = vertices->size();
  
  // NOTE FOR RUN 2 THE OLD SELECTION OF GOOD VERTICES BELOW IS DISCOURAGED
  // // Find the first vertex in the collection that passes
  // // good quality criteria
  // VertexCollection::const_iterator firstGoodVertex = vertices->end();
  // int firstGoodVertexIdx = 0;
  // for (VertexCollection::const_iterator vtx = vertices->begin(); 
  //      vtx != vertices->end(); ++vtx, ++firstGoodVertexIdx) {
  //   // Replace isFake() for miniAOD because it requires tracks and miniAOD vertices don't have tracks:
  //   // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
  //   bool isFake = vtx->isFake();
  //   if( !isAOD )
  //     isFake = (vtx->chi2()==0 && vtx->ndof()==0);
  //   // Check the goodness
  //   if ( !isFake
  // 	 &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0
  // 	 && fabs(vtx->position().Z())<=24.0) {
  //     firstGoodVertex = vtx;
  //     break;
  //   }
  // }
  
  // if ( firstGoodVertex==vertices->end() )
  //   return; // skip event if there are no good PVs
  
  // // Seems always zero. Not stored in miniAOD...?
  // pvNTracks_ = firstGoodVertex->nTracks();
  pvNTracks_ = pv.nTracks();
  

  // Get the conversions collection
  edm::Handle<reco::ConversionCollection> conversions;
  if(isAOD)
    iEvent.getByToken(conversionsToken_, conversions);
  else
    iEvent.getByToken(conversionsMiniAODToken_, conversions);

  hgcEmId_->getEventSetup(iSetup);
  hgcEmId_->getEvent(iEvent);

  // Loop over electrons
  nElectrons_ = 0;
  pt_.clear();
  etaSC_.clear();
  phiSC_.clear();
  dEtaIn_.clear();
  dEtaSeed_.clear();
  dPhiIn_.clear();
  hOverE_.clear();
  full5x5_sigmaIetaIeta_.clear();
  isoChargedHadrons_.clear();
  isoNeutralHadrons_.clear();
  isoPhotons_.clear();
  relCombIsoWithEA_.clear();
  isoChargedFromPU_.clear();
  ooEmooP_.clear();
  d0_.clear();
  dz_.clear();
  expectedMissingInnerHits_.clear();
  passConversionVeto_.clear();     
  isTrue_.clear();
  isTrackerDriven_.clear();
  isEcalDriven_.clear();
  scR9_.clear();
  eSuperClusterOverP.clear( );
  eSeedClusterOverP.clear( );
  eSeedClusterOverPout.clear( );
  eEleClusterOverPout.clear( );
  deltaEtaSuperClusterAtVtx.clear( );
  deltaEtaSeedClusterAtCalo.clear( );
  deltaEtaEleClusterAtCalo.clear( );
  deltaPhiEleClusterAtCalo.clear( );
  deltaPhiSuperClusterAtVtx.clear( );
  deltaPhiSeedClusterAtCalo.clear( );
  hOverE_hgcalSafe_.clear();
  hgcId_startPosition_.clear();
  hgcId_sigmaietaieta_.clear();
  hgcId_lengthCompatibility_.clear();

  for (size_t i = 0; i < electrons->size(); ++i){
    const auto el = electrons->ptrAt(i);
  // for (const pat::Electron &el : *electrons) {
    
    // Kinematics
    if( el->pt() < 10 ) // keep only electrons above 10 GeV
      continue;
    
    nElectrons_++;
    pt_.push_back( el->pt() );
    etaSC_.push_back( el->superCluster()->eta() );
    phiSC_.push_back( el->superCluster()->phi() );

    // Phase2
    isTrackerDriven_.push_back( el->trackerDrivenSeed() );
    isEcalDriven_.push_back( el->ecalDrivenSeed() );
    scR9_.push_back( el->r9() );
    eSuperClusterOverP.push_back( el->eSuperClusterOverP() );
    eSeedClusterOverP.push_back( el->eSeedClusterOverP() );
    eSeedClusterOverPout.push_back( el->eSeedClusterOverPout() );
    eEleClusterOverPout.push_back( el->eEleClusterOverPout() );
    deltaEtaSuperClusterAtVtx.push_back( el->deltaEtaSuperClusterTrackAtVtx() );
    deltaEtaSeedClusterAtCalo.push_back( el->deltaEtaSeedClusterTrackAtCalo() );
    deltaEtaEleClusterAtCalo.push_back( el->deltaEtaEleClusterTrackAtCalo() );
    deltaPhiEleClusterAtCalo.push_back( el->deltaPhiEleClusterTrackAtCalo() );
    deltaPhiSuperClusterAtVtx.push_back( el->deltaPhiSuperClusterTrackAtVtx() );
    deltaPhiSeedClusterAtCalo.push_back( el->deltaPhiSeedClusterTrackAtCalo() );

    if ( el->superCluster().isNonnull() && el->superCluster()->seed().isNonnull() ) {
      const reco::CaloCluster * cluster = el->superCluster()->seed().get();
      bool isHGCal = hgcEmId_->setClusterPtr(cluster);
      if ( isHGCal ) {
        hOverE_hgcalSafe_.push_back( hgcEmId_->getHadronFraction() );
        hgcId_startPosition_.push_back( std::abs(hgcEmId_->getStartPosition().z()) );
        hgcId_sigmaietaieta_.push_back( hgcEmId_->getSigmaEtaEta() );
        hgcId_lengthCompatibility_.push_back( hgcEmId_->getLengthCompatibility() );
      }
      else {
        hOverE_hgcalSafe_.push_back( el->hcalOverEcal() );
        hgcId_startPosition_.push_back( -1. );
        hgcId_sigmaietaieta_.push_back( el->full5x5_sigmaIetaIeta() );
        hgcId_lengthCompatibility_.push_back( -1. );
      }
    }
    else {
      hOverE_hgcalSafe_.push_back( -1. );
      hgcId_startPosition_.push_back( -1. );
      hgcId_sigmaietaieta_.push_back( -1. );
      hgcId_lengthCompatibility_.push_back( -1. );
    }
    
    // ID and matching
    dEtaIn_.push_back( el->deltaEtaSuperClusterTrackAtVtx() );
    // Calculation of dEtaSeed is taken from VID (by HEEP folks)
    //   https://github.com/cms-sw/cmssw/blob/CMSSW_8_1_X/RecoEgamma/ElectronIdentification/plugins/cuts/GsfEleDEtaInSeedCut.cc#L31-L32
    float dEtaSeedValue = 
      el->superCluster().isNonnull() && el->superCluster()->seed().isNonnull() 
      ? 
      el->deltaEtaSuperClusterTrackAtVtx() 
      - el->superCluster()->eta() 
      + el->superCluster()->seed()->eta() 
      : std::numeric_limits<float>::max();
    dEtaSeed_.push_back( dEtaSeedValue );
    dPhiIn_.push_back( el->deltaPhiSuperClusterTrackAtVtx() );
    hOverE_.push_back( el->hcalOverEcal() );
    full5x5_sigmaIetaIeta_.push_back( el->full5x5_sigmaIetaIeta() );
    // |1/E-1/p| = |1/E - EoverPinner/E| is computed below
    // The if protects against ecalEnergy == inf or zero
    // (always the case for miniAOD for electrons <5 GeV)
    if( el->ecalEnergy() == 0 ){
      printf("Electron energy is zero!\n");
      ooEmooP_.push_back( 1e30 );
    }else if( !std::isfinite(el->ecalEnergy())){
      printf("Electron energy is not finite!\n");
      ooEmooP_.push_back( 1e30 );
    }else{
      ooEmooP_.push_back( fabs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy() ) );
    }
    
    // Isolation
    GsfElectron::PflowIsolationVariables pfIso = el->pfIsolationVariables();
    // Compute individual PF isolations
    isoChargedHadrons_.push_back( pfIso.sumChargedHadronPt );
    isoNeutralHadrons_.push_back( pfIso.sumNeutralHadronEt );
    isoPhotons_.push_back( pfIso.sumPhotonEt );
    isoChargedFromPU_.push_back( pfIso.sumPUPt );

    // Compute combined relative PF isolation with the effective area correction for pile-up
    float abseta =  abs(el->superCluster()->eta());
    float eA = effectiveAreas_.getEffectiveArea(abseta);
    relCombIsoWithEA_.push_back( ( pfIso.sumChargedHadronPt
				   + std::max( 0.0f, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - eA*rho_) )
				 / el->pt() );

    // Impact parameter
    reco::GsfTrackRef theTrack = el->gsfTrack();
    d0_.push_back( (-1) * theTrack->dxy(pv.position() ) );
    dz_.push_back( theTrack->dz( pv.position() ) );
    // d0_.push_back( (-1) * theTrack->dxy(firstGoodVertex->position() ) );
    // dz_.push_back( theTrack->dz( firstGoodVertex->position() ) );

    // Conversion rejection
    expectedMissingInnerHits_.push_back(el->gsfTrack()->hitPattern()
					.numberOfHits(reco::HitPattern::MISSING_INNER_HITS) );

    bool passConvVeto = !ConversionTools::hasMatchedConversion(*el, 
							       conversions,
							       theBeamSpot->position());
    passConversionVeto_.push_back( (int) passConvVeto );
    
    // Match to generator level truth
    
    isTrue_.push_back( matchToTruth( el, genParticles) );
    
  }
  
  // Save this electron's info
  electronTree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
SimpleElectronNtupler::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SimpleElectronNtupler::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
SimpleElectronNtupler::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
SimpleElectronNtupler::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
SimpleElectronNtupler::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
SimpleElectronNtupler::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SimpleElectronNtupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

int SimpleElectronNtupler::matchToTruth(const edm::Ptr<reco::GsfElectron> el, 
				  const edm::Handle<edm::View<reco::GenParticle>> &prunedGenParticles){

  // 
  // Explicit loop and geometric matching method (advised by Josh Bendavid)
  //

  // Find the closest status 1 gen electron to the reco electron
  double dR = 999;
  const reco::Candidate *closestElectron = 0;
  for(size_t i=0; i<prunedGenParticles->size();i++){
    const reco::Candidate *particle = &(*prunedGenParticles)[i];
    // Drop everything that is not electron or not status 1
    if( abs(particle->pdgId()) != 11 || particle->status() != 1 )
      continue;
    //
    double dRtmp = ROOT::Math::VectorUtil::DeltaR( el->p4(), particle->p4() );
    if( dRtmp < dR ){
      dR = dRtmp;
      closestElectron = particle;
    }
  }
  // See if the closest electron (if it exists) is close enough.
  // If not, no match found.
  if( !(closestElectron != 0 && dR < 0.1) ) {
    return UNMATCHED;
  }

  // 
  int ancestorPID = -999; 
  int ancestorStatus = -999;
  findFirstNonElectronMother(closestElectron, ancestorPID, ancestorStatus);

  if( ancestorPID == -999 && ancestorStatus == -999 ){
    // No non-electron parent??? This should never happen.
    // Complain.
    printf("SimpleElectronNtupler: ERROR! Electron does not apper to have a non-electron parent\n");
    return UNMATCHED;
  }
  
  if( abs(ancestorPID) > 50 && ancestorStatus == 2 )
    return TRUE_NON_PROMPT_ELECTRON;

  if( abs(ancestorPID) == 15 && ancestorStatus == 2 )
    return TRUE_ELECTRON_FROM_TAU;

  // What remains is true prompt electrons
  return TRUE_PROMPT_ELECTRON;
}

void SimpleElectronNtupler::findFirstNonElectronMother(const reco::Candidate *particle,
						 int &ancestorPID, int &ancestorStatus){

  if( particle == 0 ){
    printf("SimpleElectronNtupler: ERROR! null candidate pointer, this should never happen\n");
    return;
  }

  // Is this the first non-electron parent? If yes, return, otherwise
  // go deeper into recursion
  if( abs(particle->pdgId()) == 11 ){
    findFirstNonElectronMother(particle->mother(0), ancestorPID, ancestorStatus);
  }else{
    ancestorPID = particle->pdgId();
    ancestorStatus = particle->status();
  }

  return;
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimpleElectronNtupler);
