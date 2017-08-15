// -*- C++ -*-
//
// Package:    EgammaWork/PhotonNtupler
// Class:      SimplePhotonNtupler
// 
/**\class SimplePhotonNtupler SimplePhotonNtupler.cc EgammaWork/PhotonNtupler/plugins/SimplePhotonNtupler.cc

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

#include "DataFormats/EgammaCandidates/interface/Photon.h"

#include "DataFormats/PatCandidates/interface/Photon.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "RecoEgamma/EgammaMCTools/interface/ElectronMCTruth.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruth.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruthFinder.h"

#include "TTree.h"
#include "Math/VectorUtil.h"

//
// class declaration
//

class SimplePhotonNtupler : public edm::EDAnalyzer {
 public:
  explicit SimplePhotonNtupler(const edm::ParameterSet&);
  ~SimplePhotonNtupler();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  enum PhotonMatchType {UNMATCHED = 0, 
			MATCHED_FROM_GUDSCB,
			MATCHED_FROM_PI0,
      MATCHED_FROM_HIGGS,
			MATCHED_FROM_OTHER_SOURCES};
  
 private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  int matchToTruth(const pat::Photon &pho, 
		   const edm::Handle<edm::View<reco::GenParticle>>  &genParticles, float& gen_pt);
  
  void findFirstNonPhotonMother(const reco::Candidate *particle,
				int &ancestorPID, int &ancestorStatus);
  
  // ----------member data ---------------------------

  // Format-independent data members
  edm::EDGetTokenT<double> rhoToken_;
  
  // AOD case data members
  edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;
  
  // MiniAOD case data members
  edm::EDGetToken photonsMiniAODToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesMiniAODToken_;
  
 
  TTree *photonTree_;
  Float_t rho_;      // the rho variable
  
  // all photon variables
  // CAUTION: whatever vectors are declared here, make sure you clear them before the loop over photons!
  Int_t nPhotons_;

  std::vector<Float_t> pt_;
  std::vector<Float_t> eta_;
  std::vector<Float_t> phi_;
  std::vector<Float_t> energy_;
  
  // Energies
  std::vector<Float_t> energy_photons_;
  std::vector<Float_t> sc_energy_;
  std::vector<Float_t> raw_sc_energy_;
  std::vector<Float_t> e3x3_;
  std::vector<Float_t> e5x5_;
  std::vector<Float_t> full5x5_e3x3_;
  std::vector<Float_t> full5x5_e5x5_;

  // Variables typically used for cut based photon ID
  std::vector<Float_t> full5x5_sigmaIetaIeta_;
  std::vector<Float_t> r9_;
  std::vector<Float_t> hOverE_;
  std::vector<Int_t> hasPixelSeed_;
  std::vector<Int_t> conversionSafeElectronVeto_;

  std::vector<Float_t> isoChargedHadrons_;
  std::vector<Float_t> isoNeutralHadrons_;
  std::vector<Float_t> isoPhotons_;

  std::vector<Float_t> isoChargedHadronsWithEA_;
  std::vector<Float_t> isoNeutralHadronsWithEA_;
  std::vector<Float_t> isoPhotonsWithEA_;

  std::vector<Float_t> isoChargedHadronsPuppi_;
  std::vector<Float_t> isoNeutralHadronsPuppi_;
  std::vector<Float_t> isoPhotonsPuppi_;

  std::vector<Int_t> isTrue_;
  std::vector<Float_t> gen_pt_;

  // gen->reco match
  int gtr_nGen_;
  std::vector<float> gtr_gen_pt_;
  std::vector<float> gtr_gen_eta_;
  std::vector<float> gtr_gen_phi_;
  std::vector<float> gtr_gen_energy_;
  std::vector<bool>  gtr_gen_conversion_;
  std::vector<float> gtr_deltaR_;
  std::vector<float> gtr_iReco_;

  // Effective area constants for all isolation types
  EffectiveAreas effAreaChHadrons_;
  EffectiveAreas effAreaNeuHadrons_;
  EffectiveAreas effAreaPhotons_;

  edm::EDGetTokenT<std::vector<SimTrack>> simTracksToken_;
  edm::EDGetTokenT<std::vector<SimVertex>> simVerticesToken_;

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
SimplePhotonNtupler::SimplePhotonNtupler(const edm::ParameterSet& iConfig):
  rhoToken_(consumes<double> (iConfig.getParameter<edm::InputTag>("rho"))),
  // Objects containing effective area constants
  effAreaChHadrons_( (iConfig.getParameter<edm::FileInPath>("effAreaChHadFile")).fullPath() ),
  effAreaNeuHadrons_( (iConfig.getParameter<edm::FileInPath>("effAreaNeuHadFile")).fullPath() ),
  effAreaPhotons_( (iConfig.getParameter<edm::FileInPath>("effAreaPhoFile")).fullPath() ),
  simTracksToken_(consumes<std::vector<SimTrack>>(iConfig.getParameter<edm::InputTag>("simTracksSrc") ) ),
  simVerticesToken_(consumes<std::vector<SimVertex>>(iConfig.getParameter<edm::InputTag>("simVerticesSrc") ) )
{

  //
  // Prepare tokens for all input collections and objects
  //
  
  
  genParticlesToken_ = mayConsume<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticles"));
    
  // MiniAOD tokens
  photonsMiniAODToken_ = mayConsume<edm::View<pat::Photon> >
    (iConfig.getParameter<edm::InputTag>
     ("photonsMiniAOD"));

  genParticlesMiniAODToken_ = mayConsume<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticlesMiniAOD"));
 
  edm::Service<TFileService> fs;
  photonTree_ = fs->make<TTree> ("PhotonTree", "Photon data");
  
  photonTree_->Branch("rho"        ,  &rho_ , "rho/F");
  photonTree_->Branch("nPho",  &nPhotons_ , "nPho/I");

  // Kinematics
  photonTree_->Branch("pt"  ,  &pt_    );
  photonTree_->Branch("eta" ,  &eta_ );
  photonTree_->Branch("phi" ,  &phi_ );
  photonTree_->Branch("energy" ,  &energy_ );

  photonTree_->Branch("energy_photons", &energy_photons_);
  photonTree_->Branch("sc_energy", &sc_energy_);
  photonTree_->Branch("raw_sc_energy", &raw_sc_energy_);
  photonTree_->Branch("e3x3", &e3x3_);
  photonTree_->Branch("e5x5", &e5x5_);
  photonTree_->Branch("full5x5_e3x3", &full5x5_e3x3_);
  photonTree_->Branch("full5x5_e5x5", &full5x5_e5x5_);

  // Variables typically used for cut based photon ID
  photonTree_->Branch("full5x5_sigmaIetaIeta"  , &full5x5_sigmaIetaIeta_);
  photonTree_->Branch("r9"  , &r9_);
  photonTree_->Branch("hOverE"                 ,  &hOverE_);
  photonTree_->Branch("hasPixelSeed"           ,  &hasPixelSeed_);
  photonTree_->Branch("conversionSafeElectronVeto"           ,  &conversionSafeElectronVeto_);

  photonTree_->Branch("isoChargedHadrons"      , &isoChargedHadrons_);
  photonTree_->Branch("isoNeutralHadrons"      , &isoNeutralHadrons_);
  photonTree_->Branch("isoPhotons"             , &isoPhotons_);

  photonTree_->Branch("isoChargedHadronsWithEA"      , &isoChargedHadronsWithEA_);
  photonTree_->Branch("isoNeutralHadronsWithEA"      , &isoNeutralHadronsWithEA_);
  photonTree_->Branch("isoPhotonsWithEA"             , &isoPhotonsWithEA_);

  photonTree_->Branch("isoChargedHadronsPuppi"      , &isoChargedHadronsPuppi_);
  photonTree_->Branch("isoNeutralHadronsPuppi"      , &isoNeutralHadronsPuppi_);
  photonTree_->Branch("isoPhotonsPuppi"             , &isoPhotonsPuppi_);

  photonTree_->Branch("isTrue"             , &isTrue_);

  photonTree_->Branch("nGen", &gtr_nGen_);
  photonTree_->Branch("gen_pt", &gtr_gen_pt_);
  photonTree_->Branch("gen_eta", &gtr_gen_eta_);
  photonTree_->Branch("gen_phi", &gtr_gen_phi_);
  photonTree_->Branch("gen_energy", &gtr_gen_energy_);
  photonTree_->Branch("gen_conversion", &gtr_gen_conversion_);
  photonTree_->Branch("gen_deltaR", &gtr_deltaR_);
  photonTree_->Branch("gen_iReco", &gtr_iReco_);
}


SimplePhotonNtupler::~SimplePhotonNtupler()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SimplePhotonNtupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;

  // Retrieve the collection of photons from the event.
  // If we fail to retrieve the collection with the standard AOD
  // name, we next look for the one with the stndard miniAOD name. 
  //   We use exactly the same handle for AOD and miniAOD formats
  // since pat::Photon objects can be recast as reco::Photon objects.
  edm::Handle<edm::View<pat::Photon> > photons;
  bool isAOD = false;
  iEvent.getByToken(photonsMiniAODToken_,photons);

  // Get generator level info
  Handle<edm::View<reco::GenParticle> > genParticles;
  if( isAOD )
    iEvent.getByToken(genParticlesToken_,genParticles);
  else
    iEvent.getByToken(genParticlesMiniAODToken_,genParticles);

  // Get rho
  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  rho_ = *rhoH;

  Handle<std::vector<SimTrack>> simTracks;
  iEvent.getByToken(simTracksToken_ ,simTracks);
  Handle<std::vector<SimVertex>> simVertices;
  iEvent.getByToken(simVerticesToken_ ,simVertices);
  std::unique_ptr<PhotonMCTruthFinder> thePhotonMCTruthFinder_(new PhotonMCTruthFinder());
  std::vector<PhotonMCTruth> mcPhotons=thePhotonMCTruthFinder_->find(*simTracks,  *simVertices);

  // Clear vectors
  nPhotons_ = 0;
  pt_.clear();
  eta_.clear();
  phi_.clear();
  energy_.clear();

  energy_photons_.clear();
  sc_energy_.clear();
  raw_sc_energy_.clear();
  e3x3_.clear();
  e5x5_.clear();
  full5x5_e3x3_.clear();
  full5x5_e5x5_.clear();

  //
  full5x5_sigmaIetaIeta_.clear();
  r9_.clear();
  hOverE_.clear();
  hasPixelSeed_.clear();
  conversionSafeElectronVeto_.clear();
  //
  isoChargedHadrons_.clear();
  isoNeutralHadrons_.clear();
  isoPhotons_.clear();
  //
  isoChargedHadronsWithEA_.clear();
  isoNeutralHadronsWithEA_.clear();
  isoPhotonsWithEA_.clear();
  //
  isoChargedHadronsPuppi_.clear();
  isoNeutralHadronsPuppi_.clear();
  isoPhotonsPuppi_.clear();
  //
  isTrue_.clear();
  gen_pt_.clear();

  // Loop over photons
  for (size_t i = 0; i < photons->size(); ++i){
    const auto pho = photons->ptrAt(i);

    // Kinematics
    if( pho->pt() < 15 ) 
      continue;
    
    nPhotons_++;

    //
    // Save photon kinematics
    //
    pt_  .push_back( pho->pt() );
    eta_ .push_back( pho->superCluster()->eta() );
    phi_ .push_back( pho->superCluster()->phi() );
    energy_ .push_back( pho->energy() );

    energy_photons_.push_back( pho->p4(pat::Photon::ecal_photons).energy() );
    sc_energy_.push_back( pho->superCluster()->energy() );
    raw_sc_energy_.push_back( pho->superCluster()->rawEnergy() );
    e3x3_.push_back( pho->e3x3() );
    e5x5_.push_back( pho->e5x5() );
    full5x5_e3x3_.push_back( pho->full5x5_e3x3() );
    full5x5_e5x5_.push_back( pho->full5x5_e5x5() );

    hOverE_                .push_back( pho->hadTowOverEm() );
    hasPixelSeed_          .push_back( (Int_t)pho->hasPixelSeed() );
    conversionSafeElectronVeto_.push_back( (Int_t)pho->passElectronVeto() );

    // Get values from ValueMaps, use the photon pointer as the key.
    // Note: starting from CMSSW 7.2.1 or so one can get any full5x5 quantity
    // directly from the photon object, there is no need for value maps anymore.
    // However 7.2.0 and prior (this includes PHYS14 MC samples) requires ValueMaps.
    full5x5_sigmaIetaIeta_ .push_back( pho->full5x5_sigmaIetaIeta() );
    r9_.push_back( pho->r9() );

    float chIso = pho->chargedHadronIso();
    float nhIso = pho->neutralHadronIso();
    float phIso = pho->photonIso();
    isoChargedHadrons_ .push_back( chIso );
    isoNeutralHadrons_ .push_back( nhIso );
    isoPhotons_        .push_back( phIso );

    float abseta = fabs( pho->superCluster()->eta());
    isoChargedHadronsWithEA_ .push_back( std::max( (float)0.0, chIso 
						   - rho_*effAreaChHadrons_.getEffectiveArea(abseta)));
    isoNeutralHadronsWithEA_ .push_back( std::max( (float)0.0, nhIso 
						   - rho_*effAreaNeuHadrons_.getEffectiveArea(abseta)));
    isoPhotonsWithEA_ .push_back( std::max( (float)0.0, phIso 
					    - rho_*effAreaPhotons_.getEffectiveArea(abseta)));

    isoChargedHadronsPuppi_.push_back(pho->puppiChargedHadronIso());
    isoNeutralHadronsPuppi_.push_back(pho->puppiNeutralHadronIso());
    isoPhotonsPuppi_.push_back(pho->puppiPhotonIso());

    // Save MC truth match
    float gen_pt = -1.;
    isTrue_.push_back( matchToTruth(*pho, genParticles, gen_pt) );
    gen_pt_.push_back(gen_pt);

   }
   
  gtr_nGen_ = 0;
  gtr_gen_pt_.clear();
  gtr_gen_eta_.clear();
  gtr_gen_phi_.clear();
  gtr_gen_energy_.clear();
  gtr_gen_conversion_.clear();
  gtr_deltaR_.clear();
  gtr_iReco_.clear();
  for(const auto& p : *genParticles) {
    if ( p.pdgId() == 22 and p.isPromptFinalState() and p.pt() > 10. ) {
      gtr_nGen_++;
      gtr_gen_pt_.push_back(p.pt());
      gtr_gen_eta_.push_back(p.eta());
      gtr_gen_phi_.push_back(p.phi());
      gtr_gen_energy_.push_back(p.energy());

      bool isConversion = false;
      for(auto& pmc : mcPhotons) {
        // auto simTrack = std::find_if(simTracks->begin(), simTracks->end(), [pmc](const SimTrack& t) { return t.trackId() == (size_t) pmc.trackId(); });
        // size_t iGenPart = simTrack->genpartIndex();
        // (but this is genParticles not prunedGenParticles)
        // This is easier than the more correct alternative, first one is always initial G4 track
        if ( reco::deltaR(pmc.fourMomentum(), p.p4()) < 0.001 ) {
          // PhotonMCTruth::vertex() is the conversion vertex (not mother vertex)
          isConversion = pmc.isAConversion() and pmc.vertex().perp() < 80.;
          break;
        }
      }
      gtr_gen_conversion_.push_back(isConversion);

      float minDr = 999.;
      int ireco = -1;
      for(int i=0; i<nPhotons_; ++i) {
        float drTest = reco::deltaR(p.eta(), p.phi(), eta_[i], phi_[i]);
        if ( drTest < minDr ) {
          minDr = drTest;
          ireco = i;
        }
      }
      gtr_deltaR_.push_back(minDr);
      gtr_iReco_.push_back(ireco);
    }
  }

  // Save the info
  photonTree_->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
SimplePhotonNtupler::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SimplePhotonNtupler::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
SimplePhotonNtupler::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
SimplePhotonNtupler::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
SimplePhotonNtupler::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
SimplePhotonNtupler::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SimplePhotonNtupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

int SimplePhotonNtupler::matchToTruth(const pat::Photon &pho, 
				   const edm::Handle<edm::View<reco::GenParticle>>  
				   &genParticles, float& gen_pt)
{
  // 
  // Explicit loop and geometric matching method 
  //

  // Find the closest status 1 gen photon to the reco photon
  double dR = 999;
  const reco::Candidate *closestPhoton = 0;
  for(size_t i=0; i<genParticles->size();i++){
    const reco::Candidate *particle = &(*genParticles)[i];
    // Drop everything that is not photon or not status 1
    if( abs(particle->pdgId()) != 22 || particle->status() != 1 )
      continue;
    //
    double dRtmp = ROOT::Math::VectorUtil::DeltaR( pho.p4(), particle->p4() );
    if( dRtmp < dR ){
      dR = dRtmp;
      closestPhoton = particle;
    }
  }
  // See if the closest photon (if it exists) is close enough.
  // If not, no match found.
  if( !(closestPhoton != 0 && dR < 0.1) ) {
    return UNMATCHED;
  }
  gen_pt = closestPhoton->pt();

  // Find ID of the parent of the found generator level photon match
  int ancestorPID = -999; 
  int ancestorStatus = -999;
  findFirstNonPhotonMother(closestPhoton, ancestorPID, ancestorStatus);

  // Allowed parens: quarks pdgId 1-5, or a gluon 21
  std::vector<int> allowedParents { -1, 1, -2, 2, -3, 3, -4, 4, -5, 5, -21, 21 };
  if( !(std::find(allowedParents.begin(), 
		 allowedParents.end(), ancestorPID)
	!= allowedParents.end()) ){
    // So it is not from g, u, d, s, c, b. Check if it is from pi0 or not. 
    if( abs(ancestorPID) == 111 )
      return MATCHED_FROM_PI0;
    else if( abs(ancestorPID) == 25 )
      return MATCHED_FROM_HIGGS;
    else
      return MATCHED_FROM_OTHER_SOURCES;
  }
  return MATCHED_FROM_GUDSCB;
   
}

void SimplePhotonNtupler::findFirstNonPhotonMother(const reco::Candidate *particle,
						int &ancestorPID, int &ancestorStatus){
  
  if( particle == 0 ){
    printf("SimplePhotonNtupler: ERROR! null candidate pointer, this should never happen\n");
    return;
  }

  // Is this the first non-photon parent? If yes, return, otherwise
  // go deeper into recursion
  if( abs(particle->pdgId()) == 22 and particle->numberOfMothers() > 0 ){
    findFirstNonPhotonMother(particle->mother(0), ancestorPID, ancestorStatus);
  }else{
    ancestorPID = particle->pdgId();
    ancestorStatus = particle->status();
  }
  
  return;
}


//define this as a plug-in
DEFINE_FWK_MODULE(SimplePhotonNtupler);
