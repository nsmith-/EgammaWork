// -*- C++ -*-
//
// Package:    EgammaWork/PhotonNtupler
// Class:      PhotonNtuplerVIDwithMVADemo
// 
/**\class PhotonNtuplerVIDwithMVADemo PhotonNtuplerVIDwithMVADemo.cc EgammaWork/PhotonNtupler/plugins/PhotonNtuplerVIDwithMVADemo.cc

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

#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "Math/VectorUtil.h"

//
// class declaration
//

class PhotonNtuplerVIDwithMVADemo : public edm::EDAnalyzer {
   public:
      explicit PhotonNtuplerVIDwithMVADemo(const edm::ParameterSet&);
      ~PhotonNtuplerVIDwithMVADemo();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      enum PhotonMatchType {UNMATCHED = 0, 
                            MATCHED_FROM_GUDSCB,
                            MATCHED_FROM_PI0,
                            MATCHED_FROM_OTHER_SOURCES};

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      int matchToTruth(const reco::Photon &pho, 
		       const edm::Handle<edm::View<reco::GenParticle>>  &genParticles);

      void findFirstNonPhotonMother(const reco::Candidate *particle,
				    int &ancestorPID, int &ancestorStatus);


      // ----------member data ---------------------------

      // Data members that are the same for AOD and miniAOD
      // ... none ...

      // AOD case data members
      edm::EDGetToken photonsToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;

      // MiniAOD case data members
      edm::EDGetToken photonsMiniAODToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesMiniAODToken_;

      // ID decisions objects
      edm::EDGetTokenT<edm::ValueMap<bool> >               phoMediumIdBoolMapToken_;
      edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > phoMediumIdFullInfoMapToken_;

      // MVA values and categories (optional)
      edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_;
      edm::EDGetTokenT<edm::ValueMap<int> > mvaCategoriesMapToken_;

  TTree *photonTree_;

  // Global info
  Int_t run_;
  Int_t lumi_;
  Int_t evtnum_;

  // all variables for the output tree
  Int_t nPhotons_;

  std::vector<Float_t> pt_;
  std::vector<Float_t> eta_;
  std::vector<Float_t> phi_;

  std::vector<Float_t> mvaValue_;
  std::vector<Int_t>   mvaCategory_;

  std::vector<Int_t> passMediumId_;

  std::vector<Int_t> isTrue_;

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
PhotonNtuplerVIDwithMVADemo::PhotonNtuplerVIDwithMVADemo(const edm::ParameterSet& iConfig):
  phoMediumIdBoolMapToken_(consumes<edm::ValueMap<bool> >
			   (iConfig.getParameter<edm::InputTag>("phoMediumIdBoolMap"))),
  phoMediumIdFullInfoMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> >
			       (iConfig.getParameter<edm::InputTag>("phoMediumIdFullInfoMap"))),
  mvaValuesMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap"))),
  mvaCategoriesMapToken_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaCategoriesMap")))
{

  //
  // Prepare tokens for all input collections and objects
  //


  // AOD tokens
  photonsToken_    = mayConsume<edm::View<reco::Photon> >
    (iConfig.getParameter<edm::InputTag>
     ("photons"));

  genParticlesToken_ = mayConsume<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticles"));

  // MiniAOD tokens
  // For photons, use the fact that pat::Photon can be cast into 
  // reco::Photon
  photonsMiniAODToken_    = mayConsume<edm::View<reco::Photon> >
    (iConfig.getParameter<edm::InputTag>
     ("photonsMiniAOD"));

  genParticlesMiniAODToken_ = mayConsume<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticlesMiniAOD"));


  edm::Service<TFileService> fs;
  photonTree_ = fs->make<TTree> ("PhotonTree", "Photon data");
  
  photonTree_->Branch("run"        ,  &run_     , "run/I");
  photonTree_->Branch("lumi"       ,  &lumi_    , "lumi/I");
  photonTree_->Branch("evtnum"     ,  &evtnum_  , "evtnum/I");

  photonTree_->Branch("nPho",  &nPhotons_ , "nPho/I");
  photonTree_->Branch("pt"  ,  &pt_    );
  photonTree_->Branch("eta" ,  &eta_ );
  photonTree_->Branch("phi" ,  &phi_ );

  photonTree_->Branch("mvaVal" ,  &mvaValue_ );
  photonTree_->Branch("mvaCat" ,  &mvaCategory_ );
  
  photonTree_->Branch("passMediumId" ,  &passMediumId_ );

  photonTree_->Branch("isTrue"             , &isTrue_);

}


PhotonNtuplerVIDwithMVADemo::~PhotonNtuplerVIDwithMVADemo()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PhotonNtuplerVIDwithMVADemo::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;

  // Save global info right away
  run_ = iEvent.id().run();
  lumi_ = iEvent.id().luminosityBlock();
  evtnum_ = iEvent.id().event();

  // Retrieve the collection of photons from the event.
  // If we fail to retrieve the collection with the standard AOD
  // name, we next look for the one with the stndard miniAOD name.
  //   We use exactly the same handle for AOD and miniAOD formats
  // since pat::Photon objects can be recast as reco::Photon objects.
  edm::Handle<edm::View<reco::Photon> > photons;
  bool isAOD = true;
  iEvent.getByToken(photonsToken_, photons);
  if( !photons.isValid() ){
    isAOD = false;
    iEvent.getByToken(photonsMiniAODToken_,photons);
  }
  
  // Get the MC collection
  Handle<edm::View<reco::GenParticle> > genParticles;
  if( isAOD )
    iEvent.getByToken(genParticlesToken_,genParticles);
  else
    iEvent.getByToken(genParticlesMiniAODToken_,genParticles);

  // Get the photon ID data from the event stream.
  // Note: this implies that the VID ID modules have been run upstream.
  // If you need more info, check with the EGM group.
  //
  // The first map simply has pass/fail for each particle
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  iEvent.getByToken(phoMediumIdBoolMapToken_,medium_id_decisions);
  //
  // The second map has the full info about the cut flow
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > medium_id_cutflow_data;
  iEvent.getByToken(phoMediumIdFullInfoMapToken_,medium_id_cutflow_data);

  // Get MVA values and categories (optional)
  edm::Handle<edm::ValueMap<float> > mvaValues;
  edm::Handle<edm::ValueMap<int> > mvaCategories;
  iEvent.getByToken(mvaValuesMapToken_,mvaValues);
  iEvent.getByToken(mvaCategoriesMapToken_,mvaCategories);

  // Clear vectors
  nPhotons_ = 0;
  pt_.clear();
  eta_.clear();
  phi_.clear();
  //
  mvaValue_.clear();
  mvaCategory_.clear();
  passMediumId_.clear();
  //
  isTrue_.clear();

  // Loop over photons
  for (size_t i = 0; i < photons->size(); ++i){
    const auto pho = photons->ptrAt(i);

    // Kinematics
    if( pho->pt() < 15 ) // keep only photons above 15 GeV
      continue;
    
    nPhotons_++;

    //
    // Save photon kinematics
    //
    pt_  .push_back( pho->pt() );
    eta_ .push_back( pho->superCluster()->eta() );
    phi_ .push_back( pho->superCluster()->phi() );

    //
    // Look up and save the ID decisions
    // 
    bool isPassMedium = (*medium_id_decisions)[pho];
    passMediumId_.push_back( (int)isPassMedium);

    vid::CutFlowResult fullCutFlowData = (*medium_id_cutflow_data)[pho];
    printf("\nDEBUG: name= %s\n", fullCutFlowData.cutFlowName().c_str());
    printf("       standard decision= %d   cutflow decision=%d\n",
	   (int) isPassMedium, (int) fullCutFlowData.cutFlowPassed());
    int ncuts = fullCutFlowData.cutFlowSize();
    for(int icut = 0; icut<ncuts; icut++){
      printf("  %d     name=%s   masked=%d  value=%f    pass=%d\n", icut,
	     fullCutFlowData.getNameAtIndex(icut).c_str(),
	     (int)fullCutFlowData.isCutMasked(icut),
	     fullCutFlowData.getValueCutUpon(icut),
	     (int)fullCutFlowData.getCutResultByIndex(icut));
    }

    mvaValue_.push_back( (*mvaValues)[pho] );
    mvaCategory_.push_back( (*mvaCategories)[pho] );

    // Save MC truth match
    isTrue_.push_back( matchToTruth(*pho, genParticles) );

   }
   
  // Save the info
  photonTree_->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
PhotonNtuplerVIDwithMVADemo::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PhotonNtuplerVIDwithMVADemo::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
PhotonNtuplerVIDwithMVADemo::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
PhotonNtuplerVIDwithMVADemo::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
PhotonNtuplerVIDwithMVADemo::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
PhotonNtuplerVIDwithMVADemo::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PhotonNtuplerVIDwithMVADemo::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

int PhotonNtuplerVIDwithMVADemo::matchToTruth(const reco::Photon &pho, 
				   const edm::Handle<edm::View<reco::GenParticle>>  
				   &genParticles)
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
    else
      return MATCHED_FROM_OTHER_SOURCES;
  }
  return MATCHED_FROM_GUDSCB;
   
}

void PhotonNtuplerVIDwithMVADemo::findFirstNonPhotonMother(const reco::Candidate *particle,
						int &ancestorPID, int &ancestorStatus){
  
  if( particle == 0 ){
    printf("PhotonNtuplerVIDDemo: ERROR! null candidate pointer, this should never happen\n");
    return;
  }

  // Is this the first non-photon parent? If yes, return, otherwise
  // go deeper into recursion
  if( abs(particle->pdgId()) == 22 ){
    findFirstNonPhotonMother(particle->mother(0), ancestorPID, ancestorStatus);
  }else{
    ancestorPID = particle->pdgId();
    ancestorStatus = particle->status();
  }
  
  return;
}

//define this as a plug-in
DEFINE_FWK_MODULE(PhotonNtuplerVIDwithMVADemo);
