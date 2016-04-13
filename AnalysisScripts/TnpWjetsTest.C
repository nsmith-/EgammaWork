#include "TLorentzVector.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TCut.h"

#include <vector>
#include <cmath>

const int nWP = 5;
const TString wpName[nWP] = 
  {"Veto", "Loose", "Medium", "Tight", "HEEP"};

const TString tagIDBranchName = "passTightId";
const TString probeIDBranchName = "passMediumId";

const TString treename = "ntupler/ElectronTree";
// const TString fname = "/afs/cern.ch/user/i/ikrav/workspace/ntuples/Spring15/WJetsToLNu_tnp_study.root";
const TString fname = "/afs/cern.ch/user/i/ikrav/workspace/ntuples/data/SingleElectron_data_2015D_tnp_study_partial.root";
//const TString fname = "../ElectronNtupler/test/electron_ntuple_mini_data100K.root";

bool verbose = false;
bool smallEventCount = false;

const float ptTagMin = 20;
const float ptProbeMin = 10;


void TnpWjetsTest(){

  // This statement below should not be needed, but in one particular node I had to
  // add it, somehow the vector header was not loaded automatically there.
  gROOT->ProcessLine("#include <vector>"); 

  // =========================================================================
  //    Configure inputs
  // =========================================================================

  //
  // Find the tree
  //
  TFile *file1 = new TFile(fname);
  if( !file1 )
    assert(0);
  TTree *tree = (TTree*)file1->Get(treename);
  if( !tree )
    assert(0);

  // Event-level variables:
  int nEle;
  // Per-electron variables
  // Kinematics
  std::vector <float> *pt = 0;     
  std::vector <float> *eta = 0;    
  std::vector <float> *phi = 0;    
  std::vector <int> *passTagId = 0;    
  std::vector <int> *passProbeId = 0;    
  std::vector <int> *charge = 0;    
  std::vector <int> *isTrue = 0;    

  // Declare branches
  TBranch *b_nEle = 0;
  TBranch *b_pt = 0;
  TBranch *b_eta = 0;
  TBranch *b_phi = 0;
  TBranch *b_passTagId = 0;
  TBranch *b_passProbeId = 0;
  TBranch *b_charge = 0;
  TBranch *b_isTrue = 0;

  // Connect variables and branches to the tree with the data
  tree->SetBranchAddress("nEle", &nEle, &b_nEle);
  tree->SetBranchAddress("pt", &pt, &b_pt);
  tree->SetBranchAddress("eta", &eta, &b_eta);  
  tree->SetBranchAddress("phi", &phi, &b_phi);  

  tree->SetBranchAddress(tagIDBranchName, &passTagId, &b_passTagId);
  tree->SetBranchAddress(probeIDBranchName, &passProbeId, &b_passProbeId);

  tree->SetBranchAddress("q", &charge, &b_charge);

  tree->SetBranchAddress("isTrue", &isTrue, &b_isTrue);

  if(verbose)
    printf("Configured inputs\n");
  // =========================================================================
  //    Configure outputs
  // =========================================================================

  TH1F *hmass = new TH1F("hmass","",120, 0, 120);
  hmass->SetDirectory(0);
  TH1F *hptProbe = new TH1F("hptProbe","",100,0,100);
  hptProbe->SetDirectory(0);

  TFile fout("tnpCandidates.root","recreate");
  TTree tout("tnpCandidates","the tree with tag and probe candidates");

  Int_t event;
  Float_t mass;
  Float_t tagPt, tagPhi, tagEta, tagQ;
  Int_t tagIsTrue, tagPassId;
  Float_t probePt, probePhi, probeEta, probeQ;
  Int_t probeIsTrue, probePassId;

  tout.Branch("event", &event, "event/I");

  tout.Branch("mass",&mass,"mass/F");

  tout.Branch("tagPt", &tagPt, "tagPt/F");
  tout.Branch("tagPhi", &tagPhi, "tagPhi/F");
  tout.Branch("tagEta", &tagEta, "tagEta/F");
  tout.Branch("tagQ", &tagQ, "tagQ/F");
  tout.Branch("tagIsTrue", &tagIsTrue, "tagIsTrue/I");
  tout.Branch("tagPassId", &tagPassId, "tagPassId/I");

  tout.Branch("probePt", &probePt, "probePt/F");
  tout.Branch("probePhi", &probePhi, "probePhi/F");
  tout.Branch("probeEta", &probeEta, "probeEta/F");
  tout.Branch("probeQ", &probeQ, "probeQ/F");
  tout.Branch("probeIsTrue", &probeIsTrue, "probeIsTrue/I");
  tout.Branch("probePassId", &probePassId, "probePassId/I");

  if(verbose)
    printf("Configured outputs\n");
  // =========================================================================
  //    Run over events
  // =========================================================================

  // 
  // Loop over events
  //
  UInt_t maxEvents = tree->GetEntries();
  if( smallEventCount )
    maxEvents = std::min((float)10000, (float)maxEvents);
  if(verbose)
    printf("Start loop over events, total events = %lld\n", 
           tree->GetEntries() );
  for(UInt_t ievent = 0; ievent < maxEvents; ievent++){

    if( ievent%100000 == 0){
      printf("."); fflush(stdout);
    }
    Long64_t tentry = tree->LoadTree(ievent);
    event = ievent;

    // Load the value of the number of the electrons in the event    
    b_nEle->GetEntry(tentry);
    if(verbose)
      printf("Event %d, number of electrons %u\n", ievent, nEle);

    // Get data for all electrons in this event, only vars of interest
    b_pt->GetEntry(tentry);
    b_eta->GetEntry(tentry);
    b_phi->GetEntry(tentry);
    b_charge->GetEntry(tentry);

    b_passTagId->GetEntry(tentry);
    b_passProbeId->GetEntry(tentry);
    b_isTrue->GetEntry(tentry);

    // Loop over dielectrons electrons. Use both possible assignments
    // for the tag and the probe
    for(int iele1 = 0; iele1 < nEle; iele1++){

      // This is the tag
      tagPt = pt->at(iele1);
      tagEta = eta->at(iele1);
      tagPhi = phi->at(iele1);
      tagQ = charge->at(iele1);
      tagIsTrue = isTrue->at(iele1);
      tagPassId = passTagId->at(iele1);
      
      if( ! (tagPt > ptTagMin ) )
	continue;
      if( ! ( fabs(tagEta)<2.5 ) )
	      // && ( fabs(tagEta)<1.4442 || fabs(tagEta)<1.566) ) )
	continue;

      // if( ! (tagIsTrue ==1) )
      // 	continue;

      for(int iele2 = 0; iele2 < nEle; iele2++){
	
	if( iele1 == iele2 )
	  continue;

	probePt = pt->at(iele2);
	probeEta = eta->at(iele2);
	probePhi = phi->at(iele2);
	probeQ = charge->at(iele2);
	probeIsTrue = isTrue->at(iele2);
	probePassId = passProbeId->at(iele2);

	// This is the probe
	if( ! ( probePt > ptProbeMin ) )
	  continue;
	if( ! ( fabs(probeEta)<2.5 ) )
		// && ( fabs( probeEta)<1.4442 || fabs(probeEta)<1.566) ) )
	  continue;

	hptProbe->Fill( probePt );
	
	TLorentzVector el1, el2;
	el1.SetPtEtaPhiM( tagPt, tagEta, tagPhi, 0.000511);
	el2.SetPtEtaPhiM( probePt, probeEta, probePhi, 0.000511);
	TLorentzVector dielectron = el1+el2;
	mass = dielectron.M();
	hmass->Fill(mass);	  

	tout.Fill();
	
      } // end loop over probes
    } // end loop over tags
    
  } // end loop over events
  
  TCanvas *c1 = new TCanvas("c1","c1",10,10,800,800);
  c1->cd();
  hptProbe->Draw();
  c1->Update();

  TCanvas *c2 = new TCanvas("c2","c2",100,10,800,800);
  c2->cd();
  hmass->Draw();
  c2->Update();

  tout.Write();
  fout.Close();
}
