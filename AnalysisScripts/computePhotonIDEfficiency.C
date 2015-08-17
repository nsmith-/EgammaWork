#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TCut.h"

#include <vector>

const int nWP = 3;
enum WpType { WP_LOOSE = 0,
              WP_MEDIUM, 
              WP_TIGHT};
const TString wpName[nWP] = 
  {"Loose", "Medium", "Tight"};

const TString treename = "ntupler/PhotonTree";
const TString fname = "../PhotonNtupler/test/photon_ntuple_mini.root";

bool verbose = false;
bool smallEventCount = false;
bool useWeights = true;

const float ptmin = 30;
const float ptmax = 200;

void computePhotonIDEfficiency(WpType wp)
{

  // This statement below should not be needed, but in one particular node I had to
  // add it, somehow the vector header was not loaded automatically there.
  gROOT->ProcessLine("#include <vector>"); 

  if(wp<0 || wp>nWP) assert(0);
  printf("Compute the efficiencies for the working point %s\n", wpName[wp].Data());

  //
  // Find the tree
  //
  TFile *file1 = new TFile(fname);
  if( !file1 )
    assert(0);
  TTree *tree = (TTree*)file1->Get(treename);
  if( !tree )
    assert(0);

  // Weight histograms

  TH2D *hSignal;
  TH2D *hBackground;
  if( tree->GetEntries() > 1000000 ){
    hSignal = new TH2D("hSignal","",200,-5,5,185,ptmin,ptmax);
    hBackground = new TH2D("hBackground","",200,-5,5,185,ptmin,ptmax);
  }else{
    // If we are running on a small sample, use courser binning for weights
    // histograms
    hSignal = new TH2D("hSignal","",20,-5,5,18,ptmin,ptmax);
    hBackground = new TH2D("hBackground","",20,-5,5,18,ptmin,ptmax);
  }
  TCut miscCut = "1"; // no other cuts, no electron veto like "hasPixelSeed==0";
  TCut sigCut  = "isTrue==1";
  TCut bgCut  = "isTrue!=1";
  if( smallEventCount ){
    tree->Draw("pt:eta>>hSignal",miscCut && sigCut, "colz", 100000);
    tree->Draw("pt:eta>>hBackground",miscCut && bgCut, "colz", 100000);
  }else{
    tree->Draw("pt:eta>>hSignal",miscCut && sigCut, "colz");
    tree->Draw("pt:eta>>hBackground",miscCut && bgCut, "colz");
  }

  TH1D *ptCheckHist = new TH1D("ptCheckHist","",185, ptmin, ptmax);
  TH1D *etaCheckHist = new TH1D("etaCheckHist","",200, -5, 5);

  // Event-level variables:
  int nPho;
  // Per-photon variables
  // Kinematics
  std::vector <float> *pt = 0;     
  std::vector <float> *eta = 0;    
  std::vector <int> *passChosenId = 0;    
  std::vector <int> *isTrue = 0;    

  // Declare branches
  TBranch *b_nPho = 0;
  TBranch *b_pt = 0;
  TBranch *b_eta = 0;
  TBranch *b_passChosenId = 0;
  TBranch *b_isTrue = 0;

  // Connect variables and branches to the tree with the data
  tree->SetBranchAddress("nPho", &nPho, &b_nPho);
  tree->SetBranchAddress("pt", &pt, &b_pt);
  tree->SetBranchAddress("eta", &eta, &b_eta);

  if( wp == WP_LOOSE )
    tree->SetBranchAddress("passLooseId", &passChosenId, &b_passChosenId);
  else if( wp == WP_MEDIUM)
    tree->SetBranchAddress("passMediumId", &passChosenId, &b_passChosenId);
  else if( wp == WP_TIGHT)
    tree->SetBranchAddress("passTightId", &passChosenId, &b_passChosenId);
  else
    assert(0);
  tree->SetBranchAddress("isTrue", &isTrue, &b_isTrue);

  // Declare event counters for efficiency calculations
  double sumSignalDenomEB = 0;
  double sumSignalNumEB   = 0;
  double sumSignalDenomEE = 0;
  double sumSignalNumEE   = 0;
  double sumBackDenomEB = 0;
  double sumBackNumEB   = 0;
  double sumBackDenomEE = 0;
  double sumBackNumEE   = 0;

  double sumSignalDenomEBErr2 = 0;
  double sumSignalNumEBErr2   = 0;
  double sumSignalDenomEEErr2 = 0;
  double sumSignalNumEEErr2   = 0;
  double sumBackDenomEBErr2 = 0;
  double sumBackNumEBErr2   = 0;
  double sumBackDenomEEErr2 = 0;
  double sumBackNumEEErr2   = 0;

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

    // Load the value of the number of the photons in the event    
    b_nPho->GetEntry(tentry);
    if(verbose)
      printf("Event %d, number of photons %u\n", ievent, nPho);

    // Get data for all photons in this event, only vars of interest
    b_pt->GetEntry(tentry);
    b_eta->GetEntry(tentry);
    b_passChosenId->GetEntry(tentry);
    b_isTrue->GetEntry(tentry);

    // Loop over photons
    for(int ipho = 0; ipho < nPho; ipho++){
      
      // Preselection
      if( !(pt->at(ipho) > ptmin && pt->at(ipho) < ptmax ) ) continue;
      if( fabs(eta->at(ipho))>1.4442 && fabs(eta->at(ipho))<1.566) continue;
      if( fabs(eta->at(ipho))>2.5 ) continue;
      
      bool isBarrel = (fabs(eta->at(ipho)) < 1.479);
      
      double weight = 1.0;
      if( useWeights ){
	// We reweight signal and background (separately) to have 
	// a flat pt and eta distribution. This step is a matter of definition.
	double binContent = 0;
	TH2D *hEvents = hSignal;
	if( isTrue->at(ipho) != 1 )
	  hEvents = hBackground;
	binContent = hEvents->GetBinContent
	  ( hEvents->FindBin( eta->at(ipho), pt->at(ipho) ) );
	if( binContent == 0 ){
	  printf("Error! Zero! pt=%f, eta=%f\n", pt->at(ipho), eta->at(ipho));
	}else{
	  weight = 1./binContent;
	}
      }

      if(verbose)
	printf(" next photon pt=%f id decision is %d\n",
	       pt->at(ipho),
	       passChosenId->at(ipho));

      // Fill signal counters
      if( isTrue->at(ipho) == 1 ) {
        ptCheckHist->Fill(pt->at(ipho), weight);
        etaCheckHist->Fill(eta->at(ipho), weight);

	if( isBarrel ) {
	  sumSignalDenomEB += weight;
	  sumSignalDenomEBErr2 += weight*weight;
	  if( passChosenId->at(ipho) ) {
	    sumSignalNumEB += weight;
	    sumSignalNumEBErr2 += weight*weight;
	  }
	} else {
	  sumSignalDenomEE += weight;
	  sumSignalDenomEEErr2 += weight*weight;          
	  if( passChosenId->at(ipho) ) {
	    sumSignalNumEE += weight;
	    sumSignalNumEEErr2 += weight*weight;
	  }
	  
	}// end barrel / endcap
      } // end if signal
      
      // Fill background counters
      if( isTrue->at(ipho) !=1 ) {
	
        if( isBarrel ) {
	  
          sumBackDenomEB += weight;
          sumBackDenomEBErr2 += weight*weight;
          if( passChosenId->at(ipho) ) {
            sumBackNumEB += weight;
            sumBackNumEBErr2 += weight*weight;
          }
        } else {
          sumBackDenomEE += weight;
          sumBackDenomEEErr2 += weight*weight;    
          if( passChosenId->at(ipho) ) {
            sumBackNumEE += weight;
            sumBackNumEEErr2 += weight*weight;
          }
	  
        }// end barrel / endcap
      } // end if background
      
    } // end loop over photons
    
    
  }// end loop over events  

  // DEBUG
  if(verbose){
    printf("barrel signal pass=%.1f   total=%.1f\n", sumSignalNumEB, sumSignalDenomEB);
    printf("endcap signal pass=%.1f   total=%.1f\n", sumSignalNumEE, sumSignalDenomEE);
  }

  printf("\nEfficiencies for the working point %s\n", wpName[wp].Data());

  // Compute signal efficiencies
  double effSignalEB = sumSignalNumEB / sumSignalDenomEB;
  double effSignalEBErr = sqrt( sumSignalDenomEBErr2 
                                * effSignalEB*(1-effSignalEB)
                                /(sumSignalDenomEB*sumSignalDenomEB) );
  printf("Signal barrel efficiency: %5.1f +- %5.1f %%\n", 
         effSignalEB*100, effSignalEBErr*100 );

  double effSignalEE = sumSignalNumEE / sumSignalDenomEE;
  double effSignalEEErr = sqrt( sumSignalDenomEEErr2 
                                * effSignalEE*(1-effSignalEE)
                                /(sumSignalDenomEE*sumSignalDenomEE) );
  printf("Signal endcap efficiency: %5.1f +- %5.1f %%\n", 
         effSignalEE*100, effSignalEEErr*100 );

  // Compute background efficiencies
  double effBackEB = sumBackNumEB / sumBackDenomEB;
  double effBackEBErr = sqrt( sumBackDenomEBErr2 
                                * effBackEB*(1-effBackEB)
                                /(sumBackDenomEB*sumBackDenomEB) );
  printf("Background barrel efficiency: %5.1f +- %5.1f %%\n", 
         effBackEB*100, effBackEBErr*100 );

  double effBackEE = sumBackNumEE / sumBackDenomEE;
  double effBackEEErr = sqrt( sumBackDenomEEErr2 
                                * effBackEE*(1-effBackEE)
                                /(sumBackDenomEE*sumBackDenomEE) );
  printf("Background endcap efficiency: %5.1f +- %5.1f %%\n", 
         effBackEE*100, effBackEEErr*100 );

  // Plot a few histograms as a check of the reweighting
  TCanvas *c2 = new TCanvas("c2","c2",10,10,800, 600);
  c2->Divide(2,1);
  c2->cd(1);
  ptCheckHist->Draw();
  c2->cd(2);
  etaCheckHist->Draw();
  
}
  
