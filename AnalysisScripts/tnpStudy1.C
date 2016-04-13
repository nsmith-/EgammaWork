#include "TTree.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH2F.h"
#include "TMath.h"

const TString fname = "tnpCandidates.root";
const TString treename = "tnpCandidates";

// const TString cutsTag   = "tagPassId==1 && tagPt>30";
const TString cutsTag   = "tagPassId==0 && tagPt>30";
const TString cutsEvent = "tagQ*probeQ>0";

const int nPtBins = 5;
const float ptBinLimits[nPtBins+1] = {
  10, 20, 30, 40, 50, 200
};

const int nEtaBins = 2;
const float etaBinLimits[nEtaBins+1] = {
  0.0, 1.479, 2.5
};

// const int nEtaBins = 4;
// const float etaBinLimits[nEtaBins+1] = {
//   0.0, 0.8, 1.479, 2.0, 2.5
// };

// Forward declarations
Double_t cmsShape( Double_t *x, Double_t *par);

// Main
void tnpStudy1(){

  //
  // Find the tree
  //
  TFile *file = new TFile(fname);
  if( !file )
    assert(0);
  TTree *tree = (TTree*)file->Get(treename);
  if( !tree )
    assert(0);

  //
  // Prepare plots and fit them
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,800,800);
  c1->Divide(nPtBins, nEtaBins);

  TF1 *func = new TF1("func",cmsShape, 40, 200, 5);
  // Parameters:
  //     norm - normalization
  //     alpha    - the position of the turn-on
  //     beta      - rate of the turn-on
  //     gamma - rate of the exponential tail decrease
  //     peak     -  controls normalization only and doesn’t affect the shape, fix it to a const
  func->SetParLimits(1, 30, 500);
  func->FixParameter(4, 72);

  TH1F *hist[nPtBins][nEtaBins];
  
  float resA[nPtBins][nEtaBins];
  float resB[nPtBins][nEtaBins];
  float resG[nPtBins][nEtaBins];
  float resAErr[nPtBins][nEtaBins];
  float resBErr[nPtBins][nEtaBins];
  float resGErr[nPtBins][nEtaBins];
  for(int ieta = 0; ieta<nEtaBins; ieta++){
    for(int ipt = 0; ipt<nPtBins; ipt++){
      
      c1->cd(ipt + 1 + nPtBins*ieta);
      TString hname = TString::Format("hist_pt%d_eta%d",ipt, ieta);
      hist[ipt][ieta] = new TH1F(hname, hname, 100, 0, 200);
      TString drawCommand = TString::Format("mass>>%s", hname.Data());

      TString cutsProbe 
	= TString::Format(" probePt>=%f && probePt<%f && abs(probeEta)>=%f && abs(probeEta)<%f && probePassId==0",
			  ptBinLimits[ipt], ptBinLimits[ipt+1], 
			  etaBinLimits[ieta], etaBinLimits[ieta+1]);
      TString cuts = cutsTag + TString(" && ") + cutsProbe + TString(" && ") + cutsEvent;
      tree->Draw(drawCommand,cuts);
      printf("The cut string is %s\n", cuts.Data());
      
      func->SetParameters(100, 60, 0.08, 0.022, 72);

      hist[ipt][ieta]->Fit("func","R");

      resA[ipt][ieta] = func->GetParameter(1);
      resB[ipt][ieta] = func->GetParameter(2);
      resG[ipt][ieta] = func->GetParameter(3);
      resAErr[ipt][ieta] = func->GetParError(1);
      resBErr[ipt][ieta] = func->GetParError(2);
      resGErr[ipt][ieta] = func->GetParError(3);
    }
  }

  // Draw fit parameters
  TGraphErrors *graphA[nEtaBins];
  TGraphErrors *graphB[nEtaBins];
  TGraphErrors *graphG[nEtaBins];
  float *arrayA;
  float *arrayB;
  float *arrayG;
  float *arrayAErr;
  float *arrayBErr;
  float *arrayGErr;
  float x[nPtBins];
  float dx[nPtBins];
  for(int ieta=0; ieta<nEtaBins; ieta++){
    arrayA = new float[nPtBins];
    arrayB = new float[nPtBins];
    arrayG = new float[nPtBins];
    arrayAErr = new float[nPtBins];
    arrayBErr = new float[nPtBins];
    arrayGErr = new float[nPtBins];
    for(int ipt=0; ipt<nPtBins; ipt++){
      arrayA[ipt] = resA[ipt][ieta];
      arrayB[ipt] = resB[ipt][ieta];
      arrayG[ipt] = resG[ipt][ieta];
      arrayAErr[ipt] = resAErr[ipt][ieta];
      arrayBErr[ipt] = resBErr[ipt][ieta];
      arrayGErr[ipt] = resGErr[ipt][ieta];
      x[ipt] = (ptBinLimits[ipt] + ptBinLimits[ipt+1])/2.0;
      dx[ipt] = (ptBinLimits[ipt+1] - ptBinLimits[ipt])/2.0;
    } // end loop over pt bins
    graphA[ieta] = new TGraphErrors(nPtBins, x, arrayA, dx, arrayAErr);
    graphB[ieta] = new TGraphErrors(nPtBins, x, arrayB, dx, arrayBErr);
    graphG[ieta] = new TGraphErrors(nPtBins, x, arrayG, dx, arrayGErr);
  } // end loop over eta bins

  TCanvas *c2 = new TCanvas("c2","c2",200,10,600, 600);
  c2->cd();

  TH2F *dummyA = new TH2F("dummyA", "", 100, 0, 200, 100, 0, 250);
  dummyA->GetXaxis()->SetTitle("probe pt");
  dummyA->GetYaxis()->SetTitle("alpha");
  dummyA->Draw();
  for(int ieta=0; ieta<nEtaBins; ieta++){
    graphA[ieta]->SetMarkerSize(1);
    graphA[ieta]->SetMarkerStyle(20+ieta);
    graphA[ieta]->Draw("PE,same");
  }

  TCanvas *c3 = new TCanvas("c3","c3",300,10,600, 600);
  c3->cd(3);
  TH2F *dummyB = new TH2F("dummyB", "", 100, 0, 200, 100, 0, 0.3);
  dummyB->GetXaxis()->SetTitle("probe pt");
  dummyB->GetYaxis()->SetTitle("beta");
  dummyB->Draw();
  for(int ieta=0; ieta<nEtaBins; ieta++){
    graphB[ieta]->SetMarkerSize(1);
    graphB[ieta]->SetMarkerStyle(20+ieta);
    graphB[ieta]->Draw("PE,same");
  }

  TCanvas *c4 = new TCanvas("c4","c4",400,10,600, 600);
  c4->cd();
  TH2F *dummyG = new TH2F("dummyG", "", 100, 0, 200, 100, 0, 0.1);
  dummyG->GetXaxis()->SetTitle("probe pt");
  dummyG->GetYaxis()->SetTitle("gamma");
  dummyG->Draw();
  for(int ieta=0; ieta<nEtaBins; ieta++){
    graphG[ieta]->SetMarkerSize(1);
    graphG[ieta]->SetMarkerStyle(20+ieta);
    graphG[ieta]->Draw("PE,same");
  }

}



Double_t cmsShape( Double_t *x, Double_t *par)
 { 

   Double_t xx = x[0];
   Double_t norm = par[0];
   Double_t alpha = par[1];
   Double_t beta = par[2];
   Double_t gamma = par[3];
   Double_t peak = par[4];
   
   Double_t erf = TMath::Erfc((alpha - xx) * beta);
   Double_t u = (xx - peak)*gamma;
   
   if(u < -70) u = 1e20;
   else if( u>70 ) u = 0;
   else u = exp(-u);   //exponential decay
   return norm*erf*u;
 } 
