#include <iostream>
#include "TFile.h"
#include "TTree.h"

using namespace std;

const bool printvars =false;

void vidValidationFast(TString fileName, TString treeName, TString idName){

  TFile *fin = new TFile(fileName);
  if( !fin ){
    cout << "Failed to open file " << fileName.Data() << endl;
    return;
  }
  TString fullTreeName = TString::Format("ntupler/%s",treeName.Data());
  TTree *tree = (TTree*)fin->Get(fullTreeName);
  if( !tree ){
    cout << "Failed to find TTree " << fullTreeName.Data() << endl;
    return;
  }

  if( printvars){
    cout << "\nFor reference, here is the content of the tree:" << endl;
    tree->Show(0);
  }

  TString idCut = TString::Format("%s==1",idName.Data());

  int nEleTot = tree->Draw("pt","","goff");
  int nElePass = tree->Draw("pt",idCut,"goff");
  
  cout << "\nQuick check for the ID " << idName.Data() << ":" << endl; 
  cout << " total electrons : " << nEleTot << endl;
  cout << " passed electrons: " << nElePass << endl;

  return;
}
