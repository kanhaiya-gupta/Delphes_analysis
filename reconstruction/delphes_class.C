#define delphes_class_cxx
#include "delphes_class.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include "TLorentzVector.h"
#include "TChain.h"
#include "TMath.h"
#include <vector>
using std::vector;

void delphes_class::Loop()
{
//   In a ROOT session, you can do:
//
//      Root > .L delphes_class.C
//      Root > delphes_class t
//
//
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

  
// creating the root file as output

  TFile outf("sig_tbh.root","RECREATE");
  std::cout << "creating output ..." <<std::endl;

// creating a tree named delphes
 
  TTree Delphes("Delphes","a simple tree for analysis");

  Int_t elec_n,  muon_n, jet_n;
  Float_t elec_pt;

  Delphes.Branch("elec_n", &elec_n , "elec_n/I");
  Delphes.Branch("muon_n", &muon_n , "muon_n/I");
  Delphes.Branch("jet_n", &jet_n , "jet_n/I");
  Delphes.Branch("elec_pt", &elec_pt, "elec_pt/F");
  


   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
 
     // Begin your code here
     elec_pt->clear();
     jet_n->clear();
     





     
    // End the code here

  Delphes.Fill();
     
   }

    // Print the total number of countd events here
   //  cout << "Analyzed a total of               : " << nEvent << " events" << endl;
    Delphes.Write("",TObject::kWriteDelete);
  //  outf.Write();
}

delphes_class::delphes_class(TTree *tree) : fChain(0) 
 {
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  if (tree == 0) {
 
      TChain* tchain = new TChain("delphes");
          tchain->Add("/home/kanhaiya/Delphes_Analysis/signal.root");
          tree = tchain;
   }
   Init(tree);
   Loop();
}

