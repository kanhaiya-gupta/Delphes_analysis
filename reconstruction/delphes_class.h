//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Apr 17 12:26:30 2021 by ROOT version 5.34/38
// from TTree delphes/a simple tree for analysis
// found on file: signal.root
//////////////////////////////////////////////////////////

#ifndef delphes_class_h
#define delphes_class_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
using std::vector;

// Fixed size dimensions of array or collections stored in the TTree if any.

class delphes_class {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           events;
   Int_t           event_number;
   vector<int>     *PID;
   vector<int>     *mother_1;
   vector<int>     *mother_2;
   vector<int>     *daughter_1;
   vector<int>     *daughter_2;
   vector<int>     *charge;
   Int_t           elec_n;
   vector<int>     *elec_charge;
   vector<float>   *elec_pt;
   vector<float>   *elec_eta;
   vector<float>   *elec_phi;
   vector<float>   *elec_sumpt_charged;
   Int_t           muon_n;
   vector<int>     *muon_charge;
   vector<float>   *muon_pt;
   vector<float>   *muon_eta;
   vector<float>   *muon_phi;
   vector<float>   *muon_sumpt_charged;
   Int_t           phot_n;
   vector<float>   *phot_pt;
   vector<float>   *phot_eta;
   vector<float>   *phot_phi;
   vector<float>   *phot_E;
   Float_t         met;
   Float_t         met_phi;
   Float_t         met_eta;
   Float_t         HT;
   Int_t           jet_n;
   vector<float>   *jet_pt;
   vector<float>   *jet_eta;
   vector<float>   *jet_phi;
   vector<float>   *jet_mass;
   vector<int>     *jet_flavor;
   vector<int>     *jet_btag;
   vector<int>     *jet_tautag;
   vector<float>   *jet_ncharged;
   Int_t           bjet_n;

   // List of branches
   TBranch        *b_events;   //!
   TBranch        *b_event_number;   //!
   TBranch        *b_PID;   //!
   TBranch        *b_mother_1;   //!
   TBranch        *b_mother_2;   //!
   TBranch        *b_daughter_1;   //!
   TBranch        *b_daughter_2;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_elec_n;   //!
   TBranch        *b_elec_charge;   //!
   TBranch        *b_elec_pt;   //!
   TBranch        *b_elec_eta;   //!
   TBranch        *b_elec_phi;   //!
   TBranch        *b_elec_sumpt_charged;   //!
   TBranch        *b_muon_n;   //!
   TBranch        *b_muon_charge;   //!
   TBranch        *b_muon_pt;   //!
   TBranch        *b_muon_eta;   //!
   TBranch        *b_muon_phi;   //!
   TBranch        *b_muon_sumpt_charged;   //!
   TBranch        *b_phot_n;   //!
   TBranch        *b_phot_pt;   //!
   TBranch        *b_phot_eta;   //!
   TBranch        *b_phot_phi;   //!
   TBranch        *b_phot_E;   //!
   TBranch        *b_met;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_met_eta;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_jet_n;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_mass;   //!
   TBranch        *b_jet_flavor;   //!
   TBranch        *b_jet_btag;   //!
   TBranch        *b_jet_tautag;   //!
   TBranch        *b_jet_ncharged;   //!
   TBranch        *b_bjet_n;   //!

   delphes_class(TTree *tree=0);
   virtual ~delphes_class();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef delphes_class_cxx

delphes_class::~delphes_class()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t delphes_class::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t delphes_class::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void delphes_class::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   PID = 0;
   mother_1 = 0;
   mother_2 = 0;
   daughter_1 = 0;
   daughter_2 = 0;
   charge = 0;
   elec_charge = 0;
   elec_pt = 0;
   elec_eta = 0;
   elec_phi = 0;
   elec_sumpt_charged = 0;
   muon_charge = 0;
   muon_pt = 0;
   muon_eta = 0;
   muon_phi = 0;
   muon_sumpt_charged = 0;
   phot_pt = 0;
   phot_eta = 0;
   phot_phi = 0;
   phot_E = 0;
   jet_pt = 0;
   jet_eta = 0;
   jet_phi = 0;
   jet_mass = 0;
   jet_flavor = 0;
   jet_btag = 0;
   jet_tautag = 0;
   jet_ncharged = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("events", &events, &b_events);
   fChain->SetBranchAddress("event_number", &event_number, &b_event_number);
   fChain->SetBranchAddress("PID", &PID, &b_PID);
   fChain->SetBranchAddress("mother_1", &mother_1, &b_mother_1);
   fChain->SetBranchAddress("mother_2", &mother_2, &b_mother_2);
   fChain->SetBranchAddress("daughter_1", &daughter_1, &b_daughter_1);
   fChain->SetBranchAddress("daughter_2", &daughter_2, &b_daughter_2);
   fChain->SetBranchAddress("charge", &charge, &b_charge);
   fChain->SetBranchAddress("elec_n", &elec_n, &b_elec_n);
   fChain->SetBranchAddress("elec_charge", &elec_charge, &b_elec_charge);
   fChain->SetBranchAddress("elec_pt", &elec_pt, &b_elec_pt);
   fChain->SetBranchAddress("elec_eta", &elec_eta, &b_elec_eta);
   fChain->SetBranchAddress("elec_phi", &elec_phi, &b_elec_phi);
   fChain->SetBranchAddress("elec_sumpt_charged", &elec_sumpt_charged, &b_elec_sumpt_charged);
   fChain->SetBranchAddress("muon_n", &muon_n, &b_muon_n);
   fChain->SetBranchAddress("muon_charge", &muon_charge, &b_muon_charge);
   fChain->SetBranchAddress("muon_pt", &muon_pt, &b_muon_pt);
   fChain->SetBranchAddress("muon_eta", &muon_eta, &b_muon_eta);
   fChain->SetBranchAddress("muon_phi", &muon_phi, &b_muon_phi);
   fChain->SetBranchAddress("muon_sumpt_charged", &muon_sumpt_charged, &b_muon_sumpt_charged);
   fChain->SetBranchAddress("phot_n", &phot_n, &b_phot_n);
   fChain->SetBranchAddress("phot_pt", &phot_pt, &b_phot_pt);
   fChain->SetBranchAddress("phot_eta", &phot_eta, &b_phot_eta);
   fChain->SetBranchAddress("phot_phi", &phot_phi, &b_phot_phi);
   fChain->SetBranchAddress("phot_E", &phot_E, &b_phot_E);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("met_eta", &met_eta, &b_met_eta);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("jet_n", &jet_n, &b_jet_n);
   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_mass", &jet_mass, &b_jet_mass);
   fChain->SetBranchAddress("jet_flavor", &jet_flavor, &b_jet_flavor);
   fChain->SetBranchAddress("jet_btag", &jet_btag, &b_jet_btag);
   fChain->SetBranchAddress("jet_tautag", &jet_tautag, &b_jet_tautag);
   fChain->SetBranchAddress("jet_ncharged", &jet_ncharged, &b_jet_ncharged);
   fChain->SetBranchAddress("bjet_n", &bjet_n, &b_bjet_n);
   Notify();
}

Bool_t delphes_class::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void delphes_class::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t delphes_class::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef delphes_class_cxx
