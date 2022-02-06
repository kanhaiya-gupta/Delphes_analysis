#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "/home/kanhaiya/Madgraph/MG5_aMC_v2_5_5/Delphes/classes/DelphesClasses.h"
#include "/home/kanhaiya/Madgraph/MG5_aMC_v2_5_5/Delphes/external/ExRootAnalysis/ExRootTreeReader.h"
#endif

#include <iostream>
#include <TClonesArray.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
using std::vector;
#include <iostream>
#include <fstream>
#include "TLorentzVector.h"



//------------------------------------------------------------------------------

void test_charged_higgs()
{
  gSystem->Load("/home/kanhaiya/Madgraph/MG5_aMC_v2_5_5/Delphes/libDelphes");

   
  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add("/home/kanhaiya/Madgraph/MG5_aMC_v2_5_5/charged_higgs_signal/Events/run_01/tag_1_delphes_events.root");
 // chain.Add("/home/kanhaiya/Madgraph/MG5_aMC_v2_5_5/ttbar_jets/Events/run_01/tag_1_delphes_events.root");

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t nEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
  TClonesArray *branchScalarHT = treeReader->UseBranch("ScalarHT");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchTower = treeReader->UseBranch("Tower");
//  TClonesArray *branchVertex = treeReader->UseBranch("Vertex");
//  TClonesArray *branchHectorHit = treeReader->UseBranch("HectorHit");
  

// creating the root file as output

  TFile outf("test_sig.root","RECREATE");
  std::cout << "creating output ..." <<std::endl;

// creating a tree named delphes
 
  TTree delphes("delphes","a simple tree for analysis");
  


// Initialize the parameters for the branches
  Event *event;
  GenParticle *particle;
  Electron *electron;
  Muon *muon;
  Photon *photon;
  MissingET *MeT;
  ScalarHT *Scalar_HT;
  Jet *jet;
  TObject *object;


  Int_t events, event_number, lep_n, phot_n, jet_n, bjet_n;
  Float_t weight, met, met_eta, met_phi, HT;
  Double_t random;
  std::vector<int> *PID, *status, *charge, *mother_1, *mother_2, *daughter_1, *daughter_2; 
  std::vector<int> *elec_charge, *muon_charge;
  std::vector<float> *lep_pt, *lep_eta, *lep_phi, *lep_sumpt_charged, *lep_ID;
  std::vector<float> *phot_pt, *phot_eta, *phot_phi, *phot_E;
  std::vector<float> *jet_pt, *jet_eta, *jet_phi, *jet_mass, *jet_ncharged;
  std::vector<int> *jet_flavor, *jet_btag, *jet_tautag;
  



// creating the branches
  
  delphes.Branch("events", &events , "events/I");
  delphes.Branch("event_number", &event_number, "event_number/I");
  delphes.Branch("PID", &PID);
//  delphes.Branch("status", &status);
  delphes.Branch("mother_1", &mother_1);
  delphes.Branch("mother_2", &mother_2);
  delphes.Branch("daughter_1", &daughter_1);
  delphes.Branch("daughter_2", &daughter_2);
  delphes.Branch("charge", &charge);
  delphes.Branch("lep_ID", &lep_ID);
  delphes.Branch("lep_pt", &lep_pt);
  delphes.Branch("lep_eta", &lep_eta);
  delphes.Branch("lep_phi", &lep_phi);
  delphes.Branch("lep_sumpt_charged", &lep_sumpt_charged);
  delphes.Branch("phot_n", &phot_n , "phot_n/I");
  delphes.Branch("phot_pt", &phot_pt);
  delphes.Branch("phot_eta", &phot_eta);
  delphes.Branch("phot_phi", &phot_phi);
  delphes.Branch("phot_E", &phot_E);
  delphes.Branch("met", &met , "met/F");
  delphes.Branch("met_phi", &met_phi , "met_phi/F");
  delphes.Branch("met_eta", &met_eta , "met_eta/F");
  delphes.Branch("HT", &HT , "HT/F");
  delphes.Branch("jet_n", &jet_n , "jet_n/I");
  delphes.Branch("jet_pt", &jet_pt);
  delphes.Branch("jet_eta", &jet_eta);
  delphes.Branch("jet_phi", &jet_phi);
  delphes.Branch("jet_mass", &jet_mass);
  delphes.Branch("jet_flavor", &jet_flavor);
  delphes.Branch("jet_btag", &jet_btag);
  delphes.Branch("jet_tautag", &jet_tautag);
  delphes.Branch("jet_ncharged", &jet_ncharged);
  delphes.Branch("bjet_n", &bjet_n , "bjet_n/I");
  
  
  
  
  // Loop over all events
  for(Int_t entry = 0; entry < nEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);


  events = 0;
   for(Int_t i = 0; i < branchEvent->GetEntries(); ++i)
    {
      
      event = (Event*) branchEvent->At(i);
  
      events++;
      event_number = event->Number;
   //   weight = event->Weight;
    
   }
    
    PID->clear();
  //  status->clear();
    mother_1->clear();
    mother_2->clear();
    daughter_1->clear();
    daughter_2->clear();
    charge->clear();

   for(Int_t i=0;i<branchParticle->GetEntriesFast(); i++){

   
       particle = (GenParticle*) branchParticle->At(i);
    
     //  PID->push_back(particle->PID);
     //  status->push_back(partice->Status);
     //  mother_1->push_back(particle->M1);
     //  mother_2->push_back(particle->M2);
     //  daughter_1->push_back(particle->D1);
     //  daughter_2->push_back(particle->D2);
     //  charge->push_back(particle->Charge); 

}

  
    lep_n = 0;
    lep_pt->clear();
    lep_eta->clear();
    lep_phi->clear();
    lep_ID->clear();
    lep_sumpt_charged->clear();

  // Loop over all electrons in event
    for(Int_t i = 0; i < branchElectron->GetEntries(); ++i)
    {
      electron = (Electron*) branchElectron->At(i);
   
      lep_n++;
      lep_pt->push_back(electron->PT);
      lep_eta->push_back(electron->Eta);
      lep_phi->push_back(electron->Phi);
      lep_ID->push_back((electron->Charge)*-11);
      lep_sumpt_charged->push_back(electron->SumPtCharged);
   
   }
    
     lep_n = 0;
     lep_pt->clear();
     lep_eta->clear();
     lep_phi->clear();
     lep_ID->clear();
     lep_sumpt_charged->clear();

   // Loop over all electrons in event
   for(Int_t i = 0; i < branchMuon->GetEntries(); ++i)
    {
      muon = (Muon*) branchMuon->At(i);
      lep_n++;
      lep_pt->push_back(muon->PT);
      lep_eta->push_back(muon->Eta);
      lep_phi->push_back(muon->Phi);
      lep_ID->push_back((muon->Charge)*-13);
      lep_sumpt_charged->push_back(muon->SumPtCharged);

          
    }

 
  
    phot_n = 0;
    phot_pt->clear();
    phot_eta->clear();
    phot_phi->clear();
    phot_E->clear();

    for(Int_t i = 0; i < branchPhoton->GetEntries(); ++i)
    {
      photon = (Photon*) branchPhoton->At(i);
      phot_n++;
      phot_pt->push_back(photon->PT);
      phot_eta->push_back(photon->Eta);
      phot_phi->push_back(photon->Phi);
      phot_E->push_back(photon->E);
      

          
    }

  
 if(branchMissingET->GetEntries() > 0)
    {
      // Take MET
      MeT = (MissingET*) branchMissingET->At(0); 
      met = MeT->MET;
      met_phi = MeT->Phi;
      met_eta = MeT->Eta;
}
      
      for(Int_t i = 0; i < branchScalarHT->GetEntries(); ++i)
    {
      Scalar_HT = (ScalarHT*) branchScalarHT->At(i);
      HT = Scalar_HT->HT;
          
    }

     jet_n = 0;
     bjet_n = 0;

     jet_pt->clear();
     jet_eta->clear();
     jet_phi->clear();
     jet_mass->clear();
     jet_flavor->clear();
     jet_btag->clear();
     jet_tautag->clear();
     jet_ncharged->clear();

  for(Int_t i = 0; i < branchJet->GetEntriesFast(); ++i)
    {
      // Take the jet
     jet = (Jet*) branchJet->At(i); 
 
      jet_n++;
    // count the b-jets
      if(jet->BTag==1&&jet->PT>40){bjet_n++;}

      jet_pt->push_back(jet->PT);
      jet_eta->push_back(jet->Eta);
      jet_phi->push_back(jet->Phi);
      jet_mass->push_back(jet->Mass);
      jet_flavor->push_back(jet->Flavor);
      jet_btag->push_back(jet->BTag);
      jet_tautag->push_back(jet->TauTag);
      jet_ncharged->push_back(jet->NCharged);

    


   }

// Fill inside the tree to all the created branches

 delphes.Fill();   
  
  }

delphes.Write("",TObject::kWriteDelete);    // TObject::kWriteDelete   to remove the old snapshot and save the latest copy
//
// save the Tree ; the file will be automatically closed when going out oth function scope

 
}

