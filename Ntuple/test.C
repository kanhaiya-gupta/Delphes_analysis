void test(){
  // Declare here for histos and so on to keep info
                TFile * Output          = new TFile("testfile.root","RECREATE");

                cout << "Define Tree and setup branches\n";
                //Store skimmed variables after a few of the cuts but not all
                TTree* Loose = new TTree("Loose","Loose");


                std::vector<float> *vpx, *vpy, *vpz;
              //  std::vector<float>* vpy;
              //  std::vector<float>* vpz;

                //SetBranches for Electrons
                Loose->Branch("vpx", &vpx);
                Loose->Branch("vpy", &vpy);
                Loose->Branch("vpz", &vpz);

                for(int i =0; i<100; i++){
                        vpx->clear();
                        vpy->clear();
                        vpz->clear();

                        for(int j =0; j<10; j++){
                                vpx->push_back(j);
                                vpy->push_back(1);
                                vpz->push_back(3);
                        }

                        Loose->Fill();
                }

                Output->Write();

}//END TEST
