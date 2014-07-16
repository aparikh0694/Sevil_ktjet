void kt_test_pico(TString fin="test.root") 
{
  
  gROOT->Reset();
  gROOT->Clear();
  
  // ktJet lib
  if (gClassTable->GetID("ktJet") < 0) {
    cout<<"Load ktJet lib ..."<<endl;
    gSystem->Load("libKtJet.so");
  }

  cout<<endl;
  cout<<" Test STAR pico Dst interface"<<endl;
  cout<<" ----------------------------"<<endl;
  cout<<endl;
  cout<<"Open STAR pico Dst file : "<<fin<<endl;
  cout<<endl;

  TFile *inFile = TFile::Open(fin);
  inFile->cd();
  TTree *outT = inFile->Get("JetTree");

  TStarJetPicoEvent *event = new TStarJetPicoEvent();
  TBranch *branch = outT->GetBranch("PicoJetTree");
  branch->SetAddress(&event);
  
  Int_t ievents = outT->GetEntries();

  cout<<"Number of events = "<<ievents<<endl;

  ktStarPico *mQA=new ktStarPico(true);
  mQA->PrintQACuts();

  for (Int_t i = 0; i<ievents; i++)
    {
      if (i>250) continue;
      
      outT->GetEvent(i);

      mQA->DoQAOnly(event);

    }

  mQA->WriteQAOutputFile("QA_test_out.root");

  //inFile->Close();

  delete mQA;

  cout<<endl;
  cout<<"Done ;-)"<<endl;
  cout<<endl;
}
