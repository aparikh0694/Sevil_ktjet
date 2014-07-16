void kt_test_hijing(Int_t nEvents=100)
{

  gROOT->Reset();
  gROOT->Clear();

  // ktJet lib
  if (gClassTable->GetID("ktJet") < 0) {
    cout<<"Load ktJet lib ..."<<endl;
    gSystem->Load("libKtJet.so");
  }

  Double_t mPtCut=0.2;
  Double_t mPtCutMin=0.2;
  Double_t mRc=0.7;

  ktHijing *h=new ktHijing("hijing","HIJING");
  h->SetJetPtRange(15,-1);
  h->SetJetEtaRange(-0.3,0.3);

  h->InitRHICdAu();
  cout<<endl;

  ktFastJet *fastJet=new ktFastJet();
  fastJet->SetFiducial(1,TMath::Pi());
  fastJet->SetPtCut(mPtCut);
  fastJet->SetInfo(true);
  fastJet->SetPrintJet(false);
  fastJet->SetFillFFArray(true);
  fastJet->SetNJetsForFF(2);
  fastJet->SetActiveAreaRepeats(1);
  fastJet->PrintInfo();	  
  
  for (int i=0;i<nEvents;i++)
    {

      if (fmod(i,nEvents/10)==0 && i>0) {cout<<i<<" HIJING events done ..."<<endl;}
      
      h->MakeMinbiasEvent();
      
      ktMuEvent *ev=new ktMuEvent(i);

      h->FillFastJet(fastJet,2);

      // Do FastJet analysis ...	  
      //fastJet->DoBkg();
      //cout<<endl;
      
      fastJet->RunFastJetSub(mRc,"kt",ev);
      fastJet->RunFastJetSub(mRc,"AntiKt",ev);

      ev->PrintJets(5);
      cout<<endl;     

      fastJet->Clear();
      delete ev;
    }

  //h->GetEtaHist()->DrawCopy();
  delete fastJet;

  cout<<endl;

  delete h;
}
