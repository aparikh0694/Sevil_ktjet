void kt_quench_test(Int_t nEvents=1)
{
  gROOT->Reset();
  gROOT->Clear();

  // ktJet lib
  if (gClassTable->GetID("ktJet") < 0) {
    cout<<"Load ktJet lib ..."<<endl;
    gSystem->Load("libKtJet.so");
  }

  cout<<endl;
  cout<<"    Test ktJetQuench'"<<endl;
  cout<<" ------------------------"<<endl;
  cout<<"  (Joern Putschke, 2008)"<<endl;
  cout<<endl;

  Int_t nEv=0;
  Double_t mPtCut=2.0;
  Double_t mPtCutMin=0.1;
  Double_t mRc=0.7;
  Double_t mRcQuench=0.7;
  Double_t mS0=0.25;

  ktPythia8 *p=new ktPythia8("/Users/putschke/pythia8100/xmldoc");
  p->SetVerbose(false);
  //p->InitPhaseSpace("28.0","32.0");
  p->InitPhaseSpace("12.0","60.0");
  p->Init("RHIC");

  ktMuEvent *ev=0;
  ktMuEvent *evq=0;

  ktJetQuench *q=new ktJetQuench();
  q->SetQuenchingFraction(mS0);
  q->SetQuenchConstant(2.0);
  q->SetVerbose(false);
  q->SetDoQA(true);
  q->SetFilledFromFJ(true);
  //q->SetQuenchedSlope(0.4); // ridge temp.
  //q->InitQuenching("constant");
  q->InitQuenching("fractional");
  //q->PrintInfo();
  cout<<endl;

  ktJetQuench *qPy=new ktJetQuench();
  qPy->SetQuenchingFraction(mS0);
  qPy->SetQuenchConstant(2.0);
  qPy->SetVerbose(false);
  qPy->SetDoQA(true);
  qPy->SetFilledFromFJ(false);
  //qPy->SetQuenchedSlope(0.4); // ridge temp.
  //qPy->InitQuenching("constant");
  qPy->InitQuenching("fractional");
  cout<<endl;

  do
    {
      
      p->RunWithPycell();
      
      //if (!(fabs(p->GetPycellEta())<(1-mRc) && p->GetPycellPhi()<(TMath::Pi()-mRc) &&  p->GetPycellPhi()>mRc && p->GetPycellEt()>28.0 &&  p->GetPycellEt()<32.0)) continue;
      if (!(fabs(p->GetPycellEta())<(1-mRc) && p->GetPycellPhi()<(TMath::Pi()-mRc) &&  p->GetPycellPhi()>mRc && p->GetPycellEt()>12.0)) continue;
      
      nEv++;

      if (nEvents>9 && fmod(nEv,(Int_t) nEvents/10)==0)
        cout<<"Event number : "<<nEv<<endl;


      ev=new ktMuEvent(nEv);
      evq=new ktMuEvent(nEv);
   
      ktFastJet *fastJet=new ktFastJet();
      fastJet->SetFiducial(1,TMath::Pi());
      fastJet->SetPtCut(mPtCutMin);
      fastJet->SetInfo(false);
      fastJet->SetPrintJet(false);
      fastJet->SetFillFFArray(true);
      fastJet->SetNJetsForFF(2);
      
      if (nEv<1)
	fastJet->PrintInfo();	  
     
      fastJet->AddPythia8Event(p,"EMCAL");
      fastJet->GetHighestJetForQuenching(q,mRcQuench,"Antikt");

      p->GetHighestJetForQuenching(qPy,mRcQuench,"EMCAL",false);

      // otherwise in Pythia8 ideal jet pt ...
      // modify Pythia8 accordingly ... (!?) or use FastJet definition ...
      qPy->SetJetPt(q->GetJetPt());

      q->DoQuenching();
      qPy->DoQuenching();

      //qPy->PrintInfo();
      q->DoDiHadron(3,4,2);
      qPy->DoDiHadron(3,4,2);

      // DEBUG:
      /*
      cout<<endl;
      q->PrintInfo();
      cout<<endl;

      cout<<endl;
      qPy->PrintInfo();
      cout<<endl;
      */
     
      //fastJet->SetPtCut(mPtCut);
      /*
      fastJet->Clear();
      fastJet->AddUnQuenchedJet(qPy);
      fastJet->RunFastJetSub(mRc,"kt",ev);
      fastJet->RunFastJetSub(mRc,"Antikt",ev);

      fastJet->Clear();
      fastJet->AddQuenchedJet(qPy);
      fastJet->RunFastJetSub(mRc,"kt",evq);    
      fastJet->RunFastJetSub(mRc,"Antikt",evq);

      cout<<endl;
      ev->PrintJets(0);
      cout<<endl;
      evq->PrintJets(0);
      //cout<<endl;
      */

      q->Clear();
      qPy->Clear();

      delete fastJet;      
      delete ev; delete evq;	
      
      //gObjectTable->Print();
      
    }
  while (nEv<nEvents);
 
  //q->DrawQAPlots("FastJet");
  qPy->DrawQAPlots("Pythia8");

  
  cout<<endl;
  //q->SaveQA("test_quench_fj.root");
  //qPy->SaveQA("test_quench_py.root");
  

  delete p;
  delete qPy;
  delete q;

  //cout<<endl;
  //cout<<anz<<" particles in test file."<<endl;

  cout<<endl;
  cout<<"Done :-)"<<endl;
  cout<<endl;

}
