void kt_track_dist(Int_t nEvents=1)
{
  gROOT->Reset();
  gROOT->Clear();

  // ktJet lib
  if (gClassTable->GetID("ktJet") < 0) {
    cout<<"Load ktJet lib ..."<<endl;
    gSystem->Load("libKtJet.so");
  }

  cout<<endl;
  cout<<" Test track distoritions"<<endl;
  cout<<" ------------------------"<<endl;
  cout<<"  (Joern Putschke, 2008)"<<endl;
  cout<<endl;

  Int_t nEv=0;
  Double_t mPtCut=2.0;
  Double_t mPtCutMin=0.1;
  Double_t mRc=0.4;

  ktPythia8 *p=new ktPythia8("/home/putschke/pythia8100/xmldoc");
  p->SetVerbose(false);
  p->InitPhaseSpace("28.0","32.0");
  //p->InitPhaseSpace("12.0","60.0");
  p->Init("RHIC");

  ktMuEvent *ev=0;
  ktMuEvent *evd=0;

  TH1D *hE=new TH1D("hE","",75,0,75);
  TH1D *hEd=new TH1D("hEd","",75,0,75);
  
  TH1D *hXi=new TH1D("hXi","",25,0,5);
  TH1D *hXid=new TH1D("hXid","",25,0,5);

  TH1D *hz=new TH1D("hz","",20,0,1);
  TH1D *hzd=new TH1D("hzd","",20,0,1);

  int norm=0;

  do
    {
      
      p->RunWithPycell();
      
      if (!(fabs(p->GetPycellEta())<(1-mRc) && p->GetPycellPhi()<(TMath::Pi()-mRc) &&  p->GetPycellPhi()>mRc && p->GetPycellEt()>28.0 &&  p->GetPycellEt()<32.0)) continue;
      //if (!(fabs(p->GetPycellEta())<(1-mRc) && p->GetPycellPhi()<(TMath::Pi()-mRc) &&  p->GetPycellPhi()>mRc && p->GetPycellEt()>12.0)) continue;
      
      nEv++;

      if (nEvents>9 && fmod(nEv,(Int_t) nEvents/10)==0)
        cout<<"Event number : "<<nEv<<endl;


      ev=new ktMuEvent(nEv);
      evd=new ktMuEvent(nEv);
   
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
      fastJet->RunFastJetSub(mRc,"AntiKt",ev);
      fastJet->FillFF("AntiKt",ev);

      ev->SetMedianPtPerArea(fastJet->GetMedianPtPerArea());	  
      ev->SetParticlePid("Track Simulation");
      ev->SetBkgPtCut(mPtCutMin);
      ev->SetRParam(mRc);

      //ev->PrintJets(10);
      //cout<<endl;

      delete fastJet;

      ktFastJet *fastJetd=new ktFastJet();
      fastJetd->SetFiducial(1,TMath::Pi());
      fastJetd->SetPtCut(mPtCutMin);
      fastJetd->SetInfo(false);
      fastJetd->SetPrintJet(false);
      fastJetd->SetFillFFArray(true);
      fastJetd->SetNJetsForFF(2);
      fastJetd->SetSpaceCharge(true);
      //fastJetd->SetTrackSmear(true);

      fastJetd->AddPythia8Event(p,"EMCAL");
      fastJetd->RunFastJetSub(mRc,"AntiKt",evd);
      fastJetd->FillFF("AntiKt",evd);

      evd->SetMedianPtPerArea(fastJetd->GetMedianPtPerArea());	  
      evd->SetParticlePid("Track Simulation");
      evd->SetBkgPtCut(mPtCutMin);
      evd->SetRParam(mRc);

      //evd->PrintJets(10);

      delete fastJetd;

      if (ev->GetNFJAntiKt()>0 && evd->GetNFJAntiKt()>0)
	{
	  norm++;

	  // Fill jet pt ...
	  hE->Fill(ev->GetFJAntiKt(0)->Pt());
	  hEd->Fill(evd->GetFJAntiKt(0)->Pt());

	  // Fill FF (highest) ...
	  j=(ktMuFastJet*) ev->GetFJAntiKt(0);
	  jd=(ktMuFastJet*) evd->GetFJAntiKt(0);

	  for (int i=0;i<j->GetNJetParticles();i++)
	    {
	      mp=(ktParticle*) j->GetJetParticle(i);
	      
	      if (mp->IsCharged())
		{
		  Double_t z=mp->Pt()/j->Pt();
		  Double_t xi=log(1/z);
		  
		  hXi->Fill(xi);
		  hz->Fill(z);
		}	      
	    }

	  for (int i=0;i<jd->GetNJetParticles();i++)
	    {
	      mp=(ktParticle*) jd->GetJetParticle(i);
	      
	      if (mp->IsCharged())
		{
		  Double_t z=mp->Pt()/jd->Pt();
		  Double_t xi=log(1/z);
		  
		  hXid->Fill(xi);
		  hzd->Fill(z);
		}	      
	    }	  
	}

      delete ev; delete evd;	
      
      //gObjectTable->Print();
      
    }
  while (nEv<nEvents);
  
  cout<<endl;
  
  TCanvas *c0=new TCanvas("c0","#0",800,600);
  hE->DrawCopy();
  hEd->SetLineColor(2);
  hEd->DrawCopy("same");

  TCanvas *c1=new TCanvas("c1","#1",800,600);
  hXi->Scale(1/(Double_t) norm*1/hXi->GetBinWidth(1));
  hXid->Scale(1/(Double_t) norm*1/hXid->GetBinWidth(1));
  
  hXi->DrawCopy();
  hXid->SetLineColor(2);
  hXid->DrawCopy("same");

  TCanvas *c2=new TCanvas("c2","#2",800,600);
  hz->Scale(1/(Double_t) norm*1/hz->GetBinWidth(1));
  hzd->Scale(1/(Double_t) norm*1/hzd->GetBinWidth(1));

  hz->DrawCopy();
  hzd->SetLineColor(2);
  hzd->DrawCopy("same");

  delete p;

  //cout<<endl;
  //cout<<anz<<" particles in test file."<<endl;

  cout<<endl;
  cout<<"Done :-)"<<endl;
  cout<<endl;

}
