void kt_fastjet_test(Int_t nEvents=10,TString mfName="test_fj_out.root")
{
  gROOT->Reset();
  gROOT->Clear();

  // ktJet lib
  if (gClassTable->GetID("ktJet") < 0) {
    cout<<"Load ktJet lib ..."<<endl;
    gSystem->Load("libKtJet.so");
  }

  cout<<endl;
  cout<<" Test FastJet Interface"<<endl;
  cout<<" -----------------------------------------"<<endl;
  cout<<"      (Salur, 2013)"<<endl;
  cout<<endl;

  Int_t nEv=0;
  Double_t mSeed=4.6;
  Double_t mPtCut=2;
  Double_t mPtCutMin=0.2;
  Double_t mRc=0.4;
  Double_t mRcFF=0.7;
  Double_t mRcPycell=0.7;

  Bool_t mMuWrite=true;

  ktMCBkg *b=new ktMCBkg("LHC","p+p",0.0);
  ktPythia8 *p=new ktPythia8("/Users/salur/software/pythia8130/xmldoc/");
  p->SetVerbose(false);
  p->InitPhaseSpace("100.0","105.0");
  p->Init("LHC");



  ktMuEvent *ev=0;
  TFile *ffout=0;

   if (mMuWrite)
    {
      ffout=new TFile(mfName,"RECREATE");
      outTree=new TTree("JSimuFJ","Jet event tree");
      outTree->Branch("event",&ev);
    }

  if (b->Init())
    {
      do
	{
	  
	  p->RunWithPycell();
	  
	  if (!(fabs(p->GetPycellEta())<(0.3) && p->GetPycellPhi2()<(TMath::Pi()-0.7) && p->GetPycellEt()>100.0 &&  p->GetPycellEt()<105.0)) continue;

	  ev=new ktMuEvent(nEv);
	  
	  nEv++;

	  ktFastJet *fastJet=new ktFastJet();
	  fastJet->SetFiducial(1,TMath::Pi());
	  fastJet->SetPtCut(mPtCut);
	  fastJet->SetInfo(true);
	  fastJet->SetPrintJet(false);
	  fastJet->SetFillFFArray(true);
	  fastJet->SetNJetsForFF(2);

	  if (nEv<1)
	    fastJet->PrintInfo();	  
	  
	  // STAR EMCal tower size: 0.05 x 0.05
	  ktGrid *grid=new ktGrid();
	  grid->SetGrid((int) 125/1,0,2*TMath::Pi(),(int) 40,-1.0,1.0); // grid2 (2xEMCAL)
	  grid->SetEmcal(0.05,2*TMath::Pi()-0.05,-0.95,0.95);
	  grid->SetSeed(mSeed);
	  grid->SetBkgPtCut(mPtCut);
	  grid->SetCone(mRc);
	  grid->SetBkgCone(mRc);
	  grid->SetFFCone(mRcFF);
	  grid->SetRMaxNN(1.0); // for kt-Jetfinder not used here ...
	  grid->SetNBkgSeeds(150);
	  grid->SetVerbose(0);
	  grid->Init();
	  
	  ktGrid *grid2=new ktGrid();
	  grid2->SetGrid((int) 125/1,0,2*TMath::Pi(),(int) 40,-1.0,1.0); // grid2 (2xEMCAL)
	  grid2->SetEmcal(0.05,2*TMath::Pi()-0.05,-0.95,0.95);
	  grid2->SetSeed(mSeed);
	  grid2->SetBkgPtCut(mPtCut);
	  grid2->SetCone(mRc);
	  grid2->SetBkgCone(mRc);
	  grid2->SetFFCone(mRcFF);
	  grid2->SetRMaxNN(1.0); // for kt-Jetfinder not used here ...
	  grid2->SetNBkgSeeds(150);
	  grid2->SetVerbose(0);
	  grid2->Init();

	  // Add Pythia jets ...
	  fastJet->AddPythia8Event(p,"EMCAL");
	  p->Fill(grid,"EMCAL");

	  // Add Background ...
	  b->Fill(grid,fastJet,0);
	 
	  // Do FastJet analysis ...	  
	  fastJet->DoBkg();
	  cout<<endl;

	  fastJet->RunFastJetSub(mRc,"kt",ev);
	  fastJet->RunFastJetSub(mRc,"AntiKt",ev);
	  //fastJet->RunSISCone(mRc,ev);
	  
	  // Do FF ...
	  /*
	  fastJet->FillFF("SISCone",ev);
	  fastJet->DoXiBkg("SISCone",ev);
	  fastJet->FillFF("AntiKt",ev);
	  fastJet->DoXiBkg("AntiKt",ev);
	  fastJet->FillFF("kt",ev);		
	  fastJet->DoXiBkg("kt",ev);
	  */
	  
	  // Do LOCone analysis ...
	  if (grid->GetSeeds()) 
	    {
	      //cout<<endl;	    
	      cout<<" ---> Ran LOCone with R = "<<mRc<<endl;
	      cout<<endl;

	      grid->DoJetfinding("coneBkg");
	      grid->CalcBkgNJetsRemove(2);
	      grid->BkgRandomCones();
	      grid->BkgXi();
	  
	      //cout<<endl;	      
	      //grid->PrintJetsNiceBkgSub();
	      
	      if (nEvents<5)
		{		  
		  TCanvas *c0=new TCanvas("c0","#0",800,600);
		  //c0->SetLogz();
		  grid->hGrid->DrawCopy("lego2");
		}	      
	    }
	  
	  TObjArray *foundJets=grid->GetJetList();     
	  ev->FillBkgInfo(grid);
	  ev->SetMedianPtPerArea(fastJet->GetMedianPtPerArea());	  
	  ev->SetParticlePid("Simulation test Au+Au");
	  ev->Fill(foundJets,grid->GetBkgPtPerCellIn(),grid->GetBkgPtPerCellOut(),grid->GetBkgPtCone(),grid->GetBkgPtCluster());      

	  ev->PrintJets(10);
	  cout<<endl;
	  	  
	   if (mMuWrite)
	     {
	       outTree->Fill();
	     }

	  delete fastJet;
	  delete grid;
	  delete grid2;
	  
	  delete ev;	

	  //gObjectTable->Print();
	  
	}
      while (nEv<nEvents);
    }
  else
    cout<<" Error :  MC bkg. not initialized !"<<endl;
  
  if (mMuWrite)
    {
       ffout->Write();
       ffout->Close();
    }
  
  delete p;
  delete b;

  //cout<<endl;
  //cout<<anz<<" particles in test file."<<endl;

  cout<<endl;
  cout<<"Done :-)"<<endl;
  cout<<endl;

}
