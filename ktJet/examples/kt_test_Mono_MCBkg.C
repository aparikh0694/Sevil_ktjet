void kt_test_Mono_MCBkg(Int_t nEvents=1,Double_t mv2=0.1,TString mfOutName="test.root")
{
  Double_t mRc=0.4;
  Double_t mPtCut=0.1;
  Double_t mSeed=5.0;
  
  /*
  ktMCBkg *b=new ktMCBkg("RHIC","Au+Au",mv2);
  ktPythia8 *p=new ktPythia8("/Users/putschke/pythia8100/xmldoc");
  // choose path to your pythia8 xmldoc installation directory !
  //p->SetRandomSeed();
  p->InitPhaseSpace("28.0","32.0");
  p->Init("RHIC");
  */

  
  ktMCBkg *b=new ktMCBkg("LHC","Pb+Pb",mv2);
  ktPythia8 *p=new ktPythia8("/home/putschke/pythia8100/xmldoc");
  // choose path to your pythia8 xmldoc installation directory !
  //p->SetRandomSeed();
  // jet et cut needed in addition at least for LHC setting !!!
  p->InitPhaseSpace("98.0","102.0");
  p->Init("LHC");
  

  Int_t nEv=0;

  //TFile *f=new TFile(fileName,"RECREATE");

  ktMuEvent *ev=0;
  ktMuEvent *ev2=0;

  // some fast out histo ...
  /*
  TH1D *hE=new TH1D("hE","Energy (Cone)",100,0,100);
  TH1D *hEp=new TH1D("hE","PYTHIA8 pthat",100,0,100);
  TH1D *hECellCorr=new TH1D("hECellCorr","Energy (Cone) cell corr",100,0,100);
  TH1D *hEConeCorr=new TH1D("hEConeCorr","Energy (Cone) cone corr",100,0,100);
  TH1D *hECellCorrIn=new TH1D("hECellCorrIn","Energy (Cone) cell corr in-plane",100,0,100);
  TH1D *hEConeCorrIn=new TH1D("hEConeCorrIn","Energy (Cone) cone corr in-plane",100,0,100);
  TH1D *hECellCorrOut=new TH1D("hECellCorrOut","Energy (Cone) cell corr out-of-plane",100,0,100);
  TH1D *hEConeCorrOut=new TH1D("hEConeCorrOut","Energy (Cone) cone corr out-of-inplane",100,0,100);
  */

  TFile *fout=new TFile(mfOutName,"RECREATE");

  TH1D *hE=new TH1D("hE","Energy (Cone)",200,0,200);
  TH1D *hEp=new TH1D("hEp","PYTHIA8 pthat",200,0,200);
  TH1D *hECellCorr=new TH1D("hECellCorr","Energy (Cone) cell corr",200,0,200);
  TH1D *hEConeCorr=new TH1D("hEConeCorr","Energy (Cone) cone corr",200,0,200);
  
  //TH1D *hECellCorrIn=new TH1D("hECellCorrIn","Energy (Cone) cell corr in-plane",200,0,200);
  //TH1D *hEConeCorrIn=new TH1D("hEConeCorrIn","Energy (Cone) cone corr in-plane",200,0,200);
  //TH1D *hECellCorrOut=new TH1D("hECellCorrOut","Energy (Cone) cell corr out-of-plane",200,0,200);
  //TH1D *hEConeCorrOut=new TH1D("hEConeCorrOut","Energy (Cone) cone corr out-of-inplane",200,0,200);
  

  /*
  TTree *outTree=new TTree("Jpt","Jet event tree");
  outTree->Branch("event",&ev);
  TTree *outTree2=new TTree("Jpppt","Jet event tree");
  outTree2->Branch("event",&ev2); 
  */

  if (b->Init())
    {
      //for (int i=0;i<nEvents;i++)
      do
	{

	  p->RunWithPycell();

	  //if (!(fabs(p->GetPycellEta())<(1.0-mRc) && p->GetPycellPhi2()<(TMath::Pi()-mRc) && p->GetPycellEt()>96.0 &&  p->GetPycellEt()<104.0)) {cout<<" *** PYTHIA8 not fullfilling the eta/phi phacespace cut !"<<endl; continue;}
	  if (!(fabs(p->GetPycellEta())<(1.0-mRc) && p->GetPycellPhi2()<(TMath::Pi()-mRc) && p->GetPycellEt()>90.0 &&  p->GetPycellEt()<110.0)) continue;
	      
	  ktGrid *grid=new ktGrid();
	  grid->SetGrid((int) 1*180,0,2*TMath::Pi(),(int) 100/2,-1.0,1.0); // grid2 (2xEMCAL)
          grid->SetEmcal(0.01,2*TMath::Pi()-0.01,-0.99,0.99);
          grid->SetSeed(mSeed);
          grid->SetBkgPtCut(mPtCut);
          grid->SetCone(mRc);
          grid->SetBkgCone(0.2);
          grid->SetFFCone(0.7); 
          grid->SetRMaxNN(1.0); // for kt-Jetfinder not used here ...
          grid->SetNBkgSeeds(250);
	  grid->Init();

	  ktGrid *grid2=new ktGrid();
	  grid2->SetGrid((int) 1*180,0,2*TMath::Pi(),(int) 100/2,-1.0,1.0); // grid2 (2xEMCAL)
          //grid->SetEmcal(MinPhiEmcal,MaxPhiEmcal,-MaxEtaEmcal,MaxEtaEmcal);
          grid2->SetSeed(mSeed);
          grid2->SetBkgPtCut(mPtCut);
          grid2->SetCone(mRc);
          grid2->SetBkgCone(0.2);
          grid2->SetFFCone(0.7); 
          grid2->SetRMaxNN(1.0); // for kt-Jetfinder not used here ...
          grid2->SetNBkgSeeds(250);
	  grid2->Init();
	  
	  nEv++;

	  p->Fill(grid,"ALL");
	  p->Fill(grid2,"ALL");

	  //b->Fill(grid);
	  // DEBUG:
	  //cout<<p->GetPycellPhi()<<" "<<p->GetPycellPhi2()<<endl;

	  b->Fill(grid,p->GetPycellPhi2());//+TMath::Pi()/2.0));
	  //b->Fill(grid);

	  cout<<endl;
	  cout<<(Int_t) grid->GetNCells()<<" cells in grid acceptance + cuts."<<endl;

	  grid->GetSeeds();
	  grid->DoJetfinding("coneBkg");
	  grid->CalcBkgNJetsRemove(2);
	  grid->BkgRandomCones();
	  
	  grid2->GetSeeds();
	  grid2->PrintSeeds();
	  grid2->DoJetfinding("coneBkg");

	  // check hard scattering values eta/phi from PYTHIA ...
	  p->HardProcessInfo();
	  p->Pycell();
	  cout<<endl;
	  grid2->PrintJetsNice();
	  grid->PrintJetsNiceBkgSub();
	  
	  /*
	  if (nEvents<5)
	    {
	      
	      TCanvas *c0=new TCanvas("c0","#0",800,600);
	      //c0->SetLogz();
	      grid->hGrid->DrawCopy("lego2");
	    }
	  */

	  TObjArray *foundJets=grid->GetJetList();     
	  ev=new ktMuEvent(nEv);
	  ev->FillBkgInfo(grid);
	  ev->FillJetFinderSettings(grid);
	  ev->SetRefMult(999);
	  ev->SetParticlePid("ALL Pb+Pb");
	  ev->Fill(foundJets,grid->GetBkgPtPerCellIn(),grid->GetBkgPtPerCellOut(),grid->GetBkgPtCone(),grid->GetBkgPtCluster());      
	  //outTree->Fill();

	   TObjArray *foundJets2=grid2->GetJetList();     
	   ev2=new ktMuEvent(nEv);
	   ev2->Fill(foundJets2);
	   ev2->FillJetFinderSettings(grid2);
	   ev2->SetParticlePid("ALL p+p");
	   ev2->SetRefMult(0);
	   //outTree2->Fill();

	  // DEBUG:
	  //cout<<grid->GetBkgPtCone()<<endl;
	  //cout<<grid->GetBkgPtPerCellIn()<<" "<<grid->GetBkgPtPerCellOut()<<endl;

	  if (ev->GetNConeJets()>0 && ev2->GetNConeJets()>0)
	    {
	      ktMuJet *temp=(ktMuJet*) ev->GetConeJet(0);
	      hECellCorr->Fill(temp->PtCellCorr());
	      hEConeCorr->Fill(temp->PtConeCorr());
	     
	      ktMuJet *temp2=(ktMuJet*) ev2->GetConeJet(0);
	      hE->Fill(temp2->Pt());
	      hEp->Fill(p->GetPtHat());
	      temp2->PrintSeed();
	    }

	  delete ev;
	  delete ev2;
	  
	  cout<<" ==== Event "<<nEv<<" done !"<<endl;

	  delete grid;
	  delete grid2;
	  

	  //cout<<" all deleted !"<<endl;

	  //gObjectTable->Print();

	}
      while (nEv<nEvents);
    }
  else
    cout<<" Error :  MC bkg. not initialized !"<<endl;

  cout<<endl;

  p->Statistics();

  if (b->GetInitOK())
    {
      b->HistNorm(nEvents);
      TCanvas *c1=new TCanvas("c1","Canvas #1",600,600);
      c1->Divide(2,2);
      c1->cd(1);b->heta->DrawCopy();
      c1->cd(2);b->hpt->SetMinimum(0.001);
      b->hpt->DrawCopy();gPad->SetLogy();
      
      c1->cd(3);b->hphi->DrawCopy();
      if (b->Getv2()>0)
	{
	  c1->cd(4);b->hRP->DrawCopy();
	}
    } 
  
  TCanvas *c2=new TCanvas("c2","Canvas #2",800,600);
  Int_t nrebin=5;
  hE->Rebin(nrebin);hECellCorr->Rebin(nrebin);
  hEConeCorr->Rebin(nrebin);hEp->Rebin(nrebin);

  hE->DrawCopy();
  hECellCorr->SetLineColor(2);
  hEConeCorr->SetLineColor(4);
  hECellCorr->DrawCopy("same");
  hEConeCorr->DrawCopy("same");
  hEp->SetLineStyle(2);
  hEp->DrawCopy("same");

  fout->Write();
  fout->Close();

  delete p;
  delete b;
}
