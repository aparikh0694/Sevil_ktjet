void kt_fastjet_hiforest(Int_t nEvents=10,TString mfName="test_fj_out.root")
{
  gROOT->Reset();
  gROOT->Clear();

  // ktJet lib
  if (gClassTable->GetID("ktJet") < 0) {
    cout<<"Load ktJet lib ..."<<endl;
    gSystem->Load("libKtJet.so");
    gSystem->Load("hiForest_h.so");
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
  char file_name1[512];
  char file_out[512];
 

  sprintf(file_name1,"/Users/salur/software/CMSJet/HydjetDrum03_HiForest_v05_merged_test02.root");
  TFile *FileA = TFile::Open(file_name1);

  sprintf(file_out,"DiJetQAv4_%i_%i.root",1,2);
  TFile* outf = new TFile(file_out,"recreate");
  
  cout<<"file name  "<<file_out<<endl;
    
  Int_t pHBHENoiseFilter;
  Int_t pCollisionEventSelection;
  Int_t phiEcalRecHitSpikeFilter;
  Int_t HLT_HIJet80_v1;

  int run;
  int evt;
  int hibin;

  int RunNum;
  int EventNum;
  
  int fNbinspt=800;
  int fBinlowpt=0;
  int fBinuppt=800;

 
  TTree* pfTree = (TTree*)FileA->Get("pfcandAnalyzer/pfTree");
  Int_t nPFpart;
  Int_t pfId[65536];
  Float_t pfPt[65536];
  Float_t pfEta[65536];
  Float_t pfPhi[65536];
  
  pfTree->SetBranchAddress("nPFpart", &nPFpart);
  pfTree->SetBranchAddress("pfId", pfId);
  pfTree->SetBranchAddress("pfPt", pfPt);
  pfTree->SetBranchAddress("pfEta", pfEta);
  pfTree->SetBranchAddress("pfPhi", pfPhi);
  

  cout<<"test part"<<nPFpart<<endl;
  
  TTree* t = (TTree*)FileA->Get("hltanalysis/HltTree");
  TTree* skimt = (TTree*)FileA->Get("skimanalysis/HltTree");
  TTree* Evtt = (TTree*)FileA->Get("hiEvtAnalyzer/HiTree");
  
  skimt->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);
  skimt->SetBranchAddress("pcollisionEventSelection",&pCollisionEventSelection);
  skimt->SetBranchAddress("phiEcalRecHitSpikeFilter",&phiEcalRecHitSpikeFilter);
  

  t->SetBranchAddress("HLT_HIJet80_v1",&HLT_HIJet80_v1);
  t->SetBranchAddress("Run",&run);
  t->SetBranchAddress("Event",&evt);
  
  Evtt->SetBranchAddress("hiBin",&hibin);
  
  int Nevents = t->GetEntries();
  //int evcts=1;
  
  cout<<" Number of events  "<<t->GetEntries()<<endl;
  
  //cout<<" hibin "<<hibin<<endl;
  int Nevents = 20;


  ktMuEvent *ev=0;
    TFile *ffout=0;
    
    if (mMuWrite)
      {
	ffout=new TFile(mfName,"RECREATE");
	outTree=new TTree("JSimuFJ","Jet event tree");
	outTree->Branch("event",&ev);
      } 
  for(int iev = 0; iev < Nevents; ++iev){

    t->GetEntry(iev);
    skimt->GetEntry(iev);
    Evtt->GetEntry(iev);
    pfTree->GetEntry(iev);
    
    if (iev%1000==0) cout <<iev<<" / "<<Nevents <<endl;
    
    //    if(pCollisionEventSelection==0 || pHBHENoiseFilter==0 ||HLT_HIJet80_v1==0 ||  phiEcalRecHitSpikeFilter ==0)continue;
    
   
    
    //some event selection is needed
    ev=new ktMuEvent(iev);
    
    nEv++;
    
    ktFastJet *fastJet=new ktFastJet();
    fastJet->SetFiducial(2.6,TMath::Pi());
    fastJet->SetPtCut(mPtCut);
    fastJet->SetInfo(true);
    fastJet->SetPrintJet(false);
    fastJet->SetFillFFArray(false);
    fastJet->SetNJetsForFF(2);
    
    if (nEv<1)
      fastJet->PrintInfo();	  
    
    // STAR EMCal tower size: 0.05 x 0.05
    ktGrid *grid=new ktGrid();
    grid->SetGrid((int) 125/1,0,2*TMath::Pi(),(int) 200,-2.5,2.5); // grid2 (2xEMCAL)
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
    
    //cout<<"  test "<<nPFpart<<endl;

    for (Int_t k = 0; k < nPFpart; k++) {
      
      TLorentzVector *p=new TLorentzVector();
      p->SetPtEtaPhiM(pfPt[k],pfEta[k],pfPhi[k],0);

      //cout<<"   "<<pfPt[k]<<"   "<<k<<endl;
      grid->Fill((TLorentzVector*) p->Clone(),true);
      //      p->Fill(grid,fastJet,0);

      //      fastJet->AddParticle((TLorentzVector*) p->Clone()); 
     
      Double_t ax, ay, az, aE;
      ax=(TLorentzVector*)p->Px();
      ay=p->Py();
      az=p->Pz();
      aE=p->E();
      // cout<<"  ay "<<ay<<endl;
      fastJet->AddParticle(ax,ay,az,aE);
      //fastJet->AddParticle((TLorentzVector*)p->Px(),(TLorentzVector*)p->Py(),(TLorentzVector*)p->Pz());
     //if (TMath::Abs(p->Eta())<3.0)
      // fastJet->AddParticle(px,py,pz,E);
      delete p; 
      
    }

  // =========================

   //  fastJet->RunFastJetSub(1.0);

   // cout<<fastJet->GetRap()<<" "<<fastJet->GetPhi()<<" "<<fastJet->GetArea()<<" "<<fastJet->GetPtCorr()<<endl;
 

   //  EMCAL Background ...
    // b->Fill(grid,fastJet,0);
   
   // Do FastJet analysis ...	  
    //    fastJet->DoBkgNew();
    fastJet->DoBkgStrip();
    cout<<endl;

    fastJet->RunFastJetSub(mRc,"kt",ev);   
    
   // cout<<fastJet->GetRap()<<" "<<fastJet->GetPhi()<<" "<<fastJet->GetArea()<<" "<<fastJet->GetPtCorr()<<endl;
   
    fastJet->RunFastJetSub(mRc,"AntiKt",ev);
    //cout<<fastJet->GetRap()<<" "<<fastJet->GetPhi()<<" "<<fastJet->GetArea()<<" "<<fastJet->GetPtCorr()<<endl;
    fastJet->PrintInfo();
  
   
  TObjArray *foundJets=grid->GetJetList();     
  ev->FillBkgInfo(grid);
  // ev->SetMedianPtPerArea(fastJet->GetMedianPtPerArea());	  
  ev->SetParticlePid("Simulation test Au+Au");
  ev->Fill(foundJets,grid->GetBkgPtPerCellIn(),grid->GetBkgPtPerCellOut(),grid->GetBkgPtCone(),grid->GetBkgPtCluster());      
  
  // ev->PrintJets(10);
  cout<<endl;
  
  if (mMuWrite)
    {
      outTree->Fill();
    }
 
  } 
  delete fastJet;
  delete grid;
  delete ev;	
  //gObjectTable->Print();
  

  while (nEv<nEvents);

  if (mMuWrite)
    {
      ffout->Write();
      ffout->Close();
    }
  
  //delete p;
  //cout<<endl;
  //cout<<anz<<" particles in test file."<<endl;
  
  cout<<endl;
  cout<<"Done :-)"<<endl;
  cout<<endl;
  
}
