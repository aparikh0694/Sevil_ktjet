void kt_test_pythia()
{
  gROOT->Reset();
  //gROOT->Clear();

  // Load the Event Generator abstraction library, Pythia 6
  // library, and the Pythia 6 interface library. 
  
   // ktJet lib
  if (gClassTable->GetID("ktJet") < 0) {
    cout<<"Load ktJet lib ..."<<endl;
    gSystem->Load("libKtJet.so");
  }

  if (gClassTable->GetID("TPythia6") < 0)
    {
      cout<<"Load PYTHIA6 libs ..."<<endl;
      
      gSystem->Load("libEG");
      gSystem->Load("libPythia6"); //change to your setup
      gSystem->Load("libEGPythia6");
    }

  //cout<<gSystem->GetLibraries()<<endl; 

  cout<<endl;
  cout<<" Test kt-like with seed on 'abstarct grid'"<<endl;
  cout<<" -----------------------------------------"<<
  cout<<endl;

 // create PYTHIA event ..
  TPythia6* pythia = new TPythia6();//TPythia6::Instance();//new TPythia6();
 //AliPythia* pythia=AliPythia::Instance();
  
  //cout<<"After pythia instance ..."<<endl;
  
  // set pt hard kinematics ... (check with ALICE ...)
  pythia->SetCKIN(3,100);
  pythia->SetCKIN(4,120);

  // swith off initial and final state radiation ...
  //pythia->SetMSTP(61,0);
  //pythia->SetMSTP(71,0);

  pythia->Initialize("CMS", "p", "p", 5500);
  
  // check PYTHIA interface which particles selected
  // also KF cuts !!!!!!!

  pythia->GenerateEvent();

  TClonesArray* particles = (TClonesArray*) pythia->GetListOfParticles();

  //pythia->PrintParticles();

  int nParticles=particles->GetEntries(); 

  cout<<endl;
  cout<<nParticles<<" particles generated in PYTHIA event."<<endl;
  cout<<" ---> Check particle code !!!!<<endl";
  cout<<endl;

  //cout<<pythia->GetNumberOfParticles()<<endl;
  //TDatabasePDG*  DataBase = new TDatabasePDG();

  // DEBUG:
  Int_t anz=0;
  Double_t px,py,pz,E;
  
  ktGrid *grid=new ktGrid();
  //grid->GridInfo();
  grid->SetGrid(180*0.25,0,2*TMath::Pi(),0.25*600,-10,10); // implement check on minEta, ...
  grid->GridInfo();
  grid->SetSeed(5.0);
  grid->SetBkgPtCut(0.0);
  grid->SetCone(1.0);
  grid->SetRMaxNN(1.0);
  grid->Init();

  ktGrid *grid2=new ktGrid();
  //grid->GridInfo();
  grid2->SetGrid(180*0.25,0,2*TMath::Pi(),0.25*600,-10,10); // implement check on minEta, ...
  //grid2->GridInfo();
  grid2->SetSeed(20.0);
  grid2->SetBkgPtCut(0.0);
  grid2->SetCone(1.0);
  grid2->SetRMaxNN(1.0);
  grid2->Init();

   ktGrid *grid3=new ktGrid();
  //grid->GridInfo();
  grid3->SetGrid(0.25*180,0,2*TMath::Pi(),0.25*600,-10,10); // implement check on minEta, ...
  //grid2->GridInfo();
  grid3->SetSeed(20.0);
  grid3->SetBkgPtCut(0.0);
  grid3->SetCone(1.0);
  grid3->SetRMaxNN(1.0);
  grid3->Init();

  // Fill grid with PYTHIA particles ...

  ktFastJet *fastJet=new ktFastJet();
  //ktFastJet *fastJet2=new ktFastJet();

  for (int n=0;n<nParticles;n++)
    { 
      TMCParticle *part=(TMCParticle*) particles->At(n);
  
      px=part->GetPx();
      py=part->GetPy();
      pz=part->GetPz();
      E=part->GetEnergy();

      Double_t mpt=sqrt(px*px+py*py);

      //Int_t m_PDG  = (TParticle*) part->GetPdgCode();

      if (part->GetKF()== 98)
	cout<<anz<<" "<<px<<" "<<py<<" "<<pz<<" "<<part->GetKF()<<endl;	

      if ((part->GetKF()>110 || part->GetKF()==22) && mpt>0) 
	{
	  //cout<<anz<<" "<<px<<" "<<py<<" "<<pz<<" "<<part->GetKF()<<endl;

	  TLorentzVector *p=new TLorentzVector();
	  p->SetPxPyPzE(px,py,pz,E);
	  
	  grid->Fill((TLorentzVector*) p->Clone(),true);
	  grid3->Fill((TLorentzVector*) p->Clone(),true);

	  Double_t myeta=p->Eta();

	  delete p; 
	  
	  if (TMath::Abs(myeta)<1)
	    fastJet->AddParticle(px,py,pz,E);
	  //fastJet2->AddParticle(px,py,pz,E);
	}

      anz++;
    }

  // fill background ...

  
  TF1 *eta=new TF1("eta","1",-10,10);
  TF1 *pt=new TF1("pt","1",0.1,1);
  TF1 *phi=new TF1("phi","1",-TMath::Pi(),TMath::Pi());
  
  Double_t mass=0.139;

  for (int i=0;i<10000;i++)
    {

      Double_t myeta=eta->GetRandom();
      Double_t mypt=pt->GetRandom();
      Double_t myphi=phi->GetRandom();
      
      TLorentzVector *p=new TLorentzVector();
      //p->SetPxPyPzE(px,py,pz,E);
      p->SetPtEtaPhiM(mypt,myeta,myphi,mass);

      fastJet->AddParticle(p->Px(),p->Py(),p->Pz(),p->E());

      grid2->Fill((TLorentzVector*) p->Clone());
      grid3->Fill((TLorentzVector*) p->Clone());
      delete p; 
      
      //fastJet->AddParticle(px,py,pz,E);
    }
 

  // fill background ... HIJING 

  //delete particles;
  //pythia->Delete();
  //delete pythia;

  
  //cout<<"UnLoad PYTHIA6 libs ..."<<endl;
  //gSystem->Unload("libEGPythia6"); //change to your setup
  //gSystem->Unload("libEG");
  
  /*
  if (gClassTable->GetID("THijing") < 0) {
    // Define your main function here
    cout<<"Load THIJING lib ..."<<endl;
    //gSystem->Load("libEG");
    gSystem->Load("~/lib/libTHijing.so");
  }
 
  //cout<<gSystem->GetLibraries()<<endl; 

  //gROOT->Reset();
  //gROOT->Clear();

  THijing* my_hijing = new THijing("my_hijing", "HIJING Simulation");
  //THeader          header;

  TClonesArray* particlesH = new TClonesArray("TParticle");

  //hijing->Print("B"); 
  my_hijing->SetSeed(0);
  my_hijing->SetParameter("ihpr2",12,1);
  //hijing->Initialize(energy, frame, proj, targ, projA, projZ, targA, targZ);
  //my_hijing->Initialize(5500,"CMS","A","A",197,79,197,79);
   my_hijing->Initialize(200,"CMS","A","A",197,79,197,79); 
  //header.total=
  my_hijing->GenerateEvent(0, 2, "final", particlesH);

  //cout<<"AHHH ..."<<endl;

  Int_t nBkg=particlesH->GetEntriesFast();

  cout<<endl;
  cout<<"# HIJING particles = "<<nBkg<<endl;

  for (Int_t l = 0; l < nBkg; l++)
        {
	  
          TParticle* partH = (TParticle*) particlesH->At(l);
          Double_t myeta=partH->Eta();
          Double_t mypt=partH->Pt();
          Double_t myphi=partH->Phi();
          Double_t mass=partH->GetMass();

          TLorentzVector *p=new TLorentzVector();
          p->SetPtEtaPhiM(mypt,myeta,myphi,mass);
	  
          //fastJet->AddParticle(p->Px(),p->Py(),p->Pz(),p->E());
	  
         grid2->Fill((TLorentzVector*) p->Clone());
         grid3->Fill((TLorentzVector*) p->Clone());
         delete p; 
       }
       
  
  // =========================

  */

  cout<<endl;
  cout<<anz<<" particles in event."<<endl;
  cout<<grid->GetNCells()<<" particles in grid acceptance + cuts."<<endl;
  cout<<grid2->GetNCells()<<" particles in grid acceptance + cuts incl. bkg."<<endl;
  cout<<endl;

  //fastJet->RunFastJet();
  //fastJet->RunFastJetSub();

  delete fastJet;

  TCanvas *c1=new TCanvas("c1","#1",800,600);
  c1->SetLogz();
  grid->hGrid->DrawCopy("lego2");

  //TCanvas *c2=new TCanvas("c2","#2",800,600);
  //c2->SetLogz();
  //grid2->hGrid->DrawCopy("lego2");

  TCanvas *c3=new TCanvas("c3","#3",800,600);
  c3->SetLogz();
  grid3->hGrid->DrawCopy("lego2");

  grid->GetSeeds();
  grid->DoJetfinding("coneBkg");

  grid3->GetSeeds();
  grid3->DoJetfinding("coneBkg");

  //grid->CleanGridAndCalcBkg();
 
  //cout<<" --> Average bkg.energy (pt) per cell per event (from Cone) = "<<grid->GetBkgEnergyPerCell()<<" ("<<grid->GetBkgPtPerCell()<<")"<<endl; // check pt value ???
  //cout<<"     (well, not very meaningful in this p+p test event ;-))"<<endl;

  //grid->DoJetfinding("kt");
 
  cout<<endl;
  grid->PrintJetsNice();
  grid3->PrintJetsNice();

  //grid->CleanGridAndCalcBkg();

  //cout<<" --> Average bkg.energy (pt) per cell per event (after kt-like) = "<<grid->GetBkgEnergyPerCell()<<" ("<<grid->GetBkgPtPerCell()<<")"<<endl; // check pt value ???

  //cout<<"UnLoad PYTHIA6 libs ..."<<endl;
  //gSystem->Unload("libEGPythia6"); //change to your setup
  //gSystem->Unload("libEG");

  delete grid;

  // including bkg
  /*
  grid2->GetSeeds();
  
  grid2->DoJetfinding("coneBkg");

  grid2->PrintJetsNice();

  grid2->CleanGridAndCalcBkg();
 
  cout<<" --> Average bkg.energy (pt) per cell per event (from Cone) = "<<grid2->GetBkgEnergyPerCell()<<" ("<<grid2->GetBkgPtPerCell()<<")"<<endl; // check pt value ???
  //cout<<"     (well, not very meaningful in this p+p test event ;-))"<<endl;
  */

  // including bkg + jet

  /*
  grid3->GetSeeds();
  
  grid3->DoJetfinding("coneBkg");
  
  grid3->PrintJetsNice();

  grid3->CleanGridAndCalcBkg();
  
  cout<<" --> Average bkg.energy (pt) per cell per event (from Cone) = "<<grid3->GetBkgEnergyPerCell()<<" ("<<grid3->GetBkgPtPerCell()<<")"<<endl; // check pt value ???
  //cout<<"     (well, not very meaningful in this p+p test event ;-))"<<endl;
  */

  /*
  //grid2->DoJetfinding("kt");
  cout<<endl;
  grid2->PrintJetsNice();
  grid2->CleanGridAndCalcBkg();
  cout<<" --> Average bkg.energy (pt) per cell per event (after kt-like) = "<<grid2->GetBkgEnergyPerCell()<<" ("<<grid2->GetBkgPtPerCell()<<")"<<endl; // check pt value ???
  */
  
  //gSystem->Unload("~/lib/libTHijing.so");
  
  delete pythia;

  cout<<endl;
  cout<<"Done :-)"<<endl;
  cout<<endl;

}
