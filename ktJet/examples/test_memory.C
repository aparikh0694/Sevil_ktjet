void test_memory()
{
  gROOT->Reset();
  gROOT->Clear();

  // ktJet lib
  if (gClassTable->GetID("ktJet") < 0) {
    cout<<"Load ktJet lib ..."<<endl;
    gSystem->Load("libKtJet.so");
  }

  gObjectTable->Print();

  cout<<endl;
  cout<<" Test kt-like with seed on 'abstarct grid'"<<endl;
  cout<<" -----------------------------------------"<<endl;
  cout<<"       (Joern Putschke, June 2006)"<<endl;
  cout<<endl;

  // handle input asci file (fastjet ...)

  //FILE *f1=fopen("single-event.dat","r");
  FILE *f1=fopen("Pythia-dijet-ptmin100-lhc-pileup-1ev.dat","r");
  Float_t px,py,pz,E;
  Float_t eta,phi;
  Int_t dat1;

  // DEBUG:
  Int_t anz=0;
  Double_t esum=0;
  //ktJetCell *g=new ktJetCell();
  //cout<<g->test()<<endl;
  //delete g;
  
  ktGrid *grid=new ktGrid();
  //grid->GridInfo();
  grid->SetGrid(0.25*180,0,2*TMath::Pi(),0.25*600,-10,10); // implement check on minEta, ...
  grid->GridInfo();
  grid->SetSeed(5.0);
  grid->SetBkgPtCut(0.0);
  grid->SetCone(0.8);
  grid->SetBkgCone(0.4);
  grid->SetNBkgSeeds(1000);
  grid->SetRMaxNN(1.4);
  //grid->SetNMaxIterations(0); // test to comapre with bkg case !!! Works !
  grid->Init();

  while(!feof(f1))
    {

      dat1 = fscanf(f1,"%f %f %f %f",&px,&py,&pz,&E);

      TLorentzVector *p=new TLorentzVector();
      p->SetPxPyPzE(px,py,pz,E);

      grid->Fill((TLorentzVector*) p->Clone(),true);
     
      delete p; 

      anz++;
    }

  /*
  TCanvas *c1=new TCanvas("c1","#1",800,600);
  c1->SetLogz();
  grid->hGrid->DrawClone("lego2");
  */

  cout<<anz<<" particles in event."<<endl;
  cout<<grid->GetNCells()<<" particles in grid acceptance + cuts."<<endl;

  grid->GetSeeds();
  //grid->PrintSeeds();
  grid->DoJetfinding("coneBkg");

  grid->CalcBkg();
  grid->BkgRandomCones();

  //grid->PrintJetsNiceBkgSub();

  // cout<<" --> Average bkg.energy per cell per event = "<<grid->GetBkgPtPerCell()<<endl;

  //grid->CalcBkgNJetsRemove(2);
  //cout<<" --> Average bkg.energy per cell per event (remove N highest jets) = "<<grid->GetBkgPtPerCell()<<endl;

  //grid->DoJetfinding("kt");

  cout<<endl;
  grid->PrintJetsNiceBkgSub();

  // PrintJetsNiceBkgSub random Jet ...
  //cout<<grid->GetSeedJet(0,0)->Pt()<<endl;
  //grid->GetSeedJet(0,0)->PrintJet();

  grid->CleanGridAndCalcBkg();

  /*
  ktJet *seedJet=new ktJet();
 
  //grid->GetSeedJet(seedJet,-0.64104159,0.39978924);
  //grid->GetSeedJet(seedJet,-0.476018,0.331965);
  grid->GetSeedJet(seedJet, 2.81747326, 3.46968714);
  seedJet->PrintJet();

  delete seedJet;
  */

  grid->DoJetfinding("cone");

  grid->PrintJetsNiceBkgSub();

  //grid->CleanGridAndCalcBkg();

  //cout<<" --> Average bkg.energy (pt) per cell per event (after kt-like) = "<<grid->GetBkgEnergyPerCell()<<" ("<<grid->GetBkgPtPerCell()<<")"<<endl; // check pt value ???

  //delete c1;
  delete grid;

  gObjectTable->Print();

  //cout<<endl;
  //cout<<"Done :-)"<<endl;
  //cout<<endl;

}
