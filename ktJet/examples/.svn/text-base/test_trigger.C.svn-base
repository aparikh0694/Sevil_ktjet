void test_trigger()
{

  gROOT->Reset();
  gROOT->Clear();

  // ktJet lib
  if (gClassTable->GetID("ktJet") < 0) {
    cout<<"Load ktJet lib ..."<<endl;
    gSystem->Load("libKtJet.so");
  }

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
  
  ktTrigger *grid=new ktTrigger(3,3,2,2,10);
  //grid->GridInfo();
  grid->SetGrid(0.25*180,0,2*TMath::Pi(),0.25*600,-10,10); // implement check on minEta, ...
  grid->GridInfo();
  grid->SetSeed(5.0);
  grid->SetBkgPtCut(0.0);
  grid->SetCone(0.8);
  grid->SetRMaxNN(0.8);
  grid->Init();

  grid->PrintTriggerInfo();
  
  int anz=0;

   while(!feof(f1))
    {

      dat1 = fscanf(f1,"%f %f %f %f",&px,&py,&pz,&E);

      TLorentzVector *p=new TLorentzVector();
      p->SetPxPyPzE(px,py,pz,E);

      grid->Fill((TLorentzVector*) p->Clone(),true);

      anz++;

      delete p; 
    }

   cout<<anz<<" particles in event."<<endl;
   cout<<grid->GetNCells()<<" particles in grid acceptance + cuts."<<endl;

   grid->RunTrigger();

   grid->PrintTriggerPatches()<<endl;

   cout<<"Number of trigger patches = "<<grid->GetNTriggerPatches()<<endl;
   cout<<"pt of trigger patch = "<<grid->GetTriggerPt()<<endl;
   cout<<"eta = "<<grid->GetTriggerEta()<<" | phi = "<<grid->GetTriggerPhi()<<endl;
   cout<<endl;

   TCanvas *c1=new TCanvas("c1","#1",800,600);
   c1->SetLogz();
   grid->hGrid->DrawClone("lego2");

   delete grid;
}

