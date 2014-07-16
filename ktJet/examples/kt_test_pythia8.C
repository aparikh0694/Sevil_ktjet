void kt_test_pythia8()
{

      ktPythia8 *p=new ktPythia8("/Users/salur/software/pythia8130/xmldoc/");
  // choose path to your pythia8 xmldoc installation directory !

  p->SetRandomSeed();
  p->InitPhaseSpace("28.0","32.0");
  p->Init("RHIC");

  //for (int i=0;i<10;i++) {

  ktGrid *grid=new ktGrid();
  grid->SetGrid(180*0.25,0,2*TMath::Pi(),0.25*600,-10,10); // implement check on minEta, ...
  grid->GridInfo();
  grid->SetSeed(5.0);
  grid->SetBkgPtCut(0.0);
  grid->SetCone(1.0);
  grid->SetRMaxNN(1.0);
  grid->Init();

  ktGrid *grid2=new ktGrid();
  grid2->SetGrid(180*0.25,0,2*TMath::Pi(),0.25*600,-10,10); // implement check on minEta, ...
  grid2->GridInfo();
  grid2->SetSeed(5.0);
  grid2->SetBkgPtCut(0.0);
  grid2->SetCone(1.0);
  grid2->SetRMaxNN(1.0);
  grid2->Init();

  // get Pythia8 event ...
  p->Run();
  p->EventInfo();
  //p->ParticleList();
  
  // fill grid with Pythia8 event ...
  // check and compare with ALICE sims (decay on/off !?)
  p->Fill(grid,"EMCAL");
  p->Fill(grid2,"ALL");

  cout<<endl;
  cout<<(Int_t) grid->GetNCells()<<" cells in grid acceptance + cuts."<<endl;
  cout<<(Int_t) grid2->GetNCells()<<" cells in grid acceptance + cuts."<<endl;

  TCanvas *c1=new TCanvas("c1","#1",800,600);
  c1->SetLogz();
  grid->hGrid->DrawCopy("lego2");

  grid->GetSeeds();
  grid->DoJetfinding("coneBkg");
  grid2->GetSeeds();
  grid2->DoJetfinding("coneBkg");

  // check hard scattering values eta/phi from PYTHIA ...
  p->HardProcessInfo();

  grid->PrintJetsNice();
  grid2->PrintJetsNice();

  delete grid;
  delete grid2;
  //}

  //p->Statistics();

  delete p;
  
}
