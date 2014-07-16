void kt_test_trigger()
{

  // RHIC central background ... (v2=0.0)
  ktMCBkg *b=new ktMCBkg("RHIC","Au+Au",0.0);
  
  // Pythia8 mono-jet
  ktPythia8 *p=new ktPythia8("/home/putschke/pythia8100/xmldoc");
  // choose path to your pythia8 xmldoc installation directory !
  //p->SetRandomSeed();
  p->InitPhaseSpace("28.0","32.0");
  p->Init("RHIC");

  ktTrigger *grid=new ktTrigger(2,2,1,1,5);
  grid->SetGrid(30*1,0,2*TMath::Pi(),10*1,-1,1); 
  grid->GridInfo();
  grid->SetSeed(5.0);
  grid->SetBkgPtCut(0.0);
  grid->SetCone(0.7);
  grid->SetRMaxNN(1.0);
  grid->Init();

  grid->PrintTriggerInfo();

  // get Pythia8 event ...
  p->RunWithPycell();
  //p->EventInfo();

  p->Fill(grid,"NEUTRAL");
  b->FillNeutral(grid);

  cout<<endl;
  cout<<(int) grid->GetNCells()<<" particles in grid acceptance + cuts."<<endl;
  
  cout<<endl;
  //p->HardProcessInfo();
  cout<<"Highest Pythia jet (PYCELL):"<<endl;
  cout<<"eta = "<<p->GetPycellEta()<<" | phi = "<<p->GetPycellPhi2()<<" | pt = "<<p->GetPycellPt()<<endl;

  grid->RunTrigger();
  
  //grid->PrintTriggerPatches();
  
  cout<<"Number of trigger patches = "<<grid->GetNTriggerPatches()<<endl;
  if (grid->GetNTriggerPatches()>0)
    {
      cout<<"pt of trigger patch = "<<grid->GetTriggerPt()<<endl;
      cout<<"eta = "<<grid->GetTriggerEta()<<" | phi = "<<grid->GetTriggerPhi()<<endl;
    }
  
  cout<<endl;

  TCanvas *c1=new TCanvas("c1","#1",800,600);
  //c1->SetLogz();
  grid->hGrid->DrawCopy("lego2");
  
  delete grid;
  delete p;
  delete b;
}
