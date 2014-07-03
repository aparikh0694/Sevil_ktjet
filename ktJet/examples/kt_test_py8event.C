void kt_test_py8event(Int_t nEvents=100)
{

  TFile *f=new TFile("test_py8event.root","RECREATE");

  ktPy8Event *myEv=0;
  
  TTree *outTree=new TTree("PYTHIA8","PYTHIA8 simulation tree");
  outTree->Branch("event",&myEv);

  ktPythia8 *p=new ktPythia8("/Users/putschke/pythia8100/xmldoc");
  // choose path to your pythia8 xmldoc installation directory !
  
  //p->SetRandomSeed();

  // RHIC:
  p->InitPhaseSpace("28.0","32.0");
  //p->InitPhaseSpace("15.0","75.0");
  p->Init("RHIC");

  // LHC:
  //p->InitPhaseSpace("98.0","102.0");
  //p->Init("LHC");

  for (int i=0;i<nEvents;i++) 
    {
      
      // get Pythia8 event ...
      p->RunWithPycell();
      //p->EventInfo();
      //p->ParticleList();
      //p->HardProcessInfo();


      myEv=new ktPy8Event(i);      
      p->FillEvent(myEv);

      outTree->Fill();
      
      // DEBUG:
      if (fmod(i,(Int_t) nEvents/10)==0 && i>0)
	{
	  cout<<"Event number : "<<i<<endl;
	  myEv->EventInfo();
	  cout<<endl;
	}

      delete myEv;
      
    }
  
  
  p->Statistics();
  cout<<endl;
  cout<<"sigma = "<<p->GetSigmaGen()<<endl;
  //Double_t mySigma=p->GetSigmaGen();
  //mySigma.Write("sigma");
  cout<<endl;

  f->Write();
  f->Close();
  
  delete p;
  
}
