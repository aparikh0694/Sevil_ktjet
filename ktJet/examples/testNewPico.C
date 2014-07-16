void testNewPico(Int_t nEv=10,TString mfName="test_NewPico.root")
{
  Bool_t mMuWrite=true;
  
  ktMuEvent *ev=0;
  TFile *ffout=0;
  TTree *outTree=0;

  if (mMuWrite)
    {
      ffout=new TFile(mfName,"RECREATE");
      outTree=new TTree("J","Jet event tree");
      outTree->Branch("event",&ev);
    }

  TStarJetPicoReader reader;
  reader.SetProcessV0s(kFALSE);

  //TStarJetPicoV0Cuts*    v0Cuts = reader.GetV0Cuts();
  //reader.SetV0Cuts(new TStarJetPicoV0Cuts);

  //TStarJetPicoTrackCuts* trackCuts = reader.GetTrackCuts();
  //reader.SetTrackCuts(new TStarJetPicoTrackCuts);

  //TStarJetPicoTowerCuts* towerCuts = reader.GetTowerCuts();
  //reader.SetTowerCuts(new TStarJetPicoTrackCuts);

  TStarJetPicoEventCuts* evCuts = reader.GetEventCuts();
  evCuts->SetTriggerSelection("ppJP"); //All, MB, HT, pp, ppHT, ppJP
  //evCuts->SetTriggerSelection("MB"); //All, MB, HT, pp, ppHT, ppJP
  evCuts->SetVertexZCut(30.);
  evCuts->SetRefMultCut(0);

 
  const char *sInputDir = "/data2/bruna/NewPicoDst/";///data1/putschke/picoTest";

  // build a chain with 5 files: 20-30 and 3 files: 100-160 -> total 8 files.
  TChain *chain = TStarJetPicoUtils::BuildChainFromDirectory(sInputDir, "JetTree", 1, 0);
  // TChain *chain = TStarJetPicoUtils::BuildChainFromDirectory(sInputDir, "JetTree", 3, 100, chain);
  reader.SetInputChain(chain);

  //reader.LoadDirectory(sInputDir);

  //reader.Init(1000);
  reader.Init(nEv);
  reader.PrintStatus();

  int mcount=0;

  ktStarPico *pico=new ktStarPico();
  pico->SetInfo(true);

  while(reader.NextEvent())
    {
      reader.PrintStatus(600); // show status every 10 min.

      // picoDst event info
      //reader.PrintEventInfo();

      // fill QA histograms
      //TStarJetPicoUtils::FillQAHistogramsFromPicoReader(&reader);

      TStarJetVectorContainer<TStarJetVector>* container = reader.GetOutputContainer();

      ktFastJet *fastJet=new ktFastJet();
      fastJet->SetFiducial(1,TMath::Pi());
      fastJet->SetPtCut(0.2);
      fastJet->SetInfo(false);
      fastJet->SetPrintJet(false);
      fastJet->SetFillFFArray(true);
      fastJet->SetNJetsForFF(2);
      
      ev=new ktMuEvent(mcount);
      
      if (mcount<1)
	fastJet->PrintInfo();	
   
      pico->Fill(0,fastJet,reader.GetEvent(),container);
      
      // Do FastJet analysis ...	  
      fastJet->DoBkg();
      cout<<endl;

      fastJet->RunFastJetSub(0.7,"kt",ev);
      fastJet->RunFastJetSub(0.7,"AntiKt",ev);

      fastJet->FillFF("AntiKt",ev);
      fastJet->DoXiBkg("AntiKt",ev);
      fastJet->FillFF("kt",ev);		
      fastJet->DoXiBkg("kt",ev);

      ev->SetMedianPtPerArea(fastJet->GetMedianPtPerArea());	  
      ev->SetParticlePid("New pico Test");
      ev->SetBkgPtCut(0.2);
      ev->SetRParam(0.7);
      ev->SetTriggerInfo(reader.GetEvent()->GetTrigObjs());

      ev->PrintJets(5);
      cout<<endl;

      //pico->PrintTriggerInfo();
      //ev->PrintTriggerInfo();

      //DEBUG:
      //cout<<ev->GetJPEta()<<" "<<ev->GetJPPhi()<<endl;
      //cout<<ev->GetHTEta()<<" "<<ev->GetHTPhi()<<endl;

      if (mMuWrite)
	{
	  //cout<<"Fill tree ..."<<endl;
	  outTree->Fill();
	}
      
      delete ev;
      delete fastJet;
      
      pico->Clear();
      
      //gObjectTable->Print();

      mcount++;
    }

  delete pico;

  reader.PrintStatus();

  if (mMuWrite)
    {
       ffout->Write();
       ffout->Close();
    }

  // save the QA histograms
  //TStarJetPicoQAHistograms *qa = TStarJetPicoQAHistograms::Instance();
  //qa->WriteHistogramsToFile("qaHisto.root");
}
