void kt_test_ana()
{

  TChain *J=0;
  // TString cName="test_NewPico.root";
  TString cName="test_fj_out.root";
    
  //J=new TChain("J");//SimuFJ");
  J=new TChain("JSimuFJ");

  J->Add(cName);

  TFile *f=new TFile("test_ana.root","RECREATE");

  cout<<"Input  = "<<cName<<endl;
  cout<<endl;
  cout<<"# rec. jet events = "<<J->GetEntries()<<endl;
  cout<<endl;
  
  Int_t Nevents=J->GetEntries();

  ktAna *mAna=new ktAna();
  mAna->SetJetAlgo("AntiKt");
  mAna->SetVerbose(true);
  //mAna->InitHistograms();
  // mAna->SetFFPtRange(5,50);
  mAna->SetDiJetPtRange(5,500);
  //mAna->SetDoBkgCorrections(false);
  mAna->InitHistograms();
  /*
  ktAna *mAna3=new ktAna();
  mAna3->SetJetAlgo("kt");
  mAna3->SetFFPtRange(5,50);
  mAna3->SetDiJetPtRange(5,50);
  mAna3->SetVerbose(false);
  mAna3->SetDoBkgCorrections(false);
  mAna3->InitHistograms();

  ktAna *mAna4=new ktAna();
  mAna4->SetJetAlgo("Cone");
  mAna4->SetVerbose(true);
  mAna4->InitHistograms();
  */

  mAna->PrintAnalysisCuts();

  for (Int_t i=0;i<Nevents;i++) 
    {
      ktMuEvent *e=new ktMuEvent();
      
      J->SetBranchAddress("event",&e);
      J->GetEntry(i);

      //cout<<i<<endl;
      int Nj=e->GetNFJAntiKt();
      cout<<"  Nj  "<<Nj<<endl;
      

      cout<<" median pt  "<<e->GetMedianPtPerArea()<<endl;
      
      for(Int_t cj=0; cj<Nj;cj++){
	ktMuFastJet* ktjet;
	ktjet=e->GetFJAntiKt(cj);
	cout<<" area "<<ktjet->Area()<<" pt  "<<ktjet->Pt()<< "  eta  "<<ktjet->Eta()<<"  phi  "<<ktjet->Phi()<<"  rho  "<<ktjet->MedianPtPerArea()<<endl;

      }
      
      
      //mAna->AnalyzeEvent(e);
      //mAna->AnalyzeEventTrig(e);
      //mAna3->AnalyzeEvent(e);

      //mAna->AnalyzeEvent(e);
      //mAna4->AnalyzeEvent(e);

      //mAna2->CompareEvent(mAna3);     

      delete e;
   }

  mAna->Finish();
  //mAna2->Finish();
  //mAna3->Finish();
  //mAna4->Finish();
  
  // Example ...
  /*
  (TH1D*) (mAna->GetXi())->DrawClone();
  (TH1D*) (mAna->GetXiBkg(2,2,1))->DrawClone("same");
  (TH1D*) (mAna->GetXiSub(23,2,1))->DrawClone("same");
  */

  //(TH1D*) (mAna2->GetZBkg())->DrawClone();
  //(TH1D*) (mAna2->GetXiBkg())->DrawClone();

  //mAna->DrawPlots();
  //mAna->DrawPlots();
  //  mAna3->DrawPlots();
  //mAna2->DrawFFPlots();
  //mAna3->DrawFFPlots();

  //mAna2->DrawFFPlotsTrig();
  //mAna3->DrawFFPlotsTrig();
  //mAna4->DrawPlots();

  f->cd();

  
  mAna->Write("AntiKtAna");
  // mAna3->Write("ktAna");
  //mAna4->Write("LOConeAna");

  f->Write();

  delete mAna;
  //delete mAna2;
  //delete mAna3;
  //delete mAna4;

  //f->Close();
}
