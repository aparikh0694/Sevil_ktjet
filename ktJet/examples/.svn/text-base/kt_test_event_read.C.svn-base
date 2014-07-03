void kt_test_event_read(TString cName="test_fj_out.root")
{

  TChain *J=0;
  
  J=new TChain("J");

  J->Add(cName);

  cout<<"Input  = "<<cName<<endl;
  cout<<endl;
  cout<<"# rec. jet events = "<<J->GetEntries()<<endl;
  cout<<endl;
  
  Int_t Nevents=J->GetEntries();

  for (Int_t i=0;i<Nevents;i++) 
    {
      ktMuEvent *e=new ktMuEvent();
      
      J->SetBranchAddress("event",&e);
      J->GetEntry(i);

      e->PrintJets(5);
      // Trigger Info

      cout<<e->GetHTEta()<<" "<<e->GetHTPhi()<<endl;
      cout<<e->GetJPEta()<<" "<<e->GetJPPhi()<<endl;

      /*
      if (i==0)
	e->GetXiBkgSISCone()->DrawCopy();
      else
	e->GetXiBkgSISCone()->DrawCopy("same");

      e->GetXiBkgAntiKt()->DrawCopy("same");
      e->GetXiBkgKt()->DrawCopy("same");
      */

      // DEBUG:
      //cout<<endl;
      //cout<<e->GetNFJkt()<<endl;

      //ktMuFastJet *temp1=(ktMuFastJet*) e->GetFJSISCone(0); //e->GetFJkt(1);
      //temp1->PrintJet();
      
      // DEBUG:
      //for (int j=0;j<temp1->GetNJetParticles();j++)
      //cout<<j<<" "<<temp1->GetJetParticle(j)->Pt()<<endl;
      //temp1->GetJetParticle(j)->PrintParticle();

      //ktMuFastJet *temp2=(ktMuFastJet*) e->GetFJAntiKt(0); //e->GetFJkt(1);
      //temp2->PrintJet();

      /*
      ktMuJet *temp3=(ktMuJet*) e->GetConeJet(0);
      temp3->PrintJet();
      */

      cout<<endl;
   }
}
