// first test (Joern Putschke)

#include "ktJet.h"
#include <Riostream.h>

ClassImp(ktJet)

ktJet::ktJet()
{

  JetCells=new TObjArray(0);
  JetCellsFF=new TObjArray(0);
  //JetCells->SetOwner(true);
  JetParticles=new TObjArray(0);
  JetParticlesFF=new TObjArray(0);
  //JetParticles->SetOwner(true);
  Jet=new TLorentzVector();
  Jet->SetPxPyPzE(0,0,0,0);
  type="Pseudo";
  NJetParticles=0;
  NJetParticlesFF=0;
  NCellsInEmcal=0;
  NCellsOutsideEmcal=0;
  NIterations=0;
  NSubJets=0;

  FFAreaScale=JetAreaScale=1.0;
  // DEBUG:
  //cout<<"Default ktJet constructor "<<endl;

}

// copy constructor  
ktJet::ktJet(ktJet *j)
{
  //JetCells=new TObjArray(0); //const TObjArray &(j->GetJetCells())); // copy ???
  //JetParticles=new TObjArray(0);

  FFAreaScale=JetAreaScale=1.0;

  JetCells=new TObjArray(0);
  JetCellsFF=new TObjArray(0);
  //JetCells->SetOwner(true);
  JetParticles=new TObjArray(0);
  //JetParticles->SetOwner(true);

  for (int i=0;i<j->NJetCells();i++)
    JetCells->AddLast((ktJetCell*) (j->GetJetCells()->At(i)));//->Clone());

  for (int i=0;i<j->NJetCellsFF();i++)
    JetCellsFF->AddLast((ktJetCell*) (j->GetJetCellsFF()->At(i)));//->Clone());

  Jet=new TLorentzVector();//j->GetJetVector());
  Jet->SetPxPyPzE(0,0,0,0);
  //Jet->SetPxPyPzE(j->GetJetVector()->Px(),j->GetJetVector()->Px(),j->GetJetVector()->Px(),j->GetJetVector()->E());
  //Finish();
  //cout<<"Copy ktJet constructor "<<endl;
}


void ktJet::Add(ktJet *j)
{
  // Adding for clustering algorithm !!! check also pair implementation for kt !!!

  // DEBUG:
  //j->PrintJet();
  //cout<<j->NJetCells()<<endl;

  // have to add stuff in particle list ? Or done in ktMuJet !!!??? Check !!!
  for (int i=0;i<j->NJetCells();i++)
    JetCells->AddLast((ktJetCell*) (j->GetJetCells()->At(i)));//->Clone());
  
  for (int i=0;i<j->NJetCellsFF();i++)
    JetCellsFF->AddLast((ktJetCell*) (j->GetJetCellsFF()->At(i)));//->Clone());

  NCellsInEmcal += j->NJetCellsInEmcal();
  NCellsOutsideEmcal += j->NJetCellsOutsideEmcal();

  // problem with finish .... check how to deal when jet added ... !
  TLorentzVector *jv=(TLorentzVector*) j->GetJetVector();
  *Jet += *jv;

  for (int i=0;i<j->GetJetCells()->GetEntriesFast();i++)
    {
      ktJetCell *myCell=(ktJetCell*) j->GetJetCells()->At(i);
      for (int k=0;k<myCell->NParticles();k++)
	{
	  JetParticles->AddLast((TLorentzVector*) (myCell->GetParticleList()->At(k)));
	}
      NJetParticles += myCell->NParticles();
    }

  for (int i=0;i<j->GetJetCellsFF()->GetEntriesFast();i++)
    {
      ktJetCell *myCell=(ktJetCell*) j->GetJetCellsFF()->At(i);

      for (int k=0;k<myCell->NFFParticles();k++)
	{	  
	  JetParticlesFF->AddLast((TLorentzVector*) (myCell->GetFFParticleList()->At(k)));
	}
      NJetParticlesFF += myCell->NFFParticles(); 
    }
}

void ktJet::ClearJet()
{
  JetCells->Clear();
  JetCellsFF->Clear();
  Jet->SetPxPyPzE(0,0,0,0);
}

void ktJet::ResetJet()
{
  //JetCells->Clear();
  Jet->SetPxPyPzE(0,0,0,0);
}

ktJet::~ktJet()
{

  //cout<<type<<endl;
  //cout<<"Jet cells part ..."<<endl;

  //JetCells->Delete();
  delete JetCells;
  delete JetCellsFF;

  //cout<<"Jet part ..."<<endl;

  //JetParticles->Delete();
  //cout<<"->Delete() done"<<endl;
  delete JetParticles;
  //cout<<"FF now ..."<<endl;
  //cout<<JetParticlesFF<<endl;
  
  if (JetParticlesFF!=0)
    {
      //cout<<JetParticlesFF->GetEntries()<<endl;
      //cout<<"delete ..."<<endl;
      delete JetParticlesFF;
    }

  //cout<<"JetParticles done."<<endl;

  delete Jet;
  //delete hJet;

  // DEBUG:
  //cout<<"Default ktJet destructor"<<endl;

}

void ktJet::AddCellFF(ktJetCell *JetCell)
{
  // check !!! just to get HIJING FF !!!????
  JetCell->SetIsFF(true);
  JetCellsFF->AddLast((ktJetCell*) JetCell);
}

void ktJet::AddCellBkg(ktJetCell *JetCell)
{
  JetCells->AddLast((ktJetCell*) JetCell);
}

void ktJet::AddCell(ktJetCell *JetCell)
{

  if (!JetCell->InJet()) // !!!
    {
      JetCell->SetInJet(true);
      JetCells->AddLast((ktJetCell*) JetCell);//->Clone());
    }
  else
    {
      //DEBUG:
      cout<<"Warning : Cell already in Jet"<<endl;
    }

}

void ktJet::AddCell(ktJetCell *JetCell,Bool_t InEmcal)
{

  if (!JetCell->InJet()) // !!!
    {
      JetCell->SetInJet(true);
      JetCells->AddLast((ktJetCell*) JetCell);//->Clone());
      if (InEmcal)
	NCellsInEmcal++;
      else
	NCellsOutsideEmcal++;
    }
  else
    {
      //DEBUG:
      cout<<"Warning : Cell already in Jet"<<endl;
    }

}

void ktJet::AddCell(ktJetCell *JetCell,Bool_t InEmcal,Bool_t mSharing)
{

  if (!JetCell->InJet() || mSharing) // !!!
    {
      JetCell->SetInJet(true);
      JetCells->AddLast((ktJetCell*) JetCell);//->Clone());
      if (InEmcal)
	NCellsInEmcal++;
      else
	NCellsOutsideEmcal++;
    }
  else
    {
      //DEBUG:
      cout<<"Warning : Cell already in Jet"<<endl;
    }

}

Double_t ktJet::Phi() const
{
  Double_t mPhi;
  mPhi=Jet->Phi();if (mPhi<0) mPhi += (2*TMath::Pi());

  return mPhi;
}

//void ktJet::Finish(Int_t Nphi, Double_t phiMin, Double_t phiMax, Int_t Neta, Double_t  etaMin, Double_t etaMax)
void ktJet::Finish()
{
  // Dummy

  //hJet=new TH2F();

   NJetParticles=0;
   Jet->SetPxPyPzE(0,0,0,0);

   for (int i=0;i<JetCells->GetEntriesFast();i++)
    {
      ktJetCell *myCell=(ktJetCell*) JetCells->At(i);
      *Jet += *(myCell->GetCellParticle());
      NJetParticles += myCell->NParticles();
      //*Jet += *mPart;
    }
}

void ktJet::CleanCellsInJet()
{
  // clean InJet flag for bkg summing

  for (int i=0;i<JetCells->GetEntries();i++)
    {
      ktJetCell *myCell=(ktJetCell*) JetCells->At(i);
      myCell->SetInJet(false);
    }
}

void ktJet::SetCellsInJet()
{
  // clean InJet flag for bkg summing

  for (int i=0;i<JetCells->GetEntries();i++)
    {
      ktJetCell *myCell=(ktJetCell*) JetCells->At(i);
      myCell->SetInJet(true);
    }
}

void ktJet::SetInCluster()
{
   for (int i=0;i<JetCells->GetEntries();i++)
    {
      ktJetCell *myCell=(ktJetCell*) JetCells->At(i);
      myCell->SetInCluster(true);
    }
}

Double_t ktJet::GetDistance(ktJet* mDjet)
{
  Double_t md=0;

  Double_t mdphi=Phi()-mDjet->Phi();
  Double_t mdeta=Eta()-mDjet->Eta();

  md=TMath::Sqrt(mdphi*mdphi+mdeta*mdeta);
  
  return md;
}


Double_t ktJet::GetDistance(TLorentzVector* mDjet)
{
  Double_t md=0;
  Double_t mPhi=0;
  
  mPhi=mDjet->Phi();if (mPhi<0) mPhi += (2*TMath::Pi());

  Double_t mdphi=Phi()-mPhi;
  Double_t mdeta=Eta()-mDjet->Eta();

  md=TMath::Sqrt(mdphi*mdphi+mdeta*mdeta);
  
  return md;
}

/*
TH1D* ktJet::DrawXi()
{
  // only DEBUG: take pt=100; // memory leak ...

  Double_t xiPt=100.0;
   
  TH1D *hJetXi=new TH1D("hJetXi","Xi dist. of Background per event",50,0,10);
  
  for (int i=0;i<JetCellsFF->GetEntriesFast();i++)
    {
      ktJetCell *myCell=(ktJetCell*) JetCellsFF->At(i);

      for (int k=0;k<myCell->NFFParticles();k++)
	{
	  TLorentzVector *lvec=(TLorentzVector*) (myCell->GetFFParticleList()->At(k));
	  Double_t mXi=TMath::Log((Double_t) xiPt/(Double_t) lvec->Pt());
	  hJetXi->Fill(mXi);
	}
    }
  
  return hJetXi;
}
*/

// inherit ktGrid (could be smater ;-))
void ktJet::Finish(Int_t Nphi, Double_t phiMin, Double_t phiMax, Int_t Neta, Double_t  etaMin, Double_t etaMax,Int_t myPhiBin, Int_t myEtaBin)
//void ktJet::Finish()
{

  // DEBUG:
  //cout<<"   Finish Jet up ... calculate properties ..."<<endl;

  NJetParticles=0;
  NJetParticlesFF=0;

  // Jet 2d histogram
  TString hName="hJet_";
  hName+=type;
  hName+="_";
  hName+=myPhiBin;
  hName+="_";
  hName+=myEtaBin;

  // DEBUG
  //cout<<hName<<endl;

  //this->SetName(hname);

  //hJet=new TH2F(hName,"Jet eta/phi (pt weighted)",Nphi,phiMin,phiMax,Neta,etaMin,etaMax);
  //hJet->SetMinimum(0.1);

  for (int i=0;i<JetCells->GetEntriesFast();i++)
    {
      ktJetCell *myCell=(ktJetCell*) JetCells->At(i);
      *Jet += *(myCell->GetCellParticle());
      //hJet->Fill(myCell->Phi(),myCell->Eta(),myCell->Pt());

      // JetParticle list

      //TObjArray *myCellPart=(TObjArray*) myCell->GetParticleList();
      
      // DEBUG:
      //cout<<myCell->NParticles()<<" "<<myCellPart->GetEntriesFast()<<endl;
      //cout<<JetParticles->GetEntriesFast()<<endl;

      for (int k=0;k<myCell->NParticles();k++)
	{
	  JetParticles->AddLast((TLorentzVector*) (myCell->GetParticleList()->At(k)));//->Clone());
	}

      NJetParticles += myCell->NParticles();

      //delete myCellPart;
      //*Jet += *mPart;
    }

  for (int i=0;i<JetCellsFF->GetEntriesFast();i++)
    {
      ktJetCell *myCell=(ktJetCell*) JetCellsFF->At(i);

      for (int k=0;k<myCell->NFFParticles();k++)
	{	  
	  //DEBUG:
	  //TLorentzVector *test=(TLorentzVector*) (myCell->GetFFParticleList()->At(k));	  
	  //cout<<test->Pt()<<" "<<test->DeltaR(*GetJetVector())<<endl;
      
	  JetParticlesFF->AddLast((TLorentzVector*) (myCell->GetFFParticleList()->At(k)));//->Clone());
	}

      NJetParticlesFF += myCell->NFFParticles(); 
    }

  // DEBUG:
  //PrintJet();
}

void ktJet::PrintJet() const
{

  //cout<<"Jet :"<<endl;
  //cout<<"-----"<<endl;
  cout<<endl;
  cout<<"# cells      = "<<JetCells->GetEntriesFast()<<endl;
  cout<<"# cells in EMCAL = "<<NCellsInEmcal<<endl;
  cout<<"# cells outside EMCAL = "<<NCellsOutsideEmcal<<endl;
  cout<<"# FF cells = "<<JetCellsFF->GetEntriesFast()<<endl;
  cout<<"# particles = "<<JetParticles->GetEntriesFast()<<endl;//NJetParticles<<endl;
  cout<<"# FF particles = "<<JetParticlesFF->GetEntriesFast()<<endl;//NJetParticles<<endl;
  cout<<"Type = "<<type<<endl;
  cout<<"Phi = "<<Phi()<<endl;
  cout<<"Eta = "<<Eta()<<endl;
  cout<<"Pt  = "<<Pt()<<endl;
  cout<<"E   = "<<E()<<endl;
  cout<<"Jet Scale = "<<JetAreaScale<<endl;
  cout<<"FF Scale = "<< FFAreaScale<<endl;
  cout<<endl;
}
