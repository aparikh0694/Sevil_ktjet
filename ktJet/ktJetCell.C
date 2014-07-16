// first test (Joern Putschke)

#include "ktJetCell.h"
#include <Riostream.h>

ClassImp(ktJetCell)

ktJetCell::ktJetCell()
{
  ParticleList=new TObjArray(0);
  ParticleList->SetOwner(kTRUE);

  FFParticleList=new TObjArray(0);
  FFParticleList->SetOwner(kTRUE);

  PIDParticleList=new TObjArray(0);
  PIDParticleList->SetOwner(kTRUE);

  PIDFFParticleList=new TObjArray(0);
  PIDFFParticleList->SetOwner(kTRUE);

  etaBin=phiBin=-10;
  //ParticleList=new TClonesArray("TLorentzVector",10);
  //TClonesArray ParticleList("TLorentzVector",10);
  CellParticle=new TLorentzVector();
  CellParticle->SetPxPyPzE(0,0,0,0);
  CellParticleFF=new TLorentzVector();
  CellParticleFF->SetPxPyPzE(0,0,0,0);
  //kt2=0;
  NNIndex=0;
  inJet=false;
  isTrack=true;
  isNN=false;
  isFF=false;
  inEmcal=false;
  inCluster=false;
  // DEBUG:
  //cout<<"Default constructor"<<endl;
 
}

ktJetCell::~ktJetCell()
{
 
  ParticleList->Delete();
  delete ParticleList;
  
  //cout<<"Delete FF list ..."<<endl;

  FFParticleList->Delete();
  delete FFParticleList;
 
  PIDParticleList->Delete();
  delete PIDParticleList;

  PIDFFParticleList->Delete();
  delete PIDFFParticleList;

  //cout<<FFParticleList->GetEntries()<<endl;
  /*
  if (FFParticleList->GetEntries()>0)
    {
      //cout<<FFParticleList->GetEntries()<<endl;
      //FFParticleList->Delete();
    }

  delete FFParticleList;
  */

  //cout<<"Cell list ..."<<endl;

  delete CellParticle;
  delete CellParticleFF;

  // DEBUG:
  //cout<<"Default ktJetCell destructor"<<endl;

}

void ktJetCell::PrintCell()
{
  cout<<"Eta/Phi Bin = "<<GetEtaBin()<<" / "<<GetPhiBin()<<endl;
  cout<<"Eta = "<<Eta()<<endl;
  cout<<"Phi = "<<Phi()<<endl;
  cout<<"E   = "<<E()<<endl;
  cout<<"Pt  = "<<Pt()<<endl;
  cout<<"Kt2 = "<<Kt2()<<endl;
  cout<<"InJet = "<<InJet()<<endl;
  cout<<"InEmcal = "<<InEmcal()<<endl;
}

Double_t ktJetCell::Phi()
{
  Double_t mPhi;
  mPhi=CellParticle->Phi();if (mPhi<0) mPhi += (2*TMath::Pi());

  return mPhi;
}

Double_t ktJetCell::Eta()
{
  Double_t mEta;
  mEta=CellParticle->PseudoRapidity();
  
  return mEta;
}

Double_t ktJetCell::PhiFF()
{
  Double_t mPhi;
  mPhi=CellParticleFF->Phi();if (mPhi<0) mPhi += (2*TMath::Pi());

  return mPhi;
}

void ktJetCell::AddParticle(TLorentzVector *mPart)
{
  //cout<<NParticles()<<endl;
  //new(ParticleList[0]) TLorentzVector(*mPart);
  ParticleList->AddLast(mPart);
  *CellParticle += *mPart;
  
  inJet=false;

  // DEBUG:
  //cout<<" = "<<CellParticle->E()<<endl;
}

void ktJetCell::AddFFParticle(TLorentzVector *mPart)
{
  //DEBUG:
  //if (mPart->Pt()<1)
  //cout<<mPart->Pt()<<endl;

  FFParticleList->AddLast(mPart);
  *CellParticleFF += *mPart;
}


void ktJetCell::AddParticlePID(ktPID *mPID)
{
  ktPID *mCopyPID=new ktPID(mPID);
  PIDParticleList->AddLast(mCopyPID);
  
  // DEBUG:
  /*
  if (mPID->GetPID()==11)
    {
      cout<<"AddParticlePID"<<endl;
      mPID->PrintPID();
    }
  */
}

void ktJetCell::AddFFParticlePID(ktPID *mPID)
{
  // Check Memory !!!
  
  ktPID *mCopyPID=new ktPID(mPID);
  // mCopyPID->PrintPID();
  PIDFFParticleList->AddLast(mCopyPID);
  // DEBUG:
  //mPID->PrintPID();
}

void ktJetCell::Clean()
{
  ParticleList->Clear();
}
