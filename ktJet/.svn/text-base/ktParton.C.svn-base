// first test (Joern Putschke)

#include "ktParton.h"
#include <Riostream.h>

ClassImp(ktParton)

ktParton::ktParton()
{

  /*
  Parton=new TLorentzVector();//j->GetJetVector());
  Parton->SetPxPyPzE(0,0,0,0);

  PartonParticles=new TObjArray(0);
  */

   // DEBUG:
  // cout<<"Default ktParton constructor "<<endl;
}

ktParton::ktParton(TLorentzVector *v)
{

  Parton=new TLorentzVector();//j->GetJetVector());
  Parton->SetPxPyPzE(v->Px(),v->Py(),v->Pz(),v->Energy());

  PartonParticles=new TObjArray(0);

   // DEBUG:
  //cout<<"Default ktParton + Lorentzvector constructor "<<endl;
}

ktParton:: ktParton(Double_t px,Double_t py,Double_t pz, Double_t e)
{

  Parton=new TLorentzVector();
  Parton->SetPxPyPzE(px,py,pz,e);
  
  PartonParticles=new TObjArray(0);

   // DEBUG:
  //cout<<"Default ktParton + Lorentzvector constructor "<<endl;

}

void ktParton::SetParton(Double_t px,Double_t py,Double_t pz, Double_t e)
{
  
  Parton->SetPxPyPzE(px,py,pz,e);

}

void ktParton::SetParton(TLorentzVector *v)
{
  
  Parton->SetPxPyPzE(v->Px(),v->Py(),v->Pz(),v->Energy());

}

// check and think on how to maybe include it pointer to AliStack !???
void ktParton::AddParticle(TLorentzVector *v)
{
  
  if (TMath::Abs(Parton->DeltaR(*v))<1.0)
    PartonParticles->AddLast(v);

}

ktParton::~ktParton()
{

  //if (PartonParticles->GetEntriesFast()>0)
  PartonParticles->Delete();
  delete PartonParticles;
  //Parton->Delete();
  delete Parton;

   // DEBUG:
  //cout<<"Default ktParton destructor "<<endl;
}

Double_t ktParton::Phi() const
{
  Double_t mPhi;
  mPhi=Parton->Phi();if (mPhi<0) mPhi += (2*TMath::Pi());

  return mPhi;
}

void ktParton::PrintParton() 
{

  //cout<<"Jet :"<<endl;
  //cout<<"-----"<<endl;
  cout<<endl;
  cout<<"# particles = "<<PartonParticles->GetEntriesFast()<<endl;
  cout<<"Phi = "<<Phi()<<endl;
  cout<<"Eta = "<<Eta()<<endl;
  cout<<"Pt  = "<<Pt()<<endl;
  cout<<"E   = "<<E()<<endl;
  cout<<endl;
}
