

#include "ktParticle.h"
#include <Riostream.h>

ClassImp(ktParticle)

ktParticle::ktParticle()
{

  Pid=0;
  inJet=false;

  // DEBUG:
  //cout<<"Default ktParticle constructor "<<endl;
}

ktParticle::ktParticle(ktParticle *p)
{
  // DEBUG:
  //p->PrintParticle();

  Pid= new ktPID(p->GetPid());
  inJet=p->IsInJet();

  SetPxPyPzE(p->Px(),p->Py(),p->Pz(),p->E());
  
}


ktParticle::ktParticle(TLorentzVector *p, ktPID mPid, Bool_t mIsInJet)
{
  Pid= new ktPID(mPid);
  inJet=mIsInJet;

  SetPxPyPzE(p->Px(),p->Py(),p->Pz(),p->E()); 
}


ktParticle::~ktParticle()
{

  delete Pid;
  // DEBUG:
  //cout<<"Default ktParticle destructor"<<endl;
}


Double_t ktParticle::GetPOverE(){

  if( !(GetPid()->GetIsMatched())) return -999;
  else return P()/GetPid()->GetETower();
  
}

void ktParticle::PrintParticle()
{
  cout<<"ktParticle information:"<<endl;
  cout<<"-----------------------"<<endl;
  cout<<"eta = "<<Eta()<<endl;
  cout<<"phi = "<<Phi()<<endl;
  cout<<"pt = "<<Pt()<<endl;
  if( Pid){
    cout<<"Pid = "<< endl;
    GetPid()->PrintPID();
  }
  cout<<"InJet = "<<IsInJet()<<endl;
  

}
