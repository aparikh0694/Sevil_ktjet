// first test (Joern Putschke)

#include "ktPID.h"
#include <Riostream.h>

ClassImp(ktPID)

ktPID::ktPID()
{
  mPID=-99;

  fNsigma[0] = -999;
  fNsigma[1] = -999;
  fNsigma[2] = -999;
  fDEdx = -999;
  mass = -999;

  eTower=-999;
  iCharge=0;

  // DEBUG:
  //cout<<"Default constructor"<<endl;
}

//_________________________________________________________________________
ktPID::ktPID(ktPID *p)
{
  mPID=p->GetPID();

  eTower=p->GetETower();
  iCharge=p->GetCharge();

  fNsigma[0] = p->GetNSigmaProton();
  fNsigma[1] = p->GetNSigmaKaon();
  fNsigma[2] = p->GetNSigmaPion();

  fDEdx = p->GetdEdx();
    
  mass = p->GetMass();

  // DEBUG:
  //cout<<"Copy constructor"<<endl;
}

//_________________________________________________________________________
ktPID::~ktPID()
{
  // DEBUG:
  //cout<<"Default ktJetCell destructor"<<endl;
}

//_________________________________________________________________________

Bool_t ktPID::GetIsTower(){
  if( eTower > 0 && GetCharge()==0) return true;
  return false;
}
//_________________________________________________________________________
Bool_t ktPID::GetIsMatched(){

  if( GetIsTrack() && GetETower()>0) return true;
  return false;
}
//________________________________________________________________________
Bool_t ktPID::IsCharged(){
  if( iCharge!=0) 
    return true;
  else
    return false;
}
//_________________________________________________________________________
void ktPID::PrintPID()
{
  //cout<<endl;
  cout<<"ktPID info:"<<endl;
  cout<<"PID = "<<GetPID()<<endl;
  cout<<"isTower = "<<GetIsTower()<<endl;
  cout<<"isTrack = "<<GetIsTrack()<<endl;
  cout<<"isMatched = "<<GetIsMatched()<<endl;
  cout<<"Mass = " <<GetMass()<<endl;
  //cout<<endl; 
}
