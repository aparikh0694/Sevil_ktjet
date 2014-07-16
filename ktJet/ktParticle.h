#ifndef ROOT_ktParticle
#define ROOT_ktParticle

#include "TBuffer.h"
#include "TLorentzVector.h"
#include "ktPID.h"

class ktParticle : public TLorentzVector
{

 protected:

  ktPID* Pid;
  Bool_t inJet;
  
 public:

  ktParticle();
  ktParticle(ktParticle *p); //copy constructor !
  ktParticle(TLorentzVector *p, ktPID mPid=0, Bool_t mIsInJet=0);
  virtual ~ktParticle();
 
  void PrintParticle();
  
  ktPID* GetPid() {return Pid;}
  Bool_t IsInJet() {return inJet;}
  Bool_t IsCharged() {return Pid->IsCharged();}

  void SetPid(ktPID *mPid) {Pid= new ktPID(mPid);}
  void SetIsInJet(Bool_t m_inJet) {inJet=m_inJet;}

  Double_t GetPOverE();

  // static int memCounter;

  ClassDef(ktParticle,2)
};

#endif
