// first test (Joern Putschke)

#ifndef ROOT_ktPID
#define ROOT_ktPID

#include "TObject.h"
#include "TMath.h"


class ktPID : public TObject 
{

 protected:
 
 Int_t mPID; // Definition:PYTHIA coding
 Float_t eTower;
 Int_t iCharge;

 Float_t fNsigma[3];
 Float_t fDEdx;

 Float_t mass;

 // if needed one can add more
 // experiment specific informations here ..
 // has to update than ktGrid->Fill(...)

 public:
 
 ktPID();
 ktPID(ktPID *p); // copy constructor

 virtual ~ktPID();

 void PrintPID();

 void SetPID(Int_t m_PID) {mPID=m_PID;}
 void SetETower(Float_t mETower) {eTower=mETower;}
 void SetCharge(Int_t Charge) {iCharge=Charge;}
 void SetdEdx(Float_t dEdx){fDEdx=dEdx;}
 void SetNSigma(Float_t *nsigma){fNsigma[0]=nsigma[0];fNsigma[1]=nsigma[1];fNsigma[2]=nsigma[2];}
 void SetNSigmaProton(Float_t nsigma){ fNsigma[0]=nsigma;}
 void SetNSigmaKaon(Float_t nsigma){fNsigma[1]=nsigma;}
 void SetNSigmaPion(Float_t nsigma){fNsigma[2]=nsigma;}
 void SetMass(Float_t value){mass = value;}

 Int_t GetPID() {return mPID;}
 Bool_t GetIsTrack() {return iCharge;}
 Bool_t GetIsTower(); 
 Bool_t GetIsMatched();
 Bool_t GetIsNeutral() {return iCharge;}
 Int_t GetCharge() {return iCharge;} 
 Bool_t IsCharged();

 Float_t GetdEdx(){return fDEdx;}
 void GetNSigma( Float_t *nsigma){nsigma[0]=fNsigma[0];nsigma[1]=fNsigma[1];nsigma[2]=fNsigma[2];}
 Float_t GetNSigmaProton(){return fNsigma[0];}
 Float_t GetNSigmaKaon(){return fNsigma[1];}
 Float_t GetNSigmaPion(){return fNsigma[2];}
 Float_t GetETower(){return eTower;}

 Float_t GetMass(){return mass;}

 ClassDef(ktPID,1)
};

#endif
