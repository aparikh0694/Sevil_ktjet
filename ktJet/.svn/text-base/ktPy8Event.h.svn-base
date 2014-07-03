// first test (Joern Putschke)

#ifndef ROOT_ktPy8Event
#define ROOT_ktPy8Event

#include "TObject.h"
#include "TBuffer.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "ktParticle.h"
#include "ktPID.h"
#include "TString.h"
#include "ktGrid.h"

class ktPy8Event : public TObject
{

 private:

  Double_t mSqrts;
  Int_t evNum;
  
  //Int_t nCharged;
  //Int_t nNeutral;
  //Int_t nAll;
  
  Double_t mptHat;
  Double_t mphiHat;
  Double_t mthetaHat;
  Double_t mWeight;

  Double_t pycellPt;
  Double_t pycellEta;
  Double_t pycellPhi;
  
  TObjArray *neutralList;
  TObjArray *chargedList;
  TObjArray *allList;
  
  
 public:
  
  ktPy8Event();
  ktPy8Event(Int_t myEvNum);
  virtual ~ktPy8Event();

  Double_t GetPtHat() {return mptHat;}
  Double_t GetPyCellEta() {return pycellEta;}
  Double_t GetPyCellPhi() {return pycellPhi;}
  Double_t GetPyCellPt() {return pycellPt;}
  Double_t GetSqrtS() {return mSqrts;}
  Double_t GetWeight() {return mWeight;}
  Double_t GetPhiHat() {return mphiHat;}
  Double_t GetThetaHat() {return mthetaHat;}
  Int_t GetEventNumber() {return evNum;}
  Int_t GetNCharged() {return chargedList->GetEntries();}
  Int_t GetNNeutral() {return neutralList->GetEntries();}
  Int_t GetNAll() {return (allList->GetEntries());}
  Int_t GetN() {return (GetNCharged()+GetNNeutral()+allList->GetEntries());}

  void SetPtHat(Double_t m_ptHat) {mptHat=m_ptHat;}
  void SetPhiHat(Double_t m_mphiHat) {mphiHat=m_mphiHat;}
  void SetThetaHat(Double_t m_mthetaHat) {mthetaHat=m_mthetaHat;}
  // only highest pt/et pycell jet ....
  void SetPyCellEta(Double_t m_pycellEta) {pycellEta=m_pycellEta;}
  void SetPyCellPhi(Double_t m_pycellPhi) {pycellPhi=m_pycellPhi;}
  void SetPyCellPt(Double_t m_pycellPt) {pycellPt=m_pycellPt;}
  void SetSqrtS(Double_t m_mSqrts) {mSqrts=m_mSqrts;}
  void SetWeight(Double_t m_mWeight) {mWeight=m_mWeight;}
  void SetEventNumber(Int_t m_evNum) {evNum=m_evNum;}

  void Add(ktParticle *v,TString mType);

  //void Fill(ktGrid* myGrid);
  void Fill(ktGrid* myGrid,TString mSel);

  void EventInfo();

  ktParticle* GetNeutral(Int_t n) {return (ktParticle*) neutralList->At(n);}
  ktParticle* GetCharged(Int_t n) {return (ktParticle*) chargedList->At(n);}
  ktParticle* GetAll(Int_t n) {return (ktParticle*) allList->At(n);}

  TObjArray* GetAll() {return allList;}
  TObjArray* GetCharged() {return allList;}
  TObjArray* GetNeutral() {return allList;}

  ClassDef(ktPy8Event,1)
};

#endif
