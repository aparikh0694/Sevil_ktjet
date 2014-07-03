// first test (Joern Putschke)

#ifndef ROOT_ktMCBkg
#define ROOT_ktMCBkg

#include "TObject.h"
#include "TBuffer.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "ktGrid.h"
#include "ktFastJet.h"

class ktMCBkg : public TObject
{

 private:
  
 TF1 *feta, *fphi, *fpt, *fRP;
 TF1 *fv2pt;

 Double_t v2,T,dNdEta;
 Double_t maxEta;
 TString mSel,mSys;
 Double_t phiRP;
 Bool_t mInitOK;
 Bool_t output;

 public:

 TH1D *heta, *hphi, *hpt, *hRP;

 ktMCBkg();
 ktMCBkg(TString m_mSel,TString m_mSys, Double_t m_v2);
 ktMCBkg(TString m_mSel,TString m_mSys, Double_t m_T, Double_t m_dNdEta, Double_t m_v2);
 virtual ~ktMCBkg();

 // function hardcoded maybe include setter for
 // all sort of functions (TF1) !!!

 Bool_t Init();

 void SetRandomSeed();
 void Fill(ktGrid* myGrid);
 void Fill(ktFastJet *myFastJet,Double_t m_RP);
 void Fill(ktGrid* myGrid, Double_t m_RP);
 void Fill(ktGrid* myGrid,ktFastJet *myFastJet,Double_t m_RP);
 void Fill(ktGrid* myGrid,TString mSel);
 void FillNeutral(ktGrid* myGrid);
 void HistNorm(Int_t nEvents);
 
 Double_t GetPhiRP() {return phiRP;};
 Bool_t GetInitOK() {return mInitOK;};

 void Setv2(Double_t mv2) {v2=mv2;};
 void SetT(Double_t mT) {T=mT;};
 void SetdNdEta(Double_t mdNdEta) {dNdEta=mdNdEta;};
 void SetCollider(TString m_mSel) {mSel=m_mSel;};
 void SetSystem(TString m_mSys) {mSys=m_mSys;};
 void SetOutput(Bool_t m_output) {output=m_output;};

 Double_t Getv2() {return v2;};
 Double_t GetT() {return T;};
 Double_t GetdNdEta() {return dNdEta;};
 TString GetCollider() {return mSel;};
 TString GetSystem() {return mSys;};
 Bool_t GetOutput() {return output;};

 void PrintInfo();

 ClassDef(ktMCBkg,0)
};

#endif
