// first test (Joern Putschke)

#ifndef ROOT_ktJetQuench
#define ROOT_ktJetQuench

#include "TObject.h"
#include "TBuffer.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "ktParticle.h"

class ktJetQuench : public TObject 
{

  private:

  TF1 *fpt;
  //TF1 *feta, *fphi, *fpt;

  TObjArray *hJetParticles;
  TObjArray *hJetParticlesNoQuench;

  Double_t jPt,jPhi,jEta;
  Double_t S0;
  Double_t C;
  Double_t T;
  Double_t R;
  Double_t minQpt,maxQpt;
  Bool_t verbose;
  Double_t dE;
  Double_t ptQuench;
  Double_t ptNoQuench;
  Int_t maxIter;
  Bool_t doQA;
  Int_t nEvents;
  Bool_t init;
  Bool_t hScale;
  Double_t minDi,maxDi,assDi;
  Int_t nTrig;
  Int_t nTrigNo;
  Bool_t fracQuench;
  TString mQMethod;
  Double_t PE;

  void ReDistribute();
  void FillQA();
  void DoDiCorr(TObjArray *mA,TH2D *mH,Bool_t quenched);

  public:

  TH1D *hXi;
  TH1D *hXiNoQuench;
  TH1D *hPt, *hPtQ, *hPtNo;
  TH1D *hJptQ, *hJptRe;
  TH2D *hdEtadPhi, *hdEtadPhiNo;
  TH1D *hdRQ;
  TH1D *hdR;
  TH2D *hdRptQ;
  TH2D *hdRpt;
  
  ktJetQuench();

  void DoQuenching();
  void DoDiHadron(Double_t minTrig=3.0, Double_t maxTrig=4.0, Double_t ptAss=2.0);
  void AddParticle(Double_t px, Double_t py, Double_t pz, Double_t E);
  void AddParticle(Double_t px, Double_t py, Double_t pz, Double_t E, Int_t mPid);
  void PrintInfo();
  void DrawQAPlots(TString mT="FastJet");
  void Clear();
  void InitQuenching(TString qMethod="fractional");
  void SaveQA(TString fName="test.root");

  Int_t GetNJetParticles() {return hJetParticles->GetEntries();}
  Int_t GetNJetParticlesNoQuench() {return hJetParticlesNoQuench->GetEntries();}

  ktParticle* GetQuenchedParticle(Int_t n) {return (ktParticle*) hJetParticles->At(n);}
  ktParticle* GetUnQuenchedParticle(Int_t n) {return (ktParticle*) hJetParticlesNoQuench->At(n);}

  Double_t GetJetPt() {return jPt;}
  Double_t GetJetPhi() {return jPhi;}
  Double_t GetJetEta() {return jEta;}
  Double_t GetSumPt();
  Double_t GetSumPtNoQuench();
  Double_t GetQuenchingFraction() {return S0;}
  Double_t GetdE() {return dE;}
  Double_t GetPtQuench() {return ptQuench;}
  Double_t GetPtNoQuench() {return ptNoQuench;}
  Double_t GetQuenchedSlope() {return T;}
  Double_t GetQuenchedRadius() {return R;}
  Double_t GetQuenchConstant() {return C;}
  Double_t GetPartonicEnergyLoss() {return PE;}

  void SetJetPt(Double_t m_jPt) {jPt=m_jPt;}
  void SetJetPhi(Double_t m_jPhi) {jPhi=m_jPhi;}
  void SetJetEta(Double_t m_jEta) {jEta=m_jEta;}
  void SetMaxIterations(Int_t m_maxIter) {maxIter=m_maxIter;}
  void SetQuenchingFraction(Double_t m_S0) {S0=m_S0;}
  void SetVerbose(Bool_t m_verbose) {verbose=m_verbose;}
  void SetDoQA(Bool_t m_doQA) {doQA=m_doQA;}
  void SetQuenchedSlope(Double_t m_T) {T=m_T;}
  void SetQuenchedRadius(Double_t m_R) {R=m_R;}
  void SetQuenchPtRange(Double_t m_min,Double_t m_max);
  void SetQuenchConstant(Double_t m_C) {C=m_C;}
  void SetPartonicEnergyLoss(Double_t m_PE) {PE=m_PE;}

  virtual ~ktJetQuench();
 
  ClassDef(ktJetQuench,1)
};

#endif
