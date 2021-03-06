

#ifndef ROOT_ktHiFastJet
#define ROOT_ktHiFastJet

#include "TObject.h"
#include "TBuffer.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "ktPythia8.h"
#include "ktMuEvent.h"
#include "ktJetQuench.h"

class ktHiFastJet : public TObject 
{

  private:
  
  Double_t mEtaMax,mPhiMax;
  Double_t mEtaAna;
  Double_t mPtAna;
  Double_t mPtCut;
  
  Double_t median_pt_per_area;
  Double_t median_pt_per_area_charged;
  Double_t median_pt_per_area_neutral;

  Double_t ghost_etamax;
  Double_t ghost_area;
  Int_t active_area_repeats;

  Bool_t info, printjet;
  Bool_t spaceCharge,smear;

  TObjArray *FFArray;
  //TObjArray *Constituents; // add if needed to store the constituents from FastJet 
  TObjArray *PhiBkgHistos;

  Bool_t fillFF;
  Bool_t bkgCalculated;
  Bool_t histFFonly;

  Double_t RcFF;
  Int_t nFF;

  //TH1D *hptBkgFJ;

  void GetJetFF(ktMuFastJet *mJ,Bool_t setInJet,Bool_t mSave);
  Bool_t IsTriggerJet(ktMuFastJet *mJ,ktMuEvent *ev); // check only the first HT trigger info
  Double_t GetdPhi(Double_t mphi,Double_t vphi);

  public:

  ktHiFastJet();
  ktHiFastJet(Int_t nXiBins,Double_t xiMax);
  virtual ~ktHiFastJet();

  void AddParticle(Double_t px, Double_t py, Double_t pz, Double_t E); // include pt,cut
  void AddParticle(Double_t px, Double_t py, Double_t pz, Double_t E, Bool_t charged);
  void AddParticle(Double_t px, Double_t py, Double_t pz, Double_t E, Bool_t charged,Int_t sign);
  void AddParticle(Double_t px, Double_t py, Double_t pz, Double_t E, Bool_t charged, ktPID* mPid);
  void AddParticle(Double_t px, Double_t py, Double_t pz, Double_t E, ktPID* mPid);
  void Init(); // Dummy !
  void Clear();

  void RunFastJetSub(Double_t Rparam=1.0,TString mAlgo="kt",ktMuEvent *ev=0); 
  void RunSISCone(Double_t Rc=1.0,ktMuEvent *ev=0);
  //void DoXiBkg(TString mAlgo="kt",Int_t nJet=2, ktMuEvent *ev=0);
  void FillFF(TString mAlgo="kt", ktMuEvent *ev=0);//, Int_t nJet=2);
  void DoXiBkg(TString mAlgo="kt",ktMuEvent *ev=0);
  void DoPhiBkg(TString mAlgo="kt",ktMuEvent *ev=0);
  void DoBkg(Double_t Rparam=0.4,TString mAlgo="kt");

  void AddPythia8Event(ktPythia8 *myp,TString mSel);
  void AddQuenchedJet(ktJetQuench *mQ=0);
  void AddUnQuenchedJet(ktJetQuench *mQ=0); 

  void ClearFFArray();

  Bool_t GetHighestJetForQuenching(ktJetQuench *mQ=0,Double_t Rparam=0.7,TString mAlgo="kt");
  Int_t GetNPhiBkgHistos() {return PhiBkgHistos->GetEntries();}
  TH2D* GetPhiBkgHisto(Int_t mN) {return (TH2D*) PhiBkgHistos->At(mN);}

  void PrintInfo();
  void PrintAreaDefinition();
  void PrintFFArray();

  void SetActiveAreaRepeats(Int_t m_active_area_repeats) {active_area_repeats=m_active_area_repeats;}
  void SetPtCut(Double_t m_mPtCut) {mPtCut=m_mPtCut;}
  void SetFiducial(Double_t mEta,Double_t mPhi);
  void SetEtaAnaFiducial(Double_t mEta) {mEtaAna=mEta;}
  void SetPtAna(Double_t mPt) {mPtAna=mPt;}
  void SetInfo(Bool_t mInfo) {info=mInfo;}
  void SetPrintJet(Bool_t mPrintjet) {printjet=mPrintjet;}
  void SetAreaDefinition(Double_t m_ghost_etamax, Int_t m_active_area_repeats, Double_t m_ghost_area);

  void SetFillFFArray(Bool_t m_fillFF) {fillFF=m_fillFF;}
  void SetRcFF(Double_t m_RcFF) {RcFF=m_RcFF;}
  void SetNJetsForFF(Int_t m_nFF) {nFF=m_nFF;}
  void SetSpaceCharge(Bool_t mSpace) {spaceCharge=mSpace;}
  void SetTrackSmear(Bool_t mSmear) {smear=mSmear;}
  void SetFFHistOnly(Bool_t mhistFFonly) {histFFonly=mhistFFonly;}

  Double_t GetPtCut() {return mPtCut;}
  Double_t GetMedianPtPerArea() {return median_pt_per_area;}
  Double_t GetMedianPtPerAreaCharged() {return median_pt_per_area_charged;}
  Double_t GetMedianPtPerAreaNeutral() {return median_pt_per_area_neutral;}
  Int_t GetNCharged() {return FFArray->GetEntriesFast();}
  Bool_t GetFillFFArray() {return fillFF;}
  Double_t GetRcFF() {return RcFF;}
  Int_t GetNJetsForFF() {return nFF;}
  TObjArray* GetFFArray(){return FFArray;}

  ClassDef(ktHiFastJet,3)
};

#endif
