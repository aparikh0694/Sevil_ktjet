// first test (Joern Putschke)

#ifndef ROOT_ktStarPico
#define ROOT_ktStarPico

#include "TObject.h"
#include "TBuffer.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TNtuple.h"
#include "TString.h"
#include "ktGrid.h"
#include "ktPID.h"
#include "ktFastJet.h"
#include "ktTriggerInfo.h"
#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetVectorContainer.h"
#include "TStarJetVector.h"
#include "TStarJetPicoTriggerInfo.h"

class ktStarPico : public TObject
{

  private:

  Int_t NPrim;
  Int_t NTowers;
  Int_t RefMult;
  Int_t nEvents;

  Double_t meanZdc;
  
  Bool_t QAOutput;
  // Bool_t ElectronCheck;
  Bool_t info;
  Bool_t SC;

  TString SCName;
  TString QAName;

  // QA histograms
  TH1D *hFitPoints, *hFitOverMax, *hDca, *hEta, *hPhi, *hPt;
  TH1D *hMatchedFitPoints, *hMatchedFitOverMax, *hMatchedDca, *hMatchedEta, *hMatchedPhi;
  TH1D *hMatchedPoverE, *hMatchedPt;

  TH1D *hTowADC, *hTowEta, *hTowPhi, *hTowE,*hTowEt;
  TH1D *hMatchedTowADC, *hMatchedTowEta, *hMatchedTowPhi, *hMatchedTowE, *hMatchedTowEt;

  TH1D *hTowIDEweighted;
  TH1D *hTowIDEweightedEt;

  TH1D *hVertexZ;
  TH1D *hRefMult;
  TH1D *hTrigId;

  TH1D *hV0MassLa;
  TH1D *hV0MassALa;
  TH1D *hV0MassK0;

  TH1D *hEplusPt;
  TH1D *hEminusPt;
  TH1D *hPlusPt;
  TH1D *hMinusPt;

  TH1D *hEplusPtW;
  TH1D *hEminusPtW;
  TH1D *hPlusPtW;
  TH1D *hMinusPtW;

  TH1D *ehEplusPtW;
  TH1D *ehEminusPtW;
  TH1D *ehPlusPtW;
  TH1D *ehMinusPtW;

  TH1D *whEplusPtW;
  TH1D *whEminusPtW;
  TH1D *whPlusPtW;
  TH1D *whMinusPtW;

  TH2D *zhEplusPtW;
  TH2D *zhEminusPtW;
  TH2D *zhPlusPtW;
  TH2D *zhMinusPtW;

  TH2D *zehEplusPtW;
  TH2D *zehEminusPtW;
  TH2D *zehPlusPtW;
  TH2D *zehMinusPtW;

  TH2D *zwhEplusPtW;
  TH2D *zwhEminusPtW;
  TH2D *zwhPlusPtW;
  TH2D *zwhMinusPtW;

  TH2D *zdchEplusPtW;
  TH2D *zdchEminusPtW;
  TH2D *zdchPlusPtW;
  TH2D *zdchMinusPtW;

  TH2D *zdcehEplusPtW;
  TH2D *zdcehEminusPtW;
  TH2D *zdcehPlusPtW;
  TH2D *zdcehMinusPtW;

  TH2D *zdcwhEplusPtW;
  TH2D *zdcwhEminusPtW;
  TH2D *zdcwhPlusPtW;
  TH2D *zdcwhMinusPtW;

  TH1D *hNevents;
  TH1D *hZDCcoin;

  //TNtuple *mNtuple;

  TClonesArray* TriggerInfoArray;

  void WriteHistograms();

  public:

  ktStarPico(Bool_t mQAOut=false,Bool_t mSC=false);
  virtual ~ktStarPico();

  Bool_t Fill(ktGrid *mGrid,ktFastJet *mFast,TStarJetPicoEvent *mEv,TStarJetVectorContainer<TStarJetVector>* container);

  void PrintQACuts();
  void PrintCuts();
  void WriteQAOutputFile(TString fName);
  void PrintEventInfo();
  void PrintTriggerInfo();

  void DoSpaceCharge(TLorentzVector *p,ktPID *pid,Float_t ZdcCoin=-1);

  void Clear();

  // Setter
  //void SetCheckForElectron(Bool_t mElectronCheck) {ElectronCheck=mElectronCheck;}
  void SetTriggerInfo(TClonesArray* trinfo);
  void SetInfo(Bool_t mInfo) {info=mInfo;}
  //void SetDoSpaceCharge(Bool_t mSC) {SC=mSC;} // is done in constructor 
  //void SetDoQA(Bool_t mQAOutput) {QAOutput=mQAOutput;} // is done in constructor 
  void SetSpaceChargeFile(TString name) {SCName=name;}
  void SetQAFile(TString name) {QAName=name;}
  void SetMeanZDC(Double_t mZdc) {meanZdc=mZdc;}

  // Getter
  Bool_t GetQAOutput() {return QAOutput;}
  Int_t GetNTowers() {return NTowers;}
  Int_t GetNPrimaries() {return NPrim;} 
  Int_t GetRefMult() {return RefMult;}
  TClonesArray* GetTriggerInfo() {return TriggerInfoArray;}

  ClassDef(ktStarPico,3)
};

#endif
