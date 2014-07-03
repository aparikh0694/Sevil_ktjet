#ifndef ROOT_ktAna
#define ROOT_ktAna

#include "TBuffer.h"
#include "TObject.h"
#include "TH1.h"
#include "TH2.h"
#include "ktMuEvent.h"
#include "ktMuJet.h"
#include "ktMuFastJet.h"

class ktAna : public TObject
{

 protected:

  Int_t nEvents;
  Double_t minPtCut;
  Double_t minFakePtCut;
  Double_t fakeMindPhi,fakeMaxdPhi;
  Double_t etaMin,etaMax;
  Double_t phiMin,phiMax;
  Double_t diPhiMin,diPhiMax;
  Double_t histPtMax;
  Double_t ffPtMin,ffPtMax;
  Double_t diJetPtMin,diJetPtMax;

  Double_t evPt,evPtDi;
  Double_t tevPt,tevPtDi;
  Double_t evArea,evAreaDi;
  Double_t tevArea,tevAreaDi;

  Int_t nJets;
  Int_t refMultCut;

  Bool_t InitHist;
  Bool_t verbose;
  Bool_t FJAna;
  Bool_t BkgCorr;

  TString JetAlgo;

  void ScaleHist(TH1D *mH);
  void ScaleHist(TH1D *mH, Double_t mS);
  void LabelHist(TH1D *mH,TString xA,TString yA);
  void LabelHist(TH2D *mH,TString xA,TString yA);

  Bool_t IsDiJet(Double_t mPhi);
  Bool_t IsFakeJet(Double_t mPhi);

  Int_t NAlgoJets(ktMuEvent *ev);
  ktMuFastJet* GetFJjet(ktMuEvent *ev,Int_t n=0);
  ktMuJet* GetJet(ktMuEvent *ev,Int_t n=0);

  TH1D* GetXiBkg(ktMuEvent *ev);
  TH1D* GetZBkg(ktMuEvent *ev);
  TH1D* GetptBkg(ktMuEvent *ev);

  void FillFFHists(TH1D *hxi,TH1D *hz,ktMuFastJet *mJ,Double_t jetPt);

  Bool_t CheckAcc(ktMuFastJet *j);
  Bool_t CheckAcc(ktMuJet *j);

  TObjArray *hArray;

 public:

  TH1D *hE,*hEDi,*hPhiDi,*hXi,*hXiBkg,*hXiDi,*hXiBkgDi;
  TH1D *hZ,*hZBkg,*hZDi,*hZBkgDi;

  TH1D *thE,*thEDi,*thPhiDi,*thXi,*thXiBkg,*thXiDi,*thXiBkgDi;
  TH1D *thZ,*thZBkg,*thZDi,*thZBkgDi;
  TH1D *thSpec, *thSpecDi;
  TH2D *thAreaPt,*thAreaPtDi;
  TH1D *thPtDiff,*thPtDiffDi;

  TH1D *hSpec, *hSpecDi;
  TH1D *hPtDiff,*hPtDiffDi;

  TH2D *hEtaPhi,*hEtaPhiDi;
  TH2D *hAreaPt,*hAreaPtDi;

  TH2D *hdPhiPt;

  TH1D *hdEdi;

  TH1D *hXiFake,*hXiBkgFake;
  TH1D *hEFake;

  TH1D *hEPtNoPt;

  Bool_t TriggerPtToRecoil;

  //TH1D *hFakeJetRatio;

  ktAna();
  ktAna(Double_t eMin, Double_t eMax, Double_t pMin, Double_t pMax, Double_t ptMin);
  virtual ~ktAna();
 
  void SetCorrectionFile(TString fName); // dummy
  void DoSpectraCorrections(); //dummy

  void AnalyzeEvent(ktMuEvent *ev);
  void AnalyzeEventPtAndNoPt(ktMuEvent *ev,ktMuEvent *ev2);
  void AnalyzeEventTrig(ktMuEvent *ev); //FJ only ! (dummy so far ...)
  void AnalyzeEvent(ktMuEvent *ev,TString m_JetAlgo);
  void InitHistograms();
  void InitHistograms(TString m_JetAlgo);
  void InitHistograms(TString m_JetAlgo,Double_t m_histPtMax);
  void PrintAnalysisCuts();
  void Finish();
  void DrawPlots();
  void DrawFFPlots();
  void DrawFFPlotsFake();
  void DrawFFPlotsTrig();
  void DrawPlotsPtNoPt();

  void CompareEvent(ktAna *e);
  void CompareEventTrig(ktAna *e);
  void DoFakeJetEstimate();

  Double_t GetdPhi(Double_t mphi,Double_t vphi); 

  TH1D* GetHist(TString hN="hE",Int_t col=1, Int_t lStyle=1, Int_t lWidth=1);
  TH1D* GetMarkerHist(TString hN="hE",Int_t mStyle=22, Int_t mCol=1, Int_t mSize=1);

  // leading and second leading 
  // xi dist.
  TH1D* GetXi(Int_t mStyle=22, Int_t mCol=1, Int_t mSize=1) {return GetMarkerHist("hXi",mStyle,mCol,mSize);}
  TH1D* GetXiBkg(Int_t col=1, Int_t lStyle=1, Int_t lWidth=1) {return GetHist("hXiBkg",col,lStyle,lWidth);}
  TH1D* GetXiSub(Int_t mStyle=22, Int_t mCol=1, Int_t mSize=1);
  TH1D* GetDiJetXi(Int_t mStyle=22, Int_t mCol=1, Int_t mSize=1) {return GetMarkerHist("hXiDi",mStyle,mCol,mSize);}
  TH1D* GetDiJetXiBkg(Int_t col=1, Int_t lStyle=1, Int_t lWidth=1) {return GetHist("hXiBkgDi",col,lStyle,lWidth);}
  TH1D* GetDiJetXiSub(Int_t mStyle=22, Int_t mCol=1, Int_t mSize=1);
  
  // fake 
  TH1D* GetXiFake(Int_t mStyle=22, Int_t mCol=1, Int_t mSize=1) {return GetMarkerHist("hXiFake",mStyle,mCol,mSize);}
  TH1D* GetXiBkgFake(Int_t col=1, Int_t lStyle=1, Int_t lWidth=1) {return GetHist("hXiBkgFake",col,lStyle,lWidth);}
  TH1D* GetXiSubFake(Int_t mStyle=22, Int_t mCol=1, Int_t mSize=1); 

  // z dist.
  TH1D* GetZ(Int_t mStyle=22, Int_t mCol=1, Int_t mSize=1) {return GetMarkerHist("hZ",mStyle,mCol,mSize);}
  TH1D* GetZBkg(Int_t col=1, Int_t lStyle=1, Int_t lWidth=1) {return GetHist("hZBkg",col,lStyle,lWidth);}
  TH1D* GetZSub(Int_t mStyle=22, Int_t mCol=1, Int_t mSize=1); // do proper z bkg. shift
  TH1D* GetDiJetZ(Int_t mStyle=22, Int_t mCol=1, Int_t mSize=1) {return GetMarkerHist("hZDi",mStyle,mCol,mSize);}
  TH1D* GetDiJetZBkg(Int_t col=1, Int_t lStyle=1, Int_t lWidth=1) {return GetHist("hZBkgDi",col,lStyle,lWidth);}
  TH1D* GetDiJetZSub(Int_t mStyle=22, Int_t mCol=1, Int_t mSize=1);
  
  // spectrum
  TH1D* GetSpectrum(Int_t mStyle=22, Int_t mCol=1, Int_t mSize=1) {return GetMarkerHist("hSpec",mStyle,mCol,mSize);}
  TH1D* GetDiJetSpectrum(Int_t mStyle=22, Int_t mCol=1, Int_t mSize=1) {return GetMarkerHist("hSpecDi",mStyle,mCol,mSize);}
  TH1D* GetSpectrumCorr(Int_t mStyle=22, Int_t mCol=1, Int_t mSize=1);
  
  // N jets 
  Int_t GetNjets(Double_t ptMin, Double_t ptMax);
  Int_t GetNDiJets(Double_t ptMin, Double_t ptMax);
  Int_t GetNFakeJets(Double_t ptMin, Double_t ptMax);
  
  // trigger and recoil
  // xi dist.
  TH1D* GetXiTrig(Int_t mStyle=22, Int_t mCol=1, Int_t mSize=1) {return GetMarkerHist("thXi",mStyle,mCol,mSize);}
  TH1D* GetXiBkgTrig(Int_t col=1, Int_t lStyle=1, Int_t lWidth=1) {return GetHist("thXiBkg",col,lStyle,lWidth);}
  TH1D* GetXiSubTrig(Int_t mStyle=22, Int_t mCol=1, Int_t mSize=1);
  TH1D* GetDiJetXiTrig(Int_t mStyle=22, Int_t mCol=1, Int_t mSize=1) {return GetMarkerHist("thXiDi",mStyle,mCol,mSize);}
  TH1D* GetDiJetXiBkgTrig(Int_t col=1, Int_t lStyle=1, Int_t lWidth=1) {return GetHist("thXiBkgDi",col,lStyle,lWidth);}
  TH1D* GetDiJetXiSubTrig(Int_t mStyle=22, Int_t mCol=1, Int_t mSize=1);
  
  // z dist.
  TH1D* GetZTrig(Int_t mStyle=22, Int_t mCol=1, Int_t mSize=1) {return GetMarkerHist("thZ",mStyle,mCol,mSize);}
  TH1D* GetZBkgTrig(Int_t col=1, Int_t lStyle=1, Int_t lWidth=1) {return GetHist("thZBkg",col,lStyle,lWidth);}
  TH1D* GetZSubTrig(Int_t mStyle=22, Int_t mCol=1, Int_t mSize=1); // do proper z bkg. shift
  TH1D* GetDiJetZTrig(Int_t mStyle=22, Int_t mCol=1, Int_t mSize=1) {return GetMarkerHist("thZDi",mStyle,mCol,mSize);}
  TH1D* GetDiJetZBkgTrig(Int_t col=1, Int_t lStyle=1, Int_t lWidth=1) {return GetHist("thZBkgDi",col,lStyle,lWidth);}
  TH1D* GetDiJetZSubTrig(Int_t mStyle=22, Int_t mCol=1, Int_t mSize=1);
  
  // spectrum
  TH1D* GetSpectrumTrig(Int_t mStyle=22, Int_t mCol=1, Int_t mSize=1) {return GetMarkerHist("thSpec",mStyle,mCol,mSize);}
  TH1D* GetDiJetSpectrumTrig(Int_t mStyle=22, Int_t mCol=1, Int_t mSize=1) {return GetMarkerHist("thSpecDi",mStyle,mCol,mSize);}
  TH1D* GetSpectrumCorrTrig(Int_t mStyle=22, Int_t mCol=1, Int_t mSize=1);
 
  TH1D* GetdEdi(Int_t col=1, Int_t lStyle=1, Int_t lWidth=1) {return GetHist("hdEdi",col,lStyle,lWidth);}
  
  //TH1D* GetFakeJetRatio(Int_t col=1, Int_t lStyle=1, Int_t lWidth=1) {return GetHist("hFakeJetRatio",col,lStyle,lWidth);}

  // N jets 
  Int_t GetNjetsTrig(Double_t ptMin, Double_t ptMax);
  Int_t GetNDiJetsTrig(Double_t ptMin, Double_t ptMax);
  
  void XiShift(TH1D* m_shift,Float_t a);
  void ZShift(TH1D* m_shift,Float_t a); //dummy !!!

  void SetMinPtCut(Double_t m_minPtCut) {minPtCut=m_minPtCut;}
  void SetMinFakePtCut(Double_t m_minFakePtCut) {minFakePtCut=m_minFakePtCut;}
  void SetFiducial(Double_t eMin, Double_t eMax, Double_t pMin, Double_t pMax);
  void SetEtaFiducial(Double_t eMin,Double_t eMax) {etaMax=eMax;etaMin=eMin;}
  void SetDijetPhiCuts(Double_t pMin, Double_t pMax);
  void SetJetAlgo(TString m_JetAlgo) {JetAlgo=m_JetAlgo;}
  void SetVerbose(Bool_t m_verbose) {verbose=m_verbose;}
  void SetHistPtMax(Double_t m_histPtMax) {histPtMax= m_histPtMax;}
  void SetDoBkgCorrections(Bool_t m_BkgCorr) {BkgCorr=m_BkgCorr;}
  void SetFFPtRange(Double_t ptMin, Double_t ptMax);
  void SetDiJetPtRange(Double_t ptMin, Double_t ptMax);
  void SetRefMultCut(Int_t m_refMultCut) {refMultCut= m_refMultCut;}

  void DoXiBkgFromPt(TH1D *hxi,TH1D *hpt,Double_t jpt);
  void DoZBkgFromPt(TH1D *hz,TH1D *hpt,Double_t jpt);

  Int_t GetNEvents() {return nEvents;}
  Double_t GetMinPtCut() {return minPtCut;}
  Double_t GetEtaMax() {return etaMax;}
  Double_t GetPhiMax() {return phiMax;}
  Double_t GetEtaMin() {return etaMin;}
  Double_t GetPhiMin() {return phiMin;}
  TString GetJetAlgo() {return JetAlgo;}
  Int_t GetRefMultCut() {return refMultCut;}

  Double_t GetEventPt() {return evPt;}
  Double_t GetEventDiJetPt() {return evPtDi;}
  Double_t GetEventArea() {return evArea;}
  Double_t GetEvenDiJetArea() {return evAreaDi;}

  Double_t GetEventPtTrig() {return tevPt;}
  Double_t GetEventDiJetPtTrig() {return tevPtDi;}
  Double_t GetEventAreaTrig() {return tevArea;}
  Double_t GetEvenDiJetAreaTrig() {return tevAreaDi;}
   
  void SetTriggerPtToRecoil(Bool_t m_TriggerPtToRecoil) {TriggerPtToRecoil=m_TriggerPtToRecoil;}

  ClassDef(ktAna,2)
};

#endif
