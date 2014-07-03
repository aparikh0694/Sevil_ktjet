// first test (Joern Putschke)

#ifndef ROOT_ktMuEvent
#define ROOT_ktMuEvent

#include "TObject.h"
#include "TBuffer.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "ktMuJet.h"
#include "ktMuFastJet.h"
#include "ktJet.h"
#include "ktParton.h"
#include "ktGrid.h"
#include "TString.h"
#include "ktUEEvent.h"
#include "ktTriggerInfo.h"

class ktFastJet;
 
class ktMuEvent : public TObject
{

  private:
  
  Int_t evNum;
  //Int_t trigger;
  Int_t refMult;
  Int_t TrigId;
  Double_t PsiRP;
  Double_t vertexZ;

  Double_t BkgPtPerCell;//!
  Double_t BkgPtPerCellAll;//!
  Double_t BkgPtPerCellIn;//!
  Double_t BkgPtPerCellAllIn;//!
  Double_t BkgPtPerCellOut;//!
  Double_t BkgPtPerCellAllOut;//!
  Double_t BkgPtCone;//!
  Double_t BkgPtCluster;//!
  Double_t BkgNCellCone;//!
  Double_t BkgNCellCluster;//!
  // maybe not per event, only once ... (check also if anything has to be double ...) !!!
  Double_t NCellsInJetArea;//!

  Double_t median_pt_per_area;

  // Add more if necessary ... ;-)

  TObjArray *ktList; //!
  TObjArray *ktBFList; //!
  TObjArray *coneList; //!
  TObjArray *allList; //!
  TObjArray *PartonList;//!

  TObjArray *fjKt;
  TObjArray *fjAntiKt;
  TObjArray *fjSISCone;

  TObjArray *ueEvents;
  //TClonesArray* TriggerInfoArray;
  //static TClonesArray *fgTrigObjs;//!
  TObjArray *TriggerInfoArray;

  Double_t BkgNCellIn;//!
  Double_t BkgNCellOut;//!
  Double_t BkgNCellInAll;//!
  Double_t BkgNCellOutAll;//!
  Double_t BkgNJetCellIn;//!
  Double_t BkgNJetCellOut;//!

  Double_t jfPt;
  Double_t jfRc;
  Double_t jfRcBkg;
  Bool_t jfEmcalSet;
  TString jfPid;

  Int_t jfNSubJets;//!
 
  TH1D *hXiBkg;//!
  TH1D *hXiBkgKt;//!
  TH1D *hXiBkgAntiKt;//!
  TH1D *hXiBkgSISCone;//!

  TH1D *hZBkg;//!
  TH1D *hZBkgKt;//!
  TH1D *hZBkgAntiKt;//!
  TH1D *hZBkgSISCone;//!

  TH1D *hptBkg;
  TH1D *hptBkgKt;
  TH1D *hptBkgAntiKt;
  TH1D *hptBkgSISCone;

  public:

  ktMuEvent();
  ktMuEvent(Int_t mEvNum);
  virtual ~ktMuEvent();
  //void Delete();

  Int_t GetEventNumber() {return evNum;}
  //Int_t GetTrigger() {return trigger;}
  Int_t GetNKtJets() {return ktList->GetEntriesFast();}
  Int_t GetNKtBFJets() {return ktBFList->GetEntriesFast();}
  Int_t GetNConeJets() {return coneList->GetEntriesFast();}
  Int_t GetNAllJets() {return allList->GetEntriesFast();}
  Int_t GetNPartons() {return PartonList->GetEntriesFast();}

  Int_t GetNFJkt() {return fjKt->GetEntriesFast();}
  Int_t GetNFJAntiKt() {return fjAntiKt->GetEntriesFast();}
  Int_t GetNFJSISCone() {return fjSISCone->GetEntriesFast();}

  Double_t GetBkgPtPerCell() {return BkgPtPerCell;}
  Double_t GetBkgPtPerCellAll() {return BkgPtPerCellAll;}
  Double_t GetBkgPtPerCellIn() {return BkgPtPerCellIn;}
  Double_t GetBkgPtPerCellAllIn() {return BkgPtPerCellAllIn;} 
  Double_t GetBkgPtCone() {return BkgPtCone;}
  Double_t GetBkgPtCluster() {return BkgPtCluster;}
  Double_t GetBkgNCellCone() {return BkgNCellCone;}
  Double_t GetBkgNCellCluster() {return BkgNCellCluster;}
  Double_t GetBkgNCellIn() {return BkgNCellIn;}
  Double_t GetBkgNCellOut() {return BkgNCellOut;}
  Double_t GetBkgNJetCellIn() {return BkgNJetCellIn;}
  Double_t GetBkgNJetCellOut() {return BkgNJetCellOut;}
  Double_t GetBkgNCellInAll() {return BkgNCellInAll;}
  Double_t GetBkgNCellOutAll() {return BkgNCellOutAll;}
  Double_t GetNCellsInJetArea() {return  NCellsInJetArea;}

  Double_t GetMedianPtPerArea() {return median_pt_per_area;}

  Int_t GetRefMult(){return refMult;}
  Int_t GetTrigId() {return TrigId;}
  Double_t GetRPAngle() {return PsiRP;}
  Double_t GetVertexZ() {return vertexZ;}
  //TClonesArray* GetTriggerInfo() {return TriggerInfoArray;}
  TObjArray* GetTriggerInfo() {return TriggerInfoArray;}

  TH1D* GetXiBkg() {return hXiBkg;}
  TH1D* GetXiBkgKt() {return hXiBkgKt;}
  TH1D* GetXiBkgAntiKt() {return hXiBkgAntiKt;}
  TH1D* GetXiBkgSISCone() {return hXiBkgSISCone;}

  TH1D* GetZBkg() {return hZBkg;}
  TH1D* GetZBkgKt() {return hZBkgKt;}
  TH1D* GetZBkgAntiKt() {return hZBkgAntiKt;}
  TH1D* GetZBkgSISCone() {return hZBkgSISCone;}

  TH1D* GetptBkg() {return hptBkg;}
  TH1D* GetptBkgKt() {return hptBkgKt;}
  TH1D* GetptBkgAntiKt() {return hptBkgAntiKt;}
  TH1D* GetptBkgSISCone() {return hptBkgSISCone;}
 
  void SetXiBkgKt(TH1D *hBkg) {hXiBkgKt=(TH1D*) hBkg->Clone();}
  void SetXiBkgAntiKt(TH1D *hBkg) {hXiBkgAntiKt=(TH1D*) hBkg->Clone();}
  void SetXiBkgSISCone(TH1D *hBkg) {hXiBkgSISCone=(TH1D*) hBkg->Clone();}

  void SetZBkgKt(TH1D *hBkg) {hZBkgKt=(TH1D*) hBkg->Clone();}
  void SetZBkgAntiKt(TH1D *hBkg) {hZBkgAntiKt=(TH1D*) hBkg->Clone();}
  void SetZBkgSISCone(TH1D *hBkg) {hZBkgSISCone=(TH1D*) hBkg->Clone();}

  void SetptBkgKt(TH1D *hBkg) {hptBkgKt=(TH1D*) hBkg->Clone();hptBkgKt->SetDirectory(0);}
  void SetptBkgAntiKt(TH1D *hBkg) {hptBkgAntiKt=(TH1D*) hBkg->Clone();hptBkgAntiKt->SetDirectory(0);}
  void SetptBkgSISCone(TH1D *hBkg) {hptBkgSISCone=(TH1D*) hBkg->Clone();hptBkgSISCone->SetDirectory(0);}

  TObjArray* GetKtJets() {return ktList;}
  TObjArray* GetKtBFJets() {return ktBFList;}
  TObjArray* GetConeJets() {return coneList;}
  TObjArray* GetAllJets() {return allList;}
  TObjArray* GetPartons() {return PartonList;}

  TObjArray* GetFJktJets() {return fjKt;}
  TObjArray* GetFJAntiKtJets() {return fjAntiKt;}
  TObjArray* GetFJSISConeJets() {return fjSISCone;}
  
  TObjArray* GetUEEvents() {return ueEvents;}


  ktMuJet* GetKtJet(Int_t n) {return (ktMuJet*) ktList->At(n);}
  ktMuJet* GetKtBFJet(Int_t n) {return (ktMuJet*) ktBFList->At(n);}
  ktMuJet* GetConeJet(Int_t n) {return (ktMuJet*) coneList->At(n);}
  ktMuJet* GetAllJet(Int_t n) {return (ktMuJet*) allList->At(n);}
  ktParton* GetParton(int n) {return (ktParton*) PartonList->At(n);}

  ktMuFastJet* GetFJkt(Int_t n) {return (ktMuFastJet*) fjKt->At(n);}
  ktMuFastJet* GetFJAntiKt(Int_t n) {return (ktMuFastJet*) fjAntiKt->At(n);}
  ktMuFastJet* GetFJSISCone(Int_t n) {return (ktMuFastJet*) fjSISCone->At(n);}
  ktUEEvent* GetUEEvent(Int_t n){return (ktUEEvent*) ueEvents->At(n);}

  // only L0 infos so far ... (extend if L2 needed)
  Double_t GetHTEta();
  Double_t GetHTPhi(); // converted to 0-2pi
  Double_t GetJPEta();
  Double_t GetJPPhi(); // 0-2pi

  void SetEventNumber(Int_t mEvNum) {evNum = mEvNum;}
  //void SetTrigger(Int_t mtrigger) {trigger = mtrigger;}

  void SetBkgPtPerCell(Double_t mBkgPtPerCell) {BkgPtPerCell=mBkgPtPerCell;}
  void SetBkgPtPerCellAll(Double_t mBkgPtPerCellAll) {BkgPtPerCellAll=mBkgPtPerCellAll;}
  void SetBkgPtPerCellIn(Double_t mBkgPtPerCellIn) {BkgPtPerCellIn=mBkgPtPerCellIn;}
  void SetBkgPtPerCellAllIn(Double_t mBkgPtPerCellAllIn) {BkgPtPerCellAllIn=mBkgPtPerCellAllIn;}
  void SetBkgPtPerCellOut(Double_t mBkgPtPerCellOut) {BkgPtPerCellOut=mBkgPtPerCellOut;}
  void SetBkgPtPerCellAllOut(Double_t mBkgPtPerCellAllOut) {BkgPtPerCellAllOut=mBkgPtPerCellAllOut;}
  void SetBkgPtCone(Double_t mBkgPtCone) {BkgPtCone=mBkgPtCone;}
  void SetBkgPtCluster(Double_t mBkgPtCluster) {BkgPtCluster=mBkgPtCluster;}
  void SetBkgNCellCone(Double_t mBkgNCellCone) {BkgNCellCone=mBkgNCellCone;}
  void SetBkgNCellCluster(Double_t mBkgNCellCluster) {BkgNCellCluster=mBkgNCellCluster;}
  void SetBkgNCellIn(Double_t mBkgNCellIn) {BkgNCellIn=mBkgNCellIn;}
  void SetBkgNCellOut(Double_t mBkgNCellOut) {BkgNCellOut=mBkgNCellOut;}
  void SetBkgNJetCellIn(Double_t mBkgNJetCellIn) {BkgNJetCellIn=mBkgNJetCellIn;}
  void SetBkgNJetCellOut(Double_t mBkgNJetCellOut) {BkgNJetCellOut=mBkgNJetCellOut;}
  void SetBkgNCellInAll(Double_t mBkgNCellInAll) {BkgNCellInAll=mBkgNCellInAll;}
  void SetBkgNCellOutAll(Double_t mBkgNCellOutAll) {BkgNCellOutAll=mBkgNCellOutAll;}
  void SetNCellsInJetArea(Double_t mNCellsInJetArea) {NCellsInJetArea=mNCellsInJetArea;}

  void SetMedianPtPerArea(Double_t m_median_pt_per_area) {median_pt_per_area=m_median_pt_per_area;}

  void SetRefMult(Int_t mRefMult) {refMult=mRefMult;}
  void SetRPAngle(Double_t mPsiRP) {PsiRP=mPsiRP;}
  void SetTrigId(Int_t mTrigId) {TrigId=mTrigId;}
  void SetVertexZ(Double_t mvertexZ) {vertexZ=mvertexZ;}
  void SetTriggerInfo(TClonesArray* trinfo);
  void SetBkgPtCut(Double_t mptCut) {jfPt=mptCut;}
  void SetRParam(Double_t mRp) {jfRc=mRp;}

  void FillBkgInfo(ktGrid *mGrid);

  void AddParton(ktParton *p) {PartonList->AddLast(p);} // just for backward comp.
  void AddParton(Double_t px,Double_t py,Double_t pz, Double_t e);

  void AddFastJet(Double_t mpt,Double_t meta,Double_t mphi,Double_t metaCorr,Double_t mphiCorr,Double_t mMedian, Double_t marea,TString mtype, Double_t mNEF,  Double_t mNEFCorr, Double_t mAerror);
  //void AddFastJet(Double_t mpt,Double_t meta,Double_t mphi,Double_t metaCorr,Double_t mphiCorr,Double_t mMedian, Double_t marea,TString mtype, Double_t mNEF, Double_t mAerror);

  void AddUEEvent( TObjArray *mJets, ktGrid *grid);
  void AddUEEvent( TObjArray *mJets, ktFastJet *fJet);

  void EventSummary();
  void PrintFoundJets();
  void PrintTriggerInfo();
  void PrintFastJets(Double_t minPtCorr=10.0);
  void PrintJets(Double_t minPtCorr=10.0);

  void FillJetFinderSettings(ktGrid *mGrid);
  void PrintJetFinderSettings();
  void SetParticlePid(TString mPid) {jfPid += mPid;}
  TString GetJetFinderSettings();

  void FillAll(TObjArray *mJets);
  void Fill(TObjArray *mJets);
  void Fill(TObjArray *mJets,Double_t cellIn,Double_t cellOut);
  void Fill(TObjArray *mJets,Double_t cellIn,Double_t cellOut, Double_t coneIn);
  void Fill(TObjArray *mJets,Double_t cellIn,Double_t cellOut, Double_t coneIn, Double_t clusterIn);

 ClassDef(ktMuEvent,8)
};

#endif
