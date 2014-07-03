// first test (Joern Putschke)

#ifndef ROOT_ktMuFastJet
#define ROOT_ktMuFastJet

#include "TObject.h"
#include "TBuffer.h"
#include "TObjArray.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "ktParticle.h"
#include "TH1.h"

class ktMuFastJet : public TObject
{

  protected:

  TObjArray *JetParticles;

  Double_t pt,eta,phi,etaCorr,phiCorr;
  Double_t median_pt_per_area;
  Double_t median_pt_per_area_charged;
  Double_t median_pt_per_area_neutral;
  Double_t area;
  Double_t areaError;
  Double_t FFAreaInAcc;
  Double_t NEF;
  Double_t NEFCorr;
  Bool_t isTrigger;

  TString type;
  TH1D *hpt;

  public:
  
  ktMuFastJet();
  ktMuFastJet(Double_t mpt,Double_t meta,Double_t mphi,Double_t mMedian, Double_t marea);
  ktMuFastJet(Double_t mpt,Double_t meta,Double_t mphi,Double_t metaCorr,Double_t mphiCorr,Double_t mMedian, Double_t marea);
  ktMuFastJet(Double_t mpt,Double_t meta,Double_t mphi,Double_t mMedian, Double_t marea,TString mtype);
  ktMuFastJet(Double_t mpt,Double_t meta,Double_t mphi,Double_t metaCorr,Double_t mphiCorr,Double_t mMedian, Double_t marea,TString mtype);

  virtual ~ktMuFastJet();

  void PrintJet();

  Bool_t IsSortable() const {return kTRUE;}

  Int_t Compare(const TObject *obj) const
  { 
    if (PtCorr() == ((ktMuFastJet*) obj)->PtCorr())
      return 0;
    else if (PtCorr() > ((ktMuFastJet*) obj)->PtCorr())
      return -1;
    else if (PtCorr() < ((ktMuFastJet*) obj)->PtCorr())
      return 1;
    else
      return -9999;
  }

  Double_t PtCorr() const {return pt-area*median_pt_per_area;}
  Double_t Pt() {return pt;}
  Double_t Area() {return area;}
  Double_t AreaError() {return areaError;}
  Double_t MedianPtPerArea() {return median_pt_per_area;}

  Double_t MedianPtPerAreaCharged() {return median_pt_per_area_charged;} // not filled needed ?
  Double_t MedianPtPerAreaNeutral() {return median_pt_per_area_neutral;} // not filled needed ?

  Double_t Eta() {return eta;}
  Double_t Phi() {return phi;}
  Double_t EtaCorr() {return etaCorr;}
  Double_t PhiCorr() {return phiCorr;}
  Double_t GetFFAreaInAcc() {return FFAreaInAcc;}
  Double_t GetNEF() {return NEF;}
  Double_t GetNEFCorr() {return NEFCorr;}
  TString GetType() {return type;}
  TH1D* GetJetPtHistogram() {return hpt;}

  Bool_t IsTriggerJet() {return isTrigger;}

  Int_t GetNJetParticles() {return JetParticles->GetEntriesFast();};//NJetParticles;}
  TObjArray* GetJetParticles() {return JetParticles;}
  //TLorentzVector* GetJetParticle(Int_t n) {return (TLorentzVector*) JetParticles->At(n);}
  ktParticle*  GetJetParticle(Int_t n) {return (ktParticle*) JetParticles->At(n);}

  // add setters ...
  void SetType(char* m_type) {type=m_type;} 
  void SetMedianPtPerArea(Double_t m_median_pt_per_area) {median_pt_per_area=m_median_pt_per_area;}
  void SetMedianPtPerAreaCharged(Double_t m_median_pt_per_area_charged) {median_pt_per_area_charged=m_median_pt_per_area_charged;}
  void SetMedianPtPerAreaNeutral(Double_t m_median_pt_per_area_neutral) {median_pt_per_area_neutral=m_median_pt_per_area_neutral;}
  void SetAreaError(Double_t m_areaError) {areaError=m_areaError;}
  void SetNEF(Double_t m_NEF) {NEF=m_NEF;}
  void SetNEFCorr(Double_t m_NEFCorr) {NEFCorr=m_NEFCorr;}
  void SetJetPtHistogram(TH1D *mhpt) {hpt=(TH1D*) mhpt->Clone();hpt->SetDirectory(0);}

  void AddJetParticle(TLorentzVector *mP) {JetParticles->AddLast(mP->Clone());}
  void AddJetParticle(ktParticle *mP);

  void SetFFAreaInAcc(Double_t mFFAreaInAcc) {FFAreaInAcc=mFFAreaInAcc;}
  void SetIsTriggerJet(Bool_t mIsTrigger) {isTrigger=mIsTrigger;}
  // Add more infos if needed !!!!

  ClassDef(ktMuFastJet,7)
};

#endif
