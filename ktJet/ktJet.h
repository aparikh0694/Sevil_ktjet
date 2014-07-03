// first test (Joern Putschke)

#ifndef ROOT_ktJet
#define ROOT_ktJet

#include "TObject.h"
#include "TBuffer.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "ktJetCell.h"
//#include "ktGrid.h"

#include <Riostream.h>

class ktJet : public TObject
{

  protected:

  TObjArray *JetCells; 
  TObjArray *JetCellsFF; 
  TObjArray *JetParticles;
  TObjArray *JetParticlesFF;
  TLorentzVector *Jet;
  TString type;

  Int_t NJetParticles;
  Int_t NJetParticlesFF;
   
  Int_t NCellsInEmcal;
  Int_t NCellsOutsideEmcal;

  Int_t NIterations;
  Int_t NSubJets;

  Double_t FFAreaScale;
  Double_t JetAreaScale;

  public:

  ktJet();
  ktJet(ktJet *j);
  virtual ~ktJet();

  Bool_t IsSortable() const {return kTRUE;}

  Int_t Compare(const TObject *obj) const
  { 
    if (Pt() == ((ktJet*) obj)->Pt())
      return 0;
    else if (Pt() > ((ktJet*) obj)->Pt())
      return -1;
    else if (Pt() < ((ktJet*) obj)->Pt())
      return 1;
    else
      return -9999;
  }

  //TH2F *hJet;

  Int_t NJetCells() {return JetCells->GetEntriesFast();}
  Int_t NJetCellsInEmcal() {return NCellsInEmcal;}
  Int_t NJetCellsOutsideEmcal() {return NCellsOutsideEmcal;}
  Int_t NJetCellsFF() {return JetCellsFF->GetEntriesFast();}
  Int_t GetNJetParticles() {return JetParticles->GetEntriesFast();} //NJetParticles;}
  Int_t GetNJetParticlesFF() {return JetParticlesFF->GetEntriesFast();}
  Int_t GetNIterations() {return NIterations;}
  TObjArray* GetJetCells() {return JetCells;}
  TObjArray* GetJetCellsFF() {return JetCellsFF;}
  TObjArray* GetJetParticles() {return JetParticles;}
  TLorentzVector*  GetJetParticle(Int_t n) {return (TLorentzVector*) JetParticles->At(n);}
  ktJetCell* GetJetCell(Int_t n) {return (ktJetCell*) JetCells->At(n);}
  ktJetCell* GetJetCellFF(Int_t n) {return (ktJetCell*) JetCellsFF->At(n);}
  TLorentzVector* GetJetVector() {return Jet;}
  Int_t GetNSubJets() {return NSubJets;}

  Double_t Phi() const; // {return CellParticle->Phi();}
  Double_t Eta() const {return Jet->PseudoRapidity();}
  Double_t Pt() const {return Jet->Pt();}
  Double_t E() const {return Jet->E();}
  Double_t Kt2() {return Pt()*Pt();};
  
  Double_t GetJetAreaScale() {return JetAreaScale;}
  Double_t GetFFAreaScale() {return FFAreaScale;}

  TString GetType() {return type;}
 
  void Add(ktJet *j);
  void AddCell(ktJetCell *JetCell);
  void AddCell(ktJetCell *JetCell,Bool_t InEmcal);
  void AddCell(ktJetCell *JetCell,Bool_t InEmcal,Bool_t mSharing);
  void AddCellBkg(ktJetCell *JetCell);
  void AddCellFF(ktJetCell *JetCell);
  void CleanCellsInJet();
  void SetCellsInJet();
  void SetNIterations(Int_t m_NIterations) {NIterations=m_NIterations;}
  void SetNSubJets(Int_t m_NSubJets) {NSubJets=m_NSubJets;}
  void SetInCluster();
  void SetJetAreaScale(Double_t mJetAreaScale) {JetAreaScale=mJetAreaScale;}
  void SetFFAreaScale(Double_t mFFAreaScale) {FFAreaScale=mFFAreaScale;}

  // only for debug ...
  //TH1D* DrawXi();
  
  Double_t GetDistance(ktJet* mDjet);
  Double_t GetDistance(TLorentzVector* mDjet);

  //void Finish();
  void Finish();
  void ClearJet();
  void ResetJet();
  void Finish(Int_t Nphi, Double_t phiMin, Double_t phiMax, Int_t Neta, Double_t  etaMin,Double_t etaMax,Int_t myPhiBin, Int_t myEtaBin);
  void SetType(char* m_type) {type=m_type;}
  void PrintJet() const;
  //virtual void Print(const Option_t*)const {PrintJet();}

  ClassDef(ktJet,1)
};

#endif
