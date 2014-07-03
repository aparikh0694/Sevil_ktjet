// first test (Joern Putschke)

#ifndef ROOT_ktMuJet
#define ROOT_ktMuJet

#include "TObject.h"
#include "TBuffer.h"
#include "TObjArray.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "ktParticle.h"
#include "ktJet.h"
#include "TString.h"

class ktMuJet : public TObject
{

  protected:

  TObjArray *JetParticles;
  TLorentzVector *Jet;
  TString type;

  Int_t NJetParticles;
  Int_t NJetCells;

  Int_t NParticles;
  Int_t NCells;

  Int_t NCellsInEmcal;
  Int_t NCellsOutsideEmcal;

  Double_t JetPtCellCorr;
  Double_t JetPtConeCorr;
  Double_t JetPtClusterCorr;

  // if necessary save seed cell !
  Double_t mSeedEta;
  Double_t mSeedPhi;
  Double_t mSeedPt;
  Double_t NeutralEnergy;
  Double_t ChargedEnergy;

  Int_t NIterations;
  Int_t NSubJets;

  Double_t FFAreaScale;
  Double_t JetAreaScale;

  Bool_t electronSeed;

  public:

  ktMuJet();
  ktMuJet(ktJet *j);
  ktMuJet(ktMuJet *j,ktMuJet *j2);
  virtual ~ktMuJet();

  Bool_t IsSortable() const {return kTRUE;}

  Int_t Compare(const TObject *obj) const
  { 
    if (PtConeCorr() == ((ktMuJet*) obj)->PtConeCorr())
      return 0;
    else if (PtConeCorr() > ((ktMuJet*) obj)->PtConeCorr())
      return -1;
    else if (PtConeCorr() < ((ktMuJet*) obj)->PtConeCorr())
      return 1;
    else
      return -9999;
  }

  Int_t GetNJetParticles() {return JetParticles->GetEntriesFast();};//NJetParticles;}
  TObjArray* GetJetParticles() {return JetParticles;}
  ktParticle*  GetJetParticle(Int_t n) {return (ktParticle*) JetParticles->At(n);}
  TLorentzVector* GetJetVector() {return Jet;}

  Double_t Phi() const; // {return CellParticle->Phi();}
  Double_t Eta() const {return Jet->PseudoRapidity();}
  Double_t Pt() const {return Jet->Pt();}
  Double_t PtCellCorr() const {return JetPtCellCorr;}
  Double_t PtConeCorr() const {return JetPtConeCorr;}
  Double_t PtClusterCorr() {return JetPtClusterCorr;}
  Double_t E() const {return Jet->E();}
  Double_t Kt2() {return Pt()*Pt();}

  Double_t SeedEta() {return mSeedEta;}
  Double_t SeedPhi() {return mSeedPhi;}
  Double_t SeedPt() {return mSeedPt;}

  Double_t GetJetAreaScale() {return JetAreaScale;}
  Double_t GetFFAreaScale() {return FFAreaScale;}

  Int_t GetNCells() const {return NJetCells;}
  Int_t GetNCellsInEmcal() {return NCellsInEmcal;}
  Int_t GetNCellsOutsideEmcal() {return NCellsOutsideEmcal;} 
  Int_t GetNIterations() const {return NIterations;}
  TString GetType() {return type;}
  Int_t GetNSubJets() {return NSubJets;}
  Bool_t IsElectronSeed() {return electronSeed;}
  Double_t GetNeutralEnergy() {return NeutralEnergy;}
  Double_t GetChargedEnergy() {return ChargedEnergy;}

  Int_t GetNCellsInJet() {return NCells;}
  Int_t GetNParticlesInJet() {return NParticles;}

  Double_t GetDistance(ktMuJet* mDjet);

  void SetType(char* m_type) {type=m_type;}
  void SetNCells(Int_t mNJetCells) {NJetCells=mNJetCells;}
  void SetNIterations(Int_t m_NIterations) {NIterations=m_NIterations;}
  void SetPtCellCorr(Double_t cellIn, Double_t cellOut);
  void SetPtConeCorr(Double_t ptConeCorr);
  void SetPtClusterCorr(Double_t ptClusterCorr);
  void SetNSubJets(Int_t m_NSubJets) {NSubJets=m_NSubJets;}

  void SetJetAreaScale(Double_t mJetAreaScale) {JetAreaScale=mJetAreaScale;}
  void SetFFAreaScale(Double_t mFFAreaScale) {FFAreaScale=mFFAreaScale;}

  void PrintJet() const;
  void PrintSeed() const;

  ClassDef(ktMuJet,4)
};

#endif
