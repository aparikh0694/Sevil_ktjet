// first test (Joern Putschke)

#ifndef ROOT_ktJetCell
#define ROOT_ktJetCell

#include "TObject.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "ktPID.h"

class ktJetCell : public TObject 
{

 protected:

  TObjArray *ParticleList;
  TObjArray *FFParticleList;

  TObjArray *PIDParticleList;
  TObjArray *PIDFFParticleList;

  //TClonesArray *ParticleList;
  //TClonesArray ParticleList;
  TLorentzVector *CellParticle; // sum of all particles in cell
  TLorentzVector *CellParticleFF;
  Int_t etaBin,phiBin;
  Bool_t inJet;
  Bool_t isTrack; // just for completeness 
  Bool_t isNN;
  Bool_t inEmcal;
  Bool_t inCluster;
  Bool_t isFF;
  Int_t NNIndex; // ?
  //Double_t kt2;

 public:

  ktJetCell();
  virtual ~ktJetCell();

  Bool_t IsSortable() const {return kTRUE;}

  Int_t Compare(const TObject *obj) const //!????
  { 
    
    if (Pt() == ((ktJetCell*) obj)->Pt())
      return 0;
    else if (Pt() > ((ktJetCell*) obj)->Pt())
      return -1;
    else if (Pt() < ((ktJetCell*) obj)->Pt())
      return 1;
    else
      return -9999;
    
  }

  void AddParticle(TLorentzVector *mPart);
  void AddFFParticle(TLorentzVector *mPart);
  void AddParticlePID(ktPID *mPID);
  void AddFFParticlePID(ktPID *mPID);

  // Getter

  TObjArray* GetParticleList() {return ParticleList;}
  TObjArray* GetFFParticleList() {return FFParticleList;}
  TObjArray* GetParticleListPID() {return PIDParticleList;}
  TObjArray* GetFFParticleListPID() {return PIDFFParticleList;}
  
  Int_t GetEtaBin() {return etaBin;}
  Int_t GetPhiBin() {return phiBin;}
  Bool_t InJet() {return inJet;}
  Bool_t InCluster() {return inCluster;}
  Bool_t InEmcal() {return inEmcal;}
  Bool_t IsNN() {return isNN;}
  Bool_t IsFF() {return isFF;}
  Bool_t GetIsTrack() {return isTrack;}
  TLorentzVector* GetCellParticle() {return CellParticle;}
  TLorentzVector* GetCellParticleFF() {return CellParticleFF;}
  void PrintCell();
  //TClonesArray* GetParticleList() {return ParticleList;};
  void Clean();
  
  Int_t NParticles() {return ParticleList->GetEntriesFast();}
  Int_t NFFParticles() {return FFParticleList->GetEntriesFast();}
  Int_t NParticlesPID(){return PIDParticleList->GetEntriesFast();}

  Double_t Phi(); // {return CellParticle->Phi();}
  Double_t Eta();// {return CellParticle->PseudoRapidity();}
  Double_t PhiFF(); // {return CellParticle->Phi();}
  Double_t EtaFF() {return CellParticleFF->PseudoRapidity();}
  Double_t Pt() const  {return CellParticle->Pt();}
  Double_t E()   {return CellParticle->E();}
  Double_t Kt2() {return Pt()*Pt();};

  // Setter

  void SetEtaBin(Int_t m_etaBin) {etaBin=m_etaBin;}
  void SetPhiBin(Int_t m_phiBin) {phiBin=m_phiBin;}
  void SetInJet(Bool_t m_InJet) {inJet=m_InJet;}
  void SetInCluster(Bool_t m_InCluster) {inCluster=m_InCluster;}
  //void SetKt2(Double_t m_kt2) {kt2=m_kt2;}
  void SetIsTrack(Bool_t m_isTrack) {isTrack=m_isTrack;}
  void SetIsNN(Bool_t m_isNN) {isNN=m_isNN;}
  void SetIsFF(Bool_t m_isFF) {isFF=m_isFF;}
  void SetInEmcal(Bool_t m_inEmcal) {inEmcal=m_inEmcal;}

  ClassDef(ktJetCell,1)
};

#endif
