// first test (Joern Putschke)

#ifndef ROOT_ktParton
#define ROOT_ktParton

#include "TObject.h"
#include "TBuffer.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TString.h"

class ktParton : public TObject
{

  protected:

  TLorentzVector *Parton;
  TObjArray *PartonParticles;

  public:

  ktParton();
  ktParton(TLorentzVector *v);
  ktParton(Double_t px,Double_t py,Double_t pz, Double_t e);

  virtual ~ktParton();

  TLorentzVector* GetParton() {return Parton;}

  Double_t Phi() const; // {return CellParticle->Phi();}
  Double_t Eta() const {return Parton->PseudoRapidity();}
  Double_t Pt() const {return Parton->Pt();}
  Double_t E() const {return Parton->E();}

  void SetParton(Double_t px,Double_t py,Double_t pz, Double_t e);
  void SetParton(TLorentzVector *v);

  void AddParticle(TLorentzVector *v);
  void PrintParton();

  ClassDef(ktParton,1)
};

#endif
