// first test (Joern Putschke)

#ifndef ROOT_ktHijing
#define ROOT_ktHijing

#include "THijing.h"
#include "TH1.h"
#include "ktFastJet.h"
#include "ktGrid.h" 

class ktHijing : public THijing
{

  private:

  TClonesArray* particles;
  Int_t nEvents;
  Float_t etaMin,etaMax;
  Bool_t etaRange;
  
  Bool_t CheckJetEtaRange();

  TH1F *etaP,*ptP;

  public:

  ktHijing(); //dummy constructor no effect ...
  ktHijing(const Char_t* name, const Char_t* title);

  void SetJetPtRange(Float_t minJetPt=15.0, Float_t maxJetPt=-1.0);
  void SetJetEtaRange(Float_t mEtaMin=-1.0, Float_t mEtaMax=1.0); 
  void SetISRandFSR(int mode=3);

  void InitRHICdAu(); // Add more or use THijing::Initialize(...)
  void InitRHICAuAu();

  UInt_t MakeSeed(int mode=0);
  void MakeMinbiasEvent();
  void FillFastJet(ktFastJet *fj=0, Int_t mSel=2);
  void FillGrid(ktGrid *mGrid=0, Int_t mSel=2); // Add filling LOCone later ...

  Int_t GetNParticles() {return particles->GetEntries();}
  Int_t GetNEvents() {return nEvents;}

  Float_t GetPtFirstParton();
  Float_t GetEtaFirstParton();
  Float_t GetPtFirstPartonFSR();
  Float_t GetEtaFirstPartonFSR();

  Float_t GetPtSecondParton();
  Float_t GetPtSecondPartonFSR();

  TH1F* GetEtaHist() {return etaP;}
  TH1F* GetPtHist() {return ptP;}

  virtual ~ktHijing();
   
  ClassDef(ktHijing,0)
};

#endif

  
