// first test (Joern Putschke)

#ifndef ROOT_ktPythia8
#define ROOT_ktPythia8

#include "TObject.h"
#include "TBuffer.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "Pythia.h"
#include "ktGrid.h"
#include "ktPy8Event.h"
#include "ktJetQuench.h"

using namespace Pythia8; 

class ktPythia8 : public TObject //, public Pythia
{

  private:
  
  Pythia *p;
  Double_t pyeta,pyphi,pyphi2,pypt,pye;
  Double_t mySqrts;

  Bool_t verbose;
  Bool_t smear;

  TF1 *ptRes;
  TF1 *eRes;
  
  public:

  ktPythia8();
  ktPythia8(char* xmlDir); 
  void SetRandomSeed(); // to be done, check Pythai8 examples ...
  void InitPhaseSpace(char* ptHatMin,char* ptHatMax);
  Bool_t ReadString(string mySetting) {return p->readString(mySetting,true);};
  void Init(TString mInit);
  void EventInfo() {p->info.list();};
  void HardProcessInfo();
  void ProcessList() {p->process.list();};
  void ParticleList() {p->event.list();};
  void Pycell();
  void Run();
  void RunWithPycell();
  void Fill(ktGrid* myGrid);
  void Fill(ktGrid* myGrid,TString mSel);
  void Fill(ktGrid* myGrid,TString mSel, Bool_t mSmear);
  void FillEvent(ktPy8Event *myEv);
  void Statistics() {p->statistics();};

  void GetHighestJetForQuenching(ktJetQuench *mQ=0,Double_t Rc=0.7,TString mSel="EMCAL", Bool_t mSmear=false);
  
  void SetVerbose(Bool_t mVerbose) {verbose=mVerbose;}
  void SetSmear(Bool_t mSmear) {smear=mSmear;}

  Double_t GetPycellEta(){return pyeta;};
  Double_t GetPycellPhi2(){return pyphi2;};
  Double_t GetPycellPhi(){return pyphi;};
  Double_t GetPycellPt(){return pypt;};
  Double_t GetPycellEt(){return pye;};
  Double_t GetPtHat() {return p->info.pTHat();};
  Double_t GetPhiHat() {return p->info.phiHat();};
  Double_t GetThetaHat() {return p->info.thetaHat();};
  Double_t GetSigmaGen() {return p->info.sigmaGen();};
  Bool_t GetSmear() {return smear;}

  Pythia* GetPythia() {return p;};

  virtual ~ktPythia8();
 
  ClassDef(ktPythia8,0)
};

#endif
