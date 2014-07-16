
// first test (Joern Putschke)

#ifndef ROOT_ktGrid
#define ROOT_ktGrid

#include "TObject.h"
//#include "TNamed.h"
#include "TBuffer.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "ktJetCell.h"
#include "ktJet.h"
#include "ktNN.h"
#include "ktPseudoJet.h"
#include "ktPID.h"

//class ktNN;

struct Index{Int_t iPhi;Int_t iEta;};

class ktGrid : public TObject 
{

 private:

  Double_t phiMin;
  Double_t phiMax;
  Double_t etaMin;
  Double_t etaMax;

  Double_t phiMinEmcal;
  Double_t phiMaxEmcal;
  Double_t etaMinEmcal;
  Double_t etaMaxEmcal;

  Bool_t EmcalSet;
  Bool_t SeedAccCut;

  Int_t Neta;
  Int_t Nphi;

  Int_t NMaxIterations;

  Double_t Seed;
  Double_t ClSeed;
  Int_t NSeeds;
  Int_t NBkgSeeds;

  Double_t RMaxNN;
  Double_t Rc;
  Double_t RcCl;
  Double_t RcClMax;
  Double_t RcB;
  Double_t Rff;
  Double_t BkgPtCut;
  Double_t minJetEnergy;

  Double_t BkgEnergyPerCell;
  Double_t BkgPtPerCell;
  Double_t BkgPtPerCellAll;
  Double_t BkgPtPerCellIn;
  Double_t BkgPtPerCellAllIn;
  Double_t BkgPtPerCellOut;
  Double_t BkgPtPerCellAllOut;
  Double_t BkgPtCone;
  Double_t BkgPtCluster;
  Double_t BkgNCellCone;
  Double_t BkgNCellCluster;

  Double_t BkgNCellIn;
  Double_t BkgNCellOut;
  Double_t BkgNCellInAll;
  Double_t BkgNCellOutAll;
  Double_t BkgNJetCellIn;
  Double_t BkgNJetCellOut;

  Bool_t RcBset;
  Bool_t Sharing;
  Bool_t XiBkg;
  
  Bool_t verbose;

  Int_t MaxNCellsInJet;
  Int_t MaxNCellsInFF;

  /*
  Int_t PhiBin(Double_t m_phi) {return (hGrid->GetXaxis()->FindBin((Double_t) m_phi)-1);}
  Int_t EtaBin(Double_t m_eta) {return (hGrid->GetYaxis()->FindBin((Double_t) m_eta)-1);}
  Double_t PhiFromBin(Int_t m_phiBin) {return (hGrid->GetXaxis()->GetBinCenter(m_phiBin+1));}
  Double_t EtaFromBin(Int_t m_etaBin) {return (hGrid->GetYaxis()->GetBinCenter(m_etaBin+1));}
  */

  ktJetCell **myGrid;
  //ktJetCell *mySeeds[100];
  TObjArray *mySeedList;
  //ktJetCell *mySeedsSort[100];
  TObjArray *Jets;
  TObjArray *ClusterJets; // only for debug reasons !!!!

  Index GridIndex[100000]; // just for test reasons (define later); 
                           // higher value needed for HIJING ;-) -> change to dynamically !!
                           // or use to speed up print grid and clear grid + bkg !

 public:

  TH2F *hGrid;
  TH1D *hXiBkg;

  Int_t NCells;

  ktGrid();
  virtual ~ktGrid();

  Int_t PhiBin(Double_t m_phi) {return (hGrid->GetXaxis()->FindBin((Double_t) m_phi)-1);}
  Int_t EtaBin(Double_t m_eta) {return (hGrid->GetYaxis()->FindBin((Double_t) m_eta)-1);}
  Double_t PhiFromBin(Int_t m_phiBin) {return (hGrid->GetXaxis()->GetBinCenter(m_phiBin+1));}
  Double_t EtaFromBin(Int_t m_etaBin) {return (hGrid->GetYaxis()->GetBinCenter(m_etaBin+1));}

  void Fill(TLorentzVector *mPart);
  void Fill(TLorentzVector *mPart, Bool_t isCharged);
  void Fill(TLorentzVector *mPart, ktPID *mPID);

  void DoJetfinding(TString option);

  void FindKtJets();
  void FindKtJetsBF();
  void FindConeJetsBkg();
  void FindConeJets();
  void FindConeCluster();

  void AddFFParticles(ktJet* mJ);

  void CleanGridAndCalcBkg();
  void CleanGrid();
  void CalcBkg();
  void CalcBkgNJetsRemove(Int_t nJR=1);
  void CalcBkgRemoveDijet();
  void BkgRandomCones();
  void BkgRandomConesFast();
  void BkgXi();
  void BkgXiScaled();

  // Getter

  void PrintIndex();
  void PrintGrid();
  void PrintSeeds();
  void PrintJets();
  void PrintJetsNice();//Float_t MinJetPtCut=0);
  void PrintJetsNiceBkgSub();
  //void PrintSeedsSort();

  void GridInfo();
  Bool_t GetSeeds();
  Bool_t GetSeeds(Double_t mRdistance);
  void GetSeeds(TObjArray* clSeedList, Double_t etaCenter, Double_t phiCenter, Double_t mRcCl);
  void GetSeeds(TObjArray* clSeedList,ktJetCell *bkgSeedCell, Double_t mRcCl);

  TObjArray* GetJetList() {return Jets;};

  // only debug reasons !!!
  TObjArray* GetClusterJetList() {return ClusterJets;}

  //ktJet* GetSeedJet(Double_t etaCenter, Double_t phiCenter);
  void GetSeedJet(ktJet* seedJet,Double_t etaCenter, Double_t phiCenter);

  // calc in constructor and return value ...
  Double_t CelldEta() {return (Double_t) TMath::Abs(etaMax-etaMin)/Neta;}
  Double_t CelldPhi() {return (Double_t) TMath::Abs(phiMax-phiMin)/Nphi;} 
  Double_t GetRc() {return Rc;}
  Double_t GetRcCl() {return RcCl;}
  Double_t GetRcClMax() {return RcClMax;}
  Double_t GetRff() {return Rff;}
  Double_t GetBkgPtCut() {return BkgPtCut;}
  Double_t GetNCells() {return NCells;}
  Double_t GetSeed() {return Seed;}
  Double_t GetClusterSeed() {return ClSeed;}
  Int_t GetNBkgSeeds() {return NBkgSeeds;}
  Double_t GetRMaxNN() {return RMaxNN;}
  ktJetCell* GetCell(Int_t mPhi,Int_t mEta) {return &myGrid[mPhi][mEta];}
  //ktJetCell* GetCell(Double_t mPhi,Double_t mEta) {return &myGrid[PhiBin(mPhi)][EtaBin(mEta)];}
  Int_t NJets() {return Jets->GetEntriesFast();}
  Int_t GetPhiConeBin() {return (Int_t) (Rc/CelldPhi())+1;}
  Int_t GetEtaConeBin() {return (Int_t) (Rc/CelldEta())+1;}
  Int_t GetPhiConeBinCl() {return (Int_t) (RcClMax/CelldPhi())+1;}
  Int_t GetEtaConeBinCl() {return (Int_t) (RcClMax/CelldEta())+1;}
  Int_t GetPhiConeBinBkg() {return (Int_t) (RcB/CelldPhi())+1;}
  Int_t GetEtaConeBinBkg() {return (Int_t) (RcB/CelldEta())+1;}
  Int_t GetPhiConeBinFF() {return (Int_t) (Rff/CelldPhi())+1;}
  Int_t GetEtaConeBinFF() {return (Int_t) (Rff/CelldEta())+1;}
  Int_t GetPhiNNConeBin() {return (Int_t) (RMaxNN/CelldPhi())+1;}
  Int_t GetEtaNNConeBin() {return (Int_t) (RMaxNN/CelldEta())+1;}
  Int_t GetEtaBinRc(Double_t mC) {return (Int_t) (mC/CelldEta())+1;}
  Int_t GetPhiBinRc(Double_t mC) {return (Int_t) (mC/CelldPhi())+1;}
  Double_t GetMinJetEnergy() {return minJetEnergy;}
  Double_t GetBkgEnergyPerCell() {return BkgEnergyPerCell;}
  Double_t GetBkgPtPerCell() {return BkgPtPerCell;}
  Double_t GetBkgPtPerCellAll() {return BkgPtPerCellAll;}
  Double_t GetBkgPtPerCellIn() {return BkgPtPerCellIn;}
  Double_t GetBkgPtPerCellAllIn() {return BkgPtPerCellAllIn;}
  Double_t GetBkgPtPerCellOut() {return BkgPtPerCellOut;}
  Double_t GetBkgPtPerCellAllOut() {return BkgPtPerCellAllOut;}  
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
  Bool_t GetSharing() {return Sharing;}
  Int_t GetNCellsInJet(Double_t mC);
  Int_t GetNCellsInJet(Double_t etaCenter, Double_t phiCenter,Double_t mC);
  Int_t GetNCellsInCluster();//TObjArray *mL);
  Double_t CalcJetArea(Double_t mC) {return TMath::Pi()*mC*mC;}
  Double_t CalcRcFromBins(Int_t mN) {return TMath::Sqrt(mN*CelldEta()*CelldPhi()/TMath::Pi());}
  Double_t CalcBinsInJetArea(Double_t mC) {return TMath::Pi()*mC*mC/(CelldEta()*CelldPhi());}
  Double_t CalcScaleRatio(Double_t mC1,Double_t mC2) {return mC1*mC1/(mC2*mC2);}
  Bool_t GetEmcalSet() {return EmcalSet;}
  Double_t GetBkgCone() {return RcB;}
  TH1D* GetXiBkg() {return hXiBkg;}

  Bool_t InEmcal(Double_t mPhi,Double_t mEta);
  
  Int_t GetNeta() {return Neta;}
  Int_t GetNphi() {return Nphi;}
  Int_t GetNMaxIterations() {return NMaxIterations;}

  // Setter

  void SetGrid(Int_t m_Nphi,Double_t m_phiMin,Double_t m_phiMax,Int_t m_Neta,Double_t m_etaMin,Double_t m_etaMax);
  void SetEmcal(Double_t m_phiMin,Double_t m_phiMax,Double_t m_etaMin,Double_t m_etaMax);
  void SetSeed(Double_t m_Seed) {Seed=m_Seed;}
  void SetClusterSeed(Double_t m_ClSeed) {ClSeed=m_ClSeed;}
  void SetRMaxNN(Double_t m_RMaxNN) {RMaxNN=m_RMaxNN;}
  void SetCone(Double_t m_Rc) {Rc=m_Rc;}
  void SetClusterCone(Double_t m_RcCl) {RcCl=m_RcCl;}
  void SetClusterConeMax(Double_t m_RcClMax) {RcClMax=m_RcClMax;}
  void SetBkgCone(Double_t m_RcB) {RcB=m_RcB;RcBset=true;}
  void SetFFCone(Double_t m_Rff) {Rff=m_Rff;}
  void SetBkgPtCut(Double_t m_BkgPtCut) {BkgPtCut=m_BkgPtCut;}
  void SetMinJetEnergy(Double_t m_minJetEnergy) {minJetEnergy=m_minJetEnergy;}
  //void SetBkgEnergyPerCell(Double_t m_BkgEnergyPerCell) {BkgEnergyPerCell=m_BkgEnergyPerCell;} // Allow only once one have average Bkg calculations
  //void SetBkgtPterCell(Double_t m_BkgPtPerCell) {BkgPtPerCell=m_BkgPtPerCell;}
  void SetNBkgSeeds(Int_t m_NBkgSeeds) {NBkgSeeds=m_NBkgSeeds;}
  void SetNMaxIterations(Int_t m_NMaxIterations) {NMaxIterations=m_NMaxIterations;}
  void SetSharing(Bool_t mSharing) {Sharing=mSharing;}
  void SetXiBkg(Bool_t mXiBkg) {XiBkg=mXiBkg;} // bad fix because of hist keys !!!????
  void SetVerbose(Bool_t mVerbose) {verbose=mVerbose;}
  void Rebin(Double_t mRebin);
  void Init();

  ClassDef(ktGrid,0)
};

#endif
