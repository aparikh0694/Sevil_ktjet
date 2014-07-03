// first test (Joern Putschke)

#ifndef ROOT_ktNN
#define ROOT_ktNN

#include "TObject.h"
#include "TBuffer.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TLorentzVector.h"
//#include "ktJetCell.h"
//#include "ktJet.h"
#include "ktGrid.h"

class ktGrid;
class ktJet;
class ktJetCell;

class ktNN : public TObject 
{

 private:

 ktGrid *NNGrid;
 //ktJetCell* NNCell;
 TObjArray *NNCells;

 public:

  ktNN();
  void Init(ktGrid *myNNGrid);
  Bool_t FindNNAroundSeed(ktJetCell *SeedCell);
  Bool_t FindNN(ktJet *NNJet); // will be a more sophisticated ansatz ;-) !
  Bool_t FindAll();

  //ktJetCell* GetNN() {return NNCell;}
  Int_t GetNNN() {return NNCells->GetEntriesFast();}
  TObjArray* GetNN() {return NNCells;}

  virtual ~ktNN();

  ClassDef(ktNN,1)
};

#endif
