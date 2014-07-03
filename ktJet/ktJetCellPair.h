// first test (Joern Putschke)

#ifndef ROOT_ktJetCellPair
#define ROOT_ktJetCellPair

#include "TObject.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "ktJetCell.h"
#include "ktPseudoJet.h"

class ktJetCellPair : public TObject //public ktJetCell // ???
{

 private:

  ktJetCell *A;
  ktJetCell *B;
  //ktJetCell *M;
  ktPseudoJet *AJ;
  ktPseudoJet *BJ;

 public:

  Double_t dij;

  ktJetCellPair();
  ktJetCellPair(ktJetCell *m_A,ktJetCell *m_B);
  ktJetCellPair(ktPseudoJet *m_AJ,ktPseudoJet *m_BJ);

  //ktJetCell* Merge() {return M;}
  ktJetCell* Merge();
  ktJetCell* GetA() {return A;}
  ktJetCell* GetB() {return B;}
  ktPseudoJet* GetAJ() {return AJ;}
  ktPseudoJet* GetBJ() {return BJ;}

  virtual ~ktJetCellPair();

  // Getter

  // Setter

 ClassDef(ktJetCellPair,1)
};

#endif
