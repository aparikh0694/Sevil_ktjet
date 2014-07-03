// first test (Joern Putschke)

#ifndef ROOT_ktPseudoJet
#define ROOT_ktPseudoJet

#include "TObject.h"
#include "TBuffer.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "ktJetCell.h"
#include "ktJet.h"

class ktPseudoJet : public ktJet
{

  private:

  public:

  ktPseudoJet();
  ktPseudoJet(ktJetCell* mJCell);
  //virtual ~ktPseudoJet();
  ~ktPseudoJet();

  void Merge(ktPseudoJet *mMerge);

  ClassDef(ktPseudoJet,1)
};

#endif

