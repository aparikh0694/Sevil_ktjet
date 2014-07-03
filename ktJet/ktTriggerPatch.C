// first test (Joern Putschke)

// If necessary the grid cells can be stored in addition
// see ktJet.h and ktJet.C for further details how to implement

#include "ktTriggerPatch.h"
#include <Riostream.h>

ClassImp(ktTriggerPatch)

ktTriggerPatch::ktTriggerPatch()
{

  //ktGrid::ktGrid();

  patchPt=0.0;
  patchEtaBin=0;
  patchPhiBin=0;
  patchEta=0;
  patchPhi=0;
  NCellsX=0;
  NCellsY=0;

   // DEBUG:
  //cout<<"Default ktTriggerPatch constructor "<<endl;
}

ktTriggerPatch::ktTriggerPatch(Double_t mPt,Int_t mEtaBin,Int_t mPhiBin,Int_t mNCellsX,Int_t mNCellsY)
{

  patchPt=mPt;
  patchEtaBin=mEtaBin;
  patchPhiBin=mPhiBin;
  patchEta=0;
  patchPhi=0;
  NCellsX=mNCellsX;
  NCellsY=mNCellsY;

}

ktTriggerPatch::~ktTriggerPatch()
{

   // DEBUG:
  //cout<<"Default ktTriggerPatch destructor ..."<<endl;
}
