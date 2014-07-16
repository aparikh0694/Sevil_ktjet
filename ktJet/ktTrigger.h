// first test (Joern Putschke)

#ifndef ROOT_ktTrigger
#define ROOT_ktTrigger

#include "ktGrid.h"
#include "ktTriggerPatch.h"

class ktTrigger : public ktGrid
{

 private:

  Int_t NCellsX;
  Int_t NCellsY;
  Int_t NStepX;
  Int_t NStepY;

  Double_t minPatchPt;

  TObjArray *TriggerPatches;

 public:

  ktTrigger();
  ktTrigger(Int_t mNCellsX,Int_t mNCellsY,Int_t mNStepX, Int_t mNStepY,Double_t mMinPatchPt);
  virtual ~ktTrigger();

  Double_t GetPatchEnergy(Int_t mStartX,Int_t mStartY);
  Double_t GetPatchEnergyCenter(Int_t mStartX,Int_t mStartY);

  ktTriggerPatch* GetTriggerPatch();

  Double_t GetTriggerPt();
  Double_t GetTriggerEta();
  Double_t GetTriggerPhi();
  Double_t GetTriggerTreshold() {return minPatchPt;}

  Int_t GetNTriggerPatches() {return TriggerPatches->GetEntries();}

  void RunTrigger();

  void SetTriggerTreshold(Double_t  mMinPatchPt) {minPatchPt=mMinPatchPt;}
  void SetPatchSize(Int_t mNCellsX,Int_t mNCellsY) {NCellsX=mNCellsX;NCellsY=mNCellsY;}
  void SetStepSize(Int_t mNStepX,Int_t mNStepY) {NStepX=mNStepX;NStepY=mNStepY;}

  void PrintTriggerPatches();
  void PrintTriggerInfo();
  
  ClassDef(ktTrigger,0)
};

#endif
