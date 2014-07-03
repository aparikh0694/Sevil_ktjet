// first test (Joern Putschke)

#ifndef ROOT_ktTriggerPatch
#define ROOT_ktTriggerPatch

#include "TObject.h"

class ktTriggerPatch : public TObject
{

 private:

  Double_t patchPt;

  Int_t patchEtaBin;
  Int_t patchPhiBin;
  Int_t NCellsX;
  Int_t NCellsY;

  Double_t patchEta;
  Double_t patchPhi;

 public:

  ktTriggerPatch();
  ktTriggerPatch(Double_t mPt,Int_t mEtaBin,Int_t mPhiBin,Int_t mNCellsX,Int_t mNCellsY);

  virtual ~ktTriggerPatch();

  Bool_t IsSortable() const {return kTRUE;}

  Int_t Compare(const TObject *obj) const
  { 
    if (Pt() == ((ktTriggerPatch*) obj)->Pt())
      return 0;
    else if (Pt() > ((ktTriggerPatch*) obj)->Pt())
      return -1;
    else if (Pt() < ((ktTriggerPatch*) obj)->Pt())
      return 1;
    else
      return -9999;
  }

  void SetPatchEtaPhi(Double_t mEta, Double_t mPhi) {patchEta=mEta;patchPhi=mPhi;}

  Double_t Pt() const {return patchPt;}

  Int_t GetEtaBin() {return patchEtaBin;}
  Int_t GetPhiBin() {return patchPhiBin;}
  Int_t GetNCellsX() {return NCellsX;}
  Int_t GetNCellsY() {return NCellsY;}

  Double_t GetEta() {return patchEta;}
  Double_t GetPhi() {return patchPhi;}

  ClassDef(ktTriggerPatch,0)
};

#endif
