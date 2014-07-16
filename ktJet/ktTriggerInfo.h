#ifndef __KTTRIGGERINFO_HH
#define __KTTRIGGERINFO_HH

#include <TObject.h>
#include <TString.h>

class ktTriggerInfo : public TObject
{
 public:
  ktTriggerInfo();
  ktTriggerInfo(const ktTriggerInfo &t);
  ktTriggerInfo(ktTriggerInfo* t);
 
  // ktTriggerInfo(TString trigtype,Int_t flag,Float_t eta, Float_t phi);
  virtual ~ktTriggerInfo();

  void    Clear(Option_t *Option = "");

  //Setters
  // void SetTriggerType(TString trigtype){fTrigType=trigtype;}
  void SetEta(Float_t eta) {fEta=eta;}
  void SetPhi(Float_t phi) {fPhi=phi;}
  void SetTriggerFlag(Int_t flag) {fTrigFlag=flag;}

  //Getters
  //TString GetTriggerType() const {return fTrigType;}
  Float_t GetEta() const {return fEta;}
  Float_t GetPhi() const {return fPhi;}
  Int_t GetTriggerFlag() const {return fTrigFlag;}
  Int_t isJPL0();
  Int_t isJPL2();
  Int_t isHTL0();
  Int_t isHTL2();
  Int_t isBBC();
  void PrintInfo();


 private:
  //TString fTrigType; //JetPatchL0,JetPatchL2,HighTowerL0,HighTowerL2
  Float_t fEta; //eta of the tower/patch in the reference frame 0,0,0 
  Float_t fPhi; //phi of the tower/patch in the reference frame 0,0,0 
  Int_t fTrigFlag;//1=HTL0, 2=JPL0, 3=HTL2, 4=JPL2

 ClassDef(ktTriggerInfo, 1)
   };

#endif
