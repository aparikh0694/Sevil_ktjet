#include "ktTriggerInfo.h"
#include "Riostream.h"
ClassImp(ktTriggerInfo)

//____________________________________________________
ktTriggerInfo::ktTriggerInfo()
 : TObject()
  // , fTrigType(0)
   , fEta(0)
   , fPhi(0)
 , fTrigFlag(0)
{}
//____________________________________________________
ktTriggerInfo::ktTriggerInfo(const ktTriggerInfo &t)
: TObject(t)
  //, fTrigType(t.fTrigType)
  , fEta(t.fEta)
  , fPhi(t.fPhi)
  , fTrigFlag(t.fTrigFlag)
{}

//____________________________________________________
ktTriggerInfo::ktTriggerInfo(ktTriggerInfo* t)
{
  fEta=t->GetEta();
  fPhi=t->GetPhi();
  fTrigFlag=t->GetTriggerFlag();
}

//____________________________________________________
// ktTriggerInfo::ktTriggerInfo(TString trigtype,Int_t flag,Float_t eta, Float_t phi){
//   fTrigType=trigtype;
//   fEta=eta;
//   fPhi=phi;
//   fTrigFlag=flag;
// }
//____________________________________________________
ktTriggerInfo::~ktTriggerInfo(){
  //destructor
  fEta=0;
  fPhi=0;
  fTrigFlag=0;
}
//____________________________________________________
Int_t ktTriggerInfo::isJPL0(){
  if(fTrigFlag==2)return 1;
  else return 0;
}
//____________________________________________________
Int_t ktTriggerInfo::isJPL2(){
  if(fTrigFlag==4)return 1;
  else return 0;
}
//____________________________________________________
Int_t ktTriggerInfo::isHTL0(){
  if(fTrigFlag==1)return 1;
  else return 0;
}
//____________________________________________________
Int_t ktTriggerInfo::isHTL2(){
  if(fTrigFlag==3)return 1;
  else return 0;
}
//____________________________________________________
Int_t ktTriggerInfo::isBBC(){
  if(fTrigFlag==5)return 1;
  else return 0;
}
//____________________________________________________
void ktTriggerInfo::Clear(Option_t */*Option*/){
  // fTrigType="";
  fEta=0;
  fPhi=0;
  fTrigFlag=0;

}
//____________________________________________________
void ktTriggerInfo::PrintInfo(){
  cout<<"=============== Trigger Info =============="<<endl;
  if(fTrigFlag==1)cout<<"HTL0;  eta="<<fEta<<"   phi="<<fPhi<<endl;
  if(fTrigFlag==2)cout<<"JPL0;  eta="<<fEta<<"   phi="<<fPhi<<endl;
  if(fTrigFlag==3)cout<<"HTL2;  eta="<<fEta<<"   phi="<<fPhi<<endl;
  if(fTrigFlag==4)cout<<"JPL2;  eta="<<fEta<<"   phi="<<fPhi<<endl;
  cout<<"=========================================="<<endl;
}
