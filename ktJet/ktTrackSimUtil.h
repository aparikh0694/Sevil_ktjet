#include "TF1.h"

// hardcoded for STAR resolution (fix later !)
TF1 *ptRes=new TF1("ptRes","0.0051+0.0111*x",0,50);
TF1 *eRes=new TF1("eRes","0.15/sqrt(x)+0.03",0,50);

// STAR residual space charge distortions ...
Double_t ptSpaceCharge(Double_t pt,Int_t sign)
{
  if (sign>0)
    return pt+pt*pt*0.00086;
  else
    return pt-pt*pt*0.00086;
}

Double_t ptSpaceCharge(Double_t pt,Int_t sign,Double_t A)
{
  if (sign>0)
    return pt+pt*pt*A;
  else
    return pt-pt*pt*A;
}
