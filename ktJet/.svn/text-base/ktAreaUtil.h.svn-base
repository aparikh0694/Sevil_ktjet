// first test (Joern Putschke)
// check if all cases are computed correctly !!!!!
// some cases not fully covered !????
// works for symmetry in eta and if input in phi is 0..2pi

#include <Riostream.h>
#include "TMath.h"

Double_t piR4(Double_t mR)
{
  return TMath::Pi()*mR*mR/4.0;
}

Double_t AreaCircle(Double_t mR,Double_t x)
{
  Double_t mA=0.0;
  Double_t sqrtRX=sqrt(mR*mR-x*x);
  
  mA +=piR4(mR);

  if (x<mR && sqrtRX>0)
    mA += 0.5*(x*sqrtRX+mR*mR*atan(x/sqrtRX));
  else
    mA -= piR4(mR);
  
  return mA;
}

Double_t AAreaCircle(Double_t mR,Double_t x)
{
  Double_t mA=0.0;
  Double_t sqrtRX=sqrt(mR*mR-x*x);

  if (x<mR && sqrtRX>0)
    mA += 0.5*(x*sqrtRX+mR*mR*atan(x/sqrtRX));
  else
    mA -= piR4(mR);
  
  return mA;
}

Double_t EtaAccTrans(Double_t eta, Double_t etaMax)
{
  return etaMax-eta;
}

Double_t PhiAccTrans(Double_t phi, Double_t phiMax)
{
  return phiMax-phi;
}

Double_t IntersectionX(Double_t mR,Double_t etaAcc)
{
  Double_t mb=-(etaAcc*etaAcc)+mR*mR;
  if (mb>0)
    return sqrt(mb);
  else
    return -99;
}

Double_t AreaInAcc(Double_t mR,Double_t m_eta, Double_t m_phi, Double_t etaMax, Double_t phiMax)
{
  Double_t mAIn=0.0;

  Double_t eta=TMath::Abs(m_eta); // symmetric in eta ...
  Double_t phi=-99.0;       // transform 0-2pi ...

  if (m_phi<TMath::Pi())
    phi=2*TMath::Pi()-m_phi;
  else
    phi=m_phi;

  Double_t etaAcc=EtaAccTrans(eta,etaMax);
  Double_t phiAcc=PhiAccTrans(phi,phiMax);

  if ((TMath::Abs(eta)+mR)<=etaMax && (TMath::Abs(phi)+mR)<=phiMax)
    return TMath::Pi()*mR*mR;

  Double_t interX=IntersectionX(mR,etaAcc);

  if (interX>-90)
    {
      Double_t mA1=AreaCircle(mR,-interX);
      Double_t mA2=0.0;
      Double_t mA3=0.0;

      if (phiAcc<mR)
	{
	  if ((2*interX)>phiAcc)
	    {
	      mA2=etaAcc*(interX+phiAcc);
	      mA3=AreaCircle(mR,phiAcc);
	    }
	  else
	    {
	      Double_t mA11=AAreaCircle(mR,phiAcc)-AAreaCircle(mR,interX);
	      // DEBUG:
	      //cout<<mA1<<" "<<mA11<<endl;

	      mA1 +=mA11;
	      mA2=etaAcc*(2*interX);	    
	      mA3=AreaCircle(mR,phiAcc);
	    }
	}
      else
	{	  
	  mA1=2*mA1;
	  mA2=etaAcc*(2*interX);
	  mA3=TMath::Pi()*mR*mR/2.0;
	}
      
      // DEBUG:
      //cout<<mA1<<" "<<mA2<<" "<<mA3<<endl;

      mAIn=mA1+mA2+mA3;
    }
  else
    mAIn=2*AreaCircle(mR,phiAcc);

  return mAIn;
}
