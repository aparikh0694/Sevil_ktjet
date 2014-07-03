// first test (Joern Putschke)
#include "ktHijing.h"
#include <Riostream.h>
#include "TParticle.h"
#include "TDatime.h"
#include "TRandom.h"
#include "TRanMar.h"
#include "TSystem.h"
#include "TLorentzVector.h"

ClassImp(ktHijing)

// -------------------------------------------------------
// Selection 0=all , 1=charged only , 2=charged+neutral
Bool_t Select_TT(Int_t m_sel,TParticle *m_p)
{
  Bool_t result=false;
  Int_t m_PDG=m_p->GetPdgCode();
  TParticlePDG *m_PDG_p=(TParticlePDG*) m_p->GetPDG(1);

  if (m_sel!=0)
    {
      if (m_sel==1)
        {
          if (m_PDG_p->Charge() !=0)
            result=true;
          else
            result=false;
        }
      else if (m_sel==2)
        {
          if (m_PDG_p->Charge() !=0)
            result=true;
          else if ((m_PDG==111 || m_PDG==22))          
            result=true;
          else
            result=false;
        }
      else if (m_sel==3)
        {
          if (m_PDG_p->Charge() !=0)
            result=true;
          else if ((m_PDG==111 || m_PDG==22))        
            result=true;
          else
            result=false;
        }
      else
        result=true;  // default : m_sel = 0 or print error message
    }
  else
    result=true;

  // DEBUG:
  //cout<<m_PDG_p->GetName()<<" "<<result<<endl;

  return result;
}


ktHijing::ktHijing() : THijing()
{
  //DEBUG:
  cout<<endl;
  cout<<"Default constructor of ktHijing dummy, no effect ..."<<endl;
  particles = new TClonesArray("TParticle");
  nEvents=0;
  etaMin=-1.0;
  etaMax=1.0;
  etaRange=false;
  ptP=new TH1F("ptP","pt distribution first parton HIJING",80,0,20);
  etaP=new TH1F("etaP","eta distribution first parton HIJING",100,-5,5);
}

ktHijing::ktHijing(const Char_t* name, const Char_t* title)
  : THijing(name,title)
{
  //DBEUG:
  //cout<<endl;
  //cout<<name<<" "<<title<<endl;
  //cout<<"Default constructor of ktHijing ..."<<endl;
  particles = new TClonesArray("TParticle");
  nEvents=0;
  etaMin=-1.0;
  etaMax=1.0;
  etaRange=false;
  ptP=new TH1F("ptP","pt distribution first parton HIJING",80,0,20);
  etaP=new TH1F("etaP","eta distribution first parton HIJING",100,-5,5);
}

ktHijing::~ktHijing()
{
  //DEBUG:
  //cout<<endl;
  //cout<<"Default destructor of ktHijing ..."<<endl;
  particles->Delete();
  delete particles;

  delete ptP;delete etaP;
}

UInt_t ktHijing::MakeSeed(int mode)
{
  switch (mode) {
  case 1: 
    {
      TDatime date;
      return date.GetDate();
    }
  case 2:
    {
      TDatime date;
      UInt_t  seed1 = (date.GetDate() / 1000000 * 100 + 
		       date.GetTime() / 1000000);
      UInt_t  seed2 = (date.GetTime() % 10000);
      TRanMar ranmar(seed1, seed2);
      Int_t   eat = (date.GetDate() - 19980101 + seed2 +
		     gSystem->GetPid()); 
      for (Int_t i = 0; i < (eat + 250000) ; i++) 
	ranmar.Rndm();
      return 2 * Int_t(ranmar.Rndm() * (TMath::Power(2,30) - 1)) + 1;
    }
    break;
  }
  return 0;
}

void ktHijing::SetJetPtRange(Float_t minJetPt, Float_t maxJetPt)
{
  SetParameter("ihpr2",3,1); //IHPR2(3)=1
  SetParameter("hipr1",9,maxJetPt);
  SetParameter("hipr1",10,minJetPt);
}

void ktHijing::InitRHICdAu()
{
  Initialize(200, "CMS", "A", "A", 2, 1, 197, 79); //dAu @ RHIC
}

void ktHijing::InitRHICAuAu()
{
  Initialize(200, "CMS", "A", "A", 197, 79, 197, 79); //AuAu @ RHIC
}

void ktHijing::SetISRandFSR(int mode)
{
  //IHPR2(2):(D=3)
  //switch for initial and final state radiation in the hard scattering.
  //=0: both initial and final radiation off;
  //=1: initial radiation on and final radiation off;
  //=2: initial radiation off and final radiation on;
  //=3: both initial and final radiation on.

  SetParameter("ihpr2",2,mode);
}

void ktHijing::MakeMinbiasEvent()
{
  if (etaRange)
    {
      do
	{
	  GenerateEvent(0, 15, "final", particles);
	}
      while (!CheckJetEtaRange());
      nEvents++;
    }
  else
    {
      GenerateEvent(0, 15, "final", particles);
      nEvents++;
    }

  // eta,pt dist. of first parton = "jet"
  etaP->Fill(GetEtaFirstParton());
  ptP->Fill(GetPtFirstParton());
}

void ktHijing::FillFastJet(ktFastJet *fj,Int_t mSel)
{
  if (fj!=0)
    {
      // DEBUG:
      //cout<<GetNParticles()<<endl;

      for (Int_t l = 0; l < GetNParticles(); l++) 
	{
	  TParticle* part = (TParticle*) particles->At(l);	
	  
	  //if (part->GetStatusCode() !=1) continue;
	  // take all final state particles (implement cuts on charged and charged +  gamma ....
	  
	  if (Select_TT(mSel,part))
	    {
	      
	      //TLorentzVector *p=new TLorentzVector();
	      //p->SetPxPyPzE(part->Px(),part->Py(),part->Pz(),part->Energy());
	      
	      // check if particle is charged (for extracting the FF function in ktMuJet !!!)
	      Bool_t isCharged=true;
	      if (!(Select_TT(1,part))) isCharged=false;
	      
	      // check why, look up hijing ...
	      if (part->Pt()>0)
		fj->AddParticle(part->Px(),part->Py(),part->Pz(),part->Energy(),isCharged);
	      
	      //delete p;
	      
	    }
	}
      
    }
  else
    cout<<" Warning : Can not fill ktFastJet: zero pointer !!!"<<endl;
}

void ktHijing::FillGrid(ktGrid *mGrid, Int_t mSel)
{
   if (mGrid!=0)
    {
    }
   else
     cout<<" Warning : Can not fill ktGrid: zero pointer !!!"<<endl;
}

void ktHijing:: SetJetEtaRange(Float_t mEtaMin,Float_t mEtaMax)
{
  etaMin=mEtaMin;
  etaMax=mEtaMax;
  etaRange=true;
}

Bool_t ktHijing::CheckJetEtaRange()
{
  Float_t jEta=GetEtaFirstParton();
  
  if (jEta>etaMin && jEta<etaMax)
    return true;
  else
    return false;
}

Float_t ktHijing::GetEtaFirstParton()
{
  Float_t px=GetParameter("HINT1",21);
  Float_t py=GetParameter("HINT1",22);
  Float_t pz=GetParameter("HINT1",23);
  Float_t costheta=pz/(TMath::Sqrt(px*px+py*py+pz*pz));

  //Float_t eta=-TMath::Log(TMath::Tan(theta/2.0));
  Float_t eta=TMath::ATanH(costheta);
  return eta;
}

Float_t ktHijing::GetEtaFirstPartonFSR()
{
  Float_t px=GetParameter("HINT1",26);
  Float_t py=GetParameter("HINT1",27);
  Float_t pz=GetParameter("HINT1",28);
  Float_t costheta=pz/(TMath::Sqrt(px*px+py*py+pz*pz));

  //Float_t eta=-TMath::Log(TMath::Tan(theta/2.0));
  Float_t eta=TMath::ATanH(costheta);
  return eta;
}

Float_t ktHijing::GetPtFirstParton()
{
  Float_t px=GetParameter("HINT1",21);
  Float_t py=GetParameter("HINT1",22);
  
  return TMath::Sqrt(px*px+py*py);
}

Float_t ktHijing::GetPtFirstPartonFSR()
{
  Float_t px=GetParameter("HINT1",26);
  Float_t py=GetParameter("HINT1",27);
  
  return TMath::Sqrt(px*px+py*py);
}


Float_t ktHijing::GetPtSecondParton()
{
  Float_t px=GetParameter("HINT1",31);
  Float_t py=GetParameter("HINT1",32);
  
  return TMath::Sqrt(px*px+py*py);
}

Float_t ktHijing::GetPtSecondPartonFSR()
{
  Float_t px=GetParameter("HINT1",36);
  Float_t py=GetParameter("HINT1",37);
  
  return TMath::Sqrt(px*px+py*py);
}
