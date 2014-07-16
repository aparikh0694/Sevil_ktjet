// first test (Joern Putschke)

#include "ktMCBkg.h"
#include <Riostream.h>
#include "TRandom.h"

const Double_t PionMass=0.139;

// v2 param. made for high-pt
// fit not realistic at low-pt
// (necessary here for bkg.) !!!

/*
TF1* getv2pt(TString minormax)
{
  Int_t mx=8;
  Int_t mx2=9;

  TF1 *mv24=new TF1("mv24","(landau(0)+landau(3))/2.0",0.0,10.5);
  mv24->SetParameters((0.2618+0.3735*mx-0.0444*mx*mx),3.1,1.5,(0.2618+0.3735*mx2-0.0444*mx2*mx2),3.1,1.5);
  
  TF1 *mv2rp=new TF1("mv2rp","(landau(0)+landau(3))/2.0",0.0,10.5);
  mv2rp->SetParameters((0.759+0.278*mx-0.0367*mx*mx),(4.11-0.0934*mx),(2.623-0.15*mx),(0.759+0.278*mx2-0.0367*mx2*mx2),(4.11-0.0934*mx2),(2.623-0.15*mx2));

  TF1 *mval=new TF1("mval","(mv24+mv2rp)/2.0",0.0,10.5);

  TF1 *m0=new TF1("m0","0",0.0,10.5);

  if (minormax=="Max")
    return mv2rp;
  else if (minormax=="Min")
    return mv24;
  else if (minormax=="Mean")
    return mval;
  else
    return m0;
}

TF1* v2phi = new TF1("v2phi","1.0 + 2.0*[0]*[1]*TMath::Cos(2.0*x)",-1.0*TMath::Pi(),TMath::Pi());
 
  if (m_TrigPt<6)
  v2phi->SetParameter(0,v2ptnew->Eval(m_TrigPt));// *m_v2scale*prot_scale);
  else
    v2phi->SetParameter(0,v2ptnew->Eval(6.0));// *m_v2scale*prot_scale);
    
  v2phi->SetParameter(1,v2ptnew->Eval(m_AssPt));// *m_v2scale*prot_scale);
  
*/

ClassImp(ktMCBkg)

ktMCBkg::ktMCBkg()
{

  feta=new TF1();
  fphi=new TF1();
  fpt=new TF1();
  fRP=new TF1();
  fv2pt=new TF1();

  maxEta=1.0;

  heta=new TH1D("heta","eta distribution",20,-maxEta,maxEta);
  //hphi=new TH1D("hphi","phi distribution",180,0,2*TMath::Pi());
  //hRP=new TH1D("hRP","RP angle distribution",180,0,2*TMath::Pi());
  hphi=new TH1D("hphi","phi distribution",90,-TMath::Pi(),TMath::Pi());
  hRP=new TH1D("hRP","RP angle distribution",90,-TMath::Pi(),TMath::Pi());
  hpt=new TH1D("hpt","1/pt dN/pt distribution",50,0,5);

  mInitOK=false;
  output=true;

  v2=0.0;
  T=0.0;
  dNdEta=0.0;
  phiRP=0.0;

  mSel="RHIC";
  mSys="Au+Au";

  // DEBUG:
  //cout<<"Default ktMCBkg constructor "<<endl;

}

ktMCBkg::ktMCBkg(TString m_mSel,TString m_mSys, Double_t m_T, Double_t m_dNdEta, Double_t m_v2)
{
  // Dummy ... (to be implemented if necessary)
  
  maxEta=1.0;

  heta=new TH1D("heta","eta distribution",20,-maxEta,maxEta);
  //hphi=new TH1D("hphi","phi distribution",180,0,2*TMath::Pi());
  //hRP=new TH1D("hRP","RP angle distribution",180,0,2*TMath::Pi());
  hphi=new TH1D("hphi","phi distribution",90,-TMath::Pi(),TMath::Pi());
  hRP=new TH1D("hRP","RP angle distribution",90,-TMath::Pi(),TMath::Pi());
  hpt=new TH1D("hpt","1/pt dN/pt distribution",50,0,5);

  feta=new TF1();
  fphi=new TF1();
  fpt=new TF1();
  fv2pt=new TF1();

  mInitOK=false;
  output=true;

  mSel=m_mSel;
  mSys=m_mSys;
  v2=m_v2;
  T=m_T;
  dNdEta=m_dNdEta;
   phiRP=0.0;

  // DEBUG:
   // cout<<"full ktMCBkg constructor "<<endl;
} 

ktMCBkg::ktMCBkg(TString m_mSel,TString m_mSys, Double_t m_v2)
{

  maxEta=1.0;

  heta=new TH1D("heta","eta distribution",20,-maxEta,maxEta);
  //hphi=new TH1D("hphi","phi distribution",180,0,2*TMath::Pi());
  //hRP=new TH1D("hRP","RP angle distribution",180,0,2*TMath::Pi());
  hphi=new TH1D("hphi","phi distribution",90,-TMath::Pi(),TMath::Pi());
  hRP=new TH1D("hRP","RP angle distribution",90,-TMath::Pi(),TMath::Pi());
  hpt=new TH1D("hpt","1/pt dN/pt distribution",50,0,5);

  mInitOK=false;
  output=true;
   
  mSel=m_mSel;
  mSys=m_mSys;
  v2=m_v2;
  phiRP=0.0;

  if (mSel=="LHC")
    {
      dNdEta=1200.0;
      T=0.350;
    }
  else
    {
      // set default then ...
      mSel="RHIC";

      if (mSys=="Cu+Cu")
	{
	  dNdEta=200.0;
	  T=0.291;
	}
      else
	{
	  mSys="Au+Au";
	  dNdEta=650.0;
	  T=0.291;
	}
    }
  

  feta=new TF1("feta","1",-maxEta,maxEta);
  //fRP=new TF1("fRP","1",0,2*TMath::Pi());
  //fphi=new TF1("fphi","1+2*[0]*cos(2*(x-[1]))",0,2*TMath::Pi());
  fRP=new TF1("fRP","1",-TMath::Pi(),TMath::Pi());
  fphi=new TF1("fphi","1+2*[0]*cos(2*(x-[1]))",-TMath::Pi(),TMath::Pi());
  fphi->SetParameter(0,v2);
  fpt=new TF1("fpt","x*exp(-x/[0])",0.1,5);
  fpt->SetParameter(0,T);

  fv2pt=new TF1();
  // DEBUG:
  //cout<<"ktMCBkg(TString m_mSel,TString m_mSys, Double_t m_v2) constructor "<<endl;
}

ktMCBkg::~ktMCBkg()
{
  delete feta; delete fphi; delete fpt; delete fRP;
  delete heta; delete hphi; delete hpt; delete hRP;
  delete fv2pt;

  // DEBUG:
  //cout<<"Default ktMCBkg destructor"<<endl;

}

void ktMCBkg::PrintInfo()
{
  cout<<endl;
  cout<<" ---> MC background seetings:"<<endl;
  cout<<" ---> ------------------------"<<endl;
  cout<<" ---> (always top energies)"<<endl;
  cout<<" ---> Collider = "<<mSel<<endl;
  cout<<" ---> System   = "<<mSys<<endl;  
  cout<<" ---> pt spectrum (thermal)"<<endl;
  cout<<" ---> T        = "<<T<<endl;
  cout<<" ---> eta spectrum (flat) |eta|<1"<<endl;
  cout<<" ---> dN/deta  = "<<dNdEta<<" from charged * 3/2"<<endl;
  cout<<" ---> phi flat or v2 modulated"<<endl;
  cout<<" ---> v2       = "<<v2<<endl;
  cout<<" ---> Assume pion mass"<<endl;// for h+-"<<endl;
  cout<<endl;
}

Bool_t ktMCBkg::Init()
{
  Bool_t m_init=false;
  if (T!=0 && dNdEta!=0)
    {
      m_init=true;
      PrintInfo();
    }
  else
    {
      m_init=false;
      cout<<" ---> Error : ktMCBkg not proper initiliazed !"<<endl;
    }
  
  mInitOK=m_init;

  return m_init;
}

void ktMCBkg::SetRandomSeed()
{
  // Dummy ... (to be implemented if necessary)
}

void ktMCBkg::FillNeutral(ktGrid* myGrid)
{
  Int_t nMax=(Int_t) (2*maxEta*dNdEta*3/2.0);
  
  if (output)
    {
      cout<<" ---> Fill ktGrid with MC bkg. event ..."<<endl;
      cout<<" ---> Number of particles = "<<nMax<<endl;
      cout<<endl;
    }
  
  Double_t eta=0;
  Double_t phi=0;
  Double_t RP=0;
  Double_t pt=0;

  RP=fRP->GetRandom();
  hRP->Fill(RP);
  //fphi->SetParameter(0,v2); // Why SetParameter slows down running !????
  fphi->SetParameter(1,RP);

  for (int i=0;i<nMax;i++)
    {
      eta=feta->GetRandom();
      pt=fpt->GetRandom();
      //fphi->SetParameter(0,v2); // easy to implement v2(pt)
      phi=fphi->GetRandom();
      
      heta->Fill(eta);hphi->Fill(phi);hpt->Fill(pt,1/pt);

      if (myGrid!=0)
	{
	  TLorentzVector *pv=new TLorentzVector();
	  pv->SetPtEtaPhiM(pt,eta,phi,PionMass); 

	  if (gRandom->Rndm()>0.66)
	    myGrid->Fill((TLorentzVector*) pv->Clone(),false);
	  
	  delete pv;
	}
    }
}

void ktMCBkg::Fill(ktGrid* myGrid)
{
  Int_t nMax=(Int_t) (2*maxEta*dNdEta*3/2.0);
  
  if (output)
    {
      cout<<" ---> Fill ktGrid with MC bkg. event ..."<<endl;
      cout<<" ---> Number of particles = "<<nMax<<endl;
      cout<<endl;
    }
  
  Double_t eta=0;
  Double_t phi=0;
  Double_t RP=0;
  Double_t pt=0;

  RP=fRP->GetRandom();
  hRP->Fill(RP);
  //fphi->SetParameter(0,v2); // Why SetParameter slows down running !????
  fphi->SetParameter(1,RP);

  for (int i=0;i<nMax;i++)
    {
      eta=feta->GetRandom();
      pt=fpt->GetRandom();
      //fphi->SetParameter(0,v2); // easy to implement v2(pt)
      phi=fphi->GetRandom();
      
      heta->Fill(eta);hphi->Fill(phi);hpt->Fill(pt,1/pt);

      if (myGrid!=0)
	{
	  TLorentzVector *pv=new TLorentzVector();
	  pv->SetPtEtaPhiM(pt,eta,phi,PionMass); 

	  if (gRandom->Rndm()<0.66)
	    myGrid->Fill((TLorentzVector*) pv->Clone(),true);
	  else
	    myGrid->Fill((TLorentzVector*) pv->Clone(),false);
	  
	  delete pv;
	}
    }
}


void ktMCBkg::Fill(ktFastJet *myFastJet,Double_t m_RP)
{
  Int_t nMax=(Int_t) (2*maxEta*dNdEta*3/2.0);
  
  if (output)
    {
      cout<<" ---> Fill FastJet with MC bkg. event ..."<<endl;
      cout<<" ---> Number of particles = "<<nMax<<endl;
      cout<<endl;
    }
  
  Double_t eta=0;
  Double_t phi=0;
  Double_t pt=0;

  Double_t RP=m_RP;

  hRP->Fill(RP);
  //fphi->SetParameter(0,v2); // Why SetParameter slows down running !????
  fphi->SetParameter(1,RP);

  for (int i=0;i<nMax;i++)
    {
      eta=feta->GetRandom();
      pt=fpt->GetRandom();

      //fphi->SetParameter(0,v2); 
      // easy to implement v2(pt)
      // but why does it slows down !?????
      
      phi=fphi->GetRandom();
      
      heta->Fill(eta);hphi->Fill(phi);hpt->Fill(pt,1/pt);

      if (myFastJet!=0)
	{
	  TLorentzVector *pv=new TLorentzVector();
	  pv->SetPtEtaPhiM(pt,eta,phi,PionMass); 

	  if (gRandom->Rndm()<0.66)	
	      myFastJet->AddParticle(pv->Px(),pv->Py(),pv->Pz(),pv->E(),true);
	  else
	      myFastJet->AddParticle(pv->Px(),pv->Py(),pv->Pz(),pv->E(),false);

	  delete pv;
	}
    }
}

void ktMCBkg::Fill(ktGrid* myGrid,ktFastJet *myFastJet,Double_t m_RP)
{
  Int_t nMax=(Int_t) (2*maxEta*dNdEta*3/2.0);
  
  if (output)
    {
      cout<<" ---> Fill ktGrid & FastJet with MC bkg. event ..."<<endl;
      cout<<" ---> Number of particles = "<<nMax<<endl;
      cout<<endl;
    }
  
  Double_t eta=0;
  Double_t phi=0;
  Double_t pt=0;

  Double_t RP=m_RP;

  hRP->Fill(RP);
  //fphi->SetParameter(0,v2); // Why SetParameter slows down running !????
  fphi->SetParameter(1,RP);

  for (int i=0;i<nMax;i++)
    {
      eta=feta->GetRandom();
      pt=fpt->GetRandom();

      //fphi->SetParameter(0,v2); 
      // easy to implement v2(pt)
      // but why does it slows down !?????
      
      phi=fphi->GetRandom();
      
      heta->Fill(eta);hphi->Fill(phi);hpt->Fill(pt,1/pt);

      if (myFastJet!=0)
	{
	  TLorentzVector *pv=new TLorentzVector();
	  pv->SetPtEtaPhiM(pt,eta,phi,PionMass); 

	  if (gRandom->Rndm()<0.66)
	    {
	      myGrid->Fill((TLorentzVector*) pv->Clone(),true);
	      myFastJet->AddParticle(pv->Px(),pv->Py(),pv->Pz(),pv->E(),true);
	    }
	  else
	    {
	      myGrid->Fill((TLorentzVector*) pv->Clone(),false);
	      myFastJet->AddParticle(pv->Px(),pv->Py(),pv->Pz(),pv->E(),false);
	    }

	  delete pv;
	}
    }
}


void ktMCBkg::Fill(ktGrid* myGrid, Double_t m_RP)
{
  Int_t nMax=(Int_t) (2*maxEta*dNdEta*3/2.0);
  
  if (output)
    {
      cout<<" ---> Fill ktGrid with MC bkg. event ..."<<endl;
      cout<<" ---> Number of particles = "<<nMax<<endl;
      cout<<endl;
    }
  
  Double_t eta=0;
  Double_t phi=0;
  Double_t pt=0;

  Double_t RP=m_RP;

  hRP->Fill(RP);
  //fphi->SetParameter(0,v2); // Why SetParameter slows down running !????
  fphi->SetParameter(1,RP);

  for (int i=0;i<nMax;i++)
    {
      eta=feta->GetRandom();
      pt=fpt->GetRandom();

      //fphi->SetParameter(0,v2); 
      // easy to implement v2(pt)
      // but why does it slows down !?????
      
      phi=fphi->GetRandom();
      
      heta->Fill(eta);hphi->Fill(phi);hpt->Fill(pt,1/pt);

      if (myGrid!=0)
	{
	  TLorentzVector *pv=new TLorentzVector();
	  pv->SetPtEtaPhiM(pt,eta,phi,PionMass); 

	  if (gRandom->Rndm()<0.66)
	    myGrid->Fill((TLorentzVector*) pv->Clone(),true);
	  else
	    myGrid->Fill((TLorentzVector*) pv->Clone(),false);
	  
	  delete pv;
	}
    }
}

void ktMCBkg::Fill(ktGrid* myGrid,TString mSel)
{
  Int_t nMax=(Int_t) (2*maxEta*dNdEta*3/2.0);
  
  if (output)
    {	
      cout<<" ---> Fill ktGrid with MC bkg. event ..."<<endl;
      cout<<" ---> Particle selection = "<<mSel<<endl;
      cout<<" ---> Number of particles = "<<nMax<<endl;
      cout<<endl;
      cout<<" Warning: mSel not implemented yet, run as ktMCBkg::Fill()"<<endl;
      cout<<endl;
    }

  Fill(myGrid);

}

void ktMCBkg::HistNorm(Int_t nEvents)
{
  Double_t norm=(Double_t) nEvents;

  heta->Scale(1/norm*1/heta->GetBinWidth(1));
  hphi->Scale(1/norm*1/hphi->GetBinWidth(1));
  hpt->Scale(1/norm*1/hpt->GetBinWidth(1));
  //hRP->Scale(1/norm*1/hRP->GetBinWidth(1));

  heta->SetMinimum(0.0);
  hphi->SetMinimum(0.0);
  hpt->SetMinimum(0.0);
  hRP->SetMinimum(0.0);
}
