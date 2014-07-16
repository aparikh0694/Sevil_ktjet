// first test (Joern Putschke)

#include "ktPy8Event.h"
#include<Riostream.h> // needed for io

ClassImp(ktPy8Event)

ktPy8Event::ktPy8Event() 
{

  neutralList=new TObjArray(0);
  chargedList=new TObjArray(0);
  allList=new TObjArray(0);
  
  evNum=-1;
  mSqrts=0;
  mWeight=1.0;

  //nCharged=nNeutral=nAll=0;
  mptHat=mthetaHat=mphiHat=pycellPt=pycellEta=pycellPhi=0.0;
  
  //DEBUG:
  //cout<<endl;
  //cout<<"Default constructor of Py8Event ..."<<endl;
}

ktPy8Event::ktPy8Event(Int_t myEvNum) 
{

  neutralList=new TObjArray(0);
  chargedList=new TObjArray(0);
  allList=new TObjArray(0);
  
  evNum=myEvNum;
  mSqrts=0;
  mWeight=1.0;

  //nCharged=nNeutral=nAll=0;
  mptHat=mthetaHat=mphiHat=pycellPt=pycellEta=pycellPhi=0.0;
  
  //DEBUG:
  //cout<<endl;
  //cout<<"Default constructor of Py8Event ..."<<endl;
}

ktPy8Event::~ktPy8Event() 
{

  chargedList->Delete();
  delete chargedList;

  neutralList->Delete();
  delete neutralList;

  allList->Delete();
  delete allList;

  //DEBUG:
  //cout<<endl;
  //cout<<"Default destructor of Py8Event ..."<<endl;
}

void ktPy8Event::EventInfo()
{
  cout<<"Event #"<<evNum<<endl;
  cout<<"# particles = "<<GetN()<<endl;
  cout<<"Sqrt(sNN) = "<<mSqrts<<endl;
  cout<<"ptHat = "<<mptHat<<endl;
  cout<<"thetaHat = "<<mthetaHat<<" | phiHat = "<<mphiHat<<endl;
  cout<<"weight = "<<mWeight<<endl;
}

void ktPy8Event::Add(ktParticle *v,TString mType)
{
  // Add TLorentzVector to Pythia8 event structure
  // Selections : charged (charged only), neutral (charged+neutral) and all (rest)
  if (mType=="charged")
    {
      chargedList->AddLast((ktParticle*) v->Clone());
    }
  else if (mType=="neutral")
    {
      neutralList->AddLast((ktParticle*) v->Clone());
    }	
  else
    {
      allList->AddLast((ktParticle*) v->Clone());
    }
  
}

void ktPy8Event::Fill(ktGrid* myGrid,TString mSel)
{
  // Fill grid with Pythia8 event structure
  // Selections : TPC (charged only), EMCAL (charged+neutral) and ALL (all final state)
  
  if (mSel=="TPC")
    {
      for (int i=0;i<chargedList->GetEntries();i++)
	{
	  myGrid->Fill((TLorentzVector*) GetCharged(i)->Clone(),true);
	}	
    }
  else if (mSel=="EMCAL")
    {
      for (int i=0;i<chargedList->GetEntries();i++)
	{
	  myGrid->Fill((TLorentzVector*) GetCharged(i)->Clone(),true);
	}
      for (int j=0;j<neutralList->GetEntries();j++)
	{
	  myGrid->Fill((TLorentzVector*) GetNeutral(j)->Clone(),false);
	}
    }
  else
    {
      for (int i=0;i<chargedList->GetEntries();i++)
	{
	  myGrid->Fill((TLorentzVector*) GetCharged(i)->Clone(),true);
	}
      for (int j=0;j<neutralList->GetEntries();j++)
	{
	  myGrid->Fill((TLorentzVector*) GetNeutral(j)->Clone(),false);
	}
      for (int k=0;k<allList->GetEntries();k++)
	{
	   myGrid->Fill((TLorentzVector*) GetAll(k)->Clone(),false);
	}
    }

}
