// first test (Joern Putschke)

#include "ktPythia8.h"
#include "TRandom.h"
#include<iostream> // needed for io
#include<sstream>  // needed for internal io
#include<vector> 

ClassImp(ktPythia8)

ktPythia8::ktPythia8() 
{
 
  p=new Pythia("/home/putschke/pythia8100/xmldoc");
  pyeta=pyphi=pyphi2=pypt=pye=mySqrts=0.0;
  verbose=true;
  smear=false;

  // hardcoded for STAR resolution (fix later !)
  ptRes=new TF1("ptRes","0.0051+0.0111*x",0,50);
  eRes=new TF1("eRes","0.15/sqrt(x)+0.03",0,50);

  //DEBUG:
  //cout<<endl;
  //cout<<"Default constructor of Pythia8 ..."<<endl;
}

ktPythia8::ktPythia8(char* xmlDir) 
{
 
  p=new Pythia(xmlDir);
  pyeta=pyphi=pyphi2=pypt=pye=mySqrts=0.0;
  verbose=true;
  smear=false;

  // hardcoded for STAR resolution (fix later !)
  ptRes=new TF1("ptRes","0.0051+0.0111*x",0,50);
  eRes=new TF1("eRes","0.15/sqrt(x)+0.03",0,50);

  //DEBUG:
  //cout<<endl;
  //cout<<"Constructor of Pythia8 ..."<<endl;
}

ktPythia8::~ktPythia8()
{
 
  delete p;
  delete ptRes; delete eRes;

  //DEBUG:
  //cout<<endl;
  //cout<<"Destructor of Pythia8 ..."<<endl;
}

void ktPythia8::SetRandomSeed()
{
  // Set Random Seed time based
 
  cout<<" ---> Random seed set (time based)"<<endl;
  cout<<endl;
  p->readString("Random:setSeed = on");
  p->readString("Random:seed = 0");
}

void ktPythia8::InitPhaseSpace(char* ptHatMin,char* ptHatMax)
{
  string pmin="PhaseSpace:pTHatMin = ";
  string pmax="PhaseSpace:pTHatMax = ";
  pmin.append(ptHatMin);
  pmax.append(ptHatMax);

  cout<<" ---> Set Pythai8 phasespce:"<<endl;
  cout<<" ---> ----------------------"<<endl;
  cout<<" ---> "<<pmin<<endl;
  cout<<" ---> "<<pmax<<endl;

  p->readString(pmin);
  p->readString(pmax);

}

void ktPythia8::Init(TString mInit)
{
  cout<<endl;
  p->readString("HardQCD:all = on");
  p->readString("HadronLevel:Decay = on");
  //p->readString("PartonLevel:ISR = off");
  //p->readString("PartonLevel:FSR = off");

  // DEBUG:
  //p->readString("Main:showAllSettings = on");

  // check with ALICE simulation, how to deal with the decays
  // switch off part of it only ... !!!????

  cout<<" ---> Pythia8 settings:"<<endl;
  cout<<" ---> -----------------"<<endl;
  cout<<" ---> HardQCD:all = on"<<endl;
  cout<<" ---> HadronLevel:Decay = on"<<endl;
  cout<<" ---> PartonLevel:ISR = on"<<endl;
  cout<<" ---> PartonLevel:FSR = on"<<endl;

  cout<<endl;

  if (mInit=="LHC")
    {
      cout<<" ---> Init Pythia8 for LHC ..."<<endl;
      //      mySqrts= 14000.;
      mySqrts=2760.;
      p->init( 2212, 2212, 2760.);
    }
  else
    {
      cout<<"Init Pythia8 for RHIC ..."<<endl;
      mySqrts= 200.;
      p->init( 2212, 2212, 200.);
    }
  
  cout<<endl;
}

void ktPythia8::Run()
{
  if (!p->next())
    {
      cout<<"--> Error in Pythia8 event ..."<<endl;
    }
}

void ktPythia8::RunWithPycell()
{
  if (!p->next())
    {
      cout<<"--> Error in Pythia8 event ..."<<endl;
    }
  else
    {
       CellJet cellJet;
       cellJet.analyze(p->event,5,0.7,1.5);
       if (cellJet.size()>0)
	 {
	   Vec4 mvec=cellJet.pMassive(0);
	   Double_t mPhi=cellJet.phiWeighted(0);if (mPhi<0) mPhi += (2*TMath::Pi());
	   pyphi=mPhi;
	   pyphi2=cellJet.phiWeighted(0);
	   pyeta=cellJet.etaWeighted(0);
	   pypt=mvec.pT();
	   pye=cellJet.eT(0);
	 }
       else
	 {
	   cout<<" Warning: PYTHIA8 Pycell did not found jet !!!"<<endl;
	   pyphi=pyphi2=pyeta=pypt=pye=0.0;
	 }
    }
}

void ktPythia8::Pycell()
{
  CellJet cellJet;
  cellJet.analyze(p->event );
  cout<<" Pycell output "<<(Int_t) cellJet.size()<<" Jets:"<<endl;
  printf(" %5s %15s %15s %15s %15s %12s %12s\n","jet #", "eta",
         "phi", "E", "pt","n particles","Type");

  for (int i = 0; i < cellJet.size(); ++i) 
    {
      Vec4 mvec=cellJet.pMassive(i);
      Double_t mPhi=cellJet.phiWeighted(i);if (mPhi<0) mPhi += (2*TMath::Pi());
      printf(" %5u %15.8f %15.8f %15.8f %15.8f %12u %12s\n",i,cellJet.etaWeighted(i),mPhi,cellJet.eT(i),mvec.pT(),(Int_t) cellJet.multiplicity(i),"Pycell");
    }
}

void ktPythia8::HardProcessInfo()
{
  // check values with jet-finder internal and ktJet !???
  // Probably have to transform in CM frame ???
  // include cellJet analysis from PYTHIA ...

  cout<<" Hard process info from Pythia8:"<<endl;
  cout<<" -------------------------------"<<endl;
  
  Double_t hp_eta=-TMath::Log(TMath::Tan(p->info.thetaHat()/2.0));
  cout<<" eta = "<<hp_eta<<", phi = "<<p->info.phiHat()<<", theta = "<<p->info.thetaHat()<<", ptHat = "<<p->info.pTHat()<<endl;
  cout<<endl;
}

void ktPythia8::Fill(ktGrid* myGrid)
{
  if (verbose)
    {
      cout<<endl;
      cout<<" ---> Fill ktGrid with Pythia8 event ..."<<endl;
      cout<<" ---> Number of final state particles = "<< p->event.size()<<endl;
      cout<<endl;
    }

  for (int i = 0; i < p->event.size(); ++i) 
    {
      if (p->event[i].isFinal())
	{

	  Double_t px=p->event[i].px();
	  Double_t py=p->event[i].py();
	  Double_t pz=p->event[i].pz();
	  Double_t E=p->event[i].e();

	  ktParticle *pv=new ktParticle();
	  pv->SetPxPyPzE(px,py,pz,E);
	  
	  ktPID* pid= new ktPID();
	  pid->SetPID(p->event[i].id());
	  pid->SetCharge((Int_t)p->event[i].charge());
	  if( pid->GetCharge() == 0) pid->SetETower(p->event[i].e());

	  
	  myGrid->Fill((TLorentzVector*) pv->Clone(),pid);
	  
	  delete pid;
	  delete pv; 
	}
    }
  
}

// for backward compapility ...
void ktPythia8::Fill(ktGrid* myGrid,TString mSel, Bool_t mSmear)
{
  smear=mSmear;
  Fill(myGrid,mSel);
}

void ktPythia8::Fill(ktGrid* myGrid,TString mSel)
{
  // Fill grid with Pythia8
  // Selections : TPC (charged only), EMCAL (charged+neutral) and ALL (all final state)

  if (verbose)
    {
      cout<<endl;
      cout<<" ---> Fill ktGrid with Pythia8 event ..."<<endl;
      cout<<" ---> Number of particles = "<< p->event.size()<<endl;
      //cout<<" ---> Number of final state particles = "<< p->info.nFinal()<<endl;
      cout<<" ---> Particle selection = "<<mSel<<endl;
    }
  
  Int_t NSel=0;

  for (int i = 0; i < p->event.size(); ++i) 
    {
      if (p->event[i].isFinal())
	{

	  Double_t px=p->event[i].px();
	  Double_t py=p->event[i].py();
	  Double_t pz=p->event[i].pz();
	  Double_t E=p->event[i].e();

	  TLorentzVector *pv=new TLorentzVector();
	  pv->SetPxPyPzE(px,py,pz,E);

	  ktPID* pid = new ktPID();
	  pid->SetCharge((Int_t)p->event[i].charge());
	  pid->SetPID(p->event[i].id());
	  if( pid->GetCharge() ==0) pid->SetETower(p->event[i].e());
	  
	  // DEBUG:
	  //if(p->event[i].id()==111) cout<<"pi0"<<endl;
	  // check if decay on ...

	  if (mSel=="TPC")
	    {
	      if (p->event[i].isCharged())
		{
		  if (!smear)
		    myGrid->Fill((TLorentzVector*) pv->Clone(),pid);
		  else
		    {
		      Double_t mdpt=gRandom->Gaus(0,ptRes->Eval(pv->Pt()));
		      Double_t mpts=pv->Pt()+mdpt*pv->Pt();
		      
		      pv->SetPtEtaPhiE(mpts,pv->Eta(),pv->Phi(),pv->E());
		      myGrid->Fill((TLorentzVector*) pv->Clone(),pid); 
		    }
		  NSel++;
		}
	    }
	  else if (mSel=="NEUTRAL")
	    {
	      if (p->event[i].id()==22 || p->event[i].id()==111)
		{
		  if (!smear)
		    myGrid->Fill((TLorentzVector*) pv->Clone(),pid);
		  else
		    {
		      // check if E or Et ...
		      Double_t mde=gRandom->Gaus(0,eRes->Eval(pv->Pt()));
		      Double_t mes=pv->Pt()+mde*pv->Pt();
		      
		      pv->SetPtEtaPhiE(mes,pv->Eta(),pv->Phi(),pv->E());
		      myGrid->Fill((TLorentzVector*) pv->Clone(),pid); 
		    }
		  NSel++;
		}
	    }
	  else if (mSel=="EMCAL")
	    {
	      //DEBUG:
	      //cout<<"in emcal ..."<<endl;

	      if (p->event[i].isCharged())
		{
		  if (!smear)
		    myGrid->Fill((TLorentzVector*) pv->Clone(),pid);
		  else
		    {
		      Double_t mdpt=gRandom->Gaus(0,ptRes->Eval(pv->Pt()));
		      Double_t mpts=pv->Pt()+mdpt*pv->Pt();
		      
		      pv->SetPtEtaPhiE(mpts,pv->Eta(),pv->Phi(),pv->E());
		      myGrid->Fill((TLorentzVector*) pv->Clone(),pid); 
		    }
		  NSel++;
		}
	      else if (p->event[i].id()==22 || p->event[i].id()==111)
		{
		  if (!smear)
		    myGrid->Fill((TLorentzVector*) pv->Clone(),pid);
		  else
		    {
		      // check if E or Et ...
		      Double_t mde=gRandom->Gaus(0,eRes->Eval(pv->Pt()));
		      Double_t mes=pv->Pt()+mde*pv->Pt();
		      
		      pv->SetPtEtaPhiE(mes,pv->Eta(),pv->Phi(),pv->E());
		      myGrid->Fill((TLorentzVector*) pv->Clone(),pid); 
		    }
		  NSel++;
		}
	    }
	  else
	    {
	       if (p->event[i].isCharged())
		 {
		   myGrid->Fill((TLorentzVector*) pv->Clone(),pid);
		 }
	       else
		 {
		   myGrid->Fill((TLorentzVector*) pv->Clone(),pid);
		 }
	       NSel++;
	    }
	  
	  delete pv; 
	  delete pid;
	}
    }

  if (verbose)
    cout<<" ---> Number of particles filled = "<<NSel<<endl;
  //cout<<endl;
}

void ktPythia8::GetHighestJetForQuenching(ktJetQuench *mQ, Double_t Rc,TString mSel, Bool_t mSmear)
{

  smear=mSmear;

  if (verbose)
    cout<<"---> Get highest jet (pythia) for quenching ..."<<endl;

  // DEBUG:
  //cout<<GetPycellPt()<<" "<<GetPycellEta()<<" "<<GetPycellPhi()<<" "<< GetPycellPhi2()<<endl;
  
  mQ->SetJetPt(GetPycellPt());
  mQ->SetJetEta(GetPycellEta());
  mQ->SetJetPhi(GetPycellPhi2());

  for (int i = 0; i < p->event.size(); ++i) 
    {
      if (p->event[i].isFinal())
	{

	  Double_t px=p->event[i].px();
	  Double_t py=p->event[i].py();
	  Double_t pz=p->event[i].pz();
	  Double_t E=p->event[i].e();

	  TLorentzVector *pv=new TLorentzVector();
	  pv->SetPxPyPzE(px,py,pz,E);
	  
	  Double_t dPhi=pv->Phi()-GetPycellPhi2();
	  Double_t dEta=pv->Eta()-GetPycellEta();
	  Double_t dR=TMath::Sqrt(dPhi*dPhi+dEta*dEta);

	  if (dR<Rc)
	    {
	      // DEBUG:
	      //cout<<pv->Phi()<" "<<pv->Eta()<<endl;

	      if (mSel=="EMCAL")
		{
		  //DEBUG:
		  //cout<<"in emcal ..."<<endl;
		  
		  if (p->event[i].isCharged())
		    {
		      if (!smear)
			mQ->AddParticle(px,py,pz,E,1);				      
		      else
			{
			  Double_t mdpt=gRandom->Gaus(0,ptRes->Eval(pv->Pt()));
			  Double_t mpts=pv->Pt()+mdpt*pv->Pt();			  
			  
			  pv->SetPtEtaPhiE(mpts,pv->Eta(),pv->Phi(),pv->E());
			  mQ->AddParticle(pv->Px(),pv->Py(),pv->Pz(),pv->E(),1);
			}
		    }
		  else if (p->event[i].id()==22 || p->event[i].id()==111)
		    {
		      if (!smear)
			mQ->AddParticle(px,py,pz,E,0);			     
		      else
			{
			  // check if E or Et ...
			  Double_t mde=gRandom->Gaus(0,eRes->Eval(pv->Pt()));
			  Double_t mes=pv->Pt()+mde*pv->Pt();
			  
			  pv->SetPtEtaPhiE(mes,pv->Eta(),pv->Phi(),pv->E());
			  mQ->AddParticle(pv->Px(),pv->Py(),pv->Pz(),pv->E(),0);
			}
		    }
		}
	      else
		{
		  if (p->event[i].isCharged())		    
		    mQ->AddParticle(px,py,pz,E,1);		  
		  else
		    mQ->AddParticle(px,py,pz,E,0);
		  
		}
	    }
	  
	  delete pv;
	}
    }
}

void ktPythia8::FillEvent(ktPy8Event *myEv)
{
  // Fill generated event in ktPy8Event ...

  myEv->SetPtHat(p->info.pTHat());
  myEv->SetPhiHat(p->info.phiHat());
  myEv->SetThetaHat(p->info.thetaHat());
  myEv->SetPyCellEta(GetPycellEta());
  myEv->SetPyCellPhi(GetPycellPhi());
  myEv->SetPyCellPt(GetPycellPt());
  myEv->SetSqrtS(mySqrts);
  myEv->SetWeight(p->info.sigmaGen());//1.0); // check how to get it e-by-e !???

  for (int i = 0; i < p->event.size(); ++i) 
    {
      if (p->event[i].isFinal())
	{
	  Double_t px=p->event[i].px();
	  Double_t py=p->event[i].py();
	  Double_t pz=p->event[i].pz();
	  Double_t E=p->event[i].e();

	  ktParticle *pv=new ktParticle();
	  pv->SetPxPyPzE(px,py,pz,E);
	  ktPID *pid= new ktPID();
	  pid->SetPID(p->event[i].id());
	  pid->SetCharge((Int_t)p->event[i].charge());
	  if( pid->GetCharge() ==0) pid->SetETower(p->event[i].e());
	  
	  if (p->event[i].isCharged())
	    myEv->Add(pv,"charged");
	  else if (p->event[i].id()==22 || p->event[i].id()==111)
	    myEv->Add(pv,"neutral");
	  else
	    myEv->Add(pv,"all");

	  delete pv;
	  delete pid;
	}
    }
}
