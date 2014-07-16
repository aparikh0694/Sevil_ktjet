
#include "ktFastJet.h"
#include "ktAreaUtil.h"
#include "ktTrackSimUtil.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceActiveArea.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/SISConePlugin.hh" // include SISCone plugin
#include "fastjet/Selector.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include <fastjet/tools/BackgroundEstimatorBase.hh>
#include "fastjet/tools/Subtractor.hh" 
#include "TRandom.h"

#include<iostream> // needed for io
#include<sstream>  // needed for internal io
#include<vector> 

using namespace std;
using namespace fastjet;

vector<fastjet::PseudoJet> input_particles;
vector<fastjet::PseudoJet> input_particles_charged;
vector<fastjet::PseudoJet> input_particles_neutral;

// extract the inclusive jets with pt > 2 GeV, sorted by pt
const double ptmin = 2.0;

//----------------------------------------------------------------------
/// a function that pretty prints a list of jets
void print_jets_sub(const fastjet::ClusterSequenceAreaBase & clust_seq, 
		    const vector<fastjet::PseudoJet> & unsorted_jets , double median_pt_per_area) {

  // sort jets into increasing pt
 vector<fastjet::PseudoJet> jets = sorted_by_pt(unsorted_jets);  
  // the corrected jets will go in here
  vector<fastjet::PseudoJet> corrected_jets(jets.size());
  
  // get median pt per unit area
  // NB pt_per_unit_area exists only in ClusterSequenceActiveArea, and
  // not in the base class ClusterSequenceWithArea.
  // This might change in a future release
  //double median_pt_per_area = clust_seq.pt_per_unit_area();

  cout << "Printing inclusive jets with pt > "<< ptmin<<" GeV\n";
  cout << "---------------------------------------\n";
  printf(" ijet     rap     phi        Pt    area  +- area   Pt corr  (rap corr phi corr Pt corr)ext\n");

  for (size_t j = 0; j < jets.size(); j++) {

    // get area of each jet
    double area     = clust_seq.area(jets[j]);
    double areaErr     = clust_seq.area_error(jets[j]);

    // "standard" correction. Subtract only the Pt
    double pt_corr  = jets[j].perp() - area*median_pt_per_area;

    // "extended" correction
    fastjet::PseudoJet area_4vect = 
                       median_pt_per_area*clust_seq.area_4vector(jets[j]);
    if (area_4vect.perp2() >= jets[j].perp2() || 
	area_4vect.E()     >= jets[j].E()) {
      // if the correction is too large, set the jet to zero
      corrected_jets[j] =  0.0 * jets[j];
    } else {
      // otherwise do an E-scheme subtraction
      corrected_jets[j] = jets[j] - area_4vect;
    }
    // NB We could also leave out the above "if": 
    // corrected_jets[j] = jets[j] - area_4vect;
    // but the result would be different, since we would not avoid
    // jets with negative Pt or energy
    
    if (pt_corr>ptmin)

    printf("%5u %7.3f %7.3f %9.3f %7.3f %9.3f %9.3f %7.3f %7.3f %9.3f\n",
	   (int) j,jets[j].rap(),jets[j].phi(),jets[j].perp(), area, areaErr, pt_corr,
     corrected_jets[j].rap(),corrected_jets[j].phi(), corrected_jets[j].perp());
  }

  cout << endl;
  cout << "median pt_over_area = " << median_pt_per_area << endl << endl;


}


void fill_jets_in_event(const fastjet::ClusterSequenceAreaBase & clust_seq, 
			const vector<fastjet::PseudoJet> & unsorted_jets , double median_pt_per_area, ktMuEvent *mev,TString myType, double median_charged, double median_neutral) //, double etaAna=1.0, double ptAna=2.0)
{

   // sort jets into increasing pt
  vector<fastjet::PseudoJet> jets = sorted_by_pt(unsorted_jets);  
  // the corrected jets will go in here
  vector<fastjet::PseudoJet> corrected_jets(jets.size());
  TLorentzVector *nN = new TLorentzVector(0,0,0,0);
  TLorentzVector *nC = new TLorentzVector(0,0,0,0);

  for (size_t j = 0; j < jets.size(); j++) {
    double area     = clust_seq.area(jets[j]);
    double areaE    = clust_seq.area_error(jets[j]);
    // "standard" correction. Subtract only the Pt
    //double pt_corr  = jets[j].perp() - area*median_pt_per_area;

    // "extended" correction
    fastjet::PseudoJet area_4vect = 
                       median_pt_per_area*clust_seq.area_4vector(jets[j]);
    if (area_4vect.perp2() >= jets[j].perp2() || 
	area_4vect.E()     >= jets[j].E()) {
      // if the correction is too large, set the jet to zero
      corrected_jets[j] =  0.0 * jets[j];
    } else {
      // otherwise do an E-scheme subtraction
      corrected_jets[j] = jets[j] - area_4vect;
    }

    //cout<<"DEBUG: "<<median_pt_per_area<<endl;

    // some variables ...
    double mNEF=-99;
    double mNEFCorr=-99;
    
    nC->SetXYZM(0,0,0,0);
    nN->SetXYZM(0,0,0,0);
    double jPt=jets[j].perp();
    double jPtCorr=jets[j].perp()-area*median_pt_per_area;

    // get jet constituents ...
    // maybe add only for the highest jets !???
    vector<fastjet::PseudoJet> constituents = clust_seq.constituents(jets[j]);
    int nCon= constituents.size();
    //DEBUG:    
    //cout<<nCon<<endl;

    for (int i=0; i < nCon; i++)
      {
	fastjet::PseudoJet mPart=constituents[i];

	if (mPart.user_index()>0){
	  nC->SetPxPyPzE(nC->Px()+mPart.px(), 
			 nC->Py()+mPart.py(),
			 nC->Pz()+mPart.pz(), nC->E()+mPart.E());
	  // cout << "Charged: " << mPart.px() << " " << mPart.py() << " " << mPart.pz() << " " << mPart.E() << " Total: "  << nC->E() << endl;
	}
	else{
	  nN->SetPxPyPzE(nN->Px()+mPart.px(),
			 nN->Py()+mPart.py(),
			 nN->Pz()+mPart.pz(), nN->E()+mPart.E());
	  // cout << "Neutral: " << mPart.px() << " " << mPart.py() << " " << mPart.pz() << " " << mPart.E() << " Total: "  << nN->E() << endl;
	}
      }
      
    if (jPt>0){
      mNEF=nN->E()/(nN->E()+nC->E());
      // cout << "NEF " << mNEF << endl;
    }
    if (jPtCorr>0)
      //  mNEFCorr=(nN-area*median_neutral)/jPtCorr;
      mNEFCorr = (nN->E()-area*median_neutral)/((nN->E()-area*median_neutral)+(nC->Pt()-area*median_charged));

    //if (TMath::Abs(jets[j].rap())<etaAna && jPtCorr>ptAna)
    mev->AddFastJet(jets[j].perp(),jets[j].rap(),jets[j].phi(),corrected_jets[j].rap(),corrected_jets[j].phi(),median_pt_per_area,area,myType,mNEF,mNEFCorr,areaE);

  }
}




//----------------------------------------------------------------------

ClassImp(ktFastJet)

ktFastJet::ktFastJet()
{

  mEtaMax=10.0; 
  mPhiMax=TMath::Pi();
  mEtaAna=1.0;
  mPtAna=0.0;

  mPtCut=0;
  median_pt_per_area=0;
  median_pt_per_area_charged=0;
  median_pt_per_area_neutral=0;

  info=true;
  printjet=true;
  spaceCharge=false;
  smear=false;

  bkgCalculated=false; //Modify!

  // default (active) area definition
  ghost_etamax = 1.0; //default 6.0
  active_area_repeats = 1; //default 3
  ghost_area    = 0.01; // default 0.01

  FFArray=new TObjArray(0);
  FFArray->SetOwner(kTRUE);
  
  PhiBkgHistos=new TObjArray(0);
  PhiBkgHistos->SetOwner(kTRUE);
 
  median_ptArray=new TArray();
  // median_ptArray->SetOwner(kTRUE);

  fillFF=false;
  histFFonly=false;
  
  
  RcFF=0.7;
  nFF=2;
  
  // Add function to set bin an xi range as variables ....
  //hXiBkgFJ=new TH1D("hXiBkgFJ","Xi dist. of Background per event FastJet",50,0,10);
  //hXiBkgFJ->SetDirectory(0); // to avoid root from owning 

  //hZBkgFJ=new TH1D("hZBkgFJ","Z dist. of Background per event FastJet",10,0,1);
  //hZBkgFJ->SetDirectory(0); // to avoid root from owning

  //hptBkgFJ=new TH1D("hptBkgFJ","pt dist. of Background per event FastJet",200,0,10);
  //hptBkgFJ->SetDirectory(0); // to avoid root from owning

}

ktFastJet::ktFastJet(Int_t nXiBins,Double_t xiMax)
{

  mEtaMax=10.0;
  mPhiMax=TMath::Pi();
  mEtaAna=1.0;
  mPtAna=0.0;

  mPtCut=0;
  median_pt_per_area=0;
  median_pt_per_area_charged=0;
  median_pt_per_area_neutral=0;
  info=true;
  printjet=true;

  bkgCalculated=false; //Modify!

  // default (active) area definition
  ghost_etamax = 1.0; //default 6.0
  active_area_repeats = 1; //default 3
  ghost_area    = 0.01; // default 0.01

  FFArray=new TObjArray(0);
  FFArray->SetOwner(kTRUE);

  PhiBkgHistos=new TObjArray(0);
  PhiBkgHistos->SetOwner(kTRUE);

  //median_ptArray=new TArray();
  //median_ptArray->SetOwner(kTRUE);

  fillFF=false;
  histFFonly=false;

  RcFF=0.7;
  nFF=2;
  
  // Add function to set bin an xi range as variables ....
  //hXiBkgFJ=new TH1D("hXiBkgFJ","Xi dist. of Background per event FastJet",nXiBins,0,xiMax);
  //hXiBkgFJ->SetDirectory(0); // to avoid root from owning 

  //hZBkgFJ=new TH1D("hZBkgFJ","Z dist. of Background per event FastJet",(Int_t) 10,0,1);
  //hZBkgFJ->SetDirectory(0); // to avoid root from owning 

  //hptBkgFJ=new TH1D("hptBkgFJ","pt dist. of Background per event FastJet",200,0,10);
  //hptBkgFJ->SetDirectory(0); // to avoid root from owning

}

ktFastJet::~ktFastJet()
{
  //DEBUG:
  //cout<<endl;
  //cout<<"Destructor of FastJet ..."<<endl;

  //delete hXiBkgFJ;
  //delete hZBkgFJ;

  //if (hptBkgFJ !=0 )
  //delete hptBkgFJ;
    
  input_particles.clear();
  input_particles_charged.clear();
  input_particles_neutral.clear();

  FFArray->Delete();
  delete FFArray;

  PhiBkgHistos->Delete();
  delete PhiBkgHistos;
 
  //median_ptArray->Delete();
  delete median_ptArray;

}

void ktFastJet::Clear()
{
  input_particles.clear();
  input_particles_charged.clear();
  input_particles_neutral.clear();

  median_pt_per_area=0;
  median_pt_per_area_charged=0;
  median_pt_per_area_neutral=0;

  bkgCalculated=false;

  PhiBkgHistos->Delete();
  PhiBkgHistos->Compress();

  FFArray->Delete();
  FFArray->Compress();
  delete median_ptArray;
  //median_ptArray->Delete();
  // median_ptArray->Compress();

  //DEBUG:
  //cout<<hptBkgFJ<<endl;

  /*
  if (hptBkgFJ !=0 )
    {
      hptBkgFJ->Delete();
      //delete hptBkgFJ;     
      hptBkgFJ=new TH1D("hptBkgFJ","pt dist. of Background per event FastJet",200,0,10);
      //hptBkgFJ->SetDirectory(0); // to avoid root from owning
    }
  */
}

void ktFastJet::Init()
{
  // Dummy !!!
}

void ktFastJet::PrintInfo()
{
  cout<<endl;
  cout<<" FastJet infos:"<<endl;
  cout<<" --------------"<<endl;
  cout<<" |EtaMax| = "<<mEtaMax<<" and  |PhiMax| = "<<mPhiMax<<endl;
  cout<<" pt,cut > "<<mPtCut<<endl;
  cout<<" Cuts for minimizing output volume:"<<endl;
  cout<<" |EtaMax| for analysis = "<<mEtaAna<<endl;
  cout<<" pt,min cut for analysis = "<<mPtAna<<endl;
  PrintAreaDefinition();
}

void ktFastJet::SetFiducial(Double_t mEta,Double_t mPhi)
{
  mEtaMax=mEta;
  mPhiMax=mPhi;
}

void ktFastJet::SetAreaDefinition(Double_t m_ghost_etamax, Int_t m_active_area_repeats, Double_t m_ghost_area)
{
  ghost_etamax=m_ghost_etamax;
  active_area_repeats=m_active_area_repeats;
  ghost_area=m_ghost_area;
}

void ktFastJet::PrintAreaDefinition()
{
  cout<<endl;
  cout<<" FastJet (active) area definition:"<<endl;
  cout<<" ---------------------------------"<<endl;
  cout<<" ghost_etamax = "<<ghost_etamax<<endl;
  cout<<" active_area_repeats = "<<active_area_repeats<<endl;
  cout<<" ghost_area = "<<ghost_area<<endl;
  cout<<" Use passive area for SISCone and median_pt_over_area determined by kt"<<endl;
  cout<<endl;
}

void ktFastJet::PrintFFArray()
{
  for (int i=0;i<FFArray->GetEntries();i++)
    {
      ktParticle* p=(ktParticle*) FFArray->At(i);
      cout<<i<<" "<<p->Eta()<<" "<<p->Phi()<<" "<<p->Pt()<<" "<<p->IsInJet()<<endl;
    }
}

void ktFastJet::ClearFFArray()
{
  // DEBUG:
  //cout<<" Clear FF array ..."<<endl;

  for (int i=0;i<FFArray->GetEntries();i++)
    {
      ktParticle* p=(ktParticle*) FFArray->At(i);
      p->SetIsInJet(false);
    }
}

// just for backward comp.
void ktFastJet::AddParticle(Double_t px, Double_t py, Double_t pz, Double_t E)
{
  TLorentzVector *pv=new TLorentzVector();
  pv->SetPxPyPzE(px,py,pz,E);

  if (TMath::Abs(pv->Eta())<mEtaMax && TMath::Abs(pv->Phi())<mPhiMax && pv->Pt()>mPtCut)
    input_particles.push_back(fastjet::PseudoJet(px,py,pz,E)); 

  delete pv;
}

void ktFastJet::AddParticle(Double_t px, Double_t py, Double_t pz, Double_t E, Bool_t charged)
{

  // DEBUG:
  //cout<<" ---> In : ktFastJet::AddParticle(Double_t px, Double_t py, Double_t pz, Double_t E, Bool_t charged)"<<endl;
  
  ktParticle *pv=new ktParticle();
  pv->SetPxPyPzE(px,py,pz,E);
  ktPID *myPID = new ktPID();
  if( charged){
    myPID->SetCharge(+1);
  }
  else
    myPID->SetCharge(0);

  pv->SetPid(myPID);
  
  if (TMath::Abs(pv->Eta())<mEtaMax && TMath::Abs(pv->Phi())<mPhiMax) // && pv->Pt()>mPtCut)
    {
      //fill all particles in array	
      if (pv->Pt()>0.2)
	{
	  FFArray->AddLast(pv);
	  
	  if (pv->Pt()>mPtCut)
	    {
	      //fastjet::PseudoJet mPart(px,py,pz,E);
	      fastjet::PseudoJet mPart(pv->Px(),pv->Py(),pv->Pz(),pv->E());
	      
	      // user user index to set charged vs.neutral
	      // can be also the index to ktParticle in FFArray
	      // for further studies ...
	      if (charged)
		mPart.set_user_index(1);
	      else
		mPart.set_user_index(0);
	      
	      input_particles.push_back(mPart); 
	      if (charged)
		input_particles_charged.push_back(mPart); 
	      else
		input_particles_neutral.push_back(mPart); 
	      
	      //fill all particles in array	  
	      //FFArray->AddLast(pv);
	    }
	}
      else
	{
	  delete pv;
	  //delete myPID;
	}
    }
  else
    {
      delete pv;
      //delete myPID;
    }

  delete myPID;
}

void ktFastJet::AddParticle(Double_t px, Double_t py, Double_t pz, Double_t E, Bool_t charged, Int_t sign)
{

  // DEBUG:
  //cout<<" ---> In : ktFastJet::AddParticle(Double_t px, Double_t py, Double_t pz, Double_t E, Bool_t charged)"<<endl;
  
  ktParticle *pv=new ktParticle();
  pv->SetPxPyPzE(px,py,pz,E);
  ktPID *myPID = new ktPID();
  if( charged){
    myPID->SetCharge(sign);
  }
  else
    myPID->SetCharge(0);

  pv->SetPid(myPID);
  
  if (smear)
    {
      if (charged)
	{
	  Double_t mdpt=gRandom->Gaus(0,ptRes->Eval(pv->Pt()));
	  Double_t mpts=pv->Pt()+mdpt*pv->Pt();			  
	  
	  pv->SetPtEtaPhiE(mpts,pv->Eta(),pv->Phi(),pv->E());
	}
      else
	{
	  // check if E or Et ...
	  Double_t mde=gRandom->Gaus(0,eRes->Eval(pv->Pt()));
	  Double_t mes=pv->Pt()+mde*pv->Pt();
	  
	  pv->SetPtEtaPhiE(mes,pv->Eta(),pv->Phi(),pv->E());
	}
    }
  
  if (spaceCharge)
    {
      Double_t mptSS=ptSpaceCharge(pv->Pt(),sign);
	  
      //DEBUG:
      //if (pv->Pt()>4)
      //cout<<pv->Pt()<<" "<<mptSS<<" "; //sign<<endl;

      pv->SetPtEtaPhiE(mptSS,pv->Eta(),pv->Phi(),pv->E());
      
      //DEBUG:
      //cout<<pv->Pt()<<" "<<sign<<endl;
    }
  
  if (TMath::Abs(pv->Eta())<mEtaMax && TMath::Abs(pv->Phi())<mPhiMax) // && pv->Pt()>mPtCut)
    {
      //fill all particles in array
      if (pv->Pt()>0.2)
	{
	  FFArray->AddLast(pv);
	  
	  if (pv->Pt()>mPtCut)
	    {
	      //fastjet::PseudoJet mPart(px,py,pz,E);
	      fastjet::PseudoJet mPart(pv->Px(),pv->Py(),pv->Pz(),pv->E());
	      
	      // user user index to set charged vs.neutral
	      // can be also the index to ktParticle in FFArray
	      // for further studies ...
	      if (charged)
		mPart.set_user_index(1);
      else
		mPart.set_user_index(0);
	      
	      input_particles.push_back(mPart); 
	      if (charged)
		input_particles_charged.push_back(mPart); 
	      else
		input_particles_neutral.push_back(mPart); 
	      
	      //fill all particles in array	  
	      //FFArray->AddLast(pv);
	    }
	}
      else
	{
	  delete pv;
	  //delete myPID;
	}
    }
  else
    {
      delete pv;
      //delete myPID;
    }

  delete myPID;
}

void ktFastJet::AddParticle(Double_t px, Double_t py, Double_t pz, Double_t E,Bool_t charged,ktPID* mPid)
{
  // just for backward comp.
  // in case of ktPID filled bool charged flag obsolete ...
  AddParticle(px,py,pz,E,mPid);
}

void ktFastJet::AddParticle(Double_t px, Double_t py, Double_t pz, Double_t E,ktPID* mPid)
{

  // DEBUG:
  //cout<<" In:  ktFastJet::AddParticle(Double_t px, Double_t py, Double_t pz, Double_t E, Bool_t charged, ktPID* mPid)"<<endl;

  ktParticle *pv=new ktParticle();
  pv->SetPxPyPzE(px,py,pz,E);
  pv->SetPid(mPid);

  // Test
  if( fabs(mPid->GetPID()) == 3122 || fabs(mPid->GetPID()) == 310 ) cout << "Seen pid is " << pv->GetPid()->GetPID() << endl;

  if (TMath::Abs(pv->Eta())<mEtaMax && TMath::Abs(pv->Phi())<mPhiMax) // && pv->Pt()>mPtCut)
    {
      //fill all particles in array
      if (pv->Pt()>0)
	{
	  FFArray->AddLast(pv);
	  
	  if (pv->Pt()>mPtCut)
	    {
	      fastjet::PseudoJet mPart(px,py,pz,E);
	      // user user index to set charged vs.neutral
	      // can be also the index to ktParticle in FFArray
	      // for further studies ...
	      if (mPid->IsCharged())
		mPart.set_user_index(1);
	      else
		mPart.set_user_index(0);
	      
	      input_particles.push_back(mPart); 
	      if (mPid->IsCharged())
		input_particles_charged.push_back(mPart); 
	      else
		input_particles_neutral.push_back(mPart); 
	      

	      // test
	       if( fabs(pv->GetPid()->GetPID())  == 3122 ||
	      	  fabs(pv->GetPid()->GetPID())  == 310 ) cout << "aDDING V0 FF Array " <<  pv->GetPid()->GetPID() << endl;
	      //fill all particles in array
	      //FFArray->AddLast(pv);
	    }
	}
      else{
	if( fabs(pv->GetPid()->GetPID())  == 3122 ||
	    fabs(pv->GetPid()->GetPID())  == 310 ) cout << "Rejecting V0 FF Array pt cut" <<  pv->GetPid()->GetPID() << endl;
	delete pv;
      }
    }
  else{
    if( fabs(pv->GetPid()->GetPID())  == 3122 ||
	fabs(pv->GetPid()->GetPID())  == 310 ) cout << "rejecting  V0 FF Array eta phi cut" <<  pv->GetPid()->GetPID() << endl;
    delete pv;
  }
}


void ktFastJet::AddQuenchedJet(ktJetQuench *mQ)
{

  if (mQ!=0)
    {
      for (int i=0;i<mQ->GetNJetParticles();i++)
	{
	  Bool_t charged=false;
	  ktParticle *p=(ktParticle*) mQ->GetQuenchedParticle(i);
	  if (p->GetPid()>0)
	    charged=true;
	  
	  AddParticle(p->Px(),p->Py(),p->Pz(),p->E(),charged,p->GetPid());
	}
    }
  else
    cout<<" ---> Warning ktJetQuench zero pointer in ktFastJet !!!'"<<endl;

}

void ktFastJet::AddUnQuenchedJet(ktJetQuench *mQ)
{

  if (mQ!=0)
    {
      for (int i=0;i<mQ->GetNJetParticlesNoQuench();i++)
	{
	  Bool_t charged=false;
	  ktParticle *p=(ktParticle*) mQ->GetUnQuenchedParticle(i);
	  if (p->GetPid()>0)
	    charged=true;
	  AddParticle(p->Px(),p->Py(),p->Pz(),p->E(),charged,p->GetPid());
	}
    }
  else
    cout<<" ---> Warning ktJetQuench zero pointer in ktFastJet !!!'"<<endl;

}

void ktFastJet::AddPythia8Event(ktPythia8 *myp,TString mSel)
{
   // Fill grid with Pythia8
  // Selections : TPC (charged only), EMCAL (charged+neutral) and ALL (all final state)

  Pythia *p=(Pythia*) myp->GetPythia();

  if (info)
    {
      cout<<endl;
      cout<<" ---> Fill FastJet with Pythia8 event ..."<<endl;
      cout<<" ---> Number of particles = "<< p->event.size()<<endl;
      //cout<<" ---> Number of final state particles = "<< p->info.nFinal()<<endl;
      cout<<" ---> Particle selection = "<<mSel<<endl;
      cout<<" ---> |EtaMax| = "<<mEtaMax<<" and  |PhiMax| = "<<mPhiMax<<endl;
      cout<<" ---> pt,cut > "<<mPtCut<<endl;
      cout<<endl;
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
	  Int_t sign=(Int_t) p->event[i].charge();

	  TLorentzVector *pv=new TLorentzVector();
	  pv->SetPxPyPzE(px,py,pz,E);
	  
	  // DEBUG:
	  //if(p->event[i].id()==111) cout<<"pi0"<<endl;
	  // check if decay on ...

	  if (TMath::Abs(pv->Eta())<mEtaMax && TMath::Abs(pv->Phi())<mPhiMax)
	    {
	      
	      if (mSel=="TPC")
		{
		  if (p->event[i].isCharged())
		    {
		      //input_particles.push_back(fastjet::PseudoJet(px,py,pz,E)); 
		      AddParticle(px,py,pz,E);
		      NSel++;
		    }
		}
	      else if (mSel=="EMCAL")
		{
		  //DEBUG:
		  //cout<<"in emcal ..."<<endl;
		  
		  if (p->event[i].isCharged())
		    {
		      //input_particles.push_back(fastjet::PseudoJet(px,py,pz,E)); 
		      AddParticle(px,py,pz,E,true,sign);
		      NSel++;
		    }
		  else if (p->event[i].id()==22 || p->event[i].id()==111)
		    {
		      //input_particles.push_back(fastjet::PseudoJet(px,py,pz,E)); 
		      AddParticle(px,py,pz,E,false,sign);
		      NSel++;
		    }
		}
	      else
		{
		  //input_particles.push_back(fastjet::PseudoJet(px,py,pz,E)); 
		  AddParticle(px,py,pz,E);
		  NSel++;
		}
	    }

	  delete pv; 
	}
    }

  if (info)
    cout<<" ---> Number of particles filled = "<<NSel<<endl;
  //cout<<endl;
}

Double_t ktFastJet::GetdPhi(Double_t mphi,Double_t vphi)
{
  if (vphi < -1*TMath::Pi()) vphi += (2*TMath::Pi());
  else if (vphi > TMath::Pi()) vphi -= (2*TMath::Pi());
  if (mphi < -1*TMath::Pi()) mphi += (2*TMath::Pi());
  else if (mphi > TMath::Pi()) mphi -= (2*TMath::Pi());
  double dphi = mphi-vphi;
  if (dphi < -1*TMath::Pi()) dphi += (2*TMath::Pi());
  else if (dphi > TMath::Pi()) dphi -= (2*TMath::Pi());

  return dphi;
}

Bool_t ktFastJet:: IsTriggerJet(ktMuFastJet *mJ,ktMuEvent *ev)
{
  Double_t jEta=mJ->Eta();
  Double_t jPhi=mJ->Phi();

  Bool_t isTrigger=false;

  if (ev->GetHTEta()>-10 && ev->GetHTPhi()>-10)
    {
      Double_t tEta=ev->GetHTEta();
      Double_t tPhi=ev->GetHTPhi();
      
      Double_t dphi=GetdPhi(jPhi,tPhi);
      Double_t deta=jEta-tEta;

      Double_t d=TMath::Sqrt(deta*deta+dphi*dphi);
      
      if (d<0.3) // check ... or set as variable ...
	isTrigger=true;
    }

  if (ev->GetJPEta()>-10 && ev->GetJPPhi()>-10)
    {
      Double_t tEta=ev->GetJPEta();
      Double_t tPhi=ev->GetJPPhi();

      Double_t dphi=GetdPhi(jPhi,tPhi);
      Double_t deta=jEta-tEta;

      //Double_t d=TMath::Sqrt(deta*deta+dphi*dphi);

      if (dphi<0.5 && deta<0.5) // check !!!
	isTrigger=true;
    }

  mJ->SetIsTriggerJet(isTrigger);

  return isTrigger;
}


void ktFastJet::GetJetFF(ktMuFastJet *mJ,Bool_t setInJet, Bool_t mSave)
{
  TH1D *hptJet=new TH1D("hptJet","pt dist. in jet FF",600,0,30);
  hptJet->SetDirectory(0);
 
  Double_t meta=mJ->Eta();
  Double_t mphi=mJ->Phi();

  // cout << "FF Array has " << FFArray->GetEntries() << endl;
  for (int i=0;i<FFArray->GetEntries();i++)
    {
      ktParticle *v=(ktParticle*) FFArray->At(i);
      Double_t veta=v->Eta();      
      Double_t vphi=v->Phi();

      Double_t dphi=GetdPhi(mphi,vphi);

      Double_t d=sqrt((meta-veta)*(meta-veta)+dphi*dphi);
     
      if (d<RcFF)
	{
	  if (setInJet)
	    v->SetIsInJet(true);
	  
	  if (mSave)
	    {
	      if (!histFFonly)
		{
		  //  cout << "For " << mJ->Pt() << "Adding " << v->Pt() << " " << d << " dphi :" << dphi << endl;
		  mJ->AddJetParticle((ktParticle*) v);
		}
	      else
		{
		  // do not use e- in charged FF
		  // check for simulations (?); should work on data (new picoDsts)
		  if (v->IsCharged() && v->GetPid()->GetPID()!=11)
		    {
		      hptJet->Fill(v->Pt());
		    }
		}
	    }	
	}
      else
	{
	  //	  cout << "For " << mJ->Pt() << " rejecting " << v->Pt() << " " << d << " dphi :" << dphi << " " << mJ->Eta() << " " << v->Eta() << " " << mJ->Phi() << " " << v->Phi() << endl;
	}
    }
  //  cout << "My jet now has " << mJ->GetNJetParticles() << " jet pt " << mJ->Pt() <<   " " << mJ->PtCorr() <<endl;
  if (histFFonly && mSave)
    mJ->SetJetPtHistogram(hptJet);
  
  delete hptJet;
  
}

void ktFastJet::FillFF(TString mAlgo, ktMuEvent *ev) //, Int_t nJet)
{
  if (!ev==0 && fillFF)
    {
      if (info)
	{
	  cout<<" ---> Fill FF for "<<mAlgo<<" for the "<<nFF<<" highest jets ..."<<endl;
  cout<<endl;

	  // DEBUG:
	  //cout<<" ---> "<<FFArray->GetEntries()<<endl;
	}
      
      Int_t nMax=0;
      Int_t nFoundJets=0;
      Int_t nFilled=0;

      if (mAlgo.Contains("Anti"))	    
	nFoundJets=ev->GetNFJAntiKt();	  
      else if (mAlgo.Contains("SIS"))	    
	nFoundJets=ev->GetNFJSISCone();	  
      else	    
	nFoundJets=ev->GetNFJkt();	  
      
      if (nFF>nFoundJets)
	nMax=nFoundJets;
      else
nMax=nFF;	      	  	  
      
      // DEBUG:
      //cout<<" ---> "<<nFoundJets<<endl;
      
      // include cuts to reduce output volume at this point !
      // list sorted by pt,corr ...
      //for (int n=0;n<nMax;n++)
      for (int n=0;n<nFoundJets;n++)
	{
	 	 
	  if (!(nFilled<nMax)) break;
	  
	  ktMuFastJet *mJ=0;	  

	  if (mAlgo.Contains("Anti"))
	    mJ=(ktMuFastJet*) ev->GetFJAntiKt(n);
	  else if (mAlgo.Contains("SIS"))
	    mJ=(ktMuFastJet*) ev->GetFJSISCone(n);
	  else		
	    mJ=(ktMuFastJet*) ev->GetFJkt(n);	 	
	  
	  //mJ->SetFFAreaInAcc(AreaInAcc(RcFF,meta,mphi,mEtaMax,(2*mPhiMax)));
	  // change max phi in area calculation arbitrarily high to remove phi border
	  mJ->SetFFAreaInAcc(AreaInAcc(RcFF,mJ->Eta(),mJ->Phi(),mEtaMax,(10*mPhiMax)));

	  Bool_t isTrig= IsTriggerJet(mJ,ev);

	  // Fill/save FF for particular jet fullfilling analysis cuts 	
	  // mark/remove the particles for FF bkg. for the nFF highest jets only !
	  // think how to handle if trigger jet is not in the NFF highest ...
	  if (TMath::Abs(mJ->Eta())<mEtaAna && mJ->PtCorr()>mPtAna)
	    {
	      if (n<nMax)
		{
		  GetJetFF(mJ,true,true);
		}
	      else{
		GetJetFF(mJ,false,true);
	      }
	      nFilled++;
	    }
	  else
	    {
	      if (n<nMax){
		GetJetFF(mJ,true,false);
	      }
	    }
	}
    }
  else
    {
      if (fillFF)
	{
	  cout<<endl;
	  cout<<" ---> Warning : FF can  not be filled. Need to fill ktMuEvent before !!!"<<endl;
	  cout<<endl;
	}
    }
}


// This needs to be generalized for PIDed particles its currently all charged only 
void ktFastJet::DoXiBkg(TString mAlgo,ktMuEvent *ev)
{

  if (info)
    cout << " Warning: you are only calc. Xi bkg for all charged particles " << endl;

  if (!ev==0 && fillFF)
    {
      if (info)
	{
	  cout<<" ---> Calculate Xi bkg for "<<mAlgo<<" by removing the "<<nFF<<" highest jets ..."<<endl;
	  cout<<endl;	
	}

      TH1D *hptBkgFJ=new TH1D("hptBkgFJ","pt dist. of Background per event FastJet",200,0,10);
      hptBkgFJ->SetDirectory(0); // to avoid root from owning

      // DEBUG:
      //PrintFFArray();
      Int_t nOut=0;

      for (int i=0;i<FFArray->GetEntries();i++)
	{
	  ktParticle *v=(ktParticle*) FFArray->At(i);
	  
	  // veto against electrons to be consistent with charged FF definition in FillFF(...)
	  if (!v->IsInJet() && v->IsCharged() && v->GetPid()->GetPID()!=11)
	    {
	      if (!v || v->Pt()>0) // check why !!!!???? happens before !!???
		{
		  hptBkgFJ->Fill(v->Pt());
		}
	      
	      nOut++;
	    }
	}
      
      // DEBUG:
      //cout<<"DBEUG: Particles outside the 2 highest jets = "<<nOut<<endl;      

      // Get area occupied by the 2 highest jets ...
      Double_t mJetArea=0;
      Double_t mBkgArea=2*mEtaMax*2*mPhiMax;

      Int_t nMax=0;
      Int_t nFoundJets=0;
      
      if (mAlgo.Contains("Anti"))	    
	nFoundJets=ev->GetNFJAntiKt();	  
      else if (mAlgo.Contains("SIS"))	    
	nFoundJets=ev->GetNFJSISCone();	  
      else	    
	nFoundJets=ev->GetNFJkt();	  
      
      if (nFF>nFoundJets)
	nMax=nFoundJets;
      else
	nMax=nFF;	      
      
      for (int n=0;n<nMax;n++)
	{
	  
	  ktMuFastJet *mJ=0;
	  
	  if (mAlgo.Contains("Anti"))
	    mJ=(ktMuFastJet*) ev->GetFJAntiKt(n);
	  else if (mAlgo.Contains("SIS"))
	    mJ=(ktMuFastJet*) ev->GetFJSISCone(n);
	  else		
	    mJ=(ktMuFastJet*) ev->GetFJkt(n);

	  mJetArea += mJ->GetFFAreaInAcc();
	  
	  // DEBUG:
	  //cout<<mJ->GetNJetParticles()<<" "<<mJ->GetFFAreaInAcc()<<endl;
	}

      mBkgArea -= mJetArea;   
      hptBkgFJ->Scale(1/mBkgArea);

      if (mAlgo.Contains("Anti"))
	{
	  ev->SetptBkgAntiKt(hptBkgFJ);
	}
      else if (mAlgo.Contains("SIS"))
	{
	  ev->SetptBkgSISCone(hptBkgFJ);
	}
      else
	{
	  ev->SetptBkgKt(hptBkgFJ);
	}
      
      delete hptBkgFJ;
    }
  else
    {
      if (fillFF)
	{
	  cout<<endl;
	  cout<<" ---> Warning : Xi bkg not calculated. Need to fill ktMuEvent before !!!"<<endl;
	  cout<<endl;				
	}
    }
  
  ClearFFArray();
}

void ktFastJet::DoPhiBkg(TString mAlgo,ktMuEvent *ev)
{

  if (!ev==0)
    {
      //if (info)
	{
	  cout<<" ---> Calculate dPhi dep. background ..."<<endl;
	  cout<<endl;	
	}

      Int_t nMax=0;
      Int_t nFoundJets=0;
      
      if (mAlgo.Contains("Anti"))	    
	nFoundJets=ev->GetNFJAntiKt();	  
      else if (mAlgo.Contains("SIS"))	    
	nFoundJets=ev->GetNFJSISCone();	  
      else	    
	nFoundJets=ev->GetNFJkt();	  
      
      if (nFF>nFoundJets)
	nMax=nFoundJets;
      else
	nMax=nFF;	      
      
      for (int n=0;n<nMax;n++)
	{	 

	  ktMuFastJet *mJ=0;
	  
	  if (mAlgo.Contains("Anti"))
	    mJ=(ktMuFastJet*) ev->GetFJAntiKt(n);
	  else if (mAlgo.Contains("SIS"))
	    mJ=(ktMuFastJet*) ev->GetFJSISCone(n);
	  else		
	    mJ=(ktMuFastJet*) ev->GetFJkt(n);
	  
	  Double_t jEta=mJ->Eta();
	  Double_t jPhi=mJ->Phi();

	  // DEBUG:
	  cout<<n<<" "<<jEta<<" "<<jPhi<<" "<<mJ->PtCorr()<<endl;

	  //TH2D *phiHist=new TH2D("phiHist","Energy dR vs. dPhi",40,0,4,45,-TMath::Pi(),TMath::Pi());
	  TH2D *phiHist=new TH2D("phiHist","Energy dR vs. dPhi",40,0,2,45,-TMath::Pi(),TMath::Pi());

	  phiHist->Sumw2();
	  phiHist->SetDirectory(0); // to avoid root from owning 
	  
	  for (int i=0;i<FFArray->GetEntries();i++)
	    {
	      ktParticle *v=(ktParticle*) FFArray->At(i);
	      // DEBUG:
	      //cout<<v->Eta()<<" "<<v->Phi()<<endl;
	      
	      Double_t veta=v->Eta(); 
	      Double_t deta=jEta-veta;

	      Double_t vphi;
	      vphi=v->Phi();if (vphi<0) vphi += (2*TMath::Pi());
	      Double_t dphi=jPhi-vphi;

	      Double_t R=sqrt(deta*deta+dphi*dphi);
	
	      //phiHist->Fill(R,dphi,v->Pt());
	      phiHist->Fill(fabs(deta),dphi,v->Pt());
	    }

	  PhiBkgHistos->AddLast(phiHist);
	}
    }
}

void ktFastJet::RunSISCone(Double_t Rc,ktMuEvent *ev)
{

  // define a pluginpoint
  fastjet::JetDefinition::Plugin * plugin;

  // allocate a new plugin
  double overlap_threshold = 0.75;
  plugin = new fastjet::SISConePlugin (Rc, overlap_threshold,0,0.0);
  // problems if protojet pt,cut > 0 with active area ....  // see manual and use passive or direct from SISCone !???

  // create a jet-definition based on the plugin
  fastjet::JetDefinition jet_def(plugin);

  // see mail from gregory !!!! check passive/active area SISCone !!!
  fastjet::ActiveAreaSpec area_spec(ghost_etamax, active_area_repeats, ghost_area);
  // run the jet clustering with the above jet definition
  
  fastjet::AreaDefinition area_def_cone(fastjet::passive_area,area_spec);
  // cout << "input particles " << input_particles.size() << endl;
  if (input_particles.size()>0)
    {
      
      fastjet::ClusterSequenceArea clust_seq(input_particles,jet_def,area_def_cone);
      
      // tell the user what was done
      if (info)
	{
	  cout << " ---> Ran " << jet_def.description() << endl;
	  cout<<endl;      
	  cout << "Strategy adopted by FastJet was "<<clust_seq.strategy_string()<<endl<<endl;
	}
      
      vector<fastjet::PseudoJet> inclusive_jets = clust_seq.inclusive_jets(ptmin);
      
      // print them out
      if (printjet)
	print_jets_sub(clust_seq, inclusive_jets, median_pt_per_area);
      
      /*
      if (ev!=0)
	{
	  if (!fillFF)
	    fill_jets_in_event(clust_seq, inclusive_jets, median_pt_per_area,ev,"SISCone");
	  else
	    fill_jets_in_event(clust_seq, inclusive_jets, median_pt_per_area,ev,"SISCone",FFArray,RcFF,nFF,mEtaMax,2*mPhiMax);
	}
      */

      fill_jets_in_event(clust_seq, inclusive_jets, median_pt_per_area,ev,"SISCone",median_pt_per_area_charged,median_pt_per_area_neutral);//,mEtaAna,mPtAna);

      // not recomended by Gregory in conncetion with SISCone !
      //fastjet::ClusterSequenceActiveArea clust_seq2(input_particles,jet_def, area_spec);
      //print_jets_sub(clust_seq2, inclusive_jets, median_pt_per_area);
    }
  else
    cout<<"  ---> Warning (ktFastJet) : Empty event !!!"<<endl;
  
  delete plugin;

}


void ktFastJet::DoBkgStrip(Double_t Rparam,TString mAlgo)
{ 
  fastjet::Strategy strategy = fastjet::Best;
  //fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme;  
  fastjet::JetAlgorithm fjAlgoBkg;
  //fastjet::GhostedAreaSpec area_spec(ghost_etamax);


  if (mAlgo.Contains("AntiKt"))
    fjAlgoBkg= fastjet::antikt_algorithm;
  else
    fjAlgoBkg= fastjet::kt_algorithm;
  
  fastjet::JetDefinition jet_def_bkg(fastjet::kt_algorithm, 0.4, strategy);
  //fastjet::ActiveAreaSpec area_spec_bkg(ghost_etamax, active_area_repeats,ghost_area);  
  fastjet::AreaDefinition area_def_bkg(active_area_explicit_ghosts, 
				       fastjet::GhostedAreaSpec(ghost_etamax));
  cout<<" input particle size  "<<input_particles.size()<<endl;
  

  fastjet::JetDefinition jet_def_rho(fastjet::kt_algorithm, Rparam, strategy);
  fastjet::ActiveAreaSpec area_spec_rho(ghost_etamax, active_area_repeats, ghost_area);


  float Rin=Rparam;
  float Rout=3*Rparam;
  

 if (input_particles.size()>0)
    {
      
      fastjet::ClusterSequenceArea clust_seq(input_particles, jet_def_bkg, area_def_bkg);
      //  fastjet::Selector selector = fastjet::SelectorAbsRapMax(5) * (!fastjet::SelectorNHardest(2);//removes the 2 hardest jets in a hardwired way
      fastjet::ClusterSequenceActiveArea clust_seq_rho(input_particles, jet_def_rho, area_spec_rho);
      // tell the user what was done
      if (info)
	{
	  cout << " ---> Ran " << jet_def_rho.description() << endl;
	  cout<<endl;      
	  cout << "Strategy adopted by FastJet was "<<clust_seq.strategy_string()<<endl<<endl;
	}     
      
      vector<fastjet::PseudoJet> inclusive_jets_rho = clust_seq_rho.inclusive_jets(ptmin);
      

      // vector<fastjet::PseudoJet> inclusivejets = clust_seq.inclusive_jets(ptmin);
      //  vector<fastjet::PseudoJet> inclusivejets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));
      //check the fastjet manual for other background subtraction methods
      fastjet::Selector selector = fastjet::SelectorDoughnut(Rin, Rout);
      fastjet::JetMedianBackgroundEstimator bkg_estimator(selector, jet_def_bkg, area_def_bkg);      
      fastjet::Subtractor subtractor(&bkg_estimator);
      bkg_estimator.set_particles(input_particles);


      vector<fastjet::PseudoJet> jets_rho = sorted_by_pt(inclusive_jets_rho);  
      // the corrected jets will go in here
      // vector<fastjet::PseudoJet> corrected_jets_rho(jets_rho.size());
  
      for (Size_t j = 0; j < jets_rho.size(); j++) {
	
      median_pt_per_area = bkg_estimator.rho(jets_rho[j]);
      cout<<"New median pt values"<< median_pt_per_area<<endl;
      median_ptArray->SetAt(median_pt_per_area,j);
      }
    }
  else
    cout<<"  ---> Warning (ktFastJet::DoBkg) : Empty event !!!"<<endl;
  
  

  if (info)
    {
      cout<<" ---> Calculate background : R = "<<Rparam<<" with "<<mAlgo<<endl;
      cout<<" ---> MedianPerPtArea         = "<<median_pt_per_area<<endl;
        }




    bkgCalculated=true;
 

}
void ktFastJet::DoBkgNew(Double_t Rparam,TString mAlgo)
{
 
  fastjet::Strategy strategy = fastjet::Best;
  //fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme;
  
  fastjet::JetAlgorithm fjAlgoBkg;

  // fastjet::JetDefinition::Plugin * plugin;
  // allocate a new plugin
  // plugin = new fastjet::SISConePlugin (Rc, overlap_threshold,0,0.0);
  // problems if protojet pt,cut > 0 with active area ....



  //fastjet::GhostedAreaSpec area_spec(ghost_etamax);
  //  virtual void ;

  if (mAlgo.Contains("AntiKt"))
    fjAlgoBkg= fastjet::antikt_algorithm;
  else
    fjAlgoBkg= fastjet::kt_algorithm;
  
  fastjet::JetDefinition jet_def_bkg(fjAlgoBkg, Rparam, strategy);
  //fastjet::ActiveAreaSpec area_spec_bkg(ghost_etamax, active_area_repeats,ghost_area);
  
  fastjet::AreaDefinition area_def_bkg(active_area_explicit_ghosts, 
				       fastjet::GhostedAreaSpec(ghost_etamax));
  cout<<" input particle size  "<<input_particles.size()<<endl;
  float Rin=Rparam;
  float Rout=3*Rparam;
  
 if (input_particles.size()>0)
    {
      fastjet::ClusterSequenceArea clust_seq(input_particles, jet_def_bkg, area_def_bkg);
      fastjet::Selector selector = fastjet::SelectorAbsRapMax(5) * (!fastjet::SelectorNHardest(2));
      // estimator=new fastjet::JetMedianBackgroundEstimator(selector, jet_def_bkg, area_def_bkg);
      // vector<fastjet::PseudoJet> inclusivejets = clust_seq.inclusive_jets(ptmin);
      //  vector<fastjet::PseudoJet> inclusivejets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));
      
      //   fastjet::Selector selector = fastjet::SelectorDoughnut(Rin, Rout);

      fastjet::JetMedianBackgroundEstimator bkg_estimator(selector, jet_def_bkg, area_def_bkg);      
      fastjet::Subtractor subtractor(&bkg_estimator);
      bkg_estimator.set_particles(input_particles);
      median_pt_per_area = bkg_estimator.rho();
      cout<<"New median pt values"<< median_pt_per_area<<endl;
    }
  else
    cout<<"  ---> Warning (ktFastJet::DoBkg) : Empty event !!!"<<endl;
  
  

  if (info)
    {
      cout<<" ---> Calculate background : R = "<<Rparam<<" with "<<mAlgo<<endl;
      cout<<" ---> MedianPerPtArea         = "<<median_pt_per_area<<endl;
        }

    bkgCalculated=true;
 

}
void ktFastJet::DoBkg(Double_t Rparam,TString mAlgo)
{
 
  fastjet::Strategy strategy = fastjet::Best;
  //fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme;
  
  fastjet::JetAlgorithm fjAlgoBkg;
  
  if (mAlgo.Contains("AntiKt"))
    fjAlgoBkg= fastjet::antikt_algorithm;
  else
    fjAlgoBkg= fastjet::kt_algorithm;
  
  fastjet::JetDefinition jet_def_bkg(fjAlgoBkg, Rparam, strategy);
  fastjet::ActiveAreaSpec area_spec_bkg(ghost_etamax, active_area_repeats, ghost_area);


   
  //  Subtractor subtractor(&bkgd_estimator);
  //  bkgd_estimator.set_particles(input_particles);
 
  // DEBUG
  //cout<<input_particles.size()<<endl;

 
  if (input_particles.size()>0)
    {
      fastjet::ClusterSequenceActiveArea clust_seq_bkg(input_particles, 
						   jet_def_bkg, area_spec_bkg);
      
      median_pt_per_area = clust_seq_bkg.pt_per_unit_area();
      
      if (input_particles_charged.size()>0)
	{
	  fastjet::ClusterSequenceActiveArea clust_seq_bkg_charged(input_particles_charged, 
							   jet_def_bkg, area_spec_bkg);
	  
	  median_pt_per_area_charged = clust_seq_bkg_charged.pt_per_unit_area();
	}
      
      if (input_particles_neutral.size()>0)
	{
	  fastjet::ClusterSequenceActiveArea clust_seq_bkg_neutral(input_particles_neutral, 
								   jet_def_bkg, area_spec_bkg);
	  
	  median_pt_per_area_neutral = clust_seq_bkg_neutral.pt_per_unit_area();
	}
      
    }
  else
    cout<<"  ---> Warning (ktFastJet3::DoBkg) : Empty event !!!"<<endl;
  
  

  if (info)
    {
      cout<<" ---> Calculate background : R = "<<Rparam<<" with "<<mAlgo<<endl;
      cout<<" ---> MedianPerPtArea         = "<<median_pt_per_area<<endl;
      cout<<" ---> MedianPerPtArea charged = "<<median_pt_per_area_charged<<endl;
      cout<<" ---> MedianPerPtArea neutral = "<<median_pt_per_area_neutral<<endl;
      cout<<" ---> Sum = "<<median_pt_per_area_charged+median_pt_per_area_neutral<<" | Ratio C/ALL = "<<median_pt_per_area_charged/median_pt_per_area<<endl;
    }
 
    bkgCalculated=true;
 
}


void ktFastJet::RunFastJetSub(Double_t Rparam,TString mAlgo,ktMuEvent *ev)
{

  //double Rparam = 1.0;
  fastjet::Strategy strategy = fastjet::Best;
  //fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme;

  fastjet::JetAlgorithm fjAlgo;

  if (mAlgo.Contains("AntiKt"))
    fjAlgo= fastjet::antikt_algorithm;
  else
    fjAlgo= fastjet::kt_algorithm;

  fastjet::JetDefinition jet_def(fjAlgo, Rparam, strategy);
  fastjet::ActiveAreaSpec area_spec(ghost_etamax, active_area_repeats, 
                                    ghost_area);

  // DEBUG:
  //cout<<input_particles.size()<<endl;

  if (input_particles.size()>0)
    {
      // DEBUG
      //cout<<input_particles.size()<<endl;

      // run the jet clustering with the above jet definition
      fastjet::ClusterSequenceActiveArea clust_seq(input_particles, 
						   jet_def, area_spec);
      
      // tell the user what was done
      if (info)
	{
	  cout << " ---> Ran " << jet_def.description() << endl;
	  cout<<endl;      
	  cout << "Strategy adopted by FastJet was "<<clust_seq.strategy_string()<<endl<<endl;
	}     
      
      vector<fastjet::PseudoJet> inclusive_jets = clust_seq.inclusive_jets(ptmin);
      /*
      if (mAlgo.Contains("kt") && !bkgCalculated)
	{
	  median_pt_per_area = clust_seq.pt_per_unit_area();
	  bkgCalculated=true;
	  
	  //DEBUG:
	  //cout<<mAlgo<<endl;
	  if (info)
	    cout<<" ---> MedianPerPtArea (standard) = "<<median_pt_per_area<<endl;
	}
      */

      // print them out
      if (printjet)
	print_jets_sub(clust_seq, inclusive_jets, median_pt_per_area);
      
      fill_jets_in_event(clust_seq, inclusive_jets, median_pt_per_area,ev,mAlgo,median_pt_per_area_charged,median_pt_per_area_neutral);//,mEtaAna,mPtAna);
      
    }
  else
    cout<<"  ---> Warning (ktFastJet) : Empty event !!!"<<endl;
}

Bool_t ktFastJet::GetHighestJetForQuenching(ktJetQuench *mQ,Double_t Rparam,TString mAlgo)
{
  
  Bool_t foundJet=false;

  if (info)
    cout<<"---> Get highest jet (pythia) for quenching ..."<<endl;

  //double Rparam = 1.0;
  fastjet::Strategy strategy = fastjet::Best;
  //fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme;

  fastjet::JetAlgorithm fjAlgo;

  if (mAlgo.Contains("AntiKt"))
    fjAlgo= fastjet::antikt_algorithm;
  else
    fjAlgo= fastjet::kt_algorithm;

  fastjet::JetDefinition jet_def(fjAlgo, Rparam, strategy);
  fastjet::ActiveAreaSpec area_spec(ghost_etamax, active_area_repeats, 
                                    ghost_area);

  // DEBUG:
  //cout<<input_particles.size()<<endl;

  if (input_particles.size()>0 && mQ!=0)
    {

      // run the jet clustering with the above jet definition
      fastjet::ClusterSequenceActiveArea clust_seq(input_particles, 
						   jet_def, area_spec);
      
      // tell the user what was done
      if (info)
	{
	  cout << " ---> Ran " << jet_def.description() << endl;
	  cout<<endl;      
	  cout << "Strategy adopted by FastJet was "<<clust_seq.strategy_string()<<endl<<endl;
	}     
      
      vector<fastjet::PseudoJet> inclusive_jets = clust_seq.inclusive_jets(ptmin);

      // fill ktJetQuench class with highest found jet ...

      // sort jets into increasing pt
      vector<fastjet::PseudoJet> jets = sorted_by_pt(inclusive_jets);  
      
      if (jets.size()>0)
	{
	  if (printjet)
	    print_jets_sub(clust_seq, inclusive_jets, 0.0);

	  // Get constituents ...
	  vector<fastjet::PseudoJet> constituents=clust_seq.constituents(jets[0]);

	  // DEBUG:
	  //cout<<"N const. of highest jet = "<<constituents.size()<<endl;

	  mQ->SetJetPt(jets[0].perp());
	  mQ->SetJetEta(jets[0].rap());
	  mQ->SetJetPhi(jets[0].phi());

	  for (int i=0;i<(int) constituents.size();i++)
	    {
	      fastjet::PseudoJet mPart=constituents[i];	     
	      mQ->AddParticle(mPart.px(),mPart.py(),mPart.pz(),mPart.E(),mPart.user_index());	      
	    }	  
	  
	  foundJet=true;
	  
	}
      else
	cout<<"  ---> Warning (ktFastJet) : No jet for quenching found !!!"<<endl; 
      
    }
  else
   cout<<"  ---> Warning (ktFastJet) : Empty event !!!"<<endl; 

  return foundJet;
}
