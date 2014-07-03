// first test (Joern Putschke)

#include "ktMuJet.h"
#include <Riostream.h>

ClassImp(ktMuJet)

ktMuJet::ktMuJet()
{

  // set default values for completeness !!!

  // DEBUG:
  //cout<<"Default ktMuJet constructor "<<endl;

}

ktMuJet::ktMuJet(ktMuJet *j,ktMuJet *j2)
{
  // Add constructor for jet clustering !

  JetParticles=new TObjArray(0);
  Jet=new TLorentzVector();
  
  TLorentzVector *jv=(TLorentzVector*) j->GetJetVector();
  TLorentzVector *j2v=(TLorentzVector*) j2->GetJetVector();

  *Jet = *jv;
  *Jet += *j2v;

  NCellsInEmcal=j->GetNCellsInEmcal()+j2->GetNCellsInEmcal();
  NCellsOutsideEmcal=j->GetNCellsOutsideEmcal()+j2->GetNCellsOutsideEmcal();

  JetPtCellCorr=j->PtCellCorr()+j2->PtCellCorr();

   // Fill FF  ....
  NJetParticles=0;
  NSubJets=0;
  NJetCells=j->GetNCells()+j2->GetNCells();

  // Modify according to ktMuJet class !!!
  /*
  for (int i=0;i<j->GetNJetCells();i++)
    {
      ktJetCell *myCell=(ktJetCell*) j->GetJetCellFF(i);
      
      //cout<<i<<" "<<myCell->NFFParticles()<<endl;
      
      for (int k=0;k<myCell->NFFParticles();k++)
	{
	  JetParticles->AddLast((TLorentzVector*) (myCell->GetFFParticleList()->At(k)->Clone()));
	} 
    }
  
  for (int i=0;i<j2->GetNJetCells();i++)
    {
      ktJetCell *myCell=(ktJetCell*) j2->GetJetCellFF(i);
      
      //cout<<i<<" "<<myCell->NFFParticles()<<endl;
      
      for (int k=0;k<myCell->NFFParticles();k++)
	{
	  JetParticles->AddLast((TLorentzVector*) (myCell->GetFFParticleList()->At(k)->Clone()));
	} 
    }
  */

  type="ConeCluster";//j->GetType();
  
}

ktMuJet::ktMuJet(ktJet *j)
{
    //DEBUG:
   //j->PrintJet();
  electronSeed=false;
  NeutralEnergy=0.0;
  ChargedEnergy=0.0;

  JetAreaScale=j->GetJetAreaScale();
  FFAreaScale=j->GetFFAreaScale();

  JetParticles=new TObjArray(0);
  Jet=new TLorentzVector();

  Jet->SetPxPyPzE(j->GetJetVector()->Px(),j->GetJetVector()->Py(),j->GetJetVector()->Pz(),j->GetJetVector()->Energy());

  NCellsInEmcal=j->NJetCellsInEmcal();
  NCellsOutsideEmcal=j->NJetCellsOutsideEmcal();

  JetPtCellCorr=j->Pt();
  JetPtConeCorr=j->Pt();//0;
  JetPtClusterCorr=0;

  // Fill FF  ....
  NJetParticles=j->GetNJetParticlesFF();
  NJetCells=j->NJetCellsFF();

  NParticles=j->GetNJetParticles();
  NCells=j->NJetCells();
  NSubJets=j->GetNSubJets();

  NIterations=j->GetNIterations();

  // fill FF related informations ...
  for (int i=0;i<NJetCells;i++)
    {
      ktJetCell *myCell=(ktJetCell*) j->GetJetCellFF(i);

      //DEBUG:
      //cout<<i<<" NFF: "<<myCell->NFFParticles()<<endl;
      // cout<<i<<" NP: "<<myCell->GetParticleList()->GetEntries()<<" NPID: "<<myCell->GetParticleListPID()->GetEntries()<<endl;
      // Add FF particle PID !!!

      for (int k=0;k<myCell->NParticlesPID();k++)
	{
	  JetParticles->AddLast((ktParticle*) (myCell->GetParticleListPID()->At(k)->Clone()));
	  //cout << "Adding ktParticle " << endl;
	} 
    }

  // jet infos neutral/charged energy (uncorrected)
  // DEBUG
  //cout<<"Fill jet neutral ..."<<endl;
  //cout<<j->NJetCells()<<endl;

  for (int i=0;i<(Int_t) j->NJetCells();i++)
    {
      ktJetCell *myCell=(ktJetCell*) j->GetJetCell(i);
      
      // DEBUG:
      //cout<<i<<" "<<myCell->GetParticleList()->GetEntries()<<" "<<myCell->GetParticleListPID()->GetEntries()<<endl;

      if (i==0)
	{
	  // fill seed info ...
	  mSeedEta=myCell->Eta();
	  mSeedPhi=myCell->Phi();
	  mSeedPt=myCell->Pt(); 
	}

      // some PID infos if filled ...
      for (int k=0;k<myCell->GetParticleListPID()->GetEntries();k++)
	{
	  ktPID *myPID=(ktPID*) myCell->GetParticleListPID()->At(k);
	  TLorentzVector *mVec=(TLorentzVector*) myCell->GetParticleList()->At(k);
	  if (myPID->GetCharge() != 0)
	    ChargedEnergy += mVec->Pt();
	  else
	    NeutralEnergy += mVec->Pt();

	  if (i==0)
	    {
	      // check if seed has electron mark as electron tagged !
	      if (myPID->GetPID() == 11)
		{
		  electronSeed=true;
		  // DEBUG:
		  //cout<<"Jet has electron seed ..."<<endl;
		  //cout<<mVec->Eta()<<" "<<mVec->Phi()<<" "<<mVec->Pt();
		}
	    }
	}
    }

  //DEBUG:
  //cout<<ChargedEnergy<<" "<<NeutralEnergy<<" "<<ChargedEnergy+NeutralEnergy<<endl;

  type=j->GetType();

}

ktMuJet::~ktMuJet()
{
  
  //cout<<"Jet part ..."<<endl;
  JetParticles->Delete();
  //cout<<"->Delete() done"<<endl;
  delete JetParticles;
 
  delete Jet;

  // DEBUG:
  //cout<<"Default ktMuJet destructor"<<endl;

}

Double_t ktMuJet::Phi() const
{
  Double_t mPhi;
  mPhi=Jet->Phi();if (mPhi<0) mPhi += (2*TMath::Pi());

  return mPhi;
}

Double_t ktMuJet::GetDistance(ktMuJet* mDjet)
{
  Double_t md=0;

  Double_t mdphi=Phi()-mDjet->Phi();
  Double_t mdeta=Eta()-mDjet->Eta();

  md=TMath::Sqrt(mdphi*mdphi+mdeta*mdeta);
  
  return md;
}

void ktMuJet::SetPtCellCorr(Double_t cellIn,Double_t cellOut)
{
  //DEBUG
  //cout<<JetPtCellCorr<<endl;
  //cout<<NCellsInEmcal<<" "<<cellIn<<endl;

  JetPtCellCorr=Pt()-NCellsInEmcal*cellIn-NCellsOutsideEmcal*cellOut;

  //DEBUG
  //cout<<JetPtCellCorr<<endl;
}

void ktMuJet::SetPtConeCorr(Double_t ptConeCorr)
{
  JetPtConeCorr=Pt()-ptConeCorr;
}

void ktMuJet::SetPtClusterCorr(Double_t ptClusterCorr)
{
  JetPtClusterCorr=Pt()-ptClusterCorr;
}


void ktMuJet::PrintJet() const
{

  //cout<<"Jet :"<<endl;
  //cout<<"-----"<<endl;
  cout<<endl;
  cout<<"# cells     = "<<GetNCells()<<endl;
  cout<<"# cells in EMCAL = "<<NCellsInEmcal<<endl;
  cout<<"# cells outside EMCAL = "<<NCellsOutsideEmcal<<endl;
  cout<<"# particles = "<<JetParticles->GetEntriesFast()<<endl;
  cout<<"# iterations = "<<GetNIterations()<<endl;
  cout<<"Type = "<<type<<endl;
  cout<<"Phi = "<<Phi()<<endl;
  cout<<"Eta = "<<Eta()<<endl;
  cout<<"Pt  = "<<Pt()<<endl;
  cout<<"Pt Cell Corr = "<<JetPtCellCorr<<endl;
  cout<<"Pt Cone Corr = "<<JetPtConeCorr<<endl;
  cout<<"Pt Cluster Corr = "<<JetPtClusterCorr<<endl;
  cout<<"E   = "<<E()<<endl;
  cout<<"Jet Scale = "<<JetAreaScale<<endl;
  cout<<"FF Scale = "<< FFAreaScale<<endl;
  cout<<endl;
}

void ktMuJet::PrintSeed() const
{
  cout<<endl;
  cout<<"Jet Seed:"<<endl;
  cout<<"pt = "<<mSeedPt<<"\t eta = "<<mSeedEta<<"\t phi = "<<mSeedPhi<<endl;
  cout<<"Jet:"<<endl;
  cout<<"pt = "<<Pt()<<"\t eta = "<<Eta()<<"\t phi = "<<Phi()<<endl;
  cout<<endl;
}
