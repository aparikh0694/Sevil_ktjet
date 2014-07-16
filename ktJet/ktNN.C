// first test (Joern Putschke)

#include "ktNN.h"
#include <Riostream.h>

ClassImp(ktNN)

ktNN::ktNN()
{  
  // DEBUG:
  //cout<<"Default ktNN constructor "<<endl;
  NNCells=new TObjArray(0);
  
}

// check function def. !? Maybe better to split find and get NN !??
/*
ktJetCell* ktNN:GetNN()
{
  return 0;
}
*/

Bool_t ktNN::FindAll()
{
  
  // At some point to speed up, really use grid index !!!!!
  
  // DEBUG:
  //cout<<"FindAll ..."<<endl;
  //cout<<NNGrid->GetNphi()<<" "<<NNGrid->GetNeta()<<endl;

  for (int k=0;k<NNGrid->GetNeta();k++)
    {
      for (int l=0;l<NNGrid->GetNphi();l++)
	{
	  
	  //DEBUG:
	  //cout<<k<<" "<<l<<" "<<endl;

	  ktJetCell *myCell=NNGrid->GetCell(l,k);
	  
	  // check bkg pt cut && already in jet && lt Rc
	  if (myCell->NParticles()<1) continue;
	  if (myCell->InJet()) continue;
	 
	  myCell->SetIsNN(true);

	  ktPseudoJet *myPJet=new ktPseudoJet(myCell);
	  NNCells->AddLast(myPJet);
	  
	}
    }
  
  if (NNCells->GetEntriesFast()<2)
    return false;
  else
    return true;
}

// Find geomtrical NN cells on grid ! (based on jet or more independent index and
// adding in ktGrid !??? (hmmm, a bit more thinking needed ;-))
Bool_t ktNN::FindNN(ktJet *NNjet)
{
  return false;
}

// Just get all cells in RMaxNN cone to do brute force kt !!!
// (not really smart :-(, implement function above)
Bool_t ktNN::FindNNAroundSeed(ktJetCell *SeedCell)
{
  // DEBUG:
  //cout<<"Find NN in area around seed ..."<<endl;
  //cout<<NNJet->NJetCells()<<endl;
  
  //DEBUG:
  //ktJetCell *SeedCell=(ktJetCell*) NNJet->GetJetCells()->At(0);
  //cout<<NNCell->Eta()<<" "<<NNCell->Phi()<<endl;

  //Int_t etaBin=SeedCell->GetEtaBin();
  //Int_t phiBin=SeedCell->GetPhiBin();
  //DEBUG:
  //cout<<etaBin<<" "<<phiBin<<endl;

  // check grid bondaries (well, should really do it nicer :-( ) ...
 
  Int_t myEtaBinMin=SeedCell->GetEtaBin()-NNGrid->GetEtaNNConeBin();
  Int_t myEtaBinMax=SeedCell->GetEtaBin()+NNGrid->GetEtaNNConeBin();
  Int_t myPhiBinMin=SeedCell->GetPhiBin()-NNGrid->GetPhiNNConeBin();
  Int_t myPhiBinMax=SeedCell->GetPhiBin()+NNGrid->GetPhiNNConeBin();
  
  // check array border (not very nice :-(, write function :-))
  if (myEtaBinMin<0) myEtaBinMin=0;
  if (myPhiBinMin<0) myPhiBinMin=0;
  if (myEtaBinMax>NNGrid->GetNeta()) myEtaBinMax=NNGrid->GetNeta();
  if (myPhiBinMax>NNGrid->GetNphi()) myPhiBinMax=NNGrid->GetNphi();
  
  // DEBUG:
  //cout<<myEtaBinMin<<" - "<<myEtaBinMax<<endl;
  //cout<<myPhiBinMin<<" - "<<myPhiBinMax<<endl;
  Int_t mCount=0;

  //for (int k=(int) (mySeeds[i]->GetEtaBin()-GetEtaConeBin());k<(int) (mySeeds[i]->GetEtaBin()+GetEtaConeBin());k++)
  for (int k=myEtaBinMin;k<myEtaBinMax;k++)
    {
      //for (int l=(int) (mySeeds[i]->GetPhiBin()-GetPhiConeBin());l<(int) (mySeeds[i]->GetPhiBin()+GetPhiConeBin());l++)
      for (int l=myPhiBinMin;l<myPhiBinMax;l++)
	{
	  
	  ktJetCell *myCell=NNGrid->GetCell(l,k);
	  
	  // check bkg pt cut && already in jet && lt Rc
	  if (myCell->NParticles()<1) continue;
	  if (myCell->InJet()) continue;
	  
	  Double_t myR=TMath::Sqrt(TMath::Abs(SeedCell->Eta()-myCell->Eta())*TMath::Abs(SeedCell->Eta()-myCell->Eta())+TMath::Abs(SeedCell->Phi()-myCell->Phi())*TMath::Abs(SeedCell->Phi()-myCell->Phi()));
	  Double_t myPt=myCell->Pt();
		  
	  if (myR>NNGrid->GetRMaxNN()) continue;
	  if (myPt<NNGrid->GetBkgPtCut()) continue;

	  myCell->SetIsNN(true);

	  ktPseudoJet *myPJet=new ktPseudoJet(myCell);
	  NNCells->AddLast(myPJet);
	  //NNCells->AddLast(myCell);

	  // DEBUG: destructor test ...
	  //delete myPJet;
	  
	  mCount++;

	}
    }

  if (NNCells->GetEntriesFast()<2)
    return false;
  else
    return true;
}

void ktNN::Init(ktGrid *myNNGrid)
{
  NNGrid=myNNGrid;
  // DEBUG:
  //NNGrid->GridInfo();
}

ktNN::~ktNN()
{
  // DEBUG:
  //cout<<"Default ktNN destructor"<<endl;
  //NNCells->Delete(); // memory leak with PseudoJets !???? :-(
  delete NNCells; // check if cells are not deleted !!!

}
