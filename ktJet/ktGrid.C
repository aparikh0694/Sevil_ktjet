// first test (Joern Putschke)

#include "ktGrid.h"
#include "ktJetCellPair.h"
#include <Riostream.h>
#include "TF2.h"
#include "TRandom.h"

ClassImp(ktGrid)

ktGrid::ktGrid()
{

  verbose=true;
  SeedAccCut=false;
  
  phiMin=0;
  phiMax=2*TMath::Pi();
  etaMin=-1;
  etaMax=1;

  phiMinEmcal=0;
  phiMaxEmcal=2*TMath::Pi();
  etaMinEmcal=-1;
  etaMaxEmcal=1;

  EmcalSet=false;
  RcBset=false;

  Neta=100;
  Nphi=180;
  NCells=0;
  NSeeds=0;

  NMaxIterations=10;

  Seed=5.0;
  ClSeed=Seed;
  RMaxNN=1.5; // should be smaller by default if really NN search implemented (and not only area ...) :-(
  Rc=1.0; // default p+p case
  RcCl=0.2;
  RcClMax=1.0; // or at some point allow full EMCal for example ...
  RcB=Rc/2.0;
  Rff=0.7;

  BkgPtCut=0.0;
  minJetEnergy=0.0;
  BkgEnergyPerCell=0.0;
  BkgPtPerCell=0.0; 
  BkgPtPerCellAll=0.0; 
  BkgPtPerCellIn=0.0; 
  BkgPtPerCellAllIn=0.0; 
  BkgPtPerCellOut=0.0; 
  BkgPtPerCellAllOut=0.0;  
  BkgPtCone=0.0;
  BkgNCellCone=0.0;
  BkgPtCluster=0.0;
  BkgNCellCluster=0.0;

  BkgNCellIn=0.0;
  BkgNCellOut=0.0;
  BkgNJetCellIn=0.0;
  BkgNJetCellOut=0.0;
  BkgNCellInAll=0.0;
  BkgNCellOutAll=0.0;

  NBkgSeeds=100;

  MaxNCellsInJet=MaxNCellsInFF=0;

  Jets=new TObjArray(0);
  Jets->SetOwner(kTRUE);

  // only debug reasons !!!!
  ClusterJets=new TObjArray(0);
  ClusterJets->SetOwner(kTRUE);

  mySeedList=new TObjArray(0);
  Sharing=false;
  XiBkg=false;
  //hXiBkg=0;
  // DEBUG:
  //cout<<"Default ktGrid constructor "<<endl;
 
}

ktGrid::~ktGrid()
{

   // DEBUG:
  //cout<<"Default ktGrid destructor ..."<<endl;
  //cout<<"hXiBkg ..."<<endl;
  //hXiBkg->Dump();

  delete hXiBkg;
  delete hGrid;

  //cout<<"Delete myGrid[i] ..."<<endl;

  //cout<<Nphi<<endl;

  for (int i=0;i<(Nphi);i++) // ????
    {
      //cout<<i<<" "<<myGrid[i]<<endl;
      //if (myGrid[i])
      delete [] myGrid[i];
    }

  delete [] myGrid;

  //delete [] mySeeds;
  
  
  //cout<<"Delete Jet TObjArray  ... "<<Jets->GetEntries()<<endl;
  Jets->Delete();//Clear();
  delete Jets;

  // only debug reasons !!!!
  //cout<<"Delete ClusterJets ..."<<endl;
  ClusterJets->Delete();//Clear();
  delete ClusterJets;

   //cout<<"SeedList"<<endl;
   
  delete mySeedList;

}

// use then grid index (check in detail !!!)
void ktGrid::FindKtJetsBF()
{

 cout<<endl;
 cout<<" --> FindKtJets (brute force) ..."<<endl;
 cout<<" --> (Compare 1:1 to fastJet)"<<endl;
 //cout<<endl; 

 // copy from FindKtJets (should do it more clever ;-))
 
 ktNN *myNN=new ktNN();
 myNN->Init(this);
 
 Double_t myMinDij=1000.0;
 
 if (myNN->FindAll())
   {
     
     // DEBUG:
     //cout<<myNN->GetNNN()<<" NN in event found ..."<<endl;

     do {
       
       TObjArray *PairList=new TObjArray(0);
       
       // generate all pairs
       Int_t mNNN=myNN->GetNNN();
       
       for (int m=0;m<mNNN;m++)
	 {
	   for (int j=m+1;j<mNNN;j++)
	     {
	       
	       ktJetCellPair *mPair=new ktJetCellPair((ktPseudoJet*) myNN->GetNN()->At(m),(ktPseudoJet*) myNN->GetNN()->At(j));
	       PairList->AddLast(mPair);
	       
	     }
	 }
       
       // FindMinimum ..
       Int_t minIndex=0;
       myMinDij=1000.0;
       
       // write function :-)
       for (int k=0;k<(PairList->GetEntriesFast());k++)
	 {
	   ktJetCellPair *myPair=(ktJetCellPair*) PairList->At(k);
		      
	   if (myPair->dij<myMinDij)
	     {
	       minIndex=k;
	       myMinDij=myPair->dij;
	     }
	 }
		  
       // find minimum diB
       
       Int_t minIndexDiB=0;
       Double_t myMinDiB=1000.0;
       
       for (int l=0;l<mNNN;l++)
	 {
	   if (((ktPseudoJet*) myNN->GetNN()->At(l))->Kt2()<myMinDiB)
	     {
	       minIndexDiB=l;
	       myMinDiB=((ktPseudoJet*) myNN->GetNN()->At(l))->Kt2();
	     }
	 }
       
       ktJetCellPair *minCellPair=(ktJetCellPair*) PairList->At(minIndex);
       
       if (myMinDij<myMinDiB) 
	 {
	   
	   minCellPair->GetAJ()->Merge(minCellPair->GetBJ());
	   
	   // construct new NN list with merged Pair ...
	   myNN->GetNN()->Remove((ktPseudoJet*) minCellPair->GetBJ());
	   myNN->GetNN()->Compress();
	   
	   // clean up !
	   
	   PairList->Delete();
	   delete PairList; // check ich CellPairs get deleted, but Cells !???
	   
	 }
       else
	 {
	   
	   ktJet *jet=new ktJet((ktPseudoJet*) myNN->GetNN()->At(minIndexDiB));
	   // :-( stupid root with hamdling hist names !!!!
	   jet->SetType("KtBF");
	   jet->Finish(Nphi,phiMin,phiMax,Neta,etaMin,etaMax,1,1);
	   //jet->PrintJet();
	   // check minimun jet energy ...
	   Jets->AddLast(jet);

	   myNN->GetNN()->Remove((ktPseudoJet*) myNN->GetNN()->At(minIndexDiB));
	   myNN->GetNN()->Compress();
	   
	   // clean up !
	   PairList->Delete();
	   delete PairList; // check ich CellPairs get deleted, but Cells !???
	   
	 }
     } while (myNN->GetNNN()>1); // || myMinDij<25); // just a test
       
     ktJet *jet=new ktJet((ktPseudoJet*) myNN->GetNN()->At(0));
     jet->SetType("KtBF");
     jet->Finish(Nphi,phiMin,phiMax,Neta,etaMin,etaMax,1,1);
     // check minimun jet energy ...
     Jets->AddLast(jet);
     
   }
 else
   {
     cout<<"Warning: No NN in event found ..."<<endl;
   }
 
 delete myNN;
 
}


void ktGrid::FindKtJets()
{
  if (verbose)
    {
      cout<<endl;
      cout<<" --> FindKtJets (using seeds > "<<Seed<<" GeV) ..."<<endl;
      cout<<" --> Max. NN search radius = "<<RMaxNN<<endl;
    }
  //cout<<endl;
  //cout<<" --> (in this case just to define area ...)"<<endl;
  //cout<<" --> Take all cells around Max. NN and do kt brute force in that area ..."<<endl;
  //cout<<" --> More sophisticated NN search will be implemented soon (hopefully) ;-)"<<endl;
  //cout<<endl;

  // DEBUG:
  //cout<<"Loop over "<<NSeeds<<" Seeds ..."<<endl;
  //cout<<endl;

  // loop over seeds 
  //for (int i=0;i<NSeeds;i++)
  for (int i=0;i<mySeedList->GetEntriesFast();i++)
    {

      //DEBUG
      //PrintSeeds();
      //cout<<endl;
      //cout<<i<<" "<<mySeeds[i]->Pt()<<" "<<mySeeds[i]->Eta()<<" "<<mySeeds[i]->Phi()<<" "<<mySeeds[i]->GetEtaBin()<<" "<<mySeeds[i]->GetPhiBin()<<" "<<mySeeds[i]->InJet()<<endl;
      
      //if (!(mySeeds[i]->InJet()))
      if (!(((ktJetCell*) mySeedList->At(i))->InJet()))
	{
	  
	  ktNN *myNN=new ktNN();
	  myNN->Init(this);

	  //ktJet *kjet=new ktJet();
	  //jet->AddCell(mySeeds[i]);
	  //mySeeds[i]->SetInJet(true);

	  Double_t myMinDij=1000.0;

	  if (myNN->FindNNAroundSeed(((ktJetCell*) mySeedList->At(i))))//(mySeeds[i]))
	    {

	      do {
		//do {
		 
		// DEBUG:
		//cout<<myNN->GetNNN()<<" NN in area found ..."<<endl;

		TObjArray *PairList=new TObjArray(0);
		  
		  // generate all pairs
		  Int_t mNNN=myNN->GetNNN();
		  
		  // DEBUG
		  //cout<<"Generate all pairs ..."<<endl;
		  
		  for (int m=0;m<mNNN;m++)
		    {
		      for (int j=m+1;j<mNNN;j++)
			{
			  // DEBUG
			  //cout<<m<<" "<<j<<endl;
			  
			  ktJetCellPair *mPair=new ktJetCellPair((ktPseudoJet*) myNN->GetNN()->At(m),(ktPseudoJet*) myNN->GetNN()->At(j));
			  PairList->AddLast(mPair);
			  
			  //DEBUG:
			  //cout<<mPair->dij<<endl;
			  
			}
		    }
		  
		  //DEBUG
		  //cout<<"#Pairs = "<<PairList->GetEntriesFast()<<endl;
		  
		  // FindMinimum ..
		  Int_t minIndex=0;
		  //ktJetCellPair *minCellPair=0;
		  myMinDij=1000.0;
		  
		  // DEBUG
		  //cout<<"Find Minimum Pair dij ..."<<endl;
		  
		  // write function :-)
		  for (int k=0;k<(PairList->GetEntriesFast());k++)
		    {
		      ktJetCellPair *myPair=(ktJetCellPair*) PairList->At(k);
		      //ktJetCellPair *myPairB=(ktJetCellPair*) PairList->At(k);
		      
		      if (myPair->dij<myMinDij)
			{
			  minIndex=k;
			  myMinDij=myPair->dij;
			}
		    }
		  
		  // DEBUG:
		  //cout<<minIndex<<" "<<myMinDij<<endl;
		  //if (myMinDij>25) break;
		  
		  // find minimum diB
		  
		  Int_t minIndexDiB=0;
		  Double_t myMinDiB=1000.0;
		  
		  for (int l=0;l<mNNN;l++)
		    {
		      if (((ktPseudoJet*) myNN->GetNN()->At(l))->Kt2()<myMinDiB)
			{
			  minIndexDiB=l;
			  myMinDiB=((ktPseudoJet*) myNN->GetNN()->At(l))->Kt2();
			}
		    }
		  
		  // DEBUG:
		  //cout<<minIndexDiB<<" "<<myMinDiB<<endl;

		  ktJetCellPair *minCellPair=(ktJetCellPair*) PairList->At(minIndex);

		  if (myMinDij<myMinDiB) 
		    {
		      //DEBUG:
		      //minCellPair->GetAJ()->PrintJet();
		      //minCellPair->GetBJ()->PrintJet();
		      
		      minCellPair->GetAJ()->Merge(minCellPair->GetBJ());

		      //DEBUG:
		      //minCellPair->GetAJ()->PrintJet();
		      //minCellPair->GetBJ()->PrintJet();

		      // construct new NN list with merged Pair ...
		      myNN->GetNN()->Remove((ktPseudoJet*) minCellPair->GetBJ());
		      myNN->GetNN()->Compress();
		  
		      // clean up !

		      PairList->Delete();
		      delete PairList; // check ich CellPairs get deleted, but Cells !???
		      
		      //DEBUG
		      //cout<<"A "<<myNN->GetNNN()<<endl;
		    }
		  else
		    {

		      //DEBUG:
		      //cout<<endl;cout<<"\t ******** Ahhhh :-) ********"<<endl;cout<<endl;

		      ktJet *jet=new ktJet((ktPseudoJet*) myNN->GetNN()->At(minIndexDiB));
		      // :-( stupid root with handling hist names !!!!
		      jet->SetType("Kt");
		      jet->Finish(Nphi,phiMin,phiMax,Neta,etaMin,etaMax,((ktJetCell*) mySeedList->At(i))->GetPhiBin(),((ktJetCell*) mySeedList->At(i))->GetEtaBin());//mySeeds[i]->GetPhiBin(),mySeeds[i]->GetEtaBin());
		      //jet->PrintJet();
		      // check minimun jet energy ...
		      Jets->AddLast(jet);

		      myNN->GetNN()->Remove((ktPseudoJet*) myNN->GetNN()->At(minIndexDiB));
		      myNN->GetNN()->Compress();
		 
		      // clean up !
		      PairList->Delete();
		      delete PairList; // check ich CellPairs get deleted, but Cells !???

		      //DEBUG
		      //cout<<"B: "<<myNN->GetNNN()<<endl;		      
		    }

		  //DEBUG
		  //cout<<"C: "<<myNN->GetNNN()<<endl;

		  //} while (myNN->GetNNN()>1);
	      } while (myNN->GetNNN()>1); // || myMinDij<25); // just a test

	      // DEBUG:
	      //cout<<" *********************** "<<myNN->GetNN()->GetEntriesFast()<<endl;
	      //((ktPseudoJet*) myNN->GetNN()->At(0))->PrintJet();

	      ktJet *jet=new ktJet((ktPseudoJet*) myNN->GetNN()->At(0));
	      jet->SetType("Kt");
	      jet->Finish(Nphi,phiMin,phiMax,Neta,etaMin,etaMax,((ktJetCell*) mySeedList->At(i))->GetPhiBin(),((ktJetCell*) mySeedList->At(i))->GetEtaBin());//mySeeds[i]->GetPhiBin(),mySeeds[i]->GetEtaBin());
	      // check minimun jet energy ...
	      Jets->AddLast(jet);

	    }
	  else
	    {
	      cout<<"Warning: No NN in area found ..."<<endl;
	    }

	  delete myNN;
	  
	}
      //else
	//cout<<" Already in Jet !"<<endl; //DEBUG
      
    }
}


void ktGrid::AddFFParticles(ktJet* mJ)
{
  cout<<"WARNING: Add FF particles to Jet ... still be implemented !!!"<<endl;
}

void ktGrid::FindConeCluster()
{
  if (verbose)
    {
      cout<<endl;
      cout<<" --> FindConeCluster ..."<<endl;
      cout<<" --> Use RcCl = "<<GetRcCl()<<" and seeds > "<<Seed<<" GeV"<<endl;
      cout<<" --> for highest cluster. If EMCal set, only clusters full in EMCal acc."<<endl;
      cout<<" --> Cluster seeds > "<<GetClusterSeed()<<" GeV in are around highest cluster of "<<endl;
      cout<<" --> RcClMax = "<<GetRcClMax()<<endl;
      cout<<" --> Use Rff = "<<GetRff()<<" for FF"<<endl;
      cout<<endl;
    }

  // standard clustering can also be done offline looping over ktMuEvents !!!!
  // check adding TLorentzVector (Energy ...) !!!!

  // CHANGE !!!!!!
  // bad fix, maybe move to new class or modifyy GetSeedJet with variable Rc !!!!
  Double_t ConeRcTemp=Rc;
  Rc=GetRcCl();
  Int_t highestCluster=0;
  Int_t mUsedSubJets=0;

  TObjArray *clList=new TObjArray(0); // think about ownership !!!! for delete and main constructor
                                      // cleaner way would had been a new class though ;-)
  // done by ClusterJets in Grid !!!! only debug !!!
  //clList->SetOwner(true);

  ktJet *jCl=new ktJet();

  for (int i=0;i<mySeedList->GetEntriesFast();i++)
    {
      
      if ((((ktJetCell*) mySeedList->At(i))->InJet())) continue;
      
      ktJet *jet=new ktJet();
      
      GetSeedJet(jet,((ktJetCell*) mySeedList->At(i))->Eta(),((ktJetCell*) mySeedList->At(i))->Phi());
      if (jet->NJetCellsOutsideEmcal()<1) 
	{
	  jet->SetType("PseudoConeCluster");
	  highestCluster=i;
	  clList->AddLast(jet);
	  mUsedSubJets++;

	  // only debug !!!
	  //jet->SetInCluster();
	  ClusterJets->AddLast(jet);
	  break;
	}
      
      delete jet;
    }

  if (mUsedSubJets>0) {

  jCl->Add((ktJet*) clList->At(0));
  jCl->SetType("PseudoConeCluster");

  // DEBUG:
  //jCl->PrintJet();

  // do clustering ...
  // get seed list around highest cluster
  
  TObjArray *clSeedList=new TObjArray(0);
  GetSeeds(clSeedList,(ktJetCell*) mySeedList->At(highestCluster),GetRcClMax());

  // DEBUG:
  //cout<<"# Cluster seeds = "<<clSeedList->GetEntries()<<endl;

  // ge sub-jet around highest cluster
  // sort seeds highest - pt ! Other criteria !??? kt-like !???
  for (int i=0;i<clSeedList->GetEntries();i++)
    {
      //DEBIUG:
      //((ktJetCell*) clSeedList->At(i))->PrintCell();

      if ((((ktJetCell*) clSeedList->At(i))->InJet())) continue;

      ktJet *myClJet=new ktJet();
      GetSeedJet(myClJet,((ktJetCell*) clSeedList->At(i))->Eta(),((ktJetCell*) clSeedList->At(i))->Phi());
      myClJet->SetType("PseudoConeCluster");
      //DEBUG:
      //myClJet->PrintJet();
      
      // only debug
      ClusterJets->AddLast(myClJet);

      clList->AddLast(myClJet);
      //jCl->Add(myClJey);
    }
    
  delete clSeedList;
  
  // DBEUG:
  //cout<<"# sub-jets = "<<(clList->GetEntries()-1)<<endl;

  // sort or do just with respect to highest cluster based on seed !!!???
  // implement both !!!
  //clList->Sort();

  // do clustering ...
  
  if (clList->GetEntries()>1)
    {
      for (int i=1;i<clList->GetEntries();i++)
	{
	  ktJet *myCluJet=(ktJet*) clList->At(i);
	  //DEBUG:
	  //myCluJet->PrintJet();
	  
	  // add distance or distance/energy weighted criteria
	  // or take just everything in seed cone area ...
	  // DEBUG:
	  //cout<<((ktJet*) clList->At(0))->GetDistance(myCluJet)<<endl;

	  // obsolete because of max. seed cone cut ...
	  // can be modified and also add a energy weighting --> kt-like !???
	  if(((ktJet*) clList->At(0))->GetDistance(myCluJet)<1.0)
	    {
	      mUsedSubJets++;
	      // DEBUG
	      //cout<<"Subjet "<<i<<" | d = "<<((ktJet*) clList->At(0))->GetDistance(myCluJet)<<endl;
	      //myCluJet->PrintJet();
	      //myCluJet->SetInCluster();
	      jCl->Add(myCluJet);
	    }
	  else
	    {
	      // DEBUG ... (check memory ... but debug only)
	      // mainly for display ...
	      ClusterJets->RemoveAt(i);
	      ClusterJets->Compress();
	    }
	}
    }
  

  // include hyprid, only if #SubCones >1 else use standard cone !!!

  jCl->SetType("ConeCluster");
  jCl->SetNSubJets(mUsedSubJets);
  Jets->AddLast(jCl);
 
  // maybe not neccessary ... 
  //Jets->Sort();

  // bad fix !!!!
   Rc=ConeRcTemp;
   Int_t mNCellsInCluster=GetNCellsInCluster();
   BkgNCellCluster=(Double_t) mNCellsInCluster;
   BkgPtCluster=mNCellsInCluster/GetNCellsInJet(GetRcCl())*GetBkgPtCone()*RcB*RcB/(Rc*Rc);

  // DEBUG:
  //cout<<" ========= "<<mNCellsInCluster<<" "<<GetNCellsInJet(0.2)<<endl;
  //cout<<" ========= Bkg. pt in cluster = "<<BkgPtCluster<<endl;//*RcB*RcB/(ConeRcTemp*ConeRcTemp)<<endl;

  //delete jCl; // just for debug
  }
  else
   Rc=ConeRcTemp; 

  delete clList; 
}

// think about how and when to clean InCluster flag ....
Int_t ktGrid::GetNCellsInCluster() //TObjArray *mL)
{
  Int_t mN=0;

  //cout<<"========= DEBUG : "<<ClusterJets->GetEntries()<<" "<<GetRcCl()<<endl;

  for (int i=0;i<ClusterJets->GetEntries();i++)
    {

      ktJet *mJet=(ktJet*) ClusterJets->At(i);
      // go back to jet seed ...
      ktJetCell *mCell=(ktJetCell*) mJet->GetJetCell(0);
      // DEBUG
      //mCell->PrintCell();

      Double_t etaCenter=mCell->Eta();
      Double_t phiCenter=mCell->Phi();
      
      // DEBUG:
      //cout<<etaCenter<<" "<<phiCenter<<endl;
      
      Int_t mEtaBin=EtaBin(etaCenter);
      Int_t mPhiBin=PhiBin(phiCenter);
      
      //cout<<mEtaBin<<" "<<mPhiBin<<endl;
      
      Int_t myEtaBinMin=0;
      Int_t myEtaBinMax=0;
      Int_t myPhiBinMin=0;
      Int_t myPhiBinMax=0;
      
      myEtaBinMin=mEtaBin-GetEtaBinRc(GetRcCl());
      myEtaBinMax=mEtaBin+GetEtaBinRc(GetRcCl());
      myPhiBinMin=mPhiBin-GetPhiBinRc(GetRcCl());
      myPhiBinMax=mPhiBin+GetPhiBinRc(GetRcCl());
      
      // check array border (not very nice :-(, write function :-))
      if (myEtaBinMin<0) myEtaBinMin=0;
      if (myPhiBinMin<0) myPhiBinMin=0;
      if (myEtaBinMax>Neta) myEtaBinMax=Neta;
      if (myPhiBinMax>Nphi) myPhiBinMax=Nphi;
      
      //DEBUG
      //cout<<myEtaBinMin<<" "<<myEtaBinMax<<endl;
      //cout<<myPhiBinMin<<" "<<myPhiBinMax<<endl;
      
      for (int k=myEtaBinMin;k<myEtaBinMax;k++)
	{
	  for (int l=myPhiBinMin;l<myPhiBinMax;l++)
	    {
	      
	      ktJetCell *myCell=GetCell(l,k);
	      if (myCell->InCluster()) continue; //{cout<<"InCluster"<<endl;continue;}

	      Double_t myR=100.0;
	      Double_t mdEta=etaCenter-EtaFromBin(k);
	      Double_t mdPhi=phiCenter-PhiFromBin(l);;
	      
	      myR=TMath::Sqrt(mdEta*mdEta+mdPhi*mdPhi);
	      
	      if (myR>GetRcCl()) continue;
	      myCell->SetInCluster(true);

	      mN++;
	    }
	}
    }

  return mN;
}

Int_t ktGrid::GetNCellsInJet(Double_t mC)
{
  Int_t mN=0;

  Double_t etaCenter=(etaMinEmcal+etaMaxEmcal)/2.0;
  Double_t phiCenter=(phiMinEmcal+phiMaxEmcal)/2.0;

  // DEBUG:
  //cout<<etaCenter<<" "<<phiCenter<<endl;

  Int_t mEtaBin=EtaBin(etaCenter);
  Int_t mPhiBin=PhiBin(phiCenter);
  
  //cout<<mEtaBin<<" "<<mPhiBin<<endl;

  Int_t myEtaBinMin=0;
  Int_t myEtaBinMax=0;
  Int_t myPhiBinMin=0;
  Int_t myPhiBinMax=0;
  
  myEtaBinMin=mEtaBin-GetEtaBinRc(mC);
  myEtaBinMax=mEtaBin+GetEtaBinRc(mC);
  myPhiBinMin=mPhiBin-GetPhiBinRc(mC);
  myPhiBinMax=mPhiBin+GetPhiBinRc(mC);
  
  // check array border (not very nice :-(, write function :-))
  if (myEtaBinMin<0) myEtaBinMin=0;
  if (myPhiBinMin<0) myPhiBinMin=0;
  if (myEtaBinMax>Neta) myEtaBinMax=Neta;
  if (myPhiBinMax>Nphi) myPhiBinMax=Nphi;

  //DEBUG
  //cout<<myEtaBinMin<<" "<<myEtaBinMax<<endl;
  //cout<<myPhiBinMin<<" "<<myPhiBinMax<<endl;

  for (int k=myEtaBinMin;k<myEtaBinMax;k++)
    {
      for (int l=myPhiBinMin;l<myPhiBinMax;l++)
	{
	  
	  //ktJetCell *myCell=GetCell(l,k);
	  Double_t myR=100.0;
	  Double_t mdEta=etaCenter-EtaFromBin(k);
	  Double_t mdPhi=phiCenter-PhiFromBin(l);;

	  myR=TMath::Sqrt(mdEta*mdEta+mdPhi*mdPhi);
	  
	  if (myR>mC) continue;

	  mN++;
	}
    }

  return mN;
}

Int_t ktGrid::GetNCellsInJet(Double_t etaCenter, Double_t phiCenter,Double_t mC)
{
  Int_t mN=0;

  //Double_t etaCenter=(etaMinEmcal+etaMaxEmcal)/2.0;
  //Double_t phiCenter=(phiMinEmcal+phiMaxEmcal)/2.0;

  // DEBUG:
  //cout<<etaCenter<<" "<<phiCenter<<endl;

  Int_t mEtaBin=EtaBin(etaCenter);
  Int_t mPhiBin=PhiBin(phiCenter);
  
  //cout<<mEtaBin<<" "<<mPhiBin<<endl;

  Int_t myEtaBinMin=0;
  Int_t myEtaBinMax=0;
  Int_t myPhiBinMin=0;
  Int_t myPhiBinMax=0;
  
  myEtaBinMin=mEtaBin-GetEtaBinRc(mC);
  myEtaBinMax=mEtaBin+GetEtaBinRc(mC);
  myPhiBinMin=mPhiBin-GetPhiBinRc(mC);
  myPhiBinMax=mPhiBin+GetPhiBinRc(mC);
  
  // check array border (not very nice :-(, write function :-))
  if (myEtaBinMin<0) myEtaBinMin=0;
  if (myPhiBinMin<0) myPhiBinMin=0;
  if (myEtaBinMax>Neta) myEtaBinMax=Neta;
  if (myPhiBinMax>Nphi) myPhiBinMax=Nphi;

  //DEBUG
  //cout<<myEtaBinMin<<" "<<myEtaBinMax<<endl;
  //cout<<myPhiBinMin<<" "<<myPhiBinMax<<endl;

  for (int k=myEtaBinMin;k<myEtaBinMax;k++)
    {
      for (int l=myPhiBinMin;l<myPhiBinMax;l++)
	{
	  
	  //ktJetCell *myCell=GetCell(l,k);
	  Double_t myR=100.0;
	  Double_t mdEta=etaCenter-EtaFromBin(k);
	  Double_t mdPhi=phiCenter-PhiFromBin(l);;

	  myR=TMath::Sqrt(mdEta*mdEta+mdPhi*mdPhi);
	  
	  if (myR>mC) continue;

	  mN++;
	}
    }

  return mN;
}

void ktGrid::GetSeedJet(ktJet* seedJet,Double_t etaCenter, Double_t phiCenter)
{
  // DEBUG:
  //seedJet->PrintJet();
  seedJet->SetType("SeedJet");
  
  Int_t mEtaRandomBin,mPhiRandomBin;
  bool newSeed=false;

  mEtaRandomBin=EtaBin(etaCenter);
  mPhiRandomBin=PhiBin(phiCenter);
  
  ktJetCell *bkgSeedCell=GetCell(mPhiRandomBin,mEtaRandomBin); 
  //bkgSeedCell->SetInJet(true);

  //DEBUG:
  //cout<<bkgSeedCell->NParticles()<<endl;

  if (bkgSeedCell->NParticles()<1)
    {
      TLorentzVector *tempPart=new TLorentzVector();
      tempPart->SetPtEtaPhiM(0.00001,etaCenter,phiCenter,0);
      bkgSeedCell->AddParticle((TLorentzVector*) tempPart->Clone());
      bkgSeedCell->SetEtaBin(mEtaRandomBin);
      bkgSeedCell->SetPhiBin(mPhiRandomBin);
      // test ...
      seedJet->AddCell(bkgSeedCell,InEmcal(bkgSeedCell->Phi(),bkgSeedCell->Eta()),GetSharing());
      delete tempPart;
      newSeed=true;
    }
  else
    {
      seedJet->AddCell(bkgSeedCell,InEmcal(bkgSeedCell->Phi(),bkgSeedCell->Eta()),GetSharing());
      //bkgSeedCell->SetInJet(true);
    }
  
  Int_t myEtaBinMin=0;
  Int_t myEtaBinMax=0;
  Int_t myPhiBinMin=0;
  Int_t myPhiBinMax=0;
  
  myEtaBinMin=bkgSeedCell->GetEtaBin()-GetEtaConeBin();
  myEtaBinMax=bkgSeedCell->GetEtaBin()+GetEtaConeBin();
  myPhiBinMin=bkgSeedCell->GetPhiBin()-GetPhiConeBin();
  myPhiBinMax=bkgSeedCell->GetPhiBin()+GetPhiConeBin();
  
  // check array border (not very nice :-(, write function :-))
  if (myEtaBinMin<0) myEtaBinMin=0;
  if (myPhiBinMin<0) myPhiBinMin=0;
  if (myEtaBinMax>Neta) myEtaBinMax=Neta;
  if (myPhiBinMax>Nphi) myPhiBinMax=Nphi;

  //DEBUG
  //cout<<myEtaBinMin<<" "<<myEtaBinMax<<endl;
  //cout<<myPhiBinMin<<" "<<myPhiBinMax<<endl;

  for (int k=myEtaBinMin;k<myEtaBinMax;k++)
    {
      for (int l=myPhiBinMin;l<myPhiBinMax;l++)
	{
	  
	  ktJetCell *myCell=GetCell(l,k);
	  
	  // check bkg pt cut && already in jet && lt Rc
	  if (myCell->NParticles()<1) continue; // && myCell->NFFParticles()<1) continue;
	  if (!GetSharing() && myCell->InJet()) continue; //{cout<<"AHH "<<endl;continue;}
	  //if (myCell->InJet()) continue;
	  
	  Double_t myR=100.0;
	  //Double_t myPt=100.0;
	  
	  myR=TMath::Sqrt(TMath::Abs(bkgSeedCell->Eta()-myCell->Eta())*TMath::Abs(bkgSeedCell->Eta()-myCell->Eta())+TMath::Abs(bkgSeedCell->Phi()-myCell->Phi())*TMath::Abs(bkgSeedCell->Phi()-myCell->Phi()));
	  //}
	  
	  if (myR>Rc) continue;
	  // if (myPt<BkgPtCut) continue;
	  
	  seedJet->AddCell(myCell,InEmcal(myCell->Phi(),myCell->Eta()),GetSharing()); 
	}
    }
  
  if (newSeed)
    {
      bkgSeedCell->Clean(); 
      // remove seedCell from TObjArray to get correct number of jet cells !?? in principle not,
      // did not get added if Nparticles <1 ... ?????
      TObjArray *tempA;
      tempA=(TObjArray*) seedJet->GetJetCells();
      tempA->RemoveAt(0);
      tempA->Compress();
			  
    }
  
  // does this work with number of celles !!!????
  //seedJet->Finish();
  seedJet->SetType("PseudoCone");
  seedJet->Finish(Nphi,phiMin,phiMax,Neta,etaMin,etaMax,bkgSeedCell->GetPhiBin(),bkgSeedCell->GetEtaBin()); 
}

void ktGrid::BkgXi()
{

  Double_t xiPt=100.0;
  Int_t mN=0;
  Int_t NBkgXi=75; //default 50

  if (verbose)
    {
      cout<<endl;
      cout<<" --> Estimate background Xi-dist using Rff = "<<GetRff()<<endl;
      cout<<" --> Use "<<NBkgXi<<" random cones with RcBkg = "<<GetRff()<<endl;
      cout<<" --> Use canonical Jet pt of pt,jet = "<<xiPt<<" GeV"<<endl;
      cout<<" --> (rescaling to real pt,jet depending on bkg. and jet-finder method)"<<endl;
      cout<<endl;
    }

  // check how to get the same area/loss then with bkground and ....
  TF2 *f2= new TF2("f2","1+0*x*y",etaMin+Rff,etaMax-Rff,phiMin+Rff,phiMax-Rff);

  // old default ALICE EMCAL settings ...
  /*
  TF2 *f2;
  if(!EmcalSet)
    f2= new TF2("f2","1+0*x*y",etaMin+RcB,etaMax-RcB,phiMin+Rc,phiMax-RcB);
  else
    f2= new TF2("f2","1+0*x*y",etaMinEmcal+RcB,etaMaxEmcal-RcB,phiMin+RcB,phiMax-RcB);
  */

   while (mN<NBkgXi)
    {
      Double_t mEtaRandom,mPhiRandom;
      Int_t mEtaRandomBin,mPhiRandomBin;
      Int_t mNFF=0;

      TH1D *hXiTemp=new TH1D("hXiTemp","Xi dist. of Background per event (temp)",50,0,10);

      f2->GetRandom2(mEtaRandom,mPhiRandom);
      mEtaRandomBin=EtaBin(mEtaRandom);
      mPhiRandomBin=PhiBin(mPhiRandom);
      
      bool m_cellsInJet=false;
      //Double_t mBkgConePt=0;
      //Int_t mNBkgCells=0;

      //  check ... !!!!  
      Int_t myEtaBinMin=0;
      Int_t myEtaBinMax=0;
      Int_t myPhiBinMin=0;
      Int_t myPhiBinMax=0;

      myEtaBinMin=mEtaRandomBin-GetEtaConeBinFF();
      myEtaBinMax=mEtaRandomBin+GetEtaConeBinFF();
      myPhiBinMin=mPhiRandomBin-GetPhiConeBinFF();
      myPhiBinMax=mPhiRandomBin+GetPhiConeBinFF();
      
      // check array border (not very nice :-(, write function :-))
      if (myEtaBinMin<0) myEtaBinMin=0;
      if (myPhiBinMin<0) myPhiBinMin=0;
      if (myEtaBinMax>Neta) myEtaBinMax=Neta;
      if (myPhiBinMax>Nphi) myPhiBinMax=Nphi;
          
      for (int k=myEtaBinMin;k<myEtaBinMax;k++)
	{
	  for (int l=myPhiBinMin;l<myPhiBinMax;l++)
	    {
	    
	      ktJetCell *myCellBkg = GetCell(l,k);

	      if (myCellBkg->InJet()) {m_cellsInJet=true;continue;} // or check with IsFF() could be acc. problem at some point ...
	      
	      //if (myCellBkg->NParticles()<1) continue; // || myCell->Pt()<0.1) continue; // && myCell->NFFParticles()<1) continue;
	      if (myCellBkg->NFFParticles()<1) continue;
	      
	      Double_t myR=100.0;
	      
	      myR=TMath::Sqrt(TMath::Abs(mEtaRandom-myCellBkg->EtaFF())*TMath::Abs(mEtaRandom-myCellBkg->EtaFF())+TMath::Abs(mPhiRandom-myCellBkg->PhiFF())*TMath::Abs(mPhiRandom-myCellBkg->PhiFF()));
	      
	      if (myR>Rff) continue;
	      
	      for (int i=0;i<myCellBkg->NFFParticles();i++)
		{
		  TLorentzVector *lvec=(TLorentzVector*) myCellBkg->GetFFParticleList()->At(i);
		  if (!lvec || lvec->Pt()>0) // check why !!!!???? happens before !!???
		    {
		      Double_t mXi=TMath::Log((Double_t) xiPt/(Double_t) lvec->Pt());
		      mNFF++;
		      hXiTemp->Fill(mXi);
		      //DEBUG:
		      //cout<<mN<<" "<<i<<" "<<mXi<<" "<<(Double_t) lvec->Pt()<<endl;
		    }
		}
	      //mNBkgCells++;

	    }
	}
   
      if (!m_cellsInJet)
	{
	  //DEBUG:
	  //cout<<mN<<" "<<mNFF<<" "<<mEtaRandom<<" "<<mPhiRandom<<endl;
	  hXiBkg->Add(hXiTemp);
	  mN++;
	}
      
      delete hXiTemp;
      
    }
   
   // check for offline consistency !!!
   hXiBkg->Scale(1/((Double_t) NBkgXi)); //*(Double_t) hXiBkg->GetBinWidth(1)));

   delete f2;
}

void ktGrid::BkgXiScaled()
{
  // Scaled meaning randomly remove clusters to
  // mimik acceptance loss for real jet ...
  // CAUTION: Works only for leading jet so far !!!!
  //          Change that hXiBkg get stored for each jet !

  Double_t xiPt=100.0;
  Int_t mN=0;
  Int_t NBkgXi=75; //default 50

  if (verbose)
    {
      cout<<endl;
      cout<<" --> Estimate background Xi-dist using Rff = "<<GetRff()<<endl;
      cout<<" --> Remove randomly cells according to leading jet GetFFAreaScale()"<<endl;
      cout<<" --> Use "<<NBkgXi<<" random cones with RcBkg = "<<GetRff()<<endl;
      cout<<" --> Use canonical Jet pt of pt,jet = "<<xiPt<<" GeV"<<endl;
      cout<<" --> (rescaling to real pt,jet depending on bkg. and jet-finder method)"<<endl;
      cout<<endl;
    }

  // check how to get the same area/loss then with bkground and ....
  TF2 *f2= new TF2("f2","1+0*x*y",etaMin+Rff,etaMax-Rff,phiMin+Rff,phiMax-Rff);

  // old default ALICE EMCAL settings ...
  /*
  TF2 *f2;
  if(!EmcalSet)
    f2= new TF2("f2","1+0*x*y",etaMin+RcB,etaMax-RcB,phiMin+Rc,phiMax-RcB);
  else
    f2= new TF2("f2","1+0*x*y",etaMinEmcal+RcB,etaMaxEmcal-RcB,phiMin+RcB,phiMax-RcB);
  */

  // get leading jet !!!
  ktJet *lJet=(ktJet*) Jets->At(0); // check if really sorted at that stage !
  Double_t mAreaScale=lJet->GetFFAreaScale();
  
  // DEBUG:
  //cout<<" ***** "<<mAreaScale<<endl;

   while (mN<NBkgXi)
    {
      Double_t mEtaRandom,mPhiRandom;
      Int_t mEtaRandomBin,mPhiRandomBin;
      Int_t mNFF=0;

      TH1D *hXiTemp=new TH1D("hXiTemp","Xi dist. of Background per event (temp)",50,0,10);

      f2->GetRandom2(mEtaRandom,mPhiRandom);
      mEtaRandomBin=EtaBin(mEtaRandom);
      mPhiRandomBin=PhiBin(mPhiRandom);
      
      bool m_cellsInJet=false;
      //Double_t mBkgConePt=0;
      //Int_t mNBkgCells=0;

      //  check ... !!!!  
      Int_t myEtaBinMin=0;
      Int_t myEtaBinMax=0;
      Int_t myPhiBinMin=0;
      Int_t myPhiBinMax=0;

      myEtaBinMin=mEtaRandomBin-GetEtaConeBinFF();
      myEtaBinMax=mEtaRandomBin+GetEtaConeBinFF();
      myPhiBinMin=mPhiRandomBin-GetPhiConeBinFF();
      myPhiBinMax=mPhiRandomBin+GetPhiConeBinFF();
      
      // check array border (not very nice :-(, write function :-))
      if (myEtaBinMin<0) myEtaBinMin=0;
      if (myPhiBinMin<0) myPhiBinMin=0;
      if (myEtaBinMax>Neta) myEtaBinMax=Neta;
      if (myPhiBinMax>Nphi) myPhiBinMax=Nphi;
          
      for (int k=myEtaBinMin;k<myEtaBinMax;k++)
	{
	  for (int l=myPhiBinMin;l<myPhiBinMax;l++)
	    {
	    
	      ktJetCell *myCellBkg = GetCell(l,k);

	      if (myCellBkg->InJet()) {m_cellsInJet=true;continue;} // or check with IsFF() could be acc. problem at some point ...
	      
	      if (gRandom->Uniform(1)>mAreaScale) continue;

	      //if (myCellBkg->NParticles()<1) continue; // || myCell->Pt()<0.1) continue; // && myCell->NFFParticles()<1) continue;
	      if (myCellBkg->NFFParticles()<1) continue;
	      
	      Double_t myR=100.0;
	      
	      myR=TMath::Sqrt(TMath::Abs(mEtaRandom-myCellBkg->EtaFF())*TMath::Abs(mEtaRandom-myCellBkg->EtaFF())+TMath::Abs(mPhiRandom-myCellBkg->PhiFF())*TMath::Abs(mPhiRandom-myCellBkg->PhiFF()));
	      
	      if (myR>Rff) continue;
	      
	      for (int i=0;i<myCellBkg->NFFParticles();i++)
		{
		  TLorentzVector *lvec=(TLorentzVector*) myCellBkg->GetFFParticleList()->At(i);
		  if (!lvec || lvec->Pt()>0) // check why !!!!???? happens before !!???
		    {
		      Double_t mXi=TMath::Log((Double_t) xiPt/(Double_t) lvec->Pt());
		      mNFF++;
		      hXiTemp->Fill(mXi);
		      //DEBUG:
		      //cout<<mN<<" "<<i<<" "<<mXi<<" "<<(Double_t) lvec->Pt()<<endl;
		    }
		}
	      //mNBkgCells++;

	    }
	}
   
      if (!m_cellsInJet)
	{
	  //DEBUG:
	  //cout<<mN<<" "<<mNFF<<" "<<mEtaRandom<<" "<<mPhiRandom<<endl;
	  hXiBkg->Add(hXiTemp);
	  mN++;
	}
      
      delete hXiTemp;
      
    }
   
   // check for offline consistency !!!
   hXiBkg->Scale(1/((Double_t) NBkgXi)); //*(Double_t) hXiBkg->GetBinWidth(1)));

   delete f2;
}

void ktGrid::BkgRandomConesFast()
{
  //cout<<endl;
  Double_t bkgJetScale=RcB*RcB/(Rc*Rc);

  if (verbose)
    {
      cout<<" --> Estimate bkg. with "<<NBkgSeeds<<" random cones in Grid accep. - RcB "<<endl;//= "<<RcB<<endl;
      cout<<" --> Use RcB = "<<RcB<<" ==> Area(BkgJet)/Area(Jet) = "<<bkgJetScale<<endl;
      //cout<<endl;
    }

  // modified for taking EMCAL acceptance into account ....
  TF2 *f2;
  
  if(!EmcalSet)
    f2= new TF2("f2","1+0*x*y",etaMin+RcB,etaMax-RcB,phiMin+Rc,phiMax-RcB);
  else
    f2= new TF2("f2","1+0*x*y",etaMinEmcal+RcB,etaMaxEmcal-RcB,phiMinEmcal+RcB,phiMaxEmcal-RcB);
 
  Double_t mBkgPtConeSum=0;
  Int_t mNCell=0;
  Int_t mN=0;

  //DEBUG:
  //Int_t mTestNcells=((ktJet*) (GetJetList()->At(0)))->NJetCells();
  //cout<<mTestNcells<<endl;

  ktJetCell *myCellBkg;

  //for (int i=0;i<NBkgSeeds;i++)
  while (mN<NBkgSeeds) // Why does it change the jets ... create new seed and clean ???
    {
      Double_t mEtaRandom,mPhiRandom;
      Int_t mEtaRandomBin,mPhiRandomBin;

      f2->GetRandom2(mEtaRandom,mPhiRandom);
      mEtaRandomBin=EtaBin(mEtaRandom);
      mPhiRandomBin=PhiBin(mPhiRandom);

      // DEBUG:
      //cout<<i<<" "<<mEtaRandom<<" "<<mPhiRandom<<endl;
      //cout<<" "<<mEtaRandomBin<<" "<<mPhiRandomBin<<endl;
      //if (GetCell(mPhiRandomBin,mEtaRandomBin)->NParticles()<1 || 
      //GetCell(mPhiRandomBin,mEtaRandomBin)->InJet()) continue; 
      //if (GetCell(mPhiRandomBin,mEtaRandomBin)->InJet()) continue; 

      //bool newSeed=false;
      bool m_cellsInJet=false;
      Double_t mBkgConePt=0;
      Int_t mNBkgCells=0;

      //  check ... !!!!  
      Int_t myEtaBinMin=0;
      Int_t myEtaBinMax=0;
      Int_t myPhiBinMin=0;
      Int_t myPhiBinMax=0;

      myEtaBinMin=mEtaRandomBin-GetEtaConeBinBkg();
      myEtaBinMax=mEtaRandomBin+GetEtaConeBinBkg();
      myPhiBinMin=mPhiRandomBin-GetPhiConeBinBkg();
      myPhiBinMax=mPhiRandomBin+GetPhiConeBinBkg();
      
      // check array border (not very nice :-(, write function :-))
      if (myEtaBinMin<0) myEtaBinMin=0;
      if (myPhiBinMin<0) myPhiBinMin=0;
      if (myEtaBinMax>Neta) myEtaBinMax=Neta;
      if (myPhiBinMax>Nphi) myPhiBinMax=Nphi;
          
      for (int k=myEtaBinMin;k<myEtaBinMax;k++)
	{
	  for (int l=myPhiBinMin;l<myPhiBinMax;l++)
	    {
	    
	      myCellBkg = GetCell(l,k);
	      // DEBUG:
	      //cout<<l<<" "<<k<<endl;
	      //myCellBkg=GetCell(10,20); // if change to l,k # jet cells change !!?????
	      /*
	      if (((ktJet*) (GetJetList()->At(0)))->NJetCells()!=mTestNcells)
		{
		  cout<<mN<<" "<<l<<" "<<k<<endl;
		  myCell->PrintCell();
		}
	      */

	      if (myCellBkg->InJet()) {m_cellsInJet=true;continue;} 
	      if (myCellBkg->NParticles()<1) continue; // || myCell->Pt()<0.1) continue; // && myCell->NFFParticles()<1) continue;
	      
	      Double_t myR=100.0;
	      //Double_t myPt=100.0;
	      
	      myR=TMath::Sqrt(TMath::Abs(mEtaRandom-myCellBkg->Eta())*TMath::Abs(mEtaRandom-myCellBkg->Eta())+TMath::Abs(mPhiRandom-myCellBkg->Phi())*TMath::Abs(mPhiRandom-myCellBkg->Phi()));
	      
	      if (myR>RcB) continue;
	      //if (myPt<BkgPtCut) continue;
	      
	      mBkgConePt += myCellBkg->Pt();
	      mNBkgCells++;

	    }
	}
   
      if (!m_cellsInJet)
	{
	  mNCell += mNBkgCells;//bkgJet->NJetCells();
	  mBkgPtConeSum += mBkgConePt;//bkgJet->Pt();
	  mN++;
	}
      
    } //while (mN<NBkgSeeds);

  // DEBUG:
  //cout<<" # of bkg jets w/o cells in jet = "<<mN<<endl;

  BkgPtCone=mBkgPtConeSum/(Double_t) mN*1/bkgJetScale;
  BkgNCellCone=mNCell/(Double_t) mN*1/bkgJetScale;

  //DEBUG:
  //cout<<BkgPtCone<<" "<<BkgPtCone/(mNCell/(Double_t) mN)<<endl;
  
  delete f2;
  
}

void ktGrid::BkgRandomCones()
{

  Double_t bkgJetScale=RcB*RcB/(Rc*Rc);

  if (verbose)
    {
      //cout<<endl;
      cout<<" --> Estimate bkg. with "<<NBkgSeeds<<" random cones in Grid accep. - RcB "<<endl;//= "<<RcB<<endl;
      cout<<" --> Use RcB = "<<RcB<<" ==> Area(BkgJet)/Area(Jet) = "<<bkgJetScale<<endl;
      //cout<<endl;
    }

  // modified for taking EMCAL acceptance into account ....
  TF2 *f2;
  
  if(!EmcalSet)
    f2= new TF2("f2","1+0*x*y",etaMin+RcB,etaMax-RcB,phiMin+Rc,phiMax-RcB);
  else
    f2= new TF2("f2","1+0*x*y",etaMinEmcal+RcB,etaMaxEmcal-RcB,phiMinEmcal+RcB,phiMaxEmcal-RcB);

  // DEBUG:
  //f2->Print();

  //f2= new TF2("f2","1+0*x*y",etaMin+Rc,etaMax-Rc,phiMin+Rc,phiMax-Rc);
 
  Double_t mBkgPtConeSum=0;
  Int_t mNCell=0;
  Int_t mN=0;

  //for (int i=0;i<NBkgSeeds;i++)
  while (mN<NBkgSeeds) // Why does it change the jets ... create new seed and clean ... ????
    {
      Double_t mEtaRandom,mPhiRandom;
      Int_t mEtaRandomBin,mPhiRandomBin;

      f2->GetRandom2(mEtaRandom,mPhiRandom);
      mEtaRandomBin=EtaBin(mEtaRandom);
      mPhiRandomBin=PhiBin(mPhiRandom);
      // DEBUG:
      //cout<<i<<" "<<mEtaRandom<<" "<<mPhiRandom<<endl;
      //cout<<" "<<mEtaRandomBin<<" "<<mPhiRandomBin<<endl;

      //if (GetCell(mPhiRandomBin,mEtaRandomBin)->NParticles()<1 || 
      //	  GetCell(mPhiRandomBin,mEtaRandomBin)->InJet()) continue; 
      
      if (GetCell(mPhiRandomBin,mEtaRandomBin)->InJet()) continue; 

      ktJet *bkgJet=new ktJet();
      ktJetCell *bkgSeedCell=GetCell(mPhiRandomBin,mEtaRandomBin);

      //bkgSeedCell->SetInJet(true);

      //DEBUG:
      //cout<<bkgSeedCell->NParticles()<<endl;
      
      // funny change in results ... !??? 
      // CHECK !!!!

      bool newSeed=false;
      bool cellsInJet=false;

      if (bkgSeedCell->NParticles()<1) 
	{
	  TLorentzVector *tempPart=new TLorentzVector();
	  //tempPart->SetPtEtaPhiM(0.00001,mEtaRandom,mPhiRandom,0);
	  tempPart->SetPtEtaPhiM(0.000001,mEtaRandom,mPhiRandom,0);
	  bkgSeedCell->AddParticle((TLorentzVector*) tempPart->Clone());
	  bkgSeedCell->SetEtaBin(mEtaRandomBin);
	  bkgSeedCell->SetPhiBin(mPhiRandomBin);
	  delete tempPart;
	  newSeed=true;
	}
      
      bkgJet->AddCellBkg(bkgSeedCell);

      //DEBUG:
      //bkgSeedCell->PrintCell();

      //  check ... !!!!  
      Int_t myEtaBinMin=0;
      Int_t myEtaBinMax=0;
      Int_t myPhiBinMin=0;
      Int_t myPhiBinMax=0;

      myEtaBinMin=bkgSeedCell->GetEtaBin()-GetEtaConeBinBkg();
      myEtaBinMax=bkgSeedCell->GetEtaBin()+GetEtaConeBinBkg();
      myPhiBinMin=bkgSeedCell->GetPhiBin()-GetPhiConeBinBkg();
      myPhiBinMax=bkgSeedCell->GetPhiBin()+GetPhiConeBinBkg();
      
      // check array border (not very nice :-(, write function :-))
      if (myEtaBinMin<0) myEtaBinMin=0;
      if (myPhiBinMin<0) myPhiBinMin=0;
      if (myEtaBinMax>Neta) myEtaBinMax=Neta;
      if (myPhiBinMax>Nphi) myPhiBinMax=Nphi;
      
      for (int k=myEtaBinMin;k<myEtaBinMax;k++)
	{
	  for (int l=myPhiBinMin;l<myPhiBinMax;l++)
	    {
	      
	      ktJetCell *myCell=GetCell(l,k);
	      
	      // check bkg pt cut && already in jet && lt Rc
	      if (myCell->NParticles()<1) continue;// || myCell->Pt()<0.1) continue; // && myCell->NFFParticles()<1) continue;
	      if (myCell->InJet()) {cellsInJet=true;continue;} //{cout<<"AHH "<<endl;continue;}
	      
	      Double_t myR=100.0;
	      Double_t myPt=100.0;
	      
	      myR=TMath::Sqrt(TMath::Abs(bkgSeedCell->Eta()-myCell->Eta())*TMath::Abs(bkgSeedCell->Eta()-myCell->Eta())+TMath::Abs(bkgSeedCell->Phi()-myCell->Phi())*TMath::Abs(bkgSeedCell->Phi()-myCell->Phi()));
	      
	      if (myR>RcB) continue;
	      if (myPt<BkgPtCut) continue;
	      
	      bkgJet->AddCellBkg(myCell);
	      
	    }
	}
      
      //cout<<"Finish ..."<<endl;
      //bkgJet->SetType("RandomConeBkg");
      bkgJet->Finish();//Nphi,phiMin,phiMax,Neta,etaMin,etaMax,bkgSeedCell->GetPhiBin(),bkgSeedCell->GetEtaBin());
      //cout<<"After Finish ..."<<endl;

      //cout<<bkgJet->Pt()<<endl;

      //cellsInJet=false;

      if (newSeed)
	{
	  if (!cellsInJet)
	    mNCell += (bkgJet->NJetCells()-1);
	  // do correct, remove new created seedCell !!!!!
	  bkgSeedCell->Clean(); // check memory !!!
	  //delete bkgSeedCell;
	  // DEBUG:
	  //cout<<bkgSeedCell->NParticles()<<endl;
	}
      else
	{
	  if (!cellsInJet)
	    mNCell += bkgJet->NJetCells();
	}
      
      if (!cellsInJet)
	{
	  mBkgPtConeSum += bkgJet->Pt();
	  mN++;
	}
      
      delete bkgJet;
      //delete bkgSeedCell;

      //cout<<"delete bkgJet"<<endl;
      
    }

  // DEBUG:
  //cout<<" # of bkg jets w/o cells in jet = "<<mN<<endl;
  //cout<<mBkgPtConeSum/(Double_t) mN<<" "<<mNCell/(Double_t) mN<<endl;

  BkgPtCone=mBkgPtConeSum/(Double_t) mN*1/bkgJetScale;
  BkgNCellCone=mNCell/(Double_t) mN*1/bkgJetScale;

  //DEBUG:
  //cout<<BkgPtCone<<" "<<BkgPtCone/(mNCell/(Double_t) mN)<<endl;
  
  delete f2;
  
}

void ktGrid::FindConeJetsBkg()
{
  if (verbose)
    {
      cout<<endl;
      cout<<" --> FindConeJets to determine Bkg per cell ..."<<endl;
      cout<<" --> Use Rc = "<<GetRc()<<" and seeds > "<<Seed<<" GeV"<<endl;
      cout<<" --> Use Rff = "<<GetRff()<<" for FF"<<endl;
      cout<<endl;
    }

  for (int i=0;i<mySeedList->GetEntriesFast();i++)
    {

      //DEBUG
      //PrintSeeds();
      //cout<<endl;
      //cout<<i<<" "<<mySeeds[i]->Pt()<<" "<<mySeeds[i]->Eta()<<" "<<mySeeds[i]->Phi()<<" "<<mySeeds[i]->GetEtaBin()<<" "<<mySeeds[i]->GetPhiBin()<<" "<<mySeeds[i]->InJet()<<endl;
      
      if (!(((ktJetCell*) mySeedList->At(i))->InJet())) //mySeeds[i]->InJet()))
	{
	  ktJet *jet=new ktJet();
	  jet->AddCell(((ktJetCell*) mySeedList->At(i)),InEmcal(((ktJetCell*) mySeedList->At(i))->Phi(),((ktJetCell*) mySeedList->At(i))->Eta()));
	  jet->AddCellFF(((ktJetCell*) mySeedList->At(i)));
     
	  ((ktJetCell*) mySeedList->At(i))->SetInJet(true);
	  
	  // collect cells < Rc around seeds

	  // DEBUG:
	  //cout<<GetPhiConeBin()<<" "<<GetEtaConeBin()<<endl;
	  //cout<<mySeeds[i]->GetEtaBin()-GetEtaConeBin()<<" "<<mySeeds[i]->GetEtaBin()+GetEtaConeBin()<<endl;

	  //  check ... !!!!  
	  Int_t myEtaBinMin=0;
	  Int_t myEtaBinMax=0;
	  Int_t myPhiBinMin=0;
	  Int_t myPhiBinMax=0;

	  if (GetRc()>=GetRff())
	    {
	      myEtaBinMin=((ktJetCell*) mySeedList->At(i))->GetEtaBin()-GetEtaConeBin();
	      myEtaBinMax=((ktJetCell*) mySeedList->At(i))->GetEtaBin()+GetEtaConeBin();
	      myPhiBinMin=((ktJetCell*) mySeedList->At(i))->GetPhiBin()-GetPhiConeBin();
	      myPhiBinMax=((ktJetCell*) mySeedList->At(i))->GetPhiBin()+GetPhiConeBin();
	    }
	  else
	    {
	      myEtaBinMin=((ktJetCell*) mySeedList->At(i))->GetEtaBin()-GetEtaConeBinFF();
	      myEtaBinMax=((ktJetCell*) mySeedList->At(i))->GetEtaBin()+GetEtaConeBinFF();
	      myPhiBinMin=((ktJetCell*) mySeedList->At(i))->GetPhiBin()-GetPhiConeBinFF();
	      myPhiBinMax=((ktJetCell*) mySeedList->At(i))->GetPhiBin()+GetPhiConeBinFF();
	    }

	  // check array border (not very nice :-(, write function :-))
	  if (myEtaBinMin<0) myEtaBinMin=0;
	  if (myPhiBinMin<0) myPhiBinMin=0;
	  if (myEtaBinMax>Neta) myEtaBinMax=Neta;
	  if (myPhiBinMax>Nphi) myPhiBinMax=Nphi;

	  // DEBUG:
	  //cout<<myEtaBinMin<<" - "<<myEtaBinMax<<endl;
	  //cout<<myPhiBinMin<<" - "<<myPhiBinMax<<endl;

	  for (int k=myEtaBinMin;k<myEtaBinMax;k++)
	    {
	      for (int l=myPhiBinMin;l<myPhiBinMax;l++)
		{

		  ktJetCell *myCell=GetCell(l,k);
		  
		  // check bkg pt cut && already in jet && lt Rc
		  if (myCell->NParticles()<1 && myCell->NFFParticles()<1) continue;
		  if (myCell->InJet()) continue;
		  
		  Double_t myR=100.0;
		  Double_t myPt=100.0;

		  if (myCell->NParticles()<1)
		    {
		      myR=TMath::Sqrt(TMath::Abs(((ktJetCell*) mySeedList->At(i))->Eta()-myCell->EtaFF())*TMath::Abs(((ktJetCell*) mySeedList->At(i))->Eta()-myCell->EtaFF())+TMath::Abs(((ktJetCell*) mySeedList->At(i))->Phi()-myCell->PhiFF())*TMath::Abs(((ktJetCell*) mySeedList->At(i))->Phi()-myCell->PhiFF()));
		    }
		  else
		    {
		      myR=TMath::Sqrt(TMath::Abs(((ktJetCell*) mySeedList->At(i))->Eta()-myCell->Eta())*TMath::Abs(((ktJetCell*) mySeedList->At(i))->Eta()-myCell->Eta())+TMath::Abs(((ktJetCell*) mySeedList->At(i))->Phi()-myCell->Phi())*TMath::Abs(((ktJetCell*) mySeedList->At(i))->Phi()-myCell->Phi()));
		    }

		  myPt=myCell->Pt();

		  // add FF cell list ... (check !!!!)
		  // Why if else structure here !???

		  if(myR<=Rff && Rff>Rc)
		    {
		      //DEBUG:  
		      //cout<<myR<<endl;
		      if (myCell->NFFParticles()>0)
			jet->AddCellFF(myCell);
		    }
		  else if(myR<=Rc && Rc>=Rff)
		    {
		      
		       if (myCell->NFFParticles()>0)
			 jet->AddCellFF(myCell);
		    }
		  
		  if (myR>Rc) continue;
		  if (myPt<BkgPtCut) continue;

		  // DEBUG:
		  //cout<<myR<<" "<<myPt<<" "<<myCell->NParticles()<<endl;

		  jet->AddCell(myCell,InEmcal(myCell->Phi(),myCell->Eta()));

		}
	      
	    }
	  
	  // DEBUG:
	  //cout<<"#Cells in current Jet = "<<jet->NJetCells()<<endl;

	  // calculate jet properties eta,phi & energy number of cells+particles ...
	  // not nice :-( (friend declaration for classes or inherit !?)

	  // areas:
	  /*
	  cout<<" *** In ktGrid area = "<<GetNCellsInJet(((ktJetCell*) mySeedList->At(i))->Eta(),((ktJetCell*) mySeedList->At(i))->Phi(),Rc)<<endl;
	  //cout<<MaxNCellsInJet<<" "<<GetNCellsInJet(Rc)<<endl;
	  cout<<" *** Jet Scale = "<<GetNCellsInJet(((ktJetCell*) mySeedList->At(i))->Eta(),((ktJetCell*) mySeedList->At(i))->Phi(),Rc)/(Double_t) MaxNCellsInJet<<endl;
	  cout<<" *** FF Scale = "<< GetNCellsInJet(((ktJetCell*) mySeedList->At(i))->Eta(),((ktJetCell*) mySeedList->At(i))->Phi(),Rff)/(Double_t) MaxNCellsInFF<<endl;
	  */

	  jet->SetType("ConeBkg");
	  // set proper area scaling due to accepatance effects
	  // especially for FF
	  jet->SetJetAreaScale(GetNCellsInJet(((ktJetCell*) mySeedList->At(i))->Eta(),((ktJetCell*) mySeedList->At(i))->Phi(),Rc)/(Double_t) MaxNCellsInJet);
	  jet->SetFFAreaScale(GetNCellsInJet(((ktJetCell*) mySeedList->At(i))->Eta(),((ktJetCell*) mySeedList->At(i))->Phi(),Rff)/(Double_t) MaxNCellsInFF);
	  jet->Finish(Nphi,phiMin,phiMax,Neta,etaMin,etaMax,((ktJetCell*) mySeedList->At(i))->GetPhiBin(),((ktJetCell*) mySeedList->At(i))->GetEtaBin());

	  //DEBUG:
	  //jet->PrintJet();	  
	  //cout<<jet->GetJetAreaScale()<<" "<<jet->GetFFAreaScale()<<endl;

	  Jets->AddLast(jet);
	}
    }

  Jets->Sort();
}

// maybe use grid index to be faster ... (could be wise ;-))
// Just a very first approach, think about a bit more ;-) !!!
// simple just count cells in jet and subtract or, subtract in grid and run iterative ..
// (at some point try both with real bkg ;-))

void ktGrid::CalcBkg()
{
  if (verbose)
    cout<<" --> Calculate Bkg. per cell"<<endl;
   
  TH1F *meanBkg=new TH1F("meanBkg","",100,0,10);
  
  Int_t NBkgCells=0;
  Int_t NBkgCellsIn=0;
  Int_t NBkgCellsOut=0;
  Int_t NBkgCellsAll=0;
  Int_t NBkgCellsAllIn=0;
  Int_t NBkgCellsAllOut=0;
  Int_t NJetCells=0; // just a test
  Int_t NJetCellsIn=0;
  Int_t NJetCellsOut=0;
  // DEBUG:
  Int_t NCluster=0;

  for (int i=0;i<Nphi;i++)
    {
      for (int j=0;j<Neta;j++)
	{
	  
	  if (myGrid[i][j].NParticles()!=0 && myGrid[i][j].Pt()!=0)
	    {
	      if (!myGrid[i][j].InJet())
		{
		  BkgEnergyPerCell +=myGrid[i][j].E();
		  BkgPtPerCell += myGrid[i][j].Pt();
		  //meanBkg->Fill(myGrid[i][j].Pt());
		  if (InEmcal(myGrid[i][j].Phi(),myGrid[i][j].Eta()))
		  //if (myGrid[i][j].InEmcal())
		    {BkgPtPerCellIn += myGrid[i][j].Pt();NBkgCellsIn++;}
		  else
		    {BkgPtPerCellOut += myGrid[i][j].Pt();NBkgCellsOut++;}
		  NBkgCells++;
		}
	      else
		{
		  if (InEmcal(myGrid[i][j].Phi(),myGrid[i][j].Eta()))
		  //if (myGrid[i][j].InEmcal())
		    {
		      NJetCellsIn++;
		    }
		 else
		   NJetCellsOut++;
		  NJetCells++;
		}
	    }
	  
	  if (myGrid[i][j].InCluster())
	     NCluster++;

	  // Please check !!!!
	  
	  //if (myGrid[i][j].InEmcal())
	  if (InEmcal(PhiFromBin(i),EtaFromBin(j)))
	    NBkgCellsAllIn++;
	  else
	    NBkgCellsAllOut++;
	  

	  NBkgCellsAll++;
       }
    }

  // DEBUG:
  //cout<<"After grid cell loop ..."<<endl;

  // get average bkg per bkg cell
  BkgEnergyPerCell=BkgEnergyPerCell/(Double_t) NBkgCells;
  BkgPtPerCellAll=BkgPtPerCell/(Double_t) (NBkgCellsAll-NJetCells);
  BkgPtPerCell=BkgPtPerCell/(Double_t) NBkgCells;

  // DEBUG:
  //cout<<"Emcal ..."<<endl;

  if(EmcalSet)
    {
      BkgPtPerCellAllIn=BkgPtPerCellIn/(Double_t) (NBkgCellsAllIn-NJetCellsIn);
      BkgPtPerCellIn=BkgPtPerCellIn/(Double_t) NBkgCellsIn;
      BkgPtPerCellAllOut=BkgPtPerCellOut/(Double_t) (NBkgCellsAllOut-NJetCellsOut);
      BkgPtPerCellOut=BkgPtPerCellOut/(Double_t) NBkgCellsOut;

      BkgNCellIn=NBkgCellsIn;
      BkgNCellOut=NBkgCellsOut;
      BkgNJetCellIn=NJetCellsIn;
      BkgNJetCellOut=NJetCellsOut;
      BkgNCellInAll=NBkgCellsAllIn;
      BkgNCellOutAll=NBkgCellsAllOut;
      
      // DEBUG:
      /*
      cout<<NCluster<<endl;
      cout<<NBkgCellsIn<<" "<<" "<<NJetCellsIn<<" "<<NBkgCellsAllIn<<" "<<NBkgCellsAllOut<<" "<<NBkgCellsAll<<endl;
      cout<<NBkgCellsIn/(Double_t) (NBkgCellsAllIn-NJetCellsIn)<<endl;
      cout<<NBkgCellsOut/(Double_t) (NBkgCellsAllOut-NJetCellsOut)<<endl;
      */
    }
  else
    {
       BkgPtPerCellIn=0.0; 
       BkgPtPerCellAllIn=0.0; 
       BkgPtPerCellOut=0;//BkgPtPerCell;//0.0; 
       BkgPtPerCellAllOut=0;//BkgPtPerCellAll;//0.0;  

       // Modify correctly for case of no Emcal ...
       BkgNCellIn=0.0;
       BkgNCellOut=0.0;
       BkgNJetCellIn=0.0;
       BkgNJetCellOut=0.0;
       BkgNCellInAll=0.0;
       BkgNCellOutAll=0.0;
    }
  delete meanBkg;
}

// Update methods form CalcBkg() ....
void ktGrid::CleanGridAndCalcBkg()
{

  // set in jet flag to 0
  // rewrite with ktJet::CleanCellsInJet();
  
  //cout<<endl;
  //cout<<" *** Clean Grid and calculate Bkg. per cell (Dummy) ***"<<endl;
  if (verbose)
    cout<<" --> Clean Grid and calculate Bkg. per cell"<<endl;
  //cout<<" --> (Ratio method not implemented yet :-( )"<<endl;
  //cout<<endl;

  CalcBkg();

  for (int i=0;i<Jets->GetEntries();i++)
    {
      ktJet *tempJet=(ktJet*) (Jets->At(i));
      tempJet->CleanCellsInJet();
    }

}

void ktGrid::CleanGrid()
{

  cout<<" --> Clean Grid "<<endl;
  for (int i=0;i<Jets->GetEntries();i++)
    {
      ktJet *tempJet=(ktJet*) (Jets->At(i));
      tempJet->CleanCellsInJet();
    }
}


void ktGrid::CalcBkgNJetsRemove(Int_t nJR)
{

  if (verbose)
    cout<<" --> Calculate Bkg. per cell (remove "<<nJR<<" highest jets)"<<endl;

  
  if (!(Jets->IsSorted())) Jets->Sort();

  if (Jets->GetEntries()>nJR && nJR>=0)
    {
      for (int i=nJR;i<Jets->GetEntries();i++)
	{
	  ktJet *tempJet=(ktJet*) (Jets->At(i));
	  tempJet->CleanCellsInJet();
	}
    }
  

  CalcBkg();

}

void ktGrid::CalcBkgRemoveDijet()
{
  cout<<"WARNING: Still to be implemented ...."<<endl;
}

void ktGrid::FindConeJets()
{

  // Calculation of Ncells incorrect  if after an iteration the new seed (center) cell is empty !!!
  // Fix asap and implement also the addition of the the FF cells and particles !!!

  if (verbose)
    {
      cout<<endl;
      cout<<" --> FindConeJets (stable Cones: varation smaller than grid size) ..."<<endl;
      cout<<" --> Use Rc = "<<GetRc()<<" and seeds > "<<Seed<<" GeV"<<endl;
      cout<<" --> Use Rff = "<<GetRff()<<" for FF"<<endl;
      cout<<" --> # Maximum iterations = "<<NMaxIterations<<endl;
      cout<<endl;
    }
  
  for (int i=0;i<mySeedList->GetEntriesFast();i++)
    {
      
      if ((((ktJetCell*) mySeedList->At(i))->InJet())) continue;

      TLorentzVector *oldJet=new TLorentzVector();
      ktJet *jet=new ktJet();

      //
      // check if really the center or the center of the grid cell, compare
      // to bkg. jet finder !!!
      //
      GetSeedJet(jet,((ktJetCell*) mySeedList->At(i))->Eta(),((ktJetCell*) mySeedList->At(i))->Phi());
      //jet->PrintJet();
      jet->CleanCellsInJet();

      oldJet->SetXYZM(jet->GetJetVector()->X(),jet->GetJetVector()->Y(),jet->GetJetVector()->Z(),jet->GetJetVector()->M());
      
      for (int k=0;k<NMaxIterations;k++)
	{
	 
	  //DEBUG:
	  if (k==(NMaxIterations-1)) cout<<"WARNING: Max. Iterations reached !!!"<<endl;

	  ktJet *jetIt=new ktJet();	  
	  GetSeedJet(jetIt,jet->Eta(),jet->Phi());
	  
	  // DEBUG:
	  //if (i==0)
	  /*
	    {
	      cout<<"Iteration = "<<k<<endl;
	      //jetIt->PrintJet();
	      cout<<"Distance = "<<jetIt->GetDistance(jet)<<endl;
	      cout<<"Distance = "<<jetIt->GetDistance(oldJet)<<endl;
	      cout<<"Seed (eta,phi) = "<<((ktJetCell*) mySeedList->At(i))->Eta()<<" "<<((ktJetCell*) mySeedList->At(i))->Phi()<<endl;
	    }
	   */

	  // check split in eta/phi; calculate distance here !!!
	  if (jetIt->GetDistance(oldJet)<CelldEta() && jetIt->GetDistance(oldJet)<CelldPhi())
	    {
	      jetIt->SetCellsInJet();
	      jetIt->SetNIterations(k);
	      //jetIt->AddFFParticles();
	      jetIt->SetType("Cone");
	      Jets->AddLast(jetIt); 
	      break;
	    }
	  
	  oldJet->SetXYZM(jetIt->GetJetVector()->X(),jetIt->GetJetVector()->Y(),jetIt->GetJetVector()->Z(),jetIt->GetJetVector()->M());
	  jetIt->CleanCellsInJet();
	  delete jetIt;

	}

      delete jet;
      delete oldJet;
    }
  
  Jets->Sort();
}

void ktGrid::DoJetfinding(TString option)
{

  if (verbose)
    {
      cout<<endl;
      cout<<"Start Jetfinding on Grid ..."<<endl;
      //cout<<endl;
    }
  
  if (GetBkgPtCut()>0)
    {
      if (verbose)
	cout<<"Bkg. cut in pt = "<<GetBkgPtCut()<<endl;
      //cout<<endl;
    }

  if (option=="kt") // || option=="")
    {
      //cout<<" --> Use kt-like Jetfinder ..."<<endl;
      FindKtJets();
    }
  else if (option=="coneBkg" || option=="")
    {
      //cout<<" --> Use Cone Jetfinder to determine Bkg per cell in Pb+Pb ..."<<endl;
      //cout<<" --> (simple: no iterations)"<<endl;
      FindConeJetsBkg();
    }
  else if (option=="cone")
    {
      //cout<<" --> Use Cone Jetfinder ..."<<endl;
      FindConeJets();
    }
  else if (option=="cluster")
    {
      FindConeCluster();
    }
  else if (option=="ktBF")
    {
      //cout<<" --> Use kt-like Jetfinder ..."<<endl;
      //cout<<" --> Max. NN search radius = "<<RMaxNN<<endl;
      FindKtJetsBF();
    }
  else
    {
      cout<<" --> Not implemented yet, exit ..."<<endl;
      //break;
    }
}

Bool_t ktGrid::GetSeeds(Double_t mRdistance)
{
	// DEBUG
	/*
	 cout<<endl;
	 cout<<"Find seeds ..."<<endl;
	 cout<<endl;
	 */
	
	Double_t mySeedPt=0;
	Double_t mySeedEta=0;
	Double_t mySeedPhi=0;
	Int_t mySeedEtaBin=0;
	Int_t mySeedPhiBin=0;
	
	Bool_t erg=false;
	Bool_t my_double=false;
	
	for (int i=0;i<NCells;i++)
	  {
	    mySeedPt=myGrid[GridIndex[i].iPhi][GridIndex[i].iEta].Pt();
	    mySeedEta=myGrid[GridIndex[i].iPhi][GridIndex[i].iEta].Eta();
	    mySeedPhi=myGrid[GridIndex[i].iPhi][GridIndex[i].iEta].Phi();
	    mySeedEtaBin=GridIndex[i].iEta;
	    mySeedPhiBin=GridIndex[i].iPhi;
	    
	    //if (mySeedPt>Seed && TMath::Abs(mySeedEta)<(etaMaxEmcal-mRdistance) &&
	    //mySeedPhi>mRdistance && mySeedPhi<(2*TMath::Pi()-mRdistance))
	    if (mySeedPt>Seed)
	      {
		//DEBUG:
		//cout<<mySeedPt<<endl;
		
		if (TMath::Abs(mySeedEta)<(etaMaxEmcal-mRdistance) &&
		    mySeedPhi>mRdistance && mySeedPhi<(2*TMath::Pi()-mRdistance))
		  erg=true;

		// check (better include index in ktJetCell or check in index array !!!)
		for (int j=0;j<mySeedList->GetEntriesFast();j++)
		  {
		    //if (mySeedPt==((ktJetCell*) mySeedList->At(j))->Pt() && mySeedEta==((ktJetCell*) mySeedList->At(j))->Eta() && mySeedPhi==((ktJetCell*) mySeedList->At(j))->Phi())
		    if (mySeedEtaBin==(Int_t) ((ktJetCell*) mySeedList->At(j))->GetEtaBin() && mySeedPhiBin==(Int_t) ((ktJetCell*) mySeedList->At(j))->GetPhiBin())
		      {my_double=true; continue;}
		  }
		
		if (!my_double)
		  {
		    //mySeeds[NSeeds]=&myGrid[GridIndex[i].iPhi][GridIndex[i].iEta];
		    mySeedList->AddLast(&myGrid[GridIndex[i].iPhi][GridIndex[i].iEta]); 
		    NSeeds++;
		  }
	      }
	    my_double=false;
	  }
	
	mySeedList->Sort();
	
	if (mySeedList->GetEntries()>0 && erg)
	  return true;
	else
	  {
	    if (verbose)
	      cout<<" Warning: No seeds with pt,seed > "<<Seed<<" and eta,phi "<<mRdistance<<" away from acc. boundaries found !"<<endl;
	    return false;
	  }
}

Bool_t ktGrid::GetSeeds()
{
  // DEBUG
  /*
  cout<<endl;
  cout<<"Find seeds ..."<<endl;
  cout<<endl;
  */

  Double_t mySeedPt=0;
  Double_t mySeedEta=0;
  Double_t mySeedPhi=0;
  Int_t mySeedEtaBin=0;
  Int_t mySeedPhiBin=0;

  Bool_t my_double=false;

  for (int i=0;i<NCells;i++)
    {
      mySeedPt=myGrid[GridIndex[i].iPhi][GridIndex[i].iEta].Pt();
      mySeedEta=myGrid[GridIndex[i].iPhi][GridIndex[i].iEta].Eta();
      mySeedPhi=myGrid[GridIndex[i].iPhi][GridIndex[i].iEta].Phi();
      mySeedEtaBin=GridIndex[i].iEta;
      mySeedPhiBin=GridIndex[i].iPhi;

      if (mySeedPt>Seed)
	{
	  //DEBUG:
	  //cout<<mySeedPt<<endl;

	  // check (better include index in ktJetCell or check in index array !!!)
	  for (int j=0;j<mySeedList->GetEntriesFast();j++)
	    {
	      //if (mySeedPt==((ktJetCell*) mySeedList->At(j))->Pt() && mySeedEta==((ktJetCell*) mySeedList->At(j))->Eta() && mySeedPhi==((ktJetCell*) mySeedList->At(j))->Phi())
	      if (mySeedEtaBin==(Int_t) ((ktJetCell*) mySeedList->At(j))->GetEtaBin() && mySeedPhiBin==(Int_t) ((ktJetCell*) mySeedList->At(j))->GetPhiBin())
		{my_double=true; continue;}
	    }

	  if (!my_double)
	    {
	      //mySeeds[NSeeds]=&myGrid[GridIndex[i].iPhi][GridIndex[i].iEta];
	      mySeedList->AddLast(&myGrid[GridIndex[i].iPhi][GridIndex[i].iEta]); 
	      NSeeds++;
	    }
	}
      my_double=false;
    }

  mySeedList->Sort();
  
  if (mySeedList->GetEntries()>0)
	  return true;
  else
	{
	  if (verbose)
		  cout<<" Warning: No seeds with pt,seed > "<<Seed<<" found !"<<endl;
	  return false;
	}

}

void ktGrid::GetSeeds(TObjArray* clSeedList, Double_t etaCenter, Double_t phiCenter, Double_t mRcCl)
{
  // Get seeds around the highest cluster ...
  cout<<"WARNING: still to be implemented !!!"<<endl;
}

void ktGrid::GetSeeds(TObjArray* clSeedList,ktJetCell *bkgSeedCell, Double_t mRcCl)
{
  // Get seeds around the highest cluster ...
  
  Int_t myEtaBinMin=0;
  Int_t myEtaBinMax=0;
  Int_t myPhiBinMin=0;
  Int_t myPhiBinMax=0;
  
  myEtaBinMin=bkgSeedCell->GetEtaBin()-GetEtaConeBinCl();
  myEtaBinMax=bkgSeedCell->GetEtaBin()+GetEtaConeBinCl();
  myPhiBinMin=bkgSeedCell->GetPhiBin()-GetPhiConeBinCl();
  myPhiBinMax=bkgSeedCell->GetPhiBin()+GetPhiConeBinCl();
  
  // check array border (not very nice :-(, write function :-))
  if (myEtaBinMin<0) myEtaBinMin=0;
  if (myPhiBinMin<0) myPhiBinMin=0;
  if (myEtaBinMax>Neta) myEtaBinMax=Neta;
  if (myPhiBinMax>Nphi) myPhiBinMax=Nphi;

  for (int k=myEtaBinMin;k<myEtaBinMax;k++)
    {
      for (int l=myPhiBinMin;l<myPhiBinMax;l++)
	{
	  
	  ktJetCell *myCell=GetCell(l,k);
	  
	  // check bkg pt cut && already in jet && lt Rc
	  if (myCell->NParticles()<1) continue; // && myCell->NFFParticles()<1) continue;
	  if (!GetSharing() && myCell->InJet()) continue;
	  if (myCell->Pt()<GetClusterSeed()) continue;
	  // check if in EMcal ... (maybe generlazie later !!!!)
	  if (!InEmcal(myCell->Phi(),myCell->Eta())) continue;
	  
	  clSeedList->AddLast(myCell);
	}
    }

  clSeedList->Sort();
}

void ktGrid::Init()
{
  // Intitialize Grid ...
  //cout<<"Intitialize Grid ..."<<endl;
  //GridInfo();

  // Grid 2d histogram
  hGrid=new TH2F("EtaPhi","Grid (pt weighted)",Nphi,phiMin,phiMax,Neta,etaMin,etaMax);
  hGrid->SetDirectory(0); // to avoid root from owning 
  hGrid->SetMinimum(0.1);

  hXiBkg=new TH1D("hXiBkg","Xi dist. of Background per event",50,0,10);
  hXiBkg->SetDirectory(0); // to avoid root from owning 
 
  if (!RcBset)
    RcB=Rc/2.0;

  myGrid=new ktJetCell* [Nphi];
  for (int i=0;i<Nphi;i++)
    myGrid[i]=new ktJetCell[Neta];
  
  // Max numbers of Cells in Jet and FF cone for current settings
  MaxNCellsInJet=GetNCellsInJet(Rc);
  MaxNCellsInFF=GetNCellsInJet(Rff);

}


void ktGrid::Fill(TLorentzVector *mPart)
{

  Double_t phi,eta,pt;
  phi=mPart->Phi();if (phi<0) phi += (2*TMath::Pi());
  eta=mPart->PseudoRapidity();
  pt=mPart->Pt();

  // Fill grid & grid histogram (QA)
  // include QA track cuts, vertex-z ... (later :-))

  if (phi>phiMin && phi<phiMax && eta>etaMin && eta<etaMax && pt>BkgPtCut)
    {
      hGrid->Fill(phi,eta,pt);
      
      myGrid[PhiBin(phi)][EtaBin(eta)].AddParticle(mPart);
      // to be backward compatible ...
      myGrid[PhiBin(phi)][EtaBin(eta)].AddFFParticle(mPart);
      myGrid[PhiBin(phi)][EtaBin(eta)].SetPhiBin(PhiBin(phi));
      myGrid[PhiBin(phi)][EtaBin(eta)].SetEtaBin(EtaBin(eta));
      
      GridIndex[NCells].iPhi=PhiBin(phi);
      GridIndex[NCells].iEta=EtaBin(eta);
      
      NCells++;
    }
  else
    delete mPart;

}

void ktGrid::Fill(TLorentzVector *mPart, Bool_t isCharged)
{

  Double_t phi,eta,pt;
  phi=mPart->Phi();if (phi<0) phi += (2*TMath::Pi());
  eta=mPart->PseudoRapidity();
  pt=mPart->Pt();

  // Fill grid & grid histogram (QA)
  // include QA track cuts, vertex-z ... (later :-))

  if (phi>phiMin && phi<phiMax && eta>etaMin && eta<etaMax)
    {
      
       if (isCharged)
	 {
	   //cout<<"charged"<<endl;
	   //cout<<mPart->Pt()<<endl;
	   myGrid[PhiBin(phi)][EtaBin(eta)].AddFFParticle((TLorentzVector*) mPart->Clone());
	   //myGrid[PhiBin(phi)][EtaBin(eta)].AddFFParticle(mPart);
	   myGrid[PhiBin(phi)][EtaBin(eta)].SetPhiBin(PhiBin(phi));
	   myGrid[PhiBin(phi)][EtaBin(eta)].SetEtaBin(EtaBin(eta));
	 }

      if (pt>BkgPtCut)
	{
	  //cout<<"all"<<endl;
	  hGrid->Fill(phi,eta,pt);
	  
	  myGrid[PhiBin(phi)][EtaBin(eta)].AddParticle(mPart);
	  myGrid[PhiBin(phi)][EtaBin(eta)].SetPhiBin(PhiBin(phi));
	  myGrid[PhiBin(phi)][EtaBin(eta)].SetEtaBin(EtaBin(eta));
	  
	  // check if in Emcal accpetance ...

	  if (EmcalSet)
	    {
	      if (InEmcal(phi,eta))
		  myGrid[PhiBin(phi)][EtaBin(eta)].SetInEmcal(true);  
	    }

	  GridIndex[NCells].iPhi=PhiBin(phi);
	  GridIndex[NCells].iEta=EtaBin(eta);
	  
	  NCells++;
	  //cout<<NCells<<endl;
	}
      else
	delete mPart;
      
      //NCells++;
    }
  else
    delete mPart;

}

// new fill rountine maintaining backwards compability !!!
void ktGrid::Fill(TLorentzVector *mPart, ktPID *mPID)
{
  // dummmy ... ktPID not filled yet

  Double_t phi,eta,pt;
  phi=mPart->Phi();if (phi<0) phi += (2*TMath::Pi());
  eta=mPart->PseudoRapidity();
  pt=mPart->Pt();
  
  // Fill grid & grid histogram (QA)
  // include QA track cuts, vertex-z ... (later :-))
  
  if (phi>phiMin && phi<phiMax && eta>etaMin && eta<etaMax)
    {
      if (mPID->GetCharge() != 0) // && mPID->GetPID()!=4)
	{
	  if (mPID->GetPID()!=11)
	    {
	      myGrid[PhiBin(phi)][EtaBin(eta)].AddFFParticle((TLorentzVector*) mPart->Clone());
	      //mPID->PrintPID();
	      myGrid[PhiBin(phi)][EtaBin(eta)].AddFFParticlePID(mPID);
	      myGrid[PhiBin(phi)][EtaBin(eta)].SetPhiBin(PhiBin(phi));
	      myGrid[PhiBin(phi)][EtaBin(eta)].SetEtaBin(EtaBin(eta));
	    }
	    // else
	    //{cout<<" ******* electron not filled ...."<<endl;}
	}

      if (pt>BkgPtCut)
	{
	  hGrid->Fill(phi,eta,pt);
	  
	  //myGrid[PhiBin(phi)][EtaBin(eta)].AddParticle((TLorentzVector*)mPart->Clone());
	  myGrid[PhiBin(phi)][EtaBin(eta)].AddParticle(mPart);
	  myGrid[PhiBin(phi)][EtaBin(eta)].AddParticlePID(mPID);
	  myGrid[PhiBin(phi)][EtaBin(eta)].SetPhiBin(PhiBin(phi));
	  myGrid[PhiBin(phi)][EtaBin(eta)].SetEtaBin(EtaBin(eta));
	  
	  // check if in Emcal accpetance ...	  
	  if (EmcalSet)
	    {
	      if (InEmcal(phi,eta))
		  myGrid[PhiBin(phi)][EtaBin(eta)].SetInEmcal(true);  
	    }

	  GridIndex[NCells].iPhi=PhiBin(phi);
	  GridIndex[NCells].iEta=EtaBin(eta);
	  
	  NCells++;
	}
      else
	delete mPart;
    }
  else
    delete mPart;
}

void ktGrid::PrintIndex()
{
  for (int i=0;i<NCells-1;i++)
    {
      cout<<i<<" "<<GridIndex[i].iPhi<<" "<<GridIndex[i].iEta<<endl;
    }

}

void ktGrid::PrintSeeds()
{

  cout<<endl;
  cout<<"Seed list ..."<<endl;
  cout<<endl;
  /*
  for (int i=0;i<NSeeds;i++)
    {
      cout<<i<<" "<<mySeeds[i]->Pt()<<" "<<mySeeds[i]->Eta()<<" "<<mySeeds[i]->Phi()<<" "<<mySeeds[i]->GetEtaBin()<<" "<<mySeeds[i]->GetPhiBin()<<" "<<mySeeds[i]->InJet()<<endl;
    }

  //cout<<endl;
  //cout<<mySeedList->GetEntriesFast()<<endl;
  cout<<endl;
  */

   for (int i=0;i<mySeedList->GetEntriesFast();i++)
    {
      cout<<i<<" "<<((ktJetCell*) mySeedList->At(i))->Pt()<<" "<<((ktJetCell*) mySeedList->At(i))->Eta()<<" "<<((ktJetCell*) mySeedList->At(i))->Phi()<<" "<<((ktJetCell*) mySeedList->At(i))->InJet()<<" "<<((ktJetCell*) mySeedList->At(i))->GetEtaBin()<<" "<<((ktJetCell*) mySeedList->At(i))->GetPhiBin()<<endl;
    }
}

// maybe use grid index to be faster ...
void ktGrid::PrintGrid()
{

  for (int i=0;i<Nphi;i++)
    {
      for (int j=0;j<Neta;j++)
	{
	  
	  if (myGrid[i][j].NParticles()!=0)
	    {
	      cout<<myGrid[i][j].Eta()<<" "<<myGrid[i][j].Phi()<<" "<<myGrid[i][j].E()<<" "<<myGrid[i][j].InJet()<<" "<<myGrid[i][j].GetEtaBin()<<" "<<myGrid[i][j].GetPhiBin()<<endl;
	      TObjArray *mList=myGrid[i][j].GetParticleList();

	      for (int n=0;n<myGrid[i][j].NParticles();n++)
		{
		  TLorentzVector *cellPart=(TLorentzVector*) mList->At(n);

		  cout<<"   "<<n<<" "<<cellPart->PseudoRapidity()<<" "<<cellPart->E()<<endl;
		}
	    }
	}
    }

}

void ktGrid::SetGrid(Int_t m_Nphi,Double_t m_phiMin,Double_t m_phiMax,Int_t m_Neta,Double_t m_etaMin,Double_t m_etaMax)
{
  Nphi=m_Nphi;
  Neta=m_Neta;
  etaMin=m_etaMin;
  etaMax=m_etaMax;
  phiMin=m_phiMin;
  phiMax=m_phiMax;
  // predefine EMcal as grid 
  etaMinEmcal=m_etaMin;
  etaMaxEmcal=m_etaMax;
  phiMinEmcal=m_phiMin;
  phiMaxEmcal=m_phiMax;
  
}

void ktGrid::Rebin(Double_t mRebin)
{
  Nphi=(Int_t) (Nphi/mRebin);
  Neta=(Int_t) (Neta/mRebin);
}

void ktGrid::SetEmcal(Double_t m_phiMin,Double_t m_phiMax,Double_t m_etaMin,Double_t m_etaMax)
{
  // define EMcal acceptance
  etaMinEmcal=m_etaMin;
  etaMaxEmcal=m_etaMax;
  phiMinEmcal=m_phiMin;
  phiMaxEmcal=m_phiMax;

  EmcalSet=true;
}

Bool_t ktGrid::InEmcal(Double_t mPhi,Double_t mEta)
{

  Bool_t result=false;

  if (TMath::Abs(mEta)<etaMaxEmcal && (mPhi>phiMinEmcal && mPhi<phiMaxEmcal))
    result=true;
  else
    result=false;

  return result;
}

// nicer layout / table (at some point :-))
void ktGrid::PrintJets()
{
  cout<<"Found "<<Jets->GetEntriesFast()<<" Jets:"<<endl;
  cout<<"-------------"<<endl;
  cout<<endl;
  for (int i=0;i<Jets->GetEntriesFast();i++)
    {
      cout<<"Jet #"<<i+1<<":"<<endl;
      ktJet *outJet=(ktJet*) Jets->At(i);
      outJet->PrintJet();
    }
  cout<<endl;
}

// done :-)
void ktGrid::PrintJetsNice() //Float_t MinJetPtCut=0)
{
  //Jets->Sort();//Jets->GetEntriesFast()+1);

  cout<<" Found "<<Jets->GetEntriesFast()<<" Jets:"<<endl;
  cout<<" Settings : Rc = "<<GetRc()<<" "<<" Pt,bkg = "<<GetBkgPtCut()<<endl;
  cout<<" -------------"<<endl;
  cout<<endl;
  printf(" %5s %15s %15s %15s %15s %12s %10s %10s\n","jet #", "eta",
         "phi", "E", "pt", "n particles","n cells","Type");
  for (int i=0;i<Jets->GetEntriesFast();i++)
    {
      ktJet *outJet=(ktJet*) Jets->At(i);
      TString outType=outJet->GetType();
      const char* outType2=outJet->GetType();

      if (outType.Contains("KtBF"))
	{
	  //DEBUG:
	  //cout<<"AHH"<<endl;
	  if (outJet->Pt()>5.0)
	    printf(" %5u %15.8f %15.8f %15.8f %15.8f %8u %11u %15s\n",i+1,outJet->Eta(),outJet->Phi(),outJet->E(),outJet->Pt(),outJet->GetNJetParticles(),outJet->NJetCells(),outType2);
	}
      else
	printf(" %5u %15.8f %15.8f %15.8f %15.8f %8u %11u %15s\n",i+1,outJet->Eta(),outJet->Phi(),outJet->E(),outJet->Pt(),outJet->GetNJetParticles(),outJet->NJetCells(),outType2);
      
    }
  cout<<endl;
}

void ktGrid::PrintJetsNiceBkgSub() //Float_t MinJetPtCut=0)
{
  cout<<" Found "<<Jets->GetEntriesFast()<<" Jets:"<<endl;
  cout<<" Settings : Rc = "<<GetRc()<<" "<<" Pt,bkg = "<<GetBkgPtCut()<<endl;
  cout<<" -------------"<<endl;
  cout<<endl;
  printf(" %5s %15s %15s %15s %15s %15s %12s %10s %10s\n","jet #", "eta",
         "phi", "pt", "pt,corr", "pt,corr,cone","n particles","n cells","Type");
  for (int i=0;i<Jets->GetEntriesFast();i++)
    {
      ktJet *outJet=(ktJet*) Jets->At(i);
      TString outType=outJet->GetType();
      const char* outType2=outJet->GetType();

      if (outType.Contains("KtBF"))
	{
	  //DEBUG:
	  //cout<<"AHH"<<endl;
	  if (outJet->Pt()>5.0)
	    printf(" %5u %15.8f %15.8f %15.8f %15.8f %8u %11u %15s\n",i+1,outJet->Eta(),outJet->Phi(),outJet->E(),outJet->Pt(),outJet->GetNJetParticles(),outJet->NJetCells(),outType2);
	}
      else
	if (EmcalSet)
	  printf(" %5u %15.8f %15.8f %15.8f %15.8f %15.8f %8u %11u %15s\n",i+1,outJet->Eta(),outJet->Phi(),outJet->Pt(),outJet->Pt()-outJet->NJetCellsInEmcal()*GetBkgPtPerCellIn()-outJet->NJetCellsOutsideEmcal()*GetBkgPtPerCellOut(),outJet->Pt()-BkgPtCone,outJet->GetNJetParticles(),outJet->NJetCells(),outType2);
	else
	  printf(" %5u %15.8f %15.8f %15.8f %15.8f %15.8f %8u %11u %15s\n",i+1,outJet->Eta(),outJet->Phi(),outJet->Pt(),outJet->Pt()-outJet->NJetCells()*GetBkgPtPerCell(),outJet->Pt()-BkgPtCone,outJet->GetNJetParticles(),outJet->NJetCells(),outType2);
      
    }
  cout<<endl;
}

void ktGrid::GridInfo()
{

  cout<<endl;
  cout<<"Grid Info:"<<endl;
  cout<<"----------"<<endl;
  cout<<"Phi : "<<Nphi<<" "<<phiMin<<" "<<phiMax<<endl;
  cout<<"Eta : "<<Neta<<" "<<etaMin<<" "<<etaMax<<endl;
  //cout<<endl;
  cout<<" => dEta x dPhi : "<<CelldEta()<<" x "<<CelldPhi()<<endl;
  cout<<endl;
  if(EmcalSet)
    {
      cout<<"Emcal acceptance:"<<endl;
      cout<<"---------------------"<<endl;
      cout<<"Phi : "<<phiMinEmcal<<" "<<phiMaxEmcal<<endl;
      cout<<"Eta : "<<etaMinEmcal<<" "<<etaMaxEmcal<<endl;
      cout<<endl;
    }
   if (GetBkgPtCut()>0)
    {
      cout<<"Use Bkg. cut in pt (Pb+Pb) = "<<GetBkgPtCut()<<endl;
      cout<<endl;
    }

}
