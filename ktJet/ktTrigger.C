// first test (Joern Putschke)

#include "ktTrigger.h"
#include "ktJetCell.h"
#include <Riostream.h>

ClassImp(ktTrigger)

ktTrigger::ktTrigger() : ktGrid()
{

  NCellsX=0;
  NCellsY=0;
  NStepX=0;
  NStepY=0;
  minPatchPt=0;
  
  TriggerPatches=new TObjArray(0);
  TriggerPatches->SetOwner(kTRUE);

   // DEBUG:
  //cout<<"Default ktTrigger constructor "<<endl;
}

ktTrigger::ktTrigger(Int_t mNCellsX,Int_t mNCellsY,Int_t mNStepX, Int_t mNStepY, Double_t mMinPatchPt) : ktGrid()
{

  NCellsX=mNCellsX;
  NCellsY=mNCellsY;
  NStepX=mNStepX;
  NStepY=mNStepY;
  minPatchPt=mMinPatchPt;

  TriggerPatches=new TObjArray(0);
  TriggerPatches->SetOwner(kTRUE);

   // DEBUG:
  //cout<<"Default ktTrigger constructor "<<endl;
}


ktTrigger::~ktTrigger()
{

   // DEBUG:
  //cout<<"Default ktTrigger destructor ..."<<endl;

  TriggerPatches->Delete();
  delete TriggerPatches;
}


void ktTrigger::RunTrigger()
{
  cout<<endl;
  cout<<"Run sliding patch trigger ..."<<endl;
  cout<<endl;
  
  // DEBUG:
  cout<<"GridSize = "<<GetNeta()<<" x "<<GetNphi()<<endl;
  cout<<endl;

  for (int i=0+(NCellsX);i<GetNphi()-(NCellsX);i +=NStepX)
    {
      for (int j=0+(NCellsY);j<GetNeta()-(NCellsY);j +=NStepY)
	{
	  // DEBUG
	  //cout<<i<<" "<<hGrid->GetXaxis()->GetBinCenter(i+1)<<" "<<hGrid->GetYaxis()->GetBinCenter(j+1)<<" "<<GetPatchEnergy(i,j)<<endl;
	  
	  Double_t mPatchE=GetPatchEnergy(i,j);
	  if (mPatchE<0.1 || mPatchE<minPatchPt) continue;

	  // check eta,phi +-1 histogram vs. grid !!!
	  ktTriggerPatch *mPatch=new ktTriggerPatch(mPatchE,j,i,NCellsX,NCellsY);
	  mPatch->SetPatchEtaPhi(hGrid->GetYaxis()->GetBinCenter(j+1)+(NCellsY/2)*CelldEta(),hGrid->GetXaxis()->GetBinCenter(i+1)+(NCellsX/2)*CelldPhi());
	  TriggerPatches->AddLast(mPatch);
	  
	}
    }

}

Double_t ktTrigger::GetPatchEnergy(Int_t mStartX,Int_t mStartY)
{

  Double_t mPatchEnergy=0;

  for (int k=mStartX;k<(mStartX+NCellsX);k++)
    {
      for (int l=mStartY;l<(mStartY+NCellsY);l++)
	{
	  // DEBUG:
	  //cout<<l<<" "<<k<<endl;

	  ktJetCell *myCell=GetCell(k,l);

	  // DEBUG:
	  //cout<<"After cell ..."<<endl;
	  //cout<<myCell->NParticles()<<endl;

	  if (myCell->NParticles()<1) continue;

	  // DEBUG:
	  //cout<<" after cell access ..."<<endl;

	  mPatchEnergy += myCell->Pt();
	}
    }
  
  return mPatchEnergy;
}

Double_t ktTrigger::GetPatchEnergyCenter(Int_t mCenterX,Int_t mCenterY)
{

  // Dummy if needed ...

  Double_t mPatchEnergy=0;

  return mPatchEnergy;
}

ktTriggerPatch* ktTrigger::GetTriggerPatch()
{

  if (TriggerPatches->IsSorted())
    return (ktTriggerPatch*) TriggerPatches->At(0);
  else
    {
      TriggerPatches->Sort();
      return (ktTriggerPatch*) TriggerPatches->At(0);
    }

}

Double_t ktTrigger::GetTriggerPt()
{
  
  if (TriggerPatches->IsSorted())
    return ((ktTriggerPatch*) (TriggerPatches->At(0)))->Pt();
  else
    {
      TriggerPatches->Sort();
      return  ((ktTriggerPatch*) (TriggerPatches->At(0)))->Pt(); 
    }
}

Double_t ktTrigger::GetTriggerEta()
{
  
   if (TriggerPatches->IsSorted())
     return ((ktTriggerPatch*) (TriggerPatches->At(0)))->GetEta();
  else
    {
      TriggerPatches->Sort();
      return  ((ktTriggerPatch*) (TriggerPatches->At(0)))->GetEta();
    }
}

Double_t ktTrigger::GetTriggerPhi()
{
  
  if (TriggerPatches->IsSorted())
    return ((ktTriggerPatch*) (TriggerPatches->At(0)))->GetPhi();
  else
    {
      TriggerPatches->Sort();
      return ((ktTriggerPatch*) (TriggerPatches->At(0)))->GetPhi();
    }
}

void ktTrigger::PrintTriggerInfo()
{

  cout<<endl;
  cout<<"Trigger information:"<<endl;
  cout<<"---------------------"<<endl;
  cout<<"Patch size (cell x cell) = "<<NCellsX<<" x "<<NCellsY<<endl;
  cout<<"Patch size (dEta x dPhi) = "<<NCellsX*CelldEta()<<" x "<<NCellsY*CelldPhi()<<endl;
  cout<<"Step size (#cells) X = "<<NStepX<<endl;
  cout<<"Step size (#cells) Y = "<<NStepY<<endl;
  cout<<"Trigger patch pt treshold = "<<minPatchPt<<endl;
  cout<<endl;

}

void ktTrigger::PrintTriggerPatches()
{
  if (!TriggerPatches->IsSorted()) TriggerPatches->Sort();

  cout<<endl;
  cout<<TriggerPatches->GetEntries()<<" trigger patches found:"<<endl;
  cout<<endl;

  for(int i=0;i<TriggerPatches->GetEntries();i++)
    {
      cout<<i<<"\t "<<((ktTriggerPatch*) (TriggerPatches->At(i)))->GetEta()<<"\t "<<((ktTriggerPatch*) (TriggerPatches->At(i)))->GetPhi()<<"\t "<<((ktTriggerPatch*) (TriggerPatches->At(i)))->Pt()<<endl;
    }
  cout<<endl;
}
