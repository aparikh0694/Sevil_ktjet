#include "ktPseudoJet.h"
#include <Riostream.h>

ClassImp(ktPseudoJet)

ktPseudoJet::ktPseudoJet()
{

  JetCells=new TObjArray(0);
  Jet=new TLorentzVector();
  Jet->SetPxPyPzE(0,0,0,0);
  type="Pseudo";

  // DEBUG:
  //cout<<"Default ktPseudoJet constructor "<<endl; 
}

ktPseudoJet::ktPseudoJet(ktJetCell* mJCell)
{
   JetCells=new TObjArray(0);
   type="Pseudo";
   JetCells->AddLast(mJCell);
   Finish();
   //cout<<"ktPseudoJet( ktJetCell* mJCell ) constructor "<<endl; 
}

ktPseudoJet::~ktPseudoJet()
{

  //if (JetCells->GetEntriesFast()>0)
  //delete JetCells;
  // DEBUG:
  //cout<<"Default ktPseudoJet destructor"<<endl;
}

void ktPseudoJet::Merge(ktPseudoJet *mMerge)
{
  // DEBUG:
  //cout<<" Merge PseudoJets ..."<<endl;
 
  for (int i=0;i<mMerge->NJetCells();i++)
    {
      JetCells->AddLast((ktJetCell*) (mMerge->GetJetCells()->At(i)));
      ((ktJetCell*) mMerge->GetJetCells()->At(i))->SetInJet(true);
    }
    //AddCell((ktJetCell*) (mMerge->GetJetCells()->At(i)));

  Finish();
}
