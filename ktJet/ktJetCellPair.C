// first test (Joern Putschke)

#include "ktJetCellPair.h"
#include <Riostream.h>

ClassImp(ktJetCellPair)

ktJetCellPair::ktJetCellPair()
{

  // DEBUG:
  dij=999;
  //cout<<"Default ktJetCellPair constructor"<<endl;
 
}

ktJetCellPair::ktJetCellPair(ktJetCell *m_A,ktJetCell *m_B)
{

  dij=999;
  // DEBUG:
  A=m_A; B=m_B;
  
  Double_t dEtaAB=A->Eta()-B->Eta();
  Double_t dPhiAB=A->Phi()-B->Phi();

  Double_t RAB=dEtaAB*dEtaAB+dPhiAB*dPhiAB;

  // include RktMax or only default 1 ???!!! (different from NN Rmax !!!!)
  dij=TMath::Min((Double_t) A->Kt2(),(Double_t) B->Kt2())*RAB;

  //cout<<"Init ktJetCellPair constructor"<<endl;
 
}

ktJetCellPair::ktJetCellPair(ktPseudoJet *m_AJ,ktPseudoJet *m_BJ)
{

  dij=999;
  // DEBUG:
  AJ=m_AJ; BJ=m_BJ;
  
  Double_t dEtaAB=AJ->Eta()-BJ->Eta();
  Double_t dPhiAB=AJ->Phi()-BJ->Phi();

  Double_t RAB=dEtaAB*dEtaAB+dPhiAB*dPhiAB;

  // include RktMax or only default 1 ???!!! (different from NN Rmax !!!!)
  dij=TMath::Min((Double_t) AJ->Kt2(),(Double_t) BJ->Kt2())*RAB;

  //cout<<"Init ktJetCellPair constructor"<<endl;
 
}

ktJetCell* ktJetCellPair::Merge()
{

  A->SetInJet(true);
  B->SetInJet(true);

  return 0;
}

ktJetCellPair::~ktJetCellPair()
{

  // DEBUG:
  //cout<<"Default ktJetCellPair destructor"<<endl;

}
