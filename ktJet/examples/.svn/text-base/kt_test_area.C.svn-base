#include "../ktAreaUtil.h"

void kt_test_area(Double_t eta, Double_t phi,Double_t mRc)
{
  cout<<endl;
  cout<<"Area test program:"<<endl;
  cout<<"------------------"<<endl;
  cout<<endl;

  ktGrid *grid=new ktGrid();
  grid->SetGrid((int) 125/1,0,2*TMath::Pi(),(int) 40,-1.0,1.0); // grid2 (2xEMCAL)
  grid->SetEmcal(0.05,2*TMath::Pi()-0.05,-0.95,0.95);
  grid->SetSeed(4.6);
  grid->SetBkgPtCut(0.1);
  grid->SetCone(mRc);
  grid->SetBkgCone(0.4);
  grid->SetFFCone(mRc);
  grid->SetRMaxNN(1.0); // for kt-Jetfinder not used here ...
  grid->SetNBkgSeeds(150);
  grid->SetVerbose(0);
  grid->Init();
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  mhgrid=(TH1D*) (grid->hGrid)->Clone();

  TCanvas *c1=new TCanvas("c1","Canvas #1",900,300);
  mhgrid->Draw();
  
  TEllipse *mc=new TEllipse(phi,eta,mRc);
  //mc->SetLineColor(2);  
  //c1->cd();
  mc->Draw("");
  TLine *l=new TLine(phi-1.5*mRc,eta,phi+1.5*mRc,eta);
  l->SetLineColor(4);
  l->Draw();

  TLine *l2=new TLine(phi,eta-1.5*mRc,phi,eta+1.5*mRc);
  l2->SetLineColor(3);
  l2->Draw();
  
  Double_t mArea=TMath::Pi()*mRc*mRc;
  cout<<"Area = pi * "<<mRc<<"^2 = "<<mArea<<endl;
  cout<<endl;

  cout<<grid->GetNCellsInJet(mRc)<<endl;
  cout<<grid->GetNCellsInJet(eta,phi,mRc)<<endl;
  //cout<<grid->CalcBinsInJetArea(mRc)<<endl;
  //cout<<grid->CalcRcFromBins(grid->GetNCellsInJet(eta,phi,mRc))<<endl;
  cout<<grid->CalcScaleRatio(mRc,grid->CalcRcFromBins(grid->GetNCellsInJet(eta,phi,mRc)))<<endl;
  //cout<<(Int_t) grid->GetNCellsInJet(mRc)/(Double_t) grid->GetNCellsInJet(eta,phi,mRc)<<endl;
  cout<<endl;
  cout<<"Area in Acceptance from grid-cell counting = "<<mArea/grid->CalcScaleRatio(mRc,grid->CalcRcFromBins(grid->GetNCellsInJet(eta,phi,mRc)))<<endl;
  cout<<endl;

  //cout<<AreaCircle(1,1)<<endl;
  //cout<<TMath::Pi()/4.0<<endl;

  //cout<<IntersectionX(1,1.2)<<endl;
  cout<<AreaInAcc(mRc,eta,phi,1,2*TMath::Pi())<<endl;
  cout<<endl;


  delete grid;

}
