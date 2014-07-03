

#include "ktAna.h"
#include <Riostream.h>
#include "TMath.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TF1.h"

//_________________________________________________________________________
Bool_t checkTrigId2(ktMuEvent *muev, TString TrigSel){

 if( TrigSel.Contains("All") || TrigSel.Contains("ALL")){
   return true;
 }
 for( int trig=0; trig<muev->GetTriggerInfo()->GetEntries(); trig++){
   ktTriggerInfo *trigInfo =  (ktTriggerInfo *)muev->GetTriggerInfo()->At(trig);

   if( (trigInfo->isJPL0() || trigInfo->isJPL2())  && TrigSel.Contains("JP"))
     return true;
   else if( (trigInfo->isHTL0() || trigInfo->isHTL2())  && TrigSel.Contains("HT"))
     return true;
 }

 return false;
}

//_________________________________________________________________________
Bool_t findTrigMatch(ktMuEvent *muev,ktMuFastJet *jet)
{
 
  if (!(checkTrigId2(muev,"JP"))) return false;
 
  int njets=0;
  Double_t phiJ, etaJ, phiT, etaT, phiDiff, etaDiff;
  
  if( muev->GetTriggerInfo()->GetEntries() == 0) 
    return false;

  phiJ=jet->Phi();
  etaJ=jet->Eta();

  for( int trig=0; trig<muev->GetTriggerInfo()->GetEntries(); trig++){
    ktTriggerInfo *trigInfo =  (ktTriggerInfo *)muev->GetTriggerInfo()->At(trig);
    
    if( trigInfo->isJPL0() || trigInfo->isJPL2()){
      etaT = trigInfo->GetEta();
      phiT = trigInfo->GetPhi();
      if( phiT < 0) phiT +=  2*TMath::Pi();
      etaDiff = etaT-etaJ;
      phiDiff = phiT-phiJ;
      
      if( TMath::Abs(etaDiff) < 0.5 &&
	  TMath::Abs(phiDiff) < TMath::Pi()/6){
	return true;
      }
      //cout << "Jet: " << phiJ << " " << etaJ << " Patch " << phiT << " " << etaT << " " << phiDiff << " " << etaDiff << endl;
    }
  }
  
 return false;
}

//_________________________________________________________________________
//_________________________________________________________________________

ClassImp(ktAna)

ktAna::ktAna()
{
  nEvents=0;
  nJets=0;
  evPt=evPtDi=-99.0;
  evArea=evAreaDi=-1.0;

  tevPt=tevPtDi=-99.0;
  tevArea=tevAreaDi=-1.0;

  InitHist=false;
  JetAlgo="kt";
  verbose=false;
  FJAna=true;
  BkgCorr=true;

  minPtCut=5.0;
  histPtMax=75;
  refMultCut=0;

  ffPtMin=30;
  ffPtMax=100;

  diJetPtMin=20;
  diJetPtMax=100;

  minFakePtCut=15.0;
  fakeMindPhi=1.2;
  fakeMaxdPhi=2.20;

  etaMax=0.6;etaMin=-0.6;
  phiMax=2*TMath::Pi();//-0.4;
  phiMin=0.0;//0.4;
  diPhiMin=2.5;
  diPhiMax=3.8;
  
  TriggerPtToRecoil=false;

   // DEBUG:
  //cout<<"Default ktAna constructor"<<endl;
}


ktAna::ktAna(Double_t eMin, Double_t eMax, Double_t pMin, Double_t pMax, Double_t ptMin)
{
  nEvents=0;
  nJets=0;
  evPt=evPtDi=-99.0;
  evArea=evAreaDi=-1.0;

  tevPt=tevPtDi=-99.0;
  tevArea=tevAreaDi=-1.0;

  InitHist=false;
  JetAlgo="kt";
  verbose=false;
  FJAna=true;
  BkgCorr=true;

  refMultCut=0;

  ffPtMin=30;
  ffPtMax=100;
  diJetPtMin=20;
  diJetPtMax=100;

  minPtCut=ptMin;
  histPtMax=75;

  etaMin=eMin; etaMax=eMax; phiMin=pMin; phiMax=pMax;

  diPhiMin=2.4;
  diPhiMax=3.8;

  // hFakeJetRatio=0;
  
}

ktAna::~ktAna()
{

  if (InitHist)
    {
      // delete histograms ...
      delete hE;delete hEDi;delete hPhiDi;
      delete hXi; delete hXiBkg; delete hXiDi; delete hXiBkgDi;
      delete hZ; delete hZBkg; delete hZDi; delete hZBkgDi;
      delete hEtaPhi;delete hEtaPhiDi;
      delete hSpec;delete hSpecDi;
      delete hAreaPt;delete hAreaPtDi;
      delete hXiFake; delete hXiBkgFake;
      //if (hEFake!=0)
      delete hEFake;

      delete thE;delete thEDi;delete thPhiDi;
      delete thXi; delete thXiBkg; delete thXiDi; delete thXiBkgDi;
      delete thZ; delete thZBkg; delete thZDi; delete thZBkgDi;
      delete thSpec;delete thSpecDi;
      delete thAreaPt;delete thAreaPtDi;
      delete thPtDiff;delete thPtDiffDi;

      delete hPtDiff;delete hPtDiffDi;
      delete hdEdi;
      delete hdPhiPt;

      //if (hFakeJetRatio!=0)
      //delete hFakeJetRatio;
      
      delete hEPtNoPt;

      delete hArray;
    }
  // DEBUG:
  //cout<<"Default ktAna destructor"<<endl;
}

Double_t ktAna::GetdPhi(Double_t mphi,Double_t vphi)
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

void ktAna::ScaleHist(TH1D *mH)
{
  mH->Scale(1.0/(mH->GetBinWidth(1)));
}

void ktAna::ScaleHist(TH1D *mH,Double_t mS)
{
  if (mS>0)
    mH->Scale(1.0/(mH->GetBinWidth(1)*mS));
}


void ktAna::LabelHist(TH1D *mH,TString xA,TString yA)
{
  mH->GetXaxis()->SetTitleSize(0.045);
  mH->GetYaxis()->SetTitleSize(0.045);
  mH->GetXaxis()->SetTitle(xA);
  mH->GetYaxis()->SetTitle(yA);
}

void ktAna::LabelHist(TH2D *mH,TString xA,TString yA)
{
  mH->GetXaxis()->SetTitleSize(0.045);
  mH->GetYaxis()->SetTitleSize(0.045);
  mH->GetXaxis()->SetTitle(xA);
  mH->GetYaxis()->SetTitle(yA);
}


void ktAna::XiShift(TH1D* m_shift,Float_t a)
{
  TH1D *erg=(TH1D*) m_shift->Clone();
  Double_t x,xnew;
  Int_t binnew;

  for (int i=0;i<m_shift->GetNbinsX();i++)
    {
      x=m_shift->GetBinCenter(i);
      xnew=x+TMath::Log(a);
      binnew=m_shift->GetXaxis()->FindBin(xnew);

      // DEBUG:
      //cout<<xnew<<" "<<x<<" "<<m_shift->GetBinContent(i)<<" "<<a<<endl;

      erg->SetBinContent(binnew,m_shift->GetBinContent(i));
      //erg->Fill(xnew,m_shift->GetBinContent(i));
    }

  for (int i=0;i<m_shift->GetNbinsX();i++)
    {
      m_shift->SetBinContent(i,erg->GetBinContent(i));
    }

  delete erg;

}

void ktAna::ZShift(TH1D* m_shift,Float_t a)
{
  // dummy ...
}

void ktAna::SetFiducial(Double_t eMin, Double_t eMax, Double_t pMin, Double_t pMax)
{
  etaMin=eMin; etaMax=eMax; phiMin=pMin; phiMax=pMax;
}

Bool_t ktAna::IsDiJet(Double_t mPhi)
{
  if (mPhi>diPhiMin && mPhi<diPhiMax)
    return true;
  else
    return false;
}

Bool_t ktAna::IsFakeJet(Double_t mPhi)
{
  Double_t minFake=fakeMindPhi;
  Double_t maxFake=fakeMaxdPhi;
  Double_t fakeScale=TMath::Abs(minFake-maxFake);
  Double_t maxFake2=TMath::Pi()*2-minFake;
  Double_t minFake2=maxFake2-fakeScale;

  if ((mPhi>minFake && mPhi<maxFake) || (mPhi>minFake2 && mPhi<maxFake2))
    return true;
  else
    return false;
}

Int_t ktAna::GetNjets(Double_t ptMin, Double_t ptMax)
{
  Int_t bMin=hE->FindBin(ptMin);
  Int_t bMax=hE->FindBin(ptMax);

  return (Int_t)(hE->Integral(bMin,bMax));
}

Int_t ktAna::GetNDiJets(Double_t ptMin, Double_t ptMax)
{
  Int_t bMin=hEDi->FindBin(ptMin);
  Int_t bMax=hEDi->FindBin(ptMax);

  return (Int_t)(hEDi->Integral(bMin,bMax));
}

Int_t ktAna::GetNFakeJets(Double_t ptMin, Double_t ptMax)
{
  Int_t bMin=hEFake->FindBin(ptMin);
  Int_t bMax=hEFake->FindBin(ptMax);

  return (Int_t)(hEFake->Integral(bMin,bMax));
}

Int_t ktAna::GetNjetsTrig(Double_t ptMin, Double_t ptMax)
{
  Int_t bMin=thE->FindBin(ptMin);
  Int_t bMax=thE->FindBin(ptMax);

  return (Int_t)(thE->Integral(bMin,bMax));
}

Int_t ktAna::GetNDiJetsTrig(Double_t ptMin, Double_t ptMax)
{
  Int_t bMin=thEDi->FindBin(ptMin);
  Int_t bMax=thEDi->FindBin(ptMax);

  return (Int_t)(thEDi->Integral(bMin,bMax));
}

void ktAna::SetDijetPhiCuts(Double_t pMin, Double_t pMax)
{
  diPhiMin=pMin;diPhiMax=pMax;
}

void ktAna::SetFFPtRange(Double_t ptMin, Double_t ptMax)
{
  ffPtMin=ptMin;ffPtMax=ptMax;
}

void ktAna::SetDiJetPtRange(Double_t ptMin, Double_t ptMax)
{
  diJetPtMin=ptMin;diJetPtMax=ptMax;
}

void ktAna::PrintAnalysisCuts()
{
  cout<<endl;
  cout<<"Used Analysis cuts:"<<endl;
  cout<<"-------------------"<<endl;
  //cout<<endl;
  //cout<<"# Events = "<<nEvents<<endl;
  cout<<"JetAlgo = "<<JetAlgo<<endl;
  cout<<"ptMin = "<<minPtCut<<endl;
  cout<<"ptMin for fake jets = "<<minFakePtCut<<endl;
  if (TriggerPtToRecoil)
    cout<<"Use Trigger jet pt,rec for recoil"<<endl;
  cout<<"Eta range = "<<etaMin<<" .. "<<etaMax<<endl;
  cout<<"Phi range = "<<phiMin<<" .. "<<phiMax<<endl;
  cout<<"DiJet Phi range = "<<diPhiMin<<" .. "<<diPhiMax<<endl;
  cout<<"Trigger Jet (FF) pt range = "<<ffPtMin<<" .. "<<ffPtMax<<endl;
  cout<<"and DiJet pt range "<<diJetPtMin<<" .. "<<diJetPtMax<<endl; 
  cout<<"Bkg corections = "<<BkgCorr<<endl;
  cout<<endl;
}

void ktAna::SetCorrectionFile(TString fName)
{
  // Dummy so far ...
  cout<<endl;
  cout<<"Efficiency & jet energy correction file = "<<fName;
  cout<<endl;
}

void ktAna:: DoSpectraCorrections()
{
  // Dummy so far ...
}

void ktAna::InitHistograms(TString m_JetAlgo)
{
  JetAlgo=m_JetAlgo;
  InitHistograms();
}

void ktAna::InitHistograms(TString m_JetAlgo,Double_t m_histPtMax)
{
  histPtMax= m_histPtMax;
  JetAlgo=m_JetAlgo;
  InitHistograms();
}


void ktAna::InitHistograms()
{
  InitHist=true;

  hArray=new TObjArray(0);

  Int_t nBins=(Int_t) histPtMax;

  hdEdi=new TH1D("hdEdi","p_{t,Jet} trig - p_{t,Jet} recoil "+JetAlgo,40,-20,20);
  LabelHist(hdEdi,"p_{t,Jet}^{rec.,trig}-p_{t,Jet}^{rec.,recoil}","# entries");
  hdEdi->SetDirectory(0);
  hArray->AddLast(hdEdi);

  hE=new TH1D("hE","p_{t,Jet} "+JetAlgo,nBins,0,histPtMax);
  LabelHist(hE,"p_{t,Jet}^{rec.}","# entries");
  hE->SetDirectory(0);
  hArray->AddLast(hE);

  hEDi=new TH1D("hEDi","p_{t,Jet} di-jet "+JetAlgo,nBins,0,histPtMax);
  LabelHist(hEDi,"p_{t,Jet}^{rec.}","# entries");
  hEDi->SetDirectory(0);
  hArray->AddLast(hEDi);

  thE=new TH1D("thE","p_{t,Jet} triggered "+JetAlgo,nBins,0,histPtMax);
  LabelHist(thE,"p_{t,Jet}^{rec.}","# entries");
  thE->SetDirectory(0);
  hArray->AddLast(thE);

  thEDi=new TH1D("thEDi","p_{t,Jet} di-jet triggered "+JetAlgo,nBins,0,histPtMax);
  LabelHist(thEDi,"p_{t,Jet}^{rec.}","# entries");
  thEDi->SetDirectory(0);
  hArray->AddLast(thEDi);

  hPtDiff=new TH1D("hPtDiff","",80,-20,20);
  LabelHist(hPtDiff,"p_{t,Jet}^{rec.}(Algo1)-p_{t,Jet}^{rec.}(Algo2)","# entries");
  hPtDiff->SetDirectory(0);
  hArray->AddLast(hPtDiff);

  hPtDiffDi=new TH1D("hPtDiffDi","",80,-20,20);
  LabelHist(hPtDiffDi,"p_{t,Jet}^{rec.}(Algo1)-p_{t,Jet}^{rec.}(Algo2)","# entries");
  hPtDiffDi->SetDirectory(0);
  hArray->AddLast(hPtDiffDi);

  thPtDiff=new TH1D("thPtDiff","",80,-20,20);
  LabelHist(thPtDiff,"p_{t,Jet}^{rec.}(Algo1)-p_{t,Jet}^{rec.}(Algo2)","# entries");
  thPtDiff->SetDirectory(0);
  hArray->AddLast(thPtDiff);

  thPtDiffDi=new TH1D("thPtDiffDi","",80,-20,20);
  LabelHist(thPtDiffDi,"p_{t,Jet}^{rec.}(Algo1)-p_{t,Jet}^{rec.}(Algo2)","# entries");
  thPtDiffDi->SetDirectory(0);
  hArray->AddLast(thPtDiffDi);

  hSpec=new TH1D("hSpec","p_{t,Jet} "+JetAlgo,nBins,0,histPtMax);
  LabelHist(hSpec,"p_{t,Jet}^{rec.}","dN/dp_{t,Jet}^{rec.}");
  hSpec->Sumw2();
  hSpec->SetDirectory(0);
  hArray->AddLast(hSpec);

  hSpecDi=new TH1D("hSpecDi","p_{t,Jet} di-jet"+JetAlgo,nBins,0,histPtMax);
  LabelHist(hSpecDi,"p_{t,Jet}^{rec.}","dN/dp_{t,Jet}^{rec.}");
  hSpecDi->Sumw2();
  hSpecDi->SetDirectory(0);
  hArray->AddLast(hSpecDi);

  thSpec=new TH1D("thSpec","p_{t,Jet} triggered "+JetAlgo,nBins,0,histPtMax);
  LabelHist(thSpec,"p_{t,Jet}^{rec.}","dN/dp_{t,Jet}^{rec.}");
  thSpec->Sumw2();
  thSpec->SetDirectory(0);
  hArray->AddLast(thSpec);

  thSpecDi=new TH1D("thSpecDi","p_{t,Jet} di-jet triggered "+JetAlgo,nBins,0,histPtMax);
  LabelHist(hSpecDi,"p_{t,Jet}^{rec.}","dN/dp_{t,Jet}^{rec.}");
  thSpecDi->Sumw2();
  thSpecDi->SetDirectory(0);
  hArray->AddLast(thSpecDi);

  hPhiDi=new TH1D("hPhiDi","#Delta#phi highest to second highest jet "+JetAlgo,180,0,2*TMath::Pi());
  LabelHist(hPhiDi,"#Delta#Phi=#phi_{leading}-#phi_{sec.}","# entries");
  hPhiDi->SetDirectory(0);
  hArray->AddLast(hPhiDi);

  hdPhiPt=new TH2D("hdPhiPt","#Delta#phi of triggered (also highest) to second highest jet vs p_{t}^{rec}(second)"+JetAlgo,180,0,2*TMath::Pi(),nBins,0,histPtMax);
  LabelHist(hdPhiPt,"#Delta#Phi=#phi_{highest&trigger}-#phi_{sec.}","p_{t}^{rec}(second)");
  hdPhiPt->SetDirectory(0);
  hArray->AddLast(hdPhiPt);

  thPhiDi=new TH1D("thPhiDi","#Delta#phi highest to second highest jet triggered "+JetAlgo,180,0,2*TMath::Pi());
  LabelHist(thPhiDi,"#Delta#Phi=#phi_{leading}-#phi_{sec.}","# entries");
  thPhiDi->SetDirectory(0);
  hArray->AddLast(thPhiDi);

  // leading and second leading 
  // Xi dist.
  hXi=new TH1D("hXi","#xi distribution Signal+Background "+JetAlgo,50,0,10);
  LabelHist(hXi,"#xi","1/N_{Jet}dN/d#xi");
  hXi->Sumw2();
  hXi->SetDirectory(0);
  hArray->AddLast(hXi);

  hXiBkg=new TH1D("hXiBkg","#xi distribution Background "+JetAlgo,50,0,10);
  LabelHist(hXiBkg,"#xi","1/N_{Jet}dN/d#xi");
  hXiBkg->Sumw2();
  hXiBkg->SetDirectory(0);
  hArray->AddLast(hXiBkg);

  hXiDi=new TH1D("hXiDi","#xi distribution di-jet Signal+Background "+JetAlgo,50,0,10);
  LabelHist(hXiDi,"#xi","1/N_{Jet}dN/d#xi");
  hXiDi->Sumw2();
  hXiDi->SetDirectory(0);
  hArray->AddLast(hXiDi);

  hXiBkgDi=new TH1D("hXiBkgDi","#xi distribution di-jet Background "+JetAlgo,50,0,10);
  LabelHist(hXiBkgDi,"#xi","1/N_{Jet}dN/d#xi");
  hXiBkgDi->Sumw2();
  hXiBkgDi->SetDirectory(0);
  hArray->AddLast(hXiBkgDi);

  // z dist..
  hZ=new TH1D("hZ","z distribution Signal+Background "+JetAlgo,10,0,1);
  LabelHist(hZ,"z","1/N_{Jet}dN/dz");
  hZ->Sumw2();
  hZ->SetDirectory(0);
  hArray->AddLast(hZ);

  hZBkg=new TH1D("hZBkg","z distribution Background "+JetAlgo,10,0,1);
  LabelHist(hZBkg,"z","1/N_{Jet}dN/dz");
  hZBkg->Sumw2();
  hZBkg->SetDirectory(0);
  hArray->AddLast(hZBkg);

  hZDi=new TH1D("hZDi","z distribution di-jet Signal+Background "+JetAlgo,10,0,1);
  LabelHist(hZDi,"z","1/N_{Jet}dN/dz");
  hZDi->Sumw2();
  hZDi->SetDirectory(0);
  hArray->AddLast(hZDi);

  hZBkgDi=new TH1D("hZBkgDi","z distribution di-jet Background "+JetAlgo,10,0,1);
  LabelHist(hZBkgDi,"z","1/N_{Jet}dN/dz");
  hZBkgDi->Sumw2();
  hZBkgDi->SetDirectory(0);
  hArray->AddLast(hZBkgDi);

  // trigger and recoil
  // Xi dist.
  thXi=new TH1D("thXi","#xi distribution Signal+Background triggered "+JetAlgo,50,0,10);
  LabelHist(thXi,"#xi","1/N_{Jet}dN/d#xi");
  thXi->Sumw2();
  thXi->SetDirectory(0);
  hArray->AddLast(thXi);

  thXiBkg=new TH1D("thXiBkg","#xi distribution Background triggered "+JetAlgo,50,0,10);
  LabelHist(thXiBkg,"#xi","1/N_{Jet}dN/d#xi");
  thXiBkg->Sumw2();
  thXiBkg->SetDirectory(0);
  hArray->AddLast(thXiBkg);

  thXiDi=new TH1D("thXiDi","#xi distribution di-jet Signal+Background triggered "+JetAlgo,50,0,10);
  LabelHist(thXiDi,"#xi","1/N_{Jet}dN/d#xi");
  thXiDi->Sumw2();
  thXiDi->SetDirectory(0);
  hArray->AddLast(thXiDi);

  thXiBkgDi=new TH1D("thXiBkgDi","#xi distribution di-jet Background triggered "+JetAlgo,50,0,10);
  LabelHist(thXiBkgDi,"#xi","1/N_{Jet}dN/d#xi");
  thXiBkgDi->Sumw2();
  thXiBkgDi->SetDirectory(0);
  hArray->AddLast(thXiBkgDi);

  // z dist..
  thZ=new TH1D("thZ","z distribution Signal+Background triggered "+JetAlgo,10,0,1);
  LabelHist(thZ,"z","1/N_{Jet}dN/dz");
  thZ->Sumw2();
  thZ->SetDirectory(0);
  hArray->AddLast(thZ);

  thZBkg=new TH1D("thZBkg","z distribution Background triggered "+JetAlgo,10,0,1);
  LabelHist(thZBkg,"z","1/N_{Jet}dN/dz");
  thZBkg->Sumw2();
  thZBkg->SetDirectory(0);
  hArray->AddLast(thZBkg);

  thZDi=new TH1D("thZDi","z distribution di-jet Signal+Background triggered "+JetAlgo,10,0,1);
  LabelHist(thZDi,"z","1/N_{Jet}dN/dz");
  thZDi->Sumw2();
  thZDi->SetDirectory(0);
  hArray->AddLast(thZDi);

  thZBkgDi=new TH1D("thZBkgDi","z distribution di-jet Background triggered "+JetAlgo,10,0,1);
  LabelHist(thZBkgDi,"z","1/N_{Jet}dN/dz");
  thZBkgDi->Sumw2();
  thZBkgDi->SetDirectory(0);
  hArray->AddLast(thZBkgDi);

  // misc. 
  hEtaPhi=new TH2D("hEtaPhi","Eta vs. Phi leading jet "+JetAlgo,50,-1,1,180,0,2*TMath::Pi());
  LabelHist(hEtaPhi,"#eta","#phi");
  hEtaPhi->SetDirectory(0);
  hArray->AddLast(hEtaPhi);

  hEtaPhiDi=new TH2D("hEtaPhiDi","Eta vs. Phi di-jet "+JetAlgo,50,-1,1,180,0,2*TMath::Pi());
  LabelHist(hEtaPhiDi,"#eta(dijet)","#phi(dijet)");
  hEtaPhiDi->SetDirectory(0);
  hArray->AddLast(hEtaPhiDi);

  hAreaPt=new TH2D("hAreaPt","Area vs.  p_{t,Jet}^{rec.} "+JetAlgo,nBins,0,histPtMax,40,0,2);
  LabelHist(hAreaPt,"p_{t,Jet}^{rec.}","area");
  hAreaPt->SetDirectory(0);
  hArray->AddLast(hAreaPt);

  hAreaPtDi=new TH2D("hAreaPtDi","Area vs.  p_{t,Jet}^{rec.} di-jet "+JetAlgo,nBins,0,histPtMax,40,0,2);
  LabelHist(hAreaPtDi,"p_{t,Jet}^{rec.}","area");
  hAreaPtDi->SetDirectory(0);
  hArray->AddLast(hAreaPtDi);

  thAreaPt=new TH2D("thAreaPt","Area vs.  p_{t,Jet}^{rec.} triggered "+JetAlgo,nBins,0,histPtMax,40,0,2);
  LabelHist(thAreaPt,"p_{t,Jet}^{rec.}","area");
  thAreaPt->SetDirectory(0);
  hArray->AddLast(thAreaPt);

  thAreaPtDi=new TH2D("thAreaPtDi","Area vs.  p_{t,Jet}^{rec.} di-jet triggered "+JetAlgo,nBins,0,histPtMax,40,0,2);
  LabelHist(thAreaPtDi,"p_{t,Jet}^{rec.}","area");
  thAreaPtDi->SetDirectory(0);
  hArray->AddLast(thAreaPtDi);

  // fake jets ...
  hXiFake=new TH1D("hXiFake","#xi distribution Signal+Background Fake "+JetAlgo,50,0,10);
  LabelHist(hXiFake,"#xi","1/N_{Jet}dN/d#xi");
  hXiFake->Sumw2();
  hXiFake->SetDirectory(0);
  hArray->AddLast(hXiFake);

  hXiBkgFake=new TH1D("hXiBkgFake","#xi distribution Background Fake "+JetAlgo,50,0,10);
  LabelHist(hXiBkgFake,"#xi","1/N_{Jet}dN/d#xi");
  hXiBkgFake->Sumw2();
  hXiBkgFake->SetDirectory(0);
  hArray->AddLast(hXiBkgFake);

  hEFake=new TH1D("hEFake","p_{t,Jet} Fake "+JetAlgo,nBins,0,histPtMax);
  LabelHist(hEFake,"p_{t,Jet}^{rec.} Fake","# entries");
  hEFake->SetDirectory(0);
  hArray->AddLast(hEFake);

  hEPtNoPt=new TH1D("hEPtNoPt","p_{t,Jet} trig pt,cut - p_{t,Jet} matched no pt cut "+JetAlgo,40,-20,20);
  LabelHist(hEPtNoPt,"p_{t,Jet}^{rec.,trig}(pt cut)-p_{t,Jet}^{rec.,matched}(no pt cut)","# entries");
  hEPtNoPt->SetDirectory(0);
  hArray->AddLast(hEPtNoPt);
}

TH1D* ktAna::GetHist(TString hN,Int_t col, Int_t lStyle, Int_t lWidth)
{

  TH1D *temp;
  temp=(TH1D*) (hArray->FindObject(hN));
  temp->SetLineColor(col);
  temp->SetLineStyle(lStyle);
  temp->SetLineWidth(lWidth);

  return temp;//->Clone();
}

TH1D* ktAna::GetMarkerHist(TString hN,Int_t mStyle, Int_t mCol, Int_t mSize)
{
  TH1D *temp;
  temp=(TH1D*) (hArray->FindObject(hN));
  temp->SetMarkerStyle(mStyle);
  temp->SetLineColor(mCol);
  temp->SetMarkerColor(mCol);
  temp->SetMarkerSize(mSize);

  return temp;
}

TH1D* ktAna::GetXiSub(Int_t mStyle, Int_t mCol, Int_t mSize)
{
  TH1D *hXiCorr;
  hXiCorr=(TH1D*) hXi->Clone();
  
  hXiCorr->SetMarkerStyle(mStyle);
  hXiCorr->SetMarkerColor(mCol);
  hXiCorr->SetLineColor(mCol);
  hXiCorr->SetMarkerSize(mSize);

  hXiCorr->Add(hXiBkg,-1.0);

  return hXiCorr;
}

TH1D* ktAna::GetXiSubFake(Int_t mStyle, Int_t mCol, Int_t mSize)
{
  TH1D *hXiCorr;
  hXiCorr=(TH1D*) hXiFake->Clone();
  
  hXiCorr->SetMarkerStyle(mStyle);
  hXiCorr->SetMarkerColor(mCol);
  hXiCorr->SetLineColor(mCol);
  hXiCorr->SetMarkerSize(mSize);

  hXiCorr->Add(hXiBkgFake,-1.0);

  return hXiCorr;
}

TH1D* ktAna::GetDiJetXiSub(Int_t mStyle, Int_t mCol, Int_t mSize)
{
  TH1D *hXiCorrDi;
  hXiCorrDi=(TH1D*) hXiDi->Clone();
  
  hXiCorrDi->SetMarkerStyle(mStyle);
  hXiCorrDi->SetMarkerColor(mCol);
  hXiCorrDi->SetLineColor(mCol);
  hXiCorrDi->SetMarkerSize(mSize);

  hXiCorrDi->Add(hXiBkgDi,-1.0);

  return hXiCorrDi;
}

TH1D* ktAna::GetXiSubTrig(Int_t mStyle, Int_t mCol, Int_t mSize)
{
  TH1D *thXiCorr;
  thXiCorr=(TH1D*) thXi->Clone();
  
  thXiCorr->SetMarkerStyle(mStyle);
  thXiCorr->SetMarkerColor(mCol);
  thXiCorr->SetLineColor(mCol);
  thXiCorr->SetMarkerSize(mSize);

  thXiCorr->Add(thXiBkg,-1.0);

  return thXiCorr;
}

TH1D* ktAna::GetDiJetXiSubTrig(Int_t mStyle, Int_t mCol, Int_t mSize)
{
  TH1D *thXiCorrDi;
  thXiCorrDi=(TH1D*) thXiDi->Clone();
  
  thXiCorrDi->SetMarkerStyle(mStyle);
  thXiCorrDi->SetMarkerColor(mCol);
  thXiCorrDi->SetLineColor(mCol);
  thXiCorrDi->SetMarkerSize(mSize);

  thXiCorrDi->Add(thXiBkgDi,-1.0);

  return thXiCorrDi;
}

TH1D* ktAna::GetZSub(Int_t mStyle, Int_t mCol, Int_t mSize)
{
  TH1D *hZCorr;
  hZCorr=(TH1D*) hZ->Clone();
  
  hZCorr->SetMarkerStyle(mStyle);
  hZCorr->SetMarkerColor(mCol);
  hZCorr->SetLineColor(mCol);
  hZCorr->SetMarkerSize(mSize);

  hZCorr->Add(hZBkg,-1.0);

  return hZCorr;
}

// do proper z bkg. shift
TH1D* ktAna::GetDiJetZSub(Int_t mStyle, Int_t mCol, Int_t mSize)
{
  TH1D *hZCorrDi;
  hZCorrDi=(TH1D*) hZDi->Clone();
  
  hZCorrDi->SetMarkerStyle(mStyle);
  hZCorrDi->SetMarkerColor(mCol);
  hZCorrDi->SetLineColor(mCol);
  hZCorrDi->SetMarkerSize(mSize);

  hZCorrDi->Add(hZBkgDi,-1.0);

  return hZCorrDi;
}

TH1D* ktAna::GetZSubTrig(Int_t mStyle, Int_t mCol, Int_t mSize)
{
  TH1D *thZCorr;
  thZCorr=(TH1D*) thZ->Clone();
  
  thZCorr->SetMarkerStyle(mStyle);
  thZCorr->SetMarkerColor(mCol);
  thZCorr->SetLineColor(mCol);
  thZCorr->SetMarkerSize(mSize);

  thZCorr->Add(thZBkg,-1.0);

  return thZCorr;
}

// do proper z bkg. shift
TH1D* ktAna::GetDiJetZSubTrig(Int_t mStyle, Int_t mCol, Int_t mSize)
{
  TH1D *thZCorrDi;
  thZCorrDi=(TH1D*) thZDi->Clone();
  
  thZCorrDi->SetMarkerStyle(mStyle);
  thZCorrDi->SetMarkerColor(mCol);
  thZCorrDi->SetLineColor(mCol);
  thZCorrDi->SetMarkerSize(mSize);

  thZCorrDi->Add(thZBkgDi,-1.0);

  return thZCorrDi;
}


TH1D* ktAna::GetSpectrumCorr(Int_t mStyle, Int_t mCol, Int_t mSize)
{
  // dummy so far ...
  return 0;
}

TH1D* ktAna::GetSpectrumCorrTrig(Int_t mStyle, Int_t mCol, Int_t mSize)
{
  // dummy so far ...
  return 0;
}

void ktAna::Finish()
{
  ScaleHist(hXi,GetNjets(ffPtMin,ffPtMax));
  ScaleHist(hXiBkg,GetNjets(ffPtMin,ffPtMax));

  ScaleHist(hXiDi,GetNDiJets(diJetPtMin,diJetPtMax));
  ScaleHist(hXiBkgDi,GetNDiJets(diJetPtMin,diJetPtMax));

  ScaleHist(hXiFake,GetNFakeJets(diJetPtMin,diJetPtMax));
  ScaleHist(hXiBkgFake,GetNFakeJets(diJetPtMin,diJetPtMax));

  ScaleHist(hZ,GetNjets(ffPtMin,ffPtMax));
  ScaleHist(hZBkg,GetNjets(ffPtMin,ffPtMax));
  ScaleHist(hZDi,GetNDiJets(diJetPtMin,diJetPtMax));
  ScaleHist(hZBkgDi,GetNDiJets(diJetPtMin,diJetPtMax));

  ScaleHist(hSpec);
  ScaleHist(hSpecDi);

  ScaleHist(thXi,GetNjetsTrig(ffPtMin,ffPtMax));
  ScaleHist(thXiBkg,GetNjetsTrig(ffPtMin,ffPtMax));
  ScaleHist(thXiDi,GetNDiJetsTrig(diJetPtMin,diJetPtMax));
  ScaleHist(thXiBkgDi,GetNDiJetsTrig(diJetPtMin,diJetPtMax));

  ScaleHist(thZ,GetNjetsTrig(ffPtMin,ffPtMax));
  ScaleHist(thZBkg,GetNjetsTrig(ffPtMin,ffPtMax));
  ScaleHist(thZDi,GetNDiJetsTrig(diJetPtMin,diJetPtMax));
  ScaleHist(thZBkgDi,GetNDiJetsTrig(diJetPtMin,diJetPtMax));

  ScaleHist(thSpec);
  ScaleHist(thSpecDi);
}

void ktAna::DoFakeJetEstimate()
{
  cout<<endl;
  cout<<"Estimate fake-jet rate via di-jets:"<<endl;
  cout<<"-----------------------------------"<<endl;
  cout<<endl;
  cout<<"Use triggered jets also highest jet in event \nwith pt,rec > "<<minFakePtCut;
  cout<<" and Eta range = "<<etaMin<<" .. "<<etaMax<<endl;
  cout<<"Plot dPhi vs. pt,rec(second) of second highest jet."<<endl;
  cout<<"Define di-jet: "<<diPhiMin<<" < |dPhi| < "<<diPhiMax<<endl;
  cout<<"Estimate fake pedestal in "<<fakeMindPhi<<" < |dPhi| < "<<fakeMaxdPhi<<endl;
  cout<<"and extrapolate uniformily over dPhi."<<endl;
  cout<<endl;

  // Rebin ...
  //hdPhiPt->Rebin2D(3,3);

  // dphiCut on away side
  Double_t minAway=diPhiMin;//TMath::Pi()-Raway;
  Double_t maxAway=diPhiMax;//TMath::Pi()+Raway;

  Double_t awayScale=TMath::Abs(minAway-maxAway);

  // fake area
  Double_t minFake=fakeMindPhi;
  Double_t maxFake=fakeMaxdPhi;

  Double_t fakeScale=TMath::Abs(minFake-maxFake);

  // fake area
  Double_t maxFake2=TMath::Pi()*2-minFake;
  Double_t minFake2=maxFake2-fakeScale;

  //DEBUG:
  //cout<<minFake<<" "<<maxFake<<endl;
  //cout<<minFake2<<" "<<maxFake2<<endl;

  // get bins ...
  Int_t minAwayBin=hdPhiPt->GetXaxis()->FindBin(minAway);
  Int_t maxAwayBin=hdPhiPt->GetXaxis()->FindBin(maxAway);
  Int_t minFakeBin=hdPhiPt->GetXaxis()->FindBin(minFake);
  Int_t maxFakeBin=hdPhiPt->GetXaxis()->FindBin(maxFake);
  Int_t minFakeBin2=hdPhiPt->GetXaxis()->FindBin(minFake2);
  Int_t maxFakeBin2=hdPhiPt->GetXaxis()->FindBin(maxFake2);

  //do projections ...

  TH1D *hAway=new TH1D();TH1D *hFake=new TH1D();TH1D *hFake2=new TH1D();
  hAway->SetDirectory(0);hFake->SetDirectory(0);hFake2->SetDirectory(0);

  hAway=(TH1D*) hdPhiPt->ProjectionY("hAway",minAwayBin,maxAwayBin,"E");
  hFake=(TH1D*) hdPhiPt->ProjectionY("hFake",minFakeBin,maxFakeBin,"E");
  hFake2=(TH1D*)hdPhiPt->ProjectionY("hFake2",minFakeBin2,maxFakeBin2,"E");

  hFake->Rebin(3);hAway->Rebin(3);hFake2->Rebin(3);

  // do fake estimate ...
  hFake->Add(hFake2);
  hFake->Scale(1/(2*fakeScale));
  //hFake->Scale(awayScale);// for di-jets !!! think a bit more !!!
  hFake->Scale(2*TMath::Pi());

  TF1 *sbcorr=new TF1("sbcorr","[0]/(2*TMath::Pi())",0,75);
  sbcorr->SetParameter(0,awayScale);
  
  //Show some plots ...
  TCanvas *cktAnaFake=new TCanvas("cktAnaFake"+JetAlgo,"ktAna Fake-Jet plots "+JetAlgo,800,800);
  cktAnaFake->Divide(2,2);
  cktAnaFake->cd(1);
  hdPhiPt->DrawClone("colz");
  cktAnaFake->cd(2);
  hdPhiPt->ProjectionX()->DrawClone();

  cktAnaFake->cd(3);
  gPad->SetLogy(1);
  hAway->SetLineColor(1);
  hFake->SetLineColor(1);
  hFake->SetLineStyle(2);
  hFake->DrawClone("h");
  hAway->DrawClone("hsame");

  cktAnaFake->cd(4);
  gPad->SetLogy(1);
  hAway->GetYaxis()->SetTitle("\"True\" Jets/\"Fake\" Jets");
  hAway->GetXaxis()->SetTitle("p_{t,rec}^{Jet} [GeV/c]");
  hAway->Divide(hFake);
  hAway->DrawClone("");
  hAway->Add(sbcorr,-1); // check if correct ...  !!!!
  hAway->SetLineColor(2);
  hAway->DrawClone("same");

  /*
  hFakeJetRatio=(TH1D*) hAway->Clone();
  hFakeJetRatio->SetTitle("hFakeJetRatio");
  hFakeJetRatio->SetName("hFakeJetRatio");
  hFakeJetRatio->SetDirectory(0);
  hArray->AddLast(hFakeJetRatio);
  */

  delete hAway;delete hFake;delete hFake2;
}

void ktAna::DrawFFPlotsFake()
{
  TCanvas *cktAnaFFFake=new TCanvas("cktAnaFFFake"+JetAlgo,"ktAna FF plots fake "+JetAlgo,600,600);
  cktAnaFFFake->Divide(2,2);
  cktAnaFFFake->cd(1);
  gPad->SetLogy(1);
  (TH1D*) (GetXi(22,1,1))->DrawClone();
  (TH1D*) (GetXiBkg(2,2,1))->DrawClone("same");
  cktAnaFFFake->cd(2);
  (TH1D*) (GetXiSub(23,2,1))->DrawClone();
  cktAnaFFFake->cd(3);
  gPad->SetLogy(1);
  (TH1D*) (GetXiFake(22,1,1))->DrawClone();
  (TH1D*) (GetXiBkgFake(2,2,1))->DrawClone("same");
  cktAnaFFFake->cd(4);
  (TH1D*) (GetXiSubFake(23,2,1))->DrawClone("");
  (TH1D*) (GetXiSub(22,4,1))->DrawClone("same");
  (TH1D*) (GetDiJetXiSub(24,1,1))->DrawClone("same");
}

void ktAna::DrawPlotsPtNoPt()
{
  TCanvas *cktAnaPtNoPt=new TCanvas("cktAnaPtNoPt"+JetAlgo,"ktAna FF plots "+JetAlgo,800,800);
  cktAnaPtNoPt->Divide(2,3);
  cktAnaPtNoPt->cd(1);
  (TH1D*)  GetdEdi()->DrawClone();
  cktAnaPtNoPt->cd(2);
  hEPtNoPt->DrawClone();
  cktAnaPtNoPt->cd(3);
  (TH1D*) (GetXi(22,1,1))->DrawClone();
  (TH1D*) (GetXiBkg(2,2,1))->DrawClone("same");
  cktAnaPtNoPt->cd(4);
  (TH1D*) (GetDiJetXi(22,1,1))->DrawClone();
  (TH1D*) (GetDiJetXiBkg(2,2,1))->DrawClone("same");
  cktAnaPtNoPt->cd(5);
  (TH1D*) (GetXiSub(23,2,1))->DrawClone("");
  cktAnaPtNoPt->cd(6);
  (TH1D*) (GetDiJetXiSub(23,2,1))->DrawClone("");
}

void ktAna::DrawFFPlots()
{
  TCanvas *cktAnaFF=new TCanvas("cktAnaFF"+JetAlgo,"ktAna FF plots "+JetAlgo,800,800);
  cktAnaFF->Divide(2,4);
  cktAnaFF->cd(1);
  (TH1D*) (GetXi(22,1,1))->DrawClone();
  (TH1D*) (GetXiBkg(2,2,1))->DrawClone("same");
  cktAnaFF->cd(2);
  gPad->SetLogy(1);
  (TH1D*) (GetZ(22,1,1))->DrawClone();
  (TH1D*) (GetZBkg(2,2,1))->DrawClone("same");
  cktAnaFF->cd(3);
  (TH1D*) (GetXiSub(23,2,1))->DrawClone("");
  cktAnaFF->cd(4);
  gPad->SetLogy(1);
  (TH1D*) (GetZSub(23,2,1))->DrawClone("");
  cktAnaFF->cd(5);
  (TH1D*) (GetDiJetXi(22,1,1))->DrawClone();
  (TH1D*) (GetDiJetXiBkg(2,2,1))->DrawClone("same");
  cktAnaFF->cd(6);
  gPad->SetLogy(1);
  (TH1D*) (GetDiJetZ(22,1,1))->DrawClone();
  (TH1D*) (GetDiJetZBkg(2,2,1))->DrawClone("same");
  cktAnaFF->cd(7);
  (TH1D*) (GetDiJetXiSub(23,2,1))->DrawClone("");
  cktAnaFF->cd(8);
  gPad->SetLogy(1);
  (TH1D*) (GetDiJetZSub(23,2,1))->DrawClone("");
  
}

void ktAna::DrawFFPlotsTrig()
{
  TCanvas *cktAnaFFtrig=new TCanvas("cktAnaFFtrig"+JetAlgo,"ktAna FF plots trig "+JetAlgo,800,800);
  cktAnaFFtrig->Divide(2,4);
  cktAnaFFtrig->cd(1);
  (TH1D*) (GetXiTrig(22,1,1))->DrawClone();
  (TH1D*) (GetXiBkgTrig(2,2,1))->DrawClone("same");
  cktAnaFFtrig->cd(2);
  gPad->SetLogy(1);
  (TH1D*) (GetZTrig(22,1,1))->DrawClone();
  (TH1D*) (GetZBkgTrig(2,2,1))->DrawClone("same");
  cktAnaFFtrig->cd(3);
  (TH1D*) (GetXiSubTrig(23,2,1))->DrawClone("");
  cktAnaFFtrig->cd(4);
  gPad->SetLogy(1);
  (TH1D*) (GetZSubTrig(23,2,1))->DrawClone("");
  cktAnaFFtrig->cd(5);
  (TH1D*) (GetDiJetXiTrig(22,1,1))->DrawClone();
  (TH1D*) (GetDiJetXiBkgTrig(2,2,1))->DrawClone("same");
  cktAnaFFtrig->cd(6);
  gPad->SetLogy(1);
  (TH1D*) (GetDiJetZTrig(22,1,1))->DrawClone();
  (TH1D*) (GetDiJetZBkgTrig(2,2,1))->DrawClone("same");
  cktAnaFFtrig->cd(7);
  (TH1D*) (GetDiJetXiSubTrig(23,2,1))->DrawClone("");
  cktAnaFFtrig->cd(8);
  gPad->SetLogy(1);
  (TH1D*) (GetDiJetZSubTrig(23,2,1))->DrawClone("");
  
}

void ktAna::DrawPlots()
{
  TCanvas *cktAna=new TCanvas("cktAna"+JetAlgo,"ktAna plots "+JetAlgo,800,800);
  cktAna->Divide(3,4);
  cktAna->cd(1);
  gPad->SetLogy(1);
  (TH1D*) (GetHist("hE",1,1,1))->DrawClone();
  (TH1D*) (GetHist("hEDi",2,1,1))->DrawClone("same");
  //(TH1D*) (GetHist("hEFake",4,1,1))->DrawClone("same");
  cktAna->cd(2);
  gPad->SetLogy(1);
  (TH1D*) (GetSpectrum(22,1,1))->DrawClone();
  (TH1D*) (GetDiJetSpectrum(23,2,1))->DrawClone("same");
  cktAna->cd(3);
  (TH1D*) (GetXi(22,1,1))->DrawClone();
  (TH1D*) (GetXiBkg(2,2,1))->DrawClone("same");
  cktAna->cd(4);
  (TH1D*) (GetXiSub(23,2,1))->DrawClone("");
  cktAna->cd(5);
  (TH1D*) (GetDiJetXi(22,1,1))->DrawClone();
  (TH1D*) (GetDiJetXiBkg(2,2,1))->DrawClone("same");
  cktAna->cd(6);
  (TH1D*) (GetDiJetXiSub(23,2,1))->DrawClone("");
  cktAna->cd(7);
  hEtaPhi->DrawClone("colz");
  cktAna->cd(8);
  hEtaPhiDi->DrawClone("colz");
  cktAna->cd(9);
  //(TH1D*) (GetHist("hPhiDi",1,1,1))->DrawClone();
  hdPhiPt->DrawClone("colz");
  cktAna->cd(10);
  hAreaPt->DrawClone("colz");
  cktAna->cd(11);
  hAreaPtDi->DrawClone("colz");
  cktAna->cd(12);
  //(TH1D*) (GetHist("hPtDiff",1,1,1))->DrawClone();
  (TH1D*)  GetdEdi()->DrawClone();
}

void ktAna::AnalyzeEvent(ktMuEvent *ev,TString m_JetAlgo)
{
  JetAlgo=m_JetAlgo;
  AnalyzeEvent(ev);
}

Int_t ktAna::NAlgoJets(ktMuEvent *ev)
{
  Int_t mNjets=0;
      
  if (JetAlgo.Contains("SIS"))
    {
      mNjets=ev->GetNFJSISCone();
    }
  else if (JetAlgo.Contains("Anti"))
    {
      mNjets=ev->GetNFJAntiKt();
    }
  else if (JetAlgo.Contains("Cone"))
    {
      mNjets=ev->GetNConeJets();
      FJAna=false;
    }
  else
    {
      mNjets=ev->GetNFJkt();
    }
  
  return mNjets;
}

TH1D* ktAna::GetXiBkg(ktMuEvent *ev)
{
  TH1D *htemp=0;
  
  if (JetAlgo.Contains("SIS"))
    {
      htemp=ev->GetXiBkgSISCone();
    }
  else if (JetAlgo.Contains("Anti"))
    {
      htemp=ev->GetXiBkgAntiKt();
    }
  else if (JetAlgo.Contains("Cone"))
    {
      htemp=ev->GetXiBkg();
    }
  else
    {
      htemp=ev->GetXiBkgKt();
    }

  return htemp;
}

TH1D* ktAna::GetZBkg(ktMuEvent *ev)
{
  TH1D *htemp=0;
  
  if (JetAlgo.Contains("SIS"))
    {
      htemp=ev->GetZBkgSISCone();
    }
  else if (JetAlgo.Contains("Anti"))
    {
      htemp=ev->GetZBkgAntiKt();
    }
  else if (JetAlgo.Contains("Cone"))
    {
      htemp=ev->GetZBkg();
    }
  else
    {
      htemp=ev->GetZBkgKt();
    }

  return htemp;
}

TH1D* ktAna::GetptBkg(ktMuEvent *ev)
{
  TH1D *htemp=0;
  
  if (JetAlgo.Contains("SIS"))
    {
      htemp=ev->GetptBkgSISCone();
    }
  else if (JetAlgo.Contains("Anti"))
    {
      htemp=ev->GetptBkgAntiKt();
    }
  else if (JetAlgo.Contains("Cone"))
    {
      htemp=ev->GetptBkg();
    }
  else
    {
      htemp=ev->GetptBkgKt();
    }

  return htemp;
}

ktMuFastJet* ktAna::GetFJjet(ktMuEvent *ev,Int_t n)
{

  ktMuFastJet *temp;

   if (JetAlgo.Contains("SIS"))
    {
      temp=ev->GetFJSISCone(n);
    }
  else if (JetAlgo.Contains("Anti"))
    {
      temp=ev->GetFJAntiKt(n);
    }
  else
    {
      temp=ev->GetFJkt(n);
    }

  return temp;
}

ktMuJet* ktAna::GetJet(ktMuEvent *ev,Int_t n)
{
  ktMuJet *temp=(ktMuJet*) ev->GetConeJet(n);
  return temp;
}

Bool_t ktAna::CheckAcc(ktMuFastJet *j)
{
  if (j->Phi()>phiMin && j->Phi()<phiMax && j->Eta()>etaMin && j->Eta()<etaMax)
    return true;
  else
    return false;
}

Bool_t ktAna::CheckAcc(ktMuJet *j)
{
  if (j->Phi()>phiMin && j->Phi()<phiMax && j->Eta()>etaMin && j->Eta()<etaMax)
    return true;
  else
    return false;
}

void ktAna::CompareEvent(ktAna *e)
{
  hPtDiff->SetTitle("p_{t,Jet}^{rec.}("+JetAlgo+")-p_{t,Jet}^{rec.}("+e->GetJetAlgo()+")");
  LabelHist(hPtDiff,"p_{t,Jet}^{rec.}("+JetAlgo+")-p_{t,Jet}^{rec.}("+e->GetJetAlgo()+")","# entries");
  hPtDiffDi->SetTitle("p_{t,Jet}^{rec.}("+JetAlgo+")-p_{t,Jet}^{rec.}("+e->GetJetAlgo()+")");
  LabelHist(hPtDiffDi,"p_{t,Jet}^{rec.}("+JetAlgo+")-p_{t,Jet}^{rec.}("+e->GetJetAlgo()+")","# entries");

  if (GetEventPt()>0 && e->GetEventPt()>0)
    {
      hPtDiff->Fill(GetEventPt()-e->GetEventPt());
      if (GetEventDiJetPt()>0 && e->GetEventDiJetPt()>0)
	hPtDiffDi->Fill(GetEventDiJetPt()-e->GetEventDiJetPt());
    }
}


void ktAna::CompareEventTrig(ktAna *e)
{
  thPtDiff->SetTitle("p_{t,Jet}^{rec.}("+JetAlgo+")-p_{t,Jet}^{rec.}("+e->GetJetAlgo()+")");
  LabelHist(thPtDiff,"p_{t,Jet}^{rec.}("+JetAlgo+")-p_{t,Jet}^{rec.}("+e->GetJetAlgo()+")","# entries");
  thPtDiffDi->SetTitle("p_{t,Jet}^{rec.}("+JetAlgo+")-p_{t,Jet}^{rec.}("+e->GetJetAlgo()+")");
  LabelHist(thPtDiffDi,"p_{t,Jet}^{rec.}("+JetAlgo+")-p_{t,Jet}^{rec.}("+e->GetJetAlgo()+")","# entries");

  if (GetEventPtTrig()>0 && e->GetEventPtTrig()>0)
    {
      thPtDiff->Fill(GetEventPtTrig()-e->GetEventPtTrig());
      if (GetEventDiJetPtTrig()>0 && e->GetEventDiJetPtTrig()>0)
	thPtDiffDi->Fill(GetEventDiJetPtTrig()-e->GetEventDiJetPtTrig());
    }
}

void ktAna::AnalyzeEventTrig(ktMuEvent *ev)
{
  if (verbose)
    cout<<"Analyze event # "<<ev->GetEventNumber()<<endl;

  nEvents++;

  if (ev!=0)
    {
      nJets=NAlgoJets(ev);
      
      if (nJets>0)
	{
	  if (verbose)
	    cout<<"# of "<<JetAlgo<<" jets = "<<nJets<<endl;

	  if (FJAna)
	    {
	      // ================
	      // FastJet analysis 
	      // ================

	      ktMuFastJet* mJ=(ktMuFastJet*) GetFJjet(ev,0);
	      // DEBUG:
	      //mJ->PrintJet();

	      TH1D *mXiBkg=0;
	      TH1D *mXiBkgDi=0; 
	      TH1D *mZBkg=0;
	      TH1D *mZBkgDi=0;  
	      TH1D *mptBkg=0;
	      TH1D *mptBkgDi=0;  
	      
	      if (mJ->GetNJetParticles()>0)
		{
		  
		  mXiBkg=new TH1D("mXiBkg","Xi dist. of Background per event FastJet",(Int_t) 50,0,10);
		  mXiBkgDi=new TH1D("mXiBkgDi","Xi dist. of Background per event FastJet Di",(Int_t) 50,0,10);	
		  mZBkg=new TH1D("mZBkg","Z dist. of Background per event FastJet",(Int_t) 10,0,1);
		  mZBkgDi=new TH1D("mZBkgDi","Z dist. of Background per event FastJet Di",(Int_t) 10,0,1);
		  
		  mptBkg=(TH1D*) GetptBkg(ev)->Clone();
		  mptBkgDi=(TH1D*) GetptBkg(ev)->Clone();
		  
		  // DEBUG:
		  //cout<<mZBkg->GetNbinsX()<<endl;
		}	   
	      
	      Double_t jetPt=0;
	      if (BkgCorr)
		jetPt=mJ->PtCorr();
	      else
		jetPt=mJ->Pt();
	      
	      if (jetPt>minPtCut && CheckAcc(mJ))
		{
		  
		}	    
	      else
		{
		  if (verbose)
		    cout<<"No "<<JetAlgo<<" jet above pt,min = "<<minPtCut<<" found !"<<endl;
		  evPtDi=evPt=-99.0;
		  evAreaDi=evArea=-1.0;
		}

	      delete mXiBkg;
	      delete mXiBkgDi;
	      delete mZBkg;
	      delete mZBkgDi;
	      delete mptBkg;
	      delete mptBkgDi;
	    }
	}
    }
}

void ktAna::DoXiBkgFromPt(TH1D *hxi,TH1D *hpt,Double_t jpt)
{
  for (int i=1;i<hpt->GetNbinsX();i++)
    {
      Double_t pt=hpt->GetBinCenter(i);
      Double_t dn=hpt->GetBinContent(i);
      Double_t z=pt/jpt;
      Double_t xi=TMath::Log(1/z);

      //DEBUG
      //if (dn>0)
      //cout<<xi<<" "<<pt<<" "<<dn<<" "<<jpt<<endl;

      //hxi->Fill(xi,TMath::Sqrt(dn));
      
      hxi->SetBinContent(hxi->FindBin(xi),hxi->GetBinContent(hxi->FindBin(xi))+dn);
      hxi->SetBinError(hxi->FindBin(xi),TMath::Sqrt(hxi->GetBinContent(hxi->FindBin(xi))));      
      
    }
}

void ktAna::DoZBkgFromPt(TH1D *hz,TH1D *hpt,Double_t jpt)
{
  for (int i=1;i<hpt->GetNbinsX();i++)
    {
      Double_t pt=hpt->GetBinCenter(i);
      Double_t dn=hpt->GetBinContent(i);
      Double_t z=pt/jpt;
      
      //DEBUG
      //if (dn>0)
      //cout<<z<<" "<<pt<<" "<<dn<<" "<<jpt<<endl;
      
      //hz->Fill(z,TMath::Sqrt(dn));
      
      hz->SetBinContent(hz->FindBin(z),hz->GetBinContent(hz->FindBin(z))+dn);
      hz->SetBinError(hz->FindBin(z),TMath::Sqrt(hz->GetBinContent(hz->FindBin(z))));
      
    }
}


void ktAna::FillFFHists(TH1D *mhXi,TH1D *mhZ,ktMuFastJet *mJ, Double_t jetPt)
{

  if (mJ->GetNJetParticles()>0)
    {
      for (int i=0;i<mJ->GetNJetParticles();i++)
	{
	  ktParticle* p=(ktParticle*) mJ->GetJetParticle(i);
	  
	  if (p->IsCharged()) // && p->GetPid()->GetPID()!=11)
	    {
	      Double_t xi=TMath::Log((Double_t) jetPt/(Double_t) p->Pt());
	      Double_t z=(Double_t) p->Pt()/(Double_t) jetPt;
	      
	      if (mhXi!=0)
		mhXi->Fill(xi);

	      if (mhZ!=0)
		mhZ->Fill(z);
	    }		
	}		  
    }
  else
    {
      TH1D *hpt=(TH1D*) mJ->GetJetPtHistogram();
      // DEBUG:
      //cout<<"Got hist = "<<hpt<<endl;
      //cout<<jetPt<<endl;

      if (hpt!=0)
	{
	  for (int i=1;i<hpt->GetNbinsX();i++)
	    {
	      Double_t pt=hpt->GetBinCenter(i);
	      Double_t dn=hpt->GetBinContent(i);
	      Double_t z=pt/jetPt;
	      Double_t xi=TMath::Log(1/z);
	      
	      /*
		hXi->Fill(xi,TMath::Sqrt(dn));
		hZ->Fill(z,TMath::Sqrt(dn)); 
	      */
	      
	      //cout<<pt<<" "<<dn<<endl;

	      if (mhXi!=0)
		{
		  mhXi->SetBinContent(mhXi->FindBin(xi),mhXi->GetBinContent(mhXi->FindBin(xi))+dn);
		  mhXi->SetBinError(mhXi->FindBin(xi),TMath::Sqrt(mhXi->GetBinContent(mhXi->FindBin(xi))));
		}
	      if (mhZ!=0)
		{
		  mhZ->SetBinContent(mhZ->FindBin(z),mhZ->GetBinContent(mhZ->FindBin(z))+dn);
		  mhZ->SetBinError(mhZ->FindBin(xi),TMath::Sqrt(mhZ->GetBinContent(mhZ->FindBin(xi))));
		}
	      
	    }     
	}
    }  
}

void ktAna::AnalyzeEventPtAndNoPt(ktMuEvent *ev,ktMuEvent *ev2)
{
  if (verbose)
    cout<<"Analyze event # "<<ev->GetEventNumber()<<endl;
  
  nEvents++;

  if (ev!=0 && ev2!=0)
    {
      Int_t nJets1=NAlgoJets(ev);
      Int_t nJets2=NAlgoJets(ev2);

      if (nJets1>0 && nJets2>0)
	{
	  if (verbose)
	    {
	      cout<<"Event # = "<<nEvents<<endl;
	      cout<<"# of "<<JetAlgo<<" jets pt,cut    = "<<nJets1<<endl;
	      cout<<"# of "<<JetAlgo<<" jets no pt,cut = "<<nJets2<<endl;
	    }
	  
	  // ================
	  // FastJet analysis 
	  // ================
	  
	  ktMuFastJet* mJ=(ktMuFastJet*) GetFJjet(ev,0);
	  
	  TH1D *mXiBkg=0;
	  TH1D *mXiBkgDi=0; 
	  TH1D *mXiBkgFake=0;
	  TH1D *mZBkg=0;
	  TH1D *mZBkgDi=0;  
	  TH1D *mptBkg=0;
	  TH1D *mptBkgFake=0;
	  TH1D *mptBkgDi=0;  
	  
	  
	  // use background from pt,cut>2 !!!!
	  // check with other approach to use pt,cut>0 for recoil then ...
	  
	  if (mJ->GetNJetParticles()>0 || GetptBkg(ev)!=0)
	    {	
	      mXiBkg=new TH1D("mXiBkg","Xi dist. of Background per event FastJet",(Int_t) 50,0,10);
	      mXiBkgDi=new TH1D("mXiBkgDi","Xi dist. of Background per event FastJet Di",(Int_t) 50,0,10);	
	      mXiBkgFake=new TH1D("mXiBkgFake","Xi dist. of Background per event FastJet Fake",(Int_t) 50,0,10);
	      mZBkg=new TH1D("mZBkg","Z dist. of Background per event FastJet",(Int_t) 10,0,1);
	      mZBkgDi=new TH1D("mZBkgDi","Z dist. of Background per event FastJet Di",(Int_t) 10,0,1);
	      
	      mptBkg=(TH1D*) GetptBkg(ev)->Clone();
	      mptBkgFake=(TH1D*) GetptBkg(ev)->Clone();
	      mptBkgDi=(TH1D*) GetptBkg(ev)->Clone();
	    }	
	  
	  Double_t jetPt=0;
	  if (BkgCorr)
	    jetPt=mJ->PtCorr();
	  else
	    jetPt=mJ->Pt();
	  
	    // check if highest jet for pt,cut>2 GeV is trigger jet
	  if (jetPt>minPtCut && CheckAcc(mJ) && (mJ->IsTriggerJet() || findTrigMatch(ev,mJ)))
	    {
	      hE->Fill(jetPt);
	      
	      hEtaPhi->Fill(mJ->Eta(),mJ->Phi());
	      hAreaPt->Fill(jetPt,mJ->Area());
	      
	      evPt=jetPt;
	      evArea=mJ->Area();
	      
	      if (jetPt>ffPtMin && jetPt<ffPtMax)
		{	
		  hSpec->Fill(jetPt);
		  
		  if (mptBkg!=0 && mptBkg->GetEntries()>0)		       
		    {
		      DoZBkgFromPt(mZBkg,mptBkg,jetPt);
		      mZBkg->Scale(mJ->GetFFAreaInAcc());			
		      hZBkg->Add(mZBkg);
		      
		      DoXiBkgFromPt(mXiBkg,mptBkg,jetPt);
		      mXiBkg->Scale(mJ->GetFFAreaInAcc());
		      hXiBkg->Add(mXiBkg);
		    }	
		  
		  FillFFHists(hXi,hZ,mJ,jetPt);
		  
		  // if found a trigger jet fullfilling all
		  // cuts for the FF, look for the recoil
		  // in the no,pt cut event !
		  // order by pt, so take the highest jet
		  // fullfilling di-jet criteria ...

		  for (int i=0;i<nJets2;i++)
		    {
		      ktMuFastJet* mJ2=(ktMuFastJet*) GetFJjet(ev2,i);

		      Double_t jetPt2=0;
		      if (BkgCorr)
			jetPt2=mJ2->PtCorr();
		      else
			jetPt2=mJ2->Pt();

		      if (jetPt2<minPtCut) break;

		      Double_t dphi=(GetdPhi(mJ->Phi(),mJ2->Phi()));
		      Double_t dphi2=(GetdPhi(mJ->Phi(),mJ2->Phi()));

		      if (dphi<0) dphi += (2*TMath::Pi());
		      Double_t dEta=mJ->Eta()-mJ2->Eta();
		      
		       // DEBUG:
		      //cout<<nEvents<<" "<<jetPt<<" "<<jetPt2<<" "<<dphi2<<endl;

		      if (TMath::Sqrt(dEta*dEta+dphi2*dphi2)>0.4)
			{
			  hPhiDi->Fill(dphi);
			  if (jetPt>minFakePtCut)
			    hdPhiPt->Fill(dphi,jetPt2);
			  
			  //DEBUG:
			  //cout<<i<<" "<<dphi<<" "<<jetPt2<<endl;
			}
		      else
			{
			  // plot pt,rec diff. pt,cut>2 and no pt cut ...
			  hEPtNoPt->Fill(jetPt-jetPt2);
			}
		      
		      // check for di-jet and pt,cuts for no-pt cut found jets
		      // with respect to tiggered jet found with pt,cut>2 GeV
		      if (jetPt2>diJetPtMin && jetPt2<diJetPtMax && CheckAcc(mJ2))
			{
			  if (IsDiJet(dphi))
			    {
			      hdEdi->Fill(jetPt-jetPt2);			    
			      
			      hEDi->Fill(jetPt2);
			      hSpecDi->Fill(jetPt2);
			      hEtaPhiDi->Fill(mJ2->Eta(),mJ2->Phi());
			      hAreaPtDi->Fill(jetPt2,mJ2->Area());
			      
			      evPtDi=jetPt2;
			      evArea=mJ2->Area();	
			      			        
			      if (TriggerPtToRecoil)
				jetPt2=jetPt;

			      if (mptBkgDi!=0 && mptBkgDi->GetEntries()>0)		       
				{
				  DoZBkgFromPt(mZBkgDi,mptBkgDi,jetPt2);
				  mZBkgDi->Scale(mJ2->GetFFAreaInAcc());
				  hZBkgDi->Add(mZBkgDi);
				  
				  DoXiBkgFromPt(mXiBkgDi,mptBkg,jetPt2);
				  mXiBkgDi->Scale(mJ2->GetFFAreaInAcc());
				  hXiBkgDi->Add(mXiBkgDi);
				}
			      
			      FillFFHists(hXiDi,hZDi,mJ2,jetPt2);
			    }
			  else if (IsFakeJet(dphi))
			    {
			      // look for fake jets then with no-pt cut ...

			      if (TriggerPtToRecoil)
				jetPt2=jetPt;

			      hEFake->Fill(jetPt2);

			      if (mptBkgDi!=0 && mptBkgDi->GetEntries()>0)		       
				{
				  DoXiBkgFromPt(mXiBkgFake,mptBkgFake,jetPt2);
				  mXiBkgFake->Scale(mJ2->GetFFAreaInAcc());
				  hXiBkgFake->Add(mXiBkgFake);
				}
			      
			      FillFFHists(hXiFake,0,mJ2,jetPt2);
			    }
			  else
			    {
			      evPtDi=-99.0;
			      evAreaDi=-1.0;
			    }	
			}
		      else
			{
			  evPtDi=-99.0;
			  evAreaDi=-1.0;
			}
		    }
		}
	    }
	  

	  delete mXiBkg;
	  delete mXiBkgDi;
	  delete mXiBkgFake;
	  delete mZBkg;
	  delete mZBkgDi;
	  delete mptBkg;
	  delete mptBkgDi;
	  delete mptBkgFake;
	  
	  // end ...
	}
      else
	{
	  if (verbose)
	    cout<<"No "<<JetAlgo<<" jet found !"<<endl;
	  
	  evPtDi=evPt=-99.0;
	  evAreaDi=evArea=-1.0;
	}
    }
  else
    {
      if (verbose) 
	cout<<"Warning : Zero pointer to ktMuEvent !"<<endl;
      
      evPtDi=evPt=-99.0;
      evAreaDi=evArea=-1.0;
    }
}

void ktAna::AnalyzeEvent(ktMuEvent *ev)
{
  if (verbose)
    cout<<"Analyze event # "<<ev->GetEventNumber()<<endl;

  nEvents++;

  if (ev!=0)
    {
      nJets=NAlgoJets(ev);

      if (nJets>0)
	{
	  if (verbose)
	    cout<<"# of "<<JetAlgo<<" jets = "<<nJets<<endl;

	  if (FJAna)
	    {
	      // ================
	      // FastJet analysis 
	      // ================

	      ktMuFastJet* mJ=(ktMuFastJet*) GetFJjet(ev,0);
	   
	      TH1D *mXiBkg=0;
	      TH1D *mXiBkgDi=0; 
	      TH1D *mXiBkgFake=0;
	      TH1D *mZBkg=0;
	      TH1D *mZBkgDi=0;  
	      TH1D *mptBkg=0;
	      TH1D *mptBkgFake=0;
	      TH1D *mptBkgDi=0;  
	      
	      if (mJ->GetNJetParticles()>0 || GetptBkg(ev)!=0)
		{	
		  mXiBkg=new TH1D("mXiBkg","Xi dist. of Background per event FastJet",(Int_t) 50,0,10);
		  mXiBkgDi=new TH1D("mXiBkgDi","Xi dist. of Background per event FastJet Di",(Int_t) 50,0,10);	
		  mXiBkgFake=new TH1D("mXiBkgFake","Xi dist. of Background per event FastJet Fake",(Int_t) 50,0,10);
		  mZBkg=new TH1D("mZBkg","Z dist. of Background per event FastJet",(Int_t) 10,0,1);
		  mZBkgDi=new TH1D("mZBkgDi","Z dist. of Background per event FastJet Di",(Int_t) 10,0,1);

		  mptBkg=(TH1D*) GetptBkg(ev)->Clone();
		  mptBkgFake=(TH1D*) GetptBkg(ev)->Clone();
		  mptBkgDi=(TH1D*) GetptBkg(ev)->Clone();
		}	   
	 
	      Double_t jetPt=0;
	      if (BkgCorr)
		jetPt=mJ->PtCorr();
	      else
		jetPt=mJ->Pt();

	      // check if IsTriggerJet correct ...
	      if (jetPt>minPtCut && CheckAcc(mJ) && (mJ->IsTriggerJet() || findTrigMatch(ev,mJ)))
		{
		  hE->Fill(jetPt);
		  
		  hEtaPhi->Fill(mJ->Eta(),mJ->Phi());
		  hAreaPt->Fill(jetPt,mJ->Area());

		  evPt=jetPt;
		  evArea=mJ->Area();

		  if (jetPt>ffPtMin && jetPt<ffPtMax)
		    {	
		      hSpec->Fill(jetPt);
	   
		      if (mptBkg!=0 && mptBkg->GetEntries()>0)		       
			{
			  DoZBkgFromPt(mZBkg,mptBkg,jetPt);
			  mZBkg->Scale(mJ->GetFFAreaInAcc());			
			  hZBkg->Add(mZBkg);

			  DoXiBkgFromPt(mXiBkg,mptBkg,jetPt);
			  mXiBkg->Scale(mJ->GetFFAreaInAcc());
			  hXiBkg->Add(mXiBkg);
			}	

		      FillFFHists(hXi,hZ,mJ,jetPt);
		 
		    }
		  
		  // -------------------
		  // Check if dijets ... 
		  // -------------------
		  
		  if (nJets>1)
		    {
		      ktMuFastJet* mJ1=(ktMuFastJet*) GetFJjet(ev,0);
		      ktMuFastJet* mJ2=(ktMuFastJet*) GetFJjet(ev,1);
		      //Double_t dphi=TMath::Abs(mJ1->Phi()-mJ2->Phi());
		      Double_t dphi=(GetdPhi(mJ1->Phi(),mJ2->Phi()));
		      if (dphi<0) dphi += (2*TMath::Pi());

		      hPhiDi->Fill(dphi);
		     
		      Double_t jetPt1=0;
		      Double_t jetPt2=0;
		      if (BkgCorr)
			{
			  jetPt1=mJ1->PtCorr();
			  jetPt2=mJ2->PtCorr();
			}
		      else
			{
			  jetPt1=mJ1->Pt();
			  jetPt2=mJ2->Pt();
			}

		      if (jetPt1>minFakePtCut)
			hdPhiPt->Fill(dphi,jetPt2);

		      //if (IsDiJet(dphi) && jetPt1>diJetPtMin && jetPt1<diJetPtMax && jetPt2>diJetPtMin && jetPt2<diJetPtMax && CheckAcc(mJ1) && CheckAcc(mJ2))
		      if (jetPt1>ffPtMin && jetPt1<ffPtMax && jetPt2>diJetPtMin && jetPt2<diJetPtMax && CheckAcc(mJ1) && CheckAcc(mJ2))
			{

			  if (IsDiJet(dphi))
			    {
					    
			      hdEdi->Fill(jetPt1-jetPt2);			    	      
			      // DEBUG:
			      //cout<<jetPt1<<" "<<jetPt2<<" "<<dphi<<endl;
			      
			      if (TriggerPtToRecoil)
				jetPt2=jetPt1;

			      hEDi->Fill(jetPt2);
			      hSpecDi->Fill(jetPt2);
			      hEtaPhiDi->Fill(mJ2->Eta(),mJ2->Phi());
			      hAreaPtDi->Fill(jetPt2,mJ2->Area());
			      
			      evPtDi=jetPt2;
			      evArea=mJ2->Area();		
			      
			      if (mptBkgDi!=0 && mptBkgDi->GetEntries()>0)		       
				{
				  DoZBkgFromPt(mZBkgDi,mptBkgDi,jetPt2);
				  mZBkgDi->Scale(mJ2->GetFFAreaInAcc());
				  hZBkgDi->Add(mZBkgDi);
				  
				  DoXiBkgFromPt(mXiBkgDi,mptBkg,jetPt2);
				  mXiBkgDi->Scale(mJ2->GetFFAreaInAcc());
				  hXiBkgDi->Add(mXiBkgDi);
				}
			      
			      FillFFHists(hXiDi,hZDi,mJ2,jetPt2);
			    }
			  else if (IsFakeJet(dphi))
			    {

			      hEFake->Fill(jetPt2);
			       
			      // evPtDi=jetPt2;
			      //evArea=mJ2->Area();
			      
			      if (mptBkgDi!=0 && mptBkgDi->GetEntries()>0)		       
				{
				  DoXiBkgFromPt(mXiBkgFake,mptBkgFake,jetPt2);
				  mXiBkgFake->Scale(mJ2->GetFFAreaInAcc());
				  hXiBkgFake->Add(mXiBkgFake);
				}

			      FillFFHists(hXiFake,0,mJ2,jetPt2);
			    }
			  else
			    {
			      evPtDi=-99.0;
			      evAreaDi=-1.0;
			    }			  
			}		 
		      else
			{
			  evPtDi=-99.0;
			  evAreaDi=-1.0;
			}
		    }
		  
		}
	      else
		{
		  if (verbose)
		    cout<<"No "<<JetAlgo<<" jet above pt,min = "<<minPtCut<<" found !"<<endl;
		  evPtDi=evPt=-99.0;
		  evAreaDi=evArea=-1.0;
		}
	      
	      delete mXiBkg;
	      delete mXiBkgDi;
	      delete mXiBkgFake;
	      delete mZBkg;
	      delete mZBkgDi;
	      delete mptBkg;
	      delete mptBkgDi;
	      delete mptBkgFake;
	    }
	  else
	    {

	      // *******************************************
	      //
	      // **** has to be updated !!!! ****
	      // **** not maintained at the moment !!! ****
	      //
	      // ================
	      // LOCone analysis 
	      // ================
	      //
	      // *******************************************
	   	      
	    }
	}
      else
	{
	  if (verbose)
	    cout<<"No "<<JetAlgo<<" jet found !"<<endl;

	  evPtDi=evPt=-99.0;
	  evAreaDi=evArea=-1.0;
	}
    }
  else
    {
      if (verbose) 
	cout<<"Warning : Zero pointer to ktMuEvent !"<<endl;

      evPtDi=evPt=-99.0;
      evAreaDi=evArea=-1.0;
    }
  
}
