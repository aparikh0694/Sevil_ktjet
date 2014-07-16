// first test (Joern Putschke)

#include "ktStarPico.h"
#include "TFile.h"
#include "TRandom.h"
#include<iostream> // needed for io
#include<sstream>  // needed for internal io
#include<vector> 
#include "TStarJetPicoUtils.h"

ClassImp(ktStarPico)

ktStarPico::ktStarPico(Bool_t mQAOut,Bool_t mSC) 
{

  NPrim=NTowers=RefMult=0;
  QAOutput=mQAOut;
  SC=mSC;
  nEvents=0;

  meanZdc=12781.7; // Mean of first 400k AuAu HT trigger Y07

  SCName="sc_test.root";
  QAName="qa_test.root";

  TriggerInfoArray = new TClonesArray("ktTriggerInfo",10);

  //ElectronCheck=true;
  info=false;

  // create QA histograms (fill them also here; compare to picoDst reader !)
  if (QAOutput)
    {
      hFitPoints=new TH1D("hFitPoints","hFitPoints",50,0,50);
      hFitOverMax=new TH1D("hFitOverMax","hFitOverMax",40,0,1.1);
      hDca=new TH1D("hDca","hDca",30,0,3);
      hEta=new TH1D("hEta","hEta",50,-1,1);
      hPhi=new TH1D("hPhi","hPhi",180,-TMath::Pi(),TMath::Pi());
      hPt=new TH1D("hPt","hPt",80,0,20);
      hFitPoints->SetDirectory(0);hFitOverMax->SetDirectory(0);hDca->SetDirectory(0);hEta->SetDirectory(0);hPhi->SetDirectory(0);hPt->SetDirectory(0);

      hMatchedFitPoints=new TH1D("hMatchedFitPoints","hMatchedFitPoints",50,0,50);
      hMatchedFitOverMax=new TH1D("hMatchedFitOverMax","hMatchedFitOverMax",40,0,1.1);
      hMatchedDca=new TH1D("hMatchedDca","hMatchedDca",30,0,3);
      hMatchedEta=new TH1D("hMatchedEta","hMatchedEta",50,-1,1);
      hMatchedPhi=new TH1D("hMatchedPhi","hMatchedPhi",180,-TMath::Pi(),TMath::Pi());
      hMatchedPt=new TH1D("hMatchedPt","hMatchedPt",80,0,20);
      hMatchedPoverE=new TH1D("hMatchedPoverE","hMatchedPoverE",80,0,2);
      hMatchedFitPoints->SetDirectory(0);hMatchedFitOverMax->SetDirectory(0);hMatchedDca->SetDirectory(0);hMatchedEta->SetDirectory(0);hMatchedPhi->SetDirectory(0);hMatchedPt->SetDirectory(0);hMatchedPoverE->SetDirectory(0);

      hTowADC=new TH1D("hTowADC","hTwoADC",250,0,250);
      hTowEta=new TH1D("hTowEta","hTowEta",50,-1,1);
      hTowPhi=new TH1D("hTowPhi","hTowPhi",180,-TMath::Pi(),TMath::Pi());
      hTowE=new TH1D("hTowE","hTowE",80,0,20);
      hTowEt=new TH1D("hTowEt","hTowEt",80,0,20);
      hTowADC->SetDirectory(0);hTowEta->SetDirectory(0);hTowPhi->SetDirectory(0);hTowE->SetDirectory(0);hTowEt->SetDirectory(0);

      hMatchedTowADC=new TH1D("hMatchedTowADC","hMatchedTwoADC",250,0,250);
      hMatchedTowEta=new TH1D("hMatchedTowEta","hMatchedTowEta",50,-1,1);
      hMatchedTowPhi=new TH1D("hMatchedTowPhi","hMatchedTowPhi",180,-TMath::Pi(),TMath::Pi());
      hMatchedTowE=new TH1D("hMatchedTowE","hMatchedTowE",80,0,20);
      hMatchedTowEt=new TH1D("hMatchedTowEt","hMatchedTowEt",80,0,20);
      hMatchedTowADC->SetDirectory(0);hMatchedTowEta->SetDirectory(0);hMatchedTowPhi->SetDirectory(0);hMatchedTowE->SetDirectory(0);hMatchedTowEt->SetDirectory(0);

      hTowIDEweighted=new TH1D("hTowIDEweighted","hTowIDEweighted",5000,0,5000);
      hTowIDEweightedEt=new TH1D("hTowIDEweightedEt","hTowIDEweightedEt",5000,0,5000);
      hTowIDEweighted->SetDirectory(0);hTowIDEweightedEt->SetDirectory(0);

      hVertexZ=new TH1D("hVertexZ","hVertexZ",100,-50,50);
      hRefMult=new TH1D("hRefMult","hRefMult",750,0,750);
      hTrigId=new TH1D(); // check range ...
      hVertexZ->SetDirectory(0);hRefMult->SetDirectory(0);hTrigId->SetDirectory(0);

      //hTowIDEweightedRunID=new TH2D("hTowIDEweightedRunID","hTowIDEweightedRunID",5000,0,5000,20000,0,20000);

      //mNtuple=new TNtuple("ht","Ntuple for HT per Run","RunId:TrigId:TowId:E:Et");
      hV0MassLa= new TH1D("hV0MassLa","Lambda Invariant Mass",100,1,1.3);
      hV0MassALa= new TH1D("hV0MassALa","Anti-Lambda Invariant Mass",100,1,1.3);
      hV0MassK0= new TH1D("hV0MassK0","K0s Invariant Mass",100,0.3,0.6);
      hV0MassK0->SetDirectory(0);hV0MassLa->SetDirectory(0);hV0MassALa->SetDirectory(0);
    }

  // space charge histograms ...
  if (SC)
    {
      hNevents=new TH1D("hNevents","Number of events",1,0,1);
      hEplusPt=new TH1D("hEplusPt","e+ pt spectrum",60,0,30);
      hEminusPt=new TH1D("hEminusPt","e- pt spectrum",60,0,30);
      hPlusPt=new TH1D("hPlusPt","h+ pt spectrum",60,0,30);
      hMinusPt=new TH1D("hMinusPt","h- pt spectrum",60,0,30);   
      hEplusPt->Sumw2();hEminusPt->Sumw2();hPlusPt->Sumw2();hMinusPt->Sumw2();
      
      hEplusPtW=new TH1D("hEplusPtW","e+ pt spectrum",60,0,30);
      hEminusPtW=new TH1D("hEminusPtW","e- pt spectrum",60,0,30);
      hPlusPtW=new TH1D("hPlusPtW","h+ pt spectrum",60,0,30);
      hMinusPtW=new TH1D("hMinusPtW","h- pt spectrum",60,0,30);
      hEplusPtW->Sumw2();hEminusPtW->Sumw2();hPlusPtW->Sumw2();hMinusPtW->Sumw2();

      ehEplusPtW=new TH1D("ehEplusPtW","e+ pt spectrum",60,0,30);
      ehEminusPtW=new TH1D("ehEminusPtW","e- pt spectrum",60,0,30);
      ehPlusPtW=new TH1D("ehPlusPtW","h+ pt spectrum",60,0,30);
      ehMinusPtW=new TH1D("ehMinusPtW","h- pt spectrum",60,0,30);
      ehEplusPtW->Sumw2();ehEminusPtW->Sumw2();ehPlusPtW->Sumw2();ehMinusPtW->Sumw2();
      
      whEplusPtW=new TH1D("whEplusPtW","e+ pt spectrum",60,0,30);
      whEminusPtW=new TH1D("whEminusPtW","e- pt spectrum",60,0,30);
      whPlusPtW=new TH1D("whPlusPtW","h+ pt spectrum",60,0,30);
      whMinusPtW=new TH1D("whMinusPtW","h- pt spectrum",60,0,30);
      whEplusPtW->Sumw2();whEminusPtW->Sumw2();whPlusPtW->Sumw2();whMinusPtW->Sumw2();

      // vs. normalized zdc coinc.
      zhEplusPtW=new TH2D("zhEplusPtW","e+ pt spectrum vs. norm. zdc coin",30,0,3,60,0,30);
      zhEminusPtW=new TH2D("zhEminusPtW","e- pt spectrum vs. norm. zdc coin",30,0,3,60,0,30);
      zhPlusPtW=new TH2D("zhPlusPtW","h+ pt spectrum vs. norm. zdc coin",30,0,3,60,0,30);
      zhMinusPtW=new TH2D("zhMinusPtW","h- pt spectrum vs. norm. zdc coin",30,0,3,60,0,30);
      zhEplusPtW->Sumw2();zhEminusPtW->Sumw2();zhPlusPtW->Sumw2();zhMinusPtW->Sumw2();

      zehEplusPtW=new TH2D("zehEplusPtW","e+ pt spectrum vs. norm. zdc coin",30,0,3,60,0,30);
      zehEminusPtW=new TH2D("zehEminusPtW","e- pt spectrum vs. norm. zdc coin",30,0,3,60,0,30);
      zehPlusPtW=new TH2D("zehPlusPtW","h+ pt spectrum vs. norm. zdc coin",30,0,3,60,0,30);
      zehMinusPtW=new TH2D("zehMinusPtW","h- pt spectrum vs. norm. zdc coin",30,0,3,60,0,30);
      zehEplusPtW->Sumw2();zehEminusPtW->Sumw2();zehPlusPtW->Sumw2();zehMinusPtW->Sumw2();
      
      zwhEplusPtW=new TH2D("zwhEplusPtW","e+ pt spectrum vs. norm. zdc coin",30,0,3,60,0,30);
      zwhEminusPtW=new TH2D("zwhEminusPtW","e- pt spectrum vs. norm. zdc coin",30,0,3,60,0,30);
      zwhPlusPtW=new TH2D("zwhPlusPtW","h+ pt spectrum vs. norm. zdc coin",30,0,3,60,0,30);
      zwhMinusPtW=new TH2D("zwhMinusPtW","h- pt spectrum vs. norm. zdc coin",30,0,3,60,0,30);
      zwhEplusPtW->Sumw2();zwhEminusPtW->Sumw2();zwhPlusPtW->Sumw2();zwhMinusPtW->Sumw2();

      // vs. zdc coinc.
      zdchEplusPtW=new TH2D("zdchEplusPtW","e+ pt spectrum vs. zdc coin",150,0,30000,60,0,30);
      zdchEminusPtW=new TH2D("zdchEminusPtW","e- pt spectrum vs. zdc coin",150,0,30000,60,0,30);
      zdchPlusPtW=new TH2D("zdchPlusPtW","h+ pt spectrum vs. zdc coin",150,0,30000,60,0,30);
      zdchMinusPtW=new TH2D("zdchMinusPtW","h- pt spectrum vs. zdc coin",150,0,30000,60,0,30);
      zdchEplusPtW->Sumw2();zdchEminusPtW->Sumw2();zdchPlusPtW->Sumw2();zdchMinusPtW->Sumw2();

      zdcehEplusPtW=new TH2D("zdcehEplusPtW","e+ pt spectrum vs. zdc coin",150,0,30000,60,0,30);
      zdcehEminusPtW=new TH2D("zdcehEminusPtW","e- pt spectrum vs. zdc coin",150,0,30000,60,0,30);
      zdcehPlusPtW=new TH2D("zdcehPlusPtW","h+ pt spectrum vs. zdc coin",150,0,30000,60,0,30);
      zdcehMinusPtW=new TH2D("zdcehMinusPtW","h- pt spectrum vs. zdc coin",150,0,30000,60,0,30);
      zdcehEplusPtW->Sumw2();zdcehEminusPtW->Sumw2();zdcehPlusPtW->Sumw2();zdcehMinusPtW->Sumw2();
      
      zdcwhEplusPtW=new TH2D("zdcwhEplusPtW","e+ pt spectrum vs. zdc coin",150,0,30000,60,0,30);
      zdcwhEminusPtW=new TH2D("zdcwhEminusPtW","e- pt spectrum vs. zdc coin",150,0,30000,60,0,30);
      zdcwhPlusPtW=new TH2D("zdcwhPlusPtW","h+ pt spectrum vs. zdc coin",150,0,30000,60,0,30);
      zdcwhMinusPtW=new TH2D("zdcwhMinusPtW","h- pt spectrum vs. zdc coin",150,0,30000,60,0,30);
      zdcwhEplusPtW->Sumw2();zdcwhEminusPtW->Sumw2();zdcwhPlusPtW->Sumw2();zdcwhMinusPtW->Sumw2();

      hZDCcoin=new TH1D("hZDCCoin","hZDCcoin",40000,0,40000);
    }

  //DEBUG:
  //cout<<endl;
  //cout<<"Default constructor of ktStarPico ..."<<endl;
}

void ktStarPico::WriteHistograms()
{
  cout<<endl;
  cout<<"Write histograms to file ..."<<endl;
  cout<<endl;
  
  if (QAOutput)
    {
      cout<<" QA file = "<<QAName<<endl;
    }
  
  if (SC)
    {
      cout<<" SC file = "<<SCName<<endl;
      TFile *fsc=new TFile(SCName,"RECREATE");
      
      hNevents->SetBinContent(1,nEvents);

      // scale by bin width
      hEplusPt->Scale(1/(Double_t) hEplusPt->GetBinWidth(1));
      hEminusPt->Scale(1/(Double_t) hEminusPt->GetBinWidth(1));
      hPlusPt->Scale(1/(Double_t) hPlusPt->GetBinWidth(1));
      hMinusPt->Scale(1/(Double_t) hMinusPt->GetBinWidth(1));

      hEplusPtW->Scale(1/(Double_t) hEplusPtW->GetBinWidth(1));
      hEminusPtW->Scale(1/(Double_t) hEminusPtW->GetBinWidth(1));
      hPlusPtW->Scale(1/(Double_t) hPlusPtW->GetBinWidth(1));
      hMinusPtW->Scale(1/(Double_t) hMinusPtW->GetBinWidth(1));

      ehEplusPtW->Scale(1/(Double_t) ehEplusPtW->GetBinWidth(1));
      ehEminusPtW->Scale(1/(Double_t) ehEminusPtW->GetBinWidth(1));
      ehPlusPtW->Scale(1/(Double_t) ehPlusPtW->GetBinWidth(1));
      ehMinusPtW->Scale(1/(Double_t) ehMinusPtW->GetBinWidth(1));

      whEplusPtW->Scale(1/(Double_t) whEplusPtW->GetBinWidth(1));
      whEminusPtW->Scale(1/(Double_t) whEminusPtW->GetBinWidth(1));
      whPlusPtW->Scale(1/(Double_t) whPlusPtW->GetBinWidth(1));
      whMinusPtW->Scale(1/(Double_t) whMinusPtW->GetBinWidth(1));
      
      // mormalize 2D histos ...

      // scale by number of events (not necessary for ratio)
      /*
      hEplusPt->Scale(1/(Double_t) nEvents);
      hEminusPt->Scale(1/(Double_t) nEvents);
      hPlusPt->Scale(1/(Double_t) nEvents);
      hMinusPt->Scale(1/(Double_t) nEvents);
      */

      hEplusPtW->Scale(1/(Double_t) nEvents);
      hEminusPtW->Scale(1/(Double_t) nEvents);
      hPlusPtW->Scale(1/(Double_t) nEvents);
      hMinusPtW->Scale(1/(Double_t) nEvents);

      ehEplusPtW->Scale(1/(Double_t) nEvents);
      ehEminusPtW->Scale(1/(Double_t) nEvents);
      ehPlusPtW->Scale(1/(Double_t) nEvents);
      ehMinusPtW->Scale(1/(Double_t) nEvents);

      whEplusPtW->Scale(1/(Double_t) nEvents);
      whEminusPtW->Scale(1/(Double_t) nEvents);
      whPlusPtW->Scale(1/(Double_t) nEvents);
      whMinusPtW->Scale(1/(Double_t) nEvents);

      zhEplusPtW->Scale(1/(Double_t) nEvents);
      zhEminusPtW->Scale(1/(Double_t) nEvents);
      zhPlusPtW->Scale(1/(Double_t) nEvents);
      zhMinusPtW->Scale(1/(Double_t) nEvents);

      zehEplusPtW->Scale(1/(Double_t) nEvents);
      zehEminusPtW->Scale(1/(Double_t) nEvents);
      zehPlusPtW->Scale(1/(Double_t) nEvents);
      zehMinusPtW->Scale(1/(Double_t) nEvents);

      zwhEplusPtW->Scale(1/(Double_t) nEvents);
      zwhEminusPtW->Scale(1/(Double_t) nEvents);
      zwhPlusPtW->Scale(1/(Double_t) nEvents);
      zwhMinusPtW->Scale(1/(Double_t) nEvents);

      zdchEplusPtW->Scale(1/(Double_t) nEvents);
      zdchEminusPtW->Scale(1/(Double_t) nEvents);
      zdchPlusPtW->Scale(1/(Double_t) nEvents);
      zdchMinusPtW->Scale(1/(Double_t) nEvents);

      zdcehEplusPtW->Scale(1/(Double_t) nEvents);
      zdcehEminusPtW->Scale(1/(Double_t) nEvents);
      zdcehPlusPtW->Scale(1/(Double_t) nEvents);
      zdcehMinusPtW->Scale(1/(Double_t) nEvents);

      zdcwhEplusPtW->Scale(1/(Double_t) nEvents);
      zdcwhEminusPtW->Scale(1/(Double_t) nEvents);
      zdcwhPlusPtW->Scale(1/(Double_t) nEvents);
      zdcwhMinusPtW->Scale(1/(Double_t) nEvents);
      

      hNevents->Write("hNevents");
      hEplusPt->Write("hEplusPt");
      hEminusPt->Write("hEminusPt");
      hPlusPt->Write("hPlusPt");
      hMinusPt->Write("hMinusPt");
      
      hEplusPtW->Write("hEplusPtW");
      hEminusPtW->Write("hEminusPtW");
      hPlusPtW->Write("hPlusPtW");
      hMinusPtW->Write("hMinusPtW");

      ehEplusPtW->Write("ehEplusPtW");
      ehEminusPtW->Write("ehEminusPtW");
      ehPlusPtW->Write("ehPlusPtW");
      ehMinusPtW->Write("ehMinusPtW");

      whEplusPtW->Write("whEplusPtW");
      whEminusPtW->Write("whEminusPtW");
      whPlusPtW->Write("whPlusPtW");
      whMinusPtW->Write("whMinusPtW");

      zhEplusPtW->Write("zhEplusPtW");
      zhEminusPtW->Write("zhEminusPtW");
      zhPlusPtW->Write("zhPlusPtW");
      zhMinusPtW->Write("zhMinusPtW");

      zehEplusPtW->Write("zehEplusPtW");
      zehEminusPtW->Write("zehEminusPtW");
      zehPlusPtW->Write("zehPlusPtW");
      zehMinusPtW->Write("zehMinusPtW");

      zwhEplusPtW->Write("zwhEplusPtW");
      zwhEminusPtW->Write("zwhEminusPtW");
      zwhPlusPtW->Write("zwhPlusPtW");
      zwhMinusPtW->Write("zwhMinusPtW");

      zdchEplusPtW->Write("zdchEplusPtW");
      zdchEminusPtW->Write("zdchEminusPtW");
      zdchPlusPtW->Write("zdchPlusPtW");
      zdchMinusPtW->Write("zdchMinusPtW");

      zdcehEplusPtW->Write("zdcehEplusPtW");
      zdcehEminusPtW->Write("zdcehEminusPtW");
      zdcehPlusPtW->Write("zdcehPlusPtW");
      zdcehMinusPtW->Write("zdcehMinusPtW");

      zdcwhEplusPtW->Write("zdcwhEplusPtW");
      zdcwhEminusPtW->Write("zdcwhEminusPtW");
      zdcwhPlusPtW->Write("zdcwhPlusPtW");
      zdcwhMinusPtW->Write("zdcwhMinusPtW");

      hZDCcoin->Write("hZDCcoin");

      fsc->Write();
      fsc->Close();
    }

  cout<<endl;
}

ktStarPico::~ktStarPico()
{

  if (QAOutput || SC)
    WriteHistograms();

  // DEBUG:
  //cout<<"Delete histograms ..."<<endl;

  if (QAOutput)
    {
      delete hFitPoints; delete hFitOverMax; delete hDca;
      delete hEta; delete hPhi; delete hPt;
      delete hMatchedFitPoints; delete hMatchedFitOverMax; delete hMatchedDca;
      delete hMatchedEta; delete hMatchedPhi;delete hMatchedPoverE;
      
      delete hTowADC; delete hTowEta; delete hTowPhi; delete hTowE; delete hTowEt;
      delete hMatchedTowADC; delete hMatchedTowEta; delete hMatchedTowPhi; delete hMatchedTowE; delete hMatchedTowEt;
      
      delete hTowIDEweighted;
      delete hTowIDEweightedEt;
      
      delete hVertexZ; delete hRefMult; delete hTrigId;
      
      delete hV0MassLa; delete hV0MassALa; delete hV0MassK0;
    }
  
  if (SC)
    {
      delete hEplusPt;delete hEminusPt; delete hPlusPt; delete hMinusPt;
      delete hEplusPtW;delete hEminusPtW; delete hPlusPtW; delete hMinusPtW;
      delete ehEplusPtW;delete ehEminusPtW; delete ehPlusPtW; delete ehMinusPtW;
      delete whEplusPtW;delete whEminusPtW; delete whPlusPtW; delete whMinusPtW;
      delete zhEplusPtW;delete zhEminusPtW; delete zhPlusPtW; delete zhMinusPtW;
      delete zehEplusPtW;delete zehEminusPtW; delete zehPlusPtW; delete zehMinusPtW;
      delete zwhEplusPtW;delete zwhEminusPtW; delete zwhPlusPtW; delete zwhMinusPtW;
      
      delete zdchEplusPtW;delete zdchEminusPtW; delete zdchPlusPtW; delete zdchMinusPtW;
      delete zdcehEplusPtW;delete zdcehEminusPtW; delete zdcehPlusPtW; delete zdcehMinusPtW;
      delete zdcwhEplusPtW;delete zdcwhEminusPtW; delete zdcwhPlusPtW; delete zdcwhMinusPtW;

      delete hNevents;
      delete hZDCcoin;
    }

  //DEBUG:
  //cout<<"TriggerInfo Array ..."<<endl;

  TriggerInfoArray->Delete();//Clear("C");
  delete TriggerInfoArray;

  //DEBUG:
  //cout<<endl;
  //cout<<"Destructor of ktStarPico ..."<<endl;
}

void ktStarPico::Clear()
{
  TriggerInfoArray->Clear("C");
}

void ktStarPico::PrintEventInfo()
{
  cout<<"Event Info:"<<endl;
  cout<<"# primaries = "<<NPrim<<endl;
  cout<<"# towers    = "<<NTowers<<endl;
  cout<<"RefMult     = "<<RefMult<<endl;
}

void ktStarPico::PrintCuts()
{
  PrintQACuts();
}

void ktStarPico::PrintQACuts()
{
  // dummy (see TStarPico ....)
  // but implement output also khere ! (to be done)
}

void ktStarPico::WriteQAOutputFile(TString fName)
{
  TFile *f=new TFile(fName,"RECREATE");

  hFitPoints->Write("hFitPoints");
  hFitOverMax->Write("hFitOverMax");
  hDca->Write("hDca");
  hEta->Write("hEta");
  hPhi->Write("hPhi");
  hPt->Write("hPt");

  hMatchedFitPoints->Write("hMatchedFitPoints");
  hMatchedFitOverMax->Write("hMatchedFitOverMax");
  hMatchedDca->Write("hMatchedDca");
  hMatchedEta->Write("hMatchedEta");
  hMatchedPhi->Write("hMatchedPhi");
  hMatchedPoverE->Write("hMatchedPoverE");
  hMatchedPt->Write("hMatchedPt");
 
  hTowADC->Write("hTowADC");
  hTowEta->Write("hTowEta");
  hTowPhi->Write("hTowPhi");
  hTowE->Write("hTowE");
  hTowEt->Write("hTowEt");
  
  hMatchedTowADC->Write("hMatchedTowADC");
  hMatchedTowEta->Write("hMatchedTowEta");
  hMatchedTowPhi->Write("hMatchedTowPhi");
  hMatchedTowE->Write("hMatchedTowE");
  hMatchedTowEt->Write("hMatchedTowEt");
  
  hTowIDEweighted->Write("hTowIDEweighted");
  hTowIDEweightedEt->Write("hTowIDEweightedEt");

  //hTowIDEweightedRunID->Write("hTowIDEweightedRunID");

  hVertexZ->Write("hVertexZ");
  hRefMult->Write("hRefMult");
  hTrigId->Write("hTrigId");

  hV0MassLa->Write("hV0MassLa");
  hV0MassALa->Write("hV0MassALa");
  hV0MassK0->Write("hV0MassK0");


  //mNtuple->Write("ht");

  f->Write();
  f->Close();

  //cout<<endl;
  cout<<"QA output saved in file = "<<fName<<"\t #Events = "<<hVertexZ->GetEntries()<<endl;
  //cout<<endl;

}

void  ktStarPico::SetTriggerInfo(TClonesArray* trinfo){

  //cout<<" ------  ktStarPico::SetTriggerInfo --------"<<endl;
  TClonesArray &trigobj = *TriggerInfoArray;
  for( int i=0; i<trinfo->GetEntriesFast(); i++)
    {
      TStarJetPicoTriggerInfo *t = (TStarJetPicoTriggerInfo *)trinfo->At(i);
      
      ktTriggerInfo *trig = new(trigobj[i]) ktTriggerInfo();
      trig->SetPhi(t->GetPhi());
      trig->SetEta(t->GetEta());
      trig->SetTriggerFlag(t->GetTriggerFlag());
      //trig->PrintInfo();
  }

  //cout<<" N. of trig obj in the event = "<<TriggerInfoArray->GetEntriesFast()<<endl;
}

void ktStarPico::PrintTriggerInfo()
{
  //DEBUG:
  cout<<" ------  ktStarPico::PrintTriggerInfo --------"<<endl;
  cout<<" N. of trig obj in the event = "<<TriggerInfoArray->GetEntriesFast()<<endl;
  for( int i=0; i<TriggerInfoArray->GetEntriesFast(); i++)
    {
      TStarJetPicoTriggerInfo *t = (TStarJetPicoTriggerInfo *)TriggerInfoArray->At(i);
      t->PrintInfo();
    }
}

void ktStarPico::DoSpaceCharge(TLorentzVector *p,ktPID *pid,Float_t ZdcCoin)
{
  // DEBUG:
  /*
   if (info)
    {
      cout<<endl;
      cout<<"Fill space-charge plots ..."<<endl;
      cout<<endl;
    }
  */

  Double_t normZdc=ZdcCoin/meanZdc;

  hZDCcoin->Fill(ZdcCoin);

   // fill spectra charged separated for h+- and e+-
   // (save also tree for further cuts only charged !?)

   if (pid->IsCharged())
     {
       if (pid->GetCharge()>0)
	 {
	   hPlusPt->Fill(p->Pt());
	   hPlusPtW->Fill(p->Pt(),1/(Double_t) p->Pt());
	   zhPlusPtW->Fill(normZdc,p->Pt(),1/(Double_t) p->Pt());
	   zdchPlusPtW->Fill(ZdcCoin,p->Pt(),1/(Double_t) p->Pt());

	   if (p->Eta()>0)
	     {
	       whPlusPtW->Fill(p->Pt(),1/(Double_t) p->Pt());
	       zwhPlusPtW->Fill(normZdc,p->Pt(),1/(Double_t) p->Pt());
	       zdcwhPlusPtW->Fill(ZdcCoin,p->Pt(),1/(Double_t) p->Pt());
	     }
	   else
	     {
	       ehPlusPtW->Fill(p->Pt(),1/(Double_t) p->Pt());
	       zehPlusPtW->Fill(normZdc,p->Pt(),1/(Double_t) p->Pt());
	       zdcehPlusPtW->Fill(ZdcCoin,p->Pt(),1/(Double_t) p->Pt());
	     }
	  
	 }
       else
	 {
	   hMinusPt->Fill(p->Pt());
	   hMinusPtW->Fill(p->Pt(),1/(Double_t) p->Pt());
	   zhMinusPtW->Fill(normZdc,p->Pt(),1/(Double_t) p->Pt());
	   zdchMinusPtW->Fill(ZdcCoin,p->Pt(),1/(Double_t) p->Pt()); 

	   if (p->Eta()>0)
	     {
	       whMinusPtW->Fill(p->Pt(),1/(Double_t) p->Pt());
	       zwhMinusPtW->Fill(normZdc,p->Pt(),1/(Double_t) p->Pt());
	       zdcwhMinusPtW->Fill(ZdcCoin,p->Pt(),1/(Double_t) p->Pt());
	     }
	   else
	     {
	       ehMinusPtW->Fill(p->Pt(),1/(Double_t) p->Pt());
	       zehMinusPtW->Fill(normZdc,p->Pt(),1/(Double_t) p->Pt());
	       zdcehMinusPtW->Fill(ZdcCoin,p->Pt(),1/(Double_t) p->Pt());
	     }
	 }

       if (pid->GetPID()==11)
	 {
	   // is electron (maybe include sperate cuts at this leevel)
	   if (pid->GetCharge()>0)
	     {
	       hEplusPt->Fill(p->Pt());
	       hEplusPtW->Fill(p->Pt(),1/(Double_t) p->Pt());
	       zhEplusPtW->Fill(normZdc,p->Pt(),1/(Double_t) p->Pt());
	       zdchEplusPtW->Fill(ZdcCoin,p->Pt(),1/(Double_t) p->Pt());
	      
	       if (p->Eta()>0)
		 {
		   whEplusPtW->Fill(p->Pt(),1/(Double_t) p->Pt()); 
		   zwhEplusPtW->Fill(normZdc,p->Pt(),1/(Double_t) p->Pt()); 
		   zdcwhEplusPtW->Fill(ZdcCoin,p->Pt(),1/(Double_t) p->Pt()); 
		 }
	       else
		 {
		   ehEplusPtW->Fill(p->Pt(),1/(Double_t) p->Pt());
		   zehEplusPtW->Fill(normZdc,p->Pt(),1/(Double_t) p->Pt());
		   zdcehEplusPtW->Fill(ZdcCoin,p->Pt(),1/(Double_t) p->Pt());
		 }
	     }
	   else
	     {
	       hEminusPt->Fill(p->Pt());
	       hEminusPtW->Fill(p->Pt(),1/(Double_t) p->Pt()); 	 
	       zhEminusPtW->Fill(normZdc,p->Pt(),1/(Double_t) p->Pt()); 
	       zdchEminusPtW->Fill(ZdcCoin,p->Pt(),1/(Double_t) p->Pt()); 

	       if (p->Eta()>0)
		 {
		   whEminusPtW->Fill(p->Pt(),1/(Double_t) p->Pt()); 
		   zwhEminusPtW->Fill(normZdc,p->Pt(),1/(Double_t) p->Pt());  
		   zdcwhEminusPtW->Fill(ZdcCoin,p->Pt(),1/(Double_t) p->Pt());  
		 }
	       else
		 {
		   ehEminusPtW->Fill(p->Pt(),1/(Double_t) p->Pt()); 	   
		   zehEminusPtW->Fill(normZdc,p->Pt(),1/(Double_t) p->Pt());
		   zdcehEminusPtW->Fill(ZdcCoin,p->Pt(),1/(Double_t) p->Pt());
		 }
	     }
	 }
     }
}

Bool_t ktStarPico::Fill(ktGrid *mGrid,ktFastJet *mFast,TStarJetPicoEvent *mEv,TStarJetVectorContainer<TStarJetVector>* container)
{
  //DEBUG:
  if (info)
    {
      cout<<endl;
      cout<<"Fill ktGrid/ktFastJet from STAR jetPico Dsts: Event# = "<<mEv->GetHeader()->GetEventId()<<" Run# = "<<mEv->GetHeader()->GetRunId()<<endl;
      cout<<"TrigID = "<<TStarJetPicoUtils::GetTriggerIdsString(mEv)<<endl;
    }
  
  // Fill general event infos ...
  NPrim=mEv->GetHeader()->GetNOfPrimaryTracks();
  NTowers=mEv->GetHeader()->GetNOfTowers();
  RefMult=mEv->GetHeader()->GetReferenceMultiplicity();

  nEvents++;
  
  // Fill trigger objects ...
  SetTriggerInfo(mEv->GetTrigObjs());

  //DEBUG:
  //cout<<mEv->GetTrigObjs()->GetEntries()<<endl;

  //DEBUG:
  // PrintEventInfo();

  // loop over container including all tower/track cuts + corrections (check this ...)
  for (int i=0;i<container->GetEntries();i++)
    {
      //get output vector ...
      TStarJetVector *v=(TStarJetVector*) container->Get(i);
      TLorentzVector *pv=new TLorentzVector();
      ktPID *mPID=new ktPID();

      // prepare lorentzvecotor ...
      Double_t mPt=sqrt(v->px()*v->px()+v->py()*v->py());
      pv->SetPtEtaPhiM(mPt,v->pseudorapidity(),v->phi_std(),0);

      // fill PID infos ...
      mPID->SetCharge(v->GetCharge());
      mPID->SetPID(v->GetPID());
      mPID->SetdEdx(v->GetFeatureD(5));
      mPID->SetMass(v->GetFeatureD(1));
      mPID->SetNSigmaPion(v->GetFeatureD(2));
      mPID->SetNSigmaKaon(v->GetFeatureD(3));
      mPID->SetNSigmaProton(v->GetFeatureD(4));
      if( v->GetType() & TStarJetVector::_MATCHED) // check the bit is set
	mPID->SetETower(v->P()/v->GetPoverE());
      
      if( v->GetType() == TStarJetVector::_TOWER)
	mPID->SetETower(v->P());
      
      //DEBUG:
      if( fabs(mPID->GetPID()) == 3122 || mPID->GetPID() == 310){
	mPID->PrintPID();
      }
      //cout<<v->GetCharge()<<" "<<mPID->GetCharge()<<endl;
      //cout<<mPID->GetMass()<<endl;
      //if (mPID->GetPID()!=0)
      //cout<<mPID->GetPID()<<" "<<mPID->GetCharge()<<endl;

      if (mFast!=0)
	{
	  mFast->AddParticle(pv->Px(),pv->Py(),pv->Pz(),pv->E(),mPID);
	}	
	      
      if ( mGrid!=0) 
	{
	  mGrid->Fill((TLorentzVector*) pv->Clone(),mPID);
	}
      
      // Fill QA histograms if needed ... (to be done)

      // Do space charge plots ...
      if (SC)
	DoSpaceCharge(pv,mPID,mEv->GetHeader()->GetZdcCoincidenceRate());

      delete pv;
      delete mPID;
    }
  
}
