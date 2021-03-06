#include "ktMuEvent.h"
//#include "ktAreaUtil.h"
#include "ktParticle.h"
#include "ktFastJet.h"

#include <Riostream.h>

ClassImp(ktMuEvent)

//TClonesArray *ktMuEvent::fgTrigObjs   = 0;


ktMuEvent::ktMuEvent()
{

  ktList=new TObjArray(0);
  ktBFList=new TObjArray(0);
  coneList=new TObjArray(0);
  allList=new TObjArray(0);
  PartonList=new TObjArray(0);

  fjKt=new TObjArray(0);
  fjAntiKt=new TObjArray(0);
  fjSISCone=new TObjArray(0);

  ueEvents = new TObjArray(0);

  //if (fgTrigObjs == 0){
  //fgTrigObjs  = new TClonesArray("ktTriggerInfo",10);
  //}
  //TriggerInfoArray=new TClonesArray("ktTriggerInfo",10);//fgTrigObjs;
  
  TriggerInfoArray=new TObjArray(0);

  evNum=-1;
  BkgPtPerCell=0.0;
  BkgPtPerCellAll=0.0;
  BkgPtPerCellIn=0.0;
  BkgPtPerCellAllIn=0.0;
  BkgPtPerCellOut=0.0;
  BkgPtPerCellAllOut=0.0;
  BkgPtCone=0.0;
  BkgPtCluster=0.0;
  BkgNCellCone=0.0;
  BkgNCellCluster=0.0;
  NCellsInJetArea=0.0;

  refMult=0;
  vertexZ=-999;
  PsiRP=-999;
  TrigId=-99;

  BkgNCellIn=0.0;
  BkgNCellOut=0.0;
  BkgNJetCellIn=0.0;
  BkgNJetCellOut=0.0;
  BkgNCellInAll=0.0;
  BkgNCellOutAll=0.0;

  jfPt=0.0;
  jfRc=0.0;
  jfRcBkg=0.0;
  jfEmcalSet=false;
  jfPid="";

  jfNSubJets=0;

  median_pt_per_area=0.0;

  // DEBUG:
  //cout<<"Default ktMuEvent constructor "<<endl;

  hXiBkg=0;
  hXiBkgKt=0;
  hXiBkgAntiKt=0;
  hXiBkgSISCone=0;

  hZBkg=0;
  hZBkgKt=0;
  hZBkgAntiKt=0;
  hZBkgSISCone=0;

  hptBkg=0;
  hptBkgKt=0;
  hptBkgAntiKt=0;
  hptBkgSISCone=0;
}


ktMuEvent::ktMuEvent(Int_t mEvNum)
{

  ktList=new TObjArray(0);
  ktBFList=new TObjArray(0);
  coneList=new TObjArray(0);
  allList=new TObjArray(0);
  PartonList=new TObjArray(0);
  
  fjKt=new TObjArray(0);
  fjAntiKt=new TObjArray(0);
  fjSISCone=new TObjArray(0);
  
  ueEvents = new TObjArray(0);
 
  //if (fgTrigObjs == 0){
  //fgTrigObjs  = new TClonesArray("ktTriggerInfo",10);
  //}
  //TriggerInfoArray=new TClonesArray("ktTriggerInfo",10);//fgTrigObjs;
  
  TriggerInfoArray=new TObjArray(0);

  evNum=mEvNum;
  BkgPtPerCell=0.0;
  BkgPtPerCellAll=0.0;
  BkgPtPerCellIn=0.0;
  BkgPtPerCellAllIn=0.0;
  BkgPtPerCellOut=0.0;
  BkgPtPerCellAllOut=0.0;
  BkgPtCone=0.0;
  BkgPtCluster=0.0;
  BkgNCellCone=0.0;
  BkgNCellCluster=0.0;
  NCellsInJetArea=0.0;

  refMult=0;
  vertexZ=-999;
  PsiRP=-999;
  TrigId=-99;

  BkgNCellIn=0.0;
  BkgNCellOut=0.0;
  BkgNJetCellIn=0.0;
  BkgNJetCellOut=0.0;
  BkgNCellInAll=0.0;
  BkgNCellOutAll=0.0;

  jfPt=0.0;
  jfRc=0.0;
  jfRcBkg=0.0;
  jfEmcalSet=false;
  jfPid="";

  jfNSubJets=0;

  hXiBkg=0;
  hXiBkgKt=0;
  hXiBkgAntiKt=0;
  hXiBkgSISCone=0;

  hZBkg=0;
  hZBkgKt=0;
  hZBkgAntiKt=0;
  hZBkgSISCone=0;

  hptBkg=0;
  hptBkgKt=0;
  hptBkgAntiKt=0;
  hptBkgSISCone=0;

  median_pt_per_area=0.0;

  // DEBUG:
  //cout<<"Default ktMuEvent + evNum constructor "<<endl;
}

ktMuEvent::~ktMuEvent()
{

  //TriggerInfoArray->Clear("C");
  TriggerInfoArray->Delete();
  //cout<<TriggerInfoArray<<endl;
  delete TriggerInfoArray;

  ktList->Delete();
  delete ktList;
 
  ktBFList->Delete();
  delete ktBFList;
 
  coneList->Delete();
  delete coneList;
 
  allList->Delete();
  delete allList;
 
  PartonList->Delete();
  delete PartonList;

  fjKt->Delete();
  delete fjKt;

  fjAntiKt->Delete();
  delete fjAntiKt;
  
  fjSISCone->Delete();
  delete fjSISCone;

  //cout<<"Hmm :"<<ueEvents->GetEntries()<<endl;
  ueEvents->Delete();
  //cout<<"Aha :"<<ueEvents->GetEntries()<<endl;
  delete ueEvents;

  delete hXiBkg;
  delete hXiBkgKt;
  delete hXiBkgAntiKt;
  delete hXiBkgSISCone;

  delete hZBkg;
  delete hZBkgKt;
  delete hZBkgAntiKt;
  delete hZBkgSISCone;

  
  delete hptBkg;
  delete hptBkgKt;
  delete hptBkgAntiKt;
  delete hptBkgSISCone;
  // DEBUG:
  //cout<<"Default ktMuEvent destructor"<<endl;
}


void ktMuEvent::Fill(TObjArray *mJets)
{

   for (int i=0;i<mJets->GetEntries();i++)
    {
      
      //ktMuJet *mJet=new ktMuJet();
      ktJet *mktJet=(ktJet*) (mJets->At(i));
      //cout<<"Before ktMuJet ..."<<endl;
      ktMuJet *mJet=new ktMuJet(mktJet);
      //cout<<"After ktMuJet ..."<<endl;

       if (mJet->GetType().Contains("Cone"))
	 //FillConeBkgJet(myOutJet);
	 coneList->AddLast(mJet);
	 //coneList->AddLast(mJets->At(i));
       else if (mJet->GetType().Contains("KtBF"))
	 //FillKtBFJet(myOutJet);
	 ktBFList->AddLast(mJet);
	 //ktBFList->AddLast(mJets->At(i));
       else
	 //FillKtJet(myOutJet);
	 ktList->AddLast(mJet);
	 //ktList->AddLast(mJets->At(i));

       //delete mJet;
    }
   
   // DEBUG:
   //cout<< GetNKtJets() << " "<< GetNConeJets() << " "<< GetNKtBFJets() <<endl;

}

void ktMuEvent::Fill(TObjArray *mJets,Double_t cellIn,Double_t cellOut)
{

   for (int i=0;i<mJets->GetEntries();i++)
    {
      
      //ktMuJet *mJet=new ktMuJet();
      ktJet *mktJet=(ktJet*) (mJets->At(i));
      //cout<<"Before ktMuJet ..."<<endl;
      ktMuJet *mJet=new ktMuJet(mktJet);
      mJet->SetPtCellCorr(cellIn,cellOut);

      //cout<<"After ktMuJet ..."<<endl;

       if (mJet->GetType().Contains("Cone"))
	 //FillConeBkgJet(myOutJet);
	 coneList->AddLast(mJet);
	 //coneList->AddLast(mJets->At(i));
       else if (mJet->GetType().Contains("KtBF"))
	 //FillKtBFJet(myOutJet);
	 ktBFList->AddLast(mJet);
	 //ktBFList->AddLast(mJets->At(i));
       else
	 //FillKtJet(myOutJet);
	 ktList->AddLast(mJet);
	 //ktList->AddLast(mJets->At(i));

       //delete mJet;
    }

   if (coneList->GetEntries()>1)
     coneList->Sort();
   
   // DEBUG:
   //cout<< GetNKtJets() << " "<< GetNConeJets() << " "<< GetNKtBFJets() <<endl;

}

void ktMuEvent::AddParton(Double_t px,Double_t py,Double_t pz, Double_t e)
{
  ktParton *mParton=new ktParton(px,py,pz,e);
  PartonList->AddLast(mParton);

  // DEBUG:
  //cout<<PartonList->GetEntries()<<endl;
  //mParton->PrintParton();

}

void ktMuEvent::Fill(TObjArray *mJets,Double_t cellIn,Double_t cellOut, Double_t coneIn)
{

   for (int i=0;i<mJets->GetEntries();i++)
    {
      
      //ktMuJet *mJet=new ktMuJet();
      ktJet *mktJet=(ktJet*) (mJets->At(i));
      //cout<<"Before ktMuJet ..."<<endl;
      ktMuJet *mJet=new ktMuJet(mktJet);
      mJet->SetPtCellCorr(cellIn,cellOut);
      mJet->SetPtConeCorr(coneIn);

      jfNSubJets=mJet->GetNSubJets();

      //cout<<"After ktMuJet ..."<<endl;

       if (mJet->GetType().Contains("Cone"))
	 //FillConeBkgJet(myOutJet);
	 coneList->AddLast(mJet);
	 //coneList->AddLast(mJets->At(i));
       else if (mJet->GetType().Contains("KtBF"))
	 //FillKtBFJet(myOutJet);
	 ktBFList->AddLast(mJet);
	 //ktBFList->AddLast(mJets->At(i));
       else
	 //FillKtJet(myOutJet);
	 ktList->AddLast(mJet);
	 //ktList->AddLast(mJets->At(i));

       //delete mJet;
    }

   // Think if sorted or not ...
   /*
   if (coneList->GetEntries()>1)
     coneList->Sort();
   */

   // DEBUG:
   //cout<< GetNKtJets() << " "<< GetNConeJets() << " "<< GetNKtBFJets() <<endl;

}

void ktMuEvent::Fill(TObjArray *mJets,Double_t cellIn,Double_t cellOut, Double_t coneIn, Double_t clusterIn)
{

   for (int i=0;i<mJets->GetEntries();i++)
    {
      
      //ktMuJet *mJet=new ktMuJet();
      ktJet *mktJet=(ktJet*) (mJets->At(i));
      //cout<<"Before ktMuJet ..."<<endl;
      ktMuJet *mJet=new ktMuJet(mktJet);
      mJet->SetPtCellCorr(cellIn,cellOut);
      mJet->SetPtConeCorr(coneIn);

      if (mJet->GetType().Contains("Cluster"))
	mJet->SetPtClusterCorr(clusterIn);

      jfNSubJets=mJet->GetNSubJets();

      //cout<<"After ktMuJet ..."<<endl;

       if (mJet->GetType().Contains("Cone"))
	 //FillConeBkgJet(myOutJet);
	 coneList->AddLast(mJet);
	 //coneList->AddLast(mJets->At(i));
       else if (mJet->GetType().Contains("KtBF"))
	 //FillKtBFJet(myOutJet);
	 ktBFList->AddLast(mJet);
	 //ktBFList->AddLast(mJets->At(i));
       else
	 //FillKtJet(myOutJet);
	 ktList->AddLast(mJet);
	 //ktList->AddLast(mJets->At(i));

       //delete mJet;
    }

   // Think if sorted or not ...
   /*
   if (coneList->GetEntries()>1)
     coneList->Sort();
   */

   // DEBUG:
   //cout<< GetNKtJets() << " "<< GetNConeJets() << " "<< GetNKtBFJets() <<endl;

}

void ktMuEvent::FillBkgInfo(ktGrid *mGrid)
{
  SetBkgPtPerCell(mGrid->GetBkgPtPerCell());
  SetBkgPtPerCellAll(mGrid->GetBkgPtPerCellAll());
  SetBkgPtPerCellIn(mGrid->GetBkgPtPerCellIn());
  SetBkgPtPerCellAllIn(mGrid->GetBkgPtPerCellAllIn());
  SetBkgPtPerCellOut(mGrid->GetBkgPtPerCellOut());
  SetBkgPtPerCellAllOut(mGrid->GetBkgPtPerCellAllOut());
  SetBkgPtCone(mGrid->GetBkgPtCone());
  SetBkgPtCluster(mGrid->GetBkgPtCluster());
  SetBkgPtCluster(mGrid-> GetBkgPtCluster());
  SetBkgNCellCone(mGrid->GetBkgNCellCone());
  SetBkgNCellCluster(mGrid->GetBkgNCellCluster());
  SetBkgNCellIn(mGrid->GetBkgNCellIn());
  SetBkgNCellOut(mGrid->GetBkgNCellOut());
  SetBkgNJetCellIn(mGrid->GetBkgNJetCellIn());
  SetBkgNJetCellOut(mGrid->GetBkgNJetCellOut());
  SetBkgNCellInAll(mGrid->GetBkgNCellInAll());
  SetBkgNCellOutAll(mGrid->GetBkgNCellOutAll());
  SetNCellsInJetArea(mGrid->GetNCellsInJet(mGrid->GetRc()));
  hXiBkg=(TH1D*) mGrid->hXiBkg->Clone();
  
}

void ktMuEvent::FillJetFinderSettings(ktGrid *mGrid)
{
  
  jfPt=mGrid->GetBkgPtCut();
  jfRc=mGrid->GetRc();
  jfRcBkg=mGrid->GetBkgCone();
  jfEmcalSet=mGrid->GetEmcalSet();
  //jfPid +=;
}


void ktMuEvent::PrintJetFinderSettings()
{
  cout<<"JetFinder Settings:"<<endl;
  cout<<"-------------------"<<endl;
  cout<<"Rc = "<<jfRc<<endl;
  cout<<"Pt cut = "<<jfPt<<endl;
  cout<<"Rc(bkg) = "<<jfRcBkg<<endl;
  cout<<"Emcal set = "<<jfEmcalSet<<endl;
  cout<<"Particle PID = "<<jfPid<<endl;
}

TString ktMuEvent::GetJetFinderSettings()
{

  TString mS="Rc=";
  mS+= (::Form("%2.2g",jfRc));
  mS+=", Pt,bkg=";
  mS+= (::Form("%2.2g",jfPt));
  mS+=", ";
  mS+=jfPid;

  return mS;
}

void ktMuEvent::PrintFoundJets()
{
  cout<<" Event# "<<evNum<<" : found "<<GetConeJets()->GetEntriesFast()<<" Jets:"<<endl;
  cout<<" JetFinder settings : "<<GetJetFinderSettings()<<endl;
  printf(" %5s %15s %15s %15s %15s %15s %15s %12s %10s %10s %10s %10s\n","jet #", "eta",
         "phi", "pt", "pt,corr", "pt,corr,cone","pt,cluster,corr","n particles","n cells","Type","n subjets","e- tag");

  for (int i=0;i<GetConeJets()->GetEntriesFast();i++)
    {
      ktMuJet *outJet=(ktMuJet*) GetConeJets()->At(i);
      const char* outType2=outJet->GetType();
      
        //DEBUG:
       //outJet->PrintJet();
      
      printf(" %5u %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f %8u %11u %15s %8u %8u \n",i+1,outJet->Eta(),outJet->Phi(),outJet->Pt(),outJet->PtCellCorr(),outJet->PtConeCorr(),outJet->PtClusterCorr(),outJet->GetNParticlesInJet(),outJet->GetNCellsInJet(),outType2,outJet->GetNSubJets(),outJet->IsElectronSeed());
      
    }
}

void ktMuEvent::PrintJets(Double_t minPtCorr)
{
  cout<<" Event# "<<evNum<<" : found #LOCone = "<<GetConeJets()->GetEntriesFast()<<" | #kt = "<<fjKt->GetEntriesFast()<<" | #AntiKt = "<<fjAntiKt->GetEntriesFast()<<" | #SISCone = "<<fjSISCone->GetEntriesFast()<<endl;
  cout<<" JetFinder seetings : "<<GetJetFinderSettings()<<endl;
  cout<<endl;

  printf(" %5s %8s %8s %8s %8s %8s %8s %8s %8s %10s\n","jet #", "eta",
	 "phi", "pt", "pt,corr","area","+-","NEF","NEFCorr","Type");
  
  for (int i=0;i<GetConeJets()->GetEntriesFast();i++)
    {
      ktMuJet *outJet=(ktMuJet*) GetConeJets()->At(i);
      const char* outType2=outJet->GetType();
      printf(" %5u %8.3f %8.3f %8.3f %8.3f %8.4f %8.4f %8.3f %8.3f %8s\n",i+1,outJet->Eta(),outJet->Phi(),outJet->Pt(),outJet->PtConeCorr(),TMath::Pi()*jfRc*jfRc,0.0,0.0,0.0,outType2);
    }

  for (int i=0;i<GetNFJkt();i++)
    {
      ktMuFastJet *outJet=(ktMuFastJet*) GetFJktJets()->At(i);
      const char* outType2=outJet->GetType();
      
      if (outJet->PtCorr()>minPtCorr)
	{
	  printf(" %5u %8.3f %8.3f %8.3f %8.3f %8.4f %8.4f %8.3f %8.3f %8s\n",i+1,outJet->Eta(),outJet->Phi(),outJet->Pt(),outJet->PtCorr(),outJet->Area(),outJet->AreaError(),outJet->GetNEF(),outJet->GetNEFCorr(),outType2);
	 }
    }
  
  for (int i=0;i<GetNFJAntiKt();i++)
     {
       ktMuFastJet *outJet=(ktMuFastJet*) GetFJAntiKtJets()->At(i);
       const char* outType2=outJet->GetType();
       
       if (outJet->PtCorr()>minPtCorr)
	 {
	   printf(" %5u %8.3f %8.3f %8.3f %8.3f %8.4f %8.4f %8.3f %8.3f %8s\n",i+1,outJet->Eta(),outJet->Phi(),outJet->Pt(),outJet->PtCorr(),outJet->Area(),outJet->AreaError(),outJet->GetNEF(),outJet->GetNEFCorr(),outType2);
	 }
     }
  
  for (int i=0;i<GetNFJSISCone();i++)
    {
      ktMuFastJet *outJet=(ktMuFastJet*) GetFJSISConeJets()->At(i);
      const char* outType2=outJet->GetType();
      
      if (outJet->PtCorr()>minPtCorr)
	 {
	   printf(" %5u %8.3f %8.3f %8.3f %8.3f %8.4f %8.4f %8.3f %8.3f %8s\n",i+1,outJet->Eta(),outJet->Phi(),outJet->Pt(),outJet->PtCorr(),outJet->Area(),outJet->AreaError(),outJet->GetNEF(),outJet->GetNEFCorr(),outType2);
	 }
    }

  
}

void ktMuEvent::PrintFastJets(Double_t minPtCorr)
{
  cout<<" Event# "<<evNum<<" : found #kt = "<<fjKt->GetEntriesFast()<<" | #AntiKt = "<<fjAntiKt->GetEntriesFast()<<" | #SISCone = "<<fjSISCone->GetEntriesFast()<<endl;
  
   printf(" %5s %8s %8s %8s %8s %8s %10s\n","jet #", "eta",
	  "phi", "pt", "pt,corr","area","Type");

   for (int i=0;i<GetNFJkt();i++)
     {
       ktMuFastJet *outJet=(ktMuFastJet*) GetFJktJets()->At(i);
       const char* outType2=outJet->GetType();
       
       if (outJet->PtCorr()>minPtCorr)
	 {
	   printf(" %5u %8.3f %8.3f %8.3f %8.3f %8.3f %8s\n",i+1,outJet->Eta(),outJet->Phi(),outJet->Pt(),outJet->PtCorr(),outJet->Area(),outType2);
	 }
     }

    for (int i=0;i<GetNFJAntiKt();i++)
     {
       ktMuFastJet *outJet=(ktMuFastJet*) GetFJAntiKtJets()->At(i);
       const char* outType2=outJet->GetType();
       
       if (outJet->PtCorr()>minPtCorr)
	 {
	   printf(" %5u %8.3f %8.3f %8.3f %8.3f %8.3f %8s\n",i+1,outJet->Eta(),outJet->Phi(),outJet->Pt(),outJet->PtCorr(),outJet->Area(),outType2);
	 }
     }

     for (int i=0;i<GetNFJSISCone();i++)
     {
       ktMuFastJet *outJet=(ktMuFastJet*) GetFJSISConeJets()->At(i);
       const char* outType2=outJet->GetType();
       
       if (outJet->PtCorr()>minPtCorr)
	 {
	   printf(" %5u %8.3f %8.3f %8.3f %8.3f %8.3f %8s\n",i+1,outJet->Eta(),outJet->Phi(),outJet->Pt(),outJet->PtCorr(),outJet->Area(),outType2);
	 }
     }

}

void ktMuEvent::AddFastJet(Double_t mpt,Double_t meta,Double_t mphi,Double_t metaCorr,Double_t mphiCorr,Double_t mMedian, Double_t marea,TString mtype, Double_t mNEF, Double_t mNEFCorr,Double_t mAerror)
{

  ktMuFastJet *mJ=new ktMuFastJet(mpt,meta,mphi,metaCorr,mphiCorr,mMedian,marea,mtype);
  mJ->SetNEF(mNEF);
  mJ->SetNEFCorr(mNEFCorr);
  mJ->SetAreaError(mAerror);

  // DEBUG:
  //mJ->PrintJet();

  if (mtype.Contains("Anti"))
    {
      fjAntiKt->AddLast(mJ);
      fjAntiKt->Sort();
    }
  else if (mtype.Contains("SIS"))
    {
      fjSISCone->AddLast(mJ);
      fjSISCone->Sort();
    }
  else
    {
      fjKt->AddLast(mJ); 
      fjKt->Sort();
    }

}


void ktMuEvent::FillAll(TObjArray *mJets)
{
  for (int i=0;i<mJets->GetEntriesFast();i++)
    {
      allList->AddLast(mJets->At(i));
    }
}

void ktMuEvent::AddUEEvent(TObjArray *mJets, ktGrid *grid){

  if( mJets->GetEntries() > 0){
    ktUEEvent* ueEvent = new ktUEEvent();
    ueEvent->Fill(mJets, grid);
    TString type = ((ktMuJet *)(mJets->At(0)))->GetType();
    ueEvent->SetJetType((char *)type.Data());
    ueEvents->AddLast(ueEvent);
  }
}


void ktMuEvent::AddUEEvent(TObjArray *mJets, ktFastJet *fJet){

  if( mJets->GetEntries() > 0){
    ktUEEvent* ueEvent = new ktUEEvent();
    ueEvent->Fill(mJets, fJet);
    ueEvent->SetJetType((char *)((ktMuFastJet *)(mJets->At(0)))->GetType().Data());
    ueEvents->AddLast(ueEvent);
  }
}


void ktMuEvent::EventSummary()
{
  cout<<"Event number : "<<evNum<<endl;
  cout<<"# Cone jets = "<<GetNConeJets()<<endl;
  cout<<"# kt jets      = "<<GetNKtJets()<<endl;
  cout<<"# ktBF jets  = "<<GetNKtBFJets()<<endl;
  cout<<"# partons   = "<<GetNPartons()<<endl;
}

/*
void  ktMuEvent::SetTriggerInfo(TClonesArray* trinfo){
  
  //DEBUG:
  //cout<<" ------  ktMuEvent::SetTriggerInfo --------"<<endl;
  //cout<<" N. of trig obj in the event = "<<TriggerInfoArray->GetEntriesFast()<<endl;

  if (trinfo->GetEntriesFast()>0)
    {
      TClonesArray &trigobj = *TriggerInfoArray;
      for( int i=0; i<trinfo->GetEntriesFast(); i++){
	ktTriggerInfo *t = (ktTriggerInfo *)trinfo->At(i);
	ktTriggerInfo *trig = new(trigobj[i]) ktTriggerInfo(*t);
	//trig->PrintInfo();
      }
    }
}
*/

void  ktMuEvent::SetTriggerInfo(TClonesArray* trinfo){
  
  //DEBUG:
  //cout<<" ------  ktMuEvent::SetTriggerInfo --------"<<endl;
  //cout<<" N. of trig obj in the event = "<<TriggerInfoArray->GetEntriesFast()<<endl;

  if (trinfo->GetEntriesFast()>0)
    {
      for( int i=0; i<trinfo->GetEntriesFast(); i++){
	ktTriggerInfo *t = (ktTriggerInfo *)trinfo->At(i);
	//ktTriggerInfo *trig = new ktTriggerInfo(*t);
	ktTriggerInfo *trig = new ktTriggerInfo(t);
	TriggerInfoArray->AddLast(trig);

	//trig->PrintInfo();
      }
    }
}

void ktMuEvent::PrintTriggerInfo()
{

  cout<<" ------  ktMuEvent::PrintTriggerInfo --------"<<endl;
  cout<<" N. of trig obj in the event = "<<GetTriggerInfo()->GetEntriesFast()<<endl;

  for( int i=0; i<(Int_t) GetTriggerInfo()->GetEntriesFast(); i++)
    {
      ktTriggerInfo *t = (ktTriggerInfo*) GetTriggerInfo()->At(i);          
      t->PrintInfo();
    }
}

Double_t ktMuEvent::GetHTEta()
{
  for( int i=0; i<(Int_t) GetTriggerInfo()->GetEntriesFast(); i++)
    {
      ktTriggerInfo *t = (ktTriggerInfo*) GetTriggerInfo()->At(i);          
      if (t->isHTL0())
	return (Double_t) t->GetEta();
    }

  return -99;
}

Double_t ktMuEvent::GetHTPhi()
{
  for( int i=0; i<(Int_t) GetTriggerInfo()->GetEntriesFast(); i++)
    {
      ktTriggerInfo *t = (ktTriggerInfo*) GetTriggerInfo()->At(i); 
      if (t->isHTL0())
	{
	  Double_t phi=(Double_t) t->GetPhi();
	  if (phi<0) phi += (2*TMath::Pi());
	  return phi;
	}
    
    }

  return -99;
}

Double_t ktMuEvent::GetJPEta()
{
  for( int i=0; i<(Int_t) GetTriggerInfo()->GetEntriesFast(); i++)
    {
      ktTriggerInfo *t = (ktTriggerInfo*) GetTriggerInfo()->At(i);          
      if (t->isJPL0())
	return (Double_t) t->GetEta();
    }
  return -99;
}

Double_t ktMuEvent::GetJPPhi()
{
  for( int i=0; i<(Int_t) GetTriggerInfo()->GetEntriesFast(); i++)
    {
      ktTriggerInfo *t = (ktTriggerInfo*) GetTriggerInfo()->At(i);          
      if (t->isJPL0())
	 return (Double_t) t->GetPhi();
    }
  return -99;
}
