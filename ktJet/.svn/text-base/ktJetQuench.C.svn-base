// first test (Joern Putschke)

#include "ktJetQuench.h"
#include "TRandom.h"
#include <Riostream.h>
#include "TCanvas.h"
#include "TFile.h"

ClassImp(ktJetQuench)

ktJetQuench::ktJetQuench() 
{

  hJetParticles=new TObjArray(0);
  hJetParticles->SetOwner(kTRUE);
  hJetParticlesNoQuench=new TObjArray(0);
  hJetParticlesNoQuench->SetOwner(kTRUE);

  jPt=jPhi=jEta=0.0;
  S0=0.25;
  T=0.291;
  C=2.0;
  R=0.7;
  minQpt=0.2;
  maxQpt=15;
  minDi=3.0;
  maxDi=4.0;
  assDi=2.0;
  nTrig=0;
  nTrigNo=0;
  mQMethod="fractinal";
  PE=5.0;

  dE=0.0;
  ptQuench=ptNoQuench=0.0;
  maxIter=100;
  nEvents=0;
  
  verbose=false;
  doQA=false;
  init=false;
  hScale=false;
  fracQuench=true;

  //DEBUG:
  //cout<<endl;
  //cout<<"Default constructor of ktJetQuench() ..."<<endl;

  //feta=0;fphi=0;fpt=0;
  fpt=0;

  hXi=new TH1D("hXi","#xi distribution quenched",50,0,10);
  hXi->SetDirectory(0);
  hXiNoQuench=new TH1D("hXiNoQuench","#xi distribution no quenching",50,0,10);
  hXiNoQuench->SetDirectory(0);
  hPt=new TH1D("hPt","Jet pt",75,0,75);
  hPt->SetDirectory(0);
  hJptQ=new TH1D("hJptQ","Jet pt quenched",75,0,75);
  hJptQ->SetDirectory(0);
  hJptRe=new TH1D("hJptRe","Jet pt redistributed",75,0,75);
  hJptRe->SetDirectory(0);
  hPtQ=new TH1D("hPtQ","Charged particle pt quenched",60,0,15);
  hPtQ->SetDirectory(0);
  hPtNo=new TH1D("hPtNo","Charged particle pt no-quenched",60,0,15);
  hPtNo->SetDirectory(0);
  hdEtadPhi=new TH2D("hdEtadPhi","dEta vs. dPhi quenched",180,-TMath::Pi(),TMath::Pi(),50,-2,2);
  hdEtadPhi->SetDirectory(0);
  hdEtadPhiNo=new TH2D("hdEtadPhiNo","dEta vs. dPhi no-quenched",180,-TMath::Pi(),TMath::Pi(),50,-2,2);
  hdEtadPhiNo->SetDirectory(0);
  hdRQ=new TH1D("hdRQ","dR for redistributed quenched particles",20,0,1);
  hdRQ->SetDirectory(0);
  hdR=new TH1D("hdR","dR for unquenched partticles",20,0,1);
  hdR->SetDirectory(0);
  hdRptQ=new TH2D("hdRptQ","dR vs. pt for redistributed quenched particles",20,0,1,60,0,15);
  hdRptQ->SetDirectory(0);
  hdRpt=new TH2D("hdRpt","dR vs. pt for unquenched partticles",20,0,1,60,0,15);
  hdRpt->SetDirectory(0);
}


void ktJetQuench::SetQuenchPtRange(Double_t m_min,Double_t m_max)
{
  minQpt=m_min;maxQpt=m_max;

}
void ktJetQuench::InitQuenching(TString qMethod)
{
   // DEBUG:
  //cout<<"Set functions ..."<<endl;
  init=true;

  /*
  feta=new TF1("feta","1",0,R);
  fphi=new TF1("fphi","1",0,R);
  */

  /*
  if (qMethod.Contains("frac"))
    fracQuench=true;
  else
    fracQuench=false;
  */

  mQMethod=qMethod;

  fpt=new TF1("fpt","x*exp(-x/[0])",minQpt,maxQpt);
  fpt->SetParameter(0,T);

  cout<<"Quenching info:"<<endl;
  cout<<"---------------"<<endl;
  if (fracQuench)
    cout<<"Quenching fraction S0 = "<<S0<<endl;
  else
    cout<<"Constant quenching C = "<<C<<" GeV"<<endl;

  cout<<"Quenching slope (thermal spectrum) T = "<<T<<endl;
  cout<<"in pt range "<<minQpt<<" .. "<<maxQpt<<" GeV"<<endl;
  cout<<"Uniform in R around jet-axis with R = "<<R<<endl;
  cout<<"Max iterations = "<<maxIter<<endl;

}

void ktJetQuench::AddParticle(Double_t px, Double_t py, Double_t pz, Double_t E)
{

  // DEBUG:
  //cout<<E<<endl;
  //E=TMath::Sqrt(px*px+py*py);
  //cout<<E<<endl;

  //ktParticle *p=new ktParticle();
  //p->SetPxPyPzE(px,py,pz,E);
  //cout<<" Add particle ..."<<endl;

  ktParticle *pNo=new ktParticle();
  pNo->SetPxPyPzE(px,py,pz,E);

  //hJetParticles->AddLast(p);
  hJetParticlesNoQuench->AddLast(pNo);

}

void ktJetQuench::AddParticle(Double_t px, Double_t py, Double_t pz, Double_t E, int mPid)
{

  // DEBUG:
  //cout<<E<<endl;
  //E=TMath::Sqrt(px*px+py*py);
  //cout<<E<<endl;

  //ktParticle *p=new ktParticle();
  //p->SetPxPyPzE(px,py,pz,E);
  //cout<<" Add particle ..."<<endl;

  ktParticle *pNo=new ktParticle();
  pNo->SetPxPyPzE(px,py,pz,E);

  ktPID *mPID = new ktPID();
  //mPID->SetPID(mPid);
  mPID->SetCharge(mPid);
  pNo->SetPid(mPID);

  //hJetParticles->AddLast(p);
  hJetParticlesNoQuench->AddLast(pNo);

}

// check sumn of 4-vector vs. sum of pt ...
Double_t ktJetQuench::GetSumPt()
{
  Double_t mSum=0.0;
  
  for (int i=0;i<GetNJetParticles();i++)
    {
      ktParticle *p=(ktParticle*) hJetParticles->At(i);
      mSum += p->Pt();
    }

  return mSum;
}

Double_t ktJetQuench::GetSumPtNoQuench()
{
  Double_t mSum=0.0;

  for (int i=0;i<GetNJetParticlesNoQuench();i++)
    {
      ktParticle *p=(ktParticle*) hJetParticlesNoQuench->At(i);
      mSum += p->Pt();
    }

  return mSum;
}

void ktJetQuench::FillQA()
{

  if (verbose)
    {
      cout<<" ---> Fill QA histograms ..."<<endl;
    }

  hPt->Fill(jPt);

  for (int i=0;i<GetNJetParticles();i++)
    {
      ktParticle *p=(ktParticle*) hJetParticles->At(i);
      Double_t xi=TMath::Log((Double_t) jPt/(Double_t) p->Pt());
      if (p->GetPid()>0)
	{
	  hXi->Fill(xi);
	  hPtQ->Fill(p->Pt(),1/p->Pt());
	}
    }

  for (int i=0;i<GetNJetParticlesNoQuench();i++)
    {
      ktParticle *p=(ktParticle*) hJetParticlesNoQuench->At(i);
      Double_t xi=TMath::Log((Double_t) jPt/(Double_t) p->Pt());
      if (p->GetPid()>0)
	{
	  hXiNoQuench->Fill(xi);
	  hPtNo->Fill(p->Pt(),1/p->Pt());
	}
    }
}

void ktJetQuench::DoDiHadron(Double_t minTrig, Double_t maxTrig, Double_t ptAss)
{

  if (verbose)
    {
      cout<<" ---> Do di-hadron analysis "<<minTrig<<" < pt,trig < "<<maxTrig<<" and ptassoc > "<<ptAss<<endl;
    }
  
  minDi=minTrig;maxDi=maxTrig;assDi=ptAss;
  
  DoDiCorr(hJetParticles,hdEtadPhi,true);
  DoDiCorr(hJetParticlesNoQuench,hdEtadPhiNo,false);

}

void ktJetQuench::DoDiCorr(TObjArray *mA,TH2D *mH,Bool_t quenched)
{
  Int_t evNtrig=0;
  TArrayI *trigIndex=new TArrayI(0);

  for (int i=0;i<(Int_t) mA->GetEntries();i++)
    {
      ktParticle *p=(ktParticle*) mA->At(i);
      // DEBUG:
      //cout<<i<<" "<<p->Pt()<<endl;
      
      if (p->Pt()>minDi && p->Pt()<maxDi && p->GetPid()>0)
	{
	  // DEBUG:	  
	  //cout<<evNtrig<<" "<<i<<endl;

	  trigIndex->Set(evNtrig+1);
	  trigIndex->AddAt(i,evNtrig);

	  evNtrig++;
	  if (quenched)
	    nTrig++;
	  else
	    nTrigNo++;
	}
    }

  // DEBUG:
  //cout<<"Found Ntrig "<<evNtrig<<" in event."<<endl;
  //cout<<trigIndex->GetSize()<<endl;

  for (int j=0;j<(Int_t) trigIndex->GetSize();j++)
    {
      Int_t mIndex=(Int_t) trigIndex->At(j);
      ktParticle *pTrig=(ktParticle*) mA->At(mIndex);

      Double_t etaTrig=pTrig->Eta();
      Double_t phiTrig=pTrig->Phi();
      Double_t ptTrig=pTrig->Pt();

      // DEBUG:
      //cout<<"   "<<pTrig->Pt()<<" "<<pTrig->Eta()<<" "<<pTrig->Phi()<<endl;

      for (int i=0;i<(Int_t) mA->GetEntries();i++)
	{
	  ktParticle *p=(ktParticle*) mA->At(i);	
	  
	  if (p->Pt()>assDi && p->Pt()<ptTrig && i!=mIndex && p->GetPid()>0)
	    {
	      // DEBUG:
	      //cout<<p->Pt()<<" "<<p->Eta()<<" "<<p->Phi()<<endl;

	      Double_t dEta=etaTrig-p->Eta();
	      Double_t dPhi=phiTrig-p->Phi();

	      mH->Fill(dPhi,dEta);

	    }
	}
    }

  delete trigIndex;
}

void ktJetQuench::DrawQAPlots(TString mT)
{
   if (nEvents>0 && doQA)
     {
       TCanvas *cQuench=new TCanvas("cQuench_"+mT,"xi distributions filled from "+mT,800,800);
       cQuench->Divide(3,3);  
       
       hScale=true;
       
       cQuench->cd(1);
       hXi->Scale(1.0/(hXi->GetBinWidth(1)*nEvents));
       hXiNoQuench->Scale(1.0/(hXiNoQuench->GetBinWidth(1)*nEvents));
       
       hXi->SetLineColor(2);
       hXi->DrawCopy();
       hXiNoQuench->DrawCopy("same");
       
       cQuench->cd(2);
       TH1D *hratio=(TH1D*) hXi->Clone();
       hratio->Divide(hXiNoQuench);
       hratio->DrawCopy();
       
       cQuench->cd(3);
       hPt->DrawCopy();
       hJptQ->SetLineColor(2);
       hJptRe->SetLineColor(4);
       hJptQ->DrawCopy("same");
       hJptRe->DrawCopy("same");

       cQuench->cd(4);
       gPad->SetLogy();
       hPtQ->Scale(1/hPtQ->GetBinWidth(1));
       hPtNo->Scale(1/hPtNo->GetBinWidth(1));
       hPtQ->SetLineColor(2);
       
       hPtQ->DrawCopy();
       hPtNo->DrawCopy("same");

       cQuench->cd(5);
       TH1D *hraa=(TH1D*) hPtQ->Clone();
       hraa->Divide(hPtNo);
       hraa->DrawCopy();
       
       /*
       cQuench->cd(6);
       hdEtadPhi->Scale(1/(Double_t) nTrig);// *1/hdEtadPhi->GetXaxis()->GetBinWidth(1)*1/hdEtadPhi->GetYaxis()->GetBinWidth(1));
       hdEtadPhi->DrawCopy("colz");
       
       cQuench->cd(7);
       hdEtadPhiNo->Scale(1/(Double_t) nTrigNo);// *1/hdEtadPhiNo->GetXaxis()->GetBinWidth(1)*1/hdEtadPhiNo->GetYaxis()->GetBinWidth(1));
       hdEtadPhiNo->DrawCopy("colz");
       */

       cQuench->cd(6);
       TH1D* hdPhi=(TH1D*) hdEtadPhi->ProjectionX();
       hdPhi->Scale(1.0/hdPhi->GetBinWidth(1));
       hdPhi->SetLineColor(2);
       TH1D* hdPhiNo=(TH1D*) hdEtadPhiNo->ProjectionX();
       hdPhiNo->Scale(1.0/hdPhiNo->GetBinWidth(1));
       hdPhiNo->DrawCopy();
       hdPhi->DrawCopy("same");
       
       cQuench->cd(7);
       hdR->DrawCopy();
       hdRQ->SetLineColor(2);
       hdRQ->DrawCopy("same");

       cQuench->cd(8);
       TH1D* hptR=(TH1D*) hdRpt->ProfileX();
       TH1D* hptRQ=(TH1D*) hdRptQ->ProfileX();
       hptRQ->SetLineColor(2);
       hptR->DrawCopy();
       hptRQ->DrawCopy("same");
       
     }
   else
     cout<<" ---> No events analyzed to show QA plots (ktJetQuench) !!!"<<endl;
}


void ktJetQuench::SaveQA(TString fName)
{
  TFile *f=new TFile(fName,"RECREATE");

  if (!hScale)
    {
      hXi->Scale(1.0/(hXi->GetBinWidth(1)*nEvents));
      hXiNoQuench->Scale(1.0/(hXiNoQuench->GetBinWidth(1)*nEvents));
      hPtQ->Scale(1/hPtQ->GetBinWidth(1));
      hPtNo->Scale(1/hPtNo->GetBinWidth(1));
    }
     
  hXi->SetLineColor(2);
  hXi->Write("hXi_Quenched");
  hXiNoQuench->Write("hXi_NoQuenched");
  hXi->Divide(hXiNoQuench);
  hXi->Write("hFFRatio");
  hPt->Write("hJetPt");
  hPtQ->SetLineColor(2);
  hPtQ->Write("hPt_Quenched");
  hPtNo->Write("hPt_NoQuenched");
  hdEtadPhi->Write("hdEtadPhi_Quenched");
  hdEtadPhiNo->Write("hdEtadPhi_NoQuenched");
  hdRQ->Write("hdR_Quenched");
  hdR->Write("hdR_UnQuenched");
  hdRptQ->Write("hdRvsPt_Quenched");
  hdRpt->Write("hdRvsPt_UnQuenched");
  hJptQ->Write("hJetPtQuenched");
  hJptRe->Write("hJetPtRedistributed");
  f->Write();
  f->Close();
  
  //cout<<endl;
  cout<<" ---> ktQuench QA plots saved in : "<<fName<<endl;
  //cout<<endl;
  
}

void ktJetQuench::DoQuenching()
{
  if (verbose)
    {
      cout<<" ---> Do quenchinhg with S0 = "<<S0<<endl;
    }

  if (init)
    if (GetNJetParticlesNoQuench()>0)
      {
	for (int i=0;i<GetNJetParticlesNoQuench();i++)
	  {
	    ktParticle *pNo=(ktParticle*) hJetParticlesNoQuench->At(i);
	    
	    ktParticle *p=new ktParticle();
	 
	    ktPID *myPID = new ktPID();
	    myPID->SetCharge(pNo->IsCharged());
	    p->SetPid(myPID);

	    if (pNo->IsCharged()>0)
	      {
		Double_t deta=pNo->Eta()-GetJetEta();
		Double_t dphi=pNo->Phi()-GetJetPhi();	 
		Double_t dR=TMath::Sqrt(deta*deta+dphi*dphi);
		hdR->Fill(dR);
		hdRpt->Fill(dR,pNo->Pt());
	      }

	    if (mQMethod.Contains("frac"))
	      {
		Double_t ptQ=pNo->Pt()*(1-S0);
		p->SetPtEtaPhiM(ptQ,pNo->Eta(),pNo->Phi(),0);	    
		hJetParticles->AddLast(p);
	      }	    
	    else if (mQMethod.Contains("const"))
	      {		
		Double_t ptQ=pNo->Pt()-C;
		if (ptQ>0)
		  {
		    p->SetPtEtaPhiM(ptQ,pNo->Eta(),pNo->Phi(),0);	    
		    hJetParticles->AddLast(p);		
		  }	
		else
		  delete p;
	      }
	    else
	      {
		p->SetPtEtaPhiM(pNo->Pt(),pNo->Eta(),pNo->Phi(),0);
		hJetParticles->AddLast(p);	
	      }

	  }
	
	if (mQMethod.Contains("part"))
	  {
	    dE=PE;
	    jPt=jPt+PE;
	    ptQuench=GetSumPt()+PE; ptNoQuench=GetSumPtNoQuench();
	  }
	else
	  {
	    ptQuench=GetSumPt(); ptNoQuench=GetSumPtNoQuench();
	    dE=ptNoQuench-ptQuench;
	  }
	
	hJptRe->Fill(dE);
	hJptQ->Fill(ptQuench);

	ReDistribute();
	
	if (doQA)
	  FillQA();
	
	nEvents++;
      }
    else
      cout<<" ---> Warning (ktJetQuench) no particles/jet !!!"<<endl;
  else
   cout<<" ---> Warning ktJetQuench not initialized do: InitQuenching() !!!"<<endl; 

}

void ktJetQuench::ReDistribute()
{
  if (verbose)
    {
      cout<<" ---> Redistribute quenched dE "<<dE<<endl;
    }

  // DEBUG:
  //cout<<endl;
  //cout<<GetSumPt()<<" "<<GetSumPt()+dE<<" "<<GetSumPtNoQuench()<<endl;
  //cout<<endl;

  Int_t niter=0;
  Double_t pt,phi,eta,mdE;
  pt=phi=eta=0;
  mdE=dE;
  
  for (int i=0;i<maxIter;i++)
    {

      if (mdE<0) continue;

      Double_t mR=gRandom->Uniform(R);
      gRandom->Circle(eta,phi,mR);
      
      Double_t dR=TMath::Sqrt(eta*eta+phi*phi);

      // DBEUG:
      //cout<<eta<<" "<<phi<<" "<<dR<<endl;

      eta=eta+jEta;
      phi=phi+jPhi;
      
      pt=fpt->GetRandom();
      
      // DEBUG:
      //cout<<pt<<" "<<eta<<" "<<phi<<" "<<mdE-pt<<endl;

      mdE=mdE-pt;

      ktParticle *p=new ktParticle();
      p->SetPtEtaPhiM(pt,eta,phi,0);
      
      if (gRandom->Rndm()>0.66){
	ktPID *myPID = new ktPID();
	myPID->SetCharge(0);
	p->SetPid(myPID);
      }
      else{
	ktPID *myPID = new ktPID();
	myPID->SetCharge(1);
	p->SetPid(myPID);
      }
      if (p->IsCharged()>0)
	{
	  hdRQ->Fill(dR);
	  hdRptQ->Fill(dR,pt);
	}

      hJetParticles->AddLast(p);     

      niter++;
    }

  ktParticle *mp=(ktParticle*) hJetParticles->At(hJetParticles->GetLast());
  mp->SetPtEtaPhiM(mp->Pt()+mdE,mp->Eta(),mp->Phi(),0);  

  // DBEUG:
  //cout<<mdE<<" "<<niter<<endl;

}

void ktJetQuench::PrintInfo()
{
  cout<<"ktJetQuench info:"<<endl;
  cout<<"-----------------"<<endl;
  cout<<"Jet (no quench) pt = "<<jPt<<endl;
  cout<<"Jet (no quench) phi = "<<jPhi<<endl;
  cout<<"Jet (no quench) eta = "<<jEta<<endl;
  cout<<"# quenched particles = "<<GetNJetParticles()<<endl;
  cout<<"# un-quenched particles = "<<GetNJetParticlesNoQuench()<<endl;
  cout<<"Quenching fraction S0 = "<<S0<<endl;
  cout<<"ptSum quenched   = "<<GetSumPt()<<endl;
  cout<<"ptSum unquenched = "<<GetSumPtNoQuench()<<endl;

}

void ktJetQuench::Clear()
{
  hJetParticles->Delete();
  hJetParticlesNoQuench->Delete();

  hJetParticles->Compress();
  hJetParticlesNoQuench->Compress();
}

ktJetQuench::~ktJetQuench()
{

  hJetParticles->Delete();
  delete hJetParticles;
  hJetParticlesNoQuench->Delete();
  delete hJetParticlesNoQuench;

  if (fpt!=0)
    {
      //delete feta;delete fphi;delete fpt;
      delete fpt;
    }

  delete hXi;delete hXiNoQuench; delete hPt;
  delete hPtQ;delete hPtNo;
  delete hdEtadPhi; delete hdEtadPhiNo;
  delete hdRQ; delete hdR;
  delete hdRptQ; delete hdRpt;
  delete hJptQ; delete hJptRe;
  //DEBUG:
  //cout<<endl;
  //cout<<"Destructor of ktJetQuench() ..."<<endl;
}
