// first test (Helen Caines)

#include "ktUEEvent.h"
#include "ktFastJet.h"
#include <Riostream.h>
#include <math.h>

//int ktUEEvent::memCounter=0;

ClassImp(ktUEEvent)

ktUEEvent::ktUEEvent()
{

  // memCounter++;
  nJet = 0;
  jetE = -999;
  jetPt = -999;
  jetPhi = -999;
  jetEta = -999;

  jetType = "";

  towerMult = 0;
  trackMult = 0;
  matchMult = 0;
  
  UEParticles = 0;
  // UEParticles = new TObjArray(0);
  //UEParticles->SetOwner(kTRUE);


  // DEBUG:
  //cout<<"Default ktUEEvent constructor "<<endl;
}

//-------------------------------------------------------------------------
ktUEEvent::~ktUEEvent()
{
  
  // memCounter--;
  //cout << "Deleting my particles !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << memCounter << endl;
   UEParticles->Delete();
  delete UEParticles;


  //cout << "Deleting my particles !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  // DEBUG:
  //cout<<"Default ktUEEvent destructor"<<endl;
}

//-------------------------------------------------------------------------
void ktUEEvent::Fill(TObjArray *mJets, ktGrid *grid)
{

  if( UEParticles == 0){
    UEParticles = new TObjArray(0);
    UEParticles->SetOwner(kTRUE); 
  }

  if( mJets->GetEntries()==0) return;

  // Find max E jet
  nJet = mJets->GetEntries();
   for (int i=0;i<nJet;i++){
     
     //     ktJet *mktJet=(ktJet*) (mJets->At(i));
     ktMuJet *mJet=  (ktMuJet *)(mJets->At(i));
     if( mJet->Pt() > jetPt){

	jetE = mJet->E();
	jetPhi = mJet->Phi();
	if( jetPhi > TMath::Pi()) jetPhi -= 2*TMath::Pi();
	jetPt = mJet->Pt();
	jetEta = mJet->Eta();
      }
    }
   
   // Fill UEvent with data outside pi/3 (and same on awayside) of leading jet

   double dPhi;
   
   double nParticles = 0;

   for( int iPhi=0; iPhi<grid->GetNphi(); iPhi++){
     for( int iEta=0; iEta<grid->GetNeta(); iEta++){
       ktJetCell *cell = grid->GetCell(iPhi,iEta);
       TObjArray* particlesPID = cell->GetParticleListPID();
       TObjArray* particles = cell->GetParticleList();
       nParticles += cell->NParticlesPID();
       if( cell->NParticlesPID() != cell->NParticles()) cout << "PANIC!!!!" << endl;
       for( int iPart=0; iPart<cell->NParticlesPID(); iPart++){
	 ktPID *pidPart = (ktPID *) particlesPID->At(iPart);
	 ktParticle *particle = (ktParticle *)particles->At(iPart);

	 dPhi = particle->Phi()-jetPhi;
	 if( dPhi > TMath::Pi()) dPhi -= 2*TMath::Pi();
	 if( dPhi < -1*TMath::Pi()) dPhi += 2*TMath::Pi();
	 if( fabs(dPhi) < TMath::Pi()/3. || fabs(dPhi) > 2.*TMath::Pi()/3.) continue;


	 if( pidPart->GetIsMatched()) matchMult;
	 else if( pidPart->GetIsTrack()) trackMult++;
	 // If it isnt a track it must be a tower
	 else towerMult++;

	 ktParticle *copyParticle = new ktParticle(particle, pidPart, kFALSE);
	 UEParticles->AddLast(copyParticle);
       }
     }
   }

   //cout << "NParticles in grid " << nParticles << endl;
}
//-------------------------------------------------------------------------
void ktUEEvent::Fill(TObjArray *mJets, ktFastJet *fastJet)
{

  if( UEParticles == 0){
    UEParticles = new TObjArray(0);
    UEParticles->SetOwner(kTRUE); 
  }

  if( mJets->GetEntries()==0) return;
  
  // Find max E jet
  nJet = mJets->GetEntries();
   for (int i=0;i<nJet;i++){     
     ktMuFastJet *mJet= (ktMuFastJet *) mJets->At(i);
     if( mJet->PtCorr() > jetPt){
       
       jetE = mJet->PtCorr();
       jetPhi = mJet->PhiCorr();
       if( jetPhi > TMath::Pi()) jetPhi -= 2*TMath::Pi();
       jetPt = mJet->PtCorr();
       jetEta = mJet->EtaCorr();
     }
   }
   
   // Fill UEvent with data outside pi/3 (and same on awayside) of leading jet
   
   double dPhi;
   TObjArray *FFArray = fastJet->GetFFArray();
   if( FFArray == 0) return;
   //cout << " FF Array has size " << FFArray->GetEntries() << endl;

   for( int iPart=0; iPart<FFArray->GetEntries(); iPart++){
     ktParticle *particle = (ktParticle *)FFArray->At(iPart);
     dPhi = particle->Phi()-jetPhi;
     if( dPhi > TMath::Pi()) dPhi -= 2*TMath::Pi();
    if( dPhi < -1*TMath::Pi()) dPhi += 2*TMath::Pi();
     //cout << "Track phi "<< particle->Phi()*180./TMath::Pi() << " dPhi is " << dPhi*180./TMath::Pi() << " and jet phi is " <<  jetPhi*180./TMath::Pi()<<endl;
    if( fabs(dPhi) < TMath::Pi()/3. || fabs(dPhi) > 2.*TMath::Pi()/3.) continue;
  
     if( particle->GetPid()->GetIsMatched()) matchMult++;
     else  if( particle->GetPid()->GetIsTrack()) trackMult++;
     // If it isnt a track it must be a tower
     else towerMult++;
     ktParticle *copyParticle = new ktParticle(particle);
     copyParticle->SetIsInJet(false);
     UEParticles->AddLast(copyParticle);
     
    
   }
}
//-------------------------------------------------------------------------
void ktUEEvent::EventSummary()
{
  cout<<endl;
  cout<<"ktUEEvent Summary follows "<< endl;
  cout<<"Jet Type is: "<<GetJetType()<<endl;
  cout<<"Number of Jets is :" << GetNJets()<<endl;
  cout<<"Jet E and Phi : " << GetJetE() << " " << GetJetPhi() << endl;
  cout<<"Underlying Event UnMatched Track Mult : "<<GetTrackMult()<<endl; 
  cout<<"Underlying Event Matched Track Mult : "<<GetMatchMult()<<endl;  
  cout<<"Underlying Event Tower Mult : "<<GetTowerMult()<<endl;
}
