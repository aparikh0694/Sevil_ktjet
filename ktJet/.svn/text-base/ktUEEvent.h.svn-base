// first test (Helen Caines)

#ifndef ROOT_ktUEEvent
#define ROOT_ktUEEvent

#include "TObject.h"
#include "TObjArray.h" 
#include "TMath.h" 
#include "TString.h"
#include "ktJet.h"
#include "ktMuJet.h" 
#include "ktGrid.h" 
#include "ktMuFastJet.h"
#include "ktParticle.h"

class ktFastJet;

class ktUEEvent : public TObject 
 { 

   private: 

   TObjArray* UEParticles;
   Int_t trackMult;
   Int_t towerMult;
   Int_t matchMult;

   Int_t nJet;
   Double_t jetE;
   Double_t jetPt;
   Double_t jetPhi;
   Double_t jetEta;


   TString jetType;

   public: 
   
   ktUEEvent(); 
   virtual ~ktUEEvent(); 
 
   //static int memCounter;

   Int_t GetNJets(){return nJet;}
   Double_t GetJetE(){return jetE;}
   Double_t GetJetPt(){return jetPt;}
   Double_t GetJetEta(){return jetEta;}
   Double_t GetJetPhi(){return jetPhi;}

   Int_t GetAllTrackMult(){return trackMult+matchMult;}
   Int_t GetTrackMult(){return trackMult;}
   Int_t GetMatchMult(){return matchMult;}
   Int_t GetTowerMult(){return towerMult;}
   TString GetJetType(){return jetType;}
   
   TObjArray* GetParticles(){return UEParticles;}
   ktParticle* GetParticle(Int_t n){return (ktParticle*)UEParticles->At(n);}
   
    void SetTowerMult(Int_t mMult){towerMult=mMult;}
   void SetTrackMult(Int_t mMult){trackMult=mMult;}
   void SetMatchMult(Int_t mMult){matchMult=mMult;}
   void SetJetType(char* type){jetType = type;}

   void Fill( TObjArray *mJets, ktGrid *grid);
   void Fill( TObjArray *mJets, ktFastJet *fastjet);
   void EventSummary(); 

  ClassDef(ktUEEvent,2) 
 }; 

 #endif 
