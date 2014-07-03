void kt_test_MCBkg(Int_t nEvents=1000,Double_t mv2=0.0,Bool_t pureMC=true)
{

  ktMCBkg *b=new ktMCBkg("RHIC","Au+Au",mv2);
 
  if (nEvents>10)
    b->SetOutput(false);

  if (b->Init())
    {
      for (int i=0;i<nEvents;i++)
	{
	  if (!pureMC)
	    {
	      ktGrid *grid=new ktGrid();
	      grid->SetGrid(180*0.25,0,2*TMath::Pi(),0.25*600,-10,10); // implement check on minEta, ...
	      //grid->GridInfo();
	      grid->SetSeed(5.0);
	      grid->SetBkgPtCut(0.0);
	      grid->SetCone(1.0);
	      grid->SetRMaxNN(1.0);
	      grid->Init();
	  
	      //b->Fill(grid);
	      b->Fill(grid,(Double_t) 0.0);
	      
	      if (nEvents<10)
		{

		  TCanvas *c0=new TCanvas("c0","#0",800,600);
		  c0->SetLogz();
		  grid->hGrid->DrawCopy("lego2");

		  cout<<endl;
		  cout<<(Int_t) grid->GetNCells()<<" cells in grid acceptance + cuts."<<endl;
		}
	      
	      delete grid;
	    }
	  else
	    {
	      ktGrid *grid=0;
	      b->Fill(grid,0.0);
	    }
	}
    }
  else
    cout<<" Error :  MC bkg. not initialized !"<<endl;
  
  if (b->GetInitOK())
    {
      b->HistNorm(nEvents);
      TCanvas *c1=new TCanvas("c1","Canvas #1",800,800);
      c1->Divide(2,2);
      c1->cd(1);b->heta->DrawCopy();
      c1->cd(2);b->hpt->SetMinimum(0.001);
      //TF1 *mfit=new TF1("mfit","[0]*exp(-x/[1])",0.1,5);mfit->SetParameters(1,0.3);
      //b->hpt->Fit(mfit,"QR");
      b->hpt->DrawCopy();gPad->SetLogy();

      c1->cd(3);b->hphi->DrawCopy();
      if (b->Getv2()>0)
	{
	  c1->cd(4);b->hRP->DrawCopy();
	}
    }

  cout<<" # MC events created = "<<nEvents<<endl;
  cout<<endl;

  delete b;
}
