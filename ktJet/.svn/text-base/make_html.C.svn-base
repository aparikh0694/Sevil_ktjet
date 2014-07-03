void make_html()
{

   // ktJet lib
  if (gClassTable->GetID("ktEvent") < 0) {
    cout<<"Load ktJet lib ..."<<endl;
    gSystem->Load("libKtJet.so");
  }
 
  cout<<"Make class documetation ROOT-Style ..."<<endl;
  THtml html;
  html.SetSourceDir("./");
  //html.SetAuthorTag("Joern Putschke");

  html.MakeClass("ktJetCell");
  html.MakeClass("ktGrid");
  html.MakeClass("ktJet");
  html.MakeClass("ktNN");
  html.MakeClass("ktJetCellPair");
  html.MakeClass("ktPseudoJet");
  html.MakeClass("ktParton");
  html.MakeClass("ktMuEvent");
  html.MakeClass("ktMuJet");
  html.MakeClass("ktFastJet");
  html.MakeClass("ktTrigger");
  html.MakeClass("ktTriggerPatch");
  html.MakeClass("ktPythia8");
  html.MakeClass("ktMCBkg");
  html.MakeClass("ktPy8Event");
  html.MakeClass("ktStarPico");
  html.MakeClass("ktPID");
  html.MakeClass("ktMuFastJet");
  html.MakeClass("ktParticle");
  html.MakeClass("ktAna");
  html.MakeClass("ktJetQuench");
  html.MakeIndex("kt*");

  //THtml htmlexp;
  //htmlexp.SetSourceDir("./");
  //cout<<"Convert examples ..."<<endl;
  //htmlexp.Convert("kt_test.C","kt-jetfinder example");
  //htmlexp.Convert("kt_test_pythia.C","kt-jetfinder example running PYTHIA6");
  //htmlexp.Convert("kt_test_pythia8.C","kt-jetfinder example running PYTHIA8");
  //htmlexp.Convert("test_trigger.C","Sliding patch trigger");
   
  //html.MakeIndex("kt*");

  cout<<"Done."<<endl;

}
