#include "TROOT.h"
#include "TRint.h"

int Error; 

extern void InitGui();  
VoidFuncPtr_t initfuncs[] = { InitGui, 0 };

TROOT root("Rint","The ROOT Interactive Interface", initfuncs);

int main(int argc, char **argv)
{

  TRint *theApp = new TRint("ROOT example",&argc, argv,initfuncs, 0,0);
  theApp->Run();

  //TApplication theApp("App", &argc, argv);
  //theApp.Run();

  return(0);

}

