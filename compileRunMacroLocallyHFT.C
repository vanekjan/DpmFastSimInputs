#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"

void compileRunMacroLocallyHFT(const char* name = "runPicoDpmAnaMakerHFT.C"){ //orig. runPicoHFMyAnaMaker.C
  Long_t ret = gROOT->ProcessLine(".L StRoot/macros/loadSharedHFLibraries.C");
  cout << ret << endl;
  loadSharedHFLibraries();
  gSystem->CompileMacro(Form("StRoot/macros/%s", name), "kfc"); //compile macro locally without running it or loading it
}
