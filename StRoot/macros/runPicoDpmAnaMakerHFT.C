
/* **************************************************
 *   Run StPicoHFMyAnaMaker in different modes
 * --------------------------------------------------
 * run as :
 *  root -l -b -q StRoot/macros/loadSharedHFLibraries.C StRoot/macros/runPicoHFMyAnaMaker.C++
 *   or
 *  root -l -b -q StRoot/macros/runPicoHFMyAnaMaker.C
 *
 * --------------------------------------------------
 *  - Different modes to use the  class
 *    - StPicoHFMaker::kAnalyze - don't write candidate trees, just fill histograms
 *        inputFile : fileList of PicoDst files or single picoDst file
 *        outputFile: baseName for outfile
 *    - StPicoHFMaker::kWrite   - write candidate trees
 *        inputFile : path to single picoDist file
 *        outputFile: baseName for outfile
 *    - StPicoHFMaker::kRead    - read candidate trees and fill histograms
 *        inputFile : fileList of PicoDst files
 *        outputFile: baseName for outfile
 *
 * --------------------------------------------------
 *  Authors:  Xin Dong        (xdong@lbl.gov)
 *            Michael Lomnitz (mrlomnitz@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *            Jochen Thaeder  (jmthader@lbl.gov)
 *
 * **************************************************
 */

#ifndef __CINT__
#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"

#include "StMaker.h"
#include "StChain.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoHFMaker/StPicoHFEvent.h"
#include "StPicoHFMaker/StHFCuts.h"

#include "macros/loadSharedHFLibraries.C"

#include <iostream>
#include <ctime>
#include <cstdio>

#include "StPicoDpmAnaMaker/StPicoDpmAnaMaker.h"

#include "StRefMultCorr/StRefMultCorr.h"
#include "StRefMultCorr/CentralityMaker.h"

using namespace std;

#else
class StChain;
#endif

StChain *chain;

void runPicoDpmAnaMakerHFT(const Char_t *inputFile="test.list", const Char_t *outputFile="outputBaseName",
    const unsigned int makerMode = 0 /*kAnalyze*/,
    const Char_t *badRunListFileName = "picoList_bad_MB.list", const Char_t *treeName = "picoHFtree",
    const Char_t *productionBasePath = "/star/data100/reco/AuAu_200_production_2016/ReversedFullField/P16ij/2016",
    const unsigned int decayChannel = 0 /* kChannel0 */) {
  // -- Check STAR Library. Please set SL_version to the original star library used in the production
  //    from http://www.star.bnl.gov/devcgi/dbProdOptionRetrv.pl
  string SL_version = "SL17d"; //originally SL16j, not available any more
  string env_SL = getenv ("STAR");
  if (env_SL.find(SL_version)==string::npos) {
    cout<<"Environment Star Library does not match the requested library in runPicoHFMyAnaMaker.C. Exiting..."<<endl;
    exit(1);
  }

#ifdef __CINT__
  gROOT->LoadMacro("loadSharedHFLibraries.C");
  loadSharedHFLibraries();
#endif


  chain = new StChain();

  // ========================================================================================
  //makerMode    = StPicoHFMaker::kAnalyze;
  // ========================================================================================

  cout << "Maker Mode    " << makerMode << endl;
  cout << "TreeName      " << treeName << endl;
  cout << "Decay Channel " << decayChannel << endl;

  TString sInputFile(inputFile);
  TString sInputListHF("");
  TString sProductionBasePath(productionBasePath);
  TString sTreeName(treeName);

  if (makerMode == StPicoHFMaker::kAnalyze) {
    if (!sInputFile.Contains(".list") && !sInputFile.Contains("picoDst.root")) {
      cout << "No input list or picoDst root file provided! Exiting..." << endl;
      exit(1);
    }
  }
  else if (makerMode == StPicoHFMaker::kWrite) {
    if (!sInputFile.Contains("picoDst.root")) {
      cout << "No input picoDst root file provided! Exiting..." << endl;
      exit(1);
    }
  }
  else if (makerMode == StPicoHFMaker::kRead) {
    if (!sInputFile.Contains(".list")) {
      cout << "No input list provided! Exiting..." << endl;
      exit(1);
    }

    // -- prepare filelist for picoDst from trees
    sInputListHF = sInputFile;
    sInputFile = "tmpPico.list";

    TString command = "sed 's|" + sTreeName + ".root|picoDst.root|g' " + sInputListHF + " > " + sInputFile;
    cout << "COMMAND : " << command << endl;
    gSystem->Exec(command.Data());

    command = "sed -i 's|^.*" + sTreeName + "|" + sProductionBasePath + "|g' " + sInputFile; // + " > " + sInputFile;
    cout << "COMMAND : " << command << endl;
    gSystem->Exec(command.Data());
  }
  else {
    cout << "Unknown makerMode! Exiting..." << endl;
    exit(1);
  }

  StPicoDstMaker* picoDstMaker = new StPicoDstMaker(StPicoDstMaker::IoRead, sInputFile, "picoDstMaker"); //SL16j: See StRoot/StPicoDstMaker/StpicodstMaker.h: 28: enum PicoIoMode {IoWrite=1, IoRead=2};
  StPicoDpmAnaMaker* picoDpmAnaMaker = new StPicoDpmAnaMaker("picoDpmAnaMaker", picoDstMaker, outputFile, sInputListHF);
  picoDpmAnaMaker->setMakerMode(makerMode);
  picoDpmAnaMaker->setDecayChannel(StPicoDpmAnaMaker::kChannel1);
  picoDpmAnaMaker->setTreeName(treeName);

  StHFCuts* hfCuts = new StHFCuts("hfBaseCuts");
  picoDpmAnaMaker->setHFBaseCuts(hfCuts);

  // ---------------------------------------------------
  // -- Set Base cuts for HF analysis

  // -- File name of bad run list
  hfCuts->setBadRunListFileName(badRunListFileName);

  // -- ADD USER CUTS HERE ----------------------------

  hfCuts->setCutVzMax(6.);
  hfCuts->setCutVzVpdVzMax(3.);


  hfCuts->addTriggerId(520001);    // VPDMB-5-p-sst
  hfCuts->addTriggerId(520011);    // VPDMB-5-p-sst
  hfCuts->addTriggerId(520021);    // VPDMB-5-p-sst
  hfCuts->addTriggerId(520031);    // VPDMB-5-p-sst
  hfCuts->addTriggerId(520041);    // VPDMB-5-p-sst
  hfCuts->addTriggerId(520051);    // VPDMB-5-p-sst

  hfCuts->addTriggerId(570002);    // VPDMB-5-nosst (production 2, nosst stream)
  hfCuts->addTriggerId(570001);    // VPDMB-5-sst (production 2, sst stream )

  hfCuts->setHFTinputsOrPIDefficiency(0); //0 - fast-sim inputs; 1 - PID efficiency; 2 - tracking efficiency sys. err.

  hfCuts->setHFTratiosOrDCAdistributions(1); //0 - DCA distributions for FastSim, 1 - HFT ratios for FastSim

  hfCuts->setCutNHitsFitMinHist(20); //for histograms, Vanek
  hfCuts->setCutRequireHFT(true);

  //check all cuts in StHFCuts after HFMaker update - see GitHub for reference
  hfCuts->setCutDca(1.5); //HFT ratios for FastSim inputs, see createQA() in StPicoDpmAnaMaker.cxx
  hfCuts->setCutDcaXy(1.0); 
  hfCuts->setCutDcaZ(1.0);

  hfCuts->setCutNHitsFitnHitsMax(0.52);

  // ---------------------------------------------------

  // -- Channel0
  picoDpmAnaMaker->setDecayMode(StPicoHFEvent::kThreeParticleDecay);

  // -- ADD USER CUTS HERE ----------------------------

  hfCuts->setCutEta(1.);


  //Single track pt
  hfCuts->setCutPtRange(0.3,50.0,StHFCuts::kPion); //used in candidates analysis
  hfCuts->setCutPtRange(0.3,50.0,StHFCuts::kKaon); //changed to 0.3

  hfCuts->setCutPtQA(0.3); //p_T used in createQA() in StPicoDpmAnaMaker.cxx

  //TPC setters
  hfCuts->setCutTPCNSigmaPion(3.0); //same as in analysis
  hfCuts->setCutTPCNSigmaKaon(2.0); //same as in analysis

  //TOF setters, need to set pt range as well
  hfCuts->setCutTOFDeltaOneOverBeta(0.03, StHFCuts::kKaon); //changed to 0.03 (same cut as used in analysis)
  hfCuts->setCutPtotRangeHybridTOF(0.3,50.0,StHFCuts::kKaon); //changed lower boundary from 0.5 to 0.3 (same cut as used in analysis)
  hfCuts->setCutTOFDeltaOneOverBeta(0.03, StHFCuts::kPion); //changed to 0.03 (same cut as used in analysis)
  hfCuts->setCutPtotRangeHybridTOF(0.3,50.0,StHFCuts::kPion); //changed lower boundary from 0.5 to 0.3 (same cut as used in analysis)

  // set refmultCorr
  StRefMultCorr* grefmultCorrUtil = CentralityMaker::instance()->getgRefMultCorr_P16id(); //new StRefMultCorr, info about Run16, SL16j in the same file as for Run14, SL16d
  picoDpmAnaMaker->setRefMutCorr(grefmultCorrUtil);
  // ========================================================================================


  chain->Init();
  cout << "chain->Init();" << endl;
  int nEvents = picoDstMaker->chain()->GetEntries();
  cout << " Total entries = " << nEvents << endl;

  for (Int_t i=0; i<nEvents; i++) {
    if(i%100000==0)
      cout << "Working on eventNumber " << i << endl;

    chain->Clear();
    int iret = chain->Make(i);

    if (iret) { cout << "Bad return code!" << iret << endl; break;}

  }

  cout << "****************************************** " << endl;
  cout << "Work done... now its time to close up shop!"<< endl;
  cout << "****************************************** " << endl;
  chain->Finish();
  cout << "****************************************** " << endl;
  cout << "total number of events  " << nEvents << endl;
  cout << "****************************************** " << endl;

  delete chain;

  // -- clean up if in read mode
  if (makerMode == StPicoHFMaker::kRead)
    gSystem->Exec(Form("rm -f %s", sInputFile.Data()));
}

