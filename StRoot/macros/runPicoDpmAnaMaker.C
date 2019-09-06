
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

#include "StPicoHFMyAnaMaker/StPicoHFMyAnaMaker.h"

#include "macros/loadSharedHFLibraries.C"

#include <iostream>
#include <ctime>
#include <cstdio>

#include "StPicoDpmAnaMaker/StPicoDpmAnaMaker.h" //kvapil

#include "StRefMultCorr/StRefMultCorr.h"
#include "StRefMultCorr/CentralityMaker.h"

using namespace std;

#else
class StChain;
#endif

StChain *chain;

void runPicoDpmAnaMaker(const Char_t *inputFile="test.list", const Char_t *outputFile="outputBaseName",
			 const unsigned int makerMode = 0 /*kAnalyze*/,
			 const Char_t *badRunListFileName = "picoList_bad_MB.list", const Char_t *treeName = "picoHFtree",
			 const Char_t *productionBasePath = "/star/data100/reco/AuAu_200_production_2016/ReversedFullField/P16ij/2016",
			 const unsigned int decayChannel = 0 /* kChannel0 */) {
  // -- Check STAR Library. Please set SL_version to the original star library used in the production
  //    from http://www.star.bnl.gov/devcgi/dbProdOptionRetrv.pl
  string SL_version = "SL16j"; //new: SL16d -> SL16j
  string env_SL = getenv ("STAR");
  if (env_SL.find(SL_version)==string::npos) {
      cout<<"Environment Star Library does not match the requested library in runPicoHFMyAnaMaker.C. Exiting..."<<endl;
      exit(1);
  }

  //Int_t nEvents = 10000000;
  //Int_t nEvents = 1000;

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
//  StPicoDstMaker* picoDstMaker = new StPicoDstMaker(2, sInputFile, "picoDstMaker"); //for local testing only
  StPicoDpmAnaMaker* picoDpmAnaMaker = new StPicoDpmAnaMaker("picoDpmAnaMaker", picoDstMaker, outputFile, sInputListHF);
  picoDpmAnaMaker->setMakerMode(makerMode);
  picoDpmAnaMaker->setDecayChannel(StPicoDpmAnaMaker::kChannel1);//kvapil
  picoDpmAnaMaker->setTreeName(treeName);
  //picoDpmAnaMaker->setMcMode(mcMode); commented kvapil

  StHFCuts* hfCuts = new StHFCuts("hfBaseCuts");
  picoDpmAnaMaker->setHFBaseCuts(hfCuts);

  // ---------------------------------------------------
  // -- Set Base cuts for HF analysis

  // -- File name of bad run list
  hfCuts->setBadRunListFileName(badRunListFileName);

  // -- ADD USER CUTS HERE ----------------------------

  hfCuts->setCutVzMax(6.);
  hfCuts->setCutVzVpdVzMax(3.);
/* SL16d triggers
  hfCuts->addTriggerId(450050);    // vpdmb-5-p-nobsmd-hlt
  hfCuts->addTriggerId(450060);    // vpdmb-5-p-nobsmd-hlt
  hfCuts->addTriggerId(450005);    // vpdmb-5-p-nobsmd
  hfCuts->addTriggerId(450015);    // vpdmb-5-p-nobsmd
  hfCuts->addTriggerId(450025);    // vpdmb-5-p-nobsmd
*/
	//SL16j triggers
  //hfCuts->addTriggerId(520802);    // VPDMB-5-p-hlt
  //hfCuts->addTriggerId(520812);    // VPDMB-5-p-hlt
  //hfCuts->addTriggerId(520822);    // VPDMB-5-p-hlt
  //hfCuts->addTriggerId(520832);    // VPDMB-5-p-hlt
  //hfCuts->addTriggerId(520842);    // VPDMB-5-p-hlt


	hfCuts->addTriggerId(520001);    // VPDMB-5-p-sst
  hfCuts->addTriggerId(520011);    // VPDMB-5-p-sst
  hfCuts->addTriggerId(520021);    // VPDMB-5-p-sst
  hfCuts->addTriggerId(520031);    // VPDMB-5-p-sst
  hfCuts->addTriggerId(520041);    // VPDMB-5-p-sst
  hfCuts->addTriggerId(520051);    // VPDMB-5-p-sst

  hfCuts->addTriggerId(570002);    // VPDMB-5-nosst (production 2, nosst stream)
  hfCuts->addTriggerId(570001);    // VPDMB-5-sst (production 2, sst stream )

  hfCuts->setHFTinputsOrPIDefficiency(0); //0 - HFT inputs; 1 - PID efficiency; 2 - tracking efficiency sys. err.

  hfCuts->setCutNHitsFitMin(15); //kvapil 20 to 15, for candidates
	hfCuts->setCutNHitsFitMinHist(20); //for histograms, Vanek
  hfCuts->setCutRequireHFT(true);
  
  //check all cuts in StHFCuts after HFMaker update - see GitHub for reference
	hfCuts->setCutDca(1.5); //for QA, see createQA() in StPicoDpmAnaMaker.cxx
	hfCuts->setCutDcaXy(1.);
	hfCuts->setCutDcaZ(1.);

  hfCuts->setCutDcaMin(0.009,StHFCuts::kPion); //federic 1aug2016
  //hfCuts->setCutDcaMin(0.01,StHFCuts::kKaon); //federic 1aug2016
  hfCuts->setCutDcaMin(0.007,StHFCuts::kKaon); //federic 3aug2016

  hfCuts->setCutNHitsFitnHitsMax(0.52);

  // ---------------------------------------------------

  // -- Channel0
  //picoHFMyAnaMaker->setDecayMode(StPicoHFEvent::kTwoParticleDecay);
  picoDpmAnaMaker->setDecayMode(StPicoHFEvent::kThreeParticleDecay); //kvapil

  // -- ADD USER CUTS HERE ----------------------------

	hfCuts->setCutEta(1.);

	hfCuts->setCutTripletdV0Max(0.022);

	//-----------SECONDARY TRIPLET CUTS----------------------------
	float dcaDaughters12Max, dcaDaughters23Max, dcaDaughters31Max;
  float decayLengthMin, decayLengthMax;
	float cosThetaMin, massMin, massMax;

  dcaDaughters12Max = 0.009;
	dcaDaughters23Max = 0.009;
	dcaDaughters31Max = 0.009;

	decayLengthMin = 0.003;
	decayLengthMax = 0.2;

	cosThetaMin = 0.997;
	massMin = 1.7;
	massMax = 2.1;

	hfCuts->setCutSecondaryTriplet(dcaDaughters12Max, dcaDaughters23Max, dcaDaughters31Max,
				 decayLengthMin, decayLengthMax,
				 cosThetaMin, massMin, massMax);
  // --- Lomnitz cuts to remove noise from ghosting

  //Single track pt
  hfCuts->setCutPtRange(0.3,50.0,StHFCuts::kPion); //used in candidates analysis
  hfCuts->setCutPtRange(0.3,50.0,StHFCuts::kKaon); //changed to 0.3

	hfCuts->setCutPtQA(0.3); //p_T used in createQA() in StPicoDpmAnaMaker.cxx
  hfCuts->setCutPtRangeQA(0.5, 4.0); //pT range for PID efficiency
  //TPC setters
  hfCuts->setCutTPCNSigmaPion(1.0); //possibly close as possible, e.g. 1.0
  hfCuts->setCutTPCNSigmaKaon(1.0); //all changed to 1.0
	hfCuts->setCutTPCNSigmaProton(1.0); //for QA and analysis

	hfCuts->setCutTPCNSigmaHadronHist(1.0, 1); //1 = pion, not used for QA
	hfCuts->setCutTPCNSigmaHadronHist(1.0, 2); //2 = kaon
  //TOF setters, need to set pt range as well
  hfCuts->setCutTOFDeltaOneOverBeta(0.03, StHFCuts::kKaon); //changed to 0.03 (same cut as used in analysis)
  hfCuts->setCutPtotRangeHybridTOF(0.3,50.0,StHFCuts::kKaon); //changed lower boundary from 0.5 to 0.3 (same cut as used in analysis)
  hfCuts->setCutTOFDeltaOneOverBeta(0.03, StHFCuts::kPion); //changed to 0.03 (same cut as used in analysis)
  hfCuts->setCutPtotRangeHybridTOF(0.3,50.0,StHFCuts::kPion); //changed lower boundary from 0.5 to 0.3 (same cut as used in analysis)


  //pion and kaon pair cuts
  hfCuts->setCutPiPairInvMassInterval(0.4, 0.6); //inv. mass interval for K0s -> pi+pi
  hfCuts->setCutKPairInvMassInterval(1.0, 1.04 ); //inv. mass interval for phi -> K+K


  // set refmultCorr
//  cout<<"test"<<endl;
  StRefMultCorr* grefmultCorrUtil = CentralityMaker::instance()->getgRefMultCorr_P16id(); //new StRefMultCorr, info about Run16, SL16j in the same file as for Run14, SL16d
  picoDpmAnaMaker->setRefMutCorr(grefmultCorrUtil);
  //cout<<"test2"<<endl;
  // ========================================================================================

  // ========================================================================================

  //clock_t start = clock(); // getting starting time
  chain->Init();
  cout << "chain->Init();" << endl;
  int nEvents = picoDstMaker->chain()->GetEntries();
  cout << " Total entries = " << nEvents << endl;
  //if(nEvents>total) nEvents = total;

  for (Int_t i=0; i<nEvents; i++) {
    if(i%100000==0) //orig 10000
      cout << "Working on eventNumber " << i << endl;

    chain->Clear();
    int iret = chain->Make(i);

    if (iret) { cout << "Bad return code!" << iret << endl; break;}

    //total++;
  }

  cout << "****************************************** " << endl;
  cout << "Work done... now its time to close up shop!"<< endl;
  cout << "****************************************** " << endl;
  chain->Finish();
  //double duration = (double) (clock() - start) / (double) CLOCKS_PER_SEC;
  cout << "****************************************** " << endl;
  cout << "total number of events  " << nEvents << endl;
  cout << "****************************************** " << endl;
 // cout << "Time needed " << duration << " s" << endl;
  cout << "****************************************** " << endl;

  delete chain;

  // -- clean up if in read mode
  if (makerMode == StPicoHFMaker::kRead)
    gSystem->Exec(Form("rm -f %s", sInputFile.Data()));
}

