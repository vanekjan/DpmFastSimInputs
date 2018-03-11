#include "StPicoDpmAnaMaker.h"
//#include "StPicoHFMaker/StHFCuts.h"
#include <iostream>

ClassImp(StPicoDpmAnaMaker)

using namespace std;

// _________________________________________________________
StPicoDpmAnaMaker::StPicoDpmAnaMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName,
				       char const* inputHFListHFtree = "") :
  StPicoHFMaker(name, picoMaker, outputBaseFileName, inputHFListHFtree),
  mDecayChannel(kChannel1), mRefmultCorrUtil(NULL),mOutFileBaseName(outputBaseFileName){


  // constructor
}

// _________________________________________________________
StPicoDpmAnaMaker::~StPicoDpmAnaMaker() {
  // destructor
}

// _________________________________________________________
int StPicoDpmAnaMaker::InitHF() {
  // -- INITIALIZE USER HISTOGRAMS ETC HERE -------------------
  //    add them to the output list mOutList which is automatically written

  // EXAMPLE //  mOutList->Add(new TH1F(...));
  // EXAMPLE //  TH1F* hist = static_cast<TH1F*>(mOutList->Last());
    cout<<"InitHF1"<<endl;
	  mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");
	  histoInit(mOutFileBaseName, true); //for createQA()
    cout<<"InitHF2"<<endl;
//	cout<<"RefMult1"<<endl;
		mRefmultCorrUtil->setVzForWeight(6, -6.0, 6.0);
    cout<<"RefMult1"<<endl;
    mRefmultCorrUtil->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_VpdnoVtx_Vpd5_Run16.txt"); //for new StRefMultCorr, Run16, SL16j
  	cout<<"RefMult2"<<endl;

   // -------------- USER VARIABLES -------------------------

	mRunNumber = 0;
  return kStOK;
}

// _________________________________________________________
void StPicoDpmAnaMaker::ClearHF(Option_t *opt="") {
  return;
}

// _________________________________________________________
int StPicoDpmAnaMaker::FinishHF() {
   if( isMakerMode() != StPicoHFMaker::kWrite )
   //ntp_Dmeson->Write(); //for candidates - probably not needed for QA
//	cout<<"FinishHF1"<<endl;
   closeFile(); // for QA
//	cout<<"FinishHF2"<<endl;
  return kStOK;
}

// _________________________________________________________
int StPicoDpmAnaMaker::MakeHF() {
  // -- process event
  //    ADD YOUR PROCESSING CODE HERE
  //    ... it is usefull to use the methods below
  //     - createCandidates()
  //     - analyzeCandidates()
//	cout<<"start"<<endl;
  std::clock_t start1 = std::clock();//kvapil
  if (isMakerMode() == StPicoHFMaker::kWrite) {
    createCandidates();
  }
  else if (isMakerMode() == StPicoHFMaker::kRead) {
    // -- the reading back of the perviously written trees happens in the background
    analyzeCandidates();
  }
  else if (isMakerMode() == StPicoHFMaker::kAnalyze) {
    //createCandidates();
    //analyzeCandidates();
    createQA();
  }


  return kStOK;
}
//_________________________________________________________
int StPicoDpmAnaMaker::createQA(){
       //int const currentRun = mPicoHFEvent->runId();
       //if(currentRun != mRunNumber)
     //  {
       //mRunNumber = currentRun;
      mRefmultCorrUtil->init(mPicoDst->event()->runId());
      if (!mRefmultCorrUtil){
         LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
         return kStWarn;
      }
      if (mRefmultCorrUtil->isBadRun(mPicoDst->event()->runId())) return kStOK;
//      cout<<"Q1"<<endl;
//      mRefmultCorrUtil->setVzForWeight(6, -6.0, 6.0);
//      mRefmultCorrUtil->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14_P16id.txt");
//	  cout<<"Q2"<<endl;
      mRefmultCorrUtil->initEvent(mPicoDst->event()->grefMult(), mPrimVtx.z(), mPicoDst->event()->ZDCx()) ;

       int const centrality = mRefmultCorrUtil->getCentralityBin9();
       const double reweight = mRefmultCorrUtil->getWeight();
       const double refmultCor = mRefmultCorrUtil->getRefMultCorr();
       //mHists->addCent(refmultCor, centrality, reweight, pVtx.z());
       UInt_t nTracks = mPicoDst->numberOfTracks();

       for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack)
       {
          StPicoTrack const* trk = mPicoDst->track(iTrack);
          if (!trk) continue;
          //StPhysicalHelixD helix = trk->helix(); //SL16d
		  		StPhysicalHelixD helix = trk->helix(mPicoDst->event()->bField()); //SL16j, Vanek
          float dca = float(helix.geometricSignedDistance(mPrimVtx));
          StThreeVectorF momentum = trk->gMom(mPrimVtx, mPicoDst->event()->bField());

         	// if (!isGoodQaTrack(trk, momentum, dca)) continue; pt, nhits, pseudorap
					if (!(mHFCuts->hasGoodPtQA(trk))) continue;
					if (!(mHFCuts->hasGoodNHitsFitMinHist(trk))) continue;
          if (!(mHFCuts->hasGoodNHitsFitnHitsMax(trk))) continue; //nHitsFit/nHitsMax
        	if (!(mHFCuts->hasGoodEta(momentum))) continue;

          StThreeVectorF dcaPoint = helix.at(helix.pathLength(mPrimVtx.x(), mPrimVtx.y()));
          float dcaZ = dcaPoint.z() - mPrimVtx.z();
          double dcaXy = helix.geometricSignedDistance(mPrimVtx.x(), mPrimVtx.y());

          bool tpcPion = false;
          bool tpcKaon = false;
          bool tpcProton = false;
		  		if(mHFCuts->hasGoodTPCnSigmaPion(trk)) tpcPion = true;
		  		if(mHFCuts->hasGoodTPCnSigmaKaon(trk)) tpcKaon = true;
	      	if(mHFCuts->hasGoodTPCnSigmaProton(trk)) tpcProton = true;
          //float hBeta = mHFCuts->getTofBetaBase(trk); //SL16d
		  		float hBeta = mHFCuts->getTofBetaBase(trk, mPicoDst->event()->bField()); //SL16j, Vanek
          bool hTofAvailable = !isnan(hBeta) && hBeta > 0;

          bool tofPion = false;
          bool tofKaon = false;
          bool tofProton = false;

	      	if(fabs(1./hBeta - sqrt(1+M_PION_PLUS*M_PION_PLUS/(momentum.mag()*momentum.mag())))<=mHFCuts->getCutTOFDeltaOneOverBeta(StHFCuts::kPion)) tofPion = true;
          if(fabs(1./hBeta - sqrt(1+M_KAON_PLUS*M_KAON_PLUS/(momentum.mag()*momentum.mag())))<=mHFCuts->getCutTOFDeltaOneOverBeta(StHFCuts::kKaon)) tofKaon = true;
		  		if(fabs(1./hBeta - sqrt(1+M_PROTON*M_PROTON/(momentum.mag()*momentum.mag())))) tofProton = true;

          bool goodPion = (hTofAvailable && tofPion && tpcPion) || (!hTofAvailable && tpcPion);//Always require TPC
          bool goodKaon = (hTofAvailable && tofKaon && tpcKaon) || (!hTofAvailable && tpcKaon);
          bool goodProton = (hTofAvailable && tofProton && tpcProton) || (!hTofAvailable && tpcProton);

          if (trk  && fabs(dca) < mHFCuts->cutDca() && trk->isHFTTrack() && (goodPion || goodKaon || goodProton)){
             addDcaPtCent(dca, dcaXy, dcaZ, goodPion, goodKaon, goodProton, momentum.perp(), centrality, momentum.pseudoRapidity(), momentum.phi(), mPrimVtx.z()); //add Dca distribution
          }
          if (trk  && fabs(dca) < mHFCuts->cutDca() && (goodPion || goodKaon || goodProton)){
             //std::cout<<"1: "<<goodPion<<" "<< goodKaon<<" "<<  goodProton<<" "<<  momentum.perp()<<" "<<  centrality<<" "<<  momentum.pseudoRapidity()<<" "<<  momentum.phi()<<" "<<  mPrimVtx.z()<<std::endl;
             addTpcDenom1(goodPion, goodKaon, goodProton, momentum.perp(), centrality, momentum.pseudoRapidity(), momentum.phi(), mPrimVtx.z()); //Dca cut on 1.5cm, add Tpc Denominator
          }
/* old, Kvapil
          if (trk && fabs(dca) < mHFCuts->cutDca() && trk->isHFTTrack() && (goodPion || goodKaon || goodProton) && fabs(dcaXy) < mHFCuts->cutDcaXy() && fabs(dcaZ) < mHFCuts->cutDcaZ()){
             addHFTNumer1(goodPion, goodKaon, goodProton, momentum.perp(), centrality,  momentum.pseudoRapidity(), momentum.phi(), mPrimVtx.z()); //Dca cut on 1.5cm, add HFT Numerator
          }
*/
          //new version, Vanek 03/10/18
          if (trk && fabs(dca) < mHFCuts->cutDca() && trk->isHFTTrack() && (goodPion || goodKaon || goodProton)){
             addHFTNumer1(goodPion, goodKaon, goodProton, momentum.perp(), centrality,  momentum.pseudoRapidity(), momentum.phi(), mPrimVtx.z()); //Dca cut on 1.5cm, add HFT Numerator
          }
       } // .. end tracks loop
   return 0;
}

// _________________________________________________________
int StPicoDpmAnaMaker::createCandidates() {
  // Creating candidates for D+- 3 body decay
  // D+- -> K+2Pi decay

  for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1) {
    StPicoTrack const *pion1 = mPicoDst->track(mIdxPicoPions[idxPion1]);
    // -- Pion selection

    for (unsigned short idxPion2 = idxPion1+1; idxPion2 < mIdxPicoPions.size(); ++idxPion2) {
      StPicoTrack const *pion2 = mPicoDst->track(mIdxPicoPions[idxPion2]);
      // -- Pion selection
      if ( !isCloseTracks(pion1,pion2,mPrimVtx, mBField)) continue;

      for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon) {
        StPicoTrack const *kaon = mPicoDst->track(mIdxPicoKaons[idxKaon]);
        // -- Kaon selection
        // -- TOF
		//if( !mHFCuts->isHybridTOFHadron(kaon, mHFCuts->getTofBetaBase(kaon), StHFCuts::kKaon) ) continue; //SL16d
        if( !mHFCuts->isHybridTOFHadron(kaon, mHFCuts->getTofBetaBase(kaon, mPicoDst->event()->bField()), StHFCuts::kKaon) ) continue; //SL16j, Vanek
        if (mIdxPicoKaons[idxKaon] == mIdxPicoPions[idxPion1]|| mIdxPicoKaons[idxKaon] == mIdxPicoPions[idxPion2] || mIdxPicoPions[idxPion1] == mIdxPicoPions[idxPion2]) continue;
				if ( !isCloseTracks(pion1,kaon,mPrimVtx, mBField)) continue;
				if ( !isCloseTracks(kaon,pion2,mPrimVtx, mBField)) continue;
        // -- Making triplet
        StHFTriplet triplet(pion1,pion2,kaon,mHFCuts->getHypotheticalMass(StHFCuts::kPion),mHFCuts->getHypotheticalMass(StHFCuts::kPion),mHFCuts->getHypotheticalMass(StHFCuts::kKaon), mIdxPicoPions[idxPion1],mIdxPicoPions[idxPion2],mIdxPicoKaons[idxKaon], mPrimVtx, mBField);
        if(mHFCuts->hasGoodTripletdV0Max(triplet)) continue;
        if (!mHFCuts->isGoodSecondaryVertexTriplet(triplet)) continue;
        mPicoHFEvent->addHFSecondaryVertexTriplet(&triplet);

      }  // for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon)
    } // for (unsigned short idxPion2 = idxPion1+1; idxPion2 < mIdxPicoPions.size(); ++idxPion2)
  } // for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1)
 return kStOK;
}

// _________________________________________________________
int StPicoDpmAnaMaker::analyzeCandidates() {

  // --- Analyze previously constructed candidates and output to ntuple
  // -- Decay channel1
  TClonesArray const * aCandidates= mPicoHFEvent->aHFSecondaryVertices();
//	cout<<"RefMultInit3"<<endl;
	  mRefmultCorrUtil->init(mPicoDst->event()->runId());
      if (!mRefmultCorrUtil){
         LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
         return kStWarn;
      }
//	cout<<"RefMultInit4"<<endl;

      if (mRefmultCorrUtil->isBadRun(mPicoDst->event()->runId())) return kStOK;
//      cout<<"A1"<<endl;
//      mRefmultCorrUtil->setVzForWeight(6, -6.0, 6.0); //monve to InitHF???
//      mRefmultCorrUtil->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14_P16id.txt");
//	  cout<<"A2"<<endl;
      mRefmultCorrUtil->initEvent(mPicoDst->event()->grefMult(), mPrimVtx.z(), mPicoDst->event()->ZDCx()) ;

       int const centrality = mRefmultCorrUtil->getCentralityBin9();
       const double reweight = mRefmultCorrUtil->getWeight();
       const double refmultCor = mRefmultCorrUtil->getRefMultCorr();


  if( mPicoHFEvent->nHFSecondaryVertices() >0 ){
    for (unsigned int idx = 0; idx <  mPicoHFEvent->nHFSecondaryVertices(); ++idx) {

      StHFTriplet const* triplet = static_cast<StHFTriplet*>(aCandidates->At(idx));
      StPicoTrack const* pion1 = mPicoDst->track(triplet->particle1Idx());
      StPicoTrack const* pion2 = mPicoDst->track(triplet->particle2Idx());
      StPicoTrack const* kaon = mPicoDst->track(triplet->particle3Idx());

      // Greates distance between tracks
      float const dcaDaughters_12 = triplet->dcaDaughters12();
      float const dcaDaughters_23 = triplet->dcaDaughters23();
      float const dcaDaughters_13 = triplet->dcaDaughters31();
      float dcaMax = dcaDaughters_12 > dcaDaughters_13 ? dcaDaughters_12 : dcaDaughters_13;
      dcaMax = dcaMax > dcaDaughters_23 ? dcaMax : dcaDaughters_23;

      //TOF ---
      float kaonBetaBase = -1;
      float pion1BetaBase = -1;
      float pion2BetaBase = -1;
/*      kaonBetaBase = mHFCuts->getTofBetaBase(kaon); //SL16d
      pion1BetaBase = mHFCuts->getTofBetaBase(pion1);
      pion2BetaBase = mHFCuts->getTofBetaBase(pion2);
*/
	  kaonBetaBase = mHFCuts->getTofBetaBase(kaon, mPicoDst->event()->bField()); //SL16j, Vanek
      pion1BetaBase = mHFCuts->getTofBetaBase(pion1, mPicoDst->event()->bField());
      pion2BetaBase = mHFCuts->getTofBetaBase(pion2, mPicoDst->event()->bField());

/* orig. Kvapil
      float kaonTOFinvbeta = fabs(1. / mHFCuts->getTofBetaBase(kaon) - sqrt(1+M_KAON_PLUS*M_KAON_PLUS/(kaon->gMom(mPrimVtx,mBField).mag()*kaon->gMom(mPrimVtx,mBField).mag())));
      float pion1TOFinvbeta = fabs(1. / mHFCuts->getTofBetaBase(pion1) - sqrt(1+M_PION_PLUS*M_PION_PLUS/(pion1->gMom(mPrimVtx,mBField).mag()*pion1->gMom(mPrimVtx,mBField).mag())));
      float pion2TOFinvbeta = fabs(1. / mHFCuts->getTofBetaBase(pion2) - sqrt(1+M_PION_PLUS*M_PION_PLUS/(pion2->gMom(mPrimVtx,mBField).mag()*pion2->gMom(mPrimVtx,mBField).mag())));
*/

//update Vanek
	  float kaonTOFinvbeta = fabs(1. / kaonBetaBase - sqrt(1+M_KAON_PLUS*M_KAON_PLUS/(kaon->gMom(mPrimVtx,mBField).mag()*kaon->gMom(mPrimVtx,mBField).mag())));
      float pion1TOFinvbeta = fabs(1. / pion1BetaBase - sqrt(1+M_PION_PLUS*M_PION_PLUS/(pion1->gMom(mPrimVtx,mBField).mag()*pion1->gMom(mPrimVtx,mBField).mag())));
      float pion2TOFinvbeta = fabs(1. / pion2BetaBase - sqrt(1+M_PION_PLUS*M_PION_PLUS/(pion2->gMom(mPrimVtx,mBField).mag()*pion2->gMom(mPrimVtx,mBField).mag())));

      // -- Flag D plus and Dminus
      flag = -99;
      if( kaon->charge()<0 && pion1->charge()>0 && pion2->charge()>0 ) flag=0; // -- D+ -> K- + 2pi+
      if( kaon->charge()>0 && pion1->charge()<0 && pion2->charge()<0 ) flag=1; // -- D- -> K+ + 2pi-

      if( kaon->charge()<0 && pion1->charge()>0 && pion2->charge()<0 ) flag=2; // -+-
      if( kaon->charge()<0 && pion1->charge()<0 && pion2->charge()>0 ) flag=2; // --+
      if( kaon->charge()>0 && pion1->charge()>0 && pion2->charge()<0 ) flag=3; // ++-
      if( kaon->charge()>0 && pion1->charge()<0 && pion2->charge()>0 ) flag=3; // +-+

      if( kaon->charge()<0 && pion1->charge()<0 && pion2->charge()<0 ) flag=4; // ---
      if( kaon->charge()>0 && pion1->charge()>0 && pion2->charge()>0 ) flag=5; // +++

      // getting centrality
      /*int const currentRun = mPicoHFEvent->runId();

      if(currentRun != mRunNumber)
      {
          // init a new run
          mRunNumber = currentRun;
          mRefmultCorrUtil->init(mRunNumber);
          mRefmultCorrUtil->setVzForWeight(6, -6.0, 6.0);
          mRefmultCorrUtil->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14.txt");
          for(Int_t i=0;i<6;i++){
            mRefmultCorrUtil->get(i, 0);
          }
      }

      mRefmultCorrUtil->initEvent(mPicoDst->event()->grefMult(), mPrimVtx.z(), mPicoDst->event()->ZDCx()) ;
      int const centrality = mRefmultCorrUtil->getCentralityBin9() ;*/

/*		Moved before for-cycle
      mRefmultCorrUtil->init(mPicoDst->event()->runId());
      if (!mRefmultCorrUtil){
         LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
         return kStWarn;
      }
      if (mRefmultCorrUtil->isBadRun(mPicoDst->event()->runId())) return kStOK;

      mRefmultCorrUtil->setVzForWeight(6, -6.0, 6.0);
      mRefmultCorrUtil->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14_P16id.txt");

      mRefmultCorrUtil->initEvent(mPicoDst->event()->grefMult(), mPrimVtx.z(), mPicoDst->event()->ZDCx()) ;

       int const centrality = mRefmultCorrUtil->getCentralityBin9();
       const double reweight = mRefmultCorrUtil->getWeight();
       const double refmultCor = mRefmultCorrUtil->getRefMultCorr();
*/

      // ---
/*      // Saving to TTree - not needed for QA
	    pi1_runId = mPicoHFEvent->runId();
      pi1_eventId = mPicoHFEvent->eventId();
      pi1_phi = pion1->gMom(mPrimVtx,mBField).phi();
      pi1_eta = pion1->gMom(mPrimVtx,mBField).pseudoRapidity();
      pi1_pt  = pion1->gPt();
      pi1_dca = triplet->particle1Dca();
      pi1_dedx = pion1->dEdx();
      pi1_nSigma = pion1->nSigmaPion();
      pi1_nHitFit = pion1->nHitsFit();
      pi1_nHitdedx = pion1->nHitsDedx();
      pi1_TOFinvbeta = pion1TOFinvbeta;
      pi1_betaBase = pion1BetaBase;

      pi2_runId = mPicoHFEvent->runId();
      pi2_eventId = mPicoHFEvent->eventId();
      pi2_phi = pion2->gMom(mPrimVtx,mBField).phi();
      pi2_eta = pion2->gMom(mPrimVtx,mBField).pseudoRapidity();
      pi2_pt = pion2->gPt();
      pi2_dca = triplet->particle2Dca();
      pi2_dedx = pion2->dEdx();
      pi2_nSigma = pion2->nSigmaPion();
      pi2_nHitFit = pion2->nHitsFit();
      pi2_nHitdedx = pion2->nHitsDedx();
      pi2_TOFinvbeta = pion2TOFinvbeta;
      pi2_betaBase = pion2BetaBase;

      k_runId = mPicoHFEvent->runId();
      k_eventId = mPicoHFEvent->eventId();
      k_phi = kaon->gMom(mPrimVtx,mBField).phi();
      k_eta = kaon->gMom(mPrimVtx,mBField).pseudoRapidity();
      k_pt = kaon->gPt();
      k_dca = triplet->particle3Dca();
      k_dedx = kaon->dEdx();
      k_nSigma = kaon->nSigmaKaon();
      k_nHitFit = kaon->nHitsFit();
      k_nHitdedx = kaon->nHitsDedx();
      k_TOFinvbeta = kaonTOFinvbeta;
      k_betaBase = kaonBetaBase;

      mdcaMax = dcaMax;
      //flag already saved on lines 554 - 564
      primVz = mPrimVtx.z();

      D_theta = triplet->pointingAngle();
      D_decayL = triplet->decayLength();
      D_phi = triplet->phi();
      D_eta = triplet->eta();
      D_pt = sqrt(pow(triplet->px(),2.0)+pow(triplet->py(),2.0));
      D_mass = triplet->m();
      D_dV0Max = triplet->dV0Max();
      mcentrality = centrality;
      refmult = mPicoDst->event()->refMult();
      mrefmultcorr = refmultCor;
      mreweight = reweight;

	  ntp_Dmeson->Fill();

*/


/*

	  int ii=0; //not used now
      float ntVar[58];

      ntVar[ii++] = mPicoHFEvent->runId();
      ntVar[ii++] = mPicoHFEvent->eventId();
      ntVar[ii++] = pion1->gMom(mPrimVtx,mBField).phi();
      ntVar[ii++] = pion1->gMom(mPrimVtx,mBField).pseudoRapidity();
      ntVar[ii++] = pion1->gPt();
      ntVar[ii++] = triplet->particle1Dca();
      ntVar[ii++] = pion1->dEdx();
      ntVar[ii++] = pion1->nSigmaPion();
      ntVar[ii++] = pion1->nHitsFit();
      ntVar[ii++] = pion1->nHitsDedx();
      ntVar[ii++] = pion1TOFinvbeta;
      ntVar[ii++] = pion1BetaBase;

      ntVar[ii++] = mPicoHFEvent->runId();
      ntVar[ii++] = mPicoHFEvent->eventId();
      ntVar[ii++] = pion2->gMom(mPrimVtx,mBField).phi();
      ntVar[ii++] = pion2->gMom(mPrimVtx,mBField).pseudoRapidity();
      ntVar[ii++] = pion2->gPt();
      ntVar[ii++] = triplet->particle2Dca();
      ntVar[ii++] = pion2->dEdx();
      ntVar[ii++] = pion2->nSigmaPion();
      ntVar[ii++] = pion2->nHitsFit();
      ntVar[ii++] = pion2->nHitsDedx();
      ntVar[ii++] = pion2TOFinvbeta;
      ntVar[ii++] = pion2BetaBase;

      ntVar[ii++] = mPicoHFEvent->runId();
      ntVar[ii++] = mPicoHFEvent->eventId();
      ntVar[ii++] = kaon->gMom(mPrimVtx,mBField).phi();
      ntVar[ii++] = kaon->gMom(mPrimVtx,mBField).pseudoRapidity();
      ntVar[ii++] = kaon->gPt();
      ntVar[ii++] = triplet->particle3Dca();
      ntVar[ii++] = kaon->dEdx();
      ntVar[ii++] = kaon->nSigmaKaon();
      ntVar[ii++] = kaon->nHitsFit();
      ntVar[ii++] = kaon->nHitsDedx();
      ntVar[ii++] = kaonTOFinvbeta;
      ntVar[ii++] = kaonBetaBase;

      ntVar[ii++] = dcaMax;

      ntVar[ii++] = flag;
      ntVar[ii++] = mPrimVtx.z();
      ntVar[ii++] = triplet->pointingAngle();
      ntVar[ii++] = triplet->decayLength();
      ntVar[ii++] = triplet->phi();
      ntVar[ii++] = triplet->eta();
      ntVar[ii++] = sqrt(pow(triplet->px(),2.0)+pow(triplet->py(),2.0));
      ntVar[ii++] = triplet->m();
      ntVar[ii++] = triplet->dV0Max();
      ntVar[ii++] = centrality;
      ntVar[ii++] = mPicoDst->event()->refMult();
      ntVar[ii++] = refmultCor;
      ntVar[ii++] = reweight;
      ntp_DMeson->Fill(ntVar);
*/

    } // for (unsigned int idx = 0; idx <  mPicoHFEvent->nHFSecondaryVertices(); ++idx) {
  }

 return kStOK;
}

// _________________________________________________________
bool StPicoDpmAnaMaker::isHadron(StPicoTrack const * const trk, int pidFlag) const {
  // -- good hadron
  return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, pidFlag));
}

// _________________________________________________________
bool StPicoDpmAnaMaker::isPion(StPicoTrack const * const trk) const {
  // -- good pion
   StThreeVectorF t = trk->pMom();
   if (fabs(t.pseudoRapidity()) > 1.) return false; //pridano fabs 1212
   //if (!mHFCuts->isHybridTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kPion) ) return false; //SL16d
	 if (!mHFCuts->isHybridTOFHadron(trk, mHFCuts->getTofBetaBase(trk, mPicoDst->event()->bField()), StHFCuts::kPion) ) return false; //SL16j, Vanek
   if (!mHFCuts->cutMinDcaToPrimVertex(trk, StPicoCutsBase::kPion)) return false;
  return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, StPicoCutsBase::kPion));
}

// _________________________________________________________
bool StPicoDpmAnaMaker::isKaon(StPicoTrack const * const trk) const {
  // -- good kaon
  StThreeVectorF t = trk->pMom();
  if (fabs(t.pseudoRapidity()) > 1.) return false;//pridano fabs 1212
	//if (!mHFCuts->isHybridTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kKaon) ) return false; //SL16d
  if (!mHFCuts->isHybridTOFHadron(trk, mHFCuts->getTofBetaBase(trk, mPicoDst->event()->bField()), StHFCuts::kKaon) ) return false; //SL16j, Vanek
  if (!mHFCuts->cutMinDcaToPrimVertex(trk, StPicoCutsBase::kKaon)) return false;
  return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, StPicoCutsBase::kKaon));
}

// _________________________________________________________
bool StPicoDpmAnaMaker::isProton(StPicoTrack const * const trk) const {
  // -- good proton
  return (mHFCuts->isGoodTrack(trk) && mHFCuts->isTPCHadron(trk, StPicoCutsBase::kProton));
}

double StPicoDpmAnaMaker::DCA(StPicoTrack const * const trk, StThreeVectorF const & vtx) const {
  // -- particle DCA
/*  StPhysicalHelixD pHelix = trk->dcaGeometry().helix(); //SL16d
  pHelix.moveOrigin(pHelix.pathLength(vtx));
  return ((pHelix.origin() - vtx).mag());
*/
	return ((trk->origin() - vtx).mag()); //SL16j, Vanek
}


bool StPicoDpmAnaMaker::isCloseTracks(StPicoTrack const * const trk1, StPicoTrack const * const trk2, StThreeVectorF const & vtx, float bField) const {
/* SL16d
  StPhysicalHelixD p1Helix = trk1->dcaGeometry().helix();
  StPhysicalHelixD p2Helix = trk2->dcaGeometry().helix();
  p1Helix.moveOrigin(p1Helix.pathLength(vtx));
  p2Helix.moveOrigin(p2Helix.pathLength(vtx));
  if( ( p1Helix.origin()-vtx ).mag()>0.2 || ( p2Helix.origin()-vtx ).mag()>0.2 ) return false;
*/
	if( ( trk1->origin()-vtx ).mag()>0.2 || ( trk2->origin()-vtx ).mag()>0.2 ) return false; //SL16j, Vanek

  //Requires loading constants
//  StThreeVectorF const p1Mom = p1Helix.momentum(bField * kilogauss); //SL16d
//  StThreeVectorF const p2Mom = p2Helix.momentum(bField * kilogauss);

	StThreeVectorF const p1Mom = trk1->gMom(); //SL16j, Vanek
  StThreeVectorF const p2Mom = trk2->gMom();
  StPhysicalHelixD const p1StraightLine(p1Mom, trk1->origin(), 0, trk1->charge());
  StPhysicalHelixD const p2StraightLine(p2Mom, trk2->origin(), 0, trk2->charge());
  //DCA
  pair<double, double> const ss = p1StraightLine.pathLengths(p2StraightLine);
  StThreeVectorF const p1AtDcaToP2 = p1StraightLine.at(ss.first);
  StThreeVectorF const p2AtDcaToP1 = p2StraightLine.at(ss.second);
  float const dca = (p1AtDcaToP2-p2AtDcaToP1).mag();
  if(dca > 0.009) return false;
// -- good pair
  return true;
}

//-----------------------------------------------------------------------------

void StPicoDpmAnaMaker::histoInit(TString fileBaseName, bool fillQaHists){
  TString m_ParticleName[m_nParticles] = {"Pion", "Kaon", "Proton"};

   float m_EtaEdgeDca[m_nEtasDca+1] = {-1.0, -0.6, -0.2, 0.2, 0.6, 1.0}; //replace bottom!!!
   float m_PhiEdgeDca[m_nPhisDca + 1] = {-3.14159, -2.80359, -2.17527, -1.54696, -0.918637, -0.290319, 0.338, 0.966319, 1.59464, 2.22296, 2.85127, 3.14159};
   float m_VzEdgeDca[m_nVzsDca + 1] = { -6.0, -3.0, 0, 3.0, 6.0};//replace bottom!!!
   float m_CentEdgeDca[m_nCentsDca + 1] = { -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5};
   float m_PtEdgeDca[m_nPtsDca + 1] = {0.3, 0.4, 0.5, 0.6,  0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 6.0, 12.0};
   float m_EtaEdgeRatio[m_nEtasRatio + 1] = { -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4 , 0.6, 0.8, 1.0}; //replace bottom!!!
   float m_PhiEdgeRatio[m_nPhisRatio + 1] = { -3.14159, -2.80359, -2.17527, -1.54696, -0.918637, -0.290319, 0.338, 0.966319, 1.59464, 2.22296, 2.85127, 3.14159};//replace bottom!!!
   float m_VzEdgeRatio[m_nVzsRatio + 1] = { -6.0, -4.0, -2.0, 0, 2.0, 4.0, 6.0};//replace bottom!!!
   float m_CentEdgeRatio[m_nCentsRatio + 1] = { -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5};
   float m_PtEdgeRatio[m_nPtsRatio + 1] =
   {
      0.3, 0.4, 0.5, 0.6 , 0.7 , 0.8 , 0.9 ,
      1. , 1.1 , 1.2 , 1.3 , 1.4 , 1.5 , 1.6 , 1.7 , 1.8 , 1.9 ,
      2. , 2.2 , 2.4 , 2.6 , 2.8 , 3.0 ,
      3.4 , 3.8 , 4.2 , 4.6 , 5.0 ,  5.5 ,
      6. , 6.5 , 7.0 , 8.0 , 9.0 , 10. , 11,  12.0
   };
  float m_DcaEdgeDca[m_nDcasDca + 1] =
   {
     -1 , -0.96 , -0.92 , -0.88 , -0.84 , -0.8 , -0.76 , -0.72 , -0.68 , -0.64 , -0.6 , -0.56 , -0.52 , -0.48 , -0.44 , -0.4 , -0.36 , -0.32 , -0.28 , -0.24 , -0.2 , -0.16 , -0.12 ,  -0.08,
     -0.078 , -0.075 , -0.072 , -0.069 , -0.066 , -0.063 , -0.06 , -0.057 , -0.054 , -0.051 , -0.048 , -0.045 , -0.042 , -0.039 , -0.036 , -0.033 , -0.03 , -0.027 , -0.024 , -0.021 , -0.018 , -0.015 , -0.012 ,
      -0.01 , -0.0096 , -0.0092 , -0.0088 , -0.0084 , -0.008 , -0.0076 , -0.0072 , -0.0068 , -0.0064 , -0.006 , -0.0056 , -0.0052 , -0.0048 , -0.0044 , -0.004 , -0.0036 , -0.0032 , -0.0028 , -0.0024 , -0.002 , -0.0016 , -0.0012 , -0.0008 , -0.0004 , 0 , 0.0004 , 0.0008 , 0.0012 , 0.0016 , 0.002 , 0.0024 , 0.0028 , 0.0032 , 0.0036 , 0.004 , 0.0044 , 0.0048 , 0.0052 , 0.0056 , 0.006 , 0.0064 , 0.0068 , 0.0072 , 0.0076 , 0.008 , 0.0084 , 0.0088 , 0.0092 , 0.0096 , 0.01 ,
      0.012 , 0.015 , 0.018 , 0.021 , 0.024 , 0.027 , 0.03 , 0.033 , 0.036 , 0.039 , 0.042 , 0.045 , 0.048 , 0.051 , 0.054 , 0.057 , 0.06 , 0.063 , 0.066 , 0.069 , 0.072 , 0.075 , 0.078 ,
      0.08 , 0.12 , 0.16 , 0.2 , 0.24 , 0.28 , 0.32 , 0.36 , 0.4 , 0.44 , 0.48 , 0.52 , 0.56 , 0.6 , 0.64 , 0.68 , 0.72 , 0.76 , 0.8 , 0.84 , 0.88 , 0.92 , 0.96 , 1
   };



   //set in private
   // for(int temp = 0;temp<m_nParticles;temp++) m_ParticleName[temp]=temp_ParticleName[temp];
	//for(int temp = 0;temp<m_nEtasDca+1;temp++) m_EtaEdgeDca[temp]=temp_EtaEdgeDca[temp];
	//for(int temp2 = 0;temp2<m_nPhisDca+1;temp2++) m_PhiEdgeDca[temp2]=temp_PhiEdgeDca[temp2];
	/*for(int temp = 0;temp<m_nVzsDca+1;temp++) m_VzEdgeDca[temp]=temp_VzEdgeDca[temp];
	for(int temp = 0;temp<m_nCentsDca+1;temp++) m_CentEdgeDca[temp]=temp_CentEdgeDca[temp];
	for(int temp = 0;temp<m_nPtsDca+1;temp++) m_PtEdgeDca[temp]=temp_PtEdgeDca[temp];
	for(int temp = 0;temp<m_nEtasRatio+1;temp++) m_EtaEdgeRatio[temp]=temp_EtaEdgeRatio[temp];
	for(int temp = 0;temp<m_nPhisRatio+1;temp++) m_PhiEdgeRatio[temp]=temp_PhiEdgeRatio[temp];
	for(int temp = 0;temp<m_nVzsRatio+1;temp++) m_VzEdgeRatio[temp]=temp_VzEdgeRatio[temp];
	for(int temp = 0;temp<m_nCentsRatio+1;temp++) m_CentEdgeRatio[temp]=temp_CentEdgeRatio[temp];
	for(int temp = 0;temp<m_nPtsRatio+1;temp++) m_PtEdgeRatio[temp]=temp_PtEdgeRatio[temp];
	for(int temp = 0;temp<m_nDcasDca+1;temp++) m_DcaEdgeDca[temp]=temp_DcaEdgeDca[temp];*/

   mFillQaHists = fillQaHists;
   //mOutFile = new TFile(fileBaseName+".hists.root", "RECREATE"); //comment Vanek - want only one output file for QA
   //mOutFile = new TFile(Form("%s.hists.root", fileBaseName.Data()), "RECREATE");
/*
   for (int iParticle = 0; iParticle < m_nParticles; iParticle++){
      for (int iEta = 0; iEta < m_nEtasDca; iEta++){
         for (int iVz = 0; iVz < m_nVzsDca; iVz++){
            for (int iCent = 0; iCent < m_nCentsDca; iCent++){
               mh3DcaXyZPtCentPartEtaVzPhi[iParticle][iEta][iVz][iCent] = NULL;
            }
         }
      }
   }

   for (int iParticle = 0; iParticle < m_nParticles; iParticle++){
      for (int iEta = 0; iEta < m_nEtasRatio; iEta++){
         for (int iVz = 0; iVz < m_nVzsRatio; iVz++){
            for (int iPhi = 0; iPhi < m_nPhisRatio; iPhi++){
               mh2Tpc1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi] = NULL;
               mh2HFT1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi] = NULL;
            }
         }
      }
   }*/



   TH1::SetDefaultSumw2();
   if (!mFillQaHists) return;

   mh1Cent         = new TH1F("mh1Cent", "EventsVsCentrality;cent;Counts", 10, -1.5, 8.5);
   mh1CentWg         = new TH1F("mh1CentWg", "EventsVsCentrality;cent;Counts", 10, -1.5, 8.5);
   mh1gRefmultCor  = new TH1F("mh1gRefmultCor", "gRefmultCor;gRefmult;Counts", 700, 0, 700);
   mh1gRefmultCorWg  = new TH1F("mh1gRefmultCorWg", "gRefmultCorWg;gRefmultCorWg;Counts", 700, 0, 700);
   mh2CentVz         = new TH2F("mh2CentVz", "CentralityVsVz;cent;Vz", 10, -1.5, 8.5, 200, -10, 10);
   mh2CentVzWg = new TH2F("mh2CentVzWg", "CentralityVsVzWg;cent;Vz", 10, -1.5, 8.5, 200, -10, 10);

   //Add some HFT ratio plots
   mh2Tpc1PtCent  = new TH2F("mh2Tpc1PtCent", "Tpc tacks;p_{T}(GeV/c);cent", 120, 0, 12, 10, -1.5, 8.5); //Dca 1.5cm
   mh2HFT1PtCent  = new TH2F("mh2HFT1PtCent", "HFT tacks;p_{T}(GeV/c);cent", 120, 0, 12, 10, -1.5, 8.5); //Dca 1.5cm
   mh2Tpc1PhiVz  = new TH2F("mh2Tpc1PhiVz", "Tpc tacks;#Phi;Vz", 100, -3.1415, 3.1415, 20, -10, 10); //Dca 1.5cm
   mh2HFT1PhiVz  = new TH2F("mh2HFT1PhiVz", "HFT tacks;#Phi;Vz", 100, -3.1415, 3.1415, 20, -10, 10); //Dca 1.5cm

   for (int iParticle = 0; iParticle < m_nParticles; iParticle++){
      for (int iEta = 0; iEta < m_nEtasRatio; iEta++){
         for (int iVz = 0; iVz < m_nVzsRatio; iVz++){
            for (int iPhi = 0; iPhi < m_nPhisRatio; iPhi++){

               mh2Tpc1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi]  = new TH2F(Form("mh2Tpc1PtCentPartEtaVzPhi_%d_%d_%d_%d", iParticle, iEta, iVz, iPhi), "mh2Tpc1PtCent_"+m_ParticleName[iParticle]+Form("_Eta%2.1f_Vz%2.1f_Phi%2.1f;p_{T}(GeV/c);cent", m_EtaEdgeRatio[iEta], m_VzEdgeRatio[iVz], m_PhiEdgeRatio[iPhi]), m_nPtsRatio, m_PtEdgeRatio, m_nCentsRatio, m_CentEdgeRatio); //Dca 1.cm
               mh2HFT1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi]  = new TH2F(Form("mh2HFT1PtCentPartEtaVzPhi_%d_%d_%d_%d", iParticle, iEta, iVz, iPhi), "mh2HFT1PtCent_"+m_ParticleName[iParticle]+Form("_Eta%2.1f_Vz%2.1f_Phi%2.1f;p_{T}(GeV/c);cent", m_EtaEdgeRatio[iEta], m_VzEdgeRatio[iVz], m_PhiEdgeRatio[iPhi]), m_nPtsRatio, m_PtEdgeRatio, m_nCentsRatio, m_CentEdgeRatio); //Dca 1.cm
            }
         }
      }
   }

   // Add some Dca, resolution
   for (int iParticle = 0; iParticle < m_nParticles; iParticle++){
      for (int iEta = 0; iEta < m_nEtasDca; iEta++){
         for (int iVz = 0; iVz < m_nVzsDca; iVz++){
            for (int iCent = 0; iCent < m_nCentsDca; iCent++){

   	    	   mh3DcaXyZPtCentPartEtaVzPhi[iParticle][iEta][iVz][iCent]  = new TH3F(Form("mh3DcaXyZPtCentPartEtaVzPhi_%d_%d_%d_%d", iParticle, iEta, iVz, iCent),"mh3DcaXyZPt_"+m_ParticleName[iParticle]+Form("_Eta%2.1f_Vz%2.1f_Cent%2.1f;p_{T}(GeV/c);DcaXy(cm);DcaZ(cm)", m_EtaEdgeDca[iEta], m_VzEdgeDca[iVz], m_CentEdgeDca[iCent]), m_nPtsDca, m_PtEdgeDca, m_nDcasDca, m_DcaEdgeDca, m_nDcasDca, m_DcaEdgeDca); //Dca 1.cm
            }
         }
      }
   }

   mh3DcaPtCent  = new TH3F("mh3DcaPtCent", "mh3DcaPtCent;p_{T}(GeV/c);cent;Dca(cm)", 120, 0, 12, 10, -1.5, 8.5, 1000, -1, 1); //Dca 1.cm
   mh3DcaXyPtCent  = new TH3F("mh3DcaXyPtCent", "mh3DcaXyPtCent;p_{T}(GeV/c);cent;DcaXy(cm)", 120, 0, 12, 10, -1.5, 8.5, 1000, -1, 1); //Dca 1.cm
   mh3DcaZPtCent  = new TH3F("mh3DcaZPtCent", "mh3DcaZPtCent;p_{T}(GeV/c);cent;DcaZ(cm)", 120, 0, 12, 10, -1.5, 8.5, 1000, -1, 1); //Dca 1.cm

}

//-----------------------------------------------------------------------
void StPicoDpmAnaMaker::addTpcDenom1(bool IsPion, bool IsKaon, bool IsProton, float pt, int centrality, float Eta, float Phi, float Vz){
   int EtaIndex = getEtaIndexRatio(Eta);
   int PhiIndex = getPhiIndexRatio(Phi);
   int VzIndex = getVzIndexRatio(Vz);
   if(EtaIndex == -1) return;
   if(PhiIndex == -1) return;
   if(VzIndex == -1) return;
   //std::cout<<"2: "<<IsPion<<" "<<IsKaon<<" "<<IsProton<<" "<<pt<<" "<<centrality<<" "<<Eta<<" "<<Phi<<" "<<Vz<<" "<<EtaIndex<<" "<<PhiIndex<<" "<<VzIndex<<std::endl;

   if (IsPion){
      mh2Tpc1PtCentPartEtaVzPhi[0][EtaIndex][VzIndex][PhiIndex]->Fill(pt, centrality);
      //if(mh2Tpc1PtCentPartEtaVzPhi[0][EtaIndex][VzIndex][PhiIndex]) std::cout<<"true"<<<<std::endl;
      //std::cout<<pt<<" "<<centrality<<std::endl;
   }
   if (IsKaon){
      mh2Tpc1PtCentPartEtaVzPhi[1][EtaIndex][VzIndex][PhiIndex]->Fill(pt, centrality);
   }
   if (IsProton){
      mh2Tpc1PtCentPartEtaVzPhi[2][EtaIndex][VzIndex][PhiIndex]->Fill(pt, centrality);
   }
   mh2Tpc1PtCent->Fill(pt, centrality);
   if (fabs(Eta) < mHFCuts->cutEta()  && pt > mHFCuts->cutPt()) mh2Tpc1PhiVz->Fill(Phi, Vz);
}
//-----------------------------------------------------------------------
void StPicoDpmAnaMaker::addHFTNumer1(bool IsPion, bool IsKaon, bool IsProton, float pt, int centrality, float Eta, float Phi, float Vz){
   int EtaIndex = getEtaIndexRatio(Eta);
   int PhiIndex = getPhiIndexRatio(Phi);
   int VzIndex = getVzIndexRatio(Vz);
   if(EtaIndex == -1) return;
   if(PhiIndex == -1) return;
   if(VzIndex == -1) return;
   if (IsPion){
      mh2HFT1PtCentPartEtaVzPhi[0][EtaIndex][VzIndex][PhiIndex]->Fill(pt, centrality);
   }
   if (IsKaon){
      mh2HFT1PtCentPartEtaVzPhi[1][EtaIndex][VzIndex][PhiIndex]->Fill(pt, centrality);
   }
   if (IsProton){
      mh2HFT1PtCentPartEtaVzPhi[2][EtaIndex][VzIndex][PhiIndex]->Fill(pt, centrality);
   }
   mh2HFT1PtCent->Fill(pt, centrality);
   if (fabs(Eta) < mHFCuts->cutEta()  && pt > mHFCuts->cutPt()) mh2HFT1PhiVz->Fill(Phi, Vz);
}
//---------------------------------------------------------------------
void StPicoDpmAnaMaker::addDcaPtCent(float dca, float dcaXy, float dcaZ, bool IsPion, bool IsKaon, bool IsProton, float pt,  int centrality, float Eta, float Phi, float Vz){
   int EtaIndex = getEtaIndexDca(Eta);
   int VzIndex = getVzIndexDca(Vz);
   if(EtaIndex == -1) return;
   if(VzIndex == -1) return;

   if (centrality < 0) return; // remove bad centrality, only keep 9 centralities
   if (IsPion){
      mh3DcaXyZPtCentPartEtaVzPhi[0][EtaIndex][VzIndex][centrality]->Fill(pt, dcaXy, dcaZ);
   }
   if (IsKaon){
      mh3DcaXyZPtCentPartEtaVzPhi[1][EtaIndex][VzIndex][centrality]->Fill(pt, dcaXy, dcaZ);
   }
   if (IsProton){
      mh3DcaXyZPtCentPartEtaVzPhi[2][EtaIndex][VzIndex][centrality]->Fill(pt, dcaXy, dcaZ);
   }
   mh3DcaPtCent->Fill(pt, centrality, dca);
   mh3DcaXyPtCent->Fill(pt, centrality, dcaXy);
   mh3DcaZPtCent->Fill(pt, centrality, dcaZ);
}
//---------------------------------------------------------------------
int StPicoDpmAnaMaker::getEtaIndexDca(float Eta){
   float EtaEdgeDca[m_nEtasDca+1] = {-1.0, -0.6, -0.2, 0.2, 0.6, 1.0};
   for (int i = 0; i < m_nEtasDca; i++){
	 if ((Eta >= EtaEdgeDca[i]) && (Eta < EtaEdgeDca[i + 1]))
         return i;
   }
   //std::cout<<"SOMETHING WENT TERRIBRU WONG"<<std::endl;
   //return m_nEtasDca -1;
   return -1;
}

//---------------------------------------------------------------------
int StPicoDpmAnaMaker::getVzIndexDca(float Vz){
  float VzEdgeDca[m_nVzsDca + 1] = { -6.0, -3.0, 0, 3.0, 6.0};
   for (int i = 0; i < m_nVzsDca; i++){
      if ((Vz >= VzEdgeDca[i]) && (Vz < VzEdgeDca[i + 1]))
         return i;
   }
//std::cout<<"SOMETHING WENT TERRIBRU WONG"<<std::endl;
   //return m_nVzsDca - 1;
   return -1;
}
//---------------------------------------------------------------------
int StPicoDpmAnaMaker::getEtaIndexRatio(float Eta){
  float EtaEdgeRatio[m_nEtasRatio + 1] = { -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4 , 0.6, 0.8, 1.0};
   for (int i = 0; i < m_nEtasRatio; i++){
      if ((Eta >= EtaEdgeRatio[i]) && (Eta < EtaEdgeRatio[i + 1]))
         return i;
   }
//std::cout<<"SOMETHING WENT TERRIBRU WONG"<<std::endl;
   //return m_nEtasRatio - 1;
   return -1;
}
//---------------------------------------------------------------------
int StPicoDpmAnaMaker::getPhiIndexRatio(float Phi){
  float PhiEdgeRatio[m_nPhisRatio + 1] = { -3.14159, -2.80359, -2.17527, -1.54696, -0.918637, -0.290319, 0.338, 0.966319, 1.59464, 2.22296, 2.85127, 3.14159};
   for (int i = 0; i < m_nPhisRatio; i++){
      if ((Phi >= PhiEdgeRatio[i]) && (Phi < PhiEdgeRatio[i + 1]))
         return i;
   }
//std::cout<<"SOMETHING WENT TERRIBRU WONG"<<std::endl;
  // return m_nPhisRatio - 1;
   return -1;
}
//---------------------------------------------------------------------
int StPicoDpmAnaMaker::getVzIndexRatio(float Vz){
  float VzEdgeRatio[m_nVzsRatio + 1] = { -6.0, -4.0, -2.0, 0, 2.0, 4.0, 6.0};
   for (int i = 0; i < m_nVzsRatio; i++) {
      if ((Vz >= VzEdgeRatio[i]) && (Vz < VzEdgeRatio[i + 1]))
         return i;
   }
//std::cout<<"SOMETHING WENT TERRIBLE WRONG"<<std::endl;
  // return m_nVzsRatio - 1;
  return -1;
}

void StPicoDpmAnaMaker::addCent(const double refmultCor, int centrality, const double reweight, const float vz)
{
   mh1gRefmultCor->Fill(refmultCor);
   mh1gRefmultCorWg->Fill(refmultCor, reweight);
   mh1Cent->Fill(centrality);
   mh1CentWg->Fill(centrality, reweight);
   mh2CentVz->Fill(centrality, vz);
   mh2CentVzWg->Fill(centrality, vz, reweight);
}

//---------------------------------------------------------------------
void StPicoDpmAnaMaker::closeFile()
{
   //mOutFile->cd(); //do not want second file this time
//	cout<<"CloseFile1"<<endl;
   mh1Cent->Write();
   mh1CentWg->Write();
   mh1gRefmultCor->Write();
   mh1gRefmultCorWg->Write();
   mh2CentVz->Write();
   mh2CentVzWg->Write();
//	cout<<"CloseFile2"<<endl;
   //HFT ratio QA
   mh2Tpc1PtCent->Write();
   mh2Tpc1PhiVz->Write();
   mh2HFT1PhiVz->Write();
   mh2HFT1PtCent->Write();
//	cout<<"CloseFile3"<<endl;
   //HFT DCA Ratio
   for (int iParticle = 0; iParticle < m_nParticles; iParticle++)
   {
      for (int iEta = 0; iEta < m_nEtasDca; iEta++)
      {
         for (int iVz = 0; iVz < m_nVzsDca; iVz++)
         {
            for (int iCent = 0; iCent < m_nCentsDca; iCent++)
            {
               mh3DcaXyZPtCentPartEtaVzPhi[iParticle][iEta][iVz][iCent]->Write();
            }
         }
      }
   }
  // std::cout<<"tuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu"<<m_nParticles<<" "<<m_nEtasRatio<<std::endl;
//cout<<"CloseFile4"<<endl;
   for (int iParticle = 0; iParticle < m_nParticles; iParticle++)
   {
      for (int iEta = 0; iEta < m_nEtasRatio; iEta++)
      {
         for (int iVz = 0; iVz < m_nVzsRatio; iVz++)
         {
            for (int iPhi = 0; iPhi < m_nPhisRatio; iPhi++)
            {
               mh2Tpc1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi]->Write();
               mh2HFT1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi]->Write();
            }
         }
      }
   }
//cout<<"CloseFile5"<<endl;
   mh3DcaPtCent->Write();
   mh3DcaXyPtCent->Write();
   mh3DcaZPtCent->Write();
//cout<<"CloseFile6"<<endl;
   // nt->Write();
   //mOutFile->Write();
   //mOutFile->Close();
   //mOutFile->Delete();
}


