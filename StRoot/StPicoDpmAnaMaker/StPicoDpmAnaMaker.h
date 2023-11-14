#ifndef StPicoDpmAnaMaker_h
#define StPicoDpmAnaMaker_h

#include "StPicoHFMaker/StPicoHFMaker.h"
#include "TNtuple.h"
#include "StRefMultCorr/StRefMultCorr.h"
#include "TH2F.h"

#include <vector>

#include "TClonesArray.h"

#include "StThreeVectorF.hh"
#include "StLorentzVectorF.hh"

#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoDstMaker.h"

#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"

#include "StPicoHFMaker/StPicoHFEvent.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "StPicoHFMaker/StHFPair.h"
#include "StPicoHFMaker/StHFTriplet.h"
#include "StBTofUtil/tofPathLength.hh"

#include "phys_constants.h"

#include "TH1F.h"
#include "TH3F.h"


#include <ctime>

/* **************************************************
 *  Sample class fo HF picoDST analysis
 * --------------------------------------------------
 * 
 *  For more info look also in the .h files in StPicoHFMaker/
 *     StPicoHFMaker/StPicoHFMaker.h      <-- Base Class for analysis
 *     StPicoHFMaker/StPicoHFEvent.h      <-- Holds candidates for one event (written to Tree)
 *     StPicoHFMaker/StHFCuts.h           <-- Cuts, can be set in run macro
 *     StPicoHFMaker/StHFPair.h           <-- Holds a pair candidate of a two body decay
 *     StPicoHFMaker/StHFTriplet.h        <-- Holds a triplet of a three body decay
 *
 *  Usage:
 *   - Implement
 *        InitHF()
 *        MakeHF()
 *        ClearHF()
 *        FinishHF()
 *
 *  - Do not ovewrite Init, Make, Clear, Finish which are inhertited from StPicoHFMaker via StMaker 

 *  - Set StHFCuts class via setHFBaseCuts(...) in run macro
 *
 *  - Set use mode of StPicoHFMaker class  via setMakerMode(...)
 *     use enum of StPicoHFMaker::eMakerMode
 *      StPicoHFMaker::kAnalyze - don't write candidate trees, just fill histograms
 *      StPicoHFMaker::kWrite   - write candidate trees
 *      StPicoHFMaker::kRead    - read candidate trees and fill histograms
 *
 *  - Set decay mode of analysis via setDecayMode(...)
 *     use enum of StPicoHFEvent::eHFEventMode (see there for more info)
 *      StPicoHFEvent::kTwoParticleDecay,
 *      StPicoHFEvent::kThreeParticleDecay
 *      StPicoHFEvent::kTwoAndTwoParticleDecay
 *
 *  - Implement these track selection methods used to fill vectors for 'good' identified particles
 *      (methods from StHFCuts utility class can/should be used)
 *       isPion
 *       isKaon
 *       isProton
 *
 *  --------------------------------------------------
 *  
 *  Initial Authors:  
 *            Xin Dong        (xdong@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *          **Jochen Thaeder  (jmthader@lbl.gov) 
 * 
 *  ** Code Maintainer
 *
 * **************************************************
 */

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StPicoHFEvent;

class StHFPair;
class StHFTriplet;
class StHFCuts;

class StRefMultCorr;

class StPicoDpmAnaMaker : public StPicoHFMaker 
{
 public:
  StPicoDpmAnaMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName,  
		       char const* inputHFListHFtree);
  virtual ~StPicoDpmAnaMaker();
  
  virtual Int_t InitHF();
  virtual Int_t MakeHF();
  virtual void  ClearHF(Option_t *opt);
  virtual Int_t FinishHF();

  virtual bool isCloseTracks(StPicoTrack const*, StPicoTrack const*,StThreeVectorF const & , float) const;
  virtual double DCA(StPicoTrack const*, StThreeVectorF const &) const;
  int createQA();
  
  // -- ADOPT DECAY CHANNELS, if wished ------------------- 
  void setDecayChannel(unsigned int u) { mDecayChannel = u; }

  enum eDecayChannel {kChannel1, kChannel2, kChannel3};

  void setRefMutCorr(StRefMultCorr* gRefMultCorr) { mRefmultCorrUtil = gRefMultCorr; }
  StRefMultCorr* getRefMultCorr() { return mRefmultCorrUtil; }

   void histoInit(TString fileBaseName,bool fillQaHists=true);
   void addTpcDenom1(bool IsPion, bool IsKaon, bool IsProton, float pt, int centrality, float reweight, float Eta, float Phi, float Vz);
   void addHFTNumer1(bool IsPion, bool IsKaon, bool IsProton, float pt, int centrality, float reweight, float Eta, float Phi, float Vz);
   void addDcaPtCent(float dca, float dcaXy, float  dcaZ, bool IsPion, bool IsKaon, bool IsProton, float pt,  int centrality, float reweight, float Eta, float Phi, float Vz);
   int getEtaIndexDca(float Eta) ;
   int getPhiIndexDca(float Phi) ;
   int getVzIndexDca(float Vz) ;
   int getEtaIndexRatio(float Eta) ;
   int getPhiIndexRatio(float Phi) ;
   int getVzIndexRatio(float Vz) ;
   void closeFile();


 protected:
  virtual bool isHadron(StPicoTrack const*, int pidFlag) const;
  virtual bool isPion(StPicoTrack const*) const;
  virtual bool isKaon(StPicoTrack const*) const;
  virtual bool isProton(StPicoTrack const*) const;

private:
  int createCandidates();
  int analyzeCandidates();
  int FillTrackEffErrData();



  // -- private members --------------------------

  unsigned int mDecayChannel;


  // -- ADD USER MEMBERS HERE ------------------- 
   //TNtuple *ntp_DMeson; //orig. Kvapil
   TTree *Pi_PID_eff; //Vanek
   TTree *K_PID_eff;

   TTree *TrackEffErr;

   StRefMultCorr* mRefmultCorrUtil;
   int mRunNumber;
       
  TString mOutFileBaseName;

  bool mFillQaHists;
   TFile* mOutFile;

   //Cuts----------------------------
   static const int m_nParticles = 2; //removed proton, 11/08/18, Vanek
   //TString m_ParticleName[m_nParticles];

   static const int m_nEtasDca = 5;
  // float m_EtaEdgeDca[m_nEtasDca+1];
   static const int m_nPhisDca = 11;
   //static float m_PhiEdgeDca[m_nPhisDca + 1];

   static const int m_nVzsDca = 4;
   //static float m_VzEdgeDca[m_nVzsDca + 1];

   static const int m_nCentsDca = 9;
   //static float m_CentEdgeDca[m_nCentsDca + 1];

   static const int m_nPtsDca = 19;
  // static float m_PtEdgeDca[m_nPtsDca + 1];

   static const int m_nEtasRatio = 10; //original 10 - alternative version 5
  // static float m_EtaEdgeRatio[m_nEtasRatio + 1];

   static const int m_nPhisRatio = 11;
   //static float m_PhiEdgeRatio[m_nPhisRatio + 1];

   static const int m_nVzsRatio = 6;
   //static float m_VzEdgeRatio[m_nVzsRatio + 1];

   static const int m_nCentsRatio = 10;
  // static float m_CentEdgeRatio[m_nCentsRatio + 1];

   static const int m_nPtsRatio = 36;
  // static float m_PtEdgeRatio[m_nPtsRatio + 1];

  static const int m_nDcasDca = 144;
 // static float m_DcaEdgeDca[m_nDcasDca + 1];
   //-----------------------------------

   TH1F* mh1Cent;
   TH1F* mh1CentWg;
   TH1F* mh1gRefmultCor;
   TH1F* mh1gRefmultCorWg;
   TH2F* mh2CentVz;
   TH2F* mh2CentVzWg;


   //HFT ratio QA
   TH2F* mh2Tpc1PtCent;
   TH2F* mh2Tpc1PhiVz;
   TH2F* mh2HFT1PtCent;
   TH2F* mh2HFT1PhiVz;
   TH2D* mh2Tpc1PtCentPartEtaVzPhi[m_nParticles][m_nEtasRatio][m_nVzsRatio][m_nPhisRatio];
   TH2D* mh2HFT1PtCentPartEtaVzPhi[m_nParticles][m_nEtasRatio][m_nVzsRatio][m_nPhisRatio];

   //HFT Dca
   TH3F* mh3DcaXyZPtCentPartEtaVzPhi[m_nParticles][m_nEtasDca][m_nVzsDca][m_nCentsDca]; //changed back to TH3F to test memory usage

   TH3F* mh3DcaPtCent;
   TH3F* mh3DcaXyPtCent;
   TH3F* mh3DcaZPtCent;

//---Variables for PID TTree---------------------------
	Float_t Pi1nSigmaTPC, K1nSigmaTPC;
  Float_t Pi1nSigmaTOF, K1nSigmaTOF; // nSigma TOF = |1/beta - 1/beta_TOF|

  Float_t Pi1Pt, K1Pt; //single particle pT
  Float_t Pi1P, K1P; //single particle p
  
  Float_t Pi2nSigmaTPC, K2nSigmaTPC;
  Float_t Pi2nSigmaTOF, K2nSigmaTOF; // nSigma TOF = |1/beta - 1/beta_TOF|

  Float_t Pi2Pt, K2Pt; //single particle pT
  Float_t Pi2P, K2P; //single particle p

  Float_t PiPairPt, KPairPt; //pT of mothers of pi/K pairs (K0s, phi)

  Float_t PiPairInvMass, KPairInvMass;

  Int_t PiPairCharge, KPairCharge; //0 = unlike-sign, 1 = like-sign

//---Variables for tracking efficiency sys. error
  Float_t Track_pt;
  Float_t Track_nHitsFit;
  Float_t Track_nHitsMax;
  Float_t Track_nSigmaTPCPi, Track_nSigmaTPCK;
  Float_t Track_nSigmaTOFPi, Track_nSigmaTOFK;
  Int_t Track_isHFT; //0 = no HFT (track not matched to the HFT), 1 = isHFTTrack

  Int_t Track_centrality;
  Float_t Track_reweight;


	

//-------------------------------------------------
  // -- ADD USER MEMBERS HERE -------------------

  ClassDef(StPicoDpmAnaMaker, 1) //set to 1
};

#endif
