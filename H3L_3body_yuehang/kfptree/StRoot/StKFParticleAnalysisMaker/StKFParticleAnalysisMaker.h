// $Id: StKFParticleAnalysisMaker.h,v 1.1 2021/04/19 22:04:11 slan Exp $
/*!
 * \class  StKFParticleAnalysisMaker
 * \author Maksym Zyzak
 * \date   2017/10/17
 * \brief  class for analysis of PicoDst
 */
#ifndef STAR_StKFParticleAnalysisMaker
#define STAR_StKFParticleAnalysisMaker
//#define __DEVT__
#ifndef StMaker_H
#include "StMaker.h"
#endif
#include "TMVA/Reader.h"
#include "EventClass/EventClass.h"

class StKFParticleInterface;
class StKFParticlePerformanceInterface;
class KFParticle;
class StPicoDst;
class StMuDst;
class TNtuple;
class TFile;
class TH2F;
class TChain;
class StRefMultCorr;
class CentralityMaker;

class StKFParticleAnalysisMaker : public StMaker {
private:
  // static const int fNNTuples = 8;
  //static const int fNNTuples = 9;
  //static const int fNNTuples = 1;
  //static const int fNNTuples = 2;
  static const int fNNTuples = 3;
  Char_t                mBeg[1];        //!
  StMuDst                          *fMuDst;
  StPicoDst                        *fPicoDst;                          //!
  StKFParticleInterface            *fStKFParticleInterface;            //!
  StKFParticlePerformanceInterface *fStKFParticlePerformanceInterface; //!
  TNtuple* fCutsNTuple[fNNTuples];
  TFile* fNTupleFile[fNNTuples];
  TFile* f;
  int fNTuplePDG[fNNTuples];
  TString fNtupleNames[fNNTuples];
  TString fNtupleCutNames[fNNTuples];
  std::vector<TString> fDaughterNames[fNNTuples];
  vector< vector<TString> > fTMVACutFile[fNNTuples];
  vector< vector<double> > fTMVACut[fNNTuples];
  vector< vector<TMVA::Reader*> > fTMVAReader[fNNTuples];
  std::vector<int> fTMVACentralityBins[fNNTuples];
  std::vector<double> fTMVAPtBins[fNNTuples];
  Char_t                mEnd[1];        //!
  std::vector<float> fTMVAParticleParameters[fNNTuples];
  int fNTrackTMVACuts;
  bool fIsPicoAnalysis;
  int fsnn;
  int fdEdXMode;
  Bool_t fStoreTmvaNTuples;
  Bool_t fProcessSignal;
  Bool_t fCollectTrackHistograms;
  Bool_t fCollectPIDHistograms;
  Bool_t fTMVAselection;
  Bool_t fStoremctree;
  
  std::vector<Int_t> mBadRuns;
  void SetBadRun();
  Bool_t IsBadRun(Int_t RunId);
  
  std::vector<int> ftrkIdV0A;
  std::vector<int> ftrkIdV0B;
  
  //tree
  TFile *mMyFile;
  EventClass *pEve;
  TTree *ptree;
  
  int beventid;
  int brunid;
  
  float bVx;
  float bVy;
  float bVr;
  float bVz;
  float bVrerr;
  float bVxerr;
  float bVyerr;
  float bVzerr;
  
  int bismc;
  int cent9;
  double refmultcor;
  double reweight;
  
  void SetHistograms();
  
  TH1F *hvtx;
  TH1F *hvtxgood;
  TH1F *hrefmult;
  TH1F *wrefmult;
  TH1F *hrefmult_wPileup;
  TH1F *hrefmult_woPileup;
  TH1F *hCentrality;
  TH1F *hCentrality_weight;
  TH2F *hvtx_xy;
  TH2F *hvtx_xy_good;
  
  static const int nHVar = 17;
  TH1F *histVar[7][nHVar];
  TH2F *histRap[2];
  
  // for check
  TH2F *hNBTofvsRM  ;
  TH2F *hNBTofvsGRM ;
  TH2F *hTofMvsRM   ;
  TH2F *hTofMvsGRM  ;
  TH2F *hRankvsGTrk ;
  TH2F *hRankvsPTrk ;
  TH2F *hRankvsNBTof;
  TH2F *hGTrkvsPTrk ;
  TH2F *hGTrkvsRM   ;
  TH2F *hGRMvsRM    ;
  TH2F *hNVpdvsRM   ;
  TH2F *hNVpdvsGRM  ;
  TH2F *hNVpdvsGTrk ;
  TH2F *hNBEMCvsRM  ;
  TH2F *hNBEMCvsGRM ;
  TH2F *hNBEMCvsGTrk;
  
  bool isPileUp;
  int refMultPrim;
  float gweight;
  int countrefmult;
  double psi_1_EPD_0;
  double psi_1_EPD_1;
  double psi_1_EPD_2;
  double psi_1_EPD_3;
  double psi_1_EPD_4;
  double psi_1_EPD_5;
  double psi_1_EPD_6;
  double psi_1_EPD;
  
  //Centrality and flow
  Bool_t fFlowAnalysis;
  TChain* fFlowChain;
  int fFlowRunId;
  int fFlowEventId;
  int fCentrality;
  std::vector<TString> fFlowFiles;
  std::map<long, int> fFlowMap;
  
  bool fRunCentralityAnalysis;
  StRefMultCorr *fRefmultCorrUtil;
  TString fCentralityFile;
  
  bool fAnalyseDsPhiPi;
  
  void GetDaughterParameters(const int iReader, int& iDaughterTrack, int& iDaughterParticle, KFParticle& particle);
  void GetParticleParameters(const int iReader, KFParticle& particle);
  long  GetUniqueEventId(const int iRun, const int iEvent) const;
  
  int GetTMVACentralityBin(int iReader, int centrality);
  int GetTMVAPtBin(int iReader, double pt);
  void SetTMVACentralityBins(int iReader, TString bins);
  void SetTMVAPtBins(int iReader, TString bins);
  void SetTMVABins(int iReader, TString centralityBins="-1:1000", TString ptBins="-1.:1000.");
  
public:
  StKFParticleAnalysisMaker(const char *name="KFParticleAnalysis");
  virtual       ~StKFParticleAnalysisMaker();
  virtual Int_t  Init();
  virtual Int_t  InitRun(Int_t runumber);
  void           BookVertexPlots();
  virtual Int_t  Make();
  virtual Int_t  Finish();
  Bool_t         Check();
  void AnalysePicoDst() { fIsPicoAnalysis = true;  }
  void AnalyseMuDst()   { fIsPicoAnalysis = false; }
  void SetDataSet(int snn) { fsnn = snn; }
  static void    PrintMem(const Char_t *opt = "");
  virtual const char *GetCVS() const {
    static const char cvs[]="Tag $Name:  $ $Id: StKFParticleAnalysisMaker.h,v 1.1 2021/04/19 22:04:11 slan Exp $ built " __DATE__ " " __TIME__ ;
    return cvs;
  }
  void Storemc() {fStoremctree = true; }
  void ProcessSignal() { fProcessSignal = true; }
  void DoNotProcessSignal() { fProcessSignal = true; }
  void StoreTMVANtuples() { fStoreTmvaNTuples = true; }
  void CollectTrackHistograms() { fCollectTrackHistograms = true; }
  void CollectPIDHistograms() { fCollectPIDHistograms = true; }
  void UseTMVA() { fTMVAselection = true; }
  void SetTMVABinsD0   (TString centralityBins, TString ptBins) { SetTMVABins(0, centralityBins, ptBins); }
  void SetTMVABinsDPlus(TString centralityBins, TString ptBins) { SetTMVABins(1, centralityBins, ptBins); }
  void SetTMVABinsDs   (TString centralityBins, TString ptBins) { SetTMVABins(2, centralityBins, ptBins); }
  void SetTMVABinsLc   (TString centralityBins, TString ptBins) { SetTMVABins(3, centralityBins, ptBins); }
  void SetTMVABinsD0KK (TString centralityBins, TString ptBins) { SetTMVABins(4, centralityBins, ptBins); }
  void SetTMVABinsD04  (TString centralityBins, TString ptBins) { SetTMVABins(5, centralityBins, ptBins); }
  void SetTMVABinsBPlus(TString centralityBins, TString ptBins) { SetTMVABins(6, centralityBins, ptBins); }
  void SetTMVABinsB0   (TString centralityBins, TString ptBins) { SetTMVABins(7, centralityBins, ptBins); }
  void SetTMVAcutsD0   (TString file, double cut, int iCentralityBin = 0, int iPtBin = 0) { fTMVACutFile[0][iCentralityBin][iPtBin] = file; fTMVACut[0][iCentralityBin][iPtBin] = cut; }
  void SetTMVAcutsDPlus(TString file, double cut, int iCentralityBin = 0, int iPtBin = 0) { fTMVACutFile[1][iCentralityBin][iPtBin] = file; fTMVACut[1][iCentralityBin][iPtBin] = cut; }
  void SetTMVAcutsDs   (TString file, double cut, int iCentralityBin = 0, int iPtBin = 0) { fTMVACutFile[2][iCentralityBin][iPtBin] = file; fTMVACut[2][iCentralityBin][iPtBin] = cut; }
  void SetTMVAcutsLc   (TString file, double cut, int iCentralityBin = 0, int iPtBin = 0) { fTMVACutFile[3][iCentralityBin][iPtBin] = file; fTMVACut[3][iCentralityBin][iPtBin] = cut; }
  void SetTMVAcutsD0KK (TString file, double cut, int iCentralityBin = 0, int iPtBin = 0) { fTMVACutFile[4][iCentralityBin][iPtBin] = file; fTMVACut[4][iCentralityBin][iPtBin] = cut; }
  void SetTMVAcutsD04  (TString file, double cut, int iCentralityBin = 0, int iPtBin = 0) { fTMVACutFile[5][iCentralityBin][iPtBin] = file; fTMVACut[5][iCentralityBin][iPtBin] = cut; }
  void SetTMVAcutsBPlus(TString file, double cut, int iCentralityBin = 0, int iPtBin = 0) { fTMVACutFile[6][iCentralityBin][iPtBin] = file; fTMVACut[6][iCentralityBin][iPtBin] = cut; }
  void SetTMVAcutsB0   (TString file, double cut, int iCentralityBin = 0, int iPtBin = 0) { fTMVACutFile[7][iCentralityBin][iPtBin] = file; fTMVACut[7][iCentralityBin][iPtBin] = cut; }
  
  void RunFlowAnalysis()         { fFlowAnalysis = true; }
  void AddFlowFile(TString file) { fFlowFiles.push_back(file); }
  
  void RunCentralityAnalysis() { fRunCentralityAnalysis = true; }
  void SetCentralityFile(TString file) { fCentralityFile = file; }
  
  void AnalyseDsPhiPi() { fAnalyseDsPhiPi = true; }
  
  ClassDef(StKFParticleAnalysisMaker,0)   //
};
#endif
// $Log: StKFParticleAnalysisMaker.h,v $
// Revision 1.1  2021/04/19 22:04:11  slan
// update v1
//
