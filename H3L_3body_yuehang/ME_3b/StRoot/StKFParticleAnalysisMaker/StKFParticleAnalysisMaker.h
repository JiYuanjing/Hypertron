// $Id: StKFParticleAnalysisMaker.h,v 1.16 2014/08/06 11:43:53 jeromel Exp $
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
#include "TH2F.h"
#include "StCuts.h"
#include <deque>
#include "TVector3.h"
// #include "StPileupUtil/StPileupUtil.h"

class StKFParticleInterface;
class StKFParticlePerformanceInterface;
class KFParticle;
class StPicoDst;
class StMuDst;
class TNtuple;
class TFile;
class TChain;
class StRefMultCorr;
// class CentralityMaker;

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

  TFile *file_out;
  TTree *lambda_tree;  
  TTree *ks_tree;  
  TTree *cascade_tree;  
  TTree *omega_tree;  
  TTree *htriton_tree; 
  TTree *htriton3_tree; 
  TTree *h4lambda_tree; 
  TTree *he3_tree;
  TTree *h_tree;

  TTree *omega_mc_tree; 
  TTree *lambda_mc_tree; 
  TTree *ks_mc_tree; 
  TTree *cascade_mc_tree;
  TTree *htriton_mc_tree;
  TTree *he3_mc_tree;

  TH1F *hvtx;
  TH1F *hvtxgood;
  TH2F *hvtx_xy;
  TH1F *hrefmult;
  TH1F *hCent;
  TH1F *hCentWt;
  TH1F *wrefmult;

  int notbadrun;
  int beventid;
  int brunid;

  bool _fill_lambda_tree;
  bool _fill_he3_tree;

  int brefmult;
  int btofmult;

  int bparticleid;
  float bparticlemass;

  float ld_chi2topo;
  float ld_chi2primary;
  float ld_chi2ndf;
  float ld_ldl;
  float ld_l;
  float ld_dl;

  float  ht_chi2;
  float  ht_NDF;
  float ht_chi2topo;
  float ht_chi2ndf;
  float ht_ldl;
  float ht_l;
  float ht_dl;

  float chi2primary_he;
  float chi2primary_proton;
  float chi2primary_pi;

  float hl_chi2;
  float hl_NDF;
  float hl_chi2topo;
  float hl_chi2ndf;
  float hl_ldl;
  float hl_l;
  float hl_dl;
  float chi2primary_h4;
  float dca_proton;
  float dca_he;
  float dca_pi;
  float dca_pion2;
  float dca_proton2;
  float dca_deuteron2;
  float om_chi2topo;
  float om_chi2ndf;
  float om_ldl;
  float om_l;
  float om_dl;
  float chi2primary_om_proton;
  float chi2primary_om_pi;
  float chi2primary_om_bach;
  float chi2primary_om_ld;
  float om_ld_chi2topo;
  float om_ld_chi2ndf;
  float om_ld_ldl;
  float om_ld_l;

  float xi_chi2topo;
  float xi_chi2ndf;
  float xi_ldl;
  float xi_l;
  float xi_dl;
  float chi2primary_xi_proton;
  float chi2primary_xi_pi;
  float chi2primary_xi_bach;
  float chi2primary_xi_ld;
  float chi2primary_d;
  float xi_ld_chi2topo;
  float xi_ld_chi2ndf;
  float xi_ld_ldl;
  float xi_ld_l;

  float bpionm2;
  float bprotonm2;
  float bdm2;

  float bpionnsigma;
  float bprotonsigma;
  int nhits_om_proton;
  int nhits_ld_proton;
  int nhits_om_pi;
  int nhits_ld_pi;
  int nhits_om_bach;

  float mass_01, mass_01_err;
  float mass_02, mass_02_err;
  float mass_12, mass_12_err;

  float v_lambda_mass_0;
  float v_lambda_ldl_0;
  float v_lambda_l_0;
  float v_lambda_chi2primary_0;
  float v_01_chi2primary;
  float v_02_chi2primary;
  float v_12_chi2primary;
float v_01_dca;
float v_02_dca;
float v_12_dca;
float v_012_dca;

float v_01_chi2ndf;
float v_02_chi2ndf;
float v_12_chi2ndf;

  int nhits_deuteron;
  float dca_deuteron; 
  float dca_pion; 
  float dedx_om_proton;
  float dedx_om_pi;
  float dedx_om_bach;
  int nhits_proton;
  int nhits_pion;
  int nhitsdedx_pion;
  float bdedx;
  float bdca;
  int bnhits;
  int nhits_he;
  int nhitsdedx_he;
  int nhits_pi;
  int nhitsdedx_pi;
  float ld_bdfvtx;
  float ld_bdfvtx2;
  float ld_bdfvtx_xy;
  float ld_bdfvtxdev_xy;
  float ld_lifetime;
  float ht_bdfvtx;
  float ht_bdfvtx2;
  float ht_lifetime;

  float v_01_pvdca;
  float v_02_pvdca;
  float v_12_pvdca;
  float v_01_x;
  float v_01_y;
  float v_01_z;
  float v_02_x;
  float v_02_y;
  float v_02_z;
  float v_12_x;
  float v_12_y;
  float v_12_z;

  float bmcpx;
  float bmcpy;
  float bmcpz;
float b0mcpx;
float b0mcpy;
float b0mcpz;
float b1mcpx;
float b1mcpy;
float b1mcpz;
float b2mcpx;
float b2mcpy;
float b2mcpz;

  float bx;
  float by;
  float bz;
  float bmcx;
  float bmcy;
  float bmcz;	
  float bpx;
  float bpy;
  float bpz;
  float bdl;
  float bzdeuteron;

  float bdpx;
  float bdpy;
  float bdpz;

  float bVx;
  float bVy; 
  float bVr;
  float bVz;
  float bVrerr;
  float bVxerr;
  float bVyerr;
  float bVzerr;

  int bbachid;
  int bpionid;
  int bprotonid;
  int ht_ndaughters;
  int bismc;       
  int bisMix;
  int cent9;
  double refmultcor;
  double reweight;

  float dedx_pi;
  float dedx_he;

  float bbachpx;
  float bbachpy;
  float bbachpz;
  float bbachmass;

  float b2mcrawpx;
  float b2mcrawpy;
  float b2mcrawpz;

  float b0mcrawpx;
  float b0mcrawpy;
  float b0mcrawpz;
 
  float b1mcrawpx;
  float b1mcrawpy;
  float b1mcrawpz;

  float px_pi;
  float py_pi;
  float pz_pi;
  float px_he;
  float py_he;
  float pz_he;

  float bpionpx;
  float bpionpy;
  float bpionpz;
  float bpionmass;
  float gweight;
  float bprotonpx;
  float bprotonpy;
  float bprotonpz;
  float bprotonmass;
  float bpl;
  float bmcrawpx;
  float bmcrawpy;
  float bmcrawpz;
  int bmcidvx;
  int bmcidvxend;
  int bmcrefmult;
  int bmcparticleid;
  float bmcrawl;
  float bmcrawpl;
  float bmcpl;
  float bmcl;
        
  double psi_1_EPD_0;
  double psi_1_EPD_1;
  double psi_1_EPD_2;
  double psi_1_EPD_3; 
  double psi_1_EPD_4; 
  double psi_1_EPD_5; 
  double psi_1_EPD_6; 
  double psi_1_EPD; 

  int countrefmult;
  int FXTMult;
  int FXTMult2;

  //Centrality and flow
  Bool_t fFlowAnalysis;
  TChain* fFlowChain;
  int fFlowRunId;
  int fFlowEventId;
  int fCentrality;
  int fMixEventBufferSize;
  std::vector<TString> fFlowFiles;
  std::map<long, int> fFlowMap;
  
  bool fRunCentralityAnalysis;
  StRefMultCorr *fRefmultCorrUtil;
  // StPileupUtil* mPileupTool;
  
  TString fCentralityFile;

  bool fMixEvent;

  // std::vector <StMixEvent> bMixEventBuffer[nCentralityBins][nEPbins];
  std::deque <StMixEvent> bMixEventBuffer[Cuts::nCentralityBins][Cuts::nEPbins][Cuts::nVtxbins];
  int EPbin, Vtxbin;
 
 //functions 
  int CheckIfBadRun(int fsnn, int brunid);
  int PassTrigger(StPicoDst* fPicoDst, int fsnn);
  bool fAnalyseDsPhiPi;

  void GetDaughterParameters(const int iReader, int& iDaughterTrack, int& iDaughterParticle, KFParticle& particle);
  void GetParticleParameters(const int iReader, KFParticle& particle);
  long  GetUniqueEventId(const int iRun, const int iEvent) const;
  
  int GetTMVACentralityBin(int iReader, int centrality);
  int GetTMVAPtBin(int iReader, double pt);
  void SetTMVACentralityBins(int iReader, TString bins);
  void SetTMVAPtBins(int iReader, TString bins);
  void SetTMVABins(int iReader, TString centralityBins="-1:1000", TString ptBins="-1.:1000.");
 
  int getEventPlaneBin(double ep)
  {
    //split ep in nEPbins for mix
    double p=3.14159265;
    double edge[Cuts::nEPbins+1];
    edge[Cuts::nEPbins]=p;
   
    for (int i=0;i<Cuts::nEPbins;i++) {edge[i]=i*p*2.0/Cuts::nEPbins-p; }
    int idx=-999;
    for (int i=0;i<Cuts::nEPbins;i++)
    {
      if (ep>=edge[i] && ep<edge[i+1]) {idx=i;break;}
    }
    if (idx>=0 && idx<Cuts::nEPbins) return idx;
    else { cout <<"event plane angle out of range! ep=" <<ep << endl; return -999;} 
  } 
  int getVtxBin(float vx, float vy, float vz);

void FillDaughterInfo(KFParticle& daughter, int trackId, int PDG, float &chi2primary, int &nhits, float &dca,float& sigma, float& bm2);
bool FillThreeDaughtersMix(KFParticle pion, KFParticle proton, KFParticle deuteron, TVector3 dV_pi, TVector3 dV_p, TVector3 dV_d );

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
    static const char cvs[]="Tag $Name:  $ $Id: StKFParticleAnalysisMaker.h,v 1.0 2017/10/07 11:43:53 mzyzak Exp $ built " __DATE__ " " __TIME__ ; 
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
  ///mix event 
  void RunMixEvent(int buffersize)       
  { 
      fMixEvent= true; 
      // fFlowAnalysis = true;
      fRunCentralityAnalysis = true;
      fMixEventBufferSize = buffersize;
  }

  Int_t Centrality(int gRefMult );
  Int_t Centrality16(int gRefMult );
  Double_t FitWeight(Double_t refMult);
  
  void RunCentralityAnalysis() { fRunCentralityAnalysis = true; }
  void SetCentralityFile(TString file) { fCentralityFile = file; }
  
  void AnalyseDsPhiPi() { fAnalyseDsPhiPi = true; }
  
  ClassDef(StKFParticleAnalysisMaker,0)   //
};
#endif
// $Log: StKFParticleAnalysisMaker.h,v $
