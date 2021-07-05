#ifndef Shift_hh
#define Shift_hh
//
#include "StEpdUtil/StEpdGeom.h"
#include "StEpdUtil/StEpdEpInfo.h"
#include "StMaker.h"
#include <string>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
//
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StPileupUtil;
class TFile;
class TTree;
class TH1;
class TH2;
class TProfile;
class TProfile2D;
class TProfile3D;

#ifndef ST_NO_NAMESPACES
using std::string;
#endif
//
//  The class declaration. It innherits from StMaker.
class Shift : public StMaker {

public:
  Shift( const Char_t *name, StPicoDstMaker *picoMaker, const Char_t *jobid );   // constructor
  virtual ~Shift(){};                                 // destructor

  virtual void Clear(Option_t *option=""); // called after every event to cleanup
  virtual Int_t  Init();                   // called once at the beginning of your job
  virtual Int_t  Make();                   // invoked for every event
  virtual Int_t  Finish();                 // called once at the end


  // My functions
  Int_t GetRunIndex( const Int_t run );
  Int_t Centrality(Int_t gRefMult);
  bool   isGoodTrack(const StPicoTrack *ptrk);
  bool   isGoodTrackSooraj(const StPicoTrack *ptrk);
  bool   isGoodEvent(const StPicoEvent *event);
  bool   isGoodTrigger(const StPicoEvent *event);

private:

  Int_t centnumber;
  Int_t runnumber;
  Int_t eventid;
  float gweight;
  double psi_1_EPD_0;
  double psi_1_EPD_1;
  double psi_1_EPD_2;
  double psi_1_EPD_3;
  double psi_1_EPD_4;
  double psi_1_EPD_5;
  double psi_1_EPD_6;
  double psi_1_EPD;

  StPicoDstMaker *mPicoDstMaker;
  StPicoDst   *mPicoDst;
  StPicoEvent *picoEvent;
  StEpdGeom *mEpdGeom;
  StEpdEpInfo *mEpdEpInfo;
  StPileupUtil* mPileupTool;
  // get recenter par
  TProfile2D *getepd_recen[2][7];
  TProfile2D *gettpc_A_recen[2];
  TProfile2D *gettpc_B_recen[2];
  //  get shift par
  TProfile3D *pp_EPDshiftpar_sin[7];
  TProfile3D *pp_EPDshiftpar_cos[7];
  TProfile3D *pp_TPCshiftpar_Asin;
  TProfile3D *pp_TPCshiftpar_Acos;
  TProfile3D *pp_TPCshiftpar_Bsin;
  TProfile3D *pp_TPCshiftpar_Bcos;

  //Fill Q vector
  TH2D *h_TPCrawep_AQxy;
  TH2D *h_TPCrawep_BQxy;
  TH2D *h_TPCrecen_AQxy;
  TH2D *h_TPCrecen_BQxy;

  //Fill histgram  EP distribution
  TH1D *h_TPCrawep_A[10];
  TH1D *h_TPCrawep_B[10];
  TH1D *h_TPCrecen_A[10];
  TH1D *h_TPCrecen_B[10];
  TH1D *h_TPCshift_A[10];
  TH1D *h_TPCshift_B[10];

  TH1D *h_EPDrawep[10][7];
  TH1D *h_EPDrecen[10][7];
  TH1D *h_EPDshift[10][7];

  TH1D *h_TPCrawep_A_B[10];
  TH1D *h_TPCrecen_A_B[10];
  TH1D *h_TPCshift_A_B[10];
  // 2D correlation
  TH2D *h_TPCrawep_AB[10];
  TH2D *h_TPCrecen_AB[10];
  TH2D *h_TPCshift_AB[10];
  TH2D *h_EPDshift_c[10][7];

  TH2D *h_TPC_A_EPD[10][7];
  TH2D *h_TPC_B_EPD[10][7];

  //resolution
  //cos term
  TProfile *p_r1_tA_tB;
  TProfile *p_r1_tA_eA;
  TProfile *p_r1_tA_eB;
  TProfile *p_r1_tA_eC;
  TProfile *p_r1_tA_eD;
  TProfile *p_r1_tA_eAB;
  TProfile *p_r1_tA_eCD;
  TProfile *p_r1_tA_eABCD;

  TProfile *p_r1_tB_eA;
  TProfile *p_r1_tB_eB;
  TProfile *p_r1_tB_eC;
  TProfile *p_r1_tB_eD;
  TProfile *p_r1_tB_eAB;
  TProfile *p_r1_tB_eCD;
  TProfile *p_r1_tB_eABCD;

  TProfile *p_r1_eA_eB;
  TProfile *p_r1_eA_eC;
  TProfile *p_r1_eA_eD;
  TProfile *p_r1_eB_eC;
  TProfile *p_r1_eB_eD;
  TProfile *p_r1_eC_eD;
  TProfile *p_r1_eAB_eCD;

  TProfile *p_r1_eA_eCD;
  TProfile *p_r1_eB_eCD;
  TProfile *p_r1_eC_eAB;
  TProfile *p_r1_eD_eAB;

  // cos term
  TProfile *p_r1_tA_tB_sin;
  TProfile *p_r1_tA_eA_sin;
  TProfile *p_r1_tA_eB_sin;
  TProfile *p_r1_tA_eC_sin;
  TProfile *p_r1_tA_eD_sin;
  TProfile *p_r1_tA_eAB_sin;
  TProfile *p_r1_tA_eCD_sin;
  TProfile *p_r1_tA_eABCD_sin;

  TProfile *p_r1_tB_eA_sin;
  TProfile *p_r1_tB_eB_sin;
  TProfile *p_r1_tB_eC_sin;
  TProfile *p_r1_tB_eD_sin;
  TProfile *p_r1_tB_eAB_sin;
  TProfile *p_r1_tB_eCD_sin;
  TProfile *p_r1_tB_eABCD_sin;

  TProfile *p_r1_eA_eB_sin;
  TProfile *p_r1_eA_eC_sin;
  TProfile *p_r1_eA_eD_sin;
  TProfile *p_r1_eB_eC_sin;
  TProfile *p_r1_eB_eD_sin;
  TProfile *p_r1_eC_eD_sin;
  TProfile *p_r1_eAB_eCD_sin;

  TTree *psi_tree;
  TFile *File;
  TString mout_shift;

  ClassDef(Shift,0);
};
#endif
