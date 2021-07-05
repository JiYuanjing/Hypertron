#include <TFile.h>
#include <TTree.h>
#include <StMessMgr.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TROOT.h>
#include <TRandom.h>
#include <StThreeVectorF.hh>
#include <StHelix.hh>
#include <TLorentzVector.h>

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoEpdHit.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "StarClassLibrary/StLorentzVectorF.hh"
#include "StEpdUtil/StEpdGeom.h"
#include "StEpdUtil/StEpdEpFinder.h"
//#include "StPileupUtil/StPileupUtil.h"
#include "phys_constants.h"

#include "../run/run.h"
#include "../run/badrun.h"
#include "Shift.h"
//#include "../run/nsigmap_mean.h"
#include "constant.h"

ClassImp(Shift)

//__________________________________________________________________________________
Shift::Shift( const Char_t *name, StPicoDstMaker *picoMaker, const Char_t *jobid ) : StMaker(name) {
  mPicoDstMaker = picoMaker;
  mPicoDst = 0;
  
  mout_shift=Form("%s_v2.root", jobid);
}

//__________________________________________________________________________________
Int_t Shift::Init() {
  cout << "Init" << endl;
  const float pi = acos(-1.0);
  mEpdGeom = new StEpdGeom();
//  mPileupTool = new StPileupUtil();
//  mPileupTool->init();
  cout << "Define Histograms" << endl;
  //===============================
  //  Define Histograms
  //===============================
  
  File=new TFile(mout_shift.Data(),"RECREATE");
  
  psi_tree = new TTree("psi_tree","psi_tree");
  psi_tree -> Branch("isPileUp",&isPileUp,"isPileUp/O");
  psi_tree -> Branch("runnumber",&runnumber,"runnumber/I");
  psi_tree -> Branch("eventid",&eventid,"eventid/I");
  psi_tree -> Branch("centnumber",&centnumber,"centnumber/I");
  psi_tree -> Branch("refMultPrim",&refMultPrim,"refMultPrim/I");
  psi_tree -> Branch("psi_1_EPD_0",&psi_1_EPD_0,"psi_1_EPD_0/D");
  psi_tree -> Branch("psi_1_EPD_1",&psi_1_EPD_1,"psi_1_EPD_1/D");
  psi_tree -> Branch("psi_1_EPD_2",&psi_1_EPD_2,"psi_1_EPD_2/D");
  psi_tree -> Branch("psi_1_EPD_3",&psi_1_EPD_3,"psi_1_EPD_3/D");
  psi_tree -> Branch("psi_1_EPD_4",&psi_1_EPD_4,"psi_1_EPD_4/D");
  psi_tree -> Branch("psi_1_EPD_5",&psi_1_EPD_5,"psi_1_EPD_5/D");
  psi_tree -> Branch("psi_1_EPD_6",&psi_1_EPD_6,"psi_1_EPD_6/D");
  psi_tree -> Branch("psi_1_EPD",&psi_1_EPD,"psi_1_EPD/D");
  psi_tree -> Branch("gweight",&gweight,"gweight/F");
  
  h_TPCrawep_AQxy = new TH2D("TPCrawep_AQxy","",100,-2,2,100,-2,2);
  h_TPCrawep_BQxy = new TH2D("TPCrawep_BQxy","",100,-2,2,100,-2,2);
  h_TPCrecen_AQxy = new TH2D("TPCrecen_AQxy","",100,-2,2,100,-2,2);
  h_TPCrecen_BQxy = new TH2D("TPCrecen_BQxy","",100,-2,2,100,-2,2);
  
  for(int i=0; i<10; i++)
  {
    h_TPCrawep_A[i] =new TH1D(Form("TPCrawep_A_cent%d",i),"",314,-pi, pi);
    h_TPCrawep_B[i] =new TH1D(Form("TPCrawep_B_cent%d",i),"",314,-pi, pi);
    h_TPCrecen_A[i] =new TH1D(Form("TPCrecen_A_cent%d",i),"",314,-pi, pi);
    h_TPCrecen_B[i] =new TH1D(Form("TPCrecen_B_cent%d",i),"",314,-pi, pi);
    h_TPCshift_A[i] =new TH1D(Form("TPCshift_A_cent%d",i),"",314,-pi, pi);
    h_TPCshift_B[i] =new TH1D(Form("TPCshift_B_cent%d",i),"",314,-pi, pi);
  }
  for(int i=0; i<10; i++)
  {
    for(int j=0; j<7; j++)
    {
      //  raw EP distribution
      h_EPDrawep[i][j] =new TH1D(Form("EPDrawep_ring%d_cent%d",j,i),"",314,-pi, pi);
      // recenter EP distribution
      h_EPDrecen[i][j] =new TH1D(Form("EPDrecen_ring%d_cent%d",j,i),"",314,-pi, pi);
      //  shift EP dsitribution
      h_EPDshift[i][j] =new TH1D(Form("EPDshift_ring%d_cent%d",j,i),"",314,-pi, pi);
      
      h_TPC_A_EPD[i][j] = new TH2D(Form("TPC_A_EPD_ring%d_cent%d",j,i),"",32,-pi, pi,32,-pi, pi);
      h_TPC_B_EPD[i][j] = new TH2D(Form("TPC_B_EPD_ring%d_cent%d",j,i),"",32,-pi, pi,32,-pi, pi);
    }
  }
  
  
  // cos term
  p_r1_tA_tB = new TProfile("p_r1_tA_tB","",9,-0.5,8.5,-1.0,1.0,"");
  
  p_r1_tA_eA = new TProfile("p_r1_tA_eA","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_tA_eB = new TProfile("p_r1_tA_eB","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_tA_eC = new TProfile("p_r1_tA_eC","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_tA_eD = new TProfile("p_r1_tA_eD","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_tA_eAB = new TProfile("p_r1_tA_eAB","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_tA_eCD = new TProfile("p_r1_tA_eCD","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_tA_eABCD = new TProfile("p_r1_tA_eABCD","",9,-0.5,8.5,-1.0,1.0,"");
  
  p_r1_tB_eA = new TProfile("p_r1_tB_eA","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_tB_eB = new TProfile("p_r1_tB_eB","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_tB_eC = new TProfile("p_r1_tB_eC","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_tB_eD = new TProfile("p_r1_tB_eD","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_tB_eAB = new TProfile("p_r1_tB_eAB","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_tB_eCD = new TProfile("p_r1_tB_eCD","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_tB_eABCD = new TProfile("p_r1_tB_eABCD","",9,-0.5,8.5,-1.0,1.0,"");
  
  p_r1_eA_eB = new TProfile("p_r1_eA_eB","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_eA_eC = new TProfile("p_r1_eA_eC","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_eA_eD = new TProfile("p_r1_eA_eD","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_eB_eC = new TProfile("p_r1_eB_eC","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_eB_eD = new TProfile("p_r1_eB_eD","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_eC_eD = new TProfile("p_r1_eC_eD","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_eAB_eCD = new TProfile("p_r1_eAB_eCD","",9,-0.5,8.5,-1.0,1.0,"");
  
  p_r1_eA_eCD = new TProfile("p_r1_eA_eCD","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_eB_eCD = new TProfile("p_r1_eB_eCD","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_eC_eAB = new TProfile("p_r1_eC_eAB","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_eD_eAB = new TProfile("p_r1_eD_eAB","",9,-0.5,8.5,-1.0,1.0,"");
  
  // sin term
  p_r1_tA_tB_sin = new TProfile("p_r1_tA_tB_sin","",9,-0.5,8.5,-1.0,1.0,"");
  
  p_r1_tA_eA_sin = new TProfile("p_r1_tA_eA_sin","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_tA_eB_sin = new TProfile("p_r1_tA_eB_sin","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_tA_eC_sin = new TProfile("p_r1_tA_eC_sin","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_tA_eD_sin = new TProfile("p_r1_tA_eD_sin","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_tA_eAB_sin = new TProfile("p_r1_tA_eAB_sin","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_tA_eCD_sin = new TProfile("p_r1_tA_eCD_sin","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_tA_eABCD_sin = new TProfile("p_r1_tA_eABCD_sin","",9,-0.5,8.5,-1.0,1.0,"");
  
  p_r1_tB_eA_sin = new TProfile("p_r1_tB_eA_sin","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_tB_eB_sin = new TProfile("p_r1_tB_eB_sin","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_tB_eC_sin = new TProfile("p_r1_tB_eC_sin","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_tB_eD_sin = new TProfile("p_r1_tB_eD_sin","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_tB_eAB_sin = new TProfile("p_r1_tB_eAB_sin","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_tB_eCD_sin = new TProfile("p_r1_tB_eCD_sin","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_tB_eABCD_sin = new TProfile("p_r1_tB_eABCD_sin","",9,-0.5,8.5,-1.0,1.0,"");
  
  p_r1_eA_eB_sin = new TProfile("p_r1_eA_eB_sin","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_eA_eC_sin = new TProfile("p_r1_eA_eC_sin","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_eA_eD_sin = new TProfile("p_r1_eA_eD_sin","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_eB_eC_sin = new TProfile("p_r1_eB_eC_sin","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_eB_eD_sin = new TProfile("p_r1_eB_eD_sin","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_eC_eD_sin = new TProfile("p_r1_eC_eD_sin","",9,-0.5,8.5,-1.0,1.0,"");
  p_r1_eAB_eCD_sin = new TProfile("p_r1_eAB_eCD_sin","",9,-0.5,8.5,-1.0,1.0,"");
  
  for(int i=0; i<10; i++)
  {
    h_TPCrawep_AB[i]=new TH2D(Form("TPCrawep_AB_cent%d",i),"",32,-pi, pi,32,-pi, pi);
    h_TPCrecen_AB[i]=new TH2D(Form("TPCrecen_AB_cent%d",i),"",32,-pi, pi,32,-pi, pi);
    h_TPCshift_AB[i]=new TH2D(Form("TPCshift_AB_cent%d",i),"",32,-pi, pi,32,-pi, pi);
    
    h_EPDshift_c[i][0]=new TH2D(Form("EPDshift_01_cent%d",i),"",32,-pi, pi,32,-pi, pi);
    h_EPDshift_c[i][1]=new TH2D(Form("EPDshift_02_cent%d",i),"",32,-pi, pi,32,-pi, pi);
    h_EPDshift_c[i][2]=new TH2D(Form("EPDshift_03_cent%d",i),"",32,-pi, pi,32,-pi, pi);
    h_EPDshift_c[i][3]=new TH2D(Form("EPDshift_12_cent%d",i),"",32,-pi, pi,32,-pi, pi);
    h_EPDshift_c[i][4]=new TH2D(Form("EPDshift_13_cent%d",i),"",32,-pi, pi,32,-pi, pi);
    h_EPDshift_c[i][5]=new TH2D(Form("EPDshift_23_cent%d",i),"",32,-pi, pi,32,-pi, pi);
    h_EPDshift_c[i][6]=new TH2D(Form("EPDshift_1234_cent%d",i),"",32,-pi, pi,32,-pi, pi);
  }
  
  
  // connecting recentering and shift parameter
  //TFile *re_file = TFile::Open("/star/u/slan/pwg/fixtarget/FXT3_85/EPDD/EPrecon/recenter/recenterpar.root");
  TFile *re_file = TFile::Open("/star/u/slan/3gev/EPDD/pikpflow/SL20d/eventplane/parameter/recenterpar.root");
  
  gettpc_A_recen[0]     = (TProfile2D*)re_file -> Get("TPCqx_A_recen");
  gettpc_A_recen[1]     = (TProfile2D*)re_file -> Get("TPCqy_A_recen");
  gettpc_B_recen[0]     = (TProfile2D*)re_file -> Get("TPCqx_B_recen");
  gettpc_B_recen[1]     = (TProfile2D*)re_file -> Get("TPCqy_B_recen");
  for(int i=0; i<7; i++){
    getepd_recen[0][i]     = (TProfile2D*)re_file -> Get(Form("EPDQx_ring%d_recen",i));
    getepd_recen[1][i]     = (TProfile2D*)re_file -> Get(Form("EPDQy_ring%d_recen",i));
  }
  
  //TFile *shift_file = TFile::Open("/star/u/slan/pwg/fixtarget/FXT3_85/EPDD/EPrecon/shiftpar/shiftpar.root");
  TFile *shift_file = TFile::Open("/star/u/slan/3gev/EPDD/pikpflow/SL20d/eventplane/parameter/shiftpar.root");
  for(int i=0; i<7; i++){
    pp_EPDshiftpar_sin[i]  = (TProfile3D*)shift_file->Get(Form("EPDshiftpar_sin_ring%d",i));
    pp_EPDshiftpar_cos[i]  = (TProfile3D*)shift_file->Get(Form("EPDshiftpar_cos_ring%d",i));
  }
  
  pp_TPCshiftpar_Asin  = (TProfile3D*)shift_file->Get("TPCshiftpar_Asin");
  pp_TPCshiftpar_Bsin  = (TProfile3D*)shift_file->Get("TPCshiftpar_Bsin");
  pp_TPCshiftpar_Acos  = (TProfile3D*)shift_file->Get("TPCshiftpar_Acos");
  pp_TPCshiftpar_Bcos  = (TProfile3D*)shift_file->Get("TPCshiftpar_Bcos");
  
  cout << "End of Histograms" << endl;
  return kStOK;
}
//__________________________________________________________________________________
void Shift::Clear(Option_t *opt)
{
  StMaker::Clear();
}

//__________________________________________________________________________________
Int_t Shift::Finish() {
  cout << "Shift::Finish()\n";
  //===============================
  //  Write Histograms
  //===============================
  
  File->cd();
  psi_tree->Write();
  
  h_TPCrawep_AQxy->Write();
  h_TPCrawep_BQxy->Write();
  h_TPCrecen_AQxy->Write();
  h_TPCrecen_BQxy->Write();
  
  for(int i=0; i<10; i++){
    for(int j=0; j<7; j++){
      //EPD EP
      h_EPDrawep[i][j]->Write();
      h_EPDrecen[i][j]->Write();
      h_EPDshift[i][j]->Write();
      
      h_TPC_A_EPD[i][j]->Write();
      h_TPC_B_EPD[i][j]->Write();
    }
    h_TPCrawep_A[i]->Write();
    h_TPCrawep_B[i]->Write();
    h_TPCrecen_A[i]->Write();
    h_TPCrecen_B[i]->Write();
    h_TPCshift_A[i]->Write();
    h_TPCshift_B[i]->Write();
    
    h_TPCrawep_AB[i]->Write();
    h_TPCrecen_AB[i]->Write();
    h_TPCshift_AB[i]->Write();
  }
  
  for(int i=0; i<10; i++){
    for(int j=0; j<7; j++){
      h_EPDshift_c[i][j]->Write();
    }
  }
  
  // cos term
  p_r1_tA_tB->Write();
  p_r1_tA_eA->Write();
  p_r1_tA_eB->Write();
  p_r1_tA_eC->Write();
  p_r1_tA_eD->Write();
  p_r1_tA_eAB->Write();
  p_r1_tA_eCD->Write();
  p_r1_tA_eABCD->Write();
  
  p_r1_tB_eA->Write();
  p_r1_tB_eB->Write();
  p_r1_tB_eC->Write();
  p_r1_tB_eD->Write();
  p_r1_tB_eAB->Write();
  p_r1_tB_eCD->Write();
  p_r1_tB_eABCD->Write();
  
  p_r1_eA_eB->Write();
  p_r1_eA_eC->Write();
  p_r1_eA_eD->Write();
  p_r1_eB_eC->Write();
  p_r1_eB_eD->Write();
  p_r1_eC_eD->Write();
  p_r1_eAB_eCD->Write();
  
  p_r1_eA_eCD->Write();
  p_r1_eB_eCD->Write();
  p_r1_eC_eAB->Write();
  p_r1_eD_eAB->Write();
  
  // sin term
  p_r1_tA_tB_sin->Write();
  p_r1_tA_eA_sin->Write();
  p_r1_tA_eB_sin->Write();
  p_r1_tA_eC_sin->Write();
  p_r1_tA_eD_sin->Write();
  p_r1_tA_eAB_sin->Write();
  p_r1_tA_eCD_sin->Write();
  p_r1_tA_eABCD_sin->Write();
  
  p_r1_tB_eA_sin->Write();
  p_r1_tB_eB_sin->Write();
  p_r1_tB_eC_sin->Write();
  p_r1_tB_eD_sin->Write();
  p_r1_tB_eAB_sin->Write();
  p_r1_tB_eCD_sin->Write();
  p_r1_tB_eABCD_sin->Write();
  
  p_r1_eA_eB_sin->Write();
  p_r1_eA_eC_sin->Write();
  p_r1_eA_eD_sin->Write();
  p_r1_eB_eC_sin->Write();
  p_r1_eB_eD_sin->Write();
  p_r1_eC_eD_sin->Write();
  p_r1_eAB_eCD_sin->Write();
  
  
  return kStOK;
}
//__________________________________________________________________________________
Int_t Shift::GetRunIndex( const Int_t run ) {
  Int_t runindex = -999;
  for(Int_t i=0; i<nrun; i++){
    if(run==numbers[i]) runindex = i;
  }
  return runindex;
}
//---------------------------------------------------------------------------------

Int_t Shift::Centrality(int gRefMult )
{
  int centrality;
  int centFull[10]={5, 9, 16, 26, 41, 60, 86, 119, 142, 195};
  if (gRefMult>=centFull[8]) centrality=8;
  else if (gRefMult>=centFull[7] && gRefMult<centFull[8] ) centrality=7;
  else if (gRefMult>=centFull[6] && gRefMult<centFull[7] ) centrality=6;
  else if (gRefMult>=centFull[5] && gRefMult<centFull[6] ) centrality=5;
  else if (gRefMult>=centFull[4] && gRefMult<centFull[5] ) centrality=4;
  else if (gRefMult>=centFull[3] && gRefMult<centFull[4] ) centrality=3;
  else if (gRefMult>=centFull[2] && gRefMult<centFull[3] ) centrality=2;
  else if (gRefMult>=centFull[1] && gRefMult<centFull[2] ) centrality=1;
  else if (gRefMult>=centFull[0] && gRefMult<centFull[1] ) centrality=0;
  else centrality = 9;
  
  return centrality;
}

//---------------------------------------------------------------------------------
//__________________________________________________________________________________
Int_t Shift::Make() {
  //Begining of Event loop
  
  //------------------------------------------------------------------
  if(!mPicoDstMaker) {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }
  mPicoDst = mPicoDstMaker->picoDst();
  if(!mPicoDst) {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    return kStWarn;
  }
  picoEvent = (StPicoEvent*)mPicoDst->event();
  if( !picoEvent ){
    LOG_WARN << " No PicoEvent! Skip! " << endm;
    return kStWarn;
  }
  //------------------------------------------------------------------
  if(!isGoodTrigger(picoEvent)) return 0;
  if(!isGoodEvent(picoEvent)) return 0;
  
  runnumber = picoEvent->runId();
  eventid = picoEvent->eventId();
  int runindex = GetRunIndex(runnumber);
  // bad run removement
  for(int ii=0; ii<24; ii++)
  {
    if(runnumber == badrun[ii]) return 0;
  }
  
  
  int countrefmult=0;
  
  const Int_t nTrack = mPicoDst->numberOfTracks();
  
  // primary track loop for determine refmult ----------------------------------------------
  
  for (Int_t itr=0;itr<nTrack;itr++) {
    const StPicoTrack *ptrk = (StPicoTrack*)mPicoDst->track(itr);
    
    if(!ptrk)  continue;
    if(!ptrk->isPrimary())  continue;  // now selecting primary tracks
    countrefmult++;
  }
  // End of track loop ---------------------------------------------
  //---------------------Centrality----------------------------------
  centnumber = Centrality(countrefmult);
  gweight = 1.0;
  if(centnumber > 8 || centnumber < 0) return kStOK;
  //   cout<<"centnumber:"<<centnumber<<" "<<countrefmult<<endl;
  //    if(centnumber > 8 || centnumber < 0) return kStOK;
  
  TVector3 pVertex = picoEvent->primaryVertex();
  
  // pileup
  //  mPileupTool->initEvent(mPicoDst);
  //  countrefmult = mPileupTool->get_refMultPrim();
  //  centnumber  = mPileupTool->get_centrality9();
  //  gweight     = mPileupTool->get_centralityWeight();
  //  refMultPrim = mPileupTool->get_refMultPrim();
  //  isPileUp    = mPileupTool->isPileupEPD();
  
  //  if(isPileUp) return 0;
  //  if(centnumber > 8 || centnumber < 0) return 0;
  
  double Qx_rawep_TPC_A=0.0, Qy_rawep_TPC_A=0.0;
  double Qx_rawep_TPC_B=0.0, Qy_rawep_TPC_B=0.0;
  double Qx_recen_TPC_A=0.0, Qy_recen_TPC_A=0.0;
  double Qx_recen_TPC_B=0.0, Qy_recen_TPC_B=0.0;
  double weight_TPC_A = 0.0, weight_TPC_B = 0.0;
  
  double psi1rawep_A_TPC=0.0, psi1rawep_B_TPC=0.0;
  double psi1recen_A_TPC=0.0, psi1recen_B_TPC=0.0;
  double psi1shift_A_TPC=0.0, psi1shift_B_TPC=0.0;
  
  // primary track loop for determine refmult ----------------------------------------------
  for (Int_t itr=0;itr<nTrack;itr++) {
    const StPicoTrack *ptrk = (StPicoTrack*)mPicoDst->track(itr);
    
    if(!ptrk)  continue;
    if(!ptrk->isPrimary())  continue;  // now selecting primary tracks
    if(!isGoodTrack(ptrk))  continue;
    
    Float_t pt  = ptrk->pMom().Perp(); // zero for global tracks
    Float_t mom = ptrk->pMom().Mag();
    Float_t phi = ptrk->pMom().Phi();
    Float_t eta = ptrk->pMom().PseudoRapidity();
    Float_t dca = ptrk->gDCA( pVertex ).Mag();
    if(! (pt > 0.2 && pt < 2.0)) continue;  // using pt[0.2,2.0]GeV/c particles to reconstruct 2nd EP
    //cout << "pt = " << pt << endl;
    if(eta > -2.0 && eta < -1.1)
    {
      Qx_rawep_TPC_A += pt*(cos(1.0*phi));
      Qy_rawep_TPC_A += pt*(sin(1.0*phi));
      Qx_recen_TPC_A += pt*(cos(1.0*phi) - gettpc_A_recen[0]->GetBinContent(centnumber+1, runindex+1));
      Qy_recen_TPC_A += pt*(sin(1.0*phi) - gettpc_A_recen[1]->GetBinContent(centnumber+1, runindex+1));
      weight_TPC_A++;
    }
    
    else if(eta > -1.0 && eta < 0)
    {
      Qx_rawep_TPC_B += pt*(cos(1.0*phi));
      Qy_rawep_TPC_B += pt*(sin(1.0*phi));
      Qx_recen_TPC_B += pt*(cos(1.0*phi) - gettpc_B_recen[0]->GetBinContent(centnumber+1, runindex+1));
      Qy_recen_TPC_B += pt*(sin(1.0*phi) - gettpc_B_recen[1]->GetBinContent(centnumber+1, runindex+1));
      weight_TPC_B++;
    }
    else continue;
  }
  if(Qy_rawep_TPC_A ==0. || Qx_rawep_TPC_A==0. || Qy_rawep_TPC_B==0. || Qx_rawep_TPC_B==0.) return 0;
  if(Qy_recen_TPC_A ==0. || Qx_recen_TPC_A==0. || Qy_recen_TPC_B==0. || Qx_recen_TPC_B==0.) return 0;
  if(weight_TPC_A == 0 || weight_TPC_B == 0) return 0;
  Qx_rawep_TPC_A /= weight_TPC_A;
  Qy_rawep_TPC_A /= weight_TPC_A;
  Qx_recen_TPC_A /= weight_TPC_A;
  Qy_recen_TPC_A /= weight_TPC_A;
  Qx_rawep_TPC_B /= weight_TPC_B;
  Qy_rawep_TPC_B /= weight_TPC_B;
  Qx_recen_TPC_B /= weight_TPC_B;
  Qy_recen_TPC_B /= weight_TPC_B;
  
  psi1rawep_A_TPC = (1.0/1.0) * atan2(Qy_rawep_TPC_A, Qx_rawep_TPC_A);
  psi1rawep_B_TPC = (1.0/1.0) * atan2(Qy_rawep_TPC_B, Qx_rawep_TPC_B);
  
  psi1recen_A_TPC = (1.0/1.0) * atan2(Qy_recen_TPC_A, Qx_recen_TPC_A);
  psi1recen_B_TPC = (1.0/1.0) * atan2(Qy_recen_TPC_B, Qx_recen_TPC_B);
  
  //if(psi1rawep_A_TPC < -pi) {psi1rawep_A_TPC += 2*pi;}
  //if(psi1rawep_B_TPC < -pi) {psi1rawep_B_TPC += 2*pi;}
  //if(psi1rawep_A_TPC >  pi) {psi1rawep_A_TPC -= 2*pi;}
  //if(psi1rawep_B_TPC >  pi) {psi1rawep_B_TPC -= 2*pi;}
  
  //if(psi1recen_A_TPC < -pi) {psi1recen_A_TPC += 2*pi;}
  //if(psi1recen_B_TPC < -pi) {psi1recen_B_TPC += 2*pi;}
  //if(psi1recen_A_TPC >  pi) {psi1recen_A_TPC -= 2*pi;}
  //if(psi1recen_B_TPC >  pi) {psi1recen_B_TPC -= 2*pi;}
  
  psi1shift_A_TPC = psi1recen_A_TPC;
  psi1shift_B_TPC = psi1recen_B_TPC;
  
  double shiftsin_TPC_A=0, shiftcos_TPC_A=0;
  double shiftsin_TPC_B=0, shiftcos_TPC_B=0;
  
  for(int i=0; i<20; i++)
  {
    psi1shift_A_TPC += (2.0/(i+1)) * (-pp_TPCshiftpar_Asin->GetBinContent(centnumber+1, i+1, runindex+1)*cos(1.0*(i+1)*psi1recen_A_TPC) + pp_TPCshiftpar_Acos->GetBinContent(centnumber+1, i+1, runindex+1)*sin(1.0*(i+1)*psi1recen_A_TPC));
    psi1shift_B_TPC += (2.0/(i+1)) * (-pp_TPCshiftpar_Bsin->GetBinContent(centnumber+1, i+1, runindex+1)*cos(1.0*(i+1)*psi1recen_B_TPC) + pp_TPCshiftpar_Bcos->GetBinContent(centnumber+1, i+1, runindex+1)*sin(1.0*(i+1)*psi1recen_B_TPC));
  }
  
  if(psi1shift_A_TPC < -pi)   psi1shift_A_TPC += 2*pi;
  if(psi1shift_B_TPC < -pi)   psi1shift_B_TPC += 2*pi;
  if(psi1shift_A_TPC >  pi)   psi1shift_A_TPC -= 2*pi;
  if(psi1shift_B_TPC >  pi)   psi1shift_B_TPC -= 2*pi;
  
  // TPC end
  
  //-----------------Get EPD information------------------------------
  Int_t nepdHits = mPicoDst->numberOfEpdHits();
  
  StPicoEpdHit *epdHit;
  TVector3 StraightLine_center;
  TVector3 StraightLine_random;
  double phi_epd_center;
  double phi_epd_random;
  double eta_epd_center;
  double eta_epd_random;
  
  double mip;
  double TileWeight           = {0};
  
  double Qx_rawep_EPD[7]={0.0}, Qy_rawep_EPD[7]={0.0};
  double Qx_recen_EPD[7]={0.0}, Qy_recen_EPD[7]={0.0};
  
  double psi1rawep_EPD[7]={0.0};
  double psi1recen_EPD[7]={0.0};
  double psi1shift_EPD[7]={0.0};
  double weight_EPD[7]={0};
  
  
  //if(nepdHits < 75) return 0;
  for(Int_t iHit=0; iHit<nepdHits; iHit++){
    epdHit = mPicoDst->epdHit(iHit);
    mip = epdHit->nMIP();
    int iring = epdHit->row() -1;//(1~16)-1 -> 0-15
    if( !epdHit) continue;
    if(epdHit->id() > 0 ) continue;                  // unique tile identifier, absolute value is 100*position+tile, sign is +1/-1 for West/East
    
    int ringgroup = mEpdEpInfo->RingGroup(iring);   //0: 0-7, 1: 8-15    0-> inner most
    if(ringgroup == -1) continue;
    
    StraightLine_center = mEpdGeom->TileCenter(epdHit->id())        - picoEvent->primaryVertex();
    StraightLine_random = mEpdGeom->RandomPointOnTile(epdHit->id()) - picoEvent->primaryVertex();
    //StraightLine_center = mEpdGeom->TileCenter(epdHit->id())       ;
    //StraightLine_random = mEpdGeom->RandomPointOnTile(epdHit->id());
    
    phi_epd_center = StraightLine_center.Phi();
    eta_epd_center = StraightLine_center.Eta();
    phi_epd_random = StraightLine_random.Phi();
    eta_epd_random = StraightLine_random.Eta();
    
    if(mip < 0.3) continue;
    TileWeight = (mip > 2) ? 2 : mip;
    
    for(int i=0; i<4; i++)     // 0: EPD-A, 1: EPD-B, 2: EPD-C, 3: EPD-D
    {
      if(i==ringgroup){
        Qx_rawep_EPD[i] += TileWeight * cos(1.0*phi_epd_center);
        Qy_rawep_EPD[i] += TileWeight * sin(1.0*phi_epd_center);
        Qx_recen_EPD[i] += TileWeight * cos(1.0*phi_epd_center) - getepd_recen[0][i]->GetBinContent(centnumber+1, runindex+1);
        Qy_recen_EPD[i] += TileWeight * sin(1.0*phi_epd_center) - getepd_recen[1][i]->GetBinContent(centnumber+1, runindex+1);
        weight_EPD[i] += TileWeight;
      }
    }
    if(ringgroup == 0 || ringgroup == 1) // 4: EPD-AB
    {
      Qx_rawep_EPD[4] += TileWeight * cos(1.0*phi_epd_center);
      Qy_rawep_EPD[4] += TileWeight * sin(1.0*phi_epd_center);
      Qx_recen_EPD[4] += TileWeight * cos(1.0*phi_epd_center) - getepd_recen[0][4]->GetBinContent(centnumber+1, runindex+1);
      Qy_recen_EPD[4] += TileWeight * sin(1.0*phi_epd_center) - getepd_recen[1][4]->GetBinContent(centnumber+1, runindex+1);
      weight_EPD[4] += TileWeight;
    }
    if(ringgroup == 2 || ringgroup == 3)  // 5: EPD-CD
    {
      Qx_rawep_EPD[5] += TileWeight * cos(1.0*phi_epd_center);
      Qy_rawep_EPD[5] += TileWeight * sin(1.0*phi_epd_center);
      Qx_recen_EPD[5] += TileWeight * cos(1.0*phi_epd_center) - getepd_recen[0][5]->GetBinContent(centnumber+1, runindex+1);
      Qy_recen_EPD[5] += TileWeight * sin(1.0*phi_epd_center) - getepd_recen[1][5]->GetBinContent(centnumber+1, runindex+1);
      weight_EPD[5] += TileWeight;
    }
    if(ringgroup == 0 || ringgroup == 1 || ringgroup == 2 || ringgroup == 3)   // 6: EPD-ABCD
    {
      Qx_rawep_EPD[6] += TileWeight * cos(1.0*phi_epd_center);
      Qy_rawep_EPD[6] += TileWeight * sin(1.0*phi_epd_center);
      Qx_recen_EPD[6] += TileWeight * cos(1.0*phi_epd_center) - getepd_recen[0][6]->GetBinContent(centnumber+1, runindex+1);
      Qy_recen_EPD[6] += TileWeight * sin(1.0*phi_epd_center) - getepd_recen[1][6]->GetBinContent(centnumber+1, runindex+1);
      weight_EPD[6] += TileWeight;
    }
  }
  
  for(int i=0; i<7; i++){
    if(Qx_rawep_EPD[i] ==0. || Qy_rawep_EPD[i] == 0. || weight_EPD[i] == 0.) return 0;
    if(Qx_recen_EPD[i] ==0. || Qy_recen_EPD[i] == 0.) return 0;
    if(weight_EPD[i] == 0) return 0;
    Qx_rawep_EPD[i] /= weight_EPD[i];
    Qy_rawep_EPD[i] /= weight_EPD[i];
    Qx_recen_EPD[i] /= weight_EPD[i];
    Qy_recen_EPD[i] /= weight_EPD[i];
  }
  
  for(int i=0; i<7; i++){
    psi1rawep_EPD[i] = atan2(Qy_rawep_EPD[i], Qx_rawep_EPD[i])/1.0;
    psi1recen_EPD[i] = atan2(Qy_recen_EPD[i], Qx_recen_EPD[i])/1.0;
    psi1shift_EPD[i] = psi1recen_EPD[i];
  }
  for(int j=0; j<7; j++)
  {
    for(int i=0; i<20; i++)
    {
      psi1shift_EPD[j] += (2.0/(i+1)) * ( -pp_EPDshiftpar_sin[j]->GetBinContent(centnumber+1, i+1, runindex+1)*cos(1.0*(i+1)*psi1recen_EPD[j]) + pp_EPDshiftpar_cos[j]->GetBinContent(centnumber+1, i+1, runindex+1)*sin(1.0*(i+1)*psi1recen_EPD[j]) );
    }
  }
  for(int i=0; i<7; i++)
  {
    if(psi1rawep_EPD[i] < -pi) {psi1rawep_EPD[i] += 2*pi;}
    if(psi1rawep_EPD[i] >  pi) {psi1rawep_EPD[i] -= 2*pi;}
    if(psi1recen_EPD[i] < -pi) {psi1recen_EPD[i] += 2*pi;}
    if(psi1recen_EPD[i] >  pi) {psi1recen_EPD[i] -= 2*pi;}
    if(psi1shift_EPD[i] < -pi) {psi1shift_EPD[i] += 2*pi;}
    if(psi1shift_EPD[i] >  pi) {psi1shift_EPD[i] -= 2*pi;}
  }
  
  psi_1_EPD_0 = psi1shift_EPD[0];
  psi_1_EPD_1 = psi1shift_EPD[1];
  psi_1_EPD_2 = psi1shift_EPD[2];
  psi_1_EPD_3 = psi1shift_EPD[3];
  psi_1_EPD_4 = psi1shift_EPD[4];
  psi_1_EPD_5 = psi1shift_EPD[5];
  psi_1_EPD_6 = psi1shift_EPD[6];
  psi_1_EPD   = psi1shift_EPD[5];
  
  psi_tree -> Fill();
  
  // EPD EP
  for(int i=0; i<7; i++)
  {
    h_EPDrawep[centnumber][i]->Fill(psi1rawep_EPD[i]);
    h_EPDrawep[9][i]         ->Fill(psi1rawep_EPD[i]);
    h_EPDrecen[centnumber][i]->Fill(psi1recen_EPD[i]);
    h_EPDrecen[9][i]         ->Fill(psi1recen_EPD[i]);
    h_EPDshift[centnumber][i]->Fill(psi1shift_EPD[i]);
    h_EPDshift[9][i]         ->Fill(psi1shift_EPD[i]);
    
    h_TPC_A_EPD[centnumber][i]->Fill(psi1shift_A_TPC, psi1shift_EPD[i]);
    h_TPC_A_EPD[9][i]         ->Fill(psi1shift_A_TPC, psi1shift_EPD[i]);
    h_TPC_B_EPD[centnumber][i]->Fill(psi1shift_B_TPC, psi1shift_EPD[i]);
    h_TPC_B_EPD[9][i]         ->Fill(psi1shift_B_TPC, psi1shift_EPD[i]);
  }
  h_EPDshift_c[centnumber][0]->Fill(psi1shift_EPD[0], psi1shift_EPD[1]);
  h_EPDshift_c[9][0]         ->Fill(psi1shift_EPD[0], psi1shift_EPD[1]);
  h_EPDshift_c[centnumber][1]->Fill(psi1shift_EPD[0], psi1shift_EPD[2]);
  h_EPDshift_c[9][1]         ->Fill(psi1shift_EPD[0], psi1shift_EPD[2]);
  h_EPDshift_c[centnumber][2]->Fill(psi1shift_EPD[0], psi1shift_EPD[3]);
  h_EPDshift_c[9][2]         ->Fill(psi1shift_EPD[0], psi1shift_EPD[3]);
  h_EPDshift_c[centnumber][3]->Fill(psi1shift_EPD[1], psi1shift_EPD[2]);
  h_EPDshift_c[9][3]         ->Fill(psi1shift_EPD[1], psi1shift_EPD[2]);
  h_EPDshift_c[centnumber][4]->Fill(psi1shift_EPD[1], psi1shift_EPD[3]);
  h_EPDshift_c[9][4]         ->Fill(psi1shift_EPD[1], psi1shift_EPD[3]);
  h_EPDshift_c[centnumber][5]->Fill(psi1shift_EPD[2], psi1shift_EPD[3]);
  h_EPDshift_c[9][5]         ->Fill(psi1shift_EPD[2], psi1shift_EPD[3]);
  
  h_EPDshift_c[centnumber][6]->Fill(psi1shift_EPD[4], psi1shift_EPD[5]);
  h_EPDshift_c[9][6]         ->Fill(psi1shift_EPD[4], psi1shift_EPD[5]);
  
  //resolution calculation using 3 sub-event method
  
  
  // cosin term
  double r1_tA_tB = cos(1.0*(psi1shift_A_TPC - psi1shift_B_TPC));
  
  // TPC 2 and EPD 6
  double r1_tA_eA = cos(1.0*(psi1shift_A_TPC - psi1shift_EPD[0]));
  double r1_tA_eB = cos(1.0*(psi1shift_A_TPC - psi1shift_EPD[1]));
  double r1_tA_eC = cos(1.0*(psi1shift_A_TPC - psi1shift_EPD[2]));
  double r1_tA_eD = cos(1.0*(psi1shift_A_TPC - psi1shift_EPD[3]));
  double r1_tA_eAB = cos(1.0*(psi1shift_A_TPC - psi1shift_EPD[4]));
  double r1_tA_eCD = cos(1.0*(psi1shift_A_TPC - psi1shift_EPD[5]));
  double r1_tA_eABCD = cos(1.0*(psi1shift_A_TPC - psi1shift_EPD[6]));
  
  double r1_tB_eA = cos(1.0*(psi1shift_B_TPC - psi1shift_EPD[0]));
  double r1_tB_eB = cos(1.0*(psi1shift_B_TPC - psi1shift_EPD[1]));
  double r1_tB_eC = cos(1.0*(psi1shift_B_TPC - psi1shift_EPD[2]));
  double r1_tB_eD = cos(1.0*(psi1shift_B_TPC - psi1shift_EPD[3]));
  double r1_tB_eAB = cos(1.0*(psi1shift_B_TPC - psi1shift_EPD[4]));
  double r1_tB_eCD = cos(1.0*(psi1shift_B_TPC - psi1shift_EPD[5]));
  double r1_tB_eABCD = cos(1.0*(psi1shift_B_TPC - psi1shift_EPD[6]));
  
  double r1_eA_eB = cos(1.0*(psi1shift_EPD[0] - psi1shift_EPD[1]));
  double r1_eA_eC = cos(1.0*(psi1shift_EPD[0] - psi1shift_EPD[2]));
  double r1_eA_eD = cos(1.0*(psi1shift_EPD[0] - psi1shift_EPD[3]));
  double r1_eB_eC = cos(1.0*(psi1shift_EPD[1] - psi1shift_EPD[2]));
  double r1_eB_eD = cos(1.0*(psi1shift_EPD[1] - psi1shift_EPD[3]));
  double r1_eC_eD = cos(1.0*(psi1shift_EPD[2] - psi1shift_EPD[3]));
  double r1_eAB_eCD = cos(1.0*(psi1shift_EPD[4] - psi1shift_EPD[5]));
  
  double r1_eA_eCD = cos(1.0*(psi1shift_EPD[0] - psi1shift_EPD[5]));
  double r1_eB_eCD = cos(1.0*(psi1shift_EPD[1] - psi1shift_EPD[5]));
  double r1_eC_eAB = cos(1.0*(psi1shift_EPD[2] - psi1shift_EPD[4]));
  double r1_eD_eAB = cos(1.0*(psi1shift_EPD[3] - psi1shift_EPD[4]));
  
  p_r1_tA_tB->Fill(centnumber, r1_tA_tB);
  
  p_r1_tA_eA->Fill(centnumber, r1_tA_eA);
  p_r1_tA_eB->Fill(centnumber, r1_tA_eB);
  p_r1_tA_eC->Fill(centnumber, r1_tA_eC);
  p_r1_tA_eD->Fill(centnumber, r1_tA_eD);
  p_r1_tA_eAB->Fill(centnumber, r1_tA_eAB);
  p_r1_tA_eCD->Fill(centnumber, r1_tA_eCD);
  p_r1_tA_eABCD->Fill(centnumber, r1_tA_eABCD);
  
  p_r1_tB_eA->Fill(centnumber, r1_tB_eA);
  p_r1_tB_eB->Fill(centnumber, r1_tB_eB);
  p_r1_tB_eC->Fill(centnumber, r1_tB_eC);
  p_r1_tB_eD->Fill(centnumber, r1_tB_eD);
  p_r1_tB_eAB->Fill(centnumber, r1_tB_eAB);
  p_r1_tB_eCD->Fill(centnumber, r1_tB_eCD);
  p_r1_tB_eABCD->Fill(centnumber, r1_tB_eABCD);
  
  p_r1_eA_eB->Fill(centnumber, r1_eA_eB);
  p_r1_eA_eC->Fill(centnumber, r1_eA_eC);
  p_r1_eA_eD->Fill(centnumber, r1_eA_eD);
  p_r1_eB_eC->Fill(centnumber, r1_eB_eC);
  p_r1_eB_eD->Fill(centnumber, r1_eB_eD);
  p_r1_eC_eD->Fill(centnumber, r1_eC_eD);
  p_r1_eAB_eCD->Fill(centnumber, r1_eAB_eCD);
  
  p_r1_eA_eCD->Fill(centnumber, r1_eA_eCD);
  p_r1_eB_eCD->Fill(centnumber, r1_eB_eCD);
  p_r1_eC_eAB->Fill(centnumber, r1_eC_eAB);
  p_r1_eD_eAB->Fill(centnumber, r1_eD_eAB);
  //cos term end
  
  // sin term
  double r1_tA_tB_sin = sin(1.0*(psi1shift_A_TPC - psi1shift_B_TPC));
  
  // TPC 2 and EPD 6
  double r1_tA_eA_sin = sin(1.0*(psi1shift_A_TPC - psi1shift_EPD[0]));
  double r1_tA_eB_sin = sin(1.0*(psi1shift_A_TPC - psi1shift_EPD[1]));
  double r1_tA_eC_sin = sin(1.0*(psi1shift_A_TPC - psi1shift_EPD[2]));
  double r1_tA_eD_sin = sin(1.0*(psi1shift_A_TPC - psi1shift_EPD[3]));
  double r1_tA_eAB_sin = sin(1.0*(psi1shift_A_TPC - psi1shift_EPD[4]));
  double r1_tA_eCD_sin = sin(1.0*(psi1shift_A_TPC - psi1shift_EPD[5]));
  double r1_tA_eABCD_sin = sin(1.0*(psi1shift_A_TPC - psi1shift_EPD[6]));
  double r1_tB_eA_sin = sin(1.0*(psi1shift_B_TPC - psi1shift_EPD[0]));
  double r1_tB_eB_sin = sin(1.0*(psi1shift_B_TPC - psi1shift_EPD[1]));
  double r1_tB_eC_sin = sin(1.0*(psi1shift_B_TPC - psi1shift_EPD[2]));
  double r1_tB_eD_sin = sin(1.0*(psi1shift_B_TPC - psi1shift_EPD[3]));
  double r1_tB_eAB_sin = sin(1.0*(psi1shift_B_TPC - psi1shift_EPD[4]));
  double r1_tB_eCD_sin = sin(1.0*(psi1shift_B_TPC - psi1shift_EPD[5]));
  double r1_tB_eABCD_sin = sin(1.0*(psi1shift_B_TPC - psi1shift_EPD[6]));
  
  double r1_eA_eB_sin = sin(1.0*(psi1shift_EPD[0] - psi1shift_EPD[1]));
  double r1_eA_eC_sin = sin(1.0*(psi1shift_EPD[0] - psi1shift_EPD[2]));
  double r1_eA_eD_sin = sin(1.0*(psi1shift_EPD[0] - psi1shift_EPD[3]));
  double r1_eB_eC_sin = sin(1.0*(psi1shift_EPD[1] - psi1shift_EPD[2]));
  double r1_eB_eD_sin = sin(1.0*(psi1shift_EPD[1] - psi1shift_EPD[3]));
  double r1_eC_eD_sin = sin(1.0*(psi1shift_EPD[2] - psi1shift_EPD[3]));
  double r1_eAB_eCD_sin = sin(1.0*(psi1shift_EPD[4] - psi1shift_EPD[5]));
  
  p_r1_tA_tB_sin->Fill(centnumber, r1_tA_tB_sin);
  
  p_r1_tA_eA_sin->Fill(centnumber, r1_tA_eA_sin);
  p_r1_tA_eB_sin->Fill(centnumber, r1_tA_eB_sin);
  p_r1_tA_eC_sin->Fill(centnumber, r1_tA_eC_sin);
  p_r1_tA_eD_sin->Fill(centnumber, r1_tA_eD_sin);
  p_r1_tA_eAB_sin->Fill(centnumber, r1_tA_eAB_sin);
  p_r1_tA_eCD_sin->Fill(centnumber, r1_tA_eCD_sin);
  p_r1_tA_eABCD_sin->Fill(centnumber, r1_tA_eABCD_sin);
  
  p_r1_tB_eA_sin->Fill(centnumber, r1_tB_eA_sin);
  p_r1_tB_eB_sin->Fill(centnumber, r1_tB_eB_sin);
  p_r1_tB_eC_sin->Fill(centnumber, r1_tB_eC_sin);
  p_r1_tB_eD_sin->Fill(centnumber, r1_tB_eD_sin);
  p_r1_tB_eAB_sin->Fill(centnumber, r1_tB_eAB_sin);
  p_r1_tB_eCD_sin->Fill(centnumber, r1_tB_eCD_sin);
  p_r1_tB_eABCD_sin->Fill(centnumber, r1_tB_eABCD_sin);
  
  p_r1_eA_eB_sin->Fill(centnumber, r1_eA_eB_sin);
  p_r1_eA_eC_sin->Fill(centnumber, r1_eA_eC_sin);
  p_r1_eA_eD_sin->Fill(centnumber, r1_eA_eD_sin);
  p_r1_eB_eC_sin->Fill(centnumber, r1_eB_eC_sin);
  p_r1_eB_eD_sin->Fill(centnumber, r1_eB_eD_sin);
  p_r1_eC_eD_sin->Fill(centnumber, r1_eC_eD_sin);
  p_r1_eAB_eCD_sin->Fill(centnumber, r1_eAB_eCD_sin);
  // sin term end
  
  
  // TPC
  // Fill Q vector
  h_TPCrawep_AQxy->Fill(Qx_rawep_TPC_A, Qy_rawep_TPC_A);
  h_TPCrawep_BQxy->Fill(Qx_rawep_TPC_B, Qy_rawep_TPC_B);
  
  h_TPCrecen_AQxy->Fill(Qx_recen_TPC_A, Qy_recen_TPC_A);
  h_TPCrecen_BQxy->Fill(Qx_recen_TPC_B, Qy_recen_TPC_B);
  
  // TPC EP
  h_TPCrawep_A[centnumber]->Fill(psi1rawep_A_TPC);
  h_TPCrawep_B[centnumber]->Fill(psi1rawep_B_TPC);
  h_TPCrawep_A[9]         ->Fill(psi1rawep_A_TPC);
  h_TPCrawep_B[9]         ->Fill(psi1rawep_B_TPC);
  h_TPCrecen_A[centnumber]->Fill(psi1recen_A_TPC);
  h_TPCrecen_B[centnumber]->Fill(psi1recen_B_TPC);
  h_TPCrecen_A[9]         ->Fill(psi1recen_A_TPC);
  h_TPCrecen_B[9]         ->Fill(psi1recen_B_TPC);
  h_TPCshift_A[centnumber]->Fill(psi1shift_A_TPC);
  h_TPCshift_B[centnumber]->Fill(psi1shift_B_TPC);
  h_TPCshift_A[9]         ->Fill(psi1shift_A_TPC);
  h_TPCshift_B[9]         ->Fill(psi1shift_B_TPC);
  h_TPCrawep_AB[centnumber]->Fill(psi1rawep_A_TPC, psi1rawep_B_TPC);
  h_TPCrawep_AB[9]         ->Fill(psi1rawep_A_TPC, psi1rawep_B_TPC);
  h_TPCrecen_AB[centnumber]->Fill(psi1recen_A_TPC, psi1recen_B_TPC);
  h_TPCrecen_AB[9]         ->Fill(psi1recen_A_TPC, psi1recen_B_TPC);
  h_TPCshift_AB[centnumber]->Fill(psi1shift_A_TPC, psi1shift_B_TPC);
  h_TPCshift_AB[9]         ->Fill(psi1shift_A_TPC, psi1shift_B_TPC);
  
  
  return kStOK;
}
//end loop

//__________________________________________________________________________________
bool Shift::isGoodEvent(const StPicoEvent *event)
{
  Float_t vx=event->primaryVertex().X();
  Float_t vy=event->primaryVertex().Y();
  Float_t vz=event->primaryVertex().Z();
  
  if(vz < 198 || vz > 202) return false;
  if(sqrt( (vx*vx) + (vy+2.0)*(vy+2.0) ) >2.0 ) return false;
  
  return true;
}

bool Shift::isGoodTrack(const StPicoTrack *ptrk) {
  const Float_t pt  = ptrk->pMom().Perp(); // zero for global tracks
  const Float_t mom = ptrk->pMom().Mag();
  const Float_t eta = ptrk->pMom().PseudoRapidity();
  const Int_t nHits = ptrk->nHits(); //TPCHits?
  const Float_t dca = ptrk->gDCA( picoEvent->primaryVertex() ).Mag();
  const Int_t nHitsFit = ptrk->nHitsFit();
  const Int_t nHitsPoss = ptrk->nHitsMax();
  const Float_t nHitsDedx = ptrk->nHitsDedx();
  const Float_t quality = (Float_t)nHitsFit/(Float_t)nHitsPoss;
  
  //if( pt < 0.06 )  return false;
  //if( mom < 0 )   return false;
  //if( eta < -2.5 || 0 < eta) return false;
  if( fabs(dca)>3.0 ) return false;
  //if( nHits < 10 )  return false;
  if( nHitsFit < 15 )  return false;
  if( quality < 0.52 )  return false;
  //if(nHitsDedx < 0 ) return false;
  
  return true;
}

bool Shift::isGoodTrackSooraj(const StPicoTrack *ptrk) {
  const Float_t pt  = ptrk->pMom().Perp(); // zero for global tracks
  const Float_t mom = ptrk->pMom().Mag();
  const Float_t eta = ptrk->pMom().PseudoRapidity();
  const Int_t nHits = ptrk->nHits(); //TPCHits?
  const Float_t dca = ptrk->gDCA( picoEvent->primaryVertex() ).Mag();
  const Int_t nHitsFit = ptrk->nHitsFit();
  const Int_t nHitsPoss = ptrk->nHitsMax();
  const Float_t nHitsDedx = ptrk->nHitsDedx();
  const Float_t quality = (Float_t)nHitsFit/(Float_t)nHitsPoss;
  
  if( fabs(dca)>3.0 ) return false;
  if( nHitsFit < 10 )  return false;
  if( quality < 0.52 )  return false;
  return true;
}
bool Shift::isGoodTrigger(const StPicoEvent* event)
{
  if(event->isTrigger(620052)) return true;
}

//__________________________________________________________________________________
