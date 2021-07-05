#ifndef _STEVENTCLASS_H_
#define _STEVENTCLASS_H_

#include <vector>
#include <iostream>
#include <cmath>

class TTree;
//class StTrackClass;

class EventClass {
  
public:
  EventClass( TTree *t );
  EventClass();
  ~EventClass();
  void InitValue();
  void InitWrite();
  void InitRead();
  
  
  // Global variables
  //// Get functions
  unsigned int get_ntrk (){ return ntrk;    }
  char  get_centid      (){ return centid;  }
  int   get_runid       (){ return runid;   }
  short get_refMult     (){ return refMult; }
  float get_vz          (){ return vz;      }
  float get_vx          (){ return vx;      }
  float get_vy          (){ return vy;      }
  float get_zdcRate     (){ return zdcRate; }
  float get_trgEff      (){ return trgEff;  }
  
  float get_refMultcor   (){ return refmultcor;  }
  float get_tofMult      (){ return ctofmult;  }
  float get_countrefmult (){ return countrefmult;  }
  bool   get_isPileUp     (){return isPileUp  ;}
  float  get_fCentrality  (){return fCentrality  ;}
  float  get_gweight      (){return gweight  ;}
  //  float get_zdcAdcE     (){ return zdcAdcE; }
  //  float get_zdcAdcW     (){ return zdcAdcW; }
  
  //  float get_ach ( int i ){ return ach [i];  }
  //  float get_achp( int i ){ return achp[i];  }
  
  double  get_psi_1_EPD  ()  {return psi_1_EPD  ;}
  double  get_psi_1_EPD_0()  {return psi_1_EPD_0;}
  double  get_psi_1_EPD_1()  {return psi_1_EPD_1;}
  double  get_psi_1_EPD_2()  {return psi_1_EPD_2;}
  double  get_psi_1_EPD_3()  {return psi_1_EPD_3;}
  double  get_psi_1_EPD_4()  {return psi_1_EPD_4;}
  double  get_psi_1_EPD_5()  {return psi_1_EPD_5;}
  double  get_psi_1_EPD_6()  {return psi_1_EPD_6;}
  
  //  float get_tpcPsi2 ( int i ){ return tpcPsi2[i]; } // [3] x [sub]
  //  float get_tpcPsi3 ( int i ){ return tpcPsi3[i]; }
  //  float get_tpcPsi4 ( int i ){ return tpcPsi4[i]; }
  //  float get_bbcPsi1 ( int i ){ return bbcPsi1[i]; } // [3]
  //  float get_bbcPsi2 ( int i ){ return bbcPsi2[i]; } // [3]
  //  float get_zdcPsi  ( int i ){ return zdcPsi [i]; } // [3]
  //  float get_tpcQxB  ( int i ){ return tpcQxB [i]; }
  //  float get_tpcQyB  ( int i ){ return tpcQyB [i]; }
  //  float get_tpcQxF  ( int i ){ return tpcQxF [i]; }
  //  float get_tpcQyF  ( int i ){ return tpcQyF [i]; }
  //  float get_tpcQwB  ( int i ){ return tpcQwB [i]; }
  //  float get_tpcQwF  ( int i ){ return tpcQwF [i]; }
  
  
  //// Set functions
  void set_ntrk    ( unsigned int val ){ ntrk    = val; }
  void set_centid  ( char  val )       { centid  = val; }
  void set_runid   ( int   val )       { runid   = val; }
  void set_refMult ( short val )       { refMult = val; }
  void set_vz      ( float val )       { vz      = val; }
  void set_vx      ( float val )       { vx      = val; }
  void set_vy      ( float val )       { vy      = val; }
  void set_zdcRate ( float val )       { zdcRate = val; }
  void set_trgEff  ( float val )       { trgEff  = val; }
  
  void set_refMultcor   (float val){ refmultcor   = val;}
  void set_tofMult      (float val){ ctofmult     = val;}
  void set_countrefmult (float val){ countrefmult     = val;}
  void set_isPileUp      (bool  val){ isPileUp     = val;}
  void set_fCentrality   (float val){ fCentrality  = val;}
  void set_gweight       (float val){ gweight      = val;}
  //  void set_zdcAdcE ( float val )       { zdcAdcE = val; }
  //  void set_zdcAdcW ( float val )       { zdcAdcW = val; }
  //  void set_ach     ( int i, float val ){ ach[i]  = val; }
  //  void set_achp    ( int i, float val ){ achp[i] = val; }
  //
  //  void set_tpcPsi2 ( int i, float val ){ tpcPsi2[i] = val; } // [3] x [sub]
  //  void set_tpcPsi3 ( int i, float val ){ tpcPsi3[i] = val; }
  //  void set_tpcPsi4 ( int i, float val ){ tpcPsi4[i] = val; }
  
  
  void  set_psi_1_EPD  (double val)  {psi_1_EPD   = val;}
  void  set_psi_1_EPD_0(double val)  {psi_1_EPD_0 = val;}
  void  set_psi_1_EPD_1(double val)  {psi_1_EPD_1 = val;}
  void  set_psi_1_EPD_2(double val)  {psi_1_EPD_2 = val;}
  void  set_psi_1_EPD_3(double val)  {psi_1_EPD_3 = val;}
  void  set_psi_1_EPD_4(double val)  {psi_1_EPD_4 = val;}
  void  set_psi_1_EPD_5(double val)  {psi_1_EPD_5 = val;}
  void  set_psi_1_EPD_6(double val)  {psi_1_EPD_6 = val;}
  
  //  void set_bbcPsi1 ( int i, float val ){ bbcPsi1[i] = val; } // [3]
  //  void set_bbcPsi2 ( int i, float val ){ bbcPsi2[i] = val; } // [3]
  //  void set_zdcPsi  ( int i, float val ){ zdcPsi [i] = val; } // [3]
  
  //  void set_tpcQxB  ( int i, float val ){ tpcQxB [i] = val; }
  //  void set_tpcQyB  ( int i, float val ){ tpcQyB [i] = val; }
  //  void set_tpcQxF  ( int i, float val ){ tpcQxF [i] = val; }
  //  void set_tpcQyF  ( int i, float val ){ tpcQyF [i] = val; }
  //  void set_tpcQwB  ( int i, float val ){ tpcQwB [i] = val; }
  //  void set_tpcQwF  ( int i, float val ){ tpcQwF [i] = val; }
  
  
  // Single track variables
  int set_crNtrk( unsigned int itrk );
  
  // set functions //
  // parent
  void set_pdg   ( int   val ){ pdg  [crNtrk] = val; }
  void set_pt    ( float val ){ pt   [crNtrk] = val; }
  void set_phi   ( float val ){ phi  [crNtrk] = val; }
  void set_eta   ( float val ){ eta  [crNtrk] = val; }
  void set_rap   ( float val ){ rap  [crNtrk] = val; }
  void set_mass  ( float val ){ mass [crNtrk] = val; }
  void set_dL    ( float val ){ dL   [crNtrk] = val; }
  void set_dca   ( float val ){ dca  [crNtrk] = val; }
  void set_dgdca ( float val ){ dgdca[crNtrk] = val; }
  //void set_dcaA  ( float val ){ dcaA [crNtrk] = val; }
  //void set_dcaB  ( float val ){ dcaB [crNtrk] = val; }
  
  void  set_SIMD_ldl     ( float val ){ SIMD_ldl         [crNtrk]= val; }
  void  set_SIMD_l       ( float val ){ SIMD_l           [crNtrk]= val; }
  void  set_SIMD_dl      ( float val ){ SIMD_dl          [crNtrk]= val; }
  void  set_SIMD_chi2topo( float val ){ SIMD_chi2topo    [crNtrk]= val; }
  void  set_SIMD_chi2ndf ( float val ){ SIMD_chi2ndf     [crNtrk]= val; }
  void  set_SIMD_dca     ( float val ){ SIMD_dca         [crNtrk]= val; }
  void  set_SIMD_dL      ( float val ){ SIMD_decaylength [crNtrk]= val; }
  void  set_SIMD_lT      ( float val ){ SIMD_lifetime    [crNtrk]= val; }
  void  set_SIMD_mass    ( float val ){ SIMD_mass        [crNtrk]= val; }
  void  set_SIMD_pt      ( float val ){ SIMD_pt          [crNtrk]= val; }
  void  set_SIMD_phi     ( float val ){ SIMD_phi         [crNtrk]= val; }
  void  set_SIMD_eta     ( float val ){ SIMD_eta         [crNtrk]= val; }
  void  set_SIMD_rap     ( float val ){ SIMD_rapidity    [crNtrk]= val; }
  void  set_SIMD_dca2D   ( float val ){ SIMD_dca2D       [crNtrk]= val; }
  void  set_SIMD_pL      ( float val ){ SIMD_pathlength  [crNtrk]= val; }
  
  // daughter baryon in parent rest frame
  void set_phDBRF    ( float val ){ phDBRF[crNtrk] = val; }
  void set_thDBRF    ( float val ){ thDBRF[crNtrk] = val; }
  
  // daughter baryon
  void set_pt_db    ( float val ){ pt_db   [crNtrk] = val; }
  void set_phi_db   ( float val ){ phi_db  [crNtrk] = val; }
  void set_eta_db   ( float val ){ eta_db  [crNtrk] = val; }
  void set_mass_db  ( float val ){ mass_db [crNtrk] = val; }
  void set_dL_db    ( float val ){ dL_db   [crNtrk] = val; }
  void set_dca_db   ( float val ){ dca_db  [crNtrk] = val; }
  void set_dgdca_db ( float val ){ dgdca_db[crNtrk] = val; }
  
  //void set_dcaC     ( float val ){ dcaC    [crNtrk] = val; }
  //void set_sinxi    ( float val ){ sinxi  [crNtrk] = val; }
  
  void  set_SIMD_A_ldl     ( float val ){ SIMD_A_ldl         [crNtrk]= val; }
  void  set_SIMD_A_l       ( float val ){ SIMD_A_l           [crNtrk]= val; }
  void  set_SIMD_A_dl      ( float val ){ SIMD_A_dl          [crNtrk]= val; }
  void  set_SIMD_A_chi2topo( float val ){ SIMD_A_chi2topo    [crNtrk]= val; }
  void  set_SIMD_A_chi2ndf ( float val ){ SIMD_A_chi2ndf     [crNtrk]= val; }
  void  set_SIMD_A_dca     ( float val ){ SIMD_A_dca         [crNtrk]= val; }
  void  set_SIMD_A_dL      ( float val ){ SIMD_A_decaylength [crNtrk]= val; }
  void  set_SIMD_A_lT      ( float val ){ SIMD_A_lifetime    [crNtrk]= val; }
  void  set_SIMD_A_mass    ( float val ){ SIMD_A_mass        [crNtrk]= val; }
  void  set_SIMD_A_pt      ( float val ){ SIMD_A_pt          [crNtrk]= val; }
  void  set_SIMD_A_phi     ( float val ){ SIMD_A_phi         [crNtrk]= val; }
  void  set_SIMD_A_eta     ( float val ){ SIMD_A_eta         [crNtrk]= val; }
  void  set_SIMD_A_rap     ( float val ){ SIMD_A_rapidity    [crNtrk]= val; }
  void  set_SIMD_A_dca2D   ( float val ){ SIMD_A_dca2D       [crNtrk]= val; }
  void  set_SIMD_A_pL      ( float val ){ SIMD_A_pathlength  [crNtrk]= val; }
  
  void set_chi2primaryA( float val ){chi2primaryA [crNtrk] = val;}
  //  void set_dcaA        ( float val ){dcaA         [crNtrk] = val;}
  void set_nhitsA      ( float val ){nhitsA       [crNtrk] = val;}
  void set_rapidityA   ( float val ){rapidityA    [crNtrk] = val;}
  void set_nsigmaA     (int i, float val ){
    if(i==0) nsigmaA0  [crNtrk] = val;
    if(i==1) nsigmaA1  [crNtrk] = val;
    if(i==2) nsigmaA2  [crNtrk] = val;
  }
  
  void set_dedxA       ( float val ){dedxA        [crNtrk] = val;}
  void set_m2A         ( float val ){m2A          [crNtrk] = val;}
  
  // grand-daughter baryon in daughter rest frame
  void set_phGDBRF    ( float val ){ phGDBRF[crNtrk] = val; }
  void set_thGDBRF    ( float val ){ thGDBRF[crNtrk] = val; }
  
  // daughters
  //void set_tofid( short val ){ tofid[crNtrk] = val; }
  void set_ptB  ( short val ){ ptB  [crNtrk] = val; }
  void set_etaB ( short val ){ etaB [crNtrk] = val; }
  void set_phiB ( short val ){ phiB [crNtrk] = val; }
  void set_ptAA ( short val ){ ptAA [crNtrk] = val; }
  void set_etaAA( short val ){ etaAA[crNtrk] = val; }
  void set_phiAA( short val ){ phiAA[crNtrk] = val; }
  void set_ptAB ( short val ){ ptAB [crNtrk] = val; }
  void set_etaAB( short val ){ etaAB[crNtrk] = val; }
  void set_phiAB( short val ){ phiAB[crNtrk] = val; }
  
  void set_chi2primaryB( float val ) {chi2primaryB[crNtrk] = val; }
  void set_dcaB        ( float val ) {dcaB        [crNtrk] = val; }
  void set_nhitsB      ( float val ) {nhitsB      [crNtrk] = val; }
  void set_rapidityB   ( float val ) {rapidityB   [crNtrk] = val; }
  void set_nsigmaB     (int i, float val ){
    if(i==0) nsigmaB0  [crNtrk] = val;
    if(i==1) nsigmaB1  [crNtrk] = val;
    if(i==2) nsigmaB2  [crNtrk] = val;
  }
  void set_dedxB       ( float val ) {dedxB       [crNtrk] = val; }
  void set_m2B         ( float val ) {m2B         [crNtrk] = val; }
  
  void set_chi2primaryAA( float val ){chi2primaryAA [crNtrk] = val;}
  void set_dcaAA        ( float val ){dcaAA         [crNtrk] = val;}
  void set_nhitsAA      ( float val ){nhitsAA       [crNtrk] = val;}
  void set_rapidityAA   ( float val ){rapidityAA    [crNtrk] = val;}
  void set_nsigmaAA     (int i, float val ){
    if(i==0) nsigmaAA0  [crNtrk] = val;
    if(i==1) nsigmaAA1  [crNtrk] = val;
    if(i==2) nsigmaAA2  [crNtrk] = val;
  }
  void set_dedxAA       ( float val ){dedxAA        [crNtrk] = val;}
  void set_m2AA         ( float val ){m2AA          [crNtrk] = val;}
  
  void set_pdgA          ( int  val ){pdg_A          [crNtrk] = val;}
  void set_pdgB          ( int  val ){pdg_B          [crNtrk] = val;}
  void set_pdgAA         ( int  val ){pdg_AA         [crNtrk] = val;}
  void set_pdgAB         ( int  val ){pdg_AB         [crNtrk] = val;}
  void set_SIMD_pdg      ( int  val ){SIMD_pdg       [crNtrk] = val;}
  void set_SIMD_A_pdg    ( int  val ){SIMD_A_pdg     [crNtrk] = val;}
  
  void set_chi2primaryAB( float val ){chi2primaryAB [crNtrk] = val;}
  void set_dcaAB        ( float val ){dcaAB         [crNtrk] = val;}
  void set_nhitsAB      ( float val ){nhitsAB       [crNtrk] = val;}
  void set_rapidityAB   ( float val ){rapidityAB    [crNtrk] = val;}
  void set_nsigmaAB     (int i, float val ){
    if(i==0) nsigmaAB0  [crNtrk] = val;
    if(i==1) nsigmaAB1  [crNtrk] = val;
    if(i==2) nsigmaAB2  [crNtrk] = val;
  }
  void set_dedxAB       ( float val ){dedxAB        [crNtrk] = val;}
  void set_m2AB         ( float val ){m2AB          [crNtrk] = val;}
  
  // get functions //
  // parent
  int   get_pdg      ( int i ){ return pdg  [i]; }
  float get_pt       ( int i ){ return pt   [i]; }
  float get_phi      ( int i ){ return phi  [i]; }
  float get_eta      ( int i ){ return eta  [i]; }
  float get_rap      ( int i ){ return rap  [i]; }
  float get_mass     ( int i ){ return mass [i]; }
  
  int get_pdgA          ( int i ){return pdg_A          [i];}
  int get_pdgB          ( int i ){return pdg_B          [i];}
  int get_pdgAA         ( int i ){return pdg_AA         [i];}
  int get_pdgAB         ( int i ){return pdg_AB         [i];}
  int get_SIMD_pdg      ( int i ){return SIMD_pdg       [i];}
  int get_SIMD_A_pdg    ( int i ){return SIMD_A_pdg     [i];}
  
  float get_dL       ( int i ){ return dL   [i]; }
  float get_dca      ( int i ){ return dca  [i]; }
  float get_dgdca    ( int i ){ return dgdca[i]; }
  //float get_dcaA     ( int i ){ return dcaA [i]; }
  //float get_dcaB     ( int i ){ return dcaB [i]; }
  
  float  get_SIMD_ldl     ( int i ){ return SIMD_ldl         [i]; }
  float  get_SIMD_l       ( int i ){ return SIMD_l           [i]; }
  float  get_SIMD_dl      ( int i ){ return SIMD_dl          [i]; }
  float  get_SIMD_chi2topo( int i ){ return SIMD_chi2topo    [i]; }
  float  get_SIMD_chi2ndf ( int i ){ return SIMD_chi2ndf     [i]; }
  float  get_SIMD_dca     ( int i ){ return SIMD_dca         [i]; }
  float  get_SIMD_dL      ( int i ){ return SIMD_decaylength [i]; }
  float  get_SIMD_lT      ( int i ){ return SIMD_lifetime    [i]; }
  float  get_SIMD_mass    ( int i ){ return SIMD_mass        [i]; }
  float  get_SIMD_pt      ( int i ){ return SIMD_pt          [i]; }
  float  get_SIMD_phi     ( int i ){ return SIMD_phi         [i]; }
  float  get_SIMD_eta     ( int i ){ return SIMD_eta         [i]; }
  float  get_SIMD_rap     ( int i ){ return SIMD_rapidity    [i]; }
  float  get_SIMD_dca2D   ( int i ){ return SIMD_dca2D       [i]; }
  float  get_SIMD_pL      ( int i ){ return SIMD_pathlength  [i]; }
  
  // v0 rest frame
  float get_phDBRF   ( int i ){ return phDBRF[i]; }
  float get_thDBRF   ( int i ){ return thDBRF[i]; }
  
  // daughter
  short get_charge   ( int i ){ return ( (pt_db[i]>0) ? +1 : -1); }
  float get_pt_db    ( int i ){ return fabs(pt_db[i]); }
  float get_phi_db   ( int i ){ return phi_db  [i]; }
  float get_eta_db   ( int i ){ return eta_db  [i]; }
  float get_mass_db  ( int i ){ return mass_db [i]; }
  float get_dL_db    ( int i ){ return dL_db   [i]; }
  float get_dca_db   ( int i ){ return dca_db  [i]; }
  float get_dgdca_db ( int i ){ return dgdca_db[i]; }
  //float get_dcaC     ( int i ){ return dcaC    [i]; }
  //float get_sinxi    ( int i ){ return sinxi  [i]; }
  float  get_SIMD_A_ldl     ( int i ){ return SIMD_A_ldl         [i]; }
  float  get_SIMD_A_l       ( int i ){ return SIMD_A_l           [i]; }
  float  get_SIMD_A_dl      ( int i ){ return SIMD_A_dl          [i]; }
  float  get_SIMD_A_chi2topo( int i ){ return SIMD_A_chi2topo    [i]; }
  float  get_SIMD_A_chi2ndf ( int i ){ return SIMD_A_chi2ndf     [i]; }
  float  get_SIMD_A_dca     ( int i ){ return SIMD_A_dca         [i]; }
  float  get_SIMD_A_dL      ( int i ){ return SIMD_A_decaylength [i]; }
  float  get_SIMD_A_lT      ( int i ){ return SIMD_A_lifetime    [i]; }
  float  get_SIMD_A_mass    ( int i ){ return SIMD_A_mass        [i]; }
  float  get_SIMD_A_pt      ( int i ){ return SIMD_A_pt          [i]; }
  float  get_SIMD_A_phi     ( int i ){ return SIMD_A_phi         [i]; }
  float  get_SIMD_A_eta     ( int i ){ return SIMD_A_eta         [i]; }
  float  get_SIMD_A_rap     ( int i ){ return SIMD_A_rapidity    [i]; }
  float  get_SIMD_A_dca2D   ( int i ){ return SIMD_A_dca2D       [i]; }
  float  get_SIMD_A_pL      ( int i ){ return SIMD_A_pathlength  [i]; }
  
  // xi rest frame
  float get_phGDBRF  ( int i ){ return phGDBRF [i]; }
  float get_thGDBRF  ( int i ){ return thGDBRF [i]; }
  
  // daughters
  float get_chi2primaryA( int i ) { return chi2primaryA[i]; }
  //  float get_dcaA        ( int i ) { return dcaA        [i]; }
  float get_nhitsA      ( int i ) { return nhitsA      [i]; }
  float get_rapidityA   ( int i ) { return rapidityA   [i]; }
  float get_nsigmaA ( int i, int j ) {
    if(i==0) return nsigmaA0[j];
    if(i==1) return nsigmaA1[j];
    if(i==2) return nsigmaA2[j];
  }
  float get_dedxA       ( int i ) { return dedxA       [i]; }
  float get_m2A         ( int i ) { return m2A         [i]; }
  
  //short get_tofid( int i ){ return tofid[i];      }
  float get_ptB  ( int i ){ return ptB[i]/1000.;  }
  float get_etaB ( int i ){ return etaB[i]/1000.; }
  float get_phiB ( int i ){ return phiB[i]/1000.; }
  float get_ptAA ( int i ){ return ptAA[i]/1000.;  }
  float get_etaAA( int i ){ return etaAA[i]/1000.; }
  float get_phiAA( int i ){ return phiAA[i]/1000.; }
  float get_ptAB ( int i ){ return ptAB[i]/1000.;  }
  float get_etaAB( int i ){ return etaAB[i]/1000.; }
  float get_phiAB( int i ){ return phiAB[i]/1000.; }
  
  float get_chi2primaryB( int i ) { return chi2primaryB[i]; }
  float get_dcaB        ( int i ) { return dcaB        [i]; }
  float get_nhitsB      ( int i ) { return nhitsB      [i]; }
  float get_rapidityB   ( int i ) { return rapidityB   [i]; }
  float get_nsigmaB ( int i, int j ) {
    if(i==0) return nsigmaB0[j];
    if(i==1) return nsigmaB1[j];
    if(i==2) return nsigmaB2[j];
  }
  float get_dedxB       ( int i ) { return dedxB       [i]; }
  float get_m2B         ( int i ) { return m2B         [i]; }
  
  float get_chi2primaryAA( int i ){ return chi2primaryAA [i];}
  float get_dcaAA        ( int i ){ return dcaAA         [i];}
  float get_nhitsAA      ( int i ){ return nhitsAA       [i];}
  float get_rapidityAA   ( int i ){ return rapidityAA    [i];}
  float get_nsigmaAA ( int i, int j ) {
    if(i==0) return nsigmaAA0[j];
    if(i==1) return nsigmaAA1[j];
    if(i==2) return nsigmaAA2[j];
  }
  float get_dedxAA       ( int i ){ return dedxAA        [i];}
  float get_m2AA         ( int i ){ return m2AA          [i];}
  
  float get_chi2primaryAB( int i ){ return chi2primaryAB [i];}
  float get_dcaAB        ( int i ){ return dcaAB         [i];}
  float get_nhitsAB      ( int i ){ return nhitsAB       [i];}
  float get_rapidityAB   ( int i ){ return rapidityAB    [i];}
  float get_nsigmaAB ( int i, int j ) {
    if(i==0) return nsigmaAB0[j];
    if(i==1) return nsigmaAB1[j];
    if(i==2) return nsigmaAB2[j];
  }
  float get_dedxAB       ( int i ){ return dedxAB        [i];}
  float get_m2AB         ( int i ){ return m2AB          [i];}
  
  void set_acmass( float val ){ acmass[crNtrk] = val;  }
  float get_acmass( int i ){ return acmass[i]; }
  
private:
  TTree *ptree;
  int   runid;      // run number
  char  centid;     // 0=75-80%,..., 15=0-5%
  short refMult;    // refMult for run11, grefMult for run14
  float vz;
  float vx;
  float vy;
  float zdcRate;    // zdc coincidence rate
  float trgEff;     // trigger weight
  
  float refmultcor ;
  float ctofmult   ;
  float countrefmult;

  bool  isPileUp   ;
  float fCentrality;
  float gweight    ;
  //		float ach [3];    // charge asymmetry (0=charge, 1=kaon, 2= proton)
  //		float achp[2];    // charge asymmetry excluidng Lambda's daughters (0=charge, 1=proton)
  //  float tpcPsi2[6]; // TPC-Psi2,  3x2 (East, West, East+West) x (|eta|>0.1, |eta|>0.5)
  //  float tpcPsi3[6]; // TPC-Psi3
  //  float tpcPsi4[6]; // TPC-Psi4
  double  psi_1_EPD  ;
  double  psi_1_EPD_0;
  double  psi_1_EPD_1;
  double  psi_1_EPD_2;
  double  psi_1_EPD_3;
  double  psi_1_EPD_4;
  double  psi_1_EPD_5;
  double  psi_1_EPD_6;
  //		float bbcPsi1[3]; // BBC-Psi1: East, West, East+West
  //		float bbcPsi2[3]; // BBC-Psi2: East, West, East+West
  //		float zdcPsi[3];  // ZDCSMD-Psi1: East, West, East+West
  //		float zdcAdcE;    // zdcAdcSum east
  //		float zdcAdcW;    // zdcAdcSum west
  
  // tpc "raw" flow vectors for -1<eta<-0.1 and 0.1<eta<1 and h+/h-
  //		float tpcQxB[8]; // backward eta, harmonics = 1, 2, 3, 4
  //		float tpcQyB[8]; // 4 harmonics x 2 charges = 8
  //		float tpcQxF[8]; // forward eta, harmonics = 1, 2, 3, 4
  //		float tpcQyF[8];
  //		float tpcQwB[4]; // 0= multiplicity, 1= sum of pT
  //		float tpcQwF[4]; // 2 weights x 2 charges = 4
  
  unsigned int ntrk;// number of v0
  unsigned int crNtrk;
  static const unsigned int ntrk_max = 1000;
  static const unsigned int kpart = 3;
  
  
  // parent info
  int   pdg  [ntrk_max]; // particle id
  float pt   [ntrk_max]; // pt
  float phi  [ntrk_max];
  float eta  [ntrk_max];
  float rap  [ntrk_max];
  float mass [ntrk_max]; // invarinat mass of parent particle
  float dL   [ntrk_max]; // decay length
  float dca  [ntrk_max]; // dca to primary vertex
  float dgdca[ntrk_max]; // dca between two daughters
  //  float dcaA   [ntrk_max]; // dca of daughter baryon from parent
  //  float dcaB   [ntrk_max]; // dca of daughter meson from parent
  //  float dcaC   [ntrk_max]; // dca of daughter meson from parent
  
  int pdg_A          [ntrk_max];
  int pdg_B          [ntrk_max];
  int pdg_AA         [ntrk_max];
  int pdg_AB         [ntrk_max];
  int SIMD_pdg       [ntrk_max];
  int SIMD_A_pdg     [ntrk_max];
  
  float SIMD_ldl         [ntrk_max];
  float SIMD_l           [ntrk_max];
  float SIMD_dl          [ntrk_max];
  float SIMD_chi2topo    [ntrk_max];
  float SIMD_chi2ndf     [ntrk_max];
  float SIMD_dca         [ntrk_max];
  float SIMD_decaylength [ntrk_max];
  float SIMD_lifetime    [ntrk_max];
  float SIMD_mass        [ntrk_max];
  float SIMD_pt          [ntrk_max];
  float SIMD_phi         [ntrk_max];
  float SIMD_eta         [ntrk_max];
  float SIMD_rapidity    [ntrk_max];
  float SIMD_dca2D       [ntrk_max];
  float SIMD_pathlength  [ntrk_max];
  
  // daughter baryon info.
  float pt_db   [ntrk_max];
  float phi_db  [ntrk_max];
  float eta_db  [ntrk_max];
  float mass_db [ntrk_max]; // (invarinat) mass of daughter baryon
  float dL_db   [ntrk_max];
  float dca_db  [ntrk_max]; // dca of daughter baryon to pvtx
  float dgdca_db[ntrk_max]; // dca between two grand daughters
  //  float dcaC    [ntrk_max]; // dca of pi from Xi
  //float sinxi  [ntrk_max];
  float SIMD_A_ldl         [ntrk_max];
  float SIMD_A_l           [ntrk_max];
  float SIMD_A_dl          [ntrk_max];
  float SIMD_A_chi2topo    [ntrk_max];
  float SIMD_A_chi2ndf     [ntrk_max];
  float SIMD_A_dca         [ntrk_max];
  float SIMD_A_decaylength [ntrk_max];
  float SIMD_A_lifetime    [ntrk_max];
  float SIMD_A_mass        [ntrk_max];
  float SIMD_A_pt          [ntrk_max];
  float SIMD_A_phi         [ntrk_max];
  float SIMD_A_eta         [ntrk_max];
  float SIMD_A_rapidity    [ntrk_max];
  float SIMD_A_dca2D       [ntrk_max];
  float SIMD_A_pathlength  [ntrk_max];
  
  float   chi2primaryA [ntrk_max];
  //  float   dcaA         [ntrk_max];
  float   nhitsA       [ntrk_max];
  float   rapidityA    [ntrk_max];
  float   nsigmaA0      [ntrk_max];
  float   nsigmaA1      [ntrk_max];
  float   nsigmaA2      [ntrk_max];
  float   dedxA        [ntrk_max];
  float   m2A          [ntrk_max];
  
  float phDBRF [ntrk_max]; // phi   of daughter baryon in parent rest frame
  float thDBRF [ntrk_max]; // theta of daughter baryon in parent rest frame
  float phGDBRF[ntrk_max]; // phi   of grand-daughter baryon in daughter baryon rest frame
  float thGDBRF[ntrk_max]; // theta of grand-daughter baryon in daughter baryon rest frame
  
  // v0 daughter info
  //short tofid[ntrk_max]; // tof matched info.
  // 1's digit for V0
  // 0=both matched, 1=both unmatched, 2=only pi matched, 3=only p matched
  // 10's digit for particle C
  // 0=unmatched, 1=matched
  short ptB  [ntrk_max]; // pt*1000 of pions or kaons with short
  short etaB [ntrk_max];
  short phiB [ntrk_max];
  short ptAA [ntrk_max]; // pt*1000 of granddaughter baryon with short
  short etaAA[ntrk_max];
  short phiAA[ntrk_max];
  short ptAB [ntrk_max]; // pt*1000 of granddaughter meson with short
  short etaAB[ntrk_max];
  short phiAB[ntrk_max];
  
  float chi2primaryB [ntrk_max];
  float dcaB         [ntrk_max];
  float nhitsB       [ntrk_max];
  float rapidityB    [ntrk_max];
  float   nsigmaB0      [ntrk_max];
  float   nsigmaB1      [ntrk_max];
  float   nsigmaB2      [ntrk_max];
  float dedxB        [ntrk_max];
  float m2B          [ntrk_max];
  
  float chi2primaryAA [ntrk_max];
  float dcaAA         [ntrk_max];
  float nhitsAA       [ntrk_max];
  float rapidityAA    [ntrk_max];
  float   nsigmaAA0      [ntrk_max];
  float   nsigmaAA1      [ntrk_max];
  float   nsigmaAA2      [ntrk_max];
  float dedxAA        [ntrk_max];
  float m2AA          [ntrk_max];
  
  float chi2primaryAB [ntrk_max];
  float dcaAB         [ntrk_max];
  float nhitsAB       [ntrk_max];
  float rapidityAB    [ntrk_max];
  float   nsigmaAB0      [ntrk_max];
  float   nsigmaAB1      [ntrk_max];
  float   nsigmaAB2      [ntrk_max];
  float dedxAB        [ntrk_max];
  float m2AB          [ntrk_max];
  
  //bool  type  [ntrk_max]; // 0= real pair,  1= rotated pair
  float acmass[ntrk_max];
  
  //ClassDef(EventClass,1);
};

#endif
