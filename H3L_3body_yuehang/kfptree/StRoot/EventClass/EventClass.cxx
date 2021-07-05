#include "EventClass.h"
#include "TTree.h"

//ClassImp(EventClass);

EventClass::EventClass(){
  
  InitValue();
}

EventClass::EventClass( TTree *t ) : ptree(t) {
  
  InitValue();
}

EventClass::~EventClass(){
}

void EventClass::InitValue(){
  
  runid   = -1;
  centid  = -1;
  refMult = 0;
  vz      = -9999.0;
  vx      = -9999.0;
  vy      = -9999.0;
  zdcRate = -9999.0;
  trgEff  = 1.0;
  
  refmultcor  = 0;
  ctofmult    = 0;
  countrefmult    = 0;
  isPileUp    = 0;
  fCentrality = -1;
  gweight     = 1.0;
  
  //	zdcAdcE = -1.0;
  //	zdcAdcW = -1.0;
  
  //	for( int k=0; k<3; k++ ) ach [k] = -9999.0;
  //	for( int k=0; k<2; k++ ) achp[k] = -9999.0;
  
  psi_1_EPD   = -9999.0;
  psi_1_EPD_0 = -9999.0;
  psi_1_EPD_1 = -9999.0;
  psi_1_EPD_2 = -9999.0;
  psi_1_EPD_3 = -9999.0;
  psi_1_EPD_4 = -9999.0;
  psi_1_EPD_5 = -9999.0;
  psi_1_EPD_6 = -9999.0;
  //	for( int isub=0; isub<6; isub++ ){
  //		tpcPsi2[isub] = -9999.0;
  //		tpcPsi3[isub] = -9999.0;
  //		tpcPsi4[isub] = -9999.0;
  //	}
  //	for( int isub=0; isub<3; isub++ ){
  //		bbcPsi1[isub] = -9999.0;
  //		bbcPsi2[isub] = -9999.0;
  //		zdcPsi [isub] = -9999.0;
  //	}
  //	for( int ih=0; ih<8; ih++ ){
  //		tpcQxB[ih] = -9999.0;
  //		tpcQyB[ih] = -9999.0;
  //		tpcQxF[ih] = -9999.0;
  //		tpcQyF[ih] = -9999.0;
  //	}
  //	for( int iw=0; iw<4; iw++ ){
  //		tpcQwB[iw] = -9999.0;
  //		tpcQwF[iw] = -9999.0;
  //	}
  
  ntrk    = 0;
  crNtrk  = 0;
  for( unsigned int itrk=0; itrk<ntrk_max; itrk++ ){
    pdg    [itrk] = -9999;
    pt     [itrk] = -9999.0;
    phi    [itrk] = -9999.0;
    eta    [itrk] = -9999.0;
    rap    [itrk] = -9999.0;
    mass   [itrk] = -9999.0;
    
    pdg_A      [itrk] = -9999;
    pdg_B      [itrk] = -9999;
    pdg_AA     [itrk] = -9999;
    pdg_AB     [itrk] = -9999;
    SIMD_pdg   [itrk] = -9999;
    SIMD_A_pdg [itrk] = -9999;
    
    dL     [itrk] = -9999.0;
    dca    [itrk] = -9999.0;
    dgdca  [itrk] = -9999.0;
    //		dcaA   [itrk] = -9999.0;
    //		dcaB   [itrk] = -9999.0;
    SIMD_ldl         [itrk] = -9999;
    SIMD_l           [itrk] = -9999;
    SIMD_dl          [itrk] = -9999;
    SIMD_chi2topo    [itrk] = -9999;
    SIMD_chi2ndf     [itrk] = -9999;
    SIMD_dca         [itrk] = -9999;
    SIMD_decaylength [itrk] = -9999;
    SIMD_lifetime    [itrk] = -9999;
    SIMD_mass        [itrk] = -9999;
    SIMD_pt          [itrk] = -9999;
    SIMD_phi         [itrk] = -9999;
    SIMD_eta         [itrk] = -9999;
    SIMD_rapidity    [itrk] = -9999;
    SIMD_dca2D       [itrk] = -9999;
    SIMD_pathlength  [itrk] = -9999;
    
    phDBRF [itrk] = -9999.0;
    thDBRF [itrk] = -9999.0;
    
    pt_db     [itrk] = -9999.0;
    phi_db    [itrk] = -9999.0;
    eta_db    [itrk] = -9999.0;
    mass_db   [itrk] = -9999.0;
    
    SIMD_A_ldl         [itrk] = -9999;
    SIMD_A_l           [itrk] = -9999;
    SIMD_A_dl          [itrk] = -9999;
    SIMD_A_chi2topo    [itrk] = -9999;
    SIMD_A_chi2ndf     [itrk] = -9999;
    SIMD_A_dca         [itrk] = -9999;
    SIMD_A_decaylength [itrk] = -9999;
    SIMD_A_lifetime    [itrk] = -9999;
    SIMD_A_mass        [itrk] = -9999;
    SIMD_A_pt          [itrk] = -9999;
    SIMD_A_phi         [itrk] = -9999;
    SIMD_A_eta         [itrk] = -9999;
    SIMD_A_rapidity    [itrk] = -9999;
    SIMD_A_dca2D       [itrk] = -9999;
    SIMD_A_pathlength  [itrk] = -9999;
    
    chi2primaryA [itrk] = -9999;
    nhitsA       [itrk] = -9999;
    rapidityA    [itrk] = -9999;
    nsigmaA0  [itrk] = -9999;
    nsigmaA1  [itrk] = -9999;
    nsigmaA2  [itrk] = -9999;
    dedxA        [itrk] = -9999;
    m2A          [itrk] = -9999;
    dL_db     [itrk] = -9999.0;
    dca_db    [itrk] = -9999.0;
    dgdca_db  [itrk] = -9999.0;
    //		dcaC     [itrk] = -9999.0;
    //		sinxi    [itrk] = -9999.0;
    
    phGDBRF    [itrk] = -9999.0;
    thGDBRF    [itrk] = -9999.0;
    
    //tofid   [itrk] = -1;
    ptB     [itrk] = -9999;
    etaB    [itrk] = -9999;
    phiB    [itrk] = -9999;
    ptAA    [itrk] = -9999;
    etaAA   [itrk] = -9999;
    phiAA   [itrk] = -9999;
    ptAB    [itrk] = -9999;
    etaAB   [itrk] = -9999;
    phiAB   [itrk] = -9999;
    
    chi2primaryB [itrk] = -9999;
    dcaB         [itrk] = -9999;
    nhitsB       [itrk] = -9999;
    rapidityB    [itrk] = -9999;
    nsigmaB0  [itrk] = -9999;
    nsigmaB1  [itrk] = -9999;
    nsigmaB2  [itrk] = -9999;
    dedxB        [itrk] = -9999;
    m2B          [itrk] = -9999;
    
    chi2primaryAA [itrk] = -9999;
    dcaAA         [itrk] = -9999;
    nhitsAA       [itrk] = -9999;
    rapidityAA    [itrk] = -9999;
    nsigmaAA0  [itrk] = -9999;
    nsigmaAA1  [itrk] = -9999;
    nsigmaAA2  [itrk] = -9999;
    dedxAA        [itrk] = -9999;
    m2AA          [itrk] = -9999;
    
    chi2primaryAB [itrk] = -9999;
    dcaAB         [itrk] = -9999;
    nhitsAB       [itrk] = -9999;
    rapidityAB    [itrk] = -9999;
    nsigmaAB0  [itrk] = -9999;
    nsigmaAB1  [itrk] = -9999;
    nsigmaAB2  [itrk] = -9999;
    dedxAB        [itrk] = -9999;
    m2AB          [itrk] = -9999;
    
    acmass[itrk] = -9999.;
  }
}

void EventClass::InitWrite(){
  
  ptree->Branch( "runid",   &runid,   "runid/I"      );
  ptree->Branch( "vz",      &vz,      "vz/F"         );
  ptree->Branch( "vx",      &vx,      "vx/F"         );
  ptree->Branch( "vy",      &vy,      "vy/F"         );
  ptree->Branch( "zdcRate", &zdcRate, "zdcRate/F"    );

  ptree->Branch( "refMult", &refMult, "refMult/S"    );
  ptree->Branch( "ctofmult",    &ctofmult   ,"ctofmult/F"   );
  ptree->Branch( "countrefmult",&countrefmult   ,"countrefmult/F"   );

  ptree->Branch( "centid",  &centid,  "centid/B"     );
  ptree->Branch( "trgEff",  &trgEff,  "trgEff/F"     );
  ptree->Branch( "refmultcor",  &refmultcor ,"refmultcor/F" );

  ptree->Branch( "ntrk",    &ntrk,    "ntrk/i"       );
  
  ptree->Branch( "isPileUp",    &isPileUp   ,"isPileUp/O"   );
  ptree->Branch( "fCentrality", &fCentrality,"fCentrality/F");
  ptree->Branch( "gweight",     &gweight    ,"gweight/F"    );
  ptree->Branch( "psi_1_EPD  ", &psi_1_EPD  , "psi_1_EPD/d"  );
  ptree->Branch( "psi_1_EPD_0", &psi_1_EPD_0, "psi_1_EPD_0/d");
  ptree->Branch( "psi_1_EPD_1", &psi_1_EPD_1, "psi_1_EPD_1/d");
  ptree->Branch( "psi_1_EPD_2", &psi_1_EPD_2, "psi_1_EPD_2/d");
  ptree->Branch( "psi_1_EPD_3", &psi_1_EPD_3, "psi_1_EPD_3/d");
  ptree->Branch( "psi_1_EPD_4", &psi_1_EPD_4, "psi_1_EPD_4/d");
  ptree->Branch( "psi_1_EPD_5", &psi_1_EPD_5, "psi_1_EPD_5/d");
  ptree->Branch( "psi_1_EPD_6", &psi_1_EPD_6, "psi_1_EPD_6/d");
  //	ptree->Branch( "tpcPsi2", tpcPsi2,  "tpcPsi2[6]/F" );
  //	ptree->Branch( "tpcPsi3", tpcPsi3,  "tpcPsi3[6]/F" );
  //	ptree->Branch( "tpcPsi4", tpcPsi4,  "tpcPsi4[6]/F" );
  //	ptree->Branch( "bbcPsi1", bbcPsi1,  "bbcPsi1[3]/F" );
  //	ptree->Branch( "bbcPsi2", bbcPsi2,  "bbcPsi2[3]/F" );
  //	ptree->Branch( "zdcPsi",  zdcPsi,   "zdcPsi[3]/F"  );
  //	ptree->Branch( "tpcQxB",  tpcQxB,   "tpcQxB[8]/F"  );
  //	ptree->Branch( "tpcQyB",  tpcQyB,   "tpcQyB[8]/F"  );
  //	ptree->Branch( "tpcQxF",  tpcQxF,   "tpcQxF[8]/F"  );
  //	ptree->Branch( "tpcQyF",  tpcQyF,   "tpcQyF[8]/F"  );
  //	ptree->Branch( "tpcQwB",  tpcQwB,   "tpcQwB[4]/F"  );
  //	ptree->Branch( "tpcQwF",  tpcQwF,   "tpcQwF[4]/F"  );
  
  // parent particle
  ptree->Branch( "pdg",   pdg,   "pdg[ntrk]/I"    );
  ptree->Branch( "pt",    pt,    "pt[ntrk]/F"     );
  ptree->Branch( "phi",   phi,   "phi[ntrk]/F"    );
  ptree->Branch( "eta",   eta,   "eta[ntrk]/F"    );
  ptree->Branch( "rap",   rap,   "rap[ntrk]/F"    );
  ptree->Branch( "mass",  mass,  "mass[ntrk]/F"   );
  
  ptree->Branch( "pdg_A"     ,pdg_A,      "pdg_A[ntrk]/I"     );
  ptree->Branch( "pdg_B"     ,pdg_B,      "pdg_B[ntrk]/I"     );
  ptree->Branch( "pdg_AA"    ,pdg_AA,     "pdg_AA[ntrk]/I"    );
  ptree->Branch( "pdg_AB"    ,pdg_AB,     "pdg_AB[ntrk]/I"    );
  ptree->Branch( "SIMD_pdg"  ,SIMD_pdg,   "SIMD_pdg[ntrk]/I"  );
  ptree->Branch( "SIMD_A_pdg",SIMD_A_pdg, "SIMD_A_pdg[ntrk]/I");
  ptree->Branch( "dL",    dL,    "dL[ntrk]/F"     );
  ptree->Branch( "dca",   dca,   "dca[ntrk]/F"    );
  ptree->Branch( "dgdca", dgdca, "dgdca[ntrk]/F"  );
  //	ptree->Branch( "dcaA",  dcaA,  "dcaA[ntrk]/F"   );
  //	ptree->Branch( "dcaB",  dcaB,  "dcaB[ntrk]/F"   );
  ptree->Branch( "SIMD_ldl",          SIMD_ldl         ,"SIMD_ldl[ntrk]/F"       );
  ptree->Branch( "SIMD_l",            SIMD_l           ,"SIMD_l[ntrk]/F"         );
  ptree->Branch( "SIMD_dl",           SIMD_dl          ,"SIMD_dl[ntrk]/F"        );
  ptree->Branch( "SIMD_chi2topo",     SIMD_chi2topo    ,"SIMD_chi2topo[ntrk]/F"  );
  ptree->Branch( "SIMD_chi2ndf",      SIMD_chi2ndf     ,"SIMD_chi2ndf[ntrk]/F"   );
  ptree->Branch( "SIMD_dca",          SIMD_dca         ,"SIMD_dca[ntrk]/F"       );
  ptree->Branch( "SIMD_decaylength",  SIMD_decaylength ,"SIMD_decaylength[ntrk]/F"       );
  ptree->Branch( "SIMD_lifetime",     SIMD_lifetime    ,"SIMD_lifetime[ntrk]/F"  );
  ptree->Branch( "SIMD_mass",         SIMD_mass        ,"SIMD_mass[ntrk]/F"      );
  ptree->Branch( "SIMD_pt",           SIMD_pt          ,"SIMD_pt[ntrk]/F"        );
  ptree->Branch( "SIMD_phi",          SIMD_phi         ,"SIMD_phi[ntrk]/F"       );
  ptree->Branch( "SIMD_eta",          SIMD_eta         ,"SIMD_eta[ntrk]/F"       );
  ptree->Branch( "SIMD_rapidity",     SIMD_rapidity    ,"SIMD_rapidity[ntrk]/F"  );
  ptree->Branch( "SIMD_dca2D",        SIMD_dca2D       ,"SIMD_dca2D[ntrk]/F"     );
  ptree->Branch( "SIMD_pathlength",   SIMD_pathlength  ,"SIMD_pathlength[ntrk]/F");
  
  // daugter baryon in parent rest frame
  ptree->Branch( "phDBRF", phDBRF, "phDBRF[ntrk]/F"  );
  ptree->Branch( "thDBRF", thDBRF, "thDBRF[ntrk]/F"  );
  
  // daughter baryon
  ptree->Branch( "pt_db",    pt_db,    "pt_db[ntrk]/F"     );
  ptree->Branch( "phi_db",   phi_db,   "phi_db[ntrk]/F"    );
  ptree->Branch( "eta_db",   eta_db,   "eta_db[ntrk]/F"    );
  ptree->Branch( "mass_db",  mass_db,  "mass_db[ntrk]/F"   );
  ptree->Branch( "dL_db",    dL_db,    "dL_db[ntrk]/F"     );
  ptree->Branch( "dca_db",   dca_db,   "dca_db[ntrk]/F"    );
  ptree->Branch( "dgdca_db", dgdca_db, "dgdca_db[ntrk]/F"  );
  //ptree->Branch( "dcaC",     dcaC,     "dcaC[ntrk]/F"      );
  //ptree->Branch( "sinxi",  sinxi,  "sinxi[ntrk]/F"   );
  
  ptree->Branch( "SIMD_A_ldl",          SIMD_A_ldl         ,"SIMD_A_ldl[ntrk]/F"       );
  ptree->Branch( "SIMD_A_l",            SIMD_A_l           ,"SIMD_A_l[ntrk]/F"         );
  ptree->Branch( "SIMD_A_dl",           SIMD_A_dl          ,"SIMD_A_dl[ntrk]/F"        );
  ptree->Branch( "SIMD_A_chi2topo",     SIMD_A_chi2topo    ,"SIMD_A_chi2topo[ntrk]/F"  );
  ptree->Branch( "SIMD_A_chi2ndf",      SIMD_A_chi2ndf     ,"SIMD_A_chi2ndf[ntrk]/F"   );
  ptree->Branch( "SIMD_A_dca",          SIMD_A_dca         ,"SIMD_A_dca[ntrk]/F"       );
  ptree->Branch( "SIMD_A_decaylength",  SIMD_A_decaylength ,"SIMD_A_decaylength[ntrk]/F"       );
  ptree->Branch( "SIMD_A_lifetime",     SIMD_A_lifetime    ,"SIMD_A_lifetime[ntrk]/F"  );
  ptree->Branch( "SIMD_A_mass",         SIMD_A_mass        ,"SIMD_A_mass[ntrk]/F"      );
  ptree->Branch( "SIMD_A_pt",           SIMD_A_pt          ,"SIMD_A_pt[ntrk]/F"        );
  ptree->Branch( "SIMD_A_phi",          SIMD_A_phi         ,"SIMD_A_phi[ntrk]/F"       );
  ptree->Branch( "SIMD_A_eta",          SIMD_A_eta         ,"SIMD_A_eta[ntrk]/F"       );
  ptree->Branch( "SIMD_A_rapidity",     SIMD_A_rapidity    ,"SIMD_A_rapidity[ntrk]/F"  );
  ptree->Branch( "SIMD_A_dca2D",        SIMD_A_dca2D       ,"SIMD_A_dca2D[ntrk]/F"     );
  ptree->Branch( "SIMD_A_pathlength",   SIMD_A_pathlength  ,"SIMD_A_pathlength[ntrk]/F");
  
  ptree->Branch( "chi2primaryA",    chi2primaryA,   "chi2primaryA[ntrk]/F"   );
  ptree->Branch( "nhitsA",          nhitsA,         "nhitsA[ntrk]/F"   );
  ptree->Branch( "rapidityA",       rapidityA,      "rapidityA[ntrk]/F"   );
  ptree->Branch( "nsigmaA0",     nsigmaA0,    "nsigmaA0[ntrk]/F"   );
  ptree->Branch( "nsigmaA1",     nsigmaA1,    "nsigmaA1[ntrk]/F"   );
  ptree->Branch( "nsigmaA2",     nsigmaA2,    "nsigmaA2[ntrk]/F"   );
  ptree->Branch( "dedxA",           dedxA,          "dedxA[ntrk]/F"   );
  ptree->Branch( "m2A",             m2A,            "m2A[ntrk]/F"   );
  
  // grand daughter in daughter baryon rest frame
  ptree->Branch( "phGDBRF",   phGDBRF, "phGDBRF[ntrk]/F"    );
  ptree->Branch( "thGDBRF",   thGDBRF, "thGDBRF[ntrk]/F"    );
  // daughters
  //ptree->Branch( "tofid",  tofid, "tofid[ntrk]/S" );
  ptree->Branch( "ptB",             ptB,            "ptB[ntrk]/S"            );
  ptree->Branch( "etaB",            etaB,           "etaB[ntrk]/S"           );
  ptree->Branch( "phiB",            phiB,           "phiB[ntrk]/S"           );
  ptree->Branch( "chi2primaryB",    chi2primaryB,   "chi2primaryB[ntrk]/F"   );
  ptree->Branch( "dcaB",            dcaB,           "dcaB[ntrk]/F"   );
  ptree->Branch( "nhitsB",          nhitsB,         "nhitsB[ntrk]/F"   );
  ptree->Branch( "rapidityB",       rapidityB,      "rapidityB[ntrk]/F"   );
  ptree->Branch( "nsigmaB0",     nsigmaB0,    "nsigmaB0[ntrk]/F"   );
  ptree->Branch( "nsigmaB1",     nsigmaB1,    "nsigmaB1[ntrk]/F"   );
  ptree->Branch( "nsigmaB2",     nsigmaB2,    "nsigmaB2[ntrk]/F"   );
  ptree->Branch( "dedxB",           dedxB,          "dedxB[ntrk]/F"   );
  ptree->Branch( "m2B",             m2B,            "m2B[ntrk]/F"   );
  
  ptree->Branch( "ptAA",             ptAA,            "ptAA[ntrk]/S"            );
  ptree->Branch( "etaAA",            etaAA,           "etaAA[ntrk]/S"           );
  ptree->Branch( "phiAA",            phiAA,           "phiAA[ntrk]/S"           );
  ptree->Branch( "chi2primaryAA",    chi2primaryAA,   "chi2primaryAA[ntrk]/F"   );
  ptree->Branch( "dcaAA",            dcaAA,           "dcaAA[ntrk]/F"   );
  ptree->Branch( "nhitsAA",          nhitsAA,         "nhitsAA[ntrk]/F"   );
  ptree->Branch( "rapidityAA",       rapidityAA,      "rapidityAA[ntrk]/F"   );
  ptree->Branch( "nsigmaAA0",     nsigmaAA0,    "nsigmaAA0[ntrk]/F"   );
  ptree->Branch( "nsigmaAA1",     nsigmaAA1,    "nsigmaAA1[ntrk]/F"   );
  ptree->Branch( "nsigmaAA2",     nsigmaAA2,    "nsigmaAA2[ntrk]/F"   );
  ptree->Branch( "dedxAA",           dedxAA,          "dedxAA[ntrk]/F"   );
  ptree->Branch( "m2AA",             m2AA,            "m2AA[ntrk]/F"   );
  
  ptree->Branch( "ptAB",             ptAB,            "ptAB[ntrk]/S"            );
  ptree->Branch( "etaAB",            etaAB,           "etaAB[ntrk]/S"           );
  ptree->Branch( "phiAB",            phiAB,           "phiAB[ntrk]/S"           );
  ptree->Branch( "chi2primaryAB",    chi2primaryAB,   "chi2primaryAB[ntrk]/F"   );
  ptree->Branch( "dcaAB",            dcaAB,           "dcaAB[ntrk]/F"   );
  ptree->Branch( "nhitsAB",          nhitsAB,         "nhitsAB[ntrk]/F"   );
  ptree->Branch( "rapidityAB",       rapidityAB,      "rapidityAB[ntrk]/F"   );
  ptree->Branch( "nsigmaAB0",     nsigmaAB0,    "nsigmaAB0[ntrk]/F"   );
  ptree->Branch( "nsigmaAB1",     nsigmaAB1,    "nsigmaAB1[ntrk]/F"   );
  ptree->Branch( "nsigmaAB2",     nsigmaAB2,    "nsigmaAB2[ntrk]/F"   );
  ptree->Branch( "dedxAB",           dedxAB,          "dedxAB[ntrk]/F"   );
  ptree->Branch( "m2AB",             m2AB,            "m2AB[ntrk]/F"   );
  
  ptree->Branch( "acmass",   acmass,  "acmass[ntrk]/F"  );
}

void EventClass::InitRead(){
  
  ptree->SetBranchAddress( "runid",   &runid    );
  ptree->SetBranchAddress( "centid",  &centid   );
  ptree->SetBranchAddress( "refMult", &refMult  );
  ptree->SetBranchAddress( "vz",      &vz       );
  ptree->SetBranchAddress( "vx",      &vx       );
  ptree->SetBranchAddress( "vy",      &vy       );
  ptree->SetBranchAddress( "zdcRate", &zdcRate  );
  ptree->SetBranchAddress( "trgEff",  &trgEff   );
  
  ptree->SetBranchAddress( "refmultcor",  &refmultcor );
  ptree->SetBranchAddress( "ctofmult",    &ctofmult   );
  ptree->SetBranchAddress( "countrefmult",&countrefmult   );
  ptree->SetBranchAddress( "isPileUp",    &isPileUp   );
  ptree->SetBranchAddress( "fCentrality", &fCentrality);
  ptree->SetBranchAddress( "gweight",     &gweight    );
  //  ptree->SetBranchAddress( "ach",     ach       );
  //  ptree->SetBranchAddress( "achp",    achp      );
  //  ptree->SetBranchAddress( "zdcAdcE", &zdcAdcE  );
  //  ptree->SetBranchAddress( "zdcAdcW", &zdcAdcW  );
  ptree->SetBranchAddress( "ntrk",    &ntrk     );
  //  ptree->SetBranchAddress( "tpcPsi2", tpcPsi2   );
  //  ptree->SetBranchAddress( "tpcPsi3", tpcPsi3   );
  //  ptree->SetBranchAddress( "tpcPsi4", tpcPsi4   );
  //  ptree->SetBranchAddress( "bbcPsi1", bbcPsi1   );
  //  ptree->SetBranchAddress( "bbcPsi2", bbcPsi2   );
  //  ptree->SetBranchAddress( "zdcPsi",  zdcPsi    );
  //  ptree->SetBranchAddress( "tpcQxB",  tpcQxB    );
  //  ptree->SetBranchAddress( "tpcQyB",  tpcQyB    );
  //  ptree->SetBranchAddress( "tpcQxF",  tpcQxF    );
  //  ptree->SetBranchAddress( "tpcQyF",  tpcQyF    );
  //  ptree->SetBranchAddress( "tpcQwB",  tpcQwB   );
  //  ptree->SetBranchAddress( "tpcQwF",  tpcQwF   );
  ptree->SetBranchAddress( "psi_1_EPD_0", &psi_1_EPD_0);
  ptree->SetBranchAddress( "psi_1_EPD_1", &psi_1_EPD_1);
  ptree->SetBranchAddress( "psi_1_EPD_2", &psi_1_EPD_2);
  ptree->SetBranchAddress( "psi_1_EPD_3", &psi_1_EPD_3);
  ptree->SetBranchAddress( "psi_1_EPD_4", &psi_1_EPD_4);
  ptree->SetBranchAddress( "psi_1_EPD_5", &psi_1_EPD_5);
  ptree->SetBranchAddress( "psi_1_EPD_6", &psi_1_EPD_6);
  //
  // parent particle
  ptree->SetBranchAddress( "pdg",   pdg     );
  ptree->SetBranchAddress( "pt",    pt      );
  ptree->SetBranchAddress( "phi",   phi     );
  ptree->SetBranchAddress( "eta",   eta     );
  ptree->SetBranchAddress( "rap",   rap     );
  ptree->SetBranchAddress( "mass",  mass    );
  ptree->SetBranchAddress( "SIMD_ldl", SIMD_ldl         );
  ptree->SetBranchAddress( "SIMD_l",   SIMD_l           );
  ptree->SetBranchAddress( "SIMD_dl",  SIMD_dl          );
  ptree->SetBranchAddress( "SIMD_chi2topo", SIMD_chi2topo    );
  ptree->SetBranchAddress( "SIMD_chi2ndf",  SIMD_chi2ndf     );
  ptree->SetBranchAddress( "SIMD_dca",      SIMD_dca         );
  ptree->SetBranchAddress( "SIMD_decaylength", SIMD_decaylength );
  ptree->SetBranchAddress( "SIMD_lifetime",    SIMD_lifetime    );
  ptree->SetBranchAddress( "SIMD_mass",        SIMD_mass        );
  ptree->SetBranchAddress( "SIMD_pt",          SIMD_pt          );
  ptree->SetBranchAddress( "SIMD_phi",         SIMD_phi         );
  ptree->SetBranchAddress( "SIMD_eta",         SIMD_eta         );
  ptree->SetBranchAddress( "SIMD_rapidity",    SIMD_rapidity    );
  ptree->SetBranchAddress( "SIMD_dca2D",       SIMD_dca2D       );
  ptree->SetBranchAddress( "SIMD_pathlength",  SIMD_pathlength  );
  
  ptree->SetBranchAddress( "pdg_A"     ,pdg_A      );
  ptree->SetBranchAddress( "pdg_B"     ,pdg_B      );
  ptree->SetBranchAddress( "pdg_AA"    ,pdg_AA     );
  ptree->SetBranchAddress( "pdg_AB"    ,pdg_AB     );
  ptree->SetBranchAddress( "SIMD_pdg"  ,SIMD_pdg   );
  ptree->SetBranchAddress( "SIMD_A_pdg",SIMD_A_pdg );
  ptree->SetBranchAddress( "dL",    dL      );
  ptree->SetBranchAddress( "dca",   dca     );
  ptree->SetBranchAddress( "dgdca", dgdca   );
  //  ptree->SetBranchAddress( "dcaA",    dcaA      );
  //  ptree->SetBranchAddress( "dcaB",    dcaB      );
  
  // daugter baryon in parent rest frame
  ptree->SetBranchAddress( "phDBRF",   phDBRF     );
  ptree->SetBranchAddress( "thDBRF",   thDBRF     );
  
  // grand daughter in daughter baryon rest frame
  ptree->SetBranchAddress( "pt_db",    pt_db      );
  ptree->SetBranchAddress( "phi_db",   phi_db     );
  ptree->SetBranchAddress( "eta_db",   eta_db     );
  ptree->SetBranchAddress( "mass_db",  mass_db    );
  ptree->SetBranchAddress( "dL_db",    dL_db      );
  ptree->SetBranchAddress( "dca_db",   dca_db     );
  ptree->SetBranchAddress( "dgdca_db", dgdca_db   );
  //ptree->SetBranchAddress( "dcaC",     dcaC       );
  //ptree->SetBranchAddress( "sinxi",  sinxi    );
  ptree->SetBranchAddress( "SIMD_A_ldl", SIMD_A_ldl         );
  ptree->SetBranchAddress( "SIMD_A_l",   SIMD_A_l           );
  ptree->SetBranchAddress( "SIMD_A_dl",  SIMD_A_dl          );
  ptree->SetBranchAddress( "SIMD_A_chi2topo", SIMD_A_chi2topo    );
  ptree->SetBranchAddress( "SIMD_A_chi2ndf",  SIMD_A_chi2ndf     );
  ptree->SetBranchAddress( "SIMD_A_dca",      SIMD_A_dca         );
  ptree->SetBranchAddress( "SIMD_A_decaylength", SIMD_A_decaylength );
  ptree->SetBranchAddress( "SIMD_A_lifetime",    SIMD_A_lifetime    );
  ptree->SetBranchAddress( "SIMD_A_mass",        SIMD_A_mass        );
  ptree->SetBranchAddress( "SIMD_A_pt",          SIMD_A_pt          );
  ptree->SetBranchAddress( "SIMD_A_phi",         SIMD_A_phi         );
  ptree->SetBranchAddress( "SIMD_A_eta",         SIMD_A_eta         );
  ptree->SetBranchAddress( "SIMD_A_rapidity",    SIMD_A_rapidity    );
  ptree->SetBranchAddress( "SIMD_A_dca2D",       SIMD_A_dca2D       );
  ptree->SetBranchAddress( "SIMD_A_pathlength",  SIMD_A_pathlength  );
  
  ptree->SetBranchAddress( "chi2primaryA",    chi2primaryA   );
  ptree->SetBranchAddress( "nhitsA",          nhitsA         );
  ptree->SetBranchAddress( "rapidityA",       rapidityA      );
  ptree->SetBranchAddress( "nsigmaA0",     nsigmaA0    );
  ptree->SetBranchAddress( "nsigmaA1",     nsigmaA1    );
  ptree->SetBranchAddress( "nsigmaA2",     nsigmaA2    );
  ptree->SetBranchAddress( "dedxA",           dedxA          );
  ptree->SetBranchAddress( "m2A",             m2A            );
  
  // grand daughter in daughter baryon rest frame
  ptree->SetBranchAddress( "phGDBRF",   phGDBRF     );
  ptree->SetBranchAddress( "thGDBRF",   thGDBRF     );
  
  // daughters
  //ptree->SetBranchAddress( "tofid",   tofid     );
  ptree->SetBranchAddress( "ptB",             ptB            );
  ptree->SetBranchAddress( "etaB",            etaB           );
  ptree->SetBranchAddress( "phiB",            phiB           );
  ptree->SetBranchAddress( "chi2primaryB",    chi2primaryB   );
  ptree->SetBranchAddress( "dcaB",            dcaB           );
  ptree->SetBranchAddress( "nhitsB",          nhitsB         );
  ptree->SetBranchAddress( "rapidityB",       rapidityB      );
  ptree->SetBranchAddress( "nsigmaB0",     nsigmaB0    );
  ptree->SetBranchAddress( "nsigmaB1",     nsigmaB1    );
  ptree->SetBranchAddress( "nsigmaB2",     nsigmaB2    );
  ptree->SetBranchAddress( "dedxB",           dedxB          );
  ptree->SetBranchAddress( "m2B",             m2B            );
  
  ptree->SetBranchAddress( "ptAB",             ptAB            );
  ptree->SetBranchAddress( "etaAB",            etaAB           );
  ptree->SetBranchAddress( "phiAB",            phiAB           );
  ptree->SetBranchAddress( "chi2primaryAB",    chi2primaryAB   );
  ptree->SetBranchAddress( "dcaAB",            dcaAB           );
  ptree->SetBranchAddress( "nhitsAB",          nhitsAB         );
  ptree->SetBranchAddress( "rapidityAB",       rapidityAB      );
  ptree->SetBranchAddress( "nsigmaAB0",     nsigmaAB0    );
  ptree->SetBranchAddress( "nsigmaAB1",     nsigmaAB1    );
  ptree->SetBranchAddress( "nsigmaAB2",     nsigmaAB2    );
  ptree->SetBranchAddress( "dedxAB",           dedxAB          );
  ptree->SetBranchAddress( "m2AB",             m2AB            );
  
  ptree->SetBranchAddress( "ptAA",             ptAA            );
  ptree->SetBranchAddress( "etaAA",            etaAA           );
  ptree->SetBranchAddress( "phiAA",            phiAA           );
  ptree->SetBranchAddress( "chi2primaryAA",    chi2primaryAA   );
  ptree->SetBranchAddress( "dcaAA",            dcaAA           );
  ptree->SetBranchAddress( "nhitsAA",          nhitsAA         );
  ptree->SetBranchAddress( "rapidityAA",       rapidityAA      );
  ptree->SetBranchAddress( "nsigmaAA0",     nsigmaAA0    );
  ptree->SetBranchAddress( "nsigmaAA1",     nsigmaAA1    );
  ptree->SetBranchAddress( "nsigmaAA2",     nsigmaAA2    );
  ptree->SetBranchAddress( "dedxAA",           dedxAA          );
  ptree->SetBranchAddress( "m2AA",             m2AA            );
  
  ptree->SetBranchAddress( "acmass",  acmass    );
}

int EventClass::set_crNtrk( unsigned int itrk ){
  
  if( itrk>=ntrk_max ) return -1;
  crNtrk = itrk;
  return 0;
}

