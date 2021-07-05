//v3: add nhits, dedx
//testing ncut=10 bdtv7
//THIS MACRO REQUIRES TMVA package, better run on rcf
//write mass/pt histograms that can be read offline
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TLorentzVector.h"

#include <fstream>
#include <sstream>


#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif


const int npttmvabins = 1;
const int nptbins = 6;
const int ncenttmvabins = 3;
const int ncentbins = 8;

double lambda_mass;
double cascade_mass;
double omega_mass;

double _weight;

TF1 *om_width_0;
TF1 *om_width_1;
TF1 *gitom;



 TH1F	*hvtx;
 TH1F	*hvtxgood;
 TH1F	*hrefmult;

 TH1F   *gvtx;
 TH1F   *gvtxgood;
 TH1F   *grefmult;

 double sgm3;
 double sgm4;

//TH2F *mkfIncMassPtUL3_0[ncenttmvabins];
//TH2F *mkfIncMassPtUL3_1[ncenttmvabins];
//TH2F *mkfIncMassPtUL4_0[ncenttmvabins];
//TH2F *mkfIncMassPtUL4_1[ncenttmvabins];
//
TH2F *mkfIncMassPtUL3_0[ncentbins];
TH2F *mkfIncMassPtUL3_1[ncentbins];
TH2F *mkfIncMassPtUL4_0[ncentbins];
TH2F *mkfIncMassPtUL4_1[ncentbins];


TH2F *h_phieta;
TH2F *h_phieta_bg;

double ptbin[nptbins+1] = {0.4,0.8,1.2,1.6,2.0,2.4,2.8};

const int ncutptbins = 3;
double ptcutbin[ncutptbins+1] = {0.4,1.2,2.0,2.8};
const int ncutcentbins = 3;
int centcutbin[ncutcentbins+1] = {0,10,40,80};

TH1F *h_om_chi2topo[ncutcentbins][ncutptbins];
TH1F *h_om_chi2ndf[ncutcentbins][ncutptbins];
TH1F *h_om_chi2topo_bg[ncutcentbins][ncutptbins];
TH1F *h_om_chi2ndf_bg[ncutcentbins][ncutptbins];

TH1F *h_om_ldl[ncutcentbins][ncutptbins];
TH1F *h_om_ld_chi2topo[ncutcentbins][ncutptbins];
TH1F *h_om_ld_chi2ndf[ncutcentbins][ncutptbins];
TH1F *h_om_ld_ldl[ncutcentbins][ncutptbins];
TH1F *h_om_ld_l[ncutcentbins][ncutptbins];
TH1F *h_om_l[ncutcentbins][ncutptbins];

TH1F *h_om_ldl_bg[ncutcentbins][ncutptbins];
TH1F *h_om_ld_chi2topo_bg[ncutcentbins][ncutptbins];
TH1F *h_om_ld_chi2ndf_bg[ncutcentbins][ncutptbins];
TH1F *h_om_ld_ldl_bg[ncutcentbins][ncutptbins];
TH1F *h_om_ld_l_bg[ncutcentbins][ncutptbins];
TH1F *h_om_l_bg[ncutcentbins][ncutptbins];

TH1F *h_nhits_om_proton[ncutcentbins][ncutptbins];
TH1F *h_nhits_om_pi[ncutcentbins][ncutptbins];
TH1F *h_nhits_om_bach[ncutcentbins][ncutptbins];
TH1F *h_dedx_om_proton[ncutcentbins][ncutptbins];
TH1F *h_dedx_om_pi[ncutcentbins][ncutptbins];
TH1F *h_dedx_om_bach[ncutcentbins][ncutptbins];

TH1F *h_nhits_om_proton_bg[ncutcentbins][ncutptbins];
TH1F *h_nhits_om_pi_bg[ncutcentbins][ncutptbins];
TH1F *h_nhits_om_bach_bg[ncutcentbins][ncutptbins];
TH1F *h_dedx_om_proton_bg[ncutcentbins][ncutptbins];
TH1F *h_dedx_om_pi_bg[ncutcentbins][ncutptbins];
TH1F *h_dedx_om_bach_bg[ncutcentbins][ncutptbins];

TH1F *h_chi2primary_om_proton[ncutcentbins][ncutptbins];
TH1F *h_chi2primary_om_pi[ncutcentbins][ncutptbins];
TH1F *h_chi2primary_om_bach[ncutcentbins][ncutptbins];

TH1F *h_chi2primary_om_proton_bg[ncutcentbins][ncutptbins];
TH1F *h_chi2primary_om_pi_bg[ncutcentbins][ncutptbins];
TH1F *h_chi2primary_om_bach_bg[ncutcentbins][ncutptbins];

TH1F *h_om_mva[ncutcentbins][ncutptbins];
TH1F *h_om_mva_bg[ncutcentbins][ncutptbins];

TH1F *h_om_mass[ncutcentbins][ncutptbins];

int centtmva_limit[ncenttmvabins+1] = {0,10,40,80};

// int cent_limit[ncentbins+1] = {0,5,10,20,30,40,50,60,70,80};
int cent_limit[ncentbins+1] = {0,5,10,20,30,40,50,60,80};
//
//
int snn=28;

int    _cut_mode       ;
//_cut_mode = 0;//loose
//_cut_mode = 1;//default
//_cut_mode = 2;//bdt

double mvaValue        ;
double delta           ;

int    brunid          ;
int    beventid        ;
//float  brefmult        ;
int  brefmult        ;

int    btofmult        ;
float  bVz             ;
int    bparticleid     ;
float  bparticlemass   ;
float  bpx             ;
float  bpy             ;
float  bpz             ;
float  bpt	           ;
int    _rotate         ;

int    pionidtruth     ;
int    protonidtruth   ;

int    _data_or_sim = 0;

/*
float ld_pv_dca_ld;
float p_pv_dca_ld;
float pi_pv_dca_ld;
float pi_p_dca_ld;
float decaylen_ld;
float ang0_ld;
*/

// float ld_chi2topo     ;
// float ld_chi2ndf      ;
// float ld_ldl          ;
// float ld_l            ;
// float ld_dl           ;
// float chi2primary_proton ;
// float chi2primary_pi     ;

// float Chi2NDF;
// float LdL;
// float Chi2Topo;
// float L;
// float chi2Primary_Pi;
// float chi2Primary_Proton;

int nhits_om_proton;
int nhits_om_pi;
int nhits_om_bach;
float dedx_om_proton;
float dedx_om_pi;
float dedx_om_bach;

float om_chi2topo;
float om_chi2ndf;
float om_ldl;
float om_ld_chi2topo;
float om_ld_chi2ndf;
float om_ld_ldl;
float om_ld_l;
float om_l;
float om_dl;
float chi2primary_om_proton;
float chi2primary_om_pi;
float chi2primary_om_bach;
float chi2primary_om_ld;

float om_chi2topo_cut;
float om_chi2ndf_cut;
float om_ldl_cut;
float om_ld_chi2topo_cut;
float om_ld_chi2ndf_cut;
float om_ld_ldl_cut;
float om_ld_l_cut;
float om_l_cut;
float om_dl_cut;
float chi2primary_om_proton_cut;
float chi2primary_om_pi_cut;
float chi2primary_om_bach_cut;
float chi2primary_om_ld_cut;

float chi2Primary_Kaon;
float chi2Primary_Proton;
float chi2Primary_Pion;
float Chi2NDF;
float LdL;
float L;
float Chi2Topo;
float Chi2NDF_ld;
float LdL_ld;
float L_ld;
float Chi2Topo_ld;

int cent9;

//standard cut
/*
float  ld_pv_dca_ld_cut = 0.8;
float  p_pv_dca_ld_cut = 0.5;
float  pi_pv_dca_ld_cut = 1.5;
float  pi_p_dca_ld_cut = 0.8;
float  decaylen_ld_cut = 3.0;
float  ang0_ld_cut = 0.0;
*/

float ld_chi2topo_cut = 5     ;
float ld_chi2ndf_cut = 10     ;
float ld_ldl_cut = 5          ;
float ld_l_cut = 5            ;
//float ld_dl           ;
float chi2primary_proton_cut = 18.6 ;
float chi2primary_pi_cut = 18.6     ;
//end standard cut
int    _data_or_mc;
Char_t path_ch[10000];
Char_t ofile[10000];

TMVA::Reader* reader[ncenttmvabins][npttmvabins];
// float bdt_ld_cut[ncenttmvabins][npttmvabins];
float bdt_ld_cut[ncenttmvabins][nptbins];


void book_histo();
void write_histo();
void process_tree_tmva_om_kf_varcut_v3(int data_or_mc=0){
    _data_or_mc = data_or_mc;//data_or_mc=0: mc, data_or_mc=1: data,
    // _cut_mode = cut_mode;

//cuts

om_chi2topo_cut = 5;
om_chi2ndf_cut = 6;
om_ldl_cut = 10;
om_ld_chi2topo_cut = 5;
om_ld_chi2ndf_cut = 10;
om_ld_ldl_cut = 5;
om_ld_l_cut = 5;
om_l_cut = 5;
chi2primary_om_proton_cut = 18.6;
chi2primary_om_pi_cut = 18.6;
chi2primary_om_bach_cut = 18.6;


lambda_mass = 1.115683;
cascade_mass = 1.32171;
omega_mass = 1.67242;


//end cuts

 // This loads the library
     // TMVA::Tools::Instance();

//CHANGE NAME TO MATCH READER>>>

//TMVA::Reader* reader = new TMVA::Reader( "!Color:!Silent" );

for(int icent=0;icent<ncenttmvabins;icent++){
   for(int ipt=0;ipt<npttmvabins;ipt++){

      reader[icent][ipt] =  new TMVA::Reader( "!Color:!Silent" );

      reader[icent][ipt]->AddVariable("chi2Primary_Kaon", &chi2Primary_Kaon);
      reader[icent][ipt]->AddVariable("chi2Primary_Proton", &chi2Primary_Proton);
      reader[icent][ipt]->AddVariable("chi2Primary_Pion", &chi2Primary_Pion);
      reader[icent][ipt]->AddVariable("Chi2NDF", &Chi2NDF);
      reader[icent][ipt]->AddVariable("LdL", &LdL);
      reader[icent][ipt]->AddVariable("L", &L);
      reader[icent][ipt]->AddVariable("Chi2Topo", &Chi2Topo);
      reader[icent][ipt]->AddVariable("Chi2NDF_ld", &Chi2NDF_ld);
      reader[icent][ipt]->AddVariable("LdL_ld", &LdL_ld);
      reader[icent][ipt]->AddVariable("L_ld", &L_ld);
      reader[icent][ipt]->AddVariable("Chi2Topo_ld", &Chi2Topo_ld);

//      reader[icent][ipt]->BookMVA( "BDT" , Form("/star/u/yhleung2/pwg/06.tmva/weights/tmva_om_weight_kf_cent_%d_pt_%d_BDT.weights.xml",icent,ipt));
//
//        reader[icent][ipt]->BookMVA( "BDT" , Form("/star/u/yhleung2/pwg/06.tmva/weights/tmva_om_weight_kf_cent_%d_pt_%d_v4_BDT.weights.xml",icent,ipt));      
//        reader[icent][ipt]->BookMVA( "BDT" , Form("/star/u/yhleung2/pwg/06.tmva/weights/tmva_om_weight_kf_cent_%d_pt_%d_v5_BDT.weights.xml",icent,ipt));      
//          reader[icent][ipt]->BookMVA( "BDT" , Form("/star/u/yhleung2/pwg/06.tmva/weights/tmva_om_weight_kf_cent_%d_pt_%d_v6_BDT.weights.xml",icent,ipt));
          reader[icent][ipt]->BookMVA( "BDT" , Form("/star/u/yhleung2/pwg/06.tmva/weights/tmva_om_weight_kf_cent_%d_pt_%d_v7_BDT.weights.xml",icent,ipt));

}//pt
}//cent

// bdt_ld_cut[0][0] = 0.0272;
// bdt_ld_cut[0][1] = -0.0502;
// bdt_ld_cut[0][2] = -0.0729;
// bdt_ld_cut[0][3] = -0.0729;
// bdt_ld_cut[0][4] = -0.0729;
// bdt_ld_cut[0][5] = -0.0729;
//
// bdt_ld_cut[1][0] = 0.0127;
// bdt_ld_cut[1][1] = -0.0793;
// bdt_ld_cut[1][2] = -0.1243;
// bdt_ld_cut[1][3] = -0.1243;
// bdt_ld_cut[1][4] = -0.1243;
// bdt_ld_cut[1][5] = -0.1243;
//
// bdt_ld_cut[2][0] = -0.0606;
// bdt_ld_cut[2][1] = -0.1185;
// bdt_ld_cut[2][2] = -0.1744;
// bdt_ld_cut[2][3] = -0.1744;
// bdt_ld_cut[2][4] = -0.1744;
// bdt_ld_cut[2][5] = -0.1744;


TChain lambda_tree("lambda_tree");
TChain cascade_tree("cascade_tree");
TChain omega_tree("omega_tree");
//TChain omega_tree("Om");

int nChains = -999999;
ifstream fin;

//fin.open("/star/u/yhleung2/pwg/02.long_rotate/taxi_training/filelist.txt");
//strcpy (ofile,"run18_27gev_tmva_train_test.root");
//fin.open("/star/u/yhleung2/pwg/02.long_rotate/taxi/filelist.txt");
//fin.open("/star/u/yhleung2/pwg/02.long_rotate/taxi/filelist_test.txt");
if(_data_or_mc==0){//mc
// fin.open("/star/u/yhleung2/pwg/05.KFParticle_rotate/taxi_om_s2b/filelist.txt");

//fin.open("/Users/yhleung/Remote/downloads/taxi_om_sig/filelist.txt");
//strcpy (ofile,"run18_27gev_om_kf_mc_cut_var.root");
//nChains = 1;

// fin.open("/Users/yhleung/Remote/downloads/taxi_om_sig/filelist_old.txt");
// strcpy (ofile,"run18_27gev_om_kf_mc_cut_var_old.root");
// nChains = 1;

//LOCAL
// fin.open("/Users/yhleung/Remote/downloads/taxi_om_sig/filelist_v1.txt");
// strcpy (ofile,"run18_27gev_om_kf_mc_cut_var_v1.root");
// nChains = 1;

//RCF
//fin.open("/star/u/yhleung2/pwg/05.KFParticle_rotate/taxi_om_sig_v3/sum/filelist.txt");
//strcpy (ofile,"run18_27gev_om_kf_mc_cut_var_v2.root");
//nChains = 1;
//
//ksigma3
//fin.open("/star/u/yhleung2/pwg/05.KFParticle_rotate/taxi_om_sig_vksigma3/filelist.txt");
//strcpy (ofile,"run18_27gev_om_kf_mc_cut_var_vksigma3.root");
//nChains = 223;
//
//ksigma3
//fin.open("/star/u/yhleung2/pwg/05.KFParticle_rotate/taxi_om_sig_vksigma3/filelist.txt");
//strcpy (ofile,"run18_27gev_om_kf_mc_cut_var_vksigma3_bdtv4.root");
//strcpy (ofile,"run18_27gev_om_kf_mc_cut_var_vksigma3_bdtv5.root");
//strcpy (ofile,"run18_27gev_om_kf_mc_cut_var_vksigma3_bdtv6.root");
//strcpy (ofile,"run18_27gev_om_kf_mc_cut_var_vksigma3_bdtv7.root");
//nChains = 223;
//
//nhits sample
fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_om_sig_vksigma3/filelist.txt");
strcpy (ofile,"run18_27gev_om_kf_mc_cut_var_vksigma3_nhits.root");
//nChains = 3;
nChains = 391;

}
if(_data_or_mc==1){//data
// fin.open("/Users/yhleung/Remote/downloads/taxi_om_bg/filelist.txt");
// strcpy (ofile,"run18_27gev_om_kf_data_cut_var.root");
// nChains = 2;



//fin.open("/star/u/yhleung2/pwg/05.KFParticle_rotate/taxi_om_bg_v3/sum/filelist.txt");
//strcpy (ofile,"run18_27gev_om_kf_data_cut_var_v2.root");
//nChains = 16;

//fin.open("/star/u/yhleung2/pwg/05.KFParticle_rotate/taxi_ob_bg_v3/sum/filelist.txt");
//strcpy (ofile,"run18_27gev_ob_kf_data_cut_var_v2.root");
//nChains = 2;
//
//
//ksigma3 smmall sample
//fin.open("/star/u/yhleung2/pwg/05.KFParticle_rotate/taxi_om_bg_vksigma3/filelist.txt");
//strcpy (ofile,"run18_27gev_om_kf_data_var_vksigma3_smallsample.root");//zaochens centrality selection
//nChains = 2;

//ksigma3 (close to full) sample
//fin.open("/star/u/yhleung2/pwg/05.KFParticle_rotate/taxi_om_bg_vksigma3_v2/filelist.txt");
//strcpy (ofile,"run18_27gev_om_kf_data_var_vksigma3_bdtv7.root");
//nChains = 7027;


//small sample for nhits only
fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_om_data_ana/filelist.txt");
strcpy (ofile,"run18_27gev_om_kf_data_var_vksigma3_nhits.root");
nChains = 1815;
//nChains = 2;

}

for(int i=0; i<nChains; i++){
      fin >> path_ch;

      omega_tree.Add(path_ch);
      cout<<path_ch<<"  "<<omega_tree.GetEntries() <<endl;


//      TFile *f0 = new TFile(path_ch,"READ");

//       gvtx = (TH1F*)f0->Get("hvtx");
//       gvtxgood = (TH1F*)f0->Get("hvtxgood");
//       grefmult = (TH1F*)f0->Get("hrefmult");
//
// if(i==0){
//  hvtx = (TH1F*)gvtx->Clone("hvtx");
//  hvtxgood = (TH1F*)gvtxgood->Clone("hvtxgood");
//  hrefmult = (TH1F*)grefmult->Clone("hrefmult");
// }else{
//  hvtx ->Add(gvtx);
//  hvtxgood ->Add(gvtxgood);
//  hrefmult ->Add(grefmult);
// }
//
}//ichain
cout << " done TChain-ing.. " << endl;

book_histo();

 gitom = new TF1("gitom","x*3.1415962*2*[0]*sqrt(1.115683*1.115683+x*x) * exp(-sqrt(1.115683*1.115683+x*x)/[1])", 0, 6);
  gitom->SetParameters(1.48179e+00,3.07372e-01);

// omega_tree.SetBranchAddress("brunid",&brunid);
// omega_tree.SetBranchAddress("beventid",&beventid);

/*
omega_tree.SetBranchAddress("refMult",&brefmult);
omega_tree.SetBranchAddress("mass",&bparticlemass);
omega_tree.SetBranchAddress("px",&bpx);
omega_tree.SetBranchAddress("py",&bpy);
omega_tree.SetBranchAddress("pz",&bpz);
*/

omega_tree.SetBranchAddress("brefmult",&brefmult);
omega_tree.SetBranchAddress("bparticlemass",&bparticlemass);
omega_tree.SetBranchAddress("bpx",&bpx);
omega_tree.SetBranchAddress("bpy",&bpy);
omega_tree.SetBranchAddress("bpz",&bpz);


omega_tree.SetBranchAddress("cent9",&cent9);
/*
omega_tree.SetBranchAddress("Chi2Topo", &om_chi2topo);
omega_tree.SetBranchAddress("Chi2NDF", &om_chi2ndf);
omega_tree.SetBranchAddress("LdL", &om_ldl);
omega_tree.SetBranchAddress("Chi2Topo_ld", &om_ld_chi2topo);
omega_tree.SetBranchAddress("Chi2NDF_ld", &om_ld_chi2ndf);
omega_tree.SetBranchAddress("LdL_ld", &om_ld_ldl);
omega_tree.SetBranchAddress("L_ld", &om_ld_l);
omega_tree.SetBranchAddress("L", &om_l);
omega_tree.SetBranchAddress("chi2Primary_Proton", &chi2primary_om_proton);
omega_tree.SetBranchAddress("chi2Primary_Pion", &chi2primary_om_pi);
omega_tree.SetBranchAddress("chi2Primary_Kaon", &chi2primary_om_bach);
*/

omega_tree.SetBranchAddress("om_chi2topo", &om_chi2topo);
omega_tree.SetBranchAddress("om_chi2ndf", &om_chi2ndf);
omega_tree.SetBranchAddress("om_ldl", &om_ldl);
omega_tree.SetBranchAddress("om_ld_chi2topo", &om_ld_chi2topo);
omega_tree.SetBranchAddress("om_ld_chi2ndf", &om_ld_chi2ndf);
omega_tree.SetBranchAddress("om_ld_ldl", &om_ld_ldl);
omega_tree.SetBranchAddress("om_ld_l", &om_ld_l);
omega_tree.SetBranchAddress("om_l", &om_l);
omega_tree.SetBranchAddress("chi2primary_om_proton", &chi2primary_om_proton);
omega_tree.SetBranchAddress("chi2primary_om_pi", &chi2primary_om_pi);
omega_tree.SetBranchAddress("chi2primary_om_bach", &chi2primary_om_bach);
 omega_tree.SetBranchAddress("nhits_om_proton", &nhits_om_proton);
 omega_tree.SetBranchAddress("nhits_om_pi", &nhits_om_pi);
 omega_tree.SetBranchAddress("nhits_om_bach", &nhits_om_bach);
 omega_tree.SetBranchAddress("dedx_om_proton", &dedx_om_proton);
 omega_tree.SetBranchAddress("dedx_om_pi", &dedx_om_pi);
 omega_tree.SetBranchAddress("dedx_om_bach", &dedx_om_bach);


// Long64_t n_lambda_Entries = lambda_tree.GetEntries();
Long64_t n_omega_Entries = omega_tree.GetEntries();

 cout<<endl<<"Start processing "<<endl<<endl;
cout<<"--------------------------------------------"<<endl;
cout<<"Total entries "<<n_omega_Entries<< endl;

om_width_0 = new TF1("om_width_0","0.0007*x+0.0008",1.6,10);
  om_width_1 = new TF1("om_width_1","0.00022*x+0.0016",0.,1.6);

 //for (Long64_t iEntry = 0; iEntry<=1e6; iEntry++){
  for (Long64_t iEntry = 0; iEntry<=n_omega_Entries; iEntry++){

	omega_tree.GetEntry(iEntry);
  if(iEntry%100000==0){
    cout<<"Processing entry:" <<iEntry<<endl;
}
int centtmva_label = -1;
int cent_label = -1;
int pt_label = -1;

int cent_cut_label = -1;
int pt_cut_label = -1;


if(snn==28){//temp
/*
    if(brefmult>6) cent_label = 7;
    // if(brefmult>6) cent_label = 8;
    // if(brefmult>12) cent_label = 7;
    if(brefmult>24) cent_label = 6;
    if(brefmult>42) cent_label = 5;
    if(brefmult>69) cent_label = 4;
    if(brefmult>107) cent_label = 3;
    if(brefmult>157) cent_label = 2;
    if(brefmult>226) cent_label = 1;
    if(brefmult>270) cent_label = 0;

    if(brefmult>6) cent_cut_label = 2;
    // if(brefmult>24) cent_cut_label = 6;
    // if(brefmult>42) cent_cut_label = 5;
    if(brefmult>69) cent_cut_label = 1;
    // if(brefmult>107) cent_cut_label = 3;
    // if(brefmult>157) cent_cut_label = 2;
    if(brefmult>226) cent_cut_label = 0;
    // if(brefmult>270) cent_cut_label = 0;
*/
    
      if(cent9==0) cent_label = 7;
      if(cent9==1) cent_label = 7;
      if(cent9==2) cent_label = 6;
      if(cent9==3) cent_label = 5;
      if(cent9==4) cent_label = 4;
      if(cent9==5) cent_label = 3;
      if(cent9==6) cent_label = 2;
      if(cent9==7) cent_label = 1;
      if(cent9==8) cent_label = 0;

      if(cent9==0) cent_cut_label = 2;
      if(cent9==1) cent_cut_label = 2;
      if(cent9==2) cent_cut_label = 2;
      if(cent9==3) cent_cut_label = 2;
      if(cent9==4) cent_cut_label = 1;
      if(cent9==5) cent_cut_label = 1;
      if(cent9==6) cent_cut_label = 1;
      if(cent9==7) cent_cut_label = 0;
      if(cent9==8) cent_cut_label = 0;

    centtmva_label = cent_cut_label;
}
  if(cent_label<0) continue;

  bpt = sqrt(bpx*bpx+bpy*bpy);

  TLorentzVector ptc(bpx,bpy,bpz,sqrt(bpx*bpx+bpy*bpy+bpz*bpz+bparticlemass*bparticlemass));
  if(fabs(ptc.Rapidity())>0.5) continue;

  if(bpt > ptcutbin[0]) pt_cut_label = 0;
  if(bpt > ptcutbin[1]) pt_cut_label = 1;
  if(bpt > ptcutbin[2]) pt_cut_label = 2;
  if(bpt > ptcutbin[3]) pt_cut_label = -1;


  if(pt_cut_label < 0) continue;


//  om_width_0 = new TF1("om_width_0","0.0007*x+0.0008",1.6,10);
//  om_width_1 = new TF1("om_width_1","0.00022*x+0.0016",0.,1.6);


//  gitom = new TF1("gitom","x*3.1415962*2*[0]*sqrt(1.115683*1.115683+x*x) * exp(-sqrt(1.115683*1.115683+x*x)/[1])", 0, 6);
//  gitom->SetParameters(1.48179e+00,3.07372e-01);


  if(bpt>1.6){
    sgm3 = om_width_0 -> Eval(bpt);
  }else{
    sgm3 = om_width_1 -> Eval(bpt);
  }

if(_data_or_mc==0){
  // _weight = gitom->Eval(sqrt(bpx*bpx + bpy*bpy));//approx.
  _weight = gitom->Eval(bpt);//approx.
}else{
  _weight=1;
}

   TLorentzVector ptc(bpx,bpy,bpz,sqrt(bpx*bpx+bpy*bpy+bpz*bpz+bparticlemass*bparticlemass));
   if(fabs(ptc.Rapidity())>0.5) continue;


  if(om_chi2topo<0) continue;



  Chi2Topo = om_chi2topo;
  Chi2NDF = om_chi2ndf;
  LdL = om_ldl;
  Chi2Topo_ld = om_ld_chi2topo;
  Chi2NDF_ld = om_ld_chi2ndf;
  LdL_ld = om_ld_ldl;
  L_ld = om_ld_l;
  L = om_l;
  chi2Primary_Proton = chi2primary_om_proton;
  chi2Primary_Pion = chi2primary_om_pi;
  chi2Primary_Kaon = chi2primary_om_bach;

//    mvaValue = reader[centtmva_label][0]->EvaluateMVA( "BDT" );
//  h_om_mass[cent_cut_label][pt_cut_label]->Fill(bparticlemass);



if(bparticlemass<omega_mass+3*sgm3 && bparticlemass>omega_mass-3*sgm3){
  if(_data_or_mc==0){

/*
  h_om_chi2topo[cent_cut_label][pt_cut_label] -> Fill(om_chi2topo,_weight);
  h_om_chi2ndf[cent_cut_label][pt_cut_label] -> Fill(om_chi2ndf,_weight);
  h_om_ldl[cent_cut_label][pt_cut_label] -> Fill(om_ldl,_weight);
  h_om_ld_chi2topo[cent_cut_label][pt_cut_label] -> Fill(om_ld_chi2topo,_weight);
  h_om_ld_chi2ndf[cent_cut_label][pt_cut_label] -> Fill(om_ld_chi2ndf,_weight);
  h_om_ld_ldl[cent_cut_label][pt_cut_label] -> Fill(om_ld_ldl,_weight);
  h_om_ld_l[cent_cut_label][pt_cut_label] -> Fill(om_ld_l,_weight);
  h_om_l[cent_cut_label][pt_cut_label] -> Fill(om_l,_weight);
  h_chi2primary_om_proton[cent_cut_label][pt_cut_label] -> Fill(chi2primary_om_proton,_weight);
  h_chi2primary_om_pi[cent_cut_label][pt_cut_label] -> Fill(chi2primary_om_pi,_weight);
  h_chi2primary_om_bach[cent_cut_label][pt_cut_label] -> Fill(chi2primary_om_bach,_weight);
  h_om_mva[cent_cut_label][pt_cut_label] -> Fill(mvaValue,_weight);

  h_phieta->Fill(ptc.Phi(), ptc.Eta(), _weight);
*/
//cout<<"nhit:"<<cent_cut_label<<" "<<pt_cut_label<<" "<<nhits_om_proton<<endl;

 h_nhits_om_proton[cent_cut_label][pt_cut_label] -> Fill(nhits_om_proton,_weight);
 h_nhits_om_pi[cent_cut_label][pt_cut_label] -> Fill(nhits_om_pi,_weight);
 h_nhits_om_bach[cent_cut_label][pt_cut_label] -> Fill(nhits_om_bach,_weight);
 h_dedx_om_proton[cent_cut_label][pt_cut_label] -> Fill(dedx_om_proton,_weight);
 h_dedx_om_pi[cent_cut_label][pt_cut_label] -> Fill(dedx_om_pi,_weight);
 h_dedx_om_bach[cent_cut_label][pt_cut_label] -> Fill(dedx_om_bach,_weight);

  }else{

/*
  h_om_chi2topo[cent_cut_label][pt_cut_label] -> Fill(om_chi2topo);
  h_om_chi2ndf[cent_cut_label][pt_cut_label] -> Fill(om_chi2ndf);
  h_om_ldl[cent_cut_label][pt_cut_label] -> Fill(om_ldl);
  h_om_ld_chi2topo[cent_cut_label][pt_cut_label] -> Fill(om_ld_chi2topo);
  h_om_ld_chi2ndf[cent_cut_label][pt_cut_label] -> Fill(om_ld_chi2ndf);
  h_om_ld_ldl[cent_cut_label][pt_cut_label] -> Fill(om_ld_ldl);
  h_om_ld_l[cent_cut_label][pt_cut_label] -> Fill(om_ld_l);
  h_om_l[cent_cut_label][pt_cut_label] -> Fill(om_l);
  h_chi2primary_om_proton[cent_cut_label][pt_cut_label] -> Fill(chi2primary_om_proton);
  h_chi2primary_om_pi[cent_cut_label][pt_cut_label] -> Fill(chi2primary_om_pi);
  h_chi2primary_om_bach[cent_cut_label][pt_cut_label] -> Fill(chi2primary_om_bach);
  h_om_mva[cent_cut_label][pt_cut_label] -> Fill(mvaValue);

  h_phieta->Fill(ptc.Phi(), ptc.Eta());
*/

//cout<<"nhit:"<<cent_cut_label<<" "<<pt_cut_label<<" "<<nhits_om_proton<<" "<<om_chi2topo<<" "<<  bparticlemass<<endl;

 h_nhits_om_proton[cent_cut_label][pt_cut_label] -> Fill(nhits_om_proton);
 h_nhits_om_pi[cent_cut_label][pt_cut_label] -> Fill(nhits_om_pi);
 h_nhits_om_bach[cent_cut_label][pt_cut_label] -> Fill(nhits_om_bach);
 h_dedx_om_proton[cent_cut_label][pt_cut_label] -> Fill(dedx_om_proton);
 h_dedx_om_pi[cent_cut_label][pt_cut_label] -> Fill(dedx_om_pi);
 h_dedx_om_bach[cent_cut_label][pt_cut_label] -> Fill(dedx_om_bach);
  }

}

if(  (bparticlemass<omega_mass+6*sgm3 && bparticlemass>omega_mass+3*sgm3) ||  (bparticlemass>omega_mass-6*sgm3 && bparticlemass<omega_mass-3*sgm3)  ){
if(_data_or_mc==1){
  h_om_chi2topo_bg[cent_cut_label][pt_cut_label] -> Fill(om_chi2topo);
  h_om_chi2ndf_bg[cent_cut_label][pt_cut_label] -> Fill(om_chi2ndf);
  h_om_ldl_bg[cent_cut_label][pt_cut_label] -> Fill(om_ldl,_weight);
  h_om_ld_chi2topo_bg[cent_cut_label][pt_cut_label] -> Fill(om_ld_chi2topo);
  h_om_ld_chi2ndf_bg[cent_cut_label][pt_cut_label] -> Fill(om_ld_chi2ndf);
  h_om_ld_ldl_bg[cent_cut_label][pt_cut_label] -> Fill(om_ld_ldl);
  h_om_ld_l_bg[cent_cut_label][pt_cut_label] -> Fill(om_ld_l);
  h_om_l_bg[cent_cut_label][pt_cut_label] -> Fill(om_l);
  h_chi2primary_om_proton_bg[cent_cut_label][pt_cut_label] -> Fill(chi2primary_om_proton);
  h_chi2primary_om_pi_bg[cent_cut_label][pt_cut_label] -> Fill(chi2primary_om_pi);
  h_chi2primary_om_bach_bg[cent_cut_label][pt_cut_label] -> Fill(chi2primary_om_bach);
  h_om_mva_bg[cent_cut_label][pt_cut_label] -> Fill(mvaValue);

  h_phieta_bg->Fill(ptc.Phi(), ptc.Eta());

 h_nhits_om_proton_bg[cent_cut_label][pt_cut_label] -> Fill(nhits_om_proton);
 h_nhits_om_pi_bg[cent_cut_label][pt_cut_label] -> Fill(nhits_om_pi);
 h_nhits_om_bach_bg[cent_cut_label][pt_cut_label] -> Fill(nhits_om_bach);
 h_dedx_om_proton_bg[cent_cut_label][pt_cut_label] -> Fill(dedx_om_proton);
 h_dedx_om_pi_bg[cent_cut_label][pt_cut_label] -> Fill(dedx_om_pi);
 h_dedx_om_bach_bg[cent_cut_label][pt_cut_label] -> Fill(dedx_om_bach);
 
    }
  }
}//entry loop


write_histo();

}

void book_histo()
{
//for(int icent=0;icent<ncenttmvabins;icent++){
for(int icent=0;icent<ncentbins;icent++){

  mkfIncMassPtUL3_0[icent] = new TH2F(Form("mkfIncMassPtUL3_0[%d]",icent), Form("mkfIncMassPtUL3_0[%d]",icent), 1000, 1, 2, 1000, 0, 10);
  mkfIncMassPtUL4_0[icent] = new TH2F(Form("mkfIncMassPtUL4_0[%d]",icent), Form("mkfIncMassPtUL4_0[%d]",icent), 1000, 1, 2, 1000, 0, 10);
  mkfIncMassPtUL3_1[icent] = new TH2F(Form("mkfIncMassPtUL3_1[%d]",icent), Form("mkfIncMassPtUL3_1[%d]",icent), 1000, 1, 2, 1000, 0, 10);
  mkfIncMassPtUL4_1[icent] = new TH2F(Form("mkfIncMassPtUL4_1[%d]",icent), Form("mkfIncMassPtUL4_1[%d]",icent), 1000, 1, 2, 1000, 0, 10);

  mkfIncMassPtUL3_0[icent]->Sumw2();
  mkfIncMassPtUL4_0[icent]->Sumw2();
  mkfIncMassPtUL3_1[icent]->Sumw2();
  mkfIncMassPtUL4_1[icent]->Sumw2();

	}

h_phieta = new TH2F("h_phieta",";#phi;#eta",100,-TMath::Pi(),TMath::Pi(), 100,-1.5,1.5);
h_phieta_bg = new TH2F("h_phieta_bg",";#phi;#eta",100,-TMath::Pi(),TMath::Pi(), 100,-1.5,1.5);

for(int icent=0;icent<ncutcentbins;icent++){
//  h_om_mass[icent] = new TH1F(Form("h_om_mass[%d]",icent), Form("h_om_mass[%d]",icent), 1000,1,2);
  for(int ipt=0;ipt<ncutptbins;ipt++){
 h_om_mass[icent][ipt] = new TH1F(Form("h_om_mass[%d][%d]",icent,ipt), Form("h_om_mass[%d][%d]",icent,ipt), 1000,1,2);

  h_om_chi2topo[icent][ipt] = new TH1F(Form("h_om_chi2topo[%d][%d]",icent,ipt), Form("h_om_chi2topo[%d][%d]",icent,ipt), 50,0,10);
  h_om_chi2ndf[icent][ipt] = new TH1F(Form("h_om_chi2ndf[%d][%d]",icent,ipt), Form("h_om_chi2ndf[%d][%d]",icent,ipt), 50,0,10);

  h_om_ldl[icent][ipt] = new TH1F(Form("h_om_ldl[%d][%d]",icent,ipt), Form("h_om_ldl[%d][%d]",icent,ipt), 100, 0, 80);
  h_om_ld_chi2topo[icent][ipt] = new TH1F(Form("h_om_ld_chi2topo[%d][%d]",icent,ipt), Form("h_om_ld_chi2topo[%d][%d]",icent,ipt), 50,0,50);
  h_om_ld_chi2ndf[icent][ipt] = new TH1F(Form("h_om_ld_chi2ndf[%d][%d]",icent,ipt), Form("h_om_ld_chi2ndf[%d][%d]",icent,ipt), 50,0,10);
  h_om_ld_ldl[icent][ipt] = new TH1F(Form("h_om_ld_ldl[%d][%d]",icent,ipt), Form("h_om_ld_ldl[%d][%d]",icent,ipt), 100, 0, 100);
  h_om_ld_l[icent][ipt] = new TH1F(Form("h_om_ld_l[%d][%d]",icent,ipt), Form("h_om_ld_l[%d][%d]",icent,ipt), 100, 0, 100);
  h_om_l[icent][ipt] = new TH1F(Form("h_om_l[%d][%d]",icent,ipt), Form("h_om_l[%d][%d]",icent,ipt), 100, 0, 50);

  h_om_mva[icent][ipt] = new TH1F(Form("h_om_mva[%d][%d]",icent,ipt), Form("h_om_mva[%d][%d]",icent,ipt), 100, -0.5,0.5);

  h_om_chi2topo_bg[icent][ipt] = new TH1F(Form("h_om_chi2topo_bg[%d][%d]",icent,ipt), Form("h_om_chi2topo_bg[%d][%d]",icent,ipt), 50,0,10);
  h_om_chi2ndf_bg[icent][ipt] = new TH1F(Form("h_om_chi2ndf_bg[%d][%d]",icent,ipt), Form("h_om_chi2ndf_bg[%d][%d]",icent,ipt), 50,0,10);

  h_om_ldl_bg[icent][ipt] = new TH1F(Form("h_om_ldl_bg[%d][%d]",icent,ipt), Form("h_om_ldl_bg[%d][%d]",icent,ipt), 100, 0, 80);
  h_om_ld_chi2topo_bg[icent][ipt] = new TH1F(Form("h_om_ld_chi2topo_bg[%d][%d]",icent,ipt), Form("h_om_ld_chi2topo_bg[%d][%d]",icent,ipt), 50,0,50);
  h_om_ld_chi2ndf_bg[icent][ipt] = new TH1F(Form("h_om_ld_chi2ndf_bg[%d][%d]",icent,ipt), Form("h_om_ld_chi2ndf_bg[%d][%d]",icent,ipt), 50,0,10);
  h_om_ld_ldl_bg[icent][ipt] = new TH1F(Form("h_om_ld_ldl_bg[%d][%d]",icent,ipt), Form("h_om_ld_ldl_bg[%d][%d]",icent,ipt), 100, 0, 100);
  h_om_ld_l_bg[icent][ipt] = new TH1F(Form("h_om_ld_l_bg[%d][%d]",icent,ipt), Form("h_om_ld_l_bg[%d][%d]",icent,ipt), 100, 0, 100);
  h_om_l_bg[icent][ipt] = new TH1F(Form("h_om_l_bg[%d][%d]",icent,ipt), Form("h_om_l_bg[%d][%d]",icent,ipt), 100, 0, 50);

  h_om_mva_bg[icent][ipt] = new TH1F(Form("h_om_mva_bg[%d][%d]",icent,ipt), Form("h_om_mva_bg[%d][%d]",icent,ipt), 100, -0.5,0.5);

  h_chi2primary_om_proton[icent][ipt] = new TH1F(Form("h_chi2primary_om_proton[%d][%d]",icent,ipt), Form("h_chi2primary_om_proton[%d][%d]",icent,ipt), 100,0,1000);
  h_chi2primary_om_pi[icent][ipt] = new TH1F(Form("h_chi2primary_om_pi[%d][%d]",icent,ipt), Form("h_chi2primary_om_pi[%d][%d]",icent,ipt), 100,0,1000);
  h_chi2primary_om_bach[icent][ipt] = new TH1F(Form("h_chi2primary_om_bach[%d][%d]",icent,ipt), Form("h_chi2primary_om_bach[%d][%d]",icent,ipt), 100,0,1000);

  h_chi2primary_om_proton_bg[icent][ipt] = new TH1F(Form("h_chi2primary_om_proton_bg[%d][%d]",icent,ipt), Form("h_chi2primary_om_proton_bg[%d][%d]",icent,ipt), 100,0,1000);
  h_chi2primary_om_pi_bg[icent][ipt] = new TH1F(Form("h_chi2primary_om_pi_bg[%d][%d]",icent,ipt), Form("h_chi2primary_om_pi_bg[%d][%d]",icent,ipt), 100,0,1000);
  h_chi2primary_om_bach_bg[icent][ipt] = new TH1F(Form("h_chi2primary_om_bach_bg[%d][%d]",icent,ipt), Form("h_chi2primary_om_bach_bg[%d][%d]",icent,ipt), 100,0,1000);


h_nhits_om_proton[icent][ipt] = new TH1F(Form("h_nhits_om_proton[%d][%d]",icent,ipt), Form("h_nhits_om_proton[%d][%d]",icent,ipt),50,0,50);
h_nhits_om_pi[icent][ipt] = new TH1F(Form("h_nhits_om_pi[%d][%d]",icent,ipt), Form("h_nhits_om_pi[%d][%d]",icent,ipt),50,0,50);
h_nhits_om_bach[icent][ipt] = new TH1F(Form("h_nhits_om_bach[%d][%d]",icent,ipt), Form("h_nhits_om_bach[%d][%d]",icent,ipt),50,0,50);
h_dedx_om_proton[icent][ipt] = new TH1F(Form("h_dedx_om_proton[%d][%d]",icent,ipt), Form("h_dedx_om_proton[%d][%d]",icent,ipt),50,0,50);
h_dedx_om_pi[icent][ipt] = new TH1F(Form("h_dedx_om_pi[%d][%d]",icent,ipt), Form("h_dedx_om_pi[%d][%d]",icent,ipt),50,0,50);
h_dedx_om_bach[icent][ipt] = new TH1F(Form("h_dedx_om_bach[%d][%d]",icent,ipt), Form("h_dedx_om_bach[%d][%d]",icent,ipt),50,0,50);

h_nhits_om_proton_bg[icent][ipt] = new TH1F(Form("h_nhits_om_proton_bg[%d][%d]",icent,ipt), Form("h_nhits_om_proton_bg[%d][%d]",icent,ipt),50,0,50);
h_nhits_om_pi_bg[icent][ipt] = new TH1F(Form("h_nhits_om_pi_bg[%d][%d]",icent,ipt), Form("h_nhits_om_pi_bg[%d][%d]",icent,ipt),50,0,50);
h_nhits_om_bach_bg[icent][ipt] = new TH1F(Form("h_nhits_om_bach_bg[%d][%d]",icent,ipt), Form("h_nhits_om_bach_bg[%d][%d]",icent,ipt),50,0,50);
h_dedx_om_proton_bg[icent][ipt] = new TH1F(Form("h_dedx_om_proton_bg[%d][%d]",icent,ipt), Form("h_dedx_om_proton_bg[%d][%d]",icent,ipt),50,0,50);
h_dedx_om_pi_bg[icent][ipt] = new TH1F(Form("h_dedx_om_pi_bg[%d][%d]",icent,ipt), Form("h_dedx_om_pi_bg[%d][%d]",icent,ipt),50,0,50);
h_dedx_om_bach_bg[icent][ipt] = new TH1F(Form("h_dedx_om_bach_bg[%d][%d]",icent,ipt), Form("h_dedx_om_bach_bg[%d][%d]",icent,ipt),50,0,50);



  }
}

}

void write_histo()
{
    TFile *outhistfile = new TFile (ofile, "RECREATE");
    outhistfile->cd();

//	for(int icent=0;icent<ncenttmvabins;icent++){
	for(int icent=0;icent<ncentbins;icent++){
    mkfIncMassPtUL3_0[icent]->Write();
    mkfIncMassPtUL4_0[icent]->Write();
    // mkfIncMassPtUL3_1[icent]->Write();
    // mkfIncMassPtUL4_1[icent]->Write();
  }

   // hvtx -> Write();
   // hvtxgood -> Write();
   // hrefmult -> Write();

   h_phieta -> Write();
   h_phieta_bg -> Write();

   for(int icent=0;icent<ncutcentbins;icent++){
//     h_om_mass[icent]->Write();
     for(int ipt=0;ipt<ncutptbins;ipt++){
          h_om_mass[icent][ipt]->Write();

     h_om_chi2topo[icent][ipt]->Write();
     h_om_chi2ndf[icent][ipt]->Write();
     h_om_ldl[icent][ipt]->Write();
     h_om_ld_chi2topo[icent][ipt]->Write();
     h_om_ld_chi2ndf[icent][ipt]->Write();
     h_om_ld_ldl[icent][ipt]->Write();
     h_om_ld_l[icent][ipt]->Write();
     h_om_l[icent][ipt]->Write();

     h_chi2primary_om_proton[icent][ipt]->Write();
     h_chi2primary_om_pi[icent][ipt]->Write();
     h_chi2primary_om_bach[icent][ipt]->Write();

     h_om_mva[icent][ipt]->Write();

h_nhits_om_proton[icent][ipt] ->Write();
h_nhits_om_pi[icent][ipt] ->Write();
h_nhits_om_bach[icent][ipt] ->Write();
h_dedx_om_proton[icent][ipt] ->Write();
h_dedx_om_pi[icent][ipt] ->Write();
h_dedx_om_bach[icent][ipt] ->Write();

     if(_data_or_mc==1){
     h_om_chi2topo_bg[icent][ipt]->Write();
     h_om_chi2ndf_bg[icent][ipt]->Write();
     h_om_ldl_bg[icent][ipt]->Write();
     h_om_ld_chi2topo_bg[icent][ipt]->Write();
     h_om_ld_chi2ndf_bg[icent][ipt]->Write();
     h_om_ld_ldl_bg[icent][ipt]->Write();
     h_om_ld_l_bg[icent][ipt]->Write();
     h_om_l_bg[icent][ipt]->Write();

     h_chi2primary_om_proton_bg[icent][ipt]->Write();
     h_chi2primary_om_pi_bg[icent][ipt]->Write();
     h_chi2primary_om_bach_bg[icent][ipt]->Write();

     h_om_mva_bg[icent][ipt]->Write();

h_nhits_om_proton_bg[icent][ipt]->Write();
h_nhits_om_pi_bg[icent][ipt] ->Write();
h_nhits_om_bach_bg[icent][ipt] ->Write();
h_dedx_om_proton_bg[icent][ipt]->Write();
h_dedx_om_pi_bg[icent][ipt] ->Write();
h_dedx_om_bach_bg[icent][ipt] ->Write();

      }

     }
   }

}
