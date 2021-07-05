//THIS MACRO REQUIRES TMVA package, better run on rcf
//write mass/pt histograms that can be read offline

//v2 tmva for corrected chi2primiary with vertex error
//testing v5 bdt: ksigma3 with default params
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"


#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif


const int npttmvabins = 1;
//const int nptbins = 6;
const int nptbins = 8;
const int ncenttmvabins = 3;
const int ncentbins = 8;

 TH1F	*hvtx;
 TH1F	*hvtxgood;
 TH1F	*hrefmult;
 TH1F	*wrefmult;

 TH1F   *gvtx;
 TH1F   *gvtxgood;
 TH1F   *grefmult;
 TH1F   *frefmult;

//TH2F *mkfIncMassPtUL3_0[ncenttmvabins];
//TH2F *mkfIncMassPtUL3_1[ncenttmvabins];
//TH2F *mkfIncMassPtUL4_0[ncenttmvabins];
//TH2F *mkfIncMassPtUL4_1[ncenttmvabins];
//
TH2F *mkfIncMassPtUL3_0[ncentbins];
TH2F *mkfIncMassPtUL3_1[ncentbins];
TH2F *mkfIncMassPtUL4_0[ncentbins];
TH2F *mkfIncMassPtUL4_1[ncentbins];

bool _applyweight;

// double ptbin[npttmvabins+1] = {0.0,1.0,2.0,3.0,5.0,10.0};
//double ptbin[nptbins+1] = {0.4,0.8,1.2,1.6,2.0,2.4,2.8};
double ptbin[nptbins+1] = {0.4,0.8,1.2,1.6,2.0,2.4,2.8,3.2,3.6};

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
int    brefmult        ;
int    btofmult        ;
float  bVz             ;
int    bparticleid     ;
float  bparticlemass   ;
float  bpx             ;
float  bpy             ;
float  bpz             ;
float  bpt	       ;
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

float ld_chi2topo     ;
float ld_chi2ndf      ;
float ld_ldl          ;
float ld_l            ;
float ld_dl           ;
float chi2primary_proton ;
float chi2primary_pi     ;

float Chi2NDF;
float LdL;
float Chi2Topo;
float L;
float chi2Primary_Pi;
float chi2Primary_Proton;

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

int notbadrun;

double reweight;
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
int    _condor;
Char_t path_ch[10000];
Char_t ofile[10000];

TMVA::Reader* reader[ncenttmvabins][npttmvabins];
// float bdt_ld_cut[ncenttmvabins][npttmvabins];
float bdt_ld_cut[ncenttmvabins][nptbins];


void book_histo();
void write_histo();
void process_tree_tmva_om_kf_v3(int condor=0, int cut_mode=0){
    _condor = condor;//condor=1: condor settings
    _cut_mode = cut_mode;

 // This loads the library
     TMVA::Tools::Instance();

//CHANGE NAME TO MATCH READER>>>

//TMVA::Reader* reader = new TMVA::Reader( "!Color:!Silent" );

for(int icent=0;icent<ncenttmvabins;icent++){
   for(int ipt=0;ipt<npttmvabins;ipt++){

      reader[icent][ipt] =  new TMVA::Reader( "!Color:!Silent" );

      // reader[icent][ipt]->AddVariable("Chi2NDF", &Chi2NDF);
      // reader[icent][ipt]->AddVariable("LdL", &LdL);
      // reader[icent][ipt]->AddVariable("Chi2Topo", &Chi2Topo);
      // reader[icent][ipt]->AddVariable("L", &L);
      // reader[icent][ipt]->AddVariable("chi2Primary_Pi", &chi2Primary_Pi);
      // reader[icent][ipt]->AddVariable("chi2Primary_Proton", &chi2Primary_Proton);

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


      // reader[icent][ipt]->BookMVA( "BDT" , Form("/star/u/yhleung2/pwg/06.tmva/weights/tmva_om_weight_kf_cent_%d_pt_%d_BDT.weights.xml",icent,ipt));
      //reader[icent][ipt]->BookMVA( "BDT" , Form("/star/u/yhleung2/pwg/06.tmva/weights/tmva_om_weight_kf_cent_%d_pt_%d_v2_BDT.weights.xml",icent,ipt));
      //reader[icent][ipt]->BookMVA( "BDT" , Form("/star/u/yhleung2/pwg/06.tmva/weights/tmva_om_weight_kf_cent_%d_pt_%d_v4_BDT.weights.xml",icent,ipt));
//      reader[icent][ipt]->BookMVA( "BDT" , Form("/star/u/yhleung2/pwg/06.tmva/weights/tmva_om_weight_kf_cent_%d_pt_%d_v5_BDT.weights.xml",icent,ipt));//ncut 11
      reader[icent][ipt]->BookMVA( "BDT" , Form("/star/u/yhleung2/pwg/06.tmva/weights/tmva_om_weight_kf_cent_%d_pt_%d_v6_BDT.weights.xml",icent,ipt));// ncut 20

// if(ipt<4 && icent==0){
//       reader[icent][ipt]->BookMVA( "BDT" , Form("/star/u/yhleung2/pwg/06.tmva/weights/tmva_ld_weight_kf_cent_%d_pt_%d_BDT.weights.xml",icent,ipt));
// }
// else if(ipt==4 && icent==0){
//       reader[icent][ipt]->BookMVA( "BDT" ,"/star/u/yhleung2/pwg/06.tmva/weights/tmva_ld_weight_kf_cent_0_pt_3_BDT.weights.xml");
// }
// else if(ipt<4 && (icent==1||icent==2)){
//       reader[icent][ipt]->BookMVA( "BDT" ,Form("/star/u/yhleung2/pwg/06.tmva/weights/tmva_ld_weight_kf_cent_1_pt_%d_BDT.weights.xml",ipt));
// }
// else if(ipt==4 && (icent==1||icent==2)){
//       reader[icent][ipt]->BookMVA( "BDT" ,"/star/u/yhleung2/pwg/06.tmva/weights/tmva_ld_weight_kf_cent_1_pt_3_BDT.weights.xml");
// }//else{
// //      reader[icent][ipt]->BookMVA( "BDT" ,"/star/u/yhleung2/pwg/06.tmva/weights/tmva_ld_weight_cent_1_pt_2_BDT.weights.xml");
// //}

}//pt
}//cent


//cuts with wrong run11 no vtx err
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


//cuts with wrong run11 with manual vtx err
/*
bdt_ld_cut[0][0] = 0.0506;
bdt_ld_cut[0][1] = -0.0197;
bdt_ld_cut[0][2] = -0.0773;
bdt_ld_cut[0][3] = -0.0773;
bdt_ld_cut[0][4] = -0.0773;
bdt_ld_cut[0][5] = -0.0773;

bdt_ld_cut[1][0] = 0.0135;
bdt_ld_cut[1][1] = -0.0724;
bdt_ld_cut[1][2] = -0.1112;
bdt_ld_cut[1][3] = -0.1112;
bdt_ld_cut[1][4] = -0.1112;
bdt_ld_cut[1][5] = -0.1112;

bdt_ld_cut[2][0] = -0.0800;
bdt_ld_cut[2][1] = -0.1618;
bdt_ld_cut[2][2] = -0.2233;
bdt_ld_cut[2][3] = -0.2233;
bdt_ld_cut[2][4] = -0.2233;
bdt_ld_cut[2][5] = -0.2233;
*/

//ksigma 2 case
/*
bdt_ld_cut[0][0] = 0.1306;
bdt_ld_cut[0][1] = 0.0802;
bdt_ld_cut[0][2] = 0.0597;
bdt_ld_cut[0][3] = 0.0315;
bdt_ld_cut[0][4] = 0.0315;
bdt_ld_cut[0][5] = 0.0315;

bdt_ld_cut[1][0] = 0.0929;
bdt_ld_cut[1][1] = 0.0493;
bdt_ld_cut[1][2] = -0.0122;
bdt_ld_cut[1][3] = -0.0122;
bdt_ld_cut[1][4] = -0.0122;
bdt_ld_cut[1][5] = -0.0146;

bdt_ld_cut[2][0] = 0.0958;
bdt_ld_cut[2][1] = -0.0226;
bdt_ld_cut[2][2] = -0.0789;
bdt_ld_cut[2][3] = -0.0789;
bdt_ld_cut[2][4] = -0.0789;
bdt_ld_cut[2][5] = -0.0789;
*/

//ksigma3 case
/*
bdt_ld_cut[0][0] = 0.1925;
bdt_ld_cut[0][1] = 0.1538;
bdt_ld_cut[0][2] = 0.1082;
bdt_ld_cut[0][3] = 0.1082;
bdt_ld_cut[0][4] = 0.1082;
bdt_ld_cut[0][5] = 0.1082;

bdt_ld_cut[1][0] = 0.1543; 
bdt_ld_cut[1][1] = 0.1054; 
bdt_ld_cut[1][2] = 0.0682;
bdt_ld_cut[1][3] = 0.0575;
bdt_ld_cut[1][4] = 0.0575;
bdt_ld_cut[1][5] = 0.0575;

bdt_ld_cut[2][0] = 0.0426; 
bdt_ld_cut[2][1] = -0.0283;
bdt_ld_cut[2][2] = -0.0283;
bdt_ld_cut[2][3] = -0.0030;
bdt_ld_cut[2][4] = -0.0030;
bdt_ld_cut[2][5] = -0.0030;
*/

//ksigma3 + v6bdt
bdt_ld_cut[0][0] = 0.1610;
bdt_ld_cut[0][1] = 0.1304;
bdt_ld_cut[0][2] = 0.0907;
bdt_ld_cut[0][3] = 0.0907;
bdt_ld_cut[0][4] = 0.0907;
bdt_ld_cut[0][5] = 0.0907;
bdt_ld_cut[0][5] = 0.0907;
bdt_ld_cut[0][6] = 0.0907;
bdt_ld_cut[0][7] = 0.0907;

bdt_ld_cut[1][0] = 0.1447; 
bdt_ld_cut[1][1] = 0.0790;
bdt_ld_cut[1][2] = 0.0462;
bdt_ld_cut[1][3] = 0.0462;
bdt_ld_cut[1][4] = 0.0462;
bdt_ld_cut[1][5] = 0.0462;
bdt_ld_cut[1][5] = 0.0462;
bdt_ld_cut[1][6] = 0.0462;
bdt_ld_cut[1][7] = 0.0462;

bdt_ld_cut[2][0] = 0.0301;
bdt_ld_cut[2][1] = 0.0021;
bdt_ld_cut[2][2] = -0.0251;
bdt_ld_cut[2][3] = -0.0251;
bdt_ld_cut[2][4] = -0.0251;
bdt_ld_cut[2][5] = -0.0251;
bdt_ld_cut[2][5] = -0.0251;
bdt_ld_cut[2][6] = -0.0251;
bdt_ld_cut[2][7] = -0.0251;

TChain lambda_tree("lambda_tree");
TChain cascade_tree("cascade_tree");
TChain omega_tree("omega_tree");
int nChains = -999999;
ifstream fin;

//fin.open("/star/u/yhleung2/pwg/02.long_rotate/taxi_training/filelist.txt");
//strcpy (ofile,"run18_27gev_tmva_train_test.root");
//fin.open("/star/u/yhleung2/pwg/02.long_rotate/taxi/filelist.txt");
//fin.open("/star/u/yhleung2/pwg/02.long_rotate/taxi/filelist_test.txt");
if(_condor==0){
// fin.open("/star/u/yhleung2/pwg/05.KFParticle_rotate/taxi_om_s2b/filelist.txt");
fin.open("/star/u/yhleung2/pwg/05.KFParticle_rotate/testfilelist.txt");
strcpy (ofile,Form("run18_27gev_om_kf_tmva_train_%d.root",_cut_mode));
//strcpy (ofile,Form("run18_27gev_tmva_train_%d_test.root",_cut_mode));
//nChains = 34;
//nChains = 5;
//nChains = 10;
nChains = 1;
}
if(_condor==1){
fin.open("filelist.txt");
strcpy (ofile,Form("ana_tree_cut_mode_%d.root",_cut_mode));
nChains = 1;
}

for(int i=0; i<nChains; i++){
      fin >> path_ch;
      // lambda_tree.Add(path_ch);
      // cout<<path_ch<<"  "<<lambda_tree.GetEntries() <<endl;

      omega_tree.Add(path_ch);
      cout<<path_ch<<"  "<<omega_tree.GetEntries() <<endl;


      TFile *f0 = new TFile(path_ch,"READ");
      gvtx = (TH1F*)f0->Get("hvtx");
      gvtxgood = (TH1F*)f0->Get("hvtxgood");
      grefmult = (TH1F*)f0->Get("hrefmult");
      frefmult = (TH1F*)f0->Get("wrefmult");

if(i==0){
 hvtx = (TH1F*)gvtx->Clone("hvtx");
 hvtxgood = (TH1F*)gvtxgood->Clone("hvtxgood");
 hrefmult = (TH1F*)grefmult->Clone("hrefmult");
 wrefmult = (TH1F*)frefmult->Clone("wrefmult");
}else{
 hvtx ->Add(gvtx);
 hvtxgood ->Add(gvtxgood);
 hrefmult ->Add(grefmult);
 wrefmult ->Add(frefmult);
}

}//ichain
cout << " done TChain-ing.. " << endl;

book_histo();
//

omega_tree.SetBranchAddress("brunid",&brunid);
omega_tree.SetBranchAddress("beventid",&beventid);
omega_tree.SetBranchAddress("brefmult",&brefmult);
omega_tree.SetBranchAddress("btofmult",&btofmult);
omega_tree.SetBranchAddress("bparticleid",&bparticleid);
omega_tree.SetBranchAddress("bparticlemass",&bparticlemass);
omega_tree.SetBranchAddress("bpx",&bpx);
omega_tree.SetBranchAddress("bpy",&bpy);
omega_tree.SetBranchAddress("bpz",&bpz);
// if(_data_or_sim==0){
// lambda_tree.SetBranchAddress("_rotate",&_rotate);
// }
//if(_data_or_sim==1){
//lambda_tree.SetBranchAddress("pionidtruth",&pionidtruth);
//lambda_tree.SetBranchAddress("protonidtruth",&protonidtruth);
//}
/*
lambda_tree.SetBranchAddress("ld_pv_dca_ld",&ld_pv_dca_ld);
lambda_tree.SetBranchAddress("p_pv_dca_ld",&p_pv_dca_ld);
lambda_tree.SetBranchAddress("pi_p_dca_ld",&pi_p_dca_ld);
lambda_tree.SetBranchAddress("pi_pv_dca_ld",&pi_pv_dca_ld);
lambda_tree.SetBranchAddress("decaylen_ld",&decaylen_ld);
lambda_tree.SetBranchAddress("ang0_ld",&ang0_ld);
*/

// lambda_tree.SetBranchAddress("ld_chi2topo",&ld_chi2topo);
// lambda_tree.SetBranchAddress("ld_chi2ndf",&ld_chi2ndf);
// lambda_tree.SetBranchAddress("ld_ldl",&ld_ldl);
// lambda_tree.SetBranchAddress("ld_l",&ld_l);
// lambda_tree.SetBranchAddress("chi2primary_proton",&chi2primary_proton);
// lambda_tree.SetBranchAddress("chi2primary_pi",&chi2primary_pi);

omega_tree.SetBranchAddress("om_chi2topo", &om_chi2topo);
omega_tree.SetBranchAddress("om_chi2ndf", &om_chi2ndf);
omega_tree.SetBranchAddress("om_ldl", &om_ldl);
omega_tree.SetBranchAddress("om_ld_chi2topo", &om_ld_chi2topo);
omega_tree.SetBranchAddress("om_ld_chi2ndf", &om_ld_chi2ndf);
omega_tree.SetBranchAddress("om_ld_ldl", &om_ld_ldl);
omega_tree.SetBranchAddress("om_ld_l", &om_ld_l);
omega_tree.SetBranchAddress("om_l", &om_l);
omega_tree.SetBranchAddress("om_dl", &om_dl);
omega_tree.SetBranchAddress("chi2primary_om_proton", &chi2primary_om_proton);
omega_tree.SetBranchAddress("chi2primary_om_pi", &chi2primary_om_pi);
omega_tree.SetBranchAddress("chi2primary_om_bach", &chi2primary_om_bach);
omega_tree.SetBranchAddress("chi2primary_om_ld", &chi2primary_om_ld);

omega_tree.SetBranchAddress("notbadrun", &notbadrun);

omega_tree.SetBranchAddress("reweight", &reweight);
omega_tree.SetBranchAddress("cent9", &cent9);



// Long64_t n_lambda_Entries = lambda_tree.GetEntries();
Long64_t n_omega_Entries = omega_tree.GetEntries();

 cout<<endl<<"Start processing "<<endl<<endl;
cout<<"--------------------------------------------"<<endl;
// cout<<"Total entries "<<n_lambda_Entries<< endl;
//" "<<n_cascade_Entries<<" "<<n_omega_Entries<<endl;

cout<<"Total entries "<<n_omega_Entries<< endl;

if(_cut_mode==1){
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
}
/*
if(_cut_mode==5){
om_chi2topo_cut = 5;
om_chi2ndf_cut = 6;
om_ldl_cut = 15;
om_ld_chi2topo_cut = 5;
om_ld_chi2ndf_cut = 10;
om_ld_ldl_cut = 5;
om_ld_l_cut = 5;
om_l_cut = 5;
chi2primary_om_proton_cut = 18.6;
chi2primary_om_pi_cut = 18.6;
chi2primary_om_bach_cut = 18.6;
}

if(_cut_mode==6){
om_chi2topo_cut = 5;
om_chi2ndf_cut = 6;
om_ldl_cut = 10;
om_ld_chi2topo_cut = 15;
om_ld_chi2ndf_cut = 10;
om_ld_ldl_cut = 5;
om_ld_l_cut = 5;
om_l_cut = 5;
chi2primary_om_proton_cut = 18.6;
chi2primary_om_pi_cut = 18.6;
chi2primary_om_bach_cut = 18.6;
}
*/

for (Long64_t iEntry = 0; iEntry<=n_omega_Entries; iEntry++){

	omega_tree.GetEntry(iEntry);
  if(iEntry%100000==0){
    cout<<"Processing entry:" <<iEntry<<endl;
}
int centtmva_label = -1;
int cent_label = -1;
int pt_label = -1;


if(notbadrun>0) continue;


if(snn==28){//temp
    // if(brefmult>6) cent_label = 7;
    // if(brefmult>24) cent_label = 6;
    // if(brefmult>42) cent_label = 5;
    // if(brefmult>69) cent_label = 4;
    // if(brefmult>107) cent_label = 3;
    // if(brefmult>157) cent_label = 2;
    // if(brefmult>226) cent_label = 1;
    // if(brefmult>270) cent_label = 0;

    if(cent9==0) cent_label = 7;
    if(cent9==1) cent_label = 7;
    if(cent9==2) cent_label = 6;
    if(cent9==3) cent_label = 5;
    if(cent9==4) cent_label = 4;
    if(cent9==5) cent_label = 3;
    if(cent9==6) cent_label = 2;
    if(cent9==7) cent_label = 1;
    if(cent9==8) cent_label = 0;

      // if(brefmult>6)   centtmva_label = 2;
      // if(brefmult>69)  centtmva_label = 1;
      // if(brefmult>226) centtmva_label = 0;

      if(cent9==0) centtmva_label = 2;
      if(cent9==1) centtmva_label = 2;
      if(cent9==2) centtmva_label = 2;
      if(cent9==3) centtmva_label = 2;
      if(cent9==4) centtmva_label = 1;
      if(cent9==5) centtmva_label = 1;
      if(cent9==6) centtmva_label = 1;
      if(cent9==7) centtmva_label = 0;
      if(cent9==8) centtmva_label = 0;

}
  if(cent_label<0) continue;//<=5refmult

  bpt = sqrt(bpx*bpx+bpy*bpy);

  if(bpt > ptbin[0]) pt_label = 0;
  if(bpt > ptbin[1]) pt_label = 1;
  if(bpt > ptbin[2]) pt_label = 2;
  if(bpt > ptbin[3]) pt_label = 3;
  if(bpt > ptbin[4]) pt_label = 4;
  if(bpt > ptbin[5]) pt_label = 5;
  if(bpt > ptbin[6]) pt_label = 6;
  if(bpt > ptbin[7]) pt_label = 7;
//  if(bpt > ptbin[8]) pt_label = 8;
  if(pt_label < 0) continue;

   _rotate = 0;
//   mvaValue = reader[cent_label][pt_label]->EvaluateMVA( "BDT" );`

  TLorentzVector ptc(bpx,bpy,bpz,sqrt(bpx*bpx+bpy*bpy+bpz*bpz+bparticlemass*bparticlemass));
  if(fabs(ptc.Rapidity())>0.5) continue;

  if(om_chi2topo<0) continue;

  _applyweight = false;
  if(brefmult<=60)  _applyweight = true;


//decide whether to use veto cut, or not. for now, lets NOT use the veto cuts

if(_cut_mode==0){


if(!_applyweight){
  if(bparticleid>0){
  mkfIncMassPtUL3_0[cent_label]->Fill(bparticlemass,sqrt(bpx*bpx+bpy*bpy));
  }else{
  mkfIncMassPtUL4_0[cent_label]->Fill(bparticlemass,sqrt(bpx*bpx+bpy*bpy));
  }
}else{
  if(bparticleid>0){
  mkfIncMassPtUL3_0[cent_label]->Fill(bparticlemass,sqrt(bpx*bpx+bpy*bpy),reweight);
  }else{
  mkfIncMassPtUL4_0[cent_label]->Fill(bparticlemass,sqrt(bpx*bpx+bpy*bpy),reweight);
  }
}

}

//if(_cut_mode==1 || _cut_mode==5 || _cut_mode==6){
if(_cut_mode==1){
  //TODO, default cuts
if(
  om_chi2topo < om_chi2topo_cut &&
  om_chi2ndf < om_chi2ndf_cut &&
  om_ldl > om_ldl_cut &&
  om_ld_chi2topo > om_ld_chi2topo_cut &&
  om_ld_chi2ndf < om_ld_chi2ndf_cut &&
  om_ld_ldl > om_ld_ldl_cut &&
  om_ld_l > om_ld_l_cut &&
  om_l > om_l_cut &&
  chi2primary_om_proton > chi2primary_om_proton_cut &&
  chi2primary_om_pi > chi2primary_om_pi_cut &&
  chi2primary_om_bach > chi2primary_om_bach_cut
){

  if(!_applyweight){
    if(bparticleid>0){
    mkfIncMassPtUL3_0[cent_label]->Fill(bparticlemass,sqrt(bpx*bpx+bpy*bpy));
    }else{
    mkfIncMassPtUL4_0[cent_label]->Fill(bparticlemass,sqrt(bpx*bpx+bpy*bpy));
    }
  }else{
    if(bparticleid>0){
    mkfIncMassPtUL3_0[cent_label]->Fill(bparticlemass,sqrt(bpx*bpx+bpy*bpy),reweight);
    }else{
    mkfIncMassPtUL4_0[cent_label]->Fill(bparticlemass,sqrt(bpx*bpx+bpy*bpy),reweight);
    }
  }

   }
}//cut mode 1
//if(_cut_mode==2 || _cut_mode==3 || _cut_mode==4){
if(_cut_mode>=2 || _cut_mode<=-2){

       //map name
       // Chi2Topo = ld_chi2topo;
       // Chi2NDF = ld_chi2ndf;
       // LdL =ld_ldl;
       // L = ld_l;
       // chi2Primary_Proton = chi2primary_proton;
       // chi2Primary_Pi = chi2primary_pi;


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

       //chi2primary_om_ld
       //end map name
/*
       if(_cut_mode==2){
         delta=0;
       }
       if(_cut_mode==3){
         delta=0.06;
       }
       if(_cut_mode==4){
         delta=0.12;
       }
*/
if(_cut_mode==2){
    delta=0;
  }
  if(_cut_mode==3){
    delta=0.04;
  }
  if(_cut_mode==4){
    delta=0.08;
  }
  if(_cut_mode==5){
    delta=0.12;
  }
  if(_cut_mode==6){
//    if(pt_label==0) delta=0.02;
//    if(pt_label==1) delta=0.06;
//    if(pt_label==2) delta=0.10;
//    if(pt_label==3) delta=0.12;
//    if(pt_label==4) delta=0.14;
//    if(pt_label==5) delta=0.16;
  delta=0.16;    
  }
  if(_cut_mode==7){
    delta=0.20;
  }
  if(_cut_mode==8){
     delta=0.24;
  }
  
  if(_cut_mode==-3){
    delta=-0.04;
  }
  if(_cut_mode==-4){
    delta=-0.08;
  }
  if(_cut_mode==-5){
    delta=-0.12;
  }

  if(_cut_mode==106){//+50% eff
    if(centtmva_label==0){
/*
    if(pt_label==0) delta=0.02;
    if(pt_label==1) delta=0.04;
    if(pt_label==2) delta=0.08;
    if(pt_label==3) delta=0.12;
    if(pt_label==4) delta=0.12;
    if(pt_label==5) delta=0.16;
    if(pt_label==6) delta=0.16;
    if(pt_label==7) delta=0.16;
*/
      delta = 0.04;
    }
  if(centtmva_label==1){
/*
    if(pt_label==0) delta=0.02;
    if(pt_label==1) delta=0.08;
    if(pt_label==2) delta=0.12;
    if(pt_label==3) delta=0.16;
    if(pt_label==4) delta=0.16;
    if(pt_label==5) delta=0.16;
    if(pt_label==6) delta=0.16;
    if(pt_label==7) delta=0.16;
*/
       delta = 0.08;
    }
  if(centtmva_label==2){
 /*
    if(pt_label==0) delta=0.12;
    if(pt_label==1) delta=0.16;
    if(pt_label==2) delta=0.20;
    if(pt_label==3) delta=0.20;
    if(pt_label==4) delta=0.20;
    if(pt_label==5) delta=0.20;
    if(pt_label==6) delta=0.20;
    if(pt_label==7) delta=0.20;
*/
      delta =0.12;
    } 
  }
  
  if(_cut_mode==-106){//-50%eff
   if(centtmva_label==0){
/*
    if(pt_label==0) delta=-0.04;
    if(pt_label==1) delta=-0.04;
    if(pt_label==2) delta=-0.08;
    if(pt_label==3) delta=-0.10;
    if(pt_label==4) delta=-0.12;
    if(pt_label==5) delta=-0.12;
    if(pt_label==6) delta=-0.16;
    if(pt_label==7) delta=-0.16;
*/   
             delta = -0.04;
    }
   if(centtmva_label==1){
/*
    if(pt_label==0) delta=-0.04;
    if(pt_label==1) delta=-0.06;
    if(pt_label==2) delta=-0.10;
    if(pt_label==3) delta=-0.12;
    if(pt_label==4) delta=-0.14;
    if(pt_label==5) delta=-0.14;
    if(pt_label==6) delta=-0.14;
    if(pt_label==7) delta=-0.14;
*/ 
             delta = -0.08;
    }
if(centtmva_label==2){
/*
    if(pt_label==0) delta=-0.12;
    if(pt_label==1) delta=-0.12;
    if(pt_label==2) delta=-0.16;
    if(pt_label==3) delta=-0.18;
    if(pt_label==4) delta=-0.20;
    if(pt_label==5) delta=-0.20;
    if(pt_label==6) delta=-0.20;
    if(pt_label==7) delta=-0.20;
*/
               delta = -0.12;
    } 
  }



       //TODO, only 1 pt bin for tmva due to statistics
       mvaValue = reader[centtmva_label][0]->EvaluateMVA( "BDT" );
//cout<<mvaValue<<" "<<bdt_ld_cut[cent_label][pt_label]<<endl;
           if(mvaValue > bdt_ld_cut[centtmva_label][pt_label] - delta){

        if(!_applyweight){
          if(bparticleid>0){
          mkfIncMassPtUL3_0[cent_label]->Fill(bparticlemass,sqrt(bpx*bpx+bpy*bpy));
          }else{
          mkfIncMassPtUL4_0[cent_label]->Fill(bparticlemass,sqrt(bpx*bpx+bpy*bpy));
          }
        }else{
          if(bparticleid>0){
          mkfIncMassPtUL3_0[cent_label]->Fill(bparticlemass,sqrt(bpx*bpx+bpy*bpy),reweight);
          }else{
          mkfIncMassPtUL4_0[cent_label]->Fill(bparticlemass,sqrt(bpx*bpx+bpy*bpy),reweight);
          }
        }


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

   hvtx -> Write();
   hvtxgood -> Write();
   hrefmult -> Write();
   wrefmult -> Write();

}
