//THIS MACRO REQUIRES TMVA package, better run on rcf
//write mass/pt histograms that can be read offline

//v2 tmva for corrected chi2primiary with vertex error
//v4: tmva, bdt scanned value for each pt-cent bin
//v5: manual cuts for very loose sample 
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
const int ncenttmvabins = 2;
//const int ncentbins = 8;
const int ncentbins = 9;

 TH1F	*hvtx;
 TH1F	*hvtxgood;
 TH1F	*hrefmult;
 TH1F	*wrefmult;

 TH1F   *gvtx;
 TH1F   *gvtxgood;
 TH1F   *grefmult;
 TH1F   *frefmult;

 TH2F  *h_dedx_p;

TH1F  *h_hl_pl;
TH1F  *h_hl_pl_bg;

 TH1F  *h_ht_mass;
 TH1F  *h_ht_mass_bg;

 TH1F  *h_ht_ldl_bg;
 TH1F  *h_ht_dl_bg;
 TH1F  *h_ht_l_bg;

 TH1F  *h_ht_ldl;
 TH1F  *h_ht_dl;
 TH1F  *h_ht_l;


//TH2F *mkfIncMassPtUL3_0[ncenttmvabins];
//TH2F *mkfIncMassPtUL3_1[ncenttmvabins];
//TH2F *mkfIncMassPtUL4_0[ncenttmvabins];
//TH2F *mkfIncMassPtUL4_1[ncenttmvabins];
//
TH2F *mkfIncMassPtUL3_0[ncentbins];
TH2F *mkfIncMassPtUL3_1[ncentbins];
TH2F *mkfIncMassPtUL4_0[ncentbins];
TH2F *mkfIncMassPtUL4_1[ncentbins];

TH2F *mkfIncMassPtUL5_0[ncentbins];
TH2F *mkfIncMassPtUL6_0[ncentbins];

TH3F *h_chi2primary_he;
TH3F *h_bnhits;
TH3F *h_bdedx;

TH2F *h_pt;
TH2F *g_pt;

TH3F *h_dca;
bool _applyweight;

// double ptbin[npttmvabins+1] = {0.0,1.0,2.0,3.0,5.0,10.0};
//double ptbin[nptbins+1] = {0.4,0.8,1.2,1.6,2.0,2.4,2.8};
double ptbin[nptbins+1] = {0.4,0.8,1.2,1.6,2.0,2.4,2.8,3.2,3.6};

int centtmva_limit[ncenttmvabins+1] = {0,10,80};

int cent_limit[ncentbins+1] = {0,5,10,20,30,40,50,60,70,80};
//int cent_limit[ncentbins+1] = {0,5,10,20,30,40,50,60,80};
bool _cent_off = true;
//
int snn=28;

int    _cut_mode       ;
//_cut_mode = 0;//loose
//_cut_mode = 1;//default
//_cut_mode = 2;//bdt

float xi_chi2topo     ;
float xi_chi2ndf      ;
float xi_ldl          ;
float xi_ld_chi2topo  ;
float xi_ld_chi2ndf   ;
float xi_ld_ldl       ;
float xi_ld_l         ;
float xi_l            ;
float xi_dl           ;
float chi2primary_xi_proton;
float chi2primary_xi_pi    ;
float chi2primary_xi_bach  ;
float chi2primary_xi_ld    ;

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
float  bp 	       ;
float  beta            ;
float  brap            ;
int    _rotate         ;
float bdca;
int bmcparticleid;
float bmcrawpx;
float bmcrawpy;
float bmcrawpz;

int bmcparticleid;
float bmcrawpx;
float bmcrawpy;
float bmcrawpz;

float bdedx;
int bnhits;

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

float chi2primary_he;
float chi2primary_pi;
float ht_ldl;
float ht_dl;
float ht_l;
float ht_chi2topo;
float ht_chi2ndf;


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

float xi_chi2topo_cut;
float xi_chi2ndf_cut;
float xi_ldl_cut;
float xi_ld_chi2topo_cut;
float xi_ld_chi2ndf_cut;
float xi_ld_ldl_cut;
float xi_ld_l_cut;
float xi_l_cut;
float xi_dl_cut;
float chi2primary_xi_proton_cut;
float chi2primary_xi_pi_cut;
float chi2primary_xi_bach_cut;
float chi2primary_xi_ld_cut;

float chi2Primary_Kaon;
float chi2Primary_Bach;
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

float ht_chi2topo_cut;
float ht_chi2ndf_cut;
float ht_ldl_cut;
float ht_l_cut;
float chi2primary_he_cut;
float chi2primary_pi_cut;

float hl_chi2topo_cut;
float hl_chi2ndf_cut;
float hl_ldl_cut;
float hl_l_cut;
float chi2primary_h4_cut;
float chi2primary_pi_cut;

double ht_mass;
double ht_width;
double hl_mass;
double hl_width;

int notbadrun;
float ycm;
double reweight;
int cent9;
/*
 ht_mass = 2.991;
 ht_width= 0.005;

 hl_mass = 3.9235;
 hl_width= 0.005;
*/

//from fit:
 ht_mass = 2.9926;
 ht_width= 0.005;

 hl_mass = 3.9239;
 hl_width= 0.005;

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

double bdt_o_cut;


float ht_chi2ndf;
float ht_ldl;
float ht_chi2topo;
float ht_dl;
float ht_l;
float chi2primary_he;
float chi2primary_pi;

float hl_chi2ndf;
float hl_ldl;
float hl_chi2topo;
float hl_dl;
float hl_l;
float chi2primary_h4;
int bismc;
TMVA::Reader* reader[ncenttmvabins][npttmvabins];
// float bdt_ld_cut[ncenttmvabins][npttmvabins];

//obsolete, used for bdt scanning
float bdt_ld_cut[ncenttmvabins][nptbins];

//used for om analysis
float bdt_om_cut[ncentbins][nptbins];
float bdt_ob_cut[ncentbins][nptbins];

//using for xi analyusis
const int ncentbdtbins = 3;
const int nptbdtbins = 7;
float bdt_xi_cut[ncentbdtbins][nptbdtbins];
double ptbdtbin[nptbdtbins+1] = {0.2,0.6,1.0,1.4,1.8,2.2,2.6,3.0};

void book_histo();
void write_histo();

void process_tree_tmva_he_kf_vmv5(int condor=0, int cut_mode=0){

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

      reader[icent][ipt]->AddVariable("chi2Primary_Bach", &chi2Primary_Bach);
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

      reader[icent][ipt]->BookMVA( "BDT" , Form("/star/u/yhleung2/pwg/06.tmva/weights/tmva_xi_weight_kf_cent_%d_pt_%d_v7_BDT.weights.xml",icent,ipt));

    }//pt
}//cent

double delta_bdtscan = -0.04;//*****This should be a minus sign

//cent00-10
bdt_xi_cut[0][0] = 0.1063;//30000
bdt_xi_cut[0][1] = 0.0388;//5069.46
bdt_xi_cut[0][2] = 0.0382;//2000
bdt_xi_cut[0][3] = -0.0208;//1300
bdt_xi_cut[0][4] = -0.0208;
bdt_xi_cut[0][5] = -0.0423;//700
bdt_xi_cut[0][6] = -0.0423;

//cent10-40
bdt_xi_cut[1][0] = 0.1109;//18868
bdt_xi_cut[1][1] = 0.0051;//1551.67
bdt_xi_cut[1][2] = -0.0285;//616
bdt_xi_cut[1][3] = -0.0476;//400
bdt_xi_cut[1][4] = -0.0476;//400
bdt_xi_cut[1][5] = -0.0980;//240
bdt_xi_cut[1][6] = -0.0980;//240

//cent40-80
bdt_xi_cut[2][0] = -0.0980;//240
bdt_xi_cut[2][1] = -0.1116;//90
bdt_xi_cut[2][2] = -0.1116;//90
bdt_xi_cut[2][3] = -0.1523;//60
bdt_xi_cut[2][4] = -0.1523;//60
bdt_xi_cut[2][5] = -0.1523;//60
bdt_xi_cut[2][6] = -0.1523;//60


TChain lambda_tree("lambda_tree");
TChain cascade_tree("cascade_tree");
TChain omega_tree("omega_tree");
TChain htriton_tree("htriton_tree");
TChain h4lambda_tree("h4lambda_tree");
TChain he3_tree("he3_tree");
TChain he3_mc_tree("he3_mc_tree");

int nChains = -999999;
ifstream fin;

if(_condor==0){
//27gev
//fin.open("/star/u/yhleung2/pwg/29.KFParticle_hy2/taxi_hy_ana_v2/filelist.txt");
//strcpy (ofile,Form("run18_27gev_ht_train_%d_test.root",_cut_mode));
//nChains = 1;
//
//fxt26
//fin.open("/star/u/yhleung2/pwg/29.KFParticle_hy2/taxi_om_anatree_26gev_loose/filelist.txt");
//strcpy (ofile,Form("run18_26fxtgev_ht_train_%d_test.root",_cut_mode));
//nChains = 1;
//
//fxt3
//fin.open("/star/u/yhleung2/pwg/29.KFParticle_hy2/taxi_om_anatree_3gev_loose/filelist.txt");
//strcpy (ofile,Form("run18_3fxtgev_ht_train_%d_test.root",_cut_mode));
//nChains = 1;
//
//hl+ht
//fxt26
//fin.open("/star/u/yhleung2/pwg/29.KFParticle_hy2/taxi_hl_anatree_26gev_loose/filelist.txt");
//strcpy (ofile,Form("run18_26fxtgev_hl_train_%d_test.root",_cut_mode));
//nChains = 1;
//fxt3
//fin.open("/star/u/yhleung2/pwg/29.KFParticle_hy2/taxi_hl_anatree_3gev_loose/filelist.txt");
//strcpy (ofile,Form("run18_3fxtgev_hl_train_%d_test.root",_cut_mode));
//nChains = 1;
//fxt3veryloose
//fin.open("/star/u/yhleung2/pwg/29.KFParticle_hy2/taxi_hl_anatree_3gev_veryloose/filelist.txt");
//strcpy (ofile,Form("run18_3fxtgev_hl_train_%d_test.root",_cut_mode));
//nChains = 1;
//27gevveryloose
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_hyl_ana/filelist.txt");
//strcpy (ofile,Form("run18_27gev_hl_train_%d_test.root",_cut_mode));
//nChains = 1;

//rotate
//fin.open("/star/u/yhleung2/pwg/29.KFParticle_hy2/taxi_om_anatree_26gev_loose_rotate/filelist.txt");
//strcpy (ofile,Form("run18_26fxtgev_ht_train_%d_test_rotate.root",_cut_mode));
//nChains = 1;
//
//fin.open("/star/u/yhleung2/pwg/29.KFParticle_hy2/taxi_om_anatree_3gev_loose_rotate/filelist.txt");
//strcpy (ofile,Form("run18_3fxtgev_ht_train_%d_test_rotate.root",_cut_mode));
//nChains = 1;
//
//hl+ht
//fxt26
//fin.open("/star/u/yhleung2/pwg/29.KFParticle_hy2/taxi_hl_anatree_26gev_loose_rotate/filelist.txt");
//strcpy (ofile,Form("run18_26fxtgev_ht_train_%d_test_rotate.root",_cut_mode));
//nChains = 1;
//fxt3
//fin.open("/star/u/yhleung2/pwg/29.KFParticle_hy2/taxi_hl_anatree_3gev_loose_rotate/filelist.txt");
//strcpy (ofile,Form("run18_3fxtgev_hl_train_%d_test_rotate.root",_cut_mode));
//nChains = 1;
//fxt3veryloose
///fin.open("/star/u/yhleung2/pwg/29.KFParticle_hy2/taxi_hl_anatree_3gev_veryloose_rotate/filelist.txt");
//strcpy (ofile,Form("run18_3fxtgev_hl_train_%d_test_rotate.root",_cut_mode));
//nChains = 1;
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_hyl_ana_rotate/filelist.txt");
//strcpy (ofile,Form("run18_27gev_hl_train_%d_test_rotate.root",_cut_mode));
//nChains = 1;

fin.open("/star/u/yhleung2/pwg/33.standard_reco/filelist_pico_he3_mc.txt");
strcpy (ofile,"run18_27gev_he_hist.root");
nChains = 1;
}
if(_condor==1){
fin.open("filelist.txt");
strcpy (ofile,Form("ana_tree_cut_mode_%d.root",_cut_mode));
nChains = 1;
}

for(int i=0; i<nChains; i++){
      fin >> path_ch;
//      htriton_tree.Add(path_ch);
//      cout<<path_ch<<" "<<htriton_tree.GetEntries() <<endl;

//      h4lambda_tree.Add(path_ch);
//      cout<<path_ch<<" "<<h4lambda_tree.GetEntries() <<endl;

he3_tree.Add(path_ch);
cout<<path_ch<<" "<<he3_tree.GetEntries() <<endl;
he3_mc_tree.Add(path_ch);
cout<<path_ch<<" "<<he3_mc_tree.GetEntries() <<endl;

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

he3_tree.SetBranchAddress("brunid",&brunid);
he3_tree.SetBranchAddress("beventid",&beventid);
he3_tree.SetBranchAddress("brefmult",&brefmult);
he3_tree.SetBranchAddress("btofmult",&btofmult);
he3_tree.SetBranchAddress("bparticleid",&bparticleid);
he3_tree.SetBranchAddress("bparticlemass",&bparticlemass);
he3_tree.SetBranchAddress("bpx",&bpx);
he3_tree.SetBranchAddress("bpy",&bpy);
he3_tree.SetBranchAddress("bpz",&bpz);
he3_tree.SetBranchAddress("bdedx", &bdedx);
he3_tree.SetBranchAddress("bnhits", &bnhits);
he3_tree.SetBranchAddress("chi2primary_he", &chi2primary_he);
he3_tree.SetBranchAddress("reweight", &reweight);
he3_tree.SetBranchAddress("cent9", &cent9);
he3_tree.SetBranchAddress("bismc", &bismc);
he3_tree.SetBranchAddress("bdca", &bdca);

he3_mc_tree.SetBranchAddress("brunid", &brunid);
he3_mc_tree.SetBranchAddress("brefmult", &brefmult);
he3_mc_tree.SetBranchAddress("bmcparticleid", &bmcparticleid);
he3_mc_tree.SetBranchAddress("bmcrawpx", &bmcrawpx);
he3_mc_tree.SetBranchAddress("bmcrawpy", &bmcrawpy);
he3_mc_tree.SetBranchAddress("bmcrawpz", &bmcrawpz);
//he3_mc_tree.SetBranchAddress("refmultcor", &refmultcor);
he3_mc_tree.SetBranchAddress("reweight", &reweight);
he3_mc_tree.SetBranchAddress("cent9", &cent9);

Long64_t n_he3_Entries = he3_tree.GetEntries();
Long64_t n_he3_mc_Entries = he3_mc_tree.GetEntries();

 cout<<endl<<"Start processing "<<endl<<endl;
cout<<"--------------------------------------------"<<endl;

cout<<"Total entries "<<n_he3_Entries<< endl;

if(_cut_mode==0){//default
ht_chi2topo_cut = 5;
ht_chi2ndf_cut = 10;
ht_ldl_cut = 5;//default
ht_l_cut = 5;//defualt
chi2primary_he_cut = 18.6;
chi2primary_pi_cut = 18.6;

hl_chi2topo_cut = 5;
hl_chi2ndf_cut = 10;
hl_ldl_cut = 5;
hl_l_cut = 5;
chi2primary_h4_cut = 18.6;
chi2primary_pi_cut = 18.6;
}
if(_cut_mode==1){//loose
ht_chi2topo_cut = 5;
ht_chi2ndf_cut = 10;
ht_ldl_cut = 3;
ht_l_cut = 1;
chi2primary_he_cut = 18.6;
chi2primary_pi_cut = 18.6;

hl_chi2topo_cut = 5;
hl_chi2ndf_cut = 10;
hl_ldl_cut = 3;
hl_l_cut = 1;
chi2primary_h4_cut = 18.6;
chi2primary_pi_cut = 18.6;
}
if(_cut_mode>=3){//very loose,scan
ht_chi2topo_cut = 5;
ht_chi2ndf_cut = 10;
ht_ldl_cut = 3;//loose
ht_l_cut = 1;//loose

hl_chi2topo_cut = 5;
hl_chi2ndf_cut = 10;
hl_ldl_cut = 3;
hl_l_cut = 1;

if(_cut_mode==5){
chi2primary_he_cut = 5.;
chi2primary_pi_cut = 5.;
chi2primary_h4_cut = 5.;
}
if(_cut_mode==4){
chi2primary_he_cut = 10.;
chi2primary_pi_cut = 10.;
chi2primary_h4_cut = 10.;
}
if(_cut_mode==3){
chi2primary_he_cut = 15.;
chi2primary_pi_cut = 15.;
chi2primary_h4_cut = 15.;
}


}

//n_he3_Entries = 100000;
for (Long64_t iEntry = 0; iEntry<=n_he3_Entries; iEntry++){

	he3_tree.GetEntry(iEntry);
  if(iEntry%100000==0){
    cout<<"Processing entry:" <<iEntry<<endl;
}
int centtmva_label = -1;
int centbdt_label = -1;
int cent_label = -1;
int pt_label = -1;


//if(notbadrun>0) continue;


if(snn==28){//temp
    if(cent9==1) cent_label = 7;
    if(cent9==2) cent_label = 6;
    if(cent9==3) cent_label = 5;
    if(cent9==4) cent_label = 4;
    if(cent9==5) cent_label = 3;
    if(cent9==6) cent_label = 2;
    if(cent9==7) cent_label = 1;
    if(cent9==8) cent_label = 0;

      if(cent9==0) centbdt_label = 2;
      if(cent9==1) centbdt_label = 2;
      if(cent9==2) centbdt_label = 2;
      if(cent9==3) centbdt_label = 2;
      if(cent9==4) centbdt_label = 1;
      if(cent9==5) centbdt_label = 1;
      if(cent9==6) centbdt_label = 1;
      if(cent9==7) centbdt_label = 0;
      if(cent9==8) centbdt_label = 0;


      if(cent9==0) centtmva_label = 1;
      if(cent9==1) centtmva_label = 1;
      if(cent9==2) centtmva_label = 1;
      if(cent9==3) centtmva_label = 1;
      if(cent9==4) centtmva_label = 1;
      if(cent9==5) centtmva_label = 1;
      if(cent9==6) centtmva_label = 1;
      if(cent9==7) centtmva_label = 0;
      if(cent9==8) centtmva_label = 0;

if(_cent_off){
cent_label = 0;
centbdt_label = 0;
centtmva_label = 0;
}

}
  if(cent_label<0) continue;//<=5refmult

  ycm = -1.045;


  bpt = sqrt(bpx*bpx+bpy*bpy);
//  bp = 0.5*sqrt(bpx*bpx+bpy*bpy+bpz*bpz);
   bp = sqrt(bpx*bpx+bpy*bpy+bpz*bpz);
   _rotate = 0;


  TLorentzVector ptc(bpx,bpy,bpz,sqrt(bpx*bpx+bpy*bpy+bpz*bpz+bparticlemass*bparticlemass));
  beta = ptc.Eta(); 
  brap = ptc.Rapidity() - ycm;

//  _applyweight = false;
//  if(brefmult<=60)  _applyweight = true;


if(bparticleid>0 && bismc==1){

h_dedx_p->Fill(bp,bdedx);
h_chi2primary_he->Fill(beta,bpt,chi2primary_he);
h_dca->Fill(beta,bpt,bdca);
if(chi2primary_he<18.6){
h_bnhits->Fill(beta,bpt,bnhits);
h_bdedx->Fill(beta,bpt,bdedx);
}
h_pt->Fill(brap,bpt);

}//cut mode 0

}//entry loop

for (Long64_t iEntry = 0; iEntry<=n_he3_mc_Entries; iEntry++){

                 he3_mc_tree.GetEntry(iEntry);
           if(iEntry%100000==0){
             cout<<"Processing mc entry:" <<iEntry<<endl;
         }

         bpt = sqrt(bmcrawpx*bmcrawpx+bmcrawpy*bmcrawpy);

         TLorentzVector ptc(bmcrawpx,bmcrawpy,bmcrawpz,sqrt(bmcrawpx*bmcrawpx+bmcrawpy*bmcrawpy+bmcrawpz*bmcrawpz+2.8094135*2.8094135));
         beta = ptc.Eta();
         brap = ptc.Rapidity() - ycm;
         g_pt->Fill(brap,bpt);

         }

write_histo();

}

void book_histo()
{

  h_dedx_p = new TH2F("h_dedx_p","",500,0,5,2000,0,100);
  h_dedx_p->Sumw2();


h_chi2primary_he = new TH3F("h_chi2primary_he","",10,-2.0,0.0,50,0,5,500,0,500);
h_bnhits = new TH3F("h_bnhits","",10,-2.0,0.0,50,0,5,70,0,70);
h_bdedx = new TH3F("h_bdedx","",10,-2.0,0.0,50,0,5,100,0,100);
h_pt = new TH2F("h_pt","",100,-1.0,1.0,100,0,10);
h_dca = new TH3F("h_dca","",10,-2.0,0.0,50,0,5,500,0,20);
g_pt = new TH2F("g_pt","",100,-1.0,1.0,100,0,10);
h_chi2primary_he ->Sumw2();
h_bnhits ->Sumw2();
h_bdedx ->Sumw2();
h_dca ->Sumw2();
h_pt->Sumw2();
g_pt->Sumw2();

}

void write_histo()
{
    TFile *outhistfile = new TFile (ofile, "RECREATE");
    outhistfile->cd();

	h_dedx_p->Write();
h_chi2primary_he->Write();
h_bnhits->Write();
h_bdedx->Write();
h_pt->Write();
g_pt->Write();
h_dca ->Write();
}
