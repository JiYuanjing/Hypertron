//THIS MACRO REQUIRES TMVA package, better run on rcf
//write mass/pt histograms that can be read offline

//v6 add countrefmult
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPaveLabel.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TNtuple.h"
#include "TText.h"
#include "TDatime.h"
#include "TRandom.h"
#include "TChain.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TRandom2.h"
#include "TVector3.h"
#include "TGraphErrors.h"

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
 TH2F *h_eta;
 TH2F *g_eta;
//TH2F *mkfIncMassPtUL3_0[ncenttmvabins];
//TH2F *mkfIncMassPtUL3_1[ncenttmvabins];
//TH2F *mkfIncMassPtUL4_0[ncenttmvabins];
//TH2F *mkfIncMassPtUL4_1[ncenttmvabins];

TH3F *h_nhits_proton;
TH3F *h_nhits_pi;

TH2F *mkfIncMassPtUL3_0[ncentbins];
TH2F *mkfIncMassPtUL3_1[ncentbins];
TH2F *mkfIncMassPtUL4_0[ncentbins];
TH2F *mkfIncMassPtUL4_1[ncentbins];

TH2F *mkfIncMassPtUL5_0[ncentbins];
TH2F *mkfIncMassPtUL6_0[ncentbins];
TH3F *h_mass_qa;
TH3F *h_mass_pl;
TH3F *h_mass_pt;
TH3F *h_chi2primary_proton;
TH3F *h_chi2primary_pi;
TH3F *h_ld_ldl;
TH3F *h_ld_l;
TH3F *h_ld_dl;
TH3F *h_ld_chi2topo;
TH3F *h_ld_chi2ndf;
TH3F *h_ld_dca_proton;
TH3F *h_ld_dca_pi;

TH3F *h_mass_pt_wgt;
TH3F *h_mass_pl_wgt;

TH3F *h_bnhits;
TH3F *h_bdedx;

TH3F *h_ld_bdfvtx;
TH3F *h_ld_bdfvtx2;
TH3F *h_ld_lifetime;
TH2F *h_lt_mc;
TH2F *g_lt_mc;

TH2F *h_pt;
TH2F *g_pt;
TH2F *g_pl;
TH2F *g_pt_wgt;
TH2F *g_pl_wgt;
TH3F *h_dca;

TH2F *h_p_mc;
TH2F *h_pt_mc;
TH2F *h_l_mc;
TH2F *h_pl_mc;
TH2F *g_p_mc;
TH2F *g_pt_mc;
TH2F *g_l_mc;
TH2F *g_pl_mc;
int cent_label;
bool _applyweight;

// double ptbin[npttmvabins+1] = {0.0,1.0,2.0,3.0,5.0,10.0};
//double ptbin[nptbins+1] = {0.4,0.8,1.2,1.6,2.0,2.4,2.8};
double ptbin[nptbins+1] = {0.4,0.8,1.2,1.6,2.0,2.4,2.8,3.2,3.6};

int centtmva_limit[ncenttmvabins+1] = {0,10,80};

int cent_limit[ncentbins+1] = {0,5,10,20,30,40,50,60,70,80};
//int cent_limit[ncentbins+1] = {0,5,10,20,30,40,50,60,80};
bool _cent_off = true;
//
//int snn=28;
int snn=3;

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
float chi2Primary_Pion;
float Chi2NDF;
float LdL;
float L;
float Chi2Topo;
float Chi2NDF_ld;
float LdL_ld;
float L_ld;
float Chi2Topo_ld;

int nhits_ld_pi;
int nhits_ld_proton;

float ht_chi2topo_cut;
float ht_chi2ndf_cut;
float ht_ldl_cut;
float ht_l_cut;
float chi2primary_he_cut;

float hl_chi2topo_cut;
float hl_chi2ndf_cut;
float hl_ldl_cut;
float hl_l_cut;
float chi2primary_h4_cut;

double ht_mass = 2.99131;
double ht_width = 0.005;
double hl_mass = 3.9239;
double hl_width = 0.005;

float ld_mass = 1.115683;

int mccountrefmult;
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

float mass_lo;
float mass_hi;

float bmcrawpl;
float hl_chi2ndf;
float hl_ldl;
float hl_chi2topo;
float hl_dl;
float hl_l;
float chi2primary_h4;
float bmcrap;
float bmcpx;
float bmcpy;
float bmcpz;
float bmcl;
float bmcpl;
float bpl;

float ld_bdfvtx;
float ld_bdfvtx2;
float ld_lifetime;

float dca_proton;
float dca_pi;

TF1 *levyfit4;
TF1 *levyfit5;
TF1 *levyfit6;
TF1 *levyfit7;
TF1 *levyfit8;

int bismc;
//TMVA::Reader* reader[ncenttmvabins][npttmvabins];
// float bdt_ld_cut[ncenttmvabins][npttmvabins];

//obsolete, used for bdt scanning
float bdt_ld_cut[ncenttmvabins][nptbins];

//used for om analysis
float bdt_om_cut[ncentbins][nptbins];
float bdt_ob_cut[ncentbins][nptbins];

float weight;
int countrefmult;

int centrality;
int centFull[ncentbins] ={15,22,32,43,57,73,92,117,133};

//using for xi analyusis
const int ncentbdtbins = 3;
const int nptbdtbins = 7;
float bdt_xi_cut[ncentbdtbins][nptbdtbins];
double ptbdtbin[nptbdtbins+1] = {0.2,0.6,1.0,1.4,1.8,2.2,2.6,3.0};
int choose_cent;

int Centrality(int aa );
void book_histo();
void write_histo();

void process_tree_tmva_ld_kf_vmv6(int condor=0, int cut_mode=0){

    _condor = condor;//condor=1: condor settings
    _cut_mode = cut_mode;

    choose_cent	= 0;
//    choose_cent = 1;
//    choose_cent = 2;


levyfit4 = new TF1("levyfit4","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.4,2.4);
levyfit5 = new TF1("levyfit5","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.2,2.4);
levyfit6 = new TF1("levyfit6","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.2,2.4);
levyfit7 = new TF1("levyfit7","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.2,2.4);
levyfit8 = new TF1("levyfit8","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.2,2.4);

levyfit4->SetParameters(1.53536e-01, -1.21064e+08, 5.25287e+07);
levyfit5->SetParameters(1.51332e-01 , -1.25927e+08, 4.86193e+07);
levyfit6->SetParameters(1.45903e-01 , -1.60954e+08, 4.74115e+07);
levyfit7->SetParameters(1.28382e-01 , -1.60610e+08, 4.92025e+07);
levyfit8->SetParameters(1.11371e-01 , -2.04689e+08, 4.22255e+07);
//

 // This loads the library
//     TMVA::Tools::Instance();

//CHANGE NAME TO MATCH READER>>>
/*
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
*/
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
TChain lambda_mc_tree("lambda_mc_tree");

int nChains = -999999;
ifstream fin;

if(_condor==0){//emb
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ld_mu/sum/filelist.txt");
//strcpy (ofile,"run18_3gev_ld_hist_mu.root");
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/filelist_ld_emb3.txt");
//strcpy (ofile,"run18_3gev_ld_hist_mu3.root");
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ld_mu3_correctedvertex/sum/filelist.txt");
//strcpy (ofile,"run18_3gev_ld_hist_mu4.root");

//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ld_mu3_correctedvertex_anderror/filelist.txt");
//strcpy (ofile,"run18_3gev_ld_hist_mu5.root");
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ld_mu3_correctedvertex_andnoerror/filelist.txt");
//strcpy (ofile,"run18_3gev_ld_hist_mu5noerr.root");

//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ld_mu3_correctedvertex_werror2/filelist.txt");
//strcpy (ofile, "run18_3gev_ld_hist_mu6.root");
//nChains = 1;
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ld_mu3_correctedvertex_werror4/filelist.txt");
//strcpy (ofile, "run18_3gev_ld_hist_mu7.root");

//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ld_mu3_correctedvertex_werror_veryloose/filelist.txt");
//strcpy (ofile, "run18_3gev_ld_hist_mu_veryloose.root");
//nChains = 1;

//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ld_mu3_correctedvertex_werror_3/sum/filelist.txt");
//strcpy (ofile, "run18_3gev_ld_hist_mu_3.root");
//nChains = 1;

//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ld_mu3_correctedvertex_w10xerror_3/sum/filelist.txt");
//strcpy (ofile, "run18_3gev_ld_hist_mu_10x_3.root");
//nChains = 1;

fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ld_mu3_correctedvertex_w5xerror_4/filelist.txt");
strcpy (ofile, Form("run18_3gev_ld_hist_mu_5x_3_newreweight_centa%d.root",choose_cent));
nChains = 1;
}
if(_condor==1){//data
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ld_pico/sum/filelist.txt");
//strcpy (ofile,"run18_3gev_ld_hist_pico.root");
//
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ld_pico_rotate/sum/filelist.txt");
//strcpy (ofile,"run18_3gev_ld_hist_pico_rotate.root");
//nChains = 1;
//
//loose cuts
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ld_pico_tree/filelist.txt");
//strcpy (ofile,"run18_3gev_ld_hist_pico_loose.root");
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ld_pico_rotate_tree/filelist.txt");
//strcpy (ofile,"run18_3gev_ld_hist_pico_rotate_loose.root");
//nChains = 1200;

//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ld_pico_tree_2/filelist.txt");
//strcpy (ofile,"run18_3gev_ld_hist_pico_loose_7.root");
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ld_pico_rotate_tree_2/filelist.txt");
//strcpy (ofile,"run18_3gev_ld_hist_pico_rotate_loose_7.root");

//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ld_pico_tree_4/filelist.txt");
//strcpy (ofile,Form("run18_3gev_ld_hist_pico_loose_4_centa%d.root",choose_cent));
fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ld_pico_rotate_tree_4/filelist.txt");
strcpy (ofile,Form("run18_3gev_ld_hist_pico_rotate_loose_4_centa%d.root",choose_cent));
nChains = 1200;
}
if(_condor==2){//debug mode
fin.open("filelist.txt");
strcpy (ofile,"testrotate.root");
nChains = 1;
}


for(int i=0; i<nChains; i++){
      fin >> path_ch;
//      htriton_tree.Add(path_ch);
//      cout<<path_ch<<" "<<htriton_tree.GetEntries() <<endl;

//      h4lambda_tree.Add(path_ch);
//      cout<<path_ch<<" "<<h4lambda_tree.GetEntries() <<endl;

lambda_tree.Add(path_ch);
cout<<path_ch<<" "<<lambda_tree.GetEntries() <<endl;
lambda_mc_tree.Add(path_ch);
cout<<path_ch<<" "<<lambda_mc_tree.GetEntries() <<endl;


//the folllwoing needs to be reimplemented
/*
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
*/
}//ichain
cout << " done TChain-ing.. " << endl;

book_histo();

//lambda_tree.SetBranchAddress("brunid",&brunid);
lambda_tree.SetBranchAddress("brefmult",&brefmult);
lambda_tree.SetBranchAddress("btofmult",&btofmult);
lambda_tree.SetBranchAddress("bparticleid",&bparticleid);
lambda_tree.SetBranchAddress("bparticlemass",&bparticlemass);
lambda_tree.SetBranchAddress("bpx",&bpx);
lambda_tree.SetBranchAddress("bpy",&bpy);
lambda_tree.SetBranchAddress("bpz",&bpz);
lambda_tree.SetBranchAddress("chi2primary_proton", &chi2primary_proton);
lambda_tree.SetBranchAddress("chi2primary_pi", &chi2primary_pi);
lambda_tree.SetBranchAddress("ld_chi2topo", &ld_chi2topo);
lambda_tree.SetBranchAddress("ld_chi2ndf", &ld_chi2ndf);
lambda_tree.SetBranchAddress("ld_ldl", &ld_ldl);
lambda_tree.SetBranchAddress("ld_l", &ld_l);
lambda_tree.SetBranchAddress("ld_dl", &ld_dl);
lambda_tree.SetBranchAddress("dca_proton",&dca_proton);
lambda_tree.SetBranchAddress("dca_pi",&dca_pi);
lambda_tree.SetBranchAddress("nhits_ld_pi",&nhits_ld_pi);
lambda_tree.SetBranchAddress("nhits_ld_proton",&nhits_ld_proton);
lambda_tree.SetBranchAddress("ld_bdfvtx",&ld_bdfvtx);
lambda_tree.SetBranchAddress("ld_bdfvtx2",&ld_bdfvtx2);
lambda_tree.SetBranchAddress("ld_lifetime",&ld_lifetime);
lambda_tree.SetBranchAddress("countrefmult",&countrefmult);
lambda_tree.SetBranchAddress("reweight", &reweight);
lambda_tree.SetBranchAddress("cent9", &cent9);
lambda_tree.SetBranchAddress("bismc", &bismc);

if(_condor==0){
lambda_tree.SetBranchAddress("bmcpx",&bmcpx);
lambda_tree.SetBranchAddress("bmcpy",&bmcpy);
lambda_tree.SetBranchAddress("bmcpz",&bmcpz);
lambda_tree.SetBranchAddress("bmcl",&bmcl);
lambda_tree.SetBranchAddress("bmcpl",&bmcpl);
lambda_tree.SetBranchAddress("bpl",&bpl);
}
//lambda_mc_tree.SetBranchAddress("brunid", &brunid);
lambda_mc_tree.SetBranchAddress("brefmult", &brefmult);
lambda_mc_tree.SetBranchAddress("bmcparticleid", &bmcparticleid);
lambda_mc_tree.SetBranchAddress("bmcrawpx", &bmcrawpx);
lambda_mc_tree.SetBranchAddress("bmcrawpy", &bmcrawpy);
lambda_mc_tree.SetBranchAddress("bmcrawpz", &bmcrawpz);
lambda_mc_tree.SetBranchAddress("bmcrawpl", &bmcrawpl);
lambda_mc_tree.SetBranchAddress("reweight", &reweight);
lambda_mc_tree.SetBranchAddress("countrefmult", &mccountrefmult);
lambda_mc_tree.SetBranchAddress("cent9", &cent9);

Long64_t n_lambda_Entries = lambda_tree.GetEntries();
Long64_t n_lambda_mc_Entries = lambda_mc_tree.GetEntries();

 cout<<endl<<"Start processing "<<endl<<endl;
cout<<"--------------------------------------------"<<endl;

cout<<"Total entries "<<n_lambda_Entries<< endl;

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

//n_lambda_Entries = 100000;
for (Long64_t iEntry = 0; iEntry<=n_lambda_Entries; iEntry++){

	lambda_tree.GetEntry(iEntry);
  if(iEntry%100000==0){
    cout<<"Processing entry:" <<iEntry<<endl;
}
int centtmva_label = -1;
int centbdt_label = -1;
 cent_label = -1;
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
}
if(snn==3){
    cent_label = Centrality(countrefmult);
}
//if(_cent_off){
//cent_label = 0;
//centbdt_label = 0;
//centtmva_label = 0;
//}


  if(cent_label<0) continue;//<=5refmult

   if(cent_label<7 && choose_cent==1) continue;//cent7,8
   if( (cent_label>6 || cent_label<2) && choose_cent==2) continue;//cent2,3,4,5,6

if(snn==3){
  ycm = -1.045;
}else{
  ycm = -999.;
}


 mass_lo = 1.105;
 mass_hi = 1.125;


  bpt = sqrt(bpx*bpx+bpy*bpy);
//  bp = 0.5*sqrt(bpx*bpx+bpy*bpy+bpz*bpz);
   bp = sqrt(bpx*bpx+bpy*bpy+bpz*bpz);
   _rotate = 0;


  TLorentzVector ptc(bpx,bpy,bpz,sqrt(bpx*bpx+bpy*bpy+bpz*bpz+bparticlemass*bparticlemass));
  beta = ptc.Eta(); 
  brap = ptc.Rapidity() - ycm;

// TLorentzVector mcptc(bmcpx,bmcpy,bmcpz,sqrt(bmcpx*bmcpx+bmcpy*bmcpy+bmcpz*bmcpz+ld_mass*ld_mass));	
// bmcrap = mcptc.Rapidity() - ycm; 

double betaxgamma;
betaxgamma = sqrt(bpx*bpx+bpy*bpy+bpz*bpz)/bparticlemass;
double ld_pl;
ld_pl = ld_l/betaxgamma;

//  _applyweight = false;
//  if(brefmult<=60)  _applyweight = true;

//cout<<"bparticleid:"<<bparticleid<<" "<<bismc<<" "<<_condor<<endl;

if(bparticleid>0 && ((bismc==1&&_condor==0) || _condor==1 || _condor==2) ){

if(bparticlemass>mass_lo && bparticlemass<mass_hi){

h_chi2primary_proton->Fill(beta,bpt,chi2primary_proton);
h_chi2primary_pi->Fill(beta,bpt,chi2primary_pi);
h_nhits_proton->Fill(beta,bpt,nhits_ld_proton);
h_nhits_pi->Fill(beta,bpt,nhits_ld_pi);
h_ld_ldl->Fill(beta,bpt,ld_ldl);
h_ld_l->Fill(beta,bpt,ld_l);
h_ld_dl->Fill(beta,bpt,ld_dl);
h_ld_chi2topo->Fill(beta,bpt,ld_chi2topo);
h_ld_chi2ndf->Fill(beta,bpt,ld_chi2ndf);
h_ld_dca_proton->Fill(beta,bpt,dca_proton);
h_ld_dca_pi->Fill(beta,bpt,dca_pi);
h_pt->Fill(brap,bpt);
h_eta->Fill(beta,bpt);
h_ld_bdfvtx->Fill(beta,bpt,ld_bdfvtx);
h_ld_bdfvtx2->Fill(beta,bpt,ld_bdfvtx2);
h_ld_lifetime->Fill(beta,bpt,ld_lifetime);

}
h_mass_qa->Fill(beta,bpt,bparticlemass);

h_mass_pt->Fill(bparticlemass,brap,bpt);
h_mass_pl->Fill(bparticlemass,brap,ld_pl);

if(bismc==1&&_condor==0){

TLorentzVector mcptc(bmcpx,bmcpy,bmcpz,sqrt(bmcpx*bmcpx+bmcpy*bmcpy+bmcpz*bmcpz+ld_mass*ld_mass));
 bmcrap = mcptc.Rapidity() - ycm;

if(fabs(bmcrap)<0.1){ weight = levyfit4->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy)); } 
else if(fabs(bmcrap)>=0.1 && fabs(bmcrap)<0.3){ weight = levyfit5->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy)); } 
else if(fabs(bmcrap)>=0.3 && fabs(bmcrap)<0.5){ weight = levyfit6->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy)); } 
else if(fabs(bmcrap)>=0.5 && fabs(bmcrap)<0.7){ weight = levyfit7->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy)); }
else{ weight = levyfit8->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy)); }

h_mass_pt_wgt->Fill(bparticlemass,brap,bpt,weight);
h_mass_pl_wgt->Fill(bparticlemass,brap,ld_pl,weight);
}


//h_dedx_p->Fill(bp,bdedx);
//h_dca->Fill(beta,bpt,bdca);
//if(chi2primary_he<18.6){
//h_bnhits->Fill(beta,bpt,bnhits);
//h_bdedx->Fill(beta,bpt,bdedx);
//}

if(_condor==0){
h_p_mc->Fill(bp,sqrt(bmcpx*bmcpx+bmcpy*bmcpy+bmcpz*bmcpz));
h_pt_mc->Fill(bpt,sqrt(bmcpx*bmcpx+bmcpy*bmcpy));
h_l_mc->Fill(ld_l,bmcl);
h_pl_mc->Fill(ld_pl,bmcpl);
h_lt_mc->Fill(ld_lifetime,bmcpl);

g_p_mc->Fill(bp,bp-sqrt(bmcpx*bmcpx+bmcpy*bmcpy+bmcpz*bmcpz));
g_pt_mc->Fill(bpt,bpt-sqrt(bmcpx*bmcpx+bmcpy*bmcpy));
g_l_mc->Fill(ld_l,ld_l-bmcl);
g_pl_mc->Fill(ld_pl,ld_pl-bmcpl);
g_lt_mc->Fill(ld_lifetime,ld_lifetime-bmcpl);
}


}//bparticleid loop

}//entry loop

for (Long64_t iEntry = 0; iEntry<=n_lambda_mc_Entries; iEntry++){

                 lambda_mc_tree.GetEntry(iEntry);
           if(iEntry%100000==0){
             cout<<"Processing mc entry:" <<iEntry<<endl;
         }

cent_label=-1;
if(snn==3){
      cent_label = Centrality(mccountrefmult);
  }
if(cent_label<0) continue;//<=5refmult
     if(cent_label<7 && choose_cent==1) continue;//cent7,8,9
     if( (cent_label>6 || cent_label<2) && choose_cent==2) continue;//cent2,3,4,5,6

         bpt = sqrt(bmcrawpx*bmcrawpx+bmcrawpy*bmcrawpy);

         TLorentzVector ptc(bmcrawpx,bmcrawpy,bmcrawpz,sqrt(bmcrawpx*bmcrawpx+bmcrawpy*bmcrawpy+bmcrawpz*bmcrawpz+ld_mass*ld_mass));
         beta = ptc.Eta();
         brap = ptc.Rapidity() - ycm;

double betaxgamma;
betaxgamma = sqrt(bmcrawpx*bmcrawpx+bmcrawpy*bmcrawpy+bmcrawpz*bmcrawpz)/ld_mass;
double ld_pl;
ld_pl = bmcrawpl/betaxgamma;

         g_pt->Fill(brap,bpt);
         g_eta->Fill(beta, bpt);
	 g_pl->Fill(brap,ld_pl);

if(fabs(brap)<0.1){ weight = levyfit4->Eval(bpt); }
else if(fabs(brap)>=0.1 && fabs(brap)<0.3){ weight = levyfit5->Eval(bpt); }
else if(fabs(brap)>=0.3 && fabs(brap)<0.5){ weight = levyfit6->Eval(bpt); }
else if(fabs(brap)>=0.5 && fabs(brap)<0.7){ weight = levyfit7->Eval(bpt); }
else{ weight = levyfit8->Eval(bpt); }
//weight = levyfit7->Eval(bpt);

g_pt_wgt->Fill(brap,bpt,weight);
g_pl_wgt->Fill(brap,ld_pl,weight);

         }

write_histo();

}

void book_histo()
{

  h_dedx_p = new TH2F("h_dedx_p","",500,0,5,2000,0,100);
  h_dedx_p->Sumw2();

h_mass_pt = new TH3F("h_mass_pt","",1000,1,2,9,-0.9,0.9,50,0,5);
h_mass_pl = new TH3F("h_mass_pl","",1000,1,2,9,-0.9,0.9,500,0,100);

h_mass_pt_wgt = new TH3F("h_mass_pt_wgt","",1000,1,2,9,-0.9,0.9,50,0,5);
h_mass_pl_wgt = new TH3F("h_mass_pl_wgt","",1000,1,2,9,-0.9,0.9,500,0,100);

h_mass_qa = new TH3F("h_mass_qa","",10,-2.0,0.0,50,0,5,1000,1,2);
h_chi2primary_proton = new TH3F("h_chi2primary_proton","",10,-2.0,0.0,50,0,5,500,0,500);
h_chi2primary_pi = new TH3F("h_chi2primary_pi","",10,-2.0,0.0,50,0,5,500,0,500);
h_ld_ldl = new TH3F("h_ld_ldl","",10,-2.0,0.0,50,0,5,1000,0,200);
h_ld_l = new TH3F("h_ld_l","",10,-2.0,0.0,50,0,5,1000,0,200);
h_ld_dl = new TH3F("h_ld_dl","",10,-2.0,0.0,50,0,5,500,0,20);
h_ld_chi2topo = new TH3F("h_ld_chi2topo","",10,-2.0,0.0,50,0,5,500,0,10);
h_ld_chi2ndf = new TH3F("h_ld_chi2ndf","",10,-2.0,0.0,50,0,5,100,0,10);
h_nhits_proton = new TH3F("h_nhits_proton","",10,-2.0,0.0,50,0,5,80,0,80);
h_nhits_pi = new TH3F("h_nhits_pi","",10,-2.0,0.0,50,0,5,80,0,80);

h_ld_dca_proton = new TH3F("h_ld_dca_proton","",10,-2.0,0.0,50,0,5,200,0,20);
h_ld_dca_pi     = new TH3F("h_ld_dca_pi","",10,-2.0,0.0,50,0,5,200,0,20);

h_ld_bdfvtx = new TH3F("h_ld_bdfvtx","",10,-2.0,0.0,50,0,5,1000,0,5);
h_ld_bdfvtx2 = new TH3F("h_ld_bdfvtx2","",10,-2.0,0.0,50,0,5,1000,0,1);
h_ld_lifetime = new TH3F("h_ld_lifetime","",10,-2.0,0.0,50,0,5,500,0,100);

h_bnhits = new TH3F("h_bnhits","",10,-2.0,0.0,50,0,5,70,0,70);
h_bdedx = new TH3F("h_bdedx","",10,-2.0,0.0,50,0,5,100,0,100);
h_pt = new TH2F("h_pt","",100,-1.0,1.0,100,0,10);
h_eta = new TH2F("h_eta","",200,-2.5,0.5,100,0,10);
h_dca = new TH3F("h_dca","",10,-2.0,0.0,50,0,5,500,0,20);
g_pt = new TH2F("g_pt","",9,-0.9,0.9,50,0,5);
//g_pt = new TH2F("g_pt","",100,-1,1,50,0,5);
g_pl = new TH2F("g_pl","",9,-0.9,0.9,500,0,100);

g_pt_wgt = new TH2F("g_pt_wgt","",9,-0.9,0.9,50,0,5);
g_pl_wgt = new TH2F("g_pl_wgt","",9,-0.9,0.9,500,0,100);

g_eta = new TH2F("g_eta","",50,-2.5,0.5,50,0,5);

h_p_mc = new TH2F("h_p_mc","",1000,0,10,1000,0,10);
h_pt_mc = new TH2F("h_pt_mc","",1000,0,10,1000,0,10);
h_l_mc = new TH2F("h_l_mc","",1000,0,100,1000,0,100);
h_pl_mc = new TH2F("h_pl_mc","",1000,0,100,1000,0,100);
h_lt_mc = new TH2F("h_lt_mc","",1000,0,100,1000,0,100);

g_p_mc = new TH2F("g_p_mc","",1000,0,10,1000,-10,10);
g_pt_mc = new TH2F("g_pt_mc","",1000,0,10,1000,-10,10);
g_l_mc = new TH2F("g_l_mc","",1000,0,100,1000,-100,100);
g_pl_mc = new TH2F("g_pl_mc","",1000,0,100,1000,-100,100);
g_lt_mc = new TH2F("g_lt_mc","",1000,0,100,1000,-100,100);

g_pt_wgt->Sumw2();
g_pl_wgt->Sumw2();
h_mass_pt_wgt->Sumw2();
h_mass_pl_wgt->Sumw2();

h_ld_bdfvtx->Sumw2();
h_ld_bdfvtx2->Sumw2();
h_ld_lifetime->Sumw2();
h_lt_mc->Sumw2();
g_lt_mc->Sumw2();

h_ld_dca_proton->Sumw2();
h_ld_dca_pi->Sumw2();
h_mass_qa->Sumw2();
h_chi2primary_proton ->Sumw2();
h_bnhits ->Sumw2();
h_bdedx ->Sumw2();
h_dca ->Sumw2();
h_pt->Sumw2();
g_pt->Sumw2();
g_eta->Sumw2();
h_mass_pl->Sumw2();

h_nhits_proton ->Sumw2();
h_nhits_pi ->Sumw2();

h_p_mc->Sumw2();
h_pt_mc->Sumw2();
h_l_mc->Sumw2();
h_pl_mc->Sumw2();

g_p_mc->Sumw2();
g_pt_mc->Sumw2();
g_l_mc->Sumw2();
g_pl_mc->Sumw2();
}

void write_histo()
{
    TFile *outhistfile = new TFile (ofile, "RECREATE");
    outhistfile->cd();

h_mass_qa->Write();
h_mass_pt->Write();
h_mass_pl->Write();
h_chi2primary_proton->Write();
h_chi2primary_pi->Write();
h_ld_ldl->Write();
h_ld_l->Write();
h_ld_dl->Write();
h_ld_chi2topo->Write();
h_ld_chi2ndf->Write();
h_ld_dca_proton->Write();
h_ld_dca_pi->Write();
h_nhits_proton->Write();
h_nhits_pi ->Write();
h_eta->Write();
h_dedx_p->Write();
h_bnhits->Write();
h_bdedx->Write();
h_pt->Write();
g_pt->Write();
g_pl->Write();
h_dca ->Write();
g_eta->Write();
h_ld_bdfvtx->Write();
h_ld_bdfvtx2->Write();
h_ld_lifetime->Write();
g_pt_wgt->Write();
g_pl_wgt->Write();
h_mass_pt_wgt->Write();
h_mass_pl_wgt->Write();


if(_condor==0){
h_p_mc -> Write();
h_pt_mc -> Write();
h_l_mc -> Write();
h_pl_mc -> Write();
g_p_mc -> Write();
g_pt_mc -> Write();
g_l_mc -> Write();
g_pl_mc -> Write();
h_lt_mc -> Write();
g_lt_mc -> Write();
}


}
int Centrality(int gRefMult )
{
    if      (gRefMult>=centFull[8]) centrality=8;
    else if (gRefMult>=centFull[7]) centrality=7;
    else if (gRefMult>=centFull[6]) centrality=6;
    else if (gRefMult>=centFull[5]) centrality=5;
    else if (gRefMult>=centFull[4]) centrality=4;
    else if (gRefMult>=centFull[3]) centrality=3;
    else if (gRefMult>=centFull[2]) centrality=2;
    else if (gRefMult>=centFull[1]) centrality=1;
    else if (gRefMult>=centFull[0]) centrality=0;
    else centrality = 9;

    return centrality;
}

