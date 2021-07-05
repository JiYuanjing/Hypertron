//v8: cut scan
//v7: add fiducial
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

 TH2F *h_eta;
 TH2F *h_rap;
 TH2F *g_eta;
//TH2F *mkfIncMassPtUL3_0[ncenttmvabins];
//TH2F *mkfIncMassPtUL3_1[ncenttmvabins];
//TH2F *mkfIncMassPtUL4_0[ncenttmvabins];
//TH2F *mkfIncMassPtUL4_1[ncenttmvabins];

TH3F *h_ht_bdt;
TH3F *h_nhits_he;
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
TH3F *h_chi2primary_he;
TH3F *h_chi2primary_pi;
TH3F *h_ht_ldl;
TH3F *h_ht_l;
TH3F *h_ht_dl;
TH3F *h_ht_chi2topo;
TH3F *h_ht_chi2ndf;
TH3F *h_ht_dca_he;
TH3F *h_ht_dca_pi;

TH3F *h_mass_pt_wgt;
TH3F *h_mass_pl_wgt;
TH1F *h_cent;

TH3F *h_bnhits;
TH3F *h_bdedx;
TH2F *g_pt_fine;
TH2F *g_pt_fine_wgt;
TH2F *g_pt_fine2;
TH2F *g_pt_fine2_wgt;

TH3F *h_ht_bdfvtx;
TH3F *h_ht_bdfvtx2;
TH3F *h_ht_lifetime;
TH2F *h_lt_mc;
TH2F *g_lt_mc;

TH2F *h_pt;
TH2F *h_r_pt;
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

TH2F *g_pt_fine_in;

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
float xi_ht_chi2topo  ;
float xi_ht_chi2ndf   ;
float xi_ht_ldl       ;
float xi_ht_l         ;
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
float ht_pv_dca_ld;
float p_pv_dca_ld;
float pi_pv_dca_ld;
float pi_p_dca_ld;
float decaylen_ld;
float ang0_ld;
*/

float ht_chi2topo     ;
float ht_chi2ndf      ;
float chi2primary_pi     ;

float chi2Primary_Pi;
float chi2Primary_Proton;

float om_chi2topo;
float om_chi2ndf;
float om_ldl;
float om_ht_chi2topo;
float om_ht_chi2ndf;
float om_ht_ldl;
float om_ht_l;
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

float dca_he_cut;
float om_chi2topo_cut;
float om_chi2ndf_cut;
float om_ldl_cut;
float om_ht_chi2topo_cut;
float om_ht_chi2ndf_cut;
float om_ht_ldl_cut;
float om_ht_l_cut;
float om_l_cut;
float om_dl_cut;
float chi2primary_om_proton_cut;
float chi2primary_om_pi_cut;
float chi2primary_om_bach_cut;
float chi2primary_om_ht_cut;

float xi_chi2topo_cut;
float xi_chi2ndf_cut;
float xi_ldl_cut;
float xi_ht_chi2topo_cut;
float xi_ht_chi2ndf_cut;
float xi_ht_ldl_cut;
float xi_ht_l_cut;
float xi_l_cut;
float xi_dl_cut;
float chi2primary_xi_proton_cut;
float chi2primary_xi_pi_cut;
float chi2primary_xi_bach_cut;
float chi2primary_xi_ht_cut;

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

int nhits_pi;
int nhits_he;

//float ht_chi2topo_cut;
//float ht_chi2ndf_cut;
//float ht_ldl_cut;
//float ht_l_cut;
//float chi2primary_he_cut;

float hl_chi2topo_cut;
float hl_chi2ndf_cut;
float hl_ldl_cut;
float hl_l_cut;
float chi2primary_h4_cut;

double ht_mass = 2.99131;
double ht_width = 0.005;
double hl_mass = 3.924;
double hl_width = 0.005;

//float ht_mass = 1.115683;

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
float  ht_pv_dca_ht_cut = 0.8;
float  p_pv_dca_ht_cut = 0.5;
float  pi_pv_dca_ht_cut = 1.5;
float  pi_p_dca_ht_cut = 0.8;
float  decaylen_ht_cut = 3.0;
float  ang0_ht_cut = 0.0;
*/

float ht_chi2topo_cut = 5     ;
float ht_chi2ndf_cut = 10     ;
float ht_ldl_cut = 5          ;
float ht_l_cut = 5            ;
//float ht_dl           ;
float chi2primary_he_cut = 18.6 ;
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

float ht_bdfvtx;
float ht_bdfvtx2;
float ht_lifetime;

float dca_he;
float dca_pi;

TFile  *fgpt_0;

TF1 *levyfit4;
TF1 *levyfit5;
TF1 *levyfit6;
TF1 *levyfit7;
TF1 *levyfit8;
TF1 *t_quadr;
TF1 *t_quadr0;
TF1 *t_quadr1;
TF1 *bolt0;
TF1 *bolt1;
TF1 *bolt2;

const int nrapbins = 9;
TH1D *g_pt_py[nrapbins];

int bismc;
//TMVA::Reader* reader[ncenttmvabins][npttmvabins];
// float bdt_ht_cut[ncenttmvabins][npttmvabins];

//obsolete, used for bdt scanning
float bdt_ht_cut[ncenttmvabins][nptbins];

//used for om analysis
float bdt_om_cut[ncentbins][nptbins];
float bdt_ob_cut[ncentbins][nptbins];

float weight;
int countrefmult;

int centrality;
int centFull[ncentbins] ={15,22,32,43,57,73,92,117,133};
//float r_cent[ncentbins+1] = {0., 0., 0.178389, 0.466883, 1.06469, 2.30611, 2.70053, 1.49454, 0.866928,0};//this is h3l
float r_cent[ncentbins+1] = {0.0693964, 0.166454, 0.178389*2.38874, 0.466883*1.09257, 1.06469*1.43903, 2.30611*0.893755, 2.70053*0.857322, 1.49454*0.814087, 0.866928*0.797996,0};
const int ncentbdtbins = 3;
const int nptbdtbins = 7;
float bdt_xi_cut[ncentbdtbins][nptbdtbins];
double ptbdtbin[nptbdtbins+1] = {0.2,0.6,1.0,1.4,1.8,2.2,2.6,3.0};
int choose_cent;
int _batch;

int _wgt;
int Centrality(int aa );
void book_histo();
void write_histo();

void process_tree_tmva_hl_kf_vmvtest(int condor=0, int cut_mode=0, int wgt=2, int batch=0){

    _condor = condor;//condor=1: condor settings
    _cut_mode = cut_mode;
    _wgt = wgt;
    _batch = batch;

    cout<<"_wgt:"<<_wgt<<endl;
    choose_cent	= 0;
//    choose_cent = 1;
//    choose_cent = 2;


levyfit4 = new TF1("levyfit4","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.4,2.4);
levyfit5 = new TF1("levyfit5","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.2,2.4);
levyfit6 = new TF1("levyfit6","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.2,2.4);
levyfit7 = new TF1("levyfit7","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.2,2.4);
levyfit8 = new TF1("levyfit8","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.2,2.4);
t_quadr = new TF1("t_quadr","[0]+[1]*x+[2]*x*x",-2,2);

t_quadr0 = new TF1("t_quadr0","[0]+[1]*x+[2]*x*x",-2,2);
t_quadr1 = new TF1("t_quadr1","[0]+[1]*x+[2]*x*x",-2,2);
//bolt0 = new TF1("bolt0", "x*sqrt(2.80941350568*2.80941350568+x*x)*exp(-sqrt(2.80941350568*2.80941350568+x*x)/[0])", 0.,5);
//bolt1 = new TF1("bolt1", "x*sqrt(3.72840015733*3.72840015733+x*x)*exp(-sqrt(3.72840015733*3.72840015733+x*x)/[0])", 0.,5);
bolt0 = new TF1("bolt0", "1e9*x*sqrt(3.924*3.924+x*x)*exp(-sqrt(3.924*3.924+x*x)/[0])", 0.,5);
bolt1 = new TF1("bolt1", "1e9*x*sqrt(3.924*3.924+x*x)*exp(-sqrt(3.924*3.924+x*x)/[0])", 0.,5);
bolt2 = new TF1("bolt2", "1e9*x*sqrt(3.924*3.924+x*x)*exp(-sqrt(3.924*3.924+x*x)/[0])", 0.,5);

levyfit4->SetParameters(1.53536e-01 , -1.21064e+08, 5.25287e+07);
levyfit5->SetParameters(1.51332e-01 , -1.25927e+08, 4.86193e+07);
levyfit6->SetParameters(1.45903e-01 , -1.60954e+08, 4.74115e+07);
levyfit7->SetParameters(1.28382e-01 , -1.60610e+08, 4.92025e+07);
levyfit8->SetParameters(1.11371e-01 , -2.04689e+08, 4.22255e+07);

//t_quadr->SetParameters(2.06539, 0.0, -1.75299);
//t_quadr0->SetParameters(3.04528e+01, 0.0, 6.71347e+01);
//t_quadr1->SetParameters(6.88265e+01, 0.0, 1.90472e+02);
        t_quadr->SetParameters(4.71255e-01,0.00000e+00,1.79971e+00*1.5);//this is from fit to 3gevv data, multiply second paramter by 1.5
        //t_quadr0->SetParameters(4.71255e-01,0.00000e+00,1.79971e+00);//this is from fit to 3gevv data,h3l
        //
        t_quadr0->SetParameters(7.03832e-01,0.,9.60834e-01);//this is from fit to 3gevv data,h4l
        //t_quadr0->SetParameters(7.03832e-01,0.,9.60834e-01*2);//this is from fit to 3gevv data,h4l

        t_quadr1->SetParameters(6.88265e+01, 0.0, 1.90472e+02);

//bolt0->SetParameter(0,0.407);
//bolt1->SetParameter(0,0.398);
	bolt0->SetParameter(0,0.13);
        bolt1->SetParameter(0,0.105);
        bolt2->SetParameter(0,0.08);

fgpt_0 = new TFile("h_gpt_ht.root","READ");
for(int irap=0;irap<nrapbins;irap++){
g_pt_py[irap] = (TH1D*)fgpt_0->Get(Form("g_pt_py[%d]",irap));
}

fgpt_0 = new TFile("hl_input_mc_fine.root","READ");
g_pt_fine_in = (TH2F*)fgpt_0->Get("g_pt_fine")->Clone("g_pt_fine_in");


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


TChain cascade_tree("cascade_tree");
TChain omega_tree("omega_tree");
TChain htriton_tree("htriton_tree");
TChain htriton_mc_tree("htriton_mc_tree");

int nChains = -999999;
ifstream fin;

if(_condor==0){//emb
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_mu/sum/filelist.txt");
//strcpy (ofile,"run18_3gev_ht_hist_mu.root");
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/filelist_ht_emb3.txt");
//strcpy (ofile,"run18_3gev_ht_hist_mu3.root");
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_mu3_correctedvertex/sum/filelist.txt");
//strcpy (ofile,"run18_3gev_ht_hist_mu4.root");

//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_mu3_correctedvertex_anderror/filelist.txt");
//strcpy (ofile,"run18_3gev_ht_hist_mu5.root");
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_mu3_correctedvertex_andnoerror/filelist.txt");
//strcpy (ofile,"run18_3gev_ht_hist_mu5noerr.root");

//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_mu3_correctedvertex_werror2/filelist.txt");
//strcpy (ofile, "run18_3gev_ht_hist_mu6.root");
//nChains = 1;
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_mu3_correctedvertex_werror4/filelist.txt");
//strcpy (ofile, "run18_3gev_ht_hist_mu7.root");

//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_mu3_correctedvertex_werror_veryloose/filelist.txt");
//strcpy (ofile, "run18_3gev_ht_hist_mu_veryloose.root");
//nChains = 1;

//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_mu3_correctedvertex_werror_3/sum/filelist.txt");
//strcpy (ofile, "run18_3gev_ht_hist_mu_3.root");
//nChains = 1;

//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_mu3_correctedvertex_w10xerror_3/sum/filelist.txt");
//strcpy (ofile, "run18_3gev_ht_hist_mu_10x_3.root");
//nChains = 1;

//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_mu3_correctedvertex_w5xerror_4/filelist.txt");
//strcpy (ofile, Form("run18_3gev_ht_hist_mu_5x_3_newreweight_centa%d.root",choose_cent));
//nChains = 1;

//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_mu3_correctedvertex_w5xerror_5/filelist.txt");
//strcpy (ofile, Form("run18_3gev_ht_hist_mu_5x_5_newreweight_centa%d.root",choose_cent));
//nChains = 326;

//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_mu/filelist.txt");
//strcpy (ofile, Form("run18_3gev_ht_hist_mu_1x_1_newreweight_centa%d.root",choose_cent));
//nChains = 90;

//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_hl_mu_4/filelist.txt");
//strcpy (ofile, Form("run18_3gev_hl_hist_mu_1x_4_newreweight_centa%d_wgt%d_cut%d.root",choose_cent,_wgt,_cut_mode));
//nChains = 262;

fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_hl_mu_newpid2_manualvtxerr/filelist.txt");
strcpy (ofile, Form("run18_3gev_hl_hist_mu_1x_8_newreweight_centa%d_wgt%d_cut%d_mve_vtest.root",choose_cent,_wgt,_cut_mode));
nChains = 1268;

}
if(_condor==1){//data
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_pico/sum/filelist.txt");
//strcpy (ofile,"run18_3gev_ht_hist_pico.root");
//
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_pico_rotate/sum/filelist.txt");
//strcpy (ofile,"run18_3gev_ht_hist_pico_rotate.root");
//nChains = 1;
//
//loose cuts
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_pico_tree/filelist.txt");
//strcpy (ofile,"run18_3gev_ht_hist_pico_loose.root");
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_pico_rotate_tree/filelist.txt");
//strcpy (ofile,"run18_3gev_ht_hist_pico_rotate_loose.root");
//nChains = 1200;

//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_pico_tree_2/filelist.txt");
//strcpy (ofile,"run18_3gev_ht_hist_pico_loose_7.root");
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_pico_rotate_tree_2/filelist.txt");
//strcpy (ofile,"run18_3gev_ht_hist_pico_rotate_loose_7.root");

//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_pico_tree_4/filelist.txt");
//strcpy (ofile,Form("run18_3gev_ht_hist_pico_loose_4_centa%d.root",choose_cent));
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_pico_rotate_tree_4/filelist.txt");
//strcpy (ofile,Form("run18_3gev_ht_hist_pico_rotate_loose_4_centa%d.root",choose_cent));
//nChains = 1200;

/*
        if(_batch==0){
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_pico_tree_newpid_0//filelist.txt");
//strcpy (ofile,Form("run18_3gev_hl_hist_pico_loose_7_centa%d_cut%d.root",choose_cent,_cut_mode));

fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_pico_tree_newpid2/filelist.txt");
strcpy (ofile,Form("run18_3gev_hl_hist_pico_loose_8_centa%d_cut%d_v9.root",choose_cent,_cut_mode));
nChains = 1222;
}
        if(_batch==1){
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_pico_tree_newpid_rotate_0/filelist.txt");
//strcpy (ofile,Form("run18_3gev_hl_hist_pico_rotate_0_loose_7_centa%d_cut%d.root",choose_cent,_cut_mode));

fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_pico_tree_newpid2_rotate_0/filelist.txt");
strcpy (ofile,Form("run18_3gev_hl_hist_pico_rotate_0_loose_8_centa%d_cut%d_v9.root",choose_cent,_cut_mode));
nChains = 1208;
}
        if(_batch==2){
fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_pico_tree_newpid2_rotate_1/filelist.txt");
strcpy (ofile,Form("run18_3gev_hl_hist_pico_rotate_1_loose_8_centa%d_cut%d_v9.root",choose_cent,_cut_mode));
nChains = 1219;
}
        if(_batch==3){
fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_pico_tree_newpid2_rotate_2/filelist.txt");
strcpy (ofile,Form("run18_3gev_hl_hist_pico_rotate_2_loose_8_centa%d_cut%d_v9.root",choose_cent,_cut_mode));
nChains = 1214;
}
        if(_batch==4){
fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_pico_tree_newpid2_rotate_3/filelist.txt");
strcpy (ofile,Form("run18_3gev_hl_hist_pico_rotate_3_loose_8_centa%d_cut%d_v9.root",choose_cent,_cut_mode));
nChains = 1220;
}
        if(_batch==5){
fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_pico_tree_newpid2_rotate_4/filelist.txt");
strcpy (ofile,Form("run18_3gev_hl_hist_pico_rotate_4_loose_8_centa%d_cut%d_v9.root",choose_cent,_cut_mode));
nChains = 1221;
}
if(_batch==6){
fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_pico_tree_newpid2_rotate_5/filelist.txt");
strcpy (ofile,Form("run18_3gev_hl_hist_pico_rotate_5_loose_8_centa%d_cut%d_v9.root",choose_cent,_cut_mode));
nChains = 1221;
}
if(_batch==7){
fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_pico_tree_newpid2_rotate_6/filelist.txt");
strcpy (ofile,Form("run18_3gev_hl_hist_pico_rotate_6_loose_8_centa%d_cut%d_v9.root",choose_cent,_cut_mode));
nChains = 1221;
}
*/
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

//      h4htriton_tree.Add(path_ch);
//      cout<<path_ch<<" "<<h4htriton_tree.GetEntries() <<endl;

htriton_tree.Add(path_ch);
cout<<path_ch<<" "<<htriton_tree.GetEntries() <<endl;
htriton_mc_tree.Add(path_ch);
cout<<path_ch<<" "<<htriton_mc_tree.GetEntries() <<endl;


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

//htriton_tree.SetBranchAddress("brunid",&brunid);
htriton_tree.SetBranchAddress("brefmult",&brefmult);
htriton_tree.SetBranchAddress("btofmult",&btofmult);
htriton_tree.SetBranchAddress("bparticleid",&bparticleid);
htriton_tree.SetBranchAddress("bparticlemass",&bparticlemass);
htriton_tree.SetBranchAddress("bpx",&bpx);
htriton_tree.SetBranchAddress("bpy",&bpy);
htriton_tree.SetBranchAddress("bpz",&bpz);
htriton_tree.SetBranchAddress("chi2primary_he", &chi2primary_he);
htriton_tree.SetBranchAddress("chi2primary_pi", &chi2primary_pi);
htriton_tree.SetBranchAddress("ht_chi2topo", &ht_chi2topo);
htriton_tree.SetBranchAddress("ht_chi2ndf", &ht_chi2ndf);
htriton_tree.SetBranchAddress("ht_ldl", &ht_ldl);
htriton_tree.SetBranchAddress("ht_l", &ht_l);
htriton_tree.SetBranchAddress("ht_dl", &ht_dl);
htriton_tree.SetBranchAddress("dca_he",&dca_he);
htriton_tree.SetBranchAddress("dca_pi",&dca_pi);
htriton_tree.SetBranchAddress("nhits_pi",&nhits_pi);
htriton_tree.SetBranchAddress("nhits_he",&nhits_he);
htriton_tree.SetBranchAddress("ht_bdfvtx",&ht_bdfvtx);
htriton_tree.SetBranchAddress("ht_bdfvtx2",&ht_bdfvtx2);
htriton_tree.SetBranchAddress("ht_lifetime",&ht_lifetime);
htriton_tree.SetBranchAddress("countrefmult",&countrefmult);
htriton_tree.SetBranchAddress("reweight", &reweight);
htriton_tree.SetBranchAddress("cent9", &cent9);
htriton_tree.SetBranchAddress("bismc", &bismc);

if(_condor==0){
htriton_tree.SetBranchAddress("bmcpx",&bmcpx);
htriton_tree.SetBranchAddress("bmcpy",&bmcpy);
htriton_tree.SetBranchAddress("bmcpz",&bmcpz);
htriton_tree.SetBranchAddress("bmcl",&bmcl);
htriton_tree.SetBranchAddress("bmcpl",&bmcpl);
htriton_tree.SetBranchAddress("bpl",&bpl);
}
//htriton_mc_tree.SetBranchAddress("brunid", &brunid);
htriton_mc_tree.SetBranchAddress("brefmult", &brefmult);
htriton_mc_tree.SetBranchAddress("bmcparticleid", &bmcparticleid);
htriton_mc_tree.SetBranchAddress("bmcrawpx", &bmcrawpx);
htriton_mc_tree.SetBranchAddress("bmcrawpy", &bmcrawpy);
htriton_mc_tree.SetBranchAddress("bmcrawpz", &bmcrawpz);
htriton_mc_tree.SetBranchAddress("bmcrawpl", &bmcrawpl);
htriton_mc_tree.SetBranchAddress("reweight", &reweight);
htriton_mc_tree.SetBranchAddress("countrefmult", &mccountrefmult);
htriton_mc_tree.SetBranchAddress("cent9", &cent9);

Long64_t n_lambda_Entries = htriton_tree.GetEntries();
Long64_t n_lambda_mc_Entries = htriton_mc_tree.GetEntries();

 cout<<endl<<"Start processing "<<endl<<endl;
cout<<"--------------------------------------------"<<endl;

cout<<"Total entries "<<n_lambda_Entries<< endl;

if( _cut_mode==1 || (_cut_mode>=11 && _cut_mode<20) ){//loose
        ht_chi2topo_cut = 5;
        ht_chi2ndf_cut = 5;
chi2primary_he_cut = 3.0;
        chi2primary_pi_cut = 3.0;

        hl_chi2topo_cut = 5;
        hl_chi2ndf_cut = 5;
chi2primary_h4_cut = 3.0;
        chi2primary_pi_cut = 3.0;
        }
if(_cut_mode==2 || _cut_mode==61 || _cut_mode==62){
        ht_chi2topo_cut = 4.;
        hl_chi2topo_cut = 4.;
        ht_chi2ndf_cut = 4;
        hl_chi2ndf_cut = 4;
        chi2primary_he_cut = 3.;
        chi2primary_pi_cut = 20.;
        chi2primary_h4_cut = 3.;
        dca_he_cut = 0.;
        }
        if(_cut_mode==3){
        ht_chi2topo_cut = 5.;
        hl_chi2topo_cut = 5.;
        ht_chi2ndf_cut = 5;
        hl_chi2ndf_cut = 5;
        chi2primary_he_cut = 3.;
        chi2primary_pi_cut = 3.;
        chi2primary_h4_cut = 3.;
        dca_he_cut = 1.0;
        }
        if(_cut_mode==4){
        ht_chi2topo_cut = 5.;
        hl_chi2topo_cut = 5.;
        ht_chi2ndf_cut = 5;
        hl_chi2ndf_cut = 5;
        chi2primary_he_cut = 10;
        chi2primary_pi_cut = 3.;
        chi2primary_h4_cut = 10.;
        }
        if(_cut_mode==5){
        ht_chi2topo_cut = 5.;
        hl_chi2topo_cut = 5.;
        ht_chi2ndf_cut = 5;
        hl_chi2ndf_cut = 5;
        chi2primary_he_cut = 20.;
        chi2primary_pi_cut = 3.;
        chi2primary_h4_cut = 20.;
        }
if(_cut_mode==6){
        ht_chi2topo_cut = 5.;
        hl_chi2topo_cut = 5.;
        ht_chi2ndf_cut = 5;
        hl_chi2ndf_cut = 5;
        chi2primary_he_cut = 30.;
        chi2primary_pi_cut = 3.;
        chi2primary_h4_cut = 30.;
        }
        if(_cut_mode==7){
        ht_chi2topo_cut = 5.;
        hl_chi2topo_cut = 5.;
        ht_chi2ndf_cut = 5;
        hl_chi2ndf_cut = 5;
        chi2primary_he_cut = 40.;
        chi2primary_pi_cut = 3.;
        chi2primary_h4_cut = 40.;
        }
        if(_cut_mode==8){
        ht_chi2topo_cut = 5.;
        hl_chi2topo_cut = 5.;
        ht_chi2ndf_cut = 5;
        hl_chi2ndf_cut = 5;
        chi2primary_he_cut = 3.;
        chi2primary_pi_cut = 3.;
        chi2primary_h4_cut = 3.;
        dca_he_cut = 1.5;
        }

        if(_cut_mode==21){
	ht_chi2topo_cut = 4.;
        hl_chi2topo_cut = 4.;
        ht_chi2ndf_cut = 4;
        hl_chi2ndf_cut = 4;
        chi2primary_he_cut = 3.;
        chi2primary_pi_cut = 3.;
        chi2primary_h4_cut = 3.;
        }
	if(_cut_mode==22){
	ht_chi2topo_cut = 4.;
        hl_chi2topo_cut = 4.;
        ht_chi2ndf_cut = 4;
        hl_chi2ndf_cut = 4;
        chi2primary_he_cut = 3.;
        chi2primary_pi_cut = 40.;
        chi2primary_h4_cut = 3.;
        }
        if(_cut_mode==31){
	ht_chi2topo_cut = 4.;
        hl_chi2topo_cut = 4.;
        ht_chi2ndf_cut = 4;
        hl_chi2ndf_cut = 4;
        chi2primary_he_cut = 10.;
        chi2primary_pi_cut = 20.;
        chi2primary_h4_cut = 3.;
        }
        if(_cut_mode==32){
	ht_chi2topo_cut = 4.;
        hl_chi2topo_cut = 4.;
        ht_chi2ndf_cut = 4;
        hl_chi2ndf_cut = 4;
        chi2primary_he_cut = 20.;
        chi2primary_pi_cut = 20.;
        chi2primary_h4_cut = 3.;
        }
        if(_cut_mode==41){
	ht_chi2topo_cut = 3.;
        hl_chi2topo_cut = 3.;
        ht_chi2ndf_cut = 4;
        hl_chi2ndf_cut = 4;
        chi2primary_he_cut = 3.;
        chi2primary_pi_cut = 20.;
        chi2primary_h4_cut = 3.;
        }
	if(_cut_mode==42){
	ht_chi2topo_cut = 5.;
        hl_chi2topo_cut = 5.;
        ht_chi2ndf_cut = 4;
        hl_chi2ndf_cut = 4;
        chi2primary_he_cut = 3.;
        chi2primary_pi_cut = 20.;
        chi2primary_h4_cut = 3.;
	}if(_cut_mode==51){
	ht_chi2topo_cut = 4.;
        hl_chi2topo_cut = 4.;
        ht_chi2ndf_cut = 3;
        hl_chi2ndf_cut = 3;
        chi2primary_he_cut = 3.;
        chi2primary_pi_cut = 20.;
        chi2primary_h4_cut = 3.;
        }
        if(_cut_mode==52){
	ht_chi2topo_cut = 4.;
        hl_chi2topo_cut = 4.;
        ht_chi2ndf_cut = 5;
        hl_chi2ndf_cut = 5;
        chi2primary_he_cut = 3.;
        chi2primary_pi_cut = 20.;
        chi2primary_h4_cut = 3.;
	}



//n_lambda_Entries = 100000;
for (Long64_t iEntry = 0; iEntry<=n_lambda_Entries; iEntry++){

	htriton_tree.GetEntry(iEntry);
  if(iEntry%100000==0){
    cout<<"Processing entry:" <<iEntry<<endl;
}
int centtmva_label = -1;
int centbdt_label = -1;
 cent_label = -1;
int pt_label = -1;




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
//ooverwrite
//ycm = -2.03;

 //mass_lo = 1.105;
 //mass_lo = 2.988;
 //mass_hi = 1.125;
 //mass_hi = 2.998;

 mass_lo = 3.919;
 mass_hi = 3.929;

   bpt = sqrt(bpx*bpx+bpy*bpy);
//  bp = 0.5*sqrt(bpx*bpx+bpy*bpy+bpz*bpz);
   bp = sqrt(bpx*bpx+bpy*bpy+bpz*bpz);
   _rotate = 0;


  TLorentzVector ptc(bpx,bpy,bpz,sqrt(bpx*bpx+bpy*bpy+bpz*bpz+bparticlemass*bparticlemass));
  beta = ptc.Eta(); 
  brap = ptc.Rapidity() - ycm;

// TLorentzVector mcptc(bmcpx,bmcpy,bmcpz,sqrt(bmcpx*bmcpx+bmcpy*bmcpy+bmcpz*bmcpz+ht_mass*ht_mass));	
// bmcrap = mcptc.Rapidity() - ycm; 

    double betaxgamma;
    betaxgamma = sqrt(bpx*bpx+bpy*bpy+bpz*bpz)/bparticlemass;
    double ht_pl;
    ht_pl = ht_l/betaxgamma;


//place cut here
if(chi2primary_he<chi2primary_he_cut || chi2primary_pi<chi2primary_pi_cut) continue;
if(ht_chi2topo>ht_chi2topo_cut) continue;
if(ht_chi2ndf>ht_chi2ndf_cut) continue;
//        if( (_cut_mode==2 || _cut_mode==3 || _cut_mode==8) && dca_he<dca_he_cut) continue;

        if(_cut_mode==61 && nhits_pi<17) continue;
        if(_cut_mode==61 && nhits_he<17) continue;

        if(_cut_mode==62 && nhits_pi<20) continue;
        if(_cut_mode==62 && nhits_he<20) continue;

	if(ht_pl>50) continue;//no signal for ht_pl>50
        if(ht_pl<2) continue;//no signal for ht_pl<2

if(bparticlemass>mass_lo && bparticlemass<mass_hi){
h_pt->Fill(brap,bpt);
h_r_pt->Fill(-brap,bpt);
}

//fiducial cut
//if( bpt>3.0) continue;
//if( brap>0.95) continue;
//if( brap>0.7  &&  bpt < 4.9*(brap-0.50)*(brap-0.50)+0.40+0.017 ) continue;
//if( brap<0.7  &&  bpt < 2.25*(brap-0.80)*(brap-0.80)+0.5-0.029+0.1195 ) continue;

if(bparticleid==3005 && ((bismc==1&&_condor==0) || _condor==1 || _condor==2) ){

h_mass_pt->Fill(bparticlemass,brap,bpt);

if( bpt>3.0) continue;
if( brap>0.95) continue;
if( brap>0.7  &&  bpt < 4.9*(brap-0.50)*(brap-0.50)+0.40+0.017 ) continue;
if( brap<0.7  &&  bpt < 2.25*(brap-0.80)*(brap-0.80)+0.5-0.029+0.1195 ) continue;

if(bparticlemass>mass_lo && bparticlemass<mass_hi && !_condor==0){




h_chi2primary_he->Fill(beta,bpt,chi2primary_he);
h_chi2primary_pi->Fill(beta,bpt,chi2primary_pi);
h_nhits_he->Fill(beta,bpt,nhits_he);
h_nhits_pi->Fill(beta,bpt,nhits_pi);
h_ht_ldl->Fill(beta,bpt,ht_ldl);
h_ht_l->Fill(beta,bpt,ht_l);
h_ht_dl->Fill(beta,bpt,ht_dl);
h_ht_chi2topo->Fill(beta,bpt,ht_chi2topo);
h_ht_chi2ndf->Fill(beta,bpt,ht_chi2ndf);
h_ht_dca_he->Fill(beta,bpt,dca_he);
h_ht_dca_pi->Fill(beta,bpt,dca_pi);
h_eta->Fill(beta,bpt);
 h_rap->Fill(-brap,bpt);
h_ht_bdfvtx->Fill(beta,bpt,ht_bdfvtx);
h_ht_bdfvtx2->Fill(beta,bpt,ht_bdfvtx2);
h_ht_lifetime->Fill(beta,bpt,ht_lifetime);

        h_cent->Fill(cent_label);
}

h_mass_qa->Fill(beta,bpt,bparticlemass);
//h_mass_pt->Fill(bparticlemass,brap,bpt);
h_mass_pl->Fill(bparticlemass,brap,ht_pl);

if(bismc==1&&_condor==0){

 TLorentzVector mcptc(bmcpx,bmcpy,bmcpz,sqrt(bmcpx*bmcpx+bmcpy*bmcpy+bmcpz*bmcpz+hl_mass*hl_mass));
 bmcrap = mcptc.Rapidity() - ycm;

     float mc_weight;
     mc_weight = 1./g_pt_fine_in->GetBinContent(g_pt_fine_in->FindBin(bmcrap,sqrt(bmcpx*bmcpx+bmcpy*bmcpy))); 
/*
if(bmcrap<-0.7){ weight = levyfit8->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy))/g_pt_py[0]->GetBinContent(g_pt_py[0]->FindBin(sqrt(bmcpx*bmcpx+bmcpy*bmcpy))); }
else if(bmcrap<-0.5){ weight = levyfit7->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy))/g_pt_py[1]->GetBinContent(g_pt_py[1]->FindBin(sqrt(bmcpx*bmcpx+bmcpy*bmcpy))); }
else if(bmcrap<-0.3){ weight = levyfit6->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy))/g_pt_py[2]->GetBinContent(g_pt_py[2]->FindBin(sqrt(bmcpx*bmcpx+bmcpy*bmcpy))); }
else if(bmcrap<-0.1){ weight = levyfit5->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy))/g_pt_py[3]->GetBinContent(g_pt_py[3]->FindBin(sqrt(bmcpx*bmcpx+bmcpy*bmcpy))); }
else if(bmcrap<0.1) { weight = levyfit4->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy))/g_pt_py[4]->GetBinContent(g_pt_py[4]->FindBin(sqrt(bmcpx*bmcpx+bmcpy*bmcpy))); }
else if(bmcrap<0.3) { weight = levyfit5->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy))/g_pt_py[5]->GetBinContent(g_pt_py[5]->FindBin(sqrt(bmcpx*bmcpx+bmcpy*bmcpy))); }
else if(bmcrap<0.5) { weight = levyfit6->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy))/g_pt_py[6]->GetBinContent(g_pt_py[6]->FindBin(sqrt(bmcpx*bmcpx+bmcpy*bmcpy))); }
else if(bmcrap<0.7) { weight = levyfit7->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy))/g_pt_py[7]->GetBinContent(g_pt_py[7]->FindBin(sqrt(bmcpx*bmcpx+bmcpy*bmcpy))); }
else { weight = levyfit8->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy))/g_pt_py[8]->GetBinContent(g_pt_py[8]->FindBin(sqrt(bmcpx*bmcpx+bmcpy*bmcpy))); }
*/

     float rap_weight;
     float pt_weight;
        float cent_weight;

        cent_weight = r_cent[cent_label];

if(_wgt==0){
        rap_weight=1.;
        pt_weight = 1.;
        }else if(_wgt==1){
        rap_weight = t_quadr->Eval(bmcrap);
        pt_weight = bolt1->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy));
        }else if(_wgt==2){
        rap_weight = t_quadr0->Eval(bmcrap);
        pt_weight =bolt0->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy));
        }else if(_wgt==3){
rap_weight = 1;
        pt_weight =bolt0->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy));
        }else if(_wgt==4){
        rap_weight = t_quadr0->Eval(bmcrap);
        pt_weight =bolt1->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy));
        }else if(_wgt==5){
        rap_weight = t_quadr0->Eval(bmcrap+2.03-1.045);
        pt_weight =bolt2->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy));
        }     

weight = rap_weight*pt_weight;
         weight *= mc_weight;
        weight *= cent_weight;

h_mass_pt_wgt->Fill(bparticlemass,brap,bpt,weight);

if( bpt>3.0) continue;
if( brap>0.95) continue;
if( brap>0.7  &&  bpt < 4.9*(brap-0.50)*(brap-0.50)+0.40+0.017 ) continue;
if( brap<0.7  &&  bpt < 2.25*(brap-0.80)*(brap-0.80)+0.5-0.029+0.1195 ) continue;

h_mass_pl_wgt->Fill(bparticlemass,brap,ht_pl,weight);

if(bparticleid==3005 && (bismc==1&&_condor==0)){
	if(bparticlemass>mass_lo && bparticlemass<mass_hi){

h_chi2primary_he->Fill(beta,bpt,chi2primary_he,weight);
h_chi2primary_pi->Fill(beta,bpt,chi2primary_pi,weight);
h_nhits_he->Fill(beta,bpt,nhits_he,weight);
h_nhits_pi->Fill(beta,bpt,nhits_pi,weight);
h_ht_ldl->Fill(beta,bpt,ht_ldl,weight);
h_ht_l->Fill(beta,bpt,ht_l,weight);
h_ht_dl->Fill(beta,bpt,ht_dl,weight);
h_ht_chi2topo->Fill(beta,bpt,ht_chi2topo,weight);
h_ht_chi2ndf->Fill(beta,bpt,ht_chi2ndf,weight);
h_ht_dca_he->Fill(beta,bpt,dca_he,weight);
h_ht_dca_pi->Fill(beta,bpt,dca_pi,weight);
        h_rap->Fill(-(brap+ycm+2.03),bpt,weight);
h_eta->Fill(beta,bpt,weight);
h_ht_bdfvtx->Fill(beta,bpt,ht_bdfvtx,weight);
h_ht_bdfvtx2->Fill(beta,bpt,ht_bdfvtx2,weight);
h_ht_lifetime->Fill(beta,bpt,ht_lifetime,weight);
        h_cent->Fill(cent_label,cent_weight);

	}
}

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
h_l_mc->Fill(ht_l,bmcl);
h_pl_mc->Fill(ht_pl,bmcpl);
h_lt_mc->Fill(ht_lifetime,bmcpl);

g_p_mc->Fill(bp,bp-sqrt(bmcpx*bmcpx+bmcpy*bmcpy+bmcpz*bmcpz));
g_pt_mc->Fill(bpt,bpt-sqrt(bmcpx*bmcpx+bmcpy*bmcpy));
g_l_mc->Fill(ht_l,ht_l-bmcl);
g_pl_mc->Fill(ht_pl,ht_pl-bmcpl);
g_lt_mc->Fill(ht_lifetime,ht_lifetime-bmcpl);
}


}//bparticleid loop

}//entry loop

for (Long64_t iEntry = 0; iEntry<=n_lambda_mc_Entries; iEntry++){

                 htriton_mc_tree.GetEntry(iEntry);
           if(iEntry%100000==0){
             cout<<"Processing mc entry:" <<iEntry<<endl;
         }
if(bmcparticleid!=3005) continue;

cent_label=-1;
if(snn==3){
      cent_label = Centrality(mccountrefmult);
  }
if(cent_label<0) continue;//<=5refmult
     if(cent_label<7 && choose_cent==1) continue;//cent7,8,9
     if( (cent_label>6 || cent_label<2) && choose_cent==2) continue;//cent2,3,4,5,6

         bpt = sqrt(bmcrawpx*bmcrawpx+bmcrawpy*bmcrawpy);

         TLorentzVector ptc(bmcrawpx,bmcrawpy,bmcrawpz,sqrt(bmcrawpx*bmcrawpx+bmcrawpy*bmcrawpy+bmcrawpz*bmcrawpz+hl_mass*hl_mass));
         beta = ptc.Eta();
         brap = ptc.Rapidity() - ycm;

	double betaxgamma;
	betaxgamma = sqrt(bmcrawpx*bmcrawpx+bmcrawpy*bmcrawpy+bmcrawpz*bmcrawpz)/hl_mass;
	double ht_pl;
	//ht_pl = bmcrawpl/betaxgamma;
	ht_pl = bmcrawpl;
        g_pt->Fill(brap,bpt);
	g_pt_fine->Fill(brap,bpt);
        g_pt_fine2->Fill(brap,bpt);

//fiducial cut
if(bpt>3.0) continue;
if(brap>0.95) continue;
if( brap>0.7  &&  bpt < 4.9*(brap-0.50)*(brap-0.50)+0.40+0.017 ) continue;
if( brap<0.7  &&  bpt < 2.25*(brap-0.80)*(brap-0.80)+0.5-0.029+0.1195 ) continue;


         g_eta->Fill(beta, bpt);
	 g_pl->Fill(brap,ht_pl);

/*
if(fabs(brap)<0.1){ weight = levyfit4->Eval(bpt); }
else if(fabs(brap)>=0.1 && fabs(brap)<0.3){ weight = levyfit5->Eval(bpt); }
else if(fabs(brap)>=0.3 && fabs(brap)<0.5){ weight = levyfit6->Eval(bpt); }
else if(fabs(brap)>=0.5 && fabs(brap)<0.7){ weight = levyfit7->Eval(bpt); }
else{ weight = levyfit8->Eval(bpt); }
//weight = levyfit7->Eval(bpt);
*/


/*
     if(brap<-0.7){ weight = levyfit8->Eval(bpt)/g_pt_py[0]->GetBinContent(g_pt_py[0]->FindBin(bpt)); }
     else if(brap<-0.5){ weight = levyfit7->Eval(bpt)/g_pt_py[1]->GetBinContent(g_pt_py[1]->FindBin(bpt)); }
     else if(brap<-0.3){ weight = levyfit6->Eval(bpt)/g_pt_py[2]->GetBinContent(g_pt_py[2]->FindBin(bpt)); }
     else if(brap<-0.1){ weight = levyfit5->Eval(bpt)/g_pt_py[3]->GetBinContent(g_pt_py[3]->FindBin(bpt)); }
     else if(brap<0.1) { weight = levyfit4->Eval(bpt)/g_pt_py[4]->GetBinContent(g_pt_py[4]->FindBin(bpt)); }
     else if(brap<0.3) { weight = levyfit5->Eval(bpt)/g_pt_py[5]->GetBinContent(g_pt_py[5]->FindBin(bpt)); }
     else if(brap<0.5) { weight = levyfit6->Eval(bpt)/g_pt_py[6]->GetBinContent(g_pt_py[6]->FindBin(bpt)); }
     else if(brap<0.7) { weight = levyfit7->Eval(bpt)/g_pt_py[7]->GetBinContent(g_pt_py[7]->FindBin(bpt)); }
     else { weight = levyfit8->Eval(bpt)/g_pt_py[8]->GetBinContent(g_pt_py[8]->FindBin(bpt)); }
*/
float mc_weight;
     mc_weight = 1./g_pt_fine_in->GetBinContent(g_pt_fine_in->FindBin(brap,bpt));
float rap_weight;
float pt_weight;
        float cent_weight;

        cent_weight = r_cent[cent_label];

if(_wgt==0){
        rap_weight=1.;
        pt_weight = 1.;
}else if(_wgt==1){
        rap_weight = t_quadr->Eval(brap);
        pt_weight = bolt1->Eval(bpt);
}else if(_wgt==2){
        rap_weight = t_quadr0->Eval(brap);
        pt_weight =bolt0->Eval(bpt);
}else if(_wgt==3){
rap_weight=1;
pt_weight =bolt0->Eval(bpt);
}else if(_wgt==4){
        rap_weight = t_quadr0->Eval(brap);
        pt_weight = bolt1->Eval(bpt);
}else if(_wgt==5){
        rap_weight = t_quadr0->Eval(brap);
        pt_weight = bolt2->Eval(bpt);
}


     weight = rap_weight*pt_weight;
     weight *= mc_weight;
        weight *= cent_weight;


g_pt_wgt->Fill(brap,bpt,weight);
g_pl_wgt->Fill(brap,ht_pl,weight);
g_pt_fine_wgt->Fill(brap,bpt,weight);
g_pt_fine2_wgt->Fill(brap,bpt,weight);

         }

write_histo();

}

void book_histo()
{

  h_dedx_p = new TH2F("h_dedx_p","",500,0,5,2000,0,100);
  h_dedx_p->Sumw2();
 
  double min_mass = 3.5;
  double max_mass = 4.5;

  h_mass_pt = new TH3F("h_mass_pt","",1000,min_mass,max_mass,10,-1.0,1.0,50,0,5);
  h_mass_pl = new TH3F("h_mass_pl","",1000,min_mass,max_mass,10,-1.0,1.0,500,0,100);

  h_mass_pt_wgt = new TH3F("h_mass_pt_wgt","",1000,min_mass,max_mass,10,-1.0,1.0,50,0,5);
  h_mass_pl_wgt = new TH3F("h_mass_pl_wgt","",1000,min_mass,max_mass,10,-1.0,1.0,500,0,100);

  h_mass_qa = new TH3F("h_mass_qa","",10,-2.0,0.0,50,0,5,1000,min_mass,max_mass);
  h_chi2primary_he = new TH3F("h_chi2primary_he","",10,-2.0,0.0,50,0,5,500,0,500);
  h_chi2primary_pi = new TH3F("h_chi2primary_pi","",10,-2.0,0.0,50,0,5,500,0,500);
  h_ht_ldl = new TH3F("h_ht_ldl","",10,-2.0,0.0,50,0,5,1000,0,200);
  h_ht_l = new TH3F("h_ht_l","",10,-2.0,0.0,50,0,5,1000,0,200);
  h_ht_dl = new TH3F("h_ht_dl","",10,-2.0,0.0,50,0,5,500,0,20);
  h_ht_chi2topo = new TH3F("h_ht_chi2topo","",10,-2.0,0.0,50,0,5,500,0,10);
  h_ht_chi2ndf = new TH3F("h_ht_chi2ndf","",10,-2.0,0.0,50,0,5,100,0,10);
  h_nhits_he = new TH3F("h_nhits_he","",10,-2.0,0.0,50,0,5,80,0,80);
  h_nhits_pi = new TH3F("h_nhits_pi","",10,-2.0,0.0,50,0,5,80,0,80);

  h_ht_dca_he = new TH3F("h_ht_dca_he","",10,-2.0,0.0,50,0,5,200,0,20);
  h_ht_dca_pi     = new TH3F("h_ht_dca_pi","",10,-2.0,0.0,50,0,5,200,0,20);

  h_ht_bdfvtx = new TH3F("h_ht_bdfvtx","",10,-2.0,0.0,50,0,5,1000,0,5);
  h_ht_bdfvtx2 = new TH3F("h_ht_bdfvtx2","",10,-2.0,0.0,50,0,5,1000,0,1);
  h_ht_lifetime = new TH3F("h_ht_lifetime","",10,-2.0,0.0,50,0,5,500,0,100);

h_ht_bdt = new TH3F("h_ht_bdt","",10,-2.0,0.0,50,0,5,200,-1.,1.);
  h_cent = new TH1F("h_cent","", 9,0,9);

  h_bnhits = new TH3F("h_bnhits","",10,-2.0,0.0,50,0,5,70,0,70);
  h_bdedx = new TH3F("h_bdedx","",10,-2.0,0.0,50,0,5,100,0,100);
  h_pt = new TH2F("h_pt","",100,-1.5,1.5,100,0,10);
  h_r_pt = new TH2F("h_r_pt","",100,-1.5,1.5,100,0,10);
  h_eta = new TH2F("h_eta","",200,-2.5,0.5,100,0,10);
  h_rap = new TH2F("h_rap","",200,-3.0,0.0,100,0,10);
  h_dca = new TH3F("h_dca","",10,-2.0,0.0,50,0,5,500,0,20);
  g_pt = new TH2F("g_pt","",10,-1.0,1.0,50,0,5);
  g_pl = new TH2F("g_pl","",10,-1.0,1.0,500,0,100);
  g_pt_fine = new TH2F("g_pt_fine","",100,-1.0,1.0,50,0,5);
  g_pt_fine2 = new TH2F("g_pt_fine2","",20,-1.0,1.0,50,0,5);
  g_pt_fine_wgt = new TH2F("g_pt_fine_wgt","",100,-1.0,1.0,50,0,5);
  g_pt_fine2_wgt = new TH2F("g_pt_fine2_wgt","",20,-1.0,1.0,50,0,5);

  g_pt_wgt = new TH2F("g_pt_wgt","",100,-1.0,1.0,50,0,5);
  g_pl_wgt = new TH2F("g_pl_wgt","",10,-1.0,1.0,500,0,100);

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
  h_mass_pt->Sumw2();
  h_mass_pl->Sumw2();

  h_ht_bdfvtx->Sumw2();
  h_ht_bdfvtx2->Sumw2();
  h_ht_lifetime->Sumw2();
  h_lt_mc->Sumw2();
  g_lt_mc->Sumw2();
h_eta->Sumw2();  
  h_ht_dca_he->Sumw2();
  h_ht_dca_pi->Sumw2();
  h_mass_qa->Sumw2();
  h_chi2primary_he ->Sumw2();
  h_chi2primary_pi ->Sumw2();
  h_bnhits ->Sumw2();
  h_bdedx ->Sumw2();
  h_dca ->Sumw2();
  h_pt->Sumw2();
  h_r_pt->Sumw2();
  g_pt->Sumw2();
  g_eta->Sumw2();
  h_mass_pl->Sumw2();
  h_ht_chi2ndf->Sumw2();
  h_ht_chi2topo->Sumw2();
  h_ht_ldl->Sumw2();
  h_ht_l->Sumw2();
h_rap->Sumw2();
  h_cent->Sumw2();

  h_nhits_he ->Sumw2();
  h_nhits_pi ->Sumw2();

  h_p_mc->Sumw2();
  h_pt_mc->Sumw2();
  h_l_mc->Sumw2();
  h_pl_mc->Sumw2();

  g_p_mc->Sumw2();
  g_pt_mc->Sumw2();
  g_l_mc->Sumw2();
  g_pl_mc->Sumw2();

  g_pt_fine->Sumw2();
  g_pt_fine_wgt->Sumw2();
  g_pt_fine2->Sumw2();
  g_pt_fine2_wgt->Sumw2();
}

void write_histo()
{
    TFile *outhistfile = new TFile (ofile, "RECREATE");
    outhistfile->cd();

h_mass_qa->Write();
h_mass_pt->Write();
h_mass_pl->Write();
h_chi2primary_he->Write();
h_chi2primary_pi->Write();
h_ht_ldl->Write();
h_ht_l->Write();
h_ht_dl->Write();
h_ht_chi2topo->Write();
h_ht_chi2ndf->Write();
h_ht_dca_he->Write();
h_ht_dca_pi->Write();
h_nhits_he->Write();
h_nhits_pi ->Write();
h_eta->Write();
  h_rap->Write();
  h_cent->Write();
h_dedx_p->Write();
h_bnhits->Write();
h_bdedx->Write();
h_pt->Write();
h_r_pt->Write();
g_pt->Write();
g_pl->Write();
h_dca ->Write();
g_eta->Write();
h_ht_bdfvtx->Write();
h_ht_bdfvtx2->Write();
h_ht_lifetime->Write();
g_pt_wgt->Write();
g_pl_wgt->Write();
h_mass_pt_wgt->Write();
h_mass_pl_wgt->Write();
h_ht_bdt->Write();//dummy
g_pt_fine_wgt->Write();
g_pt_fine->Write();
g_pt_fine2_wgt->Write();
g_pt_fine2->Write();

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

