//v15: new data, new centrality and new evt plane res
//v4: add different sys
//v3: add gweight, and eff
//v2: add resolution average method, new cent binning
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
#include "TF2.h"
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
const int ncentbins = 9;

double lambda_mass;
double cascade_mass;
double omega_mass;

double _weight;

TF1 *om_width_0;
TF1 *om_width_1;


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

float gweight;

//const int ncutptbins = 3;
//double ptcutbin[ncutptbins+1] = {0.0,1.0,2.0,3.0};
//const int ncutcentbins = 3;
//int centcutbin[ncutcentbins+1] = {0,10,40,80};
//const int ncutcentbins = 4;
//int centcutbin[ncutcentbins+1] = {0,10,20,40,60};
const int ncutcentbins = 3;
int centcutbin[ncutcentbins+1] = {0,10,40,60};

//const int ncutptbins = 4;
//double ptcutbin[ncutptbins+1] = {0.4,0.8,1.2,1.6,2.8};

//const int ncutptbins = 2;
//double ptcutbin[ncutptbins+1] = {0.4,1.0,4.0};

const int ncutptbins = 4;
double ptcutbin[ncutptbins+1] = {0.4,0.8,1.2,1.6,2.0};

//const int ncutcentbins = 9;
//int centcutbin[ncutcentbins+1] = {0,5,10,20,30,40,50,60,80};


//const int ncutptbins = 1;
//double ptcutbin[ncutptbins+1] = {0.0,5.0};
//const int ncutcentbins = 1;
//int centcutbin[ncutcentbins+1] = {0,100};

//const int nrapbins = 3;
//float rapbin[nrapbins+1] = {-0.5,0,0.5,1};
const int nrapbins = 15;
float rapbin[nrapbins+1] = {-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};

const int nrapbins_2 = 3;
float rapbin_2[nrapbins_2+1] = {-0.5,0.0,0.5,1.0};

TH2F *h_om_pt_rap[ncutcentbins];
TH2F *h_r_om_pt_rap[ncutcentbins];
TH2F *h_pt_eff_wgt_fine2;

TH1F *h_om_chi2topo[ncutcentbins][ncutptbins];
TH1F *h_om_chi2ndf[ncutcentbins][ncutptbins];
TH1F *h_om_chi2topo_bg[ncutcentbins][ncutptbins];
TH1F *h_om_chi2ndf_bg[ncutcentbins][ncutptbins];

TH1F *h_om_ldl[ncutcentbins][ncutptbins];
TH1F *h_om_ld_chi2topo[ncutcentbins][ncutptbins];

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

TH1F *h_om_mass[ncentbins][ncutptbins][nrapbins];
//TH2F *h_om_mass_phi_1o[ncutcentbins][ncutptbins][nrapbins];
//TH2F *h_om_mass_phi_2o[ncutcentbins][ncutptbins];

TH2F *h_om_mass_phi_2o[ncutcentbins][ncutptbins][nrapbins_2];
TH2F *h_om_mass_phi_1o[ncutcentbins][ncutptbins][nrapbins];

TH2F *h_om_mass_phi_2o_eff[ncutcentbins][ncutptbins][nrapbins_2];
TH2F *h_om_mass_phi_1o_eff[ncutcentbins][ncutptbins][nrapbins];

TH1F *h_om_phi_evp[ncutcentbins][ncutptbins];

int centtmva_limit[ncenttmvabins+1] = {0,10,40,80};
int cent_limit[ncentbins+1] = {0,5,10,20,30,40,50,60,70,80};

//from central to peripheral
//float _fullres_1_AB[9] = {0.287265,0.481901,0.649618,0.738704,0.740876,0.683774,0.575272,0.446576,0.346197};
//float _fullres_2_AB[9] = {0.0534926,0.156073,0.299369,0.403934,0.406805,0.336636,0.228438,0.132919,0.0783611};

//double res1_eABeCtB[9]  = {0.414801, 0.52164,  0.635928, 0.717564, 0.750364, 0.726469, 0.620534, 0.449211, 0.277282};//from shaowei
//double res12_eABeCtB[9] = {0.113904, 0.184801, 0.285307, 0.376866, 0.41955,  0.38808,  0.270056, 0.134572, 0.049776};

float _fullres_1_AB[9] = {0.277282, 0.449211, 0.620534, 0.726469, 0.750364, 0.717564,  0.635928, 0.52164, 0.414801};//from shaowei
float _fullres_2_AB[9] = {0.049776, 0.134572,  0.270056,  0.38808, 0.41955, 0.376866, 0.285307, 0.184801, 0.113904};


float _fullres_1_CD[9] = {0.466213, 0.610396, 0.652063, 0.604029, 0.498893, 0.375308, 0.263286, 0.186844, 0.146014};
float _fullres_2_CD[9] = {0.145522, 0.260327, 0.30193, 0.254339, 0.168, 0.0925498, 0.0448017, 0.0223924, 0.013635}; 
float _fullres_1[9];
float _fullres_2[9];

//float _fullres_1_AB_C[9] = {0.293989 ,  0.487802, 0.65368, 0.742256, 0.745738, 0.692386, 0.587542, 0.460459, 0.363872};
//float _fullres_2_AB_C[9] = {0.056076 ,  0.160155, 0.303633, 0.408638, 0.413294, 0.346554,  0.239265,0.14176, 0.0868215};

//float _fullres_1_AB_D[9] = {0.278407, 0.478096, 0.646629, 0.734567, 0.733228, 0.667934, 0.54375,  0.392492, 0.272611};
//float _fullres_2_AB_D[9] = {0.0501879, 0.153474, 0.296257, 0.398514, 0.396773, 0.318951, 0.202086, 0.101538, 0.0480855};


float _fullres_1_AB_C[9] = {0.267599, 0.442116, 0.615878, 0.722996, 0.745963, 0.710407, 0.622985, 0.502017, 0.388498};
float _fullres_1_AB_D[9] = {0.259867, 0.437249, 0.611714, 0.718302, 0.739225, 0.697993, 0.600531, 0.464815, 0.335032};

float _fullres_2_AB_C[9] = {0.046305, 0.130148, 0.265557, 0.383674, 0.413596, 0.368045, 0.272445, 0.17025, 0.0994078};
float _fullres_2_AB_D[9] = {0.0436283, 0.127163, 0.261577, 0.377785, 0.404622, 0.353133, 0.251089, 0.144602, 0.0732584};

//shaowei
//double res1_eABeDtB[9] = {0.335032, 0.464815, 0.600531, 0.697993, 0.739225, 0.718302, 0.611714, 0.437249, 0.259867,};
//double res1_eABeCFtB[9]={0.388498, 0.502017, 0.622985, 0.710407, 0.745963, 0.722996, 0.615878, 0.442116, 0.267599};

//double res12_eABeDtB[9]={0.0732584, 0.144602, 0.251089, 0.353133, 0.404622, 0.377785, 0.261577, 0.127163, 0.0436283};
//double res12_eABeCDtB[9]={0.0994078, 0.17025, 0.272445, 0.368045, 0.413596, 0.383674, 0.265557, 0.130148, 0.046305};


//shaowei
// resolution // from peripherial to central
// //EPD-AB resolution using EPD-C as reference EP  (default)
// double res1_eABeCtB[9]  = {0.363872, 0.460459, 0.587542, 0.692386, 0.745738, 0.742256, 0.65368, 0.487802, 0.293989};
// double res12_eABeCtB[9] = {0.0868215,0.14176,  0.239265, 0.346554, 0.413294, 0.408638, 0.303633,0.160155, 0.056076};
//
// //EPD-AB resolution using EPD-D as reference EP  (for systematic study)
// double res1_eABeDtB[9]  = {0.272611,  0.392492, 0.54375,  0.667934, 0.733228, 0.734567, 0.646629, 0.478096, 0.278407};
// double res12_eABeDtB[9] = {0.0480855, 0.101538, 0.202086, 0.318951, 0.396773, 0.398514, 0.296257, 0.153474, 0.0501879};
//
// //EPD-AB resolution using EPD-CD as reference EP  (for systematic study)
// double res1_eABeCDtB[9]  = {0.333295, 0.438144, 0.572763, 0.683997, 0.74161, 0.739785, 0.65093, 0.482731, 0.285382};
// double res12_eABeCDtB[9] = {0.0724813, 0.127709, 0.226264, 0.33689, 0.407779, 0.405361, 0.300741, 0.156644, 0.0527811};
//
//
//int snn = 28 ;
int snn = 3 ;

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

float ld_chi2topo;
float ld_chi2ndf;
float ld_ldl;
float ld_l;
float chi2primary_proton;
float chi2primary_pi;

int nhits_ld_pi;
int nhits_ld_proton;

float ld_chi2topo_cut;
float ld_chi2ndf_cut;
float ld_ldl_cut;
float ld_l_cut;
float chi2primary_proton_cut;
float chi2primary_pi_cut;

int fCentrality;
double psi_1_EPD_0;
double psi_1_EPD_1;
double psi_1_EPD_2;
double psi_1_EPD_3;
double psi_1_EPD_4;
double psi_1_EPD_5;
double psi_1_EPD;

float bphi;
float bphi_evp_1o;
float bphi_evp_2o;
float ycm;
float brap;
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

int FXTMult;
//float cent9;
int cent9;

//standard cut
/*
float ld_chi2topo_cut = 5     ;
float ld_chi2ndf_cut = 10     ;
float ld_ldl_cut = 5          ;
float ld_l_cut = 5            ;
float chi2primary_proton_cut = 18.6 ;
float chi2primary_pi_cut = 18.6     ;
*/
//end standard cut
int    _condor;
int    _batch;
int    _cut;
Char_t path_ch[10000];
Char_t ofile[10000];

float bdt_ld_cut[ncenttmvabins][nptbins];
int _evt_plane;//0 = AB, 1 = CD
bool _mb;

void book_histo();
void write_histo();
void process_tree_tmva_ld_kf_vn_v15(int condor=0,int evt_plane=0, int batch=0, int cut=1){

    _condor = condor;//condor=0: mc, condor=1: data,
    _evt_plane = evt_plane;
    _cut = cut;
    _batch = batch;
  //  _mb = true;
    _mb = false;

    if(_evt_plane==0){
      for(int icent=0;icent<9;icent++){
         _fullres_1[icent] = _fullres_1_AB[icent];
         _fullres_2[icent] = _fullres_2_AB[icent];
	if(_cut==31){
	_fullres_1[icent] = _fullres_1_AB_C[icent];
        _fullres_2[icent] = _fullres_2_AB_C[icent];
	}
	if(_cut==32){
        _fullres_1[icent] = _fullres_1_AB_D[icent];
        _fullres_2[icent] = _fullres_2_AB_D[icent];
        }
      }
    }else{
      for(int icent=0;icent<9;icent++){
         _fullres_1[icent] = _fullres_1_CD[icent];
         _fullres_2[icent] = _fullres_2_CD[icent];
      }
    }

//cuts

// TFile  *fgpt_eff = new TFile(Form("ld_vn_eff_v0_cut%d.root",_cut),"READ");
 TFile  *fgpt_eff = new TFile(Form("ld_vn_eff_v0_cut%d_v15.root",_cut),"READ");
 h_pt_eff_wgt_fine2 = (TH2F*)fgpt_eff->Get("h_pt_eff_wgt_fine2")->Clone("h_pt_eff_wgt_fine2");

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


TChain lambda_tree("lambda_tree");

int nChains = -999999;
ifstream fin;

if(condor!=0){
if(_batch==0){
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ld_ana_flow_vn_loose/filelist.txt");
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ld_ana_flow_vn_loose_newdata/filelist.txt");
fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ld_ana_flow_vn_loose_newdata_xcheck//filelist.txt");
if(!_mb){
strcpy (ofile,Form("run18_3gev_ld_kf_flow_new_weighted_evt%d_cut%d_v15.root",_evt_plane,_cut));
}else{
strcpy (ofile,Form("run18_3gev_ld_kf_flow_new_weighted_evt%d_cut%d_v15_mb.root",_evt_plane,_cut));
}
nChains = 1225;
}else{
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ld_ana_flow_vn_rotate_loose/filelist.txt");
//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ld_ana_flow_vn_rotate_loose_newdata/filelist.txt");
fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ld_ana_flow_vn_rotate_loose_newdata_xcheck/filelist.txt");
if(!_mb){
strcpy (ofile,Form("run18_3gev_ld_kf_flow_rotate_new_weighted_evt%d_cut%d_v15.root",_evt_plane,_cut));
}else{
strcpy (ofile,Form("run18_3gev_ld_kf_flow_rotate_new_weighted_evt%d_cut%d_v15_mb.root",_evt_plane,_cut));
}
nChains = 1224;
}
}else{
fin.open("filelist.txt");
strcpy (ofile,"flow_hist.root");
nChains = 1;
}


for(int i=0; i<nChains; i++){
      fin >> path_ch;
      lambda_tree.Add(path_ch);
      cout<<path_ch<<"  "<<lambda_tree.GetEntries() <<endl;

}//ichain
cout << " done TChain-ing.. " << endl;

book_histo();


lambda_tree.SetBranchAddress("brefmult",&brefmult);
lambda_tree.SetBranchAddress("bparticleid",&bparticleid);
lambda_tree.SetBranchAddress("bparticlemass",&bparticlemass);
lambda_tree.SetBranchAddress("bpx",&bpx);
lambda_tree.SetBranchAddress("bpy",&bpy);
lambda_tree.SetBranchAddress("bpz",&bpz);

lambda_tree.SetBranchAddress("nhits_ld_pi",&nhits_ld_pi);
lambda_tree.SetBranchAddress("nhits_ld_proton",&nhits_ld_proton);

lambda_tree.SetBranchAddress("cent9",&cent9);
lambda_tree.SetBranchAddress("ld_chi2topo", &ld_chi2topo);
lambda_tree.SetBranchAddress("ld_chi2ndf", &ld_chi2ndf);
lambda_tree.SetBranchAddress("ld_ldl", &ld_ldl);
lambda_tree.SetBranchAddress("ld_l", &ld_l);
lambda_tree.SetBranchAddress("chi2primary_proton", &chi2primary_proton);
lambda_tree.SetBranchAddress("chi2primary_pi", &chi2primary_pi);
lambda_tree.SetBranchAddress("fCentrality", &fCentrality);
lambda_tree.SetBranchAddress("gweight", &gweight);
lambda_tree.SetBranchAddress("FXTMult", &FXTMult);

//lambda_tree.SetBranchAddress("psi_1_EPD_0", &psi_1_EPD_0);
//lambda_tree.SetBranchAddress("psi_1_EPD_1", &psi_1_EPD_1);
//lambda_tree.SetBranchAddress("psi_1_EPD_2", &psi_1_EPD_2);
//lambda_tree.SetBranchAddress("psi_1_EPD_3", &psi_1_EPD_3);

//lambda_tree.SetBranchAddress("psi_1_EPD", &psi_1_EPD);
if(_evt_plane==0){lambda_tree.SetBranchAddress("psi_1_EPD_4", &psi_1_EPD);}
if(_evt_plane==1){lambda_tree.SetBranchAddress("psi_1_EPD_5", &psi_1_EPD);}

 Long64_t n_lambda_Entries = lambda_tree.GetEntries();


 cout<<endl<<"Start processing "<<endl<<endl;
 cout<<"--------------------------------------------"<<endl;
 cout<<"Total entries "<<n_lambda_Entries<< endl;

if(snn==3){
//ycm = -1.05;
ycm = -1.045;
}

if(_cut==1){
ld_chi2topo_cut = 5;
ld_chi2ndf_cut = 10;
ld_ldl_cut = 5;
ld_l_cut = 5;
chi2primary_proton_cut = 18.6;
chi2primary_pi_cut = 18.6;
}

if(_cut==2 || _cut==11 || _cut==12 || _cut==21 || _cut==22 || _cut==31 || _cut==32){
ld_chi2topo_cut = 5;
ld_chi2ndf_cut = 5;
ld_ldl_cut = 5;
ld_l_cut = 1;
chi2primary_proton_cut = 10.;
chi2primary_pi_cut = 10.;
}
if(_cut==3){
ld_chi2topo_cut = 5;
ld_chi2ndf_cut = 5;
ld_ldl_cut = 5;
ld_l_cut = 1;
chi2primary_proton_cut = 10.;
chi2primary_pi_cut = 5.;
}
if(_cut==4){
ld_chi2topo_cut = 5;
ld_chi2ndf_cut = 5;
ld_ldl_cut = 5;
ld_l_cut = 1;
chi2primary_proton_cut = 10.;
chi2primary_pi_cut = 15.;
}
if(_cut==5){
ld_chi2topo_cut = 5;
ld_chi2ndf_cut = 5;
ld_ldl_cut = 5;
ld_l_cut = 1;
chi2primary_proton_cut = 5.;
chi2primary_pi_cut = 10.;
}
if(_cut==6){
ld_chi2topo_cut = 5;
ld_chi2ndf_cut = 5;
ld_ldl_cut = 5;
ld_l_cut = 1;
chi2primary_proton_cut = 15.;
chi2primary_pi_cut = 10.;
}
if(_cut==7){
ld_chi2topo_cut = 4;
ld_chi2ndf_cut = 5;
ld_ldl_cut = 5;
ld_l_cut = 1;
chi2primary_proton_cut = 10.;
chi2primary_pi_cut = 10.;
}
if(_cut==8){
ld_chi2topo_cut = 3;
ld_chi2ndf_cut = 5;
ld_ldl_cut = 5;
ld_l_cut = 1;
chi2primary_proton_cut = 10.;
chi2primary_pi_cut = 10.;
}
if(_cut==9){
ld_chi2topo_cut = 5;
ld_chi2ndf_cut = 4;
ld_ldl_cut = 5;
ld_l_cut = 1;
chi2primary_proton_cut = 10.;
chi2primary_pi_cut = 10.;
}
if(_cut==10){
ld_chi2topo_cut = 5;
ld_chi2ndf_cut = 3;
ld_ldl_cut = 5;
ld_l_cut = 1;
chi2primary_proton_cut = 10.;
chi2primary_pi_cut = 10.;
}

  for (Long64_t iEntry = 0; iEntry<=n_lambda_Entries; iEntry++){

	lambda_tree.GetEntry(iEntry);
  if(iEntry%100000==0){
    cout<<"Processing entry:" <<iEntry<<endl;
}
int centtmva_label = -1;
int cent_label = -1;
int pt_label = -1;

int cent_cut_label = -1;
int pt_cut_label = -1;
int rap_label = -1;
int rap_label_2 = -1;
float eff_weight = 0.;


if(snn==3){

//use zach: https://drupal.star.bnl.gov/STAR/system/files/Sweger_3p0GeV_StandardNewest_fcv2020Nov11.pdf
//fCentrality 8 most central, 0 most peripheral
//7,8: 0-10, 4-6:10-40
if(FXTMult>195) fCentrality= 9;
else if(FXTMult>141) fCentrality=8;
else if(FXTMult>118) fCentrality=7;
else if(FXTMult>85) fCentrality=6;
else if(FXTMult>59) fCentrality=5;
else if(FXTMult>40) fCentrality=4;
else if(FXTMult>25) fCentrality=3;
else if(FXTMult>15) fCentrality=2;
else if(FXTMult>8) fCentrality=1;
else fCentrality=0;



   cent_label = fCentrality;

//0,10,40,60
if(fCentrality==7||fCentrality==8) cent_cut_label = 0;
else if(fCentrality==4||fCentrality==5||fCentrality==6) cent_cut_label = 1;
else if(fCentrality==2||fCentrality==3) cent_cut_label = 2;
else cent_cut_label = -1;

}

  if(cent_cut_label<0) continue;

	if(_mb) cent_cut_label = 0;

  if(bparticlemass>1.2) continue;


	if(chi2primary_proton<chi2primary_proton_cut || chi2primary_pi<chi2primary_pi_cut) continue;
        if(ld_chi2topo>ld_chi2topo_cut) continue;
        if(ld_chi2ndf>ld_chi2ndf_cut) continue;
	if(ld_l<ld_l_cut) continue;
	if(ld_ldl<ld_ldl_cut) continue;

	if(_cut==11 && nhits_ld_pi<17) continue;
        if(_cut==11 && nhits_ld_proton<17) continue;

        if(_cut==12 && nhits_ld_pi<20) continue;
        if(_cut==12 && nhits_ld_proton<20) continue;

  bpt = sqrt(bpx*bpx+bpy*bpy);
  bphi = atan2(bpy,bpx);
  bphi_evp_1o = bphi - psi_1_EPD; 
  bphi_evp_2o = bphi - psi_1_EPD; 

  if(bphi_evp_1o<-1.0*TMath::Pi()) bphi_evp_1o += 2*TMath::Pi();
  if(bphi_evp_1o>1.0*TMath::Pi())  bphi_evp_1o -= 2*TMath::Pi();
  bphi_evp_1o = fabs(bphi_evp_1o);
 
  if(bphi_evp_2o<-1.0*TMath::Pi()) bphi_evp_2o += 2*TMath::Pi();
  if(bphi_evp_2o>1.0*TMath::Pi())  bphi_evp_2o -= 2*TMath::Pi();
  bphi_evp_2o = fabs(bphi_evp_2o);
  if(bphi_evp_2o>0.5*TMath::Pi()) bphi_evp_2o = TMath::Pi() - bphi_evp_2o;  

  TLorentzVector ptc(bpx,bpy,bpz,sqrt(bpx*bpx+bpy*bpy+bpz*bpz+bparticlemass*bparticlemass));
  brap = ptc.Rapidity() - ycm;

for(int ipt=0; ipt<ncutptbins; ipt++){
  if(bpt > ptcutbin[ipt]) pt_cut_label = ipt;
}
  if(bpt > ptcutbin[ncutptbins]) pt_cut_label = -1;

for(int irap=0; irap<nrapbins; irap++){
  if(brap > rapbin[irap]) rap_label = irap;
}
  if(brap > rapbin[nrapbins]) rap_label = -1;

for(int irap=0; irap<nrapbins_2; irap++){
  if(brap > rapbin_2[irap]) rap_label_2 = irap;
}
  if(brap > rapbin_2[nrapbins_2]) rap_label_2 = -1;

if(bparticleid>0 && bparticlemass>lambda_mass-0.005 && bparticlemass<lambda_mass+0.005 ){
  h_om_pt_rap[cent_cut_label]->Fill(brap,bpt);
  h_r_om_pt_rap[cent_cut_label]->Fill(-brap,bpt);
}

  if(pt_cut_label < 0) continue;//is this necessary?
  eff_weight = h_pt_eff_wgt_fine2->GetBinContent(h_pt_eff_wgt_fine2->FindBin(brap,bpt));
  if(eff_weight!=0)  eff_weight = 1./eff_weight;
if(bparticleid>0){

  h_om_mass[cent_label][pt_cut_label][rap_label]->Fill(bparticlemass);
  h_om_phi_evp[cent_cut_label][pt_cut_label]->Fill(psi_1_EPD);

if(rap_label_2>=0){ 
  h_om_mass_phi_2o[cent_cut_label][pt_cut_label][rap_label_2]->Fill(bparticlemass,bphi_evp_2o,1./_fullres_2[8-fCentrality]);
  h_om_mass_phi_2o_eff[cent_cut_label][pt_cut_label][rap_label_2]->Fill(bparticlemass,bphi_evp_2o,1./_fullres_2[8-fCentrality]*eff_weight);
}

if(rap_label>=0){
  h_om_mass_phi_1o[cent_cut_label][pt_cut_label][rap_label]->Fill(bparticlemass,bphi_evp_1o,1./_fullres_1[8-fCentrality]);
  h_om_mass_phi_1o_eff[cent_cut_label][pt_cut_label][rap_label]->Fill(bparticlemass,bphi_evp_1o,1./_fullres_1[8-fCentrality]*eff_weight);
	}
}

}//entry loop


write_histo();

}

void book_histo()
{



for(int icent=0;icent<ncentbins;icent++){
  for(int ipt=0;ipt<ncutptbins;ipt++){
	for(int irap=0;irap<nrapbins;irap++){
h_om_mass[icent][ipt][irap] = new TH1F(Form("h_om_mass[%d][%d][%d]",icent,ipt,irap), Form("h_om_mass[%d][%d][%d]",icent,ipt,irap), 1000,1,2);
}
}
}

for(int icent=0;icent<ncutcentbins;icent++){

h_om_pt_rap[icent] =  new TH2F(Form("h_om_pt_rap[%d]",icent),"",500,-2,2,500,0,5);
h_r_om_pt_rap[icent] =  new TH2F(Form("h_r_om_pt_rap[%d]",icent),"",500,-2,2,500,0,5);

  for(int ipt=0;ipt<ncutptbins;ipt++){
 //h_om_mass[icent][ipt] = new TH1F(Form("h_om_mass[%d][%d]",icent,ipt), Form("h_om_mass[%d][%d]",icent,ipt), 1000,1,2);
 h_om_phi_evp[icent][ipt] = new TH1F(Form("h_om_phi_evp[%d][%d]",icent,ipt), Form("h_om_phi_evp[%d][%d]",icent,ipt), 1000,-10,10);
// h_om_mass_phi_2o[icent][ipt] = new TH2F(Form("h_om_mass_phi_evp_2o[%d][%d]",icent,ipt), Form("h_om_mass_phi_evp_2o[%d][%d]",icent,ipt), 1000,1,2, 12, 0 ,0.5*TMath::Pi());

for(int irap=0;irap<nrapbins;irap++){
  h_om_mass_phi_1o[icent][ipt][irap] = new TH2F(Form("h_om_mass_phi_evp_1o[%d][%d][%d]",icent,ipt,irap), Form("h_om_mass_phi_evp_1o[%d][%d][%d]",icent,ipt,irap), 1000,1,2, 12, 0 ,1.0*TMath::Pi());
  h_om_mass_phi_1o_eff[icent][ipt][irap] = new TH2F(Form("h_om_mass_phi_evp_1o_eff[%d][%d][%d]",icent,ipt,irap), Form("h_om_mass_phi_evp_1o_eff[%d][%d][%d]",icent,ipt,irap), 1000,1,2, 12, 0 ,1.0*TMath::Pi());
  h_om_mass_phi_1o[icent][ipt][irap]->Sumw2();
  h_om_mass_phi_1o_eff[icent][ipt][irap]->Sumw2();
}
for(int irap=0;irap<nrapbins_2;irap++){
  h_om_mass_phi_2o[icent][ipt][irap] = new TH2F(Form("h_om_mass_phi_evp_2o[%d][%d][%d]",icent,ipt,irap), Form("h_om_mass_phi_evp_2o[%d][%d][%d]",icent,ipt,irap), 1000,1,2, 12, 0 ,0.5*TMath::Pi());
  h_om_mass_phi_2o_eff[icent][ipt][irap] = new TH2F(Form("h_om_mass_phi_evp_2o_eff[%d][%d][%d]",icent,ipt,irap), Form("h_om_mass_phi_evp_2o_eff[%d][%d][%d]",icent,ipt,irap), 1000,1,2, 12, 0 ,0.5*TMath::Pi());
  h_om_mass_phi_2o[icent][ipt][irap]->Sumw2();
  h_om_mass_phi_2o_eff[icent][ipt][irap]->Sumw2();
}

  }
}

}

void write_histo()
{
    TFile *outhistfile = new TFile (ofile, "RECREATE");
    outhistfile->cd();


   // hvtx -> Write();
   // hvtxgood -> Write();
   // hrefmult -> Write();

   for(int icent=0;icent<ncentbins;icent++){
     for(int ipt=0;ipt<ncutptbins;ipt++){
       for(int irap=0;irap<nrapbins;irap++){
	h_om_mass[icent][ipt][irap]->Write();
		}
		}
		}

cout<<"ncutcentbins:"<<ncutcentbins<<" "<<ncutptbins<<endl;
   for(int icent=0;icent<ncutcentbins;icent++){
h_om_pt_rap[icent] -> Write();
h_r_om_pt_rap[icent] -> Write();
     for(int ipt=0;ipt<ncutptbins;ipt++){
//          h_om_mass[icent][ipt]->Write();
//	  h_om_mass_phi_2o[icent][ipt]->Write();
	  h_om_phi_evp[icent][ipt]->Write();
     for(int irap=0;irap<nrapbins;irap++){
	  h_om_mass_phi_1o[icent][ipt][irap]->Write();
          h_om_mass_phi_1o_eff[icent][ipt][irap]->Write();
}
for(int irap=0;irap<nrapbins_2;irap++){
          h_om_mass_phi_2o[icent][ipt][irap]->Write();
          h_om_mass_phi_2o_eff[icent][ipt][irap]->Write();
}

/*
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

     if(_condor==1){
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
*/
     }
   }

}
