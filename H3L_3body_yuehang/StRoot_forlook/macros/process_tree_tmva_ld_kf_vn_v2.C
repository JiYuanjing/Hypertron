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

//TH2F *h_om_pt_rap;

TH2F *h_phieta;
TH2F *h_phieta_bg;

double ptbin[nptbins+1] = {0.4,0.8,1.2,1.6,2.0,2.4,2.8};

//const int ncutptbins = 3;
//double ptcutbin[ncutptbins+1] = {0.0,1.0,2.0,3.0};
//const int ncutcentbins = 3;
//int centcutbin[ncutcentbins+1] = {0,10,40,80};
const int ncutcentbins = 4;
int centcutbin[ncutcentbins+1] = {0,10,20,40,60};

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
TH1F *h_om_phi_evp[ncutcentbins][ncutptbins];

int centtmva_limit[ncenttmvabins+1] = {0,10,40,80};
int cent_limit[ncentbins+1] = {0,5,10,20,30,40,50,60,70,80};

float _fullres_1_AB[9] = {0.287265,0.481901,0.649618,0.738704,0.740876,0.683774,0.575272,0.446576,0.346197};
float _fullres_2_AB[9] = {0.0534926,0.156073,0.299369,0.403934,0.406805,0.336636,0.228438,0.132919,0.0783611};
float _fullres_1_CD[9] = {0.466213, 0.610396, 0.652063, 0.604029, 0.498893, 0.375308, 0.263286, 0.186844, 0.146014};
float _fullres_2_CD[9] = {0.145522, 0.260327, 0.30193, 0.254339, 0.168, 0.0925498, 0.0448017, 0.0223924, 0.013635}; 
float _fullres_1[9];
float _fullres_2[9];

//int cent_limit[ncentbins+1] = {0,5,10,20,30,40,50,60,80};
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

float ld_chi2topo;
float ld_chi2ndf;
float ld_ldl;
float ld_l;
float chi2primary_proton;
float chi2primary_pi;
int fCentrality;
double psi_1_EPD_0;
double psi_1_EPD_1;
double psi_1_EPD_2;
double psi_1_EPD_3;
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

//float cent9;
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

//TMVA::Reader* reader[ncenttmvabins][npttmvabins];
// float bdt_ld_cut[ncenttmvabins][npttmvabins];
float bdt_ld_cut[ncenttmvabins][nptbins];
int _evt_plane;//0 = AB, 1 = CD

void book_histo();
void write_histo();
void process_tree_tmva_ld_kf_vn_v2(int condor=0,int evt_plane=0){
    _condor = condor;//condor=0: mc, condor=1: data,
    _evt_plane = evt_plane;

    if(_evt_plane==0){
      for(int icent=0;icent<9;icent++){
         _fullres_1[icent] = _fullres_1_AB[icent];
         _fullres_2[icent] = _fullres_2_AB[icent];
      }
    }else{
      for(int icent=0;icent<9;icent++){
         _fullres_1[icent] = _fullres_1_CD[icent];
         _fullres_2[icent] = _fullres_2_CD[icent];
      }
    }

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


TChain lambda_tree("lambda_tree");

int nChains = -999999;
ifstream fin;

if(condor!=0){
fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ld_ana_flow_vn_new/filelist.txt");
strcpy (ofile,Form("run18_3gev_ld_kf_flow_new_weighted_evt%d.root",_evt_plane));

//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ld_ana_flow_vn_new_rotate/filelist.txt");
//strcpy (ofile,Form("run18_3gev_ld_kf_flow_rotate_new_weighted_evt%d.root",_evt_plane));
nChains = 1222;
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

 gitom = new TF1("gitom","x*3.1415962*2*[0]*sqrt(1.115683*1.115683+x*x) * exp(-sqrt(1.115683*1.115683+x*x)/[1])", 0, 6);
  gitom->SetParameters(1.48179e+00,3.07372e-01);

lambda_tree.SetBranchAddress("brefmult",&brefmult);
lambda_tree.SetBranchAddress("bparticleid",&bparticleid);
lambda_tree.SetBranchAddress("bparticlemass",&bparticlemass);
lambda_tree.SetBranchAddress("bpx",&bpx);
lambda_tree.SetBranchAddress("bpy",&bpy);
lambda_tree.SetBranchAddress("bpz",&bpz);

lambda_tree.SetBranchAddress("cent9",&cent9);
lambda_tree.SetBranchAddress("ld_chi2topo", &ld_chi2topo);
lambda_tree.SetBranchAddress("ld_chi2ndf", &ld_chi2ndf);
lambda_tree.SetBranchAddress("ld_ldl", &ld_ldl);
lambda_tree.SetBranchAddress("ld_l", &ld_l);
lambda_tree.SetBranchAddress("chi2primary_proton", &chi2primary_proton);
lambda_tree.SetBranchAddress("chi2primary_pi", &chi2primary_pi);
lambda_tree.SetBranchAddress("fCentrality", &fCentrality);
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

if(snn==28){//temp

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

if(snn==3){
   cent_label = fCentrality;

//0,10,40,80
/*
if(fCentrality==7||fCentrality==8) cent_cut_label = 0;
else if(fCentrality==4||fCentrality==5||fCentrality==6) cent_cut_label = 1;
else if(fCentrality==0||fCentrality==1||fCentrality==2||fCentrality==3) cent_cut_label = 2;
else cent_cut_label = -1;
*/

//0,10,20,40,60
if(fCentrality==7||fCentrality==8) cent_cut_label = 0;
else if(fCentrality==6) cent_cut_label = 1;
else if(fCentrality==4||fCentrality==5) cent_cut_label = 2;
else if(fCentrality==2||fCentrality==3) cent_cut_label = 3;
else cent_cut_label = -1;
}
  if(cent_cut_label<0) continue;
  if(bparticlemass>1.2) continue;

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
//  if(fabs(ptc.Rapidity())>0.5) continue;
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
}

  if(pt_cut_label < 0) continue;//is this necessary?

if(bparticleid>0){

  h_om_mass[cent_label][pt_cut_label][rap_label]->Fill(bparticlemass);
  h_om_phi_evp[cent_cut_label][pt_cut_label]->Fill(psi_1_EPD);

if(rap_label_2>=0){ 
  //h_om_mass_phi_2o[cent_cut_label][pt_cut_label][rap_label_2]->Fill(bparticlemass,bphi_evp_2o);
  h_om_mass_phi_2o[cent_cut_label][pt_cut_label][rap_label_2]->Fill(bparticlemass,bphi_evp_2o,1./_fullres_2[8-fCentrality]);
}

if(rap_label>=0){
  //h_om_mass_phi_1o[cent_cut_label][pt_cut_label][rap_label]->Fill(bparticlemass,bphi_evp_1o);
  h_om_mass_phi_1o[cent_cut_label][pt_cut_label][rap_label]->Fill(bparticlemass,bphi_evp_1o,1./_fullres_1[8-fCentrality]);
	}
}



/*
if(bparticlemass<omega_mass+3*sgm3 && bparticlemass>omega_mass-3*sgm3){
  if(_condor==0){

 h_nhits_om_proton[cent_cut_label][pt_cut_label] -> Fill(nhits_om_proton,_weight);
 h_nhits_om_pi[cent_cut_label][pt_cut_label] -> Fill(nhits_om_pi,_weight);
 h_nhits_om_bach[cent_cut_label][pt_cut_label] -> Fill(nhits_om_bach,_weight);
 h_dedx_om_proton[cent_cut_label][pt_cut_label] -> Fill(dedx_om_proton,_weight);
 h_dedx_om_pi[cent_cut_label][pt_cut_label] -> Fill(dedx_om_pi,_weight);
 h_dedx_om_bach[cent_cut_label][pt_cut_label] -> Fill(dedx_om_bach,_weight);

  }else{

 h_nhits_om_proton[cent_cut_label][pt_cut_label] -> Fill(nhits_om_proton);
 h_nhits_om_pi[cent_cut_label][pt_cut_label] -> Fill(nhits_om_pi);
 h_nhits_om_bach[cent_cut_label][pt_cut_label] -> Fill(nhits_om_bach);
 h_dedx_om_proton[cent_cut_label][pt_cut_label] -> Fill(dedx_om_proton);
 h_dedx_om_pi[cent_cut_label][pt_cut_label] -> Fill(dedx_om_pi);
 h_dedx_om_bach[cent_cut_label][pt_cut_label] -> Fill(dedx_om_bach);
  }
}
*/
/*
if(  (bparticlemass<omega_mass+6*sgm3 && bparticlemass>omega_mass+3*sgm3) ||  (bparticlemass>omega_mass-6*sgm3 && bparticlemass<omega_mass-3*sgm3)  ){
if(_condor==1){
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
*/
}//entry loop


write_histo();

}

void book_histo()
{


//h_om_pt_rap =  new TH2F("h_om_pt_rap","h_om_pt_rap",500,-2,2,500,0,5);

for(int icent=0;icent<ncentbins;icent++){
  for(int ipt=0;ipt<ncutptbins;ipt++){
	for(int irap=0;irap<nrapbins;irap++){
h_om_mass[icent][ipt][irap] = new TH1F(Form("h_om_mass[%d][%d][%d]",icent,ipt,irap), Form("h_om_mass[%d][%d][%d]",icent,ipt,irap), 1000,1,2);
}
}
}

for(int icent=0;icent<ncutcentbins;icent++){

h_om_pt_rap[icent] =  new TH2F(Form("h_om_pt_rap[%d]",icent),"",500,-2,2,500,0,5);

  for(int ipt=0;ipt<ncutptbins;ipt++){
 //h_om_mass[icent][ipt] = new TH1F(Form("h_om_mass[%d][%d]",icent,ipt), Form("h_om_mass[%d][%d]",icent,ipt), 1000,1,2);
 h_om_phi_evp[icent][ipt] = new TH1F(Form("h_om_phi_evp[%d][%d]",icent,ipt), Form("h_om_phi_evp[%d][%d]",icent,ipt), 1000,-10,10);
// h_om_mass_phi_2o[icent][ipt] = new TH2F(Form("h_om_mass_phi_evp_2o[%d][%d]",icent,ipt), Form("h_om_mass_phi_evp_2o[%d][%d]",icent,ipt), 1000,1,2, 12, 0 ,0.5*TMath::Pi());

for(int irap=0;irap<nrapbins;irap++){
  h_om_mass_phi_1o[icent][ipt][irap] = new TH2F(Form("h_om_mass_phi_evp_1o[%d][%d][%d]",icent,ipt,irap), Form("h_om_mass_phi_evp_1o[%d][%d][%d]",icent,ipt,irap), 1000,1,2, 12, 0 ,1.0*TMath::Pi());
  h_om_mass_phi_1o[icent][ipt][irap]->Sumw2();
}
for(int irap=0;irap<nrapbins_2;irap++){
  h_om_mass_phi_2o[icent][ipt][irap] = new TH2F(Form("h_om_mass_phi_evp_2o[%d][%d][%d]",icent,ipt,irap), Form("h_om_mass_phi_evp_2o[%d][%d][%d]",icent,ipt,irap), 1000,1,2, 12, 0 ,0.5*TMath::Pi());
  h_om_mass_phi_2o[icent][ipt][irap]->Sumw2();
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
     for(int ipt=0;ipt<ncutptbins;ipt++){
//          h_om_mass[icent][ipt]->Write();
//	  h_om_mass_phi_2o[icent][ipt]->Write();
	  h_om_phi_evp[icent][ipt]->Write();
     for(int irap=0;irap<nrapbins;irap++){
	  h_om_mass_phi_1o[icent][ipt][irap]->Write();
}
for(int irap=0;irap<nrapbins_2;irap++){
          h_om_mass_phi_2o[icent][ipt][irap]->Write();
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
