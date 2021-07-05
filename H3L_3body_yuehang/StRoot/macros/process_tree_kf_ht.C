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

const int ncentbins = 9;
int    brunid          ;
int    beventid        ;
int    brefmult        ;
int    btofmult        ;
int    bismc           ;
float  bVz             ;
int    bparticleid     ;
float  bparticlemass   ;
float  bpx             ;
float  bpy             ;
float  bpz             ;
float  bpl             ;

float ht_chi2ndf        ;
float ht_ldl            ;
float ht_chi2topo       ;
float refMult           ;
float pt                ;

float ht_mass;
float sgm3;
int centrality;
float chi2primary_pi;
float chi2primary_he;
float ht_l;
float ht_pl;
float ht_chi2topo_cut = 5     ;
float ht_chi2ndf_cut = 10     ;
float ht_ldl_cut = 5          ;
float ht_l_cut = 5            ;
float chi2primary_he_cut = 18.6 ;
float chi2primary_pi_cut = 18.6 ;
float bmcpx;
float bmcpy;
float bmcpz;
float bmcrap;
float ycm;
float  bpt             ;
float  bp              ;
float  beta            ;
float  brap            ;
int countrefmult;
TFile  *fgpt_0;
TH2F *g_pt_fine_in;
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

double weight;
int cent_label;

Char_t path_ch[10000];
Char_t ofile[10000];

TFile *outhistfile;// = new TFile (ofile, "RECREATE");

const int nptbins = 1;
TTree *htriton_tree_pt[nptbins];
double ptbin[nptbins+1] =    {6.0,50.0};

bool _data_or_sim;
float r_cent[ncentbins+1] = {0., 0., 0.178389, 0.466883, 1.06469, 2.30611, 2.70053, 1.49454, 0.866928,0};//h3l fitting
int centFull[ncentbins] ={15,22,32,43,57,73,92,117,133};
void book_tree();
void write_tree();
int Centrality(int aa );
void process_tree_kf_ht(int aa);
void process_tree_kf_ht(int data_or_sim){

	_data_or_sim = data_or_sim;//0:data, comb
	//_data_or_sim =1;//1:sim, emb

	levyfit4 = new TF1("levyfit4","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.4,2.4);
        levyfit5 = new TF1("levyfit5","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.2,2.4);
        levyfit6 = new TF1("levyfit6","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.2,2.4);
        levyfit7 = new TF1("levyfit7","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.2,2.4);
        levyfit8 = new TF1("levyfit8","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.2,2.4);
        t_quadr = new TF1("t_quadr","[0]+[1]*x+[2]*x*x",-2,2);

        t_quadr0 = new TF1("t_quadr0","[0]+[1]*x+[2]*x*x",-2,2);
        t_quadr1 = new TF1("t_quadr1","[0]+[1]*x+[2]*x*x",-2,2);
//        bolt0 = new TF1("bolt0", "x*sqrt(2.80941350568*2.80941350568+x*x)*exp(-sqrt(2.80941350568*2.80941350568+x*x)/[0])", 0.,5);
//        bolt1 = new TF1("bolt1", "x*sqrt(3.72840015733*3.72840015733+x*x)*exp(-sqrt(3.72840015733*3.72840015733+x*x)/[0])", 0.,5);
	bolt0 = new TF1("bolt0", "1e9*x*sqrt(2.990*2.990+x*x)*exp(-sqrt(2.990*2.990+x*x)/[0])", 0.,5);
        bolt1 = new TF1("bolt1", "1e9*x*sqrt(2.990*2.990+x*x)*exp(-sqrt(2.990*2.990+x*x)/[0])", 0.,5);
        bolt2 = new TF1("bolt2", "1e9*x*sqrt(2.990*2.990+x*x)*exp(-sqrt(2.990*2.990+x*x)/[0])", 0.,5);


        levyfit4->SetParameters(1.53536e-01 , -1.21064e+08, 5.25287e+07);
        levyfit5->SetParameters(1.51332e-01 , -1.25927e+08, 4.86193e+07);
        levyfit6->SetParameters(1.45903e-01 , -1.60954e+08, 4.74115e+07);
        levyfit7->SetParameters(1.28382e-01 , -1.60610e+08, 4.92025e+07);
        levyfit8->SetParameters(1.11371e-01 , -2.04689e+08, 4.22255e+07);
        //t_quadr->SetParameters(2.06539, 0.0, -1.75299);
        //t_quadr0->SetParameters(3.04528e+01, 0.0, 6.71347e+01);
        //t_quadr1->SetParameters(6.88265e+01, 0.0, 1.90472e+02);
        //bolt0->SetParameter(0,0.407);
	//bolt1->SetParameter(0,0.398);

        t_quadr->SetParameters(4.71255e-01,0.00000e+00,1.79971e+00*1.5);//this is from fit to 3gevv data, multiply second paramter by 1.5
        t_quadr0->SetParameters(4.71255e-01,0.00000e+00,1.79971e+00);//this is from fit to 3gevv data
	t_quadr1->SetParameters(6.88265e+01, 0.0, 1.90472e+02);
        bolt0->SetParameter(0,0.40);
        bolt1->SetParameter(0,0.20);
        bolt2->SetParameter(0,0.15);

        fgpt_0 = new TFile("ht_input_mc_fine.root","READ");
        g_pt_fine_in = (TH2F*)fgpt_0->Get("g_pt_fine")->Clone("g_pt_fine_in");
	

	ht_mass = 2.9926;
	sgm3 = 0.00156;
	ycm = -1.045;

	ht_chi2topo_cut = 5;
        ht_chi2ndf_cut = 5;
	chi2primary_he_cut = 3.0;
        chi2primary_pi_cut = 3.0;

	TChain htriton_tree("htriton_tree");
	int nChains = -999999;
	ifstream fin;

	if(_data_or_sim==0){//sim, emb
		fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_mu_newpid2_manualvtxerr/filelist.txt");
		strcpy (ofile,"run18_3gev_tmva_ht_sig.root");
		nChains = 1107;
	}
	if(_data_or_sim==1){
		//rotation
		fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_pico_tree_newpid2_rotate_0/filelist.txt");
		strcpy (ofile,"run18_3gev_tmva_ht_bkg.root");
		nChains = 1208;

		//sideband
		//fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_pico_tree_newpid_0/filelist.txt");
                //strcpy (ofile,"run18_3gev_tmva_ht_sideband_bkg.root");
                //nChains = 1205;
	}

	for(int i=0; i<nChains; i++){
		fin >> path_ch;
		htriton_tree.Add(path_ch);
		cout<<path_ch<<"  "<<htriton_tree.GetEntries() <<endl;
	}

		cout << " done TChain-ing.. " << endl;


	htriton_tree.SetBranchAddress("bpx",&bpx);
        htriton_tree.SetBranchAddress("bpy",&bpy);
        htriton_tree.SetBranchAddress("bpz",&bpz);
	htriton_tree.SetBranchAddress("chi2primary_pi",&chi2primary_pi);
	htriton_tree.SetBranchAddress("chi2primary_he",&chi2primary_he);
	htriton_tree.SetBranchAddress("ht_l",&ht_l);
	htriton_tree.SetBranchAddress("ht_chi2ndf",&ht_chi2ndf);
	htriton_tree.SetBranchAddress("ht_ldl",&ht_ldl);
	htriton_tree.SetBranchAddress("ht_chi2topo",&ht_chi2topo);
	//htriton_tree.SetBranchAddress("refMult",&refMult);
	htriton_tree.SetBranchAddress("bparticlemass",&bparticlemass);
	htriton_tree.SetBranchAddress("bparticleid",&bparticleid);
	htriton_tree.SetBranchAddress("bpl",&bpl);
	htriton_tree.SetBranchAddress("bismc",&bismc);
        htriton_tree.SetBranchAddress("countrefmult",&countrefmult);
	if(_data_or_sim==0){
	htriton_tree.SetBranchAddress("bmcpx",&bmcpx);
	htriton_tree.SetBranchAddress("bmcpy",&bmcpy);
	htriton_tree.SetBranchAddress("bmcpz",&bmcpz);
	}

	outhistfile = new TFile (ofile, "RECREATE");
        outhistfile->cd();
	book_tree();

	Long64_t n_htriton_Entries = htriton_tree.GetEntries();

	cout<<endl<<"Start processing "<<endl<<endl;
	cout<<"--------------------------------------------"<<endl;
	cout<<"Total entries "<<n_htriton_Entries<< endl;


	for (Long64_t iEntry = 0; iEntry<=n_htriton_Entries; iEntry++)
	{
		htriton_tree.GetEntry(iEntry);
		//if(iEntry%100000==0){
		if(iEntry%10000==0){
		cout<<"Processing entry:" <<iEntry<<endl;
		}

	int pt_label = -1;

	for(int ipt=0;ipt<nptbins;ipt++){
		if(bpl > ptbin[ipt]) pt_label = ipt;
		}
	if(bpl > ptbin[nptbins]) continue;
        if(bpl<6) continue;//no signal for ht_pl<6

	//pre cuts
	if(chi2primary_he<chi2primary_he_cut || chi2primary_pi<chi2primary_pi_cut) continue;
        if(ht_chi2topo>ht_chi2topo_cut) continue;
        if(ht_chi2ndf>ht_chi2ndf_cut) continue;
	//fiductial cuts
	bpt = sqrt(bpx*bpx+bpy*bpy);
        TLorentzVector ptc(bpx,bpy,bpz,sqrt(bpx*bpx+bpy*bpy+bpz*bpz+bparticlemass*bparticlemass));
	brap = ptc.Rapidity() - ycm;
	if(bpt>3.0) continue;
        if(brap>0.95) continue;
        if( brap>0.55  &&  bpt < 3.9*(brap-0.45)*(brap-0.45)+0.65) continue;
        if( brap<0.55  &&  bpt < 1.9*(brap-0.65)*(brap-0.65)+0.67) continue;
	
	if(pt_label<0) continue;

	if(ht_chi2topo<0) continue;
        int cent_label = -1;
        cent_label = Centrality(countrefmult);          
	if(cent_label<0) continue;//<=5refmult
	if(bparticleid!=3004)  continue;

	if(_data_or_sim==0 && bismc==1){
		TLorentzVector mcptc(bmcpx,bmcpy,bmcpz,sqrt(bmcpx*bmcpx+bmcpy*bmcpy+bmcpz*bmcpz+ht_mass*ht_mass));
	        bmcrap = mcptc.Rapidity() - ycm;

		float mc_weight;
	        mc_weight = 1./g_pt_fine_in->GetBinContent(g_pt_fine_in->FindBin(bmcrap,sqrt(bmcpx*bmcpx+bmcpy*bmcpy)));	

		float rap_weight;
	        float pt_weight;
                float cent_weight;

                cent_weight = r_cent[cent_label];
		rap_weight = t_quadr0->Eval(bmcrap);
	        pt_weight =bolt1->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy));
		weight = rap_weight*pt_weight;
	        weight *= mc_weight;	
                weight *= cent_weight;
		htriton_tree_pt[pt_label]->Fill();
		}

	if(_data_or_sim==1){
		//if(bparticlemass > ht_mass-6*sgm3 && bparticlemass < ht_mass+6*sgm3){
		if( (bparticlemass > ht_mass+3*sgm3 && bparticlemass < ht_mass+6*sgm3) || (bparticlemass > ht_mass-6*sgm3 && bparticlemass < ht_mass-3*sgm3)){
		weight = 1;
		htriton_tree_pt[pt_label]->Fill();
			}
		}	

	}//iEntry

	write_tree();

}
void book_tree(){

	for(int ipt=0;ipt<nptbins;ipt++){

		htriton_tree_pt[ipt] = new TTree(Form("htriton_tree_pt[%d]",ipt), "");

		htriton_tree_pt[ipt]->Branch("ht_chi2ndf",&ht_chi2ndf,"ht_chi2ndf/F");
		htriton_tree_pt[ipt]->Branch("ht_chi2topo",&ht_chi2topo,"ht_chi2topo/F");

		htriton_tree_pt[ipt]->Branch("chi2primary_pi",&chi2primary_pi,"chi2primary_pi/F");
		htriton_tree_pt[ipt]->Branch("chi2primary_he",&chi2primary_he,"chi2primary_he/F");

		htriton_tree_pt[ipt]->Branch("ht_l",&ht_l,"ht_l/F");
		htriton_tree_pt[ipt]->Branch("ht_ldl",&ht_ldl,"ht_ldl/F");

		//htriton_tree_pt[ipt]->Branch("refMult",&refMult,"refMult/F");
		htriton_tree_pt[ipt]->Branch("bparticlemass",&bparticlemass,"bparticlemass/F");
		htriton_tree_pt[ipt]->Branch("bpl",&bpl,"bpl/F");

		htriton_tree_pt[ipt]->Branch("weight", &weight, "weight/D");

		}
}//book tree

void write_tree(){

//	TFile *outhistfile = new TFile (ofile, "RECREATE");
//	outhistfile->cd();

	for(int ipt=0;ipt<nptbins;ipt++){
		htriton_tree_pt[ipt]->Write();
	}
	outhistfile->Close();
	delete outhistfile;

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



