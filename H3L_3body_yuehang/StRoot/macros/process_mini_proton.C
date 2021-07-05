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
#include "TMVA/Reader.h"
#endif


const int MAXv0 = 10000;

const int npttmvabins = 1;
//const int nptbins = 6;
const int nptbins = 7;
const int ncenttmvabins = 2;
//const int ncentbins = 8;
const int ncentbins = 9;

TH2F *h_match_pt_pt;
TH2F *h_split_pt_pt;
TH2F *h_match_common;
TH2F *h_split_common;

 TH1F	*hvtx;
 TH1F	*hvtxgood;
 TH1F	*hrefmult;
 TH1F	*wrefmult;

 TH1F   *gvtx;
 TH1F   *gvtxgood;
 TH1F   *grefmult;
 TH1F   *frefmult;

 TH2F  *h_dedx_p;

 TH3F *h_match_nhitsfit;
 TH3F *h_match_nhitsdedx;
 TH3F *h_match_nhitsmax;

 TH3F *h_split_nhitsfit;
 TH3F *h_split_nhitsdedx;
 TH3F *h_split_nhitsmax;

TH1F  *h_hl_pl;
TH1F  *h_hl_pl_bg;

 TH1F  *h_ht_mass;
 TH1F  *h_ht_mass_bg;

 TH1F  *h_ht_ldl_bg;
 TH1F  *h_ht_dl_bg;
 TH1F  *h_ht_l_bg;

 TH2F *h_eta;
 TH2F *h_match_eta;
 TH2F *g_match_eta;

TH2F *h_match_eta_mc;
TH2F *h_split_eta_mc;

 TH2F *h_split_eta;
 TH2F *g_split_eta;
 TH2F *h_phi;
 TH2F *h_rap;
 TH2F *g_eta;

TH3F *h_nhits_he;
TH3F *h_nhits_pi;

TH2F *h_match_mc;;
TH2F *h_split_mc;

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
TH3F *h_ht_bdt;
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




//int snn=28;
int snn=3;

int    _cut_mode       ;



int    brunid          ;
int    beventid        ;
int    brefmult        ;
int    btofmult        ;
float  bVz             ;
float  bVzerr	;
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
int    _data_or_sim = 0;



/*
*/

/*
*/
Char_t path_ch[10000];
Char_t ofile[10000];




int mNMatchedPair;
int mNMcTrack;
int mNSplitPair;

float mMatchedPairsmPtMc[MAXv0];
float mMcTracksmPtMc[MAXv0];
float mMcTracksmEtaMc[MAXv0];
float mMcTracksmPhiMc[MAXv0];
UShort_t mMcTracksmGeantId[MAXv0];

float mMatchedPairsmPtGl[MAXv0];
float mMatchedPairsmEtaGl[MAXv0];
float mMatchedPairsmPhiGl[MAXv0];
Short_t mMatchedPairsmNHitMc[MAXv0];
Short_t mMatchedPairsmFitPts[MAXv0];
Short_t mMatchedPairsmDedxPts[MAXv0];
Short_t mMatchedPairsmNPossible[MAXv0];
float mSplitPairsmPhiMc[MAXv0];
Short_t mSplitPairsmNHitMc[MAXv0];
float mSplitPairsmPtGl[MAXv0];
float mSplitPairsmEtaGl[MAXv0];
float mSplitPairsmPhiGl[MAXv0];
Short_t mSplitPairsmFitPts[MAXv0];
Short_t mSplitPairsmDedxPts[MAXv0];
Short_t mSplitPairsmNPossible[MAXv0];
float mMatchedPairsmEtaMc[MAXv0];
float mSplitPairsmPtMc[MAXv0];
float mSplitPairsmEtaMc[MAXv0];

Short_t mMatchedPairsmNCommonHit[MAXv0];
Short_t mSplitPairsmNCommonHit[MAXv0];

TH2F *h_match_max;
TH2F *h_split_max;

int choose_cent;
int _batch;
int _wgt;
void book_histo();
void write_histo();

int _condor;
int _cuts;
void process_mini_proton(int condor=0){

	_condor = condor;

	TChain StMiniMcTree("StMiniMcTree");

	int nChains = -999999;
	ifstream fin;

	if(_condor==0){
	fin.open("/star/u/yhleung2/pwg/33.standard_reco/mini_x_proton.list");
	strcpy (ofile, "proton_x.root");
	nChains = 1031;
	}
		
	if(_condor==1){
        fin.open("/star/u/yhleung2/pwg/33.standard_reco/mini_proton.list");
        strcpy (ofile, "proton_1.root");
        nChains = 708;
        }


	for(int i=0; i<nChains; i++){
      	fin >> path_ch;

	StMiniMcTree.Add(path_ch);
	cout<<path_ch<<" "<<StMiniMcTree.GetEntries() <<endl;

	}//ichain

	cout << " done TChain-ing.. " << endl;

	book_histo();


	StMiniMcTree.SetBranchAddress("mNMcTrack",&mNMcTrack);
	StMiniMcTree.SetBranchAddress("mMcTracks.mPtMc",&mMcTracksmPtMc);
	StMiniMcTree.SetBranchAddress("mMcTracks.mEtaMc",&mMcTracksmEtaMc);
	StMiniMcTree.SetBranchAddress("mMcTracks.mPhiMc",&mMcTracksmPhiMc);
	StMiniMcTree.SetBranchAddress("mMcTracks.mGeantId",mMcTracksmGeantId);

	StMiniMcTree.SetBranchAddress("mNMatchedPair",&mNMatchedPair);

	StMiniMcTree.SetBranchAddress("mMatchedPairs.mPtMc",&mMatchedPairsmPtMc);
        StMiniMcTree.SetBranchAddress("mMatchedPairs.mEtaMc",&mMatchedPairsmEtaMc);

        StMiniMcTree.SetBranchAddress("mMatchedPairs.mPtGl",&mMatchedPairsmPtGl);
	StMiniMcTree.SetBranchAddress("mMatchedPairs.mEtaGl",&mMatchedPairsmEtaGl);
	StMiniMcTree.SetBranchAddress("mMatchedPairs.mPhiGl",&mMatchedPairsmPhiGl);
	StMiniMcTree.SetBranchAddress("mMatchedPairs.mNHitMc",&mMatchedPairsmNHitMc);
	StMiniMcTree.SetBranchAddress("mMatchedPairs.mFitPts",&mMatchedPairsmFitPts);
	StMiniMcTree.SetBranchAddress("mMatchedPairs.mDedxPts",&mMatchedPairsmDedxPts);
	StMiniMcTree.SetBranchAddress("mMatchedPairs.mNPossible",&mMatchedPairsmNPossible);

        StMiniMcTree.SetBranchAddress("mMatchedPairs.mNCommonHit",&mMatchedPairsmNCommonHit);

	StMiniMcTree.SetBranchAddress("mNSplitPair",&mNSplitPair);
	StMiniMcTree.SetBranchAddress("mSplitPairs.mPtMc",&mSplitPairsmPtMc);
	StMiniMcTree.SetBranchAddress("mSplitPairs.mEtaMc",&mSplitPairsmEtaMc);
	StMiniMcTree.SetBranchAddress("mSplitPairs.mPhiMc",&mSplitPairsmPhiMc);
	StMiniMcTree.SetBranchAddress("mSplitPairs.mNHitMc",&mSplitPairsmNHitMc);

	StMiniMcTree.SetBranchAddress("mSplitPairs.mNCommonHit",&mSplitPairsmNCommonHit);

	StMiniMcTree.SetBranchAddress("mSplitPairs.mPtGl",&mSplitPairsmPtGl);
	StMiniMcTree.SetBranchAddress("mSplitPairs.mEtaGl",&mSplitPairsmEtaGl);
	StMiniMcTree.SetBranchAddress("mSplitPairs.mPhiGl",&mSplitPairsmPhiGl);
	StMiniMcTree.SetBranchAddress("mSplitPairs.mFitPts",&mSplitPairsmFitPts);
	StMiniMcTree.SetBranchAddress("mSplitPairs.mDedxPts",&mSplitPairsmDedxPts);
	StMiniMcTree.SetBranchAddress("mSplitPairs.mNPossible",&mSplitPairsmNPossible);



/*
	StMiniMcTree.SetBranchAddress("brefmult",&brefmult);
	StMiniMcTree.SetBranchAddress("btofmult",&btofmult);
	StMiniMcTree.SetBranchAddress("bparticleid",&bparticleid);
	StMiniMcTree.SetBranchAddress("bparticlemass",&bparticlemass);
	StMiniMcTree.SetBranchAddress("bpx",&bpx);
	StMiniMcTree.SetBranchAddress("bpy",&bpy);
	StMiniMcTree.SetBranchAddress("bpz",&bpz);
	StMiniMcTree.SetBranchAddress("chi2primary_he", &chi2primary_he);
	StMiniMcTree.SetBranchAddress("chi2primary_pi", &chi2primary_pi);
	StMiniMcTree.SetBranchAddress("ht_chi2topo", &ht_chi2topo);
	StMiniMcTree.SetBranchAddress("ht_chi2ndf", &ht_chi2ndf);
	StMiniMcTree.SetBranchAddress("ht_ldl", &ht_ldl);
	StMiniMcTree.SetBranchAddress("ht_l", &ht_l);
	StMiniMcTree.SetBranchAddress("ht_dl", &ht_dl);
	StMiniMcTree.SetBranchAddress("dca_he",&dca_he);
	StMiniMcTree.SetBranchAddress("dca_pi",&dca_pi);
	StMiniMcTree.SetBranchAddress("nhits_pi",&nhits_pi);
	StMiniMcTree.SetBranchAddress("nhits_he",&nhits_he);
	StMiniMcTree.SetBranchAddress("ht_bdfvtx",&ht_bdfvtx);
	StMiniMcTree.SetBranchAddress("ht_bdfvtx2",&ht_bdfvtx2);
	StMiniMcTree.SetBranchAddress("ht_lifetime",&ht_lifetime);
	StMiniMcTree.SetBranchAddress("countrefmult",&countrefmult);
	StMiniMcTree.SetBranchAddress("reweight", &reweight);
	StMiniMcTree.SetBranchAddress("cent9", &cent9);
	StMiniMcTree.SetBranchAddress("bismc", &bismc);
*/

	Long64_t n_lambda_Entries = StMiniMcTree.GetEntries();

	cout<<endl<<"Start processing "<<endl<<endl;
	cout<<"--------------------------------------------"<<endl;

	cout<<"Total entries "<<n_lambda_Entries<< endl;

	for (Long64_t iEntry = 0; iEntry<n_lambda_Entries; iEntry++){

	StMiniMcTree.GetEntry(iEntry);
        if(iEntry%100000==0){
        cout<<"Processing entry:" <<iEntry<<endl;
        }

	//cout<<"mNMatchedPair:"<<mNMatchedPair<<endl;

	for(int i=0;i<mNMcTrack;i++){
		if(mMcTracksmGeantId[i]!=14) continue;
		//if(mMcTracksmGeantId[i]!=49) continue;

		h_eta->Fill(mMcTracksmEtaMc[i],mMcTracksmPtMc[i]);
		h_phi->Fill(mMcTracksmPhiMc[i],mMcTracksmPtMc[i]);
	}

	for(int i=0;i<mNMatchedPair;i++){
	//cout<<"mMatchedPairsmPtMc:"<<i<<":"<<mMatchedPairsmPtMc[i]<<endl;

	
		h_match_eta->Fill(mMatchedPairsmEtaGl[i],mMatchedPairsmPtGl[i]);
		h_match_eta_mc->Fill(mMatchedPairsmEtaMc[i],mMatchedPairsmPtMc[i]);

		h_match_nhitsfit->Fill(mMatchedPairsmEtaGl[i],mMatchedPairsmPtGl[i],mMatchedPairsmFitPts[i]);
		h_match_nhitsdedx->Fill(mMatchedPairsmEtaGl[i],mMatchedPairsmPtGl[i],mMatchedPairsmDedxPts[i]);
		h_match_nhitsmax->Fill(mMatchedPairsmEtaGl[i],mMatchedPairsmPtGl[i],mMatchedPairsmNPossible[i]);	

		if(mMatchedPairsmFitPts[i]>=15 && mMatchedPairsmDedxPts[i]>=5 && float(mMatchedPairsmFitPts[i])/float(mMatchedPairsmNPossible[i])>0.52){
		g_match_eta->Fill(mMatchedPairsmEtaGl[i],mMatchedPairsmPtGl[i]);
		}

//		cout<<"mMatchedPairsmNCommonHit[i]:"<<mMatchedPairsmNCommonHit[i]<<" "<<mMatchedPairsmFitPts[i]<<endl;

		h_match_pt_pt->Fill(mMatchedPairsmPtMc[i],mMatchedPairsmPtGl[i]);
		h_match_common->Fill(mMatchedPairsmNCommonHit[i], mMatchedPairsmFitPts[i]);
		h_match_mc->Fill(mMatchedPairsmNHitMc[i],mMatchedPairsmFitPts[i]);
                h_match_max->Fill(mMatchedPairsmNHitMc[i],mMatchedPairsmNCommonHit[i]);
	}


        for(int i=0;i<mNSplitPair;i++){
		h_split_eta->Fill(mSplitPairsmEtaGl[i],mSplitPairsmPtGl[i]);
		h_split_eta_mc->Fill(mSplitPairsmEtaMc[i],mSplitPairsmPtMc[i]);

                h_split_nhitsfit->Fill(mSplitPairsmEtaGl[i],mSplitPairsmPtGl[i],mSplitPairsmFitPts[i]);
                h_split_nhitsdedx->Fill(mSplitPairsmEtaGl[i],mSplitPairsmPtGl[i],mSplitPairsmDedxPts[i]);
                h_split_nhitsmax->Fill(mSplitPairsmEtaGl[i],mSplitPairsmPtGl[i],mSplitPairsmNPossible[i]);

		if(mSplitPairsmFitPts[i]>=15 && mSplitPairsmDedxPts[i]>=5 && float(mSplitPairsmFitPts[i])/float(mSplitPairsmNPossible[i])>0.52){
		g_split_eta->Fill(mSplitPairsmEtaGl[i],mSplitPairsmPtGl[i]);
		}

		h_split_pt_pt->Fill(mSplitPairsmPtMc[i],mSplitPairsmPtGl[i]);
                h_split_common->Fill(mSplitPairsmNCommonHit[i], mSplitPairsmFitPts[i]);
		h_split_mc->Fill(mSplitPairsmNHitMc[i],mSplitPairsmFitPts[i]);
                h_split_max->Fill(mSplitPairsmNHitMc[i],mSplitPairsmNCommonHit[i]);
        }






   	//bpt = sqrt(bpx*bpx+bpy*bpy);
  	//bp = sqrt(bpx*bpx+bpy*bpy+bpz*bpz);


	//TLorentzVector ptc(bpx,bpy,bpz,sqrt(bpx*bpx+bpy*bpy+bpz*bpz+bparticlemass*bparticlemass));
  	//beta = ptc.Eta(); 
  	//brap = ptc.Rapidity() - ycm;

	//h_pt->Fill(brap,bpt);


	}//entry loop


	write_histo();

}

void book_histo()
{

  h_dedx_p = new TH2F("h_dedx_p","",500,0,5,2000,0,100);
  h_dedx_p->Sumw2();

  h_mass_pt = new TH3F("h_mass_pt","",1000,2.5,3.5,100,-1.0,1.0,50,0,5);
  h_mass_pl = new TH3F("h_mass_pl","",1000,2.5,3.5,10,-1.0,1.0,500,0,100);

  h_mass_pt_wgt = new TH3F("h_mass_pt_wgt","",1000,2.5,3.5,100,-1.0,1.0,50,0,5);
  h_mass_pl_wgt = new TH3F("h_mass_pl_wgt","",1000,2.5,3.5,10,-1.0,1.0,500,0,100);

  h_mass_qa = new TH3F("h_mass_qa","",10,-2.0,0.0,50,0,5,1000,2.5,3.5);
  h_chi2primary_he = new TH3F("h_chi2primary_he","",10,-2.0,0.0,50,0,5,2000,0,2000);
  h_chi2primary_pi = new TH3F("h_chi2primary_pi","",10,-2.0,0.0,50,0,5,20000,0,20000);
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

  h_eta = new TH2F("h_eta","",100,-2.,0.,50,0,5);
  g_eta = new TH2F("g_eta","",100,-2.,0.,50,0,5);
  h_match_eta = new TH2F("h_match_eta","",100,-2.,0.,50,0,5);
  g_match_eta = new TH2F("g_match_eta","",100,-2.,0.,50,0,5);
  h_split_eta = new TH2F("h_split_eta","",100,-2.,0.,50,0,5);
  g_split_eta = new TH2F("g_split_eta","",100,-2.,0.,50,0,5);

  h_match_eta_mc = new TH2F("h_match_eta_mc","",100,-2.,0.,50,0,5);
  h_split_eta_mc = new TH2F("h_split_eta_mc","",100,-2.,0.,50,0,5);

  h_phi = new TH2F("h_phi","",200,-6.5,6.5,100,0,10);
  h_rap = new TH2F("h_rap","",200,-1.5,1.5,100,0,10);
  h_dca = new TH3F("h_dca","",10,-2.0,0.0,50,0,5,500,0,20);
  g_pt = new TH2F("g_pt","",100,-1.5,1.5,50,0,5);
  g_pl = new TH2F("g_pl","",10,-1.0,1.0,500,0,100);
  g_pt_fine = new TH2F("g_pt_fine","",100,-1.0,1.0,50,0,5);
  g_pt_fine2 = new TH2F("g_pt_fine2","",20,-1.0,1.0,50,0,5);
  g_pt_fine_wgt = new TH2F("g_pt_fine_wgt","",100,-1.0,1.0,50,0,5);
  g_pt_fine2_wgt = new TH2F("g_pt_fine2_wgt","",20,-1.0,1.0,50,0,5);

  g_pt_wgt = new TH2F("g_pt_wgt","",100,-1.5,1.5,50,0,5);
  g_pl_wgt = new TH2F("g_pl_wgt","",10,-1.0,1.0,500,0,100);


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

  h_match_nhitsfit = new TH3F("h_match_nhitsfit","",10,-2.0,0.0,50,0,5,60,0,60);
  h_match_nhitsdedx = new TH3F("h_match_nhitsdedx","",10,-2.0,0.0,50,0,5,60,0,60);
  h_match_nhitsmax = new TH3F("h_match_nhitsmax","",10,-2.0,0.0,50,0,5,60,0,60);

  h_match_pt_pt = new TH2F("h_match_pt_pt","",100,0,5,100,0,5);
  h_split_pt_pt = new TH2F("h_split_pt_pt","",100,0,5,100,0,5);

  h_match_common = new TH2F("h_match_common","",50,0,50,50,0,50);
  h_split_common = new TH2F("h_split_common","",50,0,50,50,0,50);

h_match_mc = new TH2F("h_match_mc","",50,0,50,50,0,50);
h_split_mc = new TH2F("h_split_mc","",50,0,50,50,0,50);

  h_match_pt_pt->Sumw2();
  h_split_pt_pt->Sumw2();
  h_match_common->Sumw2();
  h_split_common->Sumw2();

h_match_max= new TH2F("h_match_max","",50,0,50,50,0,50);
h_split_max= new TH2F("h_split_max","",50,0,50,50,0,50);

h_match_max->Sumw2();
h_split_max->Sumw2();

  h_match_nhitsfit->Sumw2();
  h_match_nhitsdedx->Sumw2();
  h_match_nhitsmax->Sumw2();

  h_split_nhitsfit = new TH3F("h_split_nhitsfit","",10,-2.0,0.0,50,0,5,60,0,60);
  h_split_nhitsdedx = new TH3F("h_split_nhitsdedx","",10,-2.0,0.0,50,0,5,60,0,60);
  h_split_nhitsmax = new TH3F("h_split_nhitsmax","",10,-2.0,0.0,50,0,5,60,0,60);

  h_split_nhitsfit->Sumw2();
  h_split_nhitsdedx->Sumw2();
  h_split_nhitsmax->Sumw2();

  g_pt_wgt->Sumw2();
  g_pl_wgt->Sumw2();
  h_mass_pt_wgt->Sumw2();
  h_mass_pl_wgt->Sumw2();
  h_mass_pt->Sumw2();
  h_mass_pl->Sumw2();
h_match_mc->Sumw2();
h_split_mc->Sumw2();

  h_ht_bdfvtx->Sumw2();
  h_ht_bdfvtx2->Sumw2();
  h_ht_lifetime->Sumw2();
  h_lt_mc->Sumw2();
  g_lt_mc->Sumw2();

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
  h_ht_bdt->Sumw2();
  h_nhits_he ->Sumw2();
  h_nhits_pi ->Sumw2();
  h_eta->Sumw2();
  h_match_eta->Sumw2();
  g_match_eta->Sumw2();
  h_split_eta->Sumw2();
  g_split_eta->Sumw2();

  h_match_eta_mc->Sumw2();
  h_split_eta_mc->Sumw2();

  h_phi->Sumw2();
  h_rap->Sumw2();
  h_cent->Sumw2();

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
h_match_mc->Write();
h_split_mc->Write();
h_eta->Write();
h_match_eta->Write();
g_match_eta->Write();
h_split_eta->Write();
g_split_eta->Write();

h_match_eta_mc->Write();
h_split_eta_mc->Write();

h_match_max->Write();
h_split_max->Write();
h_phi->Write();

h_match_nhitsfit->Write();
h_match_nhitsdedx->Write();
h_match_nhitsmax->Write();

h_split_nhitsfit->Write();
h_split_nhitsdedx->Write();
h_split_nhitsmax->Write();
h_match_pt_pt->Write();
h_split_pt_pt->Write();
h_match_common->Write();
h_split_common->Write();
}

