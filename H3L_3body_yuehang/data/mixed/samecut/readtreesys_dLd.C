#include "tree.h"
#include "Hists.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TStopwatch.h"
#include "TString.h"

#include <fstream>

// void readtree(TString mInputlist="H3L3b_tree_mc.root", int mode = 0, TString outfile="fout_H3L.root")
void readtreesys_dLd(TString mInputlist="Lambda_tree_mc.root", int const mode = 1,   TString outfile="fout_Lambda.root", int const mcState=1, int const isMix=0, int centLow=3, int centHigh=8, bool applycorr=0)
{
  bool fillQAplots=1;
  TStopwatch time;
  time.Start();
  double snn = 3;
  double ycm;
  if(snn==3) //target
  {
    ycm = -1.045;
  }else{
    ycm = -999.;
  }
  double const mass_ht = 2.99131;
  double const ht_width = 0.005;
  double const mass_hl = 3.9239;
  double const hl_width = 0.005;
  double const mass_ld = 1.11568;
  double const mass_d = 1.8756;
  double const mass_p = 0.93827;
  double const mass_pi = 0.13957;
  TGraph* gdLdCorr = new TGraph("corr.txt");

TF1* fH3Ldydpt[3];
double par[3][3]={{0.267422, 43782.6, 2.99339}, {0.200836, 92736.2,2.99339}, {0.136789, 136882, 2.99339} };
for (int irap=0;irap<3;irap++) {
      fH3Ldydpt[irap] = new TF1(Form("fH3Ldydpt[%d]",irap), "x*[1]*exp(-(sqrt([2]*[2]+x*x)-[2])/[0])", 0,5);
      fH3Ldydpt[irap]->SetParameters(par[irap]);
}

  //weighting function for MC
  if (mcState!=0){
    cout <<"starting book mc functions!" << endl;
    levyfit4 = new TF1("levyfit4","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.4,2.4);
    levyfit5 = new TF1("levyfit5","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.2,2.4);
    levyfit6 = new TF1("levyfit6","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.2,2.4);
    levyfit7 = new TF1("levyfit7","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.2,2.4);
    levyfit8 = new TF1("levyfit8","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.2,2.4);
    t_quadr = new TF1("t_quadr","[0]+[1]*x+[2]*x*x",-2,2);
    t_quadr0 = new TF1("t_quadr0","[0]+[1]*x+[2]*x*x",-2,2);
    t_quadr1 = new TF1("t_quadr1","[0]+[1]*x+[2]*x*x",-2,2);
    bolt0 = new TF1("bolt0", "1e9*x*exp(-(sqrt(2.990*2.990+x*x)-2.990)/[0])", 0.,5);
    bolt1 = new TF1("bolt1", "1e9*x*exp(-(sqrt(2.990*2.990+x*x)-2.990)/[0])", 0.,5);
    bolt2 = new TF1("bolt2", "1e9*x*exp(-(sqrt(2.990*2.990+x*x)-2.990)/[0])", 0.,5);


    if (mode==0 && mcState==1) {
      levyfit4->SetParameters(1.53536e-01 , -1.21064e+08, 5.25287e+07);
      levyfit5->SetParameters(1.51332e-01 , -1.25927e+08, 4.86193e+07);
      levyfit6->SetParameters(1.45903e-01 , -1.60954e+08, 4.74115e+07);
      levyfit7->SetParameters(1.28382e-01 , -1.60610e+08, 4.92025e+07);
      levyfit8->SetParameters(1.11371e-01 , -2.04689e+08, 4.22255e+07);

      t_quadr->SetParameters(1.25929e+00,0,-1.80963e+00-0.5);
      t_quadr0->SetParameters(1.25929e+00,0,-1.80963e+00);
      t_quadr1->SetParameters(1.25929e+00,0,-1.80963e+00+0.5);

      bolt0->SetParameter(0,0.34);
      bolt1->SetParameter(0,0.27);
      bolt2->SetParameter(0,0.20);
    }

    if (mode==1 || (mode ==0 && mcState==-20)) //MC Lambda or (MC Lambda)+d
    {
      levyfit4->SetParameters(0.141849, -1.21226e+08, 1.879e+07);
      levyfit5->SetParameters(0.145031, -1.19088e+08, 1.43524e+07);
      levyfit6->SetParameters(0.141403, -2.20994e+08, 1.2987e+07);
      levyfit7->SetParameters(0.128239, -1.60612e+08, 1.17238e+07);
      levyfit8->SetParameters(0.114106, -1.35632e+08, 8.78861e+06);
      t_quadr->SetParameters(1.15426e+00, 0 ,-8.55891e-01 ); 
    }

    ////////uniform  mc pt y distribution/////////
 if (mode==0 && mcState==1)
   {
       fgpt_0= new TFile("fMC_H3L_wt.root","READ");
       g_pt_fine_in = (TH2F*)fgpt_0->Get("hPhase")->Clone("g_pt_fine_in");
     }
 if (mode==0 && mcState==1 && mInputlist.Contains("phase")) {
     fgpt_0= new TFile("fH3L_phase_phase_wt.root","READ");  //phase space decay
     g_pt_fine_in = (TH2F*)fgpt_0->Get("hPhase")->Clone("g_pt_fine_in");
   }
 // if (mode==1) fgpt_0 = new TFile("ld_input_mc_fine_v20.root","READ");//more stats
   // if (mode==1 || (mcState==-20 && mode==0)) fgpt_0 = new TFile("flambda_phase_wt.root","READ");//more stats
   if (mode==1 || (mcState==-20 && mode==0)) {
       fgpt_0 = new TFile("ld_input_mc_fine_v20.root","READ");//more stats
       g_pt_fine_in = (TH2F*)fgpt_0->Get("g_pt_fine")->Clone("g_pt_fine_in");
   }

    cout <<"finish book mc weighting functions." << endl;
  }
  ///////////////////////end of MC weighting////////////////

  cout <<"start read tree!" << endl;
  TString treename;
  if (mode==1) treename = "lambda_tree";
  if (mode==0) treename = "htriton3_tree";
  TChain htriton3_tree(treename.Data()); 

  // TH1F* hrefmult  = new TH1F("hrefmult_tot", "refmult; hrefmult; N_{evt}", 600,0,600);
  TH1F* hcent  = new TH1F("hcent", "cent; centrality; N_{evt}", 9,-0.5,8.5);
  TH1F* hcentwt  = new TH1F("hcentwt", "centwt; centrality; N_{evt}", 9,-0.5,8.5);
  // hrefmult->SetDirectory(0);

  if (mInputlist.Contains(".root"))
  {
    htriton3_tree.Add(mInputlist.Data());
    TFile *ftmp = new TFile(mInputlist.Data());
    // TH1F* htmp = (TH1F*)ftmp->Get("hrefmult");
    // hrefmult->Add(htmp);
    // ftmp->Close();
    TH1F* htmp = (TH1F*)ftmp->Get("hCent");
    TH1F* htmp2 = (TH1F*)ftmp->Get("hCentWt");
    hcent->Add(htmp);
    hcentwt->Add(htmp2);

  }
  else
  {
    int nfile = 0;
    char tmp[2000];
    ifstream readlists;
    readlists.open(mInputlist.Data());
    while (readlists.good()){
      readlists.getline(tmp,2000);
      TFile *ftmp = new TFile(tmp);
      if (!ftmp||!(ftmp->IsOpen())||!(ftmp->GetNkeys())) {
        cout<<"Could not open this file: "<< tmp  <<endl;
        continue;
      }
      else {
        if(nfile%1==0) cout<<"read in "<<nfile<<"th file: "<< tmp <<endl;
        htriton3_tree.Add(tmp);
        TH1F* htmp = (TH1F*)ftmp->Get("hCent");
        TH1F* htmp2 = (TH1F*)ftmp->Get("hCentWt");
        hcent->Add(htmp);
        hcentwt->Add(htmp2);
        nfile++;
        ftmp->Close();
      }
    }
  }
  cout <<"finish read tree" << endl;

  // double NTotEvents = hrefmult->Integral(); //if have weight, use hwref or hCentWt

  htriton3_tree.SetBranchAddress("bismc", &bismc);
  htriton3_tree.SetBranchAddress("bparticlemass",&bparticlemass);
  htriton3_tree.SetBranchAddress("bparticleid",&bparticleid);
  htriton3_tree.SetBranchAddress("bpx",&bpx);
  htriton3_tree.SetBranchAddress("bpy",&bpy);
  htriton3_tree.SetBranchAddress("bpz",&bpz);

  htriton3_tree.SetBranchAddress("dca_proton",&dca_proton);
  htriton3_tree.SetBranchAddress("chi2primary_proton", &chi2primary_proton);
  htriton3_tree.SetBranchAddress("chi2primary_pi", &chi2primary_pi);
  htriton3_tree.SetBranchAddress("reweight", &reweight);
  htriton3_tree.SetBranchAddress("cent9", &cent9);

  if (mode==0){
    htriton3_tree.SetBranchAddress("chi2primary_d", &chi2primary_d);
    htriton3_tree.SetBranchAddress("bzdeuteron", &bzdeuteron);
    htriton3_tree.SetBranchAddress("bpionnsigma", & bpionnsigma);
    htriton3_tree.SetBranchAddress("bprotonsigma", & bprotonsigma);

    htriton3_tree.SetBranchAddress("ht_chi2topo", &ht_chi2topo);
    htriton3_tree.SetBranchAddress("ht_chi2ndf", &ht_chi2ndf);
    htriton3_tree.SetBranchAddress("ht_NDF", &ht_NDF);
    htriton3_tree.SetBranchAddress("ht_ldl", &ht_ldl);
    htriton3_tree.SetBranchAddress("ht_l", &ht_l);
    htriton3_tree.SetBranchAddress("ht_dl", &ht_dl);
    htriton3_tree.SetBranchAddress("dca_proton",&dca_proton);
    htriton3_tree.SetBranchAddress("dca_pion",&dca_pion);
    htriton3_tree.SetBranchAddress("dca_deuteron",&dca_deuteron);
    htriton3_tree.SetBranchAddress("nhits_pion",&nhits_pion);
    htriton3_tree.SetBranchAddress("nhits_deuteron",&nhits_deuteron);
    htriton3_tree.SetBranchAddress("nhits_proton",&nhits_proton);
    //htriton3_tree.SetBranchAddress("ht_bdfvtx",&ht_bdfvtx);
    //htriton3_tree.SetBranchAddress("ht_bdfvtx2",&ht_bdfvtx2);
    //htriton3_tree.SetBranchAddress("ht_lifetime",&ht_lifetime);
    htriton3_tree.SetBranchAddress("countrefmult",&countrefmult);
    htriton3_tree.SetBranchAddress("reweight", &reweight);
    htriton3_tree.SetBranchAddress("cent9", &cent9);
    htriton3_tree.SetBranchAddress("bismc", &bismc);

    htriton3_tree.SetBranchAddress("bdpx", &bdpx);
    htriton3_tree.SetBranchAddress("bdpy", &bdpy);
    htriton3_tree.SetBranchAddress("bdpz", &bdpz);
    htriton3_tree.SetBranchAddress("bpionpx", &bpionpx);
    htriton3_tree.SetBranchAddress("bpionpy", &bpionpy);
    htriton3_tree.SetBranchAddress("bpionpz", &bpionpz);
    htriton3_tree.SetBranchAddress("bprotonpx", &bprotonpx);
    htriton3_tree.SetBranchAddress("bprotonpy", &bprotonpy);
    htriton3_tree.SetBranchAddress("bprotonpz", &bprotonpz);

    htriton3_tree.SetBranchAddress("v_01_pvdca", &v_01_pvdca); //pair distance from PV 
    htriton3_tree.SetBranchAddress("v_12_dca", &v_12_dca); //pair 
    htriton3_tree.SetBranchAddress("v_01_chi2primary", &v_01_chi2primary); //
    htriton3_tree.SetBranchAddress("v_01_chi2ndf", &v_01_chi2ndf); //
    htriton3_tree.SetBranchAddress("mass_01", &mass_01); // same as below
    /* htriton3_tree.SetBranchAddress("v_lambda_mass_0", &v_lambda_mass_0); */
    htriton3_tree.SetBranchAddress("v_lambda_ldl_0", &v_lambda_ldl_0);  // will add later
    htriton3_tree.SetBranchAddress("v_lambda_l_0", &v_lambda_l_0);  // will add later
    htriton3_tree.SetBranchAddress("v_012_dca", &v_012_dca);  // will add later
    htriton3_tree.SetBranchAddress("v_12_x", &v_12_x);  // will add later
    htriton3_tree.SetBranchAddress("v_12_y", &v_12_y);  // will add later
    htriton3_tree.SetBranchAddress("v_12_z", &v_12_z);  // will add later

    htriton3_tree.SetBranchAddress("v_02_x", &v_02_x);  // will add later
    htriton3_tree.SetBranchAddress("v_02_y", &v_02_y);  // will add later
    htriton3_tree.SetBranchAddress("v_02_z", &v_02_z);  // will add later


    /* htriton3_tree.SetBranchAddress("bdpx", &bdpx); */
    /* htriton3_tree.SetBranchAddress("bdpy", &bdpy); */
    /* htriton3_tree.SetBranchAddress("bdpz", &bdpz); */
    htriton3_tree.SetBranchAddress("bpionpx", &bpionpx);
    htriton3_tree.SetBranchAddress("bpionpy", &bpionpy);
    htriton3_tree.SetBranchAddress("bpionpz", &bpionpz);
    htriton3_tree.SetBranchAddress("bprotonpx", &bprotonpx);
    htriton3_tree.SetBranchAddress("bprotonpy", &bprotonpy);
    htriton3_tree.SetBranchAddress("bprotonpz", &bprotonpz);
    htriton3_tree.SetBranchAddress("b0mcpx",&b0mcpx);
    htriton3_tree.SetBranchAddress("b0mcpy",&b0mcpy);
    htriton3_tree.SetBranchAddress("b0mcpz",&b0mcpz);
    htriton3_tree.SetBranchAddress("b1mcpx",&b1mcpx);
    htriton3_tree.SetBranchAddress("b1mcpy",&b1mcpy);
    htriton3_tree.SetBranchAddress("b1mcpz",&b1mcpz);
    if (isMix>-1) htriton3_tree.SetBranchAddress("bisMix",&bisMix);
    else bisMix=isMix;
  }
  htriton3_tree.SetBranchAddress("bmcpx",&bmcpx);
  htriton3_tree.SetBranchAddress("bmcpy",&bmcpy);                                               
  htriton3_tree.SetBranchAddress("bmcpz",&bmcpz);
  htriton3_tree.SetBranchAddress("bmcl",&bmcl);
  htriton3_tree.SetBranchAddress("bmcpl",&bmcpl);

  if (mode ==1)
  {
    htriton3_tree.SetBranchAddress("ld_chi2ndf",&v_01_chi2ndf);
    htriton3_tree.SetBranchAddress("ld_chi2primary",&v_01_chi2primary);
    htriton3_tree.SetBranchAddress("ld_ldl",&v_lambda_ldl_0);
    htriton3_tree.SetBranchAddress("ld_l",&v_lambda_l_0);
    htriton3_tree.SetBranchAddress("dca_pi",&dca_pion);
  }

  TH1F* hcharge= new TH1F("hcharge","hcharge",2,-1.5,1.5);
  ////Background QA histrograms////////
  //(ppi) QA
  // TH2F* hptppimass = new TH2F("hptppimass","hptppimass;p_{T};mass",100,0,10,400,1.05,1.15);
  // hptppimass->Sumw2();
  // TH2F* hptppichi2prim= new TH2F("hptppichi2prim","hptppichi2prim;p_{T};(p#pi) #chi^{2}_{prim}",100,0,10,50,0,40);
  // hptppichi2prim->Sumw2();
  // TH2F* hptppil= new TH2F("hptppil","hptppil;p_{T};(p#pi) l",100,0,10,50,0,40);
  // hptppil->Sumw2();
  // TH2F* hptppildl= new TH2F("hptppildl","hptppildl;p_{T};(p#pi) ldl",100,0,10,50,0,25);
  // hptppildl->Sumw2();
  // TH2F* hptppichi2ndf= new TH2F("hptppichi2ndf","hptppichi2ndf;p_{T};(p#pi) #chi^{2}_{ndf}",100,0,5,80,0,0.5);
  // hptppichi2ndf->Sumw2();
  // TH2F* hptpichi2prim= new TH2F("hptpichi2prim","hptpichi2prim;p_{T};#pi #chi^{2}_{prim}",100,0,10,50,0,40);
  // hptpichi2prim->Sumw2();
  // TH2F* hptpchi2prim= new TH2F("hptpchi2prim","hptpchi2prim;p_{T};p #chi^{2}_{prim}",100,0,10,50,0,40);
  // hptpchi2prim->Sumw2();
  // TH2F* hptpidca = new TH2F("hptpidca","hptpidca;p_{T};DCA",100,0,10,100,0,15);
  // hptpidca->Sumw2();
  // TH2F* hptpdca = new TH2F("hptpdca","hptpdca;p_{T};DCA",100,0,10,100,0,10);
  // hptpdca->Sumw2();
  // TH2F* hptsumdca = new TH2F("hptsumdca","hptsumdca;p_{T};DCA",100,0,10,100,0,20);
  // hptsumdca->Sumw2();

  // TH2F* hptppimassSB = new TH2F("hptppimassSB","hptppimass;p_{T};mass",100,0,10,400,1.05,1.15);
  // hptppimassSB->Sumw2();
  // TH2F* hptppichi2primSB= new TH2F("hptppichi2primSB","hptppichi2prim;p_{T};(p#pi) #chi^{2}_{prim}",100,0,10,50,0,40);
  // hptppichi2primSB->Sumw2();
  // TH2F* hptppilSB= new TH2F("hptppilSB","hptppil;p_{T};(p#pi) l",100,0,10,50,0,40);
  // hptppilSB->Sumw2();
  // TH2F* hptppildlSB= new TH2F("hptppildlSB","hptppildl;p_{T};(p#pi) ldl",100,0,10,50,0,25);
  // hptppildlSB->Sumw2();
  // TH2F* hptppichi2ndfSB= new TH2F("hptppichi2ndfSB","hptppichi2ndf;p_{T};(p#pi) #chi^{2}_{ndf}",100,0,5,80,0,0.5);
  // hptppichi2ndfSB->Sumw2();

  // TH2F* hptppimassSig = new TH2F("hptppimassSig","hptppimass;p_{T};mass",100,0,10,400,1.05,1.15);
  // hptppimassSig->Sumw2();
  // TH2F* hptppichi2primSig= new TH2F("hptppichi2primSig","hptppichi2prim;p_{T};(p#pi) #chi^{2}_{prim}",100,0,10,50,0,40);
  // hptppichi2primSig->Sumw2();
  // TH2F* hptppilSig= new TH2F("hptppilSig","hptppil;p_{T};(p#pi) l",100,0,10,50,0,40);
  // hptppilSig->Sumw2();
  // TH2F* hptppildlSig= new TH2F("hptppildlSig","hptppildl;p_{T};(p#pi) ldl",100,0,10,50,0,25);
  // hptppildlSig->Sumw2();
  // TH2F* hptppichi2ndfSig= new TH2F("hptppichi2ndfSig","hptppichi2ndf;p_{T};(p#pi) #chi^{2}_{ndf}",100,0,5,80,0,0.5);
  // hptppichi2ndfSig->Sumw2();

  /* TH2F* hPhase = new TH2F("hPhase","hPhase;y;pt",100,-1,1,250,0,5); */

  TH2F* hptH3Lmass  = new TH2F("hptH3Lmass","hptH3Lmass;p_{T};H3L mass",100,0,5,200,2.95,3.05);
  hptH3Lmass->Sumw2();
  TH3F* hH3LMassPtY= new TH3F("hH3LMassPtY","hH3LMassPtY;p_{T};H3L mass;Rapidity",100,0,5,200,2.95,3.05, 100, -1.,0);
  hH3LMassPtY->Sumw2();
  TH3F* hH3LKstarPtY= new TH3F("hH3LMassPtY","hH3LMassPtY;p_{T};H3L kstar;Rapidity",100,0,5,300, 0, 150, 100, -1.,0);
  hH3LKstarPtY->Sumw2();


  TH3F* hH3LMassPtCent= new TH3F("hH3LMassPtCent","hH3LMassPtCent;p_{T};H3L mass;Centrality",100,0,5,200,2.95,3.05, 9, -0.5,8.5);
  hH3LMassPtCent->Sumw2();
  TH3F* hH3LMassPtY_5_40= new TH3F("hH3LMassPtY_5_40","hH3LMassPtY;p_{T};H3L mass;Rapidity",100,0,5,200,2.95,3.05, 60, -1.5,0);
  hH3LMassPtY_5_40->Sumw2();

  TH2F* hH3LptProtonPt = new TH2F("hH3LptProtonPt","hH3LptProtonPt; p_{T};", 100, 0, 5, 80, 0, 4);
  hH3LptProtonPt->Sumw2();
  TH2F* hH3LptPionPt = new TH2F("hH3LptPionPt","hH3LptPionPt; p_{T};", 100, 0, 5, 60, 0, 3);
  hH3LptPionPt->Sumw2();

  //topological variable for H3L
  TH2F* hptH3L_l= new TH2F("hptH3L_l","hptH3L_l;p_{T};l",100,0,5,500,0,100);
  hptH3L_l->Sumw2();
  TH2F* hptH3L_ldl= new TH2F("hptH3L_ldl","hptH3L_ldl;p_{T};ldl",100,0,5,250,0,50);
  hptH3L_ldl->Sumw2();
  TH2F* hptH3L_dchi2prim= new TH2F("hptH3L_dchi2prim","hptH3L_dchi2pri;p_{T};d chi2prim",100,0,5,100,0,20);
  hptH3L_dchi2prim->Sumw2();
  TH2F* hptH3L_pichi2prim= new TH2F("hptH3L_pichi2prim","hptH3L_pichi2pri;p_{T};pi chi2prim",100,0,5,100,0,20);
  hptH3L_pichi2prim->Sumw2();
  TH2F* hptH3L_pchi2prim= new TH2F("hptH3L_pchi2prim","hptH3L_pchi2pri;p_{T};p chi2prim",100,0,5,100,0,20);
  hptH3L_pchi2prim->Sumw2();
  TH2F* hptH3L_chi2topo= new TH2F("hptH3L_chi2topo","hptH3L_chi2topo;p_{T};chi2topo",100,0,5,100,0,10);
  hptH3L_chi2topo->Sumw2();
  TH2F* hptH3L_chi2ndf= new TH2F("hptH3L_chi2ndf","hptH3L_chi2topo;p_{T};chi2ndf",100,0,5,100,0,10);
  hptH3L_chi2ndf->Sumw2();
  TH2F* hptH3L_dDca= new TH2F("hptH3L_dDca","hptH3L_dDca;p_{T};dDca",100,0,5,100,0,5);
  hptH3L_dDca->Sumw2();
  TH2F* hptH3L_piDca= new TH2F("hptH3L_piDca","hptH3L_piDca;p_{T};piDca",100,0,5,100,0,20);
  hptH3L_piDca->Sumw2();
  TH2F* hptH3L_pDca= new TH2F("hptH3L_pDca","hptH3L_pDca;p_{T};pDca",100,0,5,100,0,10);
  hptH3L_pDca->Sumw2();
  TH2F* hptH3L_dpDca= new TH2F("hptH3L_dpDca","hptH3L_dpDca;p_{T};dpDca",100,0,5,80,0,4);
  hptH3L_dpDca->Sumw2();

  //signal region

  TH2F* hptH3L_lSig= new TH2F("hptH3L_lSig","hptH3L_l;p_{T};l",100,0,5,500,0,100);
  hptH3L_lSig->Sumw2();
  TH2F* hptH3L_ldlSig = new TH2F("hptH3L_ldlSig","hptH3L_ldl;p_{T};ldl",100,0,5,250,0,50);
  hptH3L_ldlSig->Sumw2();
  TH2F* hptH3L_dchi2primSig = new TH2F("hptH3L_dchi2primSig","hptH3L_dchi2pri;p_{T};d chi2prim",100,0,5,100,0,20);
  hptH3L_dchi2primSig->Sumw2();
  TH2F* hptH3L_pichi2primSig = new TH2F("hptH3L_pichi2primSig","hptH3L_pichi2pri;p_{T};pi chi2prim",100,0,5,100,0,20);
  hptH3L_pichi2primSig->Sumw2();
  TH2F* hptH3L_pchi2primSig= new TH2F("hptH3L_pchi2primSig","hptH3L_pchi2pri;p_{T};p chi2prim",100,0,5,100,0,20);
  hptH3L_pchi2primSig->Sumw2();
  TH2F* hptH3L_chi2topoSig= new TH2F("hptH3L_chi2topoSig","hptH3L_chi2topo;p_{T};chi2topo",100,0,5,100,0,10);
  hptH3L_chi2topoSig->Sumw2();
  TH2F* hptH3L_chi2ndfSig= new TH2F("hptH3L_chi2ndfSig","hptH3L_chi2ndf;p_{T};chi2ndf",100,0,5,100,0,10);
  hptH3L_chi2ndfSig->Sumw2();
  TH2F* hptH3L_dDcaSig= new TH2F("hptH3L_dDcaSig","hptH3L_dDca;p_{T};dDca",100,0,5,100,0,5);
  hptH3L_dDcaSig->Sumw2();
  TH2F* hptH3L_piDcaSig= new TH2F("hptH3L_piDcaSig","hptH3L_piDca;p_{T};piDca",100,0,5,100,0,20);
  hptH3L_piDcaSig->Sumw2();
  TH2F* hptH3L_pDcaSig = new TH2F("hptH3L_pDcaSig","hptH3L_pDca;p_{T};pDca",100,0,5,100,0,10);
  hptH3L_pDcaSig->Sumw2();
  TH2F* hptH3L_dpDcaSig = new TH2F("hptH3L_dpDcaSig","hptH3L_dpDca;p_{T};dpDca",100,0,5,80,0,4);
  hptH3L_dpDcaSig->Sumw2();
  //SB region
  TH2F* hptH3L_lSBL= new TH2F("hptH3L_lSBL","hptH3L_l;p_{T};l",100,0,5,500,0,100);
  hptH3L_lSBL->Sumw2();
  TH2F* hptH3L_ldlSBL = new TH2F("hptH3L_ldlSBL","hptH3L_ldl;p_{T};ldl",100,0,5,250,0,50);
  hptH3L_ldlSBL->Sumw2();
  TH2F* hptH3L_dchi2primSBL = new TH2F("hptH3L_dchi2primSBL","hptH3L_dchi2pri;p_{T};d chi2prim",100,0,5,100,0,20);
  hptH3L_dchi2primSBL->Sumw2();
  TH2F* hptH3L_pichi2primSBL = new TH2F("hptH3L_pichi2primSBL","hptH3L_pichi2pri;p_{T};pi chi2prim",100,0,5,100,0,20);
  hptH3L_pichi2primSBL->Sumw2();
  TH2F* hptH3L_pchi2primSBL= new TH2F("hptH3L_pchi2primSBL","hptH3L_pchi2pri;p_{T};p chi2prim",100,0,5,100,0,20);
  hptH3L_pchi2primSBL->Sumw2();
  TH2F* hptH3L_chi2topoSBL= new TH2F("hptH3L_chi2topoSBL","hptH3L_chi2topo;p_{T};chi2topo",100,0,5,100,0,10);
  hptH3L_chi2topoSBL->Sumw2();
  TH2F* hptH3L_chi2ndfSBL= new TH2F("hptH3L_chi2ndfSBL","hptH3L_chi2ndfSBL;p_{T};chi2ndf",100,0,5,100,0,10);
  hptH3L_chi2ndfSBL->Sumw2();
  TH2F* hptH3L_dDcaSBL= new TH2F("hptH3L_dDcaSBL","hptH3L_dDca;p_{T};dDca",100,0,5,100,0,5);
  hptH3L_dDcaSBL->Sumw2();
  TH2F* hptH3L_piDcaSBL= new TH2F("hptH3L_piDcaSBL","hptH3L_piDca;p_{T};piDca",100,0,5,100,0,20);
  hptH3L_piDcaSBL->Sumw2();
  TH2F* hptH3L_pDcaSBL = new TH2F("hptH3L_pDcaSBL","hptH3L_pDca;p_{T};pDca",100,0,5,100,0,10);
  hptH3L_pDcaSBL->Sumw2();
  TH2F* hptH3L_dpDcaSBL = new TH2F("hptH3L_dpDcaSBL","hptH3L_dpDca;p_{T};dpDca",100,0,5,80,0,4);
  hptH3L_dpDcaSBL->Sumw2();

  TH2F* hptH3L_lSBR= new TH2F("hptH3L_lSBR","hptH3L_l;p_{T};l",100,0,5,500,0,100);
  hptH3L_lSBR->Sumw2();
  TH2F* hptH3L_ldlSBR = new TH2F("hptH3L_ldlSBR","hptH3L_ldl;p_{T};ldl",100,0,5,250,0,50);
  hptH3L_ldlSBR->Sumw2();
  TH2F* hptH3L_dchi2primSBR = new TH2F("hptH3L_dchi2primSBR","hptH3L_dchi2pri;p_{T};d chi2prim",100,0,5,100,0,20);
  hptH3L_dchi2primSBR->Sumw2();
  TH2F* hptH3L_pichi2primSBR = new TH2F("hptH3L_pichi2primSBR","hptH3L_pichi2pri;p_{T};pi chi2prim",100,0,5,100,0,20);
  hptH3L_pichi2primSBR->Sumw2();
  TH2F* hptH3L_pchi2primSBR= new TH2F("hptH3L_pchi2primSBR","hptH3L_pchi2pri;p_{T};p chi2prim",100,0,5,100,0,20);
  hptH3L_pchi2primSBR->Sumw2();
  TH2F* hptH3L_chi2topoSBR= new TH2F("hptH3L_chi2topoSBR","hptH3L_chi2topo;p_{T};chi2topo",100,0,5,100,0,10);
  hptH3L_chi2topoSBR->Sumw2();
  TH2F* hptH3L_chi2ndfSBR= new TH2F("hptH3L_chi2ndfSBR","hptH3L_chi2ndfSBR;p_{T};chi2ndf",100,0,5,100,0,10);
  hptH3L_chi2ndfSBR->Sumw2();
  TH2F* hptH3L_dDcaSBR= new TH2F("hptH3L_dDcaSBR","hptH3L_dDca;p_{T};dDca",100,0,5,100,0,5);
  hptH3L_dDcaSBR->Sumw2();
  TH2F* hptH3L_piDcaSBR= new TH2F("hptH3L_piDcaSBR","hptH3L_piDca;p_{T};piDca",100,0,5,100,0,20);
  hptH3L_piDcaSBR->Sumw2();
  TH2F* hptH3L_pDcaSBR = new TH2F("hptH3L_pDcaSBR","hptH3L_pDca;p_{T};pDca",100,0,5,100,0,10);
  hptH3L_pDcaSBR->Sumw2();
  TH2F* hptH3L_dpDcaSBR = new TH2F("hptH3L_dpDcaSBR","hptH3L_dpDca;p_{T};dpDca",100,0,5,80,0,4);
  hptH3L_dpDcaSBR->Sumw2();

  TH3F* h3H3L_v02xySig= new TH3F("h3H3L_v02xySig","h3H3L_v02xy;d Dca;v02;pi Dca",100,0,10,500,0,20, 100, 0, 10);
  h3H3L_v02xySig->Sumw2();

  TH3F* h3H3L_v02xySBR= new TH3F("h3H3L_v02xySBR","h3H3L_v02xy;d Dca;v02;pi Dca",100,0,10,500,0,20, 100, 0, 10);
  h3H3L_v02xySBR->Sumw2();

  TH3F* h3H3L_v02xySig2= new TH3F("h3H3L_v02xySig2","h3H3L_v02xy;chi2topo;v02;chi2ndf",100,0,10,500,0,20, 100, 0, 10);
  h3H3L_v02xySig2->Sumw2();
  TH3F* h3H3L_v02xySig3= new TH3F("h3H3L_v02xySig3","h3H3L_v02xy;chi2topo;v02;d dca",100,0,10,500,0,20, 100, 0, 10);
  h3H3L_v02xySig3->Sumw2();
  TH3F* h3H3L_v02xySig4= new TH3F("h3H3L_v02xySig4","h3H3L_v02xy;chi2topo;v02;pi Dca",100,0,10,500,0,20, 100, 0, 10);
  h3H3L_v02xySig4->Sumw2();
  TH3F* h3H3L_v02xySig5= new TH3F("h3H3L_v02xySig5","h3H3L_v02xy;chi2ndf;v02;d Dca",100,0,10,500,0,20, 100, 0, 10);
  h3H3L_v02xySig4->Sumw2();
  TH3F* h3H3L_v02xySig6= new TH3F("h3H3L_v02xySig6","h3H3L_v02xy;chi2ndf;v02;pi Dca",100,0,10,500,0,20, 100, 0, 10);
  h3H3L_v02xySig4->Sumw2();

  TH3F* h3H3L_v02xySBR2= new TH3F("h3H3L_v02xySBR2","h3H3L_v02xy;chi2topo;v02;chi2ndf",100,0,10,500,0,20, 100, 0, 10);
  h3H3L_v02xySBR2->Sumw2();
  TH3F* h3H3L_v02xySBR3= new TH3F("h3H3L_v02xySBR3","h3H3L_v02xy;chi2topo;v02;pi dca",100,0,10,500,0,20, 100, 0, 10);
  h3H3L_v02xySBR3->Sumw2();
  TH3F* h3H3L_v02xySBR4= new TH3F("h3H3L_v02xySBR4","h3H3L_v02xy;chi2topo;v02;d Dca",100,0,10,500,0,20, 100, 0, 10);
  h3H3L_v02xySBR4->Sumw2();
  TH3F* h3H3L_v02xySBR5= new TH3F("h3H3L_v02xySBR5","h3H3L_v02xy;chi2ndf;v02;d Dca",100,0,10,500,0,20, 100, 0, 10);
  h3H3L_v02xySBR5->Sumw2();
  TH3F* h3H3L_v02xySBR6= new TH3F("h3H3L_v02xySBR6","h3H3L_v02xy;chi2ndf;v02;pi Dca",100,0,10,500,0,20, 100, 0, 10);
  h3H3L_v02xySBR6->Sumw2();


  TH3F* h3H3L_v12xySig= new TH3F("h3H3L_v12xySig","h3H3L_v12xy;d Dca;v12;p Dca",100,0,10,500,0,20, 100, 0, 10);
  h3H3L_v12xySig->Sumw2();

  TH3F* h3H3L_v12xySBR= new TH3F("h3H3L_v12xySBR","h3H3L_v12xy;d Dca;v12;p Dca",100,0,10,500,0,20, 100, 0, 10);
  h3H3L_v12xySBR->Sumw2();

  TH3F* h3H3L_v12xySig2= new TH3F("h3H3L_v12xySig2","h3H3L_v12xy;chi2topo;v12;chi2ndf",100,0,10,500,0,20, 100, 0, 10);
  h3H3L_v12xySig2->Sumw2();
  TH3F* h3H3L_v12xySig3= new TH3F("h3H3L_v12xySig3","h3H3L_v12xy;chi2topo;v12;p dca",100,0,10,500,0,20, 100, 0, 10);
  h3H3L_v12xySig3->Sumw2();
  TH3F* h3H3L_v12xySig4= new TH3F("h3H3L_v12xySig4","h3H3L_v12xy;chi2topo;v12;d Dca",100,0,10,500,0,20, 100, 0, 10);
  h3H3L_v12xySig4->Sumw2();
  TH3F* h3H3L_v12xySig5= new TH3F("h3H3L_v12xySig5","h3H3L_v12xy;chi2ndf;v12;p Dca",100,0,10,500,0,20, 100, 0, 10);
  h3H3L_v12xySig4->Sumw2();
  TH3F* h3H3L_v12xySig6= new TH3F("h3H3L_v12xySig6","h3H3L_v12xy;chi2ndf;v12;d Dca",100,0,10,500,0,20, 100, 0, 10);
  h3H3L_v12xySig4->Sumw2();

  TH3F* h3H3L_v12xySBR2= new TH3F("h3H3L_v12xySBR2","h3H3L_v12xy;chi2topo;v12;chi2ndf",100,0,10,500,0,20, 100, 0, 10);
  h3H3L_v12xySBR2->Sumw2();
  TH3F* h3H3L_v12xySBR3= new TH3F("h3H3L_v12xySBR3","h3H3L_v12xy;chi2topo;v12;p dca",100,0,10,500,0,20, 100, 0, 10);
  h3H3L_v12xySBR3->Sumw2();
  TH3F* h3H3L_v12xySBR4= new TH3F("h3H3L_v12xySBR4","h3H3L_v12xy;chi2topo;v12;d Dca",100,0,10,500,0,20, 100, 0, 10);
  h3H3L_v12xySBR4->Sumw2();
  TH3F* h3H3L_v12xySBR5= new TH3F("h3H3L_v12xySBR5","h3H3L_v12xy;chi2ndf;v12;d Dca",100,0,10,500,0,20, 100, 0, 10);
  h3H3L_v12xySBR5->Sumw2();
  TH3F* h3H3L_v12xySBR6= new TH3F("h3H3L_v12xySBR6","h3H3L_v12xy;chi2ndf;v12;p Dca",100,0,10,500,0,20, 100, 0, 10);
  h3H3L_v12xySBR6->Sumw2();

  TH2F* hptH3L_ppimass = new TH2F("hptH3L_ppimass","hptH3L_ppimass;p_{T};mass",100,0,5,400,1.05,1.15);
  hptH3L_ppimass->Sumw2();
  TH2F* hptH3L_ppichi2prim= new TH2F("hptH3L_ppichi2prim","hptH3L_ppichi2prim;p_{T};(p#pi) #chi^{2}_{prim}",100,0,5,50,0,40);
  hptH3L_ppichi2prim->Sumw2();
  TH2F* hptH3L_ppil= new TH2F("hptH3L_ppil","hptH3L_ppil;p_{T};(p#pi) l",100,0,5,50,0,40);
  hptH3L_ppil->Sumw2();
  TH2F* hptH3L_ppildl= new TH2F("hptH3L_ppildl","hptH3L_ppildl;p_{T};(p#pi) ldl",100,0,5,50,0,25);
  hptH3L_ppildl->Sumw2();
  TH2F* hptH3L_ppichi2ndf= new TH2F("hptH3L_ppichi2ndf","hptH3L_ppichi2ndf;p_{T};(p#pi) #chi^{2}_{ndf}",100,0,5,80,0,0.5);
  hptH3L_ppichi2ndf->Sumw2();
  TH2F* hptH3L_ppi_d_DCA= new TH2F("hptH3L_ppi_d_DCA","hptH3L_ppi_d_DCA;p_{T};(p#pi-d) DCA",100,0,5,80,0,8);
  hptH3L_ppi_d_DCA->Sumw2();

  TH2F* hptH3L_ppimassSBL = new TH2F("hptH3L_ppimassSBL","hptH3L_ppimass;p_{T};mass",100,0,5,400,1.05,1.15);
  hptH3L_ppimassSBL->Sumw2();
  TH2F* hptH3L_ppichi2primSBL= new TH2F("hptH3L_ppichi2primSBL","hptH3L_ppichi2prim;p_{T};(p#pi) #chi^{2}_{prim}",100,0,5,50,0,40);
  hptH3L_ppichi2primSBL->Sumw2();
  TH2F* hptH3L_ppilSBL= new TH2F("hptH3L_ppilSBL","hptH3L_ppil;p_{T};(p#pi) l",100,0,5,50,0,40);
  hptH3L_ppilSBL->Sumw2();
  TH2F* hptH3L_ppildlSBL= new TH2F("hptH3L_ppildlSBL","hptH3L_ppildl;p_{T};(p#pi) ldl",100,0,5,50,0,25);
  hptH3L_ppildlSBL->Sumw2();
  TH2F* hptH3L_ppichi2ndfSBL= new TH2F("hptH3L_ppichi2ndfSBL","hptH3L_ppichi2ndf;p_{T};(p#pi) #chi^{2}_{ndf}",100,0,5,80,0,0.5);
  hptH3L_ppichi2ndfSBL->Sumw2();
  TH2F* hptH3L_ppi_d_DCASBL= new TH2F("hptH3L_ppi_d_DCASBL","hptH3L_ppi_d_DCA;p_{T};(p#pi-d) DCA",100,0,5,80,0,8);
  hptH3L_ppi_d_DCASBL->Sumw2();

  TH2F* hptH3L_ppimassSBR = new TH2F("hptH3L_ppimassSBR","hptH3L_ppimass;p_{T};mass",100,0,5,400,1.05,1.15);
  hptH3L_ppimassSBR->Sumw2();
  TH2F* hptH3L_ppichi2primSBR= new TH2F("hptH3L_ppichi2primSBR","hptH3L_ppichi2prim;p_{T};(p#pi) #chi^{2}_{prim}",100,0,5,50,0,40);
  hptH3L_ppichi2primSBR->Sumw2();
  TH2F* hptH3L_ppilSBR= new TH2F("hptH3L_ppilSBR","hptH3L_ppil;p_{T};(p#pi) l",100,0,5,50,0,40);
  hptH3L_ppilSBR->Sumw2();
  TH2F* hptH3L_ppildlSBR= new TH2F("hptH3L_ppildlSBR","hptH3L_ppildl;p_{T};(p#pi) ldl",100,0,5,50,0,25);
  hptH3L_ppildlSBR->Sumw2();
  TH2F* hptH3L_ppichi2ndfSBR= new TH2F("hptH3L_ppichi2ndfSBR","hptH3L_ppichi2ndf;p_{T};(p#pi) #chi^{2}_{ndf}",100,0,5,80,0,0.5);
  hptH3L_ppichi2ndfSBR->Sumw2();
  TH2F* hptH3L_ppi_d_DCASBR= new TH2F("hptH3L_ppi_d_DCASBR","hptH3L_ppi_d_DCA;p_{T};(p#pi-d) DCA",100,0,5,80,0,8);
  hptH3L_ppi_d_DCASBR->Sumw2();

  TH2F* hptH3L_ppimassSig = new TH2F("hptH3L_ppimassSig","hptH3L_ppimass;p_{T};mass",100,0,5,400,1.05,1.15);
  hptH3L_ppimassSig->Sumw2();
  TH2F* hptH3L_ppichi2primSig= new TH2F("hptH3L_ppichi2primSig","hptH3L_ppichi2prim;p_{T};(p#pi) #chi^{2}_{prim}",100,0,5,50,0,40);
  hptH3L_ppichi2primSig->Sumw2();
  TH2F* hptH3L_ppilSig= new TH2F("hptH3L_ppilSig","hptH3L_ppil;p_{T};(p#pi) l",100,0,5,50,0,40);
  hptH3L_ppilSig->Sumw2();
  TH2F* hptH3L_ppildlSig= new TH2F("hptH3L_ppildlSig","hptH3L_ppildl;p_{T};(p#pi) ldl",100,0,5,50,0,25);
  hptH3L_ppildlSig->Sumw2();
  TH2F* hptH3L_ppichi2ndfSig= new TH2F("hptH3L_ppichi2ndfSig","hptH3L_ppichi2ndf;p_{T};(p#pi) #chi^{2}_{ndf}",100,0,5,80,0,0.5);
  hptH3L_ppichi2ndfSig->Sumw2();
  TH2F* hptH3L_ppi_d_DCASig= new TH2F("hptH3L_ppi_d_DCASig","hptH3L_ppi_d_DCASig;p_{T};(p#pi-d) DCA",100,0,5,80,0,8);
  hptH3L_ppi_d_DCASig->Sumw2();

  TH3F* h3H3L_lSBL= new TH3F("h3H3L_lSBL","hptH3L_l;p_{T};l;y",100,0,5,500,0,100,100, -1., 0);
  h3H3L_lSBL->Sumw2();
  TH3F* h3H3L_ldlSBL = new TH3F("h3H3L_ldlSBL","hptH3L_ldl;p_{T};ldl;y",100,0,5,250,0,50, 100, -1., 0);
  h3H3L_ldlSBL->Sumw2();
  TH3F* h3H3L_chi2topoSBL= new TH3F("h3H3L_chi2topoSBL","hptH3L_chi2topo;p_{T};chi2topo;y",100,0,5,100,0,10, 100, -1., 0);
  h3H3L_chi2topoSBL->Sumw2();
  TH3F* h3H3L_chi2ndfSBL= new TH3F("h3H3L_chi2ndfSBL","hptH3L_chi2ndfSBL;p_{T};chi2ndf;y",100,0,5,100,0,10, 100, -1., 0);
  h3H3L_chi2ndfSBL->Sumw2();

  TH3F* h3H3L_lSBR= new TH3F("h3H3L_lSBR","hptH3L_l;p_{T};l;y",100,0,5,500,0,100,100, -1., 0);
  h3H3L_lSBR->Sumw2();
  TH3F* h3H3L_ldlSBR = new TH3F("h3H3L_ldlSBR","hptH3L_ldl;p_{T};ldl;y",100,0,5,250,0,50, 100, -1., 0);
  h3H3L_ldlSBR->Sumw2();
  TH3F* h3H3L_chi2topoSBR= new TH3F("h3H3L_chi2topoSBR","hptH3L_chi2topo;p_{T};chi2topo;y",100,0,5,100,0,10, 100, -1., 0);
  h3H3L_chi2topoSBR->Sumw2();
  TH3F* h3H3L_chi2ndfSBR= new TH3F("h3H3L_chi2ndfSBR","hptH3L_chi2ndfSBR;p_{T};chi2ndf;y",100,0,5,100,0,10, 100, -1., 0);
  h3H3L_chi2ndfSBR->Sumw2();

  TH3F* h3H3L_lSig= new TH3F("h3H3L_lSig","hptH3L_l;p_{T};l;y",100,0,5,500,0,100, 100, -1., 0);
  h3H3L_lSig->Sumw2();
  TH3F* h3H3L_ldlSig = new TH3F("h3H3L_ldlSig","hptH3L_ldl;p_{T};ldl;y",100,0,5,250,0,50,100, -1., 0);
  h3H3L_ldlSig->Sumw2();
  TH3F* h3H3L_chi2ndfSig= new TH3F("h3H3L_chi2ndfSig","hptH3L_chi2ndf;p_{T};chi2ndf;y",100,0,5,100,0,10, 100, -1., 0);
  h3H3L_chi2ndfSig->Sumw2();
  TH3F* h3H3L_chi2topoSig= new TH3F("h3H3L_chi2topoSig","hptH3L_chi2topo;p_{T};chi2topo;y",100,0,5,100,0,10, 100, -1., 0);
  h3H3L_chi2topoSig->Sumw2();

  TH3F* h3H3L_l= new TH3F("h3H3L_l","hptH3L_l;p_{T};l;y",100,0,5,500,0,100, 100, -1., 0);
  h3H3L_l->Sumw2();
  TH3F* h3H3L_ldl = new TH3F("h3H3L_ldl","hptH3L_ldl;p_{T};ldl;y",100,0,5,250,0,50,100, -1., 0);
  h3H3L_ldl->Sumw2();
  TH3F* h3H3L_chi2ndf= new TH3F("h3H3L_chi2ndf","hptH3L_chi2ndf;p_{T};chi2ndf;y",100,0,5,100,0,10, 100, -1., 0);
  h3H3L_chi2ndf->Sumw2();

  
  TH3F* h3H3L_chi2topo= new TH3F("h3H3L_chi2topo","hptH3L_chi2topo;p_{T};chi2topo;y",100,0,5,100,0,10, 100, -1., 0);
  h3H3L_chi2topo->Sumw2();


  TH3F* hH3LPPiMassPt = new TH3F("hH3LPPiMassPt","hH3LPPiMassPt;p_{T};H3L mass;PiP",100,0,5,200,2.95,3.05, 200,1.05,1.15);
  hH3LPPiMassPt->Sumw2();

  ////////////scan cut//////////
  TH3F* hH3LMassPtYTopoCut[7];
  TH3F* h3H3L_chi2ndfSig_toposcan[7];
  TH3F* h3H3L_chi2ndf_toposcan[7];
  TH3F* h3H3L_lSig_toposcan[7];
  TH3F* h3H3L_l_toposcan[7];
  TH3F* h3H3L_ldlSig_toposcan[7];
  TH3F* h3H3L_ldl_toposcan[7];

  for (int i=0;i<7;i++)
  {
      hH3LMassPtYTopoCut[i] = new TH3F(Form("hH3LMassPtYTopoCut%d",i),"hH3LMassPtY;p_{T};H3L mass;Rapidity",100,0,5,200,2.95,3.05, 100, -1.,0);
      hH3LMassPtYTopoCut[i]->Sumw2();
      h3H3L_chi2ndfSig_toposcan[i]= new TH3F(Form("h3H3L_chi2ndfTopoCut%dSig",i),"hptH3L_chi2ndf;p_{T};chi2ndf;y",100,0,5,100,0,10, 100, -1, 0);
      h3H3L_chi2ndfSig_toposcan[i]->Sumw2();
      h3H3L_chi2ndf_toposcan[i]= new TH3F(Form("h3H3L_chi2ndfTopoCut%d",i),"hptH3L_chi2ndf;p_{T};chi2ndf;y",100,0,5,100,0,10,100,-1,0);
      h3H3L_chi2ndf_toposcan[i]->Sumw2();
      h3H3L_ldlSig_toposcan[i]=new TH3F(Form("h3H3L_ldlTopoCut%dSig", i),"hptH3L_ldl;p_{T};ldl;y",100,0,5,250,0,50,100, -1., 0); 
      h3H3L_ldlSig_toposcan[i]->Sumw2();
      h3H3L_ldl_toposcan[i]=new TH3F(Form("h3H3L_ldlTopoCut%d", i),"hptH3L_ldl;p_{T};ldl;y",100,0,5,250,0,50,100, -1., 0); 
      h3H3L_ldl_toposcan[i]->Sumw2();
      h3H3L_lSig_toposcan[i]= new TH3F(Form("h3H3L_lTopoCut%dSig", i),"hptH3L_l;p_{T};l;y",100,0,5,500,0,100, 100, -1., 0);
      h3H3L_lSig_toposcan[i]->Sumw2();
      h3H3L_l_toposcan[i]= new TH3F(Form("h3H3L_lTopoCut%d", i),"hptH3L_l;p_{T};l;y",100,0,5,500,0,100, 100, -1., 0);
      h3H3L_l_toposcan[i]->Sumw2();
  }
  TH3F* hH3LMassPtYNDFCut[10];
  TH3F* h3H3L_chi2topoSig_ndfscan[10];
  TH3F* h3H3L_chi2topo_ndfscan[10];
  TH3F* h3H3L_lSig_ndfscan[10];
  TH3F* h3H3L_l_ndfscan[10];
  TH3F* h3H3L_ldlSig_ndfscan[10];
  TH3F* h3H3L_ldl_ndfscan[10];

  for (int i=0;i<10;i++)
  {
      hH3LMassPtYNDFCut[i] = new TH3F(Form("hH3LMassPtYNDFCut%d",i),"hH3LMassPtY;p_{T};H3L mass;Rapidity",100,0,5,200,2.95,3.05, 100, -1.,0);
      hH3LMassPtYNDFCut[i]->Sumw2();
      h3H3L_chi2topoSig_ndfscan[i]= new TH3F(Form("h3H3L_chi2topoNDFCut%dSig",i),"hptH3L_chi2topo;p_{T};chi2topo;y",100,0,5,120,0,15, 100, -1, 0);
      h3H3L_chi2topoSig_ndfscan[i]->Sumw2();
      h3H3L_chi2topo_ndfscan[i]= new TH3F(Form("h3H3L_chi2topoNDFCut%d",i),"hptH3L_chi2topo;p_{T};chi2topo;y",100,0,5,120,0,15,100,-1,0);
      h3H3L_chi2topo_ndfscan[i]->Sumw2();
      h3H3L_ldl_ndfscan[i]=new TH3F(Form("h3H3L_ldlNDFCut%d", i),"hptH3L_ldl;p_{T};ldl;y",100,0,5,250,0,50,100, -1., 0); 
      h3H3L_ldl_ndfscan[i]->Sumw2();
      h3H3L_ldlSig_ndfscan[i]=new TH3F(Form("h3H3L_ldlNDFCut%dSig", i),"hptH3L_ldl;p_{T};ldl;y",100,0,5,250,0,50,100, -1., 0); 
      h3H3L_ldlSig_ndfscan[i]->Sumw2();
      h3H3L_l_ndfscan[i]= new TH3F(Form("h3H3L_lNDFCut%d", i),"hptH3L_l;p_{T};l;y",100,0,5,500,0,100, 100, -1., 0);
      h3H3L_l_ndfscan[i]->Sumw2();
      h3H3L_lSig_ndfscan[i]= new TH3F(Form("h3H3L_lNDFCut%dSig", i),"hptH3L_l;p_{T};l;y",100,0,5,500,0,100, 100, -1., 0);
      h3H3L_lSig_ndfscan[i]->Sumw2();
  }
  int const nsys=16;
  TH3F* hH3LMassPtYSysCut[nsys];
  TH3F* h3H3L_chi2ndfSig_sysscan[nsys];
  TH3F* h3H3L_chi2topoSig_sysscan[nsys];
  TH3F* h3H3L_ldlSig_sysscan[nsys];
  TH3F* h3H3L_lSig_sysscan[nsys];
  TH3F* h3H3L_chi2ndf_sysscan[nsys];
  TH3F* h3H3L_chi2topo_sysscan[nsys];
  TH3F* h3H3L_l_sysscan[nsys];
  TH3F* h3H3L_ldl_sysscan[nsys];
  for (int i=0;i<nsys;i++)
  {
      hH3LMassPtYSysCut[i] = new TH3F(Form("hH3LMassPtYSysCut%d",i),"hH3LMassPtY;p_{T};H3L mass;Rapidity",100,0,5,200,2.95,3.05, 100, -1.,0);
      hH3LMassPtYSysCut[i]->Sumw2();
      h3H3L_chi2ndfSig_sysscan[i]= new TH3F(Form("h3H3L_chi2ndfSysCut%dSig",i),"hptH3L_chi2ndf;p_{T};chi2ndf;y",100,0,5, 100,0,10, 100, -1, 0);
      h3H3L_chi2ndfSig_sysscan[i]->Sumw2();
      h3H3L_chi2ndf_sysscan[i]= new TH3F(Form("h3H3L_chi2ndfSysCut%d",i),"hptH3L_chi2ndf;p_{T};chi2ndf;y",100,0,5, 100,0,10, 100, -1, 0);
      h3H3L_chi2ndf_sysscan[i]->Sumw2();

      h3H3L_chi2topo_sysscan[i]= new TH3F(Form("h3H3L_chi2topoSysCut%d",i),"hptH3L_chi2topo;p_{T};chi2topo;y",100,0,5,120,0,15,100,-1,0);
      h3H3L_chi2topo_sysscan[i]->Sumw2();
      h3H3L_chi2topoSig_sysscan[i]= new TH3F(Form("h3H3L_chi2topoSysCut%dSig",i),"hptH3L_chi2topo;p_{T};chi2topo;y",100,0,5,120,0,15,100,-1,0);
      h3H3L_chi2topoSig_sysscan[i]->Sumw2();
      
      h3H3L_lSig_sysscan[i] = new TH3F(Form("h3H3L_lSysCut%dSig", i),"hptH3L_l;p_{T};l;y",100,0,5,500,0,100, 100, -1., 0);
      h3H3L_lSig_sysscan[i]->Sumw2();
      h3H3L_l_sysscan[i] = new TH3F(Form("h3H3L_lSysCut%d", i),"hptH3L_l;p_{T};l;y",100,0,5,500,0,100, 100, -1., 0);
      h3H3L_l_sysscan[i]->Sumw2();
      
      h3H3L_ldlSig_sysscan[i] = new TH3F(Form("h3H3L_ldlSysCut%dSig", i),"hptH3L_ldl;p_{T};ldl;y",100,0,5,250,0,50,100, -1., 0);
      h3H3L_ldlSig_sysscan[i]->Sumw2();
      h3H3L_ldl_sysscan[i] = new TH3F(Form("h3H3L_ldlSysCut%d", i),"hptH3L_ldl;p_{T};ldl;y",100,0,5,250,0,50,100, -1., 0);
      h3H3L_ldl_sysscan[i]->Sumw2();
 }
  ////////////////////////////////////////////////////////////////

  Long64_t n_Entries = htriton3_tree.GetEntries();
  cout <<"start process "<< n_Entries<<" events" << endl;

  for (int i=0;i<n_Entries;i++)
  {
    htriton3_tree.GetEntry(i); 
    if (i%100000==0) cout <<"read "<<i<<" events!" << endl;
    if ( !(bismc == mcState)  ) continue;
    if (bisMix!=isMix) continue;
    if (cent9<centLow || cent9>centHigh) continue; //remove 50-80%

    double ptweight = 1; //reserved
    double rapweight = 1;
    /* mcMotherPt =sqrt(bmcpx*bmcpx+bmcpy*bmcpy) ; */

    /* cout << bmcrap<< endl; */
    /* double mcweight = 1./g_pt_fine_in->GetBinContent(g_pt_fine_in->FindBin(bmcrap,sqrt(bmcpx*bmcpx+bmcpy*bmcpy)));; */
    double mcweight = 1;
    double corrweight = 1;
    if (mcState!=0) {
      TLorentzVector mcptc;
      if (mode==0 && mcState==1 ) 
        mcptc.SetXYZM(bmcpx,bmcpy,bmcpz,mass_ht);
      else if (mode==0 && mcState==-20 ) {
        TLorentzVector pion;
        pion.SetXYZM(b0mcpx, b0mcpy, b0mcpz, mass_pi );
        // pion.SetXYZM(bpionpx, bpionpy, bpionpz, mass_pi );
        TLorentzVector proton;
        proton.SetXYZM( b1mcpx, b1mcpy, b1mcpz, mass_p );
        // proton.SetXYZM( bprotonpx, bprotonpy, bprotonpz, mass_p );
        mcptc = pion+proton; //using lambda pt
        TLorentzVector deu;
        deu.SetXYZM(bdpx, bdpy, bdpz, mass_d);
        /* cout <<mcptc.Pt() << endl; */
        TLorentzVector H3LR = mcptc + deu;
        TLorentzVector Qvect = (mcptc-deu);
        
        if (applycorr) {
          double Pinv = H3LR.Mag();
          double Q1 = (mass_ld*mass_ld-mass_d*mass_d)/Pinv;
          double Q=sqrt(Q1*Q1-Qvect.Mag2());
          double kstar = Q/2.0;
          if (kstar<0.075 && kstar>0.003) corrweight=gdLdCorr->Eval(kstar*1e3);
          else if (kstar<0.003) corrweight=gdLdCorr->Eval(3);
          else if (kstar>0.075) corrweight=1;
          if (corrweight<1) corrweight=1;
        }
      }
      else if (mode==1 && mcState==1) 
        mcptc.SetXYZM(bmcpx,bmcpy,bmcpz,mass_ld);
      double bmcrap = -1*(mcptc.Rapidity() - ycm);
      /* hPhase->Fill(bmcrap, mcptc.Pt()); */
      double mcMotherPt = mcptc.Pt();
      // rapweight = t_quadr->Eval(bmcrap);
      rapweight = 1;
      double mccounts=0;
       if ((mode==0&&mcState==-20) || (mode==1 && mcState==1)) 
         mccounts  = g_pt_fine_in->GetBinContent(g_pt_fine_in->GetXaxis()->FindBin(-1*bmcrap), g_pt_fine_in->GetYaxis()->FindBin(mcMotherPt));
       else if (mode==0 && mcState==1) mccounts =  g_pt_fine_in->GetBinContent(g_pt_fine_in->GetXaxis()->FindBin(bmcrap),g_pt_fine_in->GetYaxis()->FindBin(mcMotherPt));
      if (mccounts>0) mcweight = 1./mccounts;
      
      // if (mode==0 && mcState ==1) ptweight = bolt1->Eval(mcMotherPt);
      if (mode==0 && mcState ==1) {
          // ptweight = fH3Ldydpt[1]->Eval(mcMotherPt);
          if (bmcrap<0 && bmcrap>=-0.25) ptweight = fH3Ldydpt[0]->Eval(mcMotherPt);
          if (bmcrap<-0.25 && bmcrap>=-0.5) ptweight = fH3Ldydpt[1]->Eval(mcMotherPt);
          if (bmcrap<-0.5 && bmcrap>=-0.75) ptweight = fH3Ldydpt[2]->Eval(mcMotherPt);
          if (bmcrap<-0.75) ptweight = fH3Ldydpt[2]->Eval(mcMotherPt);
      }
      else if (mode==1 || mcState==-20){
        if (bmcrap<-0.7){ ptweight = levyfit8->Eval(mcMotherPt); }
        else if (bmcrap<-0.5){ ptweight = levyfit7->Eval(mcMotherPt); }
        else if (bmcrap<-0.3){ ptweight = levyfit6->Eval(mcMotherPt); }
        else if (bmcrap<-0.1){ ptweight = levyfit5->Eval(mcMotherPt); }
        else if (bmcrap<0.1) { ptweight = levyfit4->Eval(mcMotherPt); }
        else if (bmcrap<0.3) { ptweight = levyfit5->Eval(mcMotherPt); }
        else if (bmcrap<0.5) { ptweight = levyfit6->Eval(mcMotherPt); }
        else if (bmcrap<0.7) { ptweight = levyfit7->Eval(mcMotherPt); }
        else { ptweight = levyfit8->Eval(mcMotherPt); }
      }
    }
    /* cout << mcMotherPt <<" "<<mcweight<<" "<<ptweight<<" "<<rapweight<<endl; */
    gweight=reweight; //centrality weight
    // gweight=1; //centrality weight

    double weight = ptweight*rapweight*mcweight*gweight*corrweight;
    // if (mcState==0) weight = reweight; 
    if (mcState==0) weight = 1; 
    double ppi_pt = sqrt(bpx*bpx+bpy*bpy); // lambda case 
    if (mode==0) ppi_pt= sqrt(bpionpx*bpionpx+bpionpy*bpionpy+bprotonpy*bprotonpy+bprotonpx*bprotonpx);
    if (mode==0) {
      if (bparticlemass<2.95 || bparticlemass > 3.05) continue;
      if (mass_01>1.15 || mass_01< 1.05) continue;
    }
    hcharge->Fill(bparticleid>0?1:-1);
    if (bparticleid<0) continue; //currently only look at particle

    //compare H3L->ppi and Lambda->ppi
    if ((mode==0 && bismc==mcState) || (mode==1 && bismc==mcState) ){ 
      // cout << reweight << " "<<ppi_pt<< " "<<dca_pion<<" "<<chi2primary_proton<<endl;
      // hptpdca->Fill( ppi_pt,dca_proton, weight);
      // hptpidca->Fill( ppi_pt,dca_pion, weight);
      // hptsumdca->Fill( ppi_pt,dca_pion+dca_proton, weight);
      // hptpchi2prim->Fill( ppi_pt,chi2primary_proton, weight);
      // hptpichi2prim->Fill( ppi_pt,chi2primary_pi, weight);
      // hptppichi2ndf->Fill( ppi_pt, v_01_chi2ndf, weight);
      // hptppichi2prim->Fill( ppi_pt, v_01_chi2primary, weight);
      // hptppil->Fill( ppi_pt, v_lambda_l_0, weight);
      // hptppildl->Fill( ppi_pt, v_lambda_ldl_0, weight);
      // hptppimass->Fill( ppi_pt, mass_01, weight);
      if (mode==0) 
      {

      }
      else if (mode==1) {
        // hptppimass->Fill( ppi_pt, bparticlemass, weight);
      }
    }

    //compare H3L or  background
    if (mode==0 ) {
      TLorentzVector H3L;
      H3L.SetXYZM(bpx, bpy, bpz, bparticlemass );
      TLorentzVector proton;
      proton.SetXYZM(bprotonpx, bprotonpy, bprotonpz, mass_p );
      TLorentzVector pion;
      pion.SetXYZM(bpionpx, bpionpy, bpionpz, mass_pi );
      TLorentzVector deuteron;
      deuteron.SetXYZM(bdpx, bdpy, bdpz, mass_d );
      TLorentzVector Ld = pion+proton;
      TLorentzVector Qvect = (Ld-deuteron);
      double Pinv = H3L.Mag();
      double Q1 = (mass_ld*mass_ld-mass_d*mass_d)/Pinv;
      double Q=sqrt(Q1*Q1-Qvect.Mag2());
      double kstar = Q/2.0;

      double H3LpT = sqrt(bpx*bpx+bpy*bpy);
      double H3Ly = -1*(H3L.Rapidity() - ycm);
      double p_pt = sqrt(bprotonpx*bprotonpx+bprotonpy*bprotonpy);
      double p_p = sqrt(bprotonpx*bprotonpx+bprotonpy*bprotonpy+bprotonpz*bprotonpz);
      double p_d = sqrt(bdpx*bdpx+bdpy*bdpy+bdpz*bdpz);
      double pi_pt = sqrt(bpionpx*bpionpx+bpionpy*bpionpy);
      double d_pt = sqrt(bdpx*bdpx+bdpy*bdpy);

      double dpvtx = sqrt( v_12_x*v_12_x+v_12_y*v_12_y );
      // bool dcaCut = (dca_deuteron<1) && pi_pt>0.15 && p_pt>0.15 && d_pt>0.15;
      bool dcaCut = 1;
      bool nhitscut = nhits_proton>15 && nhits_deuteron>15 && nhits_pion>15;
      bool MERecCut = ht_l<200 && ht_l>1 && ht_ldl>3 && ht_chi2ndf<10 && ht_chi2topo<10;
      //default cut
      // double lcut=1, ldlcut=1, chi2topocut = 10, pdcut=5, ppcut=5, chi2ndfcut=5, chi2_dcut=0, chi2_picut=1, chi2_pcut=1;
      double lcut=8, ldlcut=5, chi2topocut = 3, pdcut=5, ppcut=5, chi2ndfcut=3.5, chi2_dcut=0, chi2_picut=10, chi2_pcut=5;
      // double lcut=8, ldlcut=5, chi2topocut = 2.5, pdcut=3, ppcut=2, chi2ndfcut=1.5, chi2_dcut=0, chi2_picut=10, chi2_pcut=5;
      // double lcut=8, ldlcut=5, chi2topocut = 2., pdcut=3, ppcut=2, chi2ndfcut=3.5, chi2_dcut=0, chi2_picut=10, chi2_pcut=5;
      bool passTopoCuts =  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut;
      bool passscantopo=  ht_l >lcut && ht_ldl>ldlcut  && fabs(p_d)<pdcut && fabs(p_p)<ppcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut;
      bool Topo_l= ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut;
      bool Topo_ldl=  ht_l >lcut && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut;
      bool Topo_chi2topo = ht_l >lcut && ht_ldl>ldlcut  && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut;
      bool Topo_chi2ndf =  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut;
      bool Topo_chi2d=  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut;
      bool Topo_chi2pi=  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_proton>chi2_pcut && MERecCut && dcaCut;
      bool Topo_chi2p=  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && MERecCut && dcaCut;

     int const nmax = 20;
      bool passTopoCutsSys[nmax];
      bool passTopoCutsSys_chi2ndf[nmax];
      bool passTopoCutsSys_chi2topo[nmax];
      bool passTopoCutsSys_l[nmax];
      bool passTopoCutsSys_ldl[nmax];
      double lcut1 = 6, lcut2 = 10, ldlcut1 = 4, ldlcut2 = 6;
      // double  chi2topocut1 = 1., chi2topocut2 = 2;
      // double  chi2topocut1 = 2, chi2topocut2 = 3;
      double  chi2topocut1 = 1.5, chi2topocut2 = 2.5;
      // double pdcut1=2.5, pdcut2=3.5;
      // double ppcut1=1.5, ppcut2=2.5;
      double chi2ndfcut1=3., chi2ndfcut2=4;
      // double chi2_dcut1=3;
      double chi2_picut1=8, chi2_picut2=12; 
      double chi2_pcut1=4, chi2_pcut2=6;
      bool dcaCut2 = (dca_deuteron<1) && pi_pt>0.1 && p_pt>0.1 && d_pt>0.1;
      bool dcaCut1 = (dca_deuteron<1) && pi_pt>0.2 && p_pt>0.2 && d_pt>0.2;
      bool nhitscut1 = nhits_pion>20 && nhits_proton>20 && nhits_deuteron>20;
      bool nhitscut2 = nhits_pion>17 && nhits_proton>17 && nhits_deuteron>17;
      passTopoCutsSys[0] =  ht_l >lcut1 && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys[1] =  ht_l >lcut2 && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys[2] =  ht_l >lcut && ht_ldl>ldlcut1  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys[3] =  ht_l >lcut && ht_ldl>ldlcut2  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys[4] =  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut1 && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys[5] =  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut2 && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys[6] =  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut1 && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys[7] =  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut2 && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys[8] =  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut1 && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys[9] =  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut2 && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys[10] =  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut1 && MERecCut && dcaCut ;
      passTopoCutsSys[11] =  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut2 && MERecCut && dcaCut ;
      passTopoCutsSys[12] =  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut && nhitscut1;
      passTopoCutsSys[13] =  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut && nhitscut2;
      passTopoCutsSys[14] =  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut1 && nhitscut;
      passTopoCutsSys[15] =  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut2 && nhitscut;

      passTopoCutsSys_chi2ndf[0] =  ht_l >lcut1 && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut &&  chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_chi2ndf[1] =  ht_l >lcut2 && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut &&  chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_chi2ndf[2] =  ht_l >lcut && ht_ldl>ldlcut1  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut &&  chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_chi2ndf[3] =  ht_l >lcut && ht_ldl>ldlcut2  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut &&  chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_chi2ndf[4] =  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut1 && fabs(p_d)<pdcut && fabs(p_p)<ppcut &&  chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_chi2ndf[5] =  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut2 && fabs(p_d)<pdcut && fabs(p_p)<ppcut &&  chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_chi2ndf[6] =  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_chi2ndf[7] =  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_chi2ndf[8] =  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut &&  chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut1 && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_chi2ndf[9] =  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut &&  chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut2 && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_chi2ndf[10] =  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut &&  chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut1 && MERecCut && dcaCut ;
      passTopoCutsSys_chi2ndf[11] =  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut &&  chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut2 && MERecCut && dcaCut ;
      passTopoCutsSys_chi2ndf[12] =  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut &&  chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut && nhitscut1;
      passTopoCutsSys_chi2ndf[13] =  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut &&  chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut && nhitscut2;
      passTopoCutsSys_chi2ndf[14] =  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut &&  chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut1 && nhitscut;
      passTopoCutsSys_chi2ndf[15] =  ht_l >lcut && ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut &&  chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut2 && nhitscut;

      passTopoCutsSys_chi2topo[0] =  ht_l >lcut1 && ht_ldl>ldlcut    && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_chi2topo[1] =  ht_l >lcut2 && ht_ldl>ldlcut    && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_chi2topo[2] =  ht_l >lcut && ht_ldl>ldlcut1    && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_chi2topo[3] =  ht_l >lcut && ht_ldl>ldlcut2    && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_chi2topo[4] =  ht_l >lcut && ht_ldl>ldlcut    && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_chi2topo[5] =  ht_l >lcut && ht_ldl>ldlcut    && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_chi2topo[6] =  ht_l >lcut && ht_ldl>ldlcut    && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut1 && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_chi2topo[7] =  ht_l >lcut && ht_ldl>ldlcut    && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut2 && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_chi2topo[8] =  ht_l >lcut && ht_ldl>ldlcut    && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut1 && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_chi2topo[9] =  ht_l >lcut && ht_ldl>ldlcut    && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut2 && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_chi2topo[10] =  ht_l >lcut && ht_ldl>ldlcut    && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut1 && MERecCut && dcaCut ;
      passTopoCutsSys_chi2topo[11] =  ht_l >lcut && ht_ldl>ldlcut    && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut2 && MERecCut && dcaCut ;
      passTopoCutsSys_chi2topo[12] =  ht_l >lcut && ht_ldl>ldlcut    && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut && nhitscut1;
      passTopoCutsSys_chi2topo[13] =  ht_l >lcut && ht_ldl>ldlcut    && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut && nhitscut2;
      passTopoCutsSys_chi2topo[14] =  ht_l >lcut && ht_ldl>ldlcut    && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut1 && nhitscut;
      passTopoCutsSys_chi2topo[15] =  ht_l >lcut && ht_ldl>ldlcut    && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut2 && nhitscut;

      passTopoCutsSys_l[0] =  ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_l[1] =  ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_l[2] =  ht_ldl>ldlcut1  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_l[3] =  ht_ldl>ldlcut2  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_l[4] =  ht_ldl>ldlcut  && ht_chi2topo<chi2topocut1 && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_l[5] =  ht_ldl>ldlcut  && ht_chi2topo<chi2topocut2 && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_l[6] =  ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut1 && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_l[7] =  ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut2 && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_l[8] =  ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut1 && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_l[9] =  ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut2 && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_l[10] =  ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut1 && MERecCut && dcaCut ;
      passTopoCutsSys_l[11] =  ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut2 && MERecCut && dcaCut ;
      passTopoCutsSys_l[12] =  ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut && nhitscut1;
      passTopoCutsSys_l[13] =  ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut && nhitscut2;
      passTopoCutsSys_l[14] =  ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut1 && nhitscut;
      passTopoCutsSys_l[15] =  ht_ldl>ldlcut  && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut2 && nhitscut;

      passTopoCutsSys_ldl[0] =  ht_l >lcut1 &&  ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_ldl[1] =  ht_l >lcut2 &&  ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_ldl[2] =  ht_l >lcut &&  ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_ldl[3] =  ht_l >lcut &&  ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_ldl[4] =  ht_l >lcut && ht_chi2topo<chi2topocut1 && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_ldl[5] =  ht_l >lcut && ht_chi2topo<chi2topocut2 && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_ldl[6] =  ht_l >lcut && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut1 && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_ldl[7] =  ht_l >lcut && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut2 && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_ldl[8] =  ht_l >lcut && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut1 && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_ldl[9] =  ht_l >lcut && ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut2 && chi2primary_proton>chi2_pcut && MERecCut && dcaCut ;
      passTopoCutsSys_ldl[10] =  ht_l >lcut &&  ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut1 && MERecCut && dcaCut ;
      passTopoCutsSys_ldl[11] =  ht_l >lcut &&  ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut2 && MERecCut && dcaCut ;
      passTopoCutsSys_ldl[12] =  ht_l >lcut &&  ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut && nhitscut1;
      passTopoCutsSys_ldl[13] =  ht_l >lcut &&  ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut && nhitscut2;
      passTopoCutsSys_ldl[14] =  ht_l >lcut &&  ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut1 && nhitscut;
      passTopoCutsSys_ldl[15] =  ht_l >lcut &&  ht_chi2topo<chi2topocut && fabs(p_d)<pdcut && fabs(p_p)<ppcut && ht_chi2ndf<chi2ndfcut && chi2primary_d>chi2_dcut  && chi2primary_pi>chi2_picut && chi2primary_proton>chi2_pcut && MERecCut && dcaCut2 && nhitscut;

      // bool MERecCut = ht_l<200 && ht_l>1 && ht_ldl>3 && ht_chi2ndf<10 && ht_chi2topo<10;
      // test secondary particles
      // bool passTopoCuts =  ht_l >8 && ht_ldl>5  && ht_chi2topo<3 && fabs(p_d)<3 && fabs(p_p)<2 && ht_chi2ndf>5  && chi2primary_d>0 && chi2primary_pi>10 && chi2primary_proton>5 && MERecCut;
      // bool Topo_l =        ht_ldl>5  && ht_chi2topo<3 && ht_chi2ndf>5  && chi2primary_d>0 && chi2primary_pi>10 && chi2primary_proton>5 && fabs(p_d)<3 && fabs(p_p)<2 && MERecCut;
      // bool Topo_ldl =      ht_l>8  && ht_chi2topo<3 && ht_chi2ndf>5  && chi2primary_d>0 && chi2primary_pi>10 && chi2primary_proton>5 && fabs(p_d)<3 && fabs(p_p)<2 && MERecCut;
      // bool Topo_chi2topo = ht_l >8 && ht_ldl>5   && ht_chi2ndf>5  && chi2primary_d>0 && chi2primary_pi>10 && chi2primary_proton>5 && fabs(p_d)<3 && fabs(p_p)<2 && MERecCut;
      // bool Topo_chi2ndf =  ht_l >8 && ht_ldl>5  && ht_chi2topo<3 &&  chi2primary_d>0 && chi2primary_pi>10 && chi2primary_proton>5 && fabs(p_d)<3 && fabs(p_p)<2 && MERecCut;
      // bool Topo_chi2d=     ht_l >8 && ht_ldl>5  && ht_chi2topo<3 && ht_chi2ndf>5  && chi2primary_pi>10 && chi2primary_proton>5 && fabs(p_d)<3 && fabs(p_p)<2 && MERecCut;
      // bool Topo_chi2pi =   ht_l >8 && ht_ldl>5  && ht_chi2topo<3 && ht_chi2ndf>5  && chi2primary_d>0 && chi2primary_proton>5 && fabs(p_d)<3 && fabs(p_p)<2 && MERecCut;
      // bool Topo_chi2p =    ht_l >8 && ht_ldl>5  && ht_chi2topo<3 && ht_chi2ndf>5  && chi2primary_d>0 && chi2primary_pi>10 && fabs(p_d)<3 && fabs(p_p)<2 && MERecCut;

      //all require them match topo cuts
      // full region
      if (bparticlemass>2.95 && bparticlemass<3.05){
        if (Topo_chi2topo) {
          hptH3L_chi2topo ->Fill( H3LpT, ht_chi2topo, weight);
          h3H3L_chi2topo ->Fill( H3LpT, ht_chi2topo,H3Ly,  weight);
        }
        if (Topo_chi2ndf) {
          hptH3L_chi2ndf ->Fill( H3LpT, ht_chi2ndf, weight);
          h3H3L_chi2ndf ->Fill( H3LpT, ht_chi2ndf, H3Ly, weight);
        }
        if (Topo_l) {
          hptH3L_l ->Fill( H3LpT, ht_l, weight);
          h3H3L_l ->Fill( H3LpT, ht_l, H3Ly, weight);
        }
        if (Topo_ldl) {
          hptH3L_ldl ->Fill( H3LpT, ht_ldl, weight);
          h3H3L_ldl ->Fill( H3LpT, ht_ldl, H3Ly, weight);
        }

        if (Topo_chi2p) hptH3L_pchi2prim ->Fill( H3LpT, chi2primary_proton, weight);
        if (Topo_chi2pi) hptH3L_pichi2prim ->Fill( H3LpT, chi2primary_pi, weight);
        if (Topo_chi2d) hptH3L_dchi2prim ->Fill( H3LpT, chi2primary_d, weight);
        if (passTopoCuts) {
          hptH3L_dDca ->Fill( H3LpT, dca_deuteron, weight);
          hptH3L_pDca ->Fill( H3LpT, dca_proton, weight);
          hptH3L_piDca ->Fill( H3LpT, dca_pion, weight);
          hptH3L_dpDca ->Fill( H3LpT, v_12_dca, weight);
          hptH3L_ppi_d_DCA->Fill( H3LpT, v_012_dca, weight);
        }
      }

      //look at sideband region
      if ((bparticlemass<2.98 && bparticlemass>2.97 ) )
      {
        if (Topo_chi2topo)  {
          h3H3L_chi2topoSBL->Fill( H3LpT, ht_chi2topo, H3Ly, weight);
          hptH3L_chi2topoSBL->Fill( H3LpT, ht_chi2topo, weight);
        }
        if (Topo_chi2ndf) {
          hptH3L_chi2ndfSBL->Fill( H3LpT, ht_chi2ndf, weight);
          h3H3L_chi2ndfSBL->Fill( H3LpT, ht_chi2ndf, H3Ly, weight);
        }
        if (Topo_l) {
          hptH3L_lSBL->Fill( H3LpT, ht_l, weight);
          h3H3L_lSBL->Fill( H3LpT, ht_l, H3Ly, weight);
        }
        if (Topo_ldl) {
          hptH3L_ldlSBL->Fill( H3LpT, ht_ldl, weight);
          h3H3L_ldlSBL->Fill( H3LpT, ht_ldl, H3Ly, weight);
        }
        
        if (ht_l >8 && ht_ldl>5 && fabs(p_d)<3 && fabs(p_p)<2  && chi2primary_d>0 && chi2primary_pi>10 && chi2primary_proton>5 && MERecCut){ 
           h3H3L_v12xySBR2->Fill(ht_chi2topo, sqrt(v_12_x*v_12_x+v_12_y*v_12_y), ht_chi2ndf, weight);
           h3H3L_v02xySBR2->Fill(ht_chi2topo, sqrt(v_02_x*v_02_x+v_02_y*v_02_y), ht_chi2ndf, weight);
        }

        if (Topo_chi2p) hptH3L_pchi2primSBL->Fill( H3LpT, chi2primary_proton, weight);
        if (Topo_chi2pi) hptH3L_pichi2primSBL->Fill( H3LpT, chi2primary_pi, weight);
        if (Topo_chi2d) hptH3L_dchi2primSBL->Fill( H3LpT, chi2primary_d, weight);
        if (passTopoCuts) {
          hptH3L_dDcaSBL->Fill( H3LpT, dca_deuteron, weight);
          hptH3L_pDcaSBL->Fill( H3LpT, dca_proton, weight);
          hptH3L_piDcaSBL->Fill( H3LpT, dca_pion, weight);
          hptH3L_dpDcaSBL->Fill( H3LpT, v_12_dca, weight);
          hptH3L_ppi_d_DCASBL->Fill( H3LpT, v_012_dca, weight);
        }
      }
      else if (( bparticlemass<3.02 && bparticlemass > 3.01 ))
      {
        if (Topo_chi2topo)  {
          hptH3L_chi2topoSBR->Fill( H3LpT, ht_chi2topo, weight);
          h3H3L_chi2topoSBR->Fill( H3LpT, ht_chi2topo, H3Ly, weight);
        }
        if (ht_l >8 && ht_ldl>5 && fabs(p_d)<3 && fabs(p_p)<2  && chi2primary_d>0 && chi2primary_pi>10 && chi2primary_proton>5 && MERecCut){ 
           h3H3L_v12xySBR2->Fill(ht_chi2topo, sqrt(v_12_x*v_12_x+v_12_y*v_12_y), ht_chi2ndf, weight);
           h3H3L_v02xySBR2->Fill(ht_chi2topo, sqrt(v_02_x*v_02_x+v_02_y*v_02_y), ht_chi2ndf, weight);
        }
        if (ht_l >8 && ht_ldl>5 && fabs(p_d)<3 && fabs(p_p)<2  && chi2primary_d>0 && chi2primary_pi>10 && chi2primary_proton>5 && ht_chi2ndf<3.5 && MERecCut){ 
           h3H3L_v12xySBR3->Fill(ht_chi2topo, sqrt(v_12_x*v_12_x+v_12_y*v_12_y), dca_proton, weight);
           h3H3L_v12xySBR4->Fill(ht_chi2topo, sqrt(v_12_x*v_12_x+v_12_y*v_12_y), dca_deuteron, weight);
           h3H3L_v02xySBR3->Fill(ht_chi2topo, sqrt(v_02_x*v_02_x+v_02_y*v_02_y), dca_pion, weight);
           h3H3L_v02xySBR4->Fill(ht_chi2topo, sqrt(v_02_x*v_02_x+v_02_y*v_02_y), dca_deuteron, weight);
        }
        if (ht_l >8 && ht_ldl>5 && fabs(p_d)<3 && fabs(p_p)<2  && chi2primary_d>0 && chi2primary_pi>10 && chi2primary_proton>5 && ht_chi2topo<3. && MERecCut){ 
           h3H3L_v12xySBR5->Fill(ht_chi2ndf, sqrt(v_12_x*v_12_x+v_12_y*v_12_y), dca_deuteron, weight);
           h3H3L_v12xySBR6->Fill(ht_chi2ndf, sqrt(v_12_x*v_12_x+v_12_y*v_12_y), dca_proton, weight);
           h3H3L_v02xySBR5->Fill(ht_chi2ndf, sqrt(v_02_x*v_02_x+v_02_y*v_02_y), dca_deuteron, weight);
           h3H3L_v02xySBR6->Fill(ht_chi2ndf, sqrt(v_02_x*v_02_x+v_02_y*v_02_y), dca_pion, weight);
        }
        if (Topo_chi2ndf) {
          hptH3L_chi2ndfSBR->Fill( H3LpT, ht_chi2ndf, weight);
          h3H3L_chi2ndfSBR->Fill( H3LpT, ht_chi2ndf, H3Ly, weight);
        }
        if (Topo_l) {
          hptH3L_lSBR->Fill( H3LpT, ht_l, weight);
          h3H3L_lSBR->Fill( H3LpT, ht_l, H3Ly, weight);
        }
        if (Topo_ldl) {
          hptH3L_ldlSBR->Fill( H3LpT, ht_ldl, weight);
          h3H3L_ldlSBR->Fill( H3LpT, ht_ldl, H3Ly, weight);
        }

        if (Topo_chi2p) hptH3L_pchi2primSBR->Fill( H3LpT, chi2primary_proton, weight);
        if (Topo_chi2pi) hptH3L_pichi2primSBR->Fill( H3LpT, chi2primary_pi, weight);
        if (Topo_chi2d) hptH3L_dchi2primSBR->Fill( H3LpT, chi2primary_d, weight);
        if (passTopoCuts) {
          hptH3L_dDcaSBR->Fill( H3LpT, dca_deuteron, weight);
          hptH3L_pDcaSBR->Fill( H3LpT, dca_proton, weight);
          hptH3L_piDcaSBR->Fill( H3LpT, dca_pion, weight);
          hptH3L_dpDcaSBR->Fill( H3LpT, v_12_dca, weight);
          hptH3L_ppi_d_DCASBR->Fill( H3LpT, v_012_dca, weight);
          h3H3L_v12xySBR->Fill(dca_deuteron, sqrt(v_12_x*v_12_x+v_12_y*v_12_y), dca_proton, weight);
        }
      }
      if (bparticlemass>2.989 && bparticlemass<2.995)
      {
        if (Topo_chi2topo) {
          hptH3L_chi2topoSig->Fill( H3LpT, ht_chi2topo, weight);
          h3H3L_chi2topoSig->Fill( H3LpT, ht_chi2topo,H3Ly, weight);
        }
        if (Topo_chi2ndf) {
          hptH3L_chi2ndfSig->Fill( H3LpT, ht_chi2ndf, weight);
          h3H3L_chi2ndfSig->Fill( H3LpT, ht_chi2ndf, H3Ly, weight);
        }
        if (Topo_l) {
          hptH3L_lSig->Fill( H3LpT, ht_l, weight);
          h3H3L_lSig->Fill( H3LpT, ht_l, H3Ly, weight);
        }

        if (Topo_ldl) {
          hptH3L_ldlSig->Fill( H3LpT, ht_ldl, weight);
          h3H3L_ldlSig->Fill( H3LpT, ht_ldl, H3Ly, weight);
        }

        if (Topo_chi2p) hptH3L_pchi2primSig->Fill( H3LpT, chi2primary_proton, weight);
        if (Topo_chi2pi) hptH3L_pichi2primSig->Fill( H3LpT, chi2primary_pi, weight);
        if (Topo_chi2d) hptH3L_dchi2primSig->Fill( H3LpT, chi2primary_d, weight);
        if (passTopoCuts){
          hptH3L_dDcaSig->Fill( H3LpT, dca_deuteron, weight);
          hptH3L_pDcaSig->Fill( H3LpT, dca_proton, weight);
          hptH3L_piDcaSig->Fill( H3LpT, dca_pion, weight);
          hptH3L_dpDcaSig->Fill( H3LpT, v_12_dca, weight);
          hptH3L_ppi_d_DCASig->Fill( H3LpT, v_012_dca, weight);
          // h3H3L_v12xySig->Fill(dca_deuteron, sqrt(v_12_x*v_12_x+v_12_y*v_12_y), dca_proton, weight);
          // h3H3L_v02xySig->Fill(dca_deuteron, sqrt(v_02_x*v_02_x+v_02_y*v_02_y), dca_pion, weight);
        }
        // if (ht_l >8 && ht_ldl>5 && fabs(p_d)<3 && fabs(p_p)<2  && chi2primary_d>0 && chi2primary_pi>10 && chi2primary_proton>5 && MERecCut){ 
        //    h3H3L_v12xySig2->Fill(ht_chi2topo, sqrt(v_12_x*v_12_x+v_12_y*v_12_y), ht_chi2ndf, weight);
        //    h3H3L_v02xySig2->Fill(ht_chi2topo, sqrt(v_02_x*v_02_x+v_02_y*v_02_y), ht_chi2ndf, weight);
        // }
        // if (ht_l >8 && ht_ldl>5 && fabs(p_d)<3 && fabs(p_p)<2  && chi2primary_d>0 && chi2primary_pi>10 && chi2primary_proton>5 && ht_chi2ndf<3.5 && MERecCut){ 
        //    h3H3L_v12xySig3->Fill(ht_chi2topo, sqrt(v_12_x*v_12_x+v_12_y*v_12_y), dca_proton, weight);
        //    h3H3L_v12xySig4->Fill(ht_chi2topo, sqrt(v_12_x*v_12_x+v_12_y*v_12_y), dca_deuteron, weight);
        //    h3H3L_v02xySig3->Fill(ht_chi2topo, sqrt(v_02_x*v_02_x+v_02_y*v_02_y), dca_deuteron, weight);
        //    h3H3L_v02xySig4->Fill(ht_chi2topo, sqrt(v_02_x*v_02_x+v_02_y*v_02_y), dca_pion, weight);
        // }
        // if (ht_l >8 && ht_ldl>5 && fabs(p_d)<3 && fabs(p_p)<2  && chi2primary_d>0 && chi2primary_pi>10 && chi2primary_proton>5 && ht_chi2topo<3. && MERecCut){ 
        //    h3H3L_v12xySig5->Fill(ht_chi2ndf, sqrt(v_12_x*v_12_x+v_12_y*v_12_y), dca_deuteron, weight);
        //    h3H3L_v12xySig6->Fill(ht_chi2ndf, sqrt(v_12_x*v_12_x+v_12_y*v_12_y), dca_proton, weight);
        //    h3H3L_v02xySig5->Fill(ht_chi2ndf, sqrt(v_02_x*v_02_x+v_02_y*v_02_y), dca_deuteron, weight);
        //    h3H3L_v02xySig6->Fill(ht_chi2ndf, sqrt(v_02_x*v_02_x+v_02_y*v_02_y), dca_pion, weight);
        // }
      }

      if (passTopoCuts) 
      {
        hptH3Lmass->Fill(H3LpT, bparticlemass, weight); 
        hH3LMassPtY->Fill(H3LpT, bparticlemass, H3Ly, weight);
        if (fabs(H3Ly)<1) hH3LMassPtCent->Fill(H3LpT, bparticlemass, cent9, weight);
        if (cent9<=7 && cent9 >=4 && fabs(p_d)<3 && fabs(p_p)<2  ) 
          hH3LMassPtY_5_40->Fill(H3LpT, bparticlemass, H3Ly, weight);
        if (bparticlemass<2.995  && bparticlemass >2.989) 
        {
          // hptppimassSig->Fill( ppi_pt, mass_01, weight);
          // // hptpdcaSig->Fill( ppi_pt,dca_proton, weight);
          // // hptpidcaSig->Fill( ppi_pt,dca_pion, weight);
          // // hptsumdcaSig->Fill( ppi_pt,dca_pion+dca_proton, weight);
          // // hptpchi2primSig->Fill( ppi_pt,chi2primary_proton, weight);
          // // hptpichi2primSig->Fill( ppi_pt,chi2primary_pi, weight);
          // hptppichi2ndfSig->Fill( ppi_pt, v_01_chi2ndf, weight);
          // hptppichi2primSig->Fill( ppi_pt, v_01_chi2primary, weight);
          // hptppilSig->Fill( ppi_pt, v_lambda_l_0, weight);
          // hptppildlSig->Fill( ppi_pt, v_lambda_ldl_0, weight);
          hptH3L_ppimassSig->Fill( H3LpT, mass_01, weight);
          hptH3L_ppichi2ndfSig->Fill( H3LpT, v_01_chi2ndf, weight);
          hptH3L_ppichi2primSig->Fill( H3LpT, v_01_chi2primary, weight);
          hptH3L_ppilSig->Fill( H3LpT, v_lambda_l_0, weight);
          hptH3L_ppildlSig->Fill( H3LpT, v_lambda_ldl_0, weight);
        }
        if ( ( bparticlemass<3.02 && bparticlemass > 3. ))
        {
          // hptppimassSB->Fill( ppi_pt, mass_01, weight);
          // // hptpdcaSB->Fill( ppi_pt,dca_proton, weight);
          // // hptpidcaSB->Fill( ppi_pt,dca_pion, weight);
          // // hptsumdcaSB->Fill( ppi_pt,dca_pion+dca_proton, weight);
          // // hptpchi2primSB->Fill( ppi_pt,chi2primary_proton, weight);
          // // hptpichi2primSB->Fill( ppi_pt,chi2primary_pi, weight);
          // hptppichi2ndfSB->Fill( ppi_pt, v_01_chi2ndf, weight);
          // hptppichi2primSB->Fill( ppi_pt, v_01_chi2primary, weight);
          // hptppilSB->Fill( ppi_pt, v_lambda_l_0, weight);
          // hptppildlSB->Fill( ppi_pt, v_lambda_ldl_0, weight);
          hptH3L_ppimassSBR->Fill( H3LpT, mass_01, weight);
          hptH3L_ppichi2ndfSBR->Fill( H3LpT, v_01_chi2ndf, weight);
          hptH3L_ppichi2primSBR->Fill( H3LpT, v_01_chi2primary, weight);
          hptH3L_ppilSBR->Fill( H3LpT, v_lambda_l_0, weight);
          hptH3L_ppildlSBR->Fill( H3LpT, v_lambda_ldl_0, weight);
        }
        else if ((bparticlemass<2.98 && bparticlemass>2.97 ) )
        {
          hptH3L_ppimassSBL->Fill( H3LpT, mass_01, weight);
          hptH3L_ppichi2ndfSBL->Fill( H3LpT, v_01_chi2ndf, weight);
          hptH3L_ppichi2primSBL->Fill( H3LpT, v_01_chi2primary, weight);
          hptH3L_ppilSBL->Fill( H3LpT, v_lambda_l_0, weight);
          hptH3L_ppildlSBL->Fill( H3LpT, v_lambda_ldl_0, weight);
          hptH3L_ppi_d_DCASBL->Fill( H3LpT, v_012_dca, weight);
        }
        hH3LptProtonPt->Fill(H3LpT, p_pt , weight );
        hH3LptPionPt->Fill(H3LpT, pi_pt, weight);       
        hH3LPPiMassPt->Fill(H3LpT, bparticlemass, mass_01, weight);
      }

      if (passscantopo)
      {
        double cutstopo[7]={0.2,0.6,1,1.5,2,2.5,3};
        for (int i=0;i<7;i++) {
          if (ht_chi2topo<cutstopo[i]){ 
             if (ht_chi2ndf<chi2ndfcut) hH3LMassPtYTopoCut[i]->Fill(H3LpT, bparticlemass, H3Ly, weight);
             if (bparticlemass>2.989 && bparticlemass<2.995) {
               h3H3L_chi2ndfSig_toposcan[i]->Fill( H3LpT, ht_chi2ndf, H3Ly, weight);
               h3H3L_lSig_toposcan[i]->Fill( H3LpT, ht_l, H3Ly, weight);
               h3H3L_ldlSig_toposcan[i]->Fill( H3LpT, ht_ldl, H3Ly, weight);

             }
             if (bparticlemass>2.97 && bparticlemass<3.02) {
               h3H3L_chi2ndf_toposcan[i]->Fill( H3LpT, ht_chi2ndf, H3Ly, weight);
               h3H3L_l_toposcan[i]->Fill( H3LpT, ht_l, H3Ly, weight);
               h3H3L_ldl_toposcan[i]->Fill( H3LpT, ht_ldl, H3Ly, weight);
             }
          }
        }

        double cutsndf[10]={1, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5};
        for (int i=0;i<10;i++) {
          if (ht_chi2ndf<cutsndf[i]){ 
             if (ht_chi2topo<chi2topocut) hH3LMassPtYNDFCut[i]->Fill(H3LpT, bparticlemass, H3Ly, weight);
             if (bparticlemass>2.989 && bparticlemass<2.995) {
               h3H3L_chi2topoSig_ndfscan[i]->Fill( H3LpT, ht_chi2topo, H3Ly, weight);
               h3H3L_lSig_ndfscan[i]->Fill( H3LpT, ht_l, H3Ly, weight);
               h3H3L_ldlSig_ndfscan[i]->Fill( H3LpT, ht_ldl, H3Ly, weight);

             }
             if (bparticlemass>2.97 && bparticlemass<3.02) {
               h3H3L_chi2topo_ndfscan[i]->Fill( H3LpT, ht_chi2topo, H3Ly, weight);
               h3H3L_l_ndfscan[i]->Fill( H3LpT, ht_l, H3Ly, weight);
               h3H3L_ldl_ndfscan[i]->Fill( H3LpT, ht_ldl, H3Ly, weight);

             }
          }
        }
      }

      for (int isys=0;isys<nsys;isys++)
      {
        if (bparticlemass>2.989 && bparticlemass<2.995){
          if (passTopoCutsSys_chi2ndf[isys]) h3H3L_chi2ndfSig_sysscan[isys]->Fill( H3LpT, ht_chi2ndf, H3Ly, weight);
          if (passTopoCutsSys_chi2topo[isys]) h3H3L_chi2topoSig_sysscan[isys]->Fill( H3LpT, ht_chi2topo,H3Ly, weight);
          if (passTopoCutsSys_l[isys]) h3H3L_lSig_sysscan[isys]->Fill( H3LpT, ht_l, H3Ly, weight);
          if (passTopoCutsSys_ldl[isys]) h3H3L_ldlSig_sysscan[isys]->Fill( H3LpT, ht_ldl,H3Ly, weight);
        }
        if (bparticlemass>2.97 && bparticlemass<3.02){
          if (passTopoCutsSys_chi2ndf[isys]) h3H3L_chi2ndf_sysscan[isys]->Fill( H3LpT, ht_chi2ndf, H3Ly, weight);
          if (passTopoCutsSys_chi2topo[isys]) h3H3L_chi2topo_sysscan[isys]->Fill( H3LpT, ht_chi2topo,H3Ly, weight);
          if (passTopoCutsSys_l[isys]) h3H3L_l_sysscan[isys]->Fill( H3LpT, ht_l,H3Ly, weight);
          if (passTopoCutsSys_ldl[isys]) h3H3L_ldl_sysscan[isys]->Fill( H3LpT, ht_ldl,H3Ly, weight);
        }
        if (passTopoCutsSys[isys]) hH3LMassPtYSysCut[isys]->Fill(H3LpT, bparticlemass, H3Ly, weight);
      }
    }
  }

  TFile* fout = new TFile(outfile.Data(),"recreate");

  ///////////////////
  // hptppimass->Write(); 
  // hptppichi2ndf->Write();
  // hptpichi2prim->Write();
  // hptpchi2prim->Write();
  // hptppichi2prim->Write();
  // hptppil->Write();
  // hptppildl->Write();
  // hptpidca->Write();
  // hptpdca->Write();
  // hptsumdca->Write();
  // //////////////////
  /* hPhase->Write(); */
  hcharge->Write();

  if (mode==0)
  {
    //////////////////////
    // hptppimassSB->Write(); 
    // hptppichi2ndfSB->Write();
    // hptppichi2primSB->Write();
    // hptppilSB->Write();
    // hptppildlSB->Write();
    //
    // hptppimassSig->Write(); 
    // hptppichi2ndfSig->Write();
    // hptppichi2primSig->Write();
    // hptppilSig->Write();
    // hptppildlSig->Write();
    // //////////////
    if (fillQAplots){ 
 
    hptH3L_ppimassSBR->Write(); 
    hptH3L_ppichi2ndfSBR->Write();
    hptH3L_ppichi2primSBR->Write();
    hptH3L_ppilSBR->Write();
    hptH3L_ppildlSBR->Write();
    hptH3L_ppi_d_DCASBR->Write();

    hptH3L_ppimassSBL->Write(); 
    hptH3L_ppichi2ndfSBL->Write();
    hptH3L_ppichi2primSBL->Write();
    hptH3L_ppilSBL->Write();
    hptH3L_ppildlSBL->Write();
    hptH3L_ppi_d_DCASBL->Write();

    hptH3L_ppimassSig->Write(); 
    hptH3L_ppichi2ndfSig->Write();
    hptH3L_ppichi2primSig->Write();
    hptH3L_ppilSig->Write();
    hptH3L_ppildlSig->Write();
    hptH3L_ppi_d_DCASig->Write();

    hptH3Lmass->Write();
    hptH3L_dDca->Write();
    hptH3L_piDca->Write();
    hptH3L_pDca->Write();
    hptH3L_dpDca->Write();
    hptH3L_chi2ndf->Write();
    hptH3L_chi2topo->Write();
    hptH3L_l->Write();
    hptH3L_ldl->Write();
    hptH3L_pichi2prim->Write();
    hptH3L_pchi2prim->Write();
    hptH3L_dchi2prim->Write();

    hptH3L_ppimass->Write(); 
    hptH3L_ppichi2ndf->Write();
    hptH3L_ppichi2prim->Write();
    hptH3L_ppil->Write();
    hptH3L_ppildl->Write();
    hptH3L_ppi_d_DCA->Write();

    hptH3L_dDcaSBR->Write();
    hptH3L_piDcaSBR->Write();
    hptH3L_pDcaSBR->Write();
    hptH3L_dpDcaSBR->Write();
    hptH3L_chi2ndfSBR->Write();
    hptH3L_chi2topoSBR->Write();
    hptH3L_lSBR->Write();
    hptH3L_ldlSBR->Write();
    hptH3L_pichi2primSBR->Write();
    hptH3L_pchi2primSBR->Write();
    hptH3L_dchi2primSBR->Write();

    hptH3L_dDcaSBL->Write();
    hptH3L_piDcaSBL->Write();
    hptH3L_pDcaSBL->Write();
    hptH3L_dpDcaSBL->Write();
    hptH3L_chi2ndfSBL->Write();
    hptH3L_chi2topoSBL->Write();
    hptH3L_lSBL->Write();
    hptH3L_ldlSBL->Write();
    hptH3L_pichi2primSBL->Write();
    hptH3L_pchi2primSBL->Write();
    hptH3L_dchi2primSBL->Write();

    hptH3L_dDcaSig->Write();
    hptH3L_piDcaSig->Write();
    hptH3L_pDcaSig->Write();
    hptH3L_dpDcaSig->Write();
    hptH3L_chi2ndfSig->Write();
    hptH3L_chi2topoSig->Write();
    hptH3L_lSig->Write();
    hptH3L_ldlSig->Write();
    hptH3L_pichi2primSig->Write();
    hptH3L_pchi2primSig->Write();
    hptH3L_dchi2primSig->Write();

    hH3LptPionPt->Write();
    hH3LptProtonPt->Write();
    hH3LPPiMassPt->Write();


    // h3H3L_v02xySig->Write();
    // h3H3L_v02xySig2->Write();
    // h3H3L_v02xySig3->Write();
    // h3H3L_v02xySig4->Write();
    // h3H3L_v02xySig5->Write();
    // h3H3L_v02xySig6->Write();
    // h3H3L_v02xySBR->Write();
    // h3H3L_v02xySBR2->Write();
    // h3H3L_v02xySBR3->Write();
    // h3H3L_v02xySBR4->Write();
    // h3H3L_v02xySBR5->Write();
    // h3H3L_v02xySBR6->Write();
    //
    // h3H3L_v12xySig->Write();
    // h3H3L_v12xySig2->Write();
    // h3H3L_v12xySig3->Write();
    // h3H3L_v12xySig4->Write();
    // h3H3L_v12xySig5->Write();
    // h3H3L_v12xySig6->Write();
    // h3H3L_v12xySBR->Write();
    // h3H3L_v12xySBR2->Write();
    // h3H3L_v12xySBR3->Write();
    // h3H3L_v12xySBR4->Write();
    // h3H3L_v12xySBR5->Write();
    // h3H3L_v12xySBR6->Write();

    }


    h3H3L_chi2ndf->Write();
    h3H3L_chi2topo->Write();
    h3H3L_l->Write();
    h3H3L_ldl->Write();

    h3H3L_chi2ndfSig->Write();
    h3H3L_chi2topoSig->Write();
    h3H3L_lSig->Write();
    h3H3L_ldlSig->Write();
    
    h3H3L_chi2ndfSBR->Write();
    h3H3L_chi2topoSBR->Write();
    h3H3L_lSBR->Write();
    h3H3L_ldlSBR->Write();

    h3H3L_chi2ndfSBL->Write();
    h3H3L_chi2topoSBL->Write();
    h3H3L_lSBL->Write();
    h3H3L_ldlSBL->Write();


    hH3LMassPtY->Write();
    for (int i=0;i<10;i++) {
      hH3LMassPtYNDFCut[i]->Write();
      h3H3L_chi2topoSig_ndfscan[i]->Write();
      h3H3L_chi2topo_ndfscan[i]->Write();
      h3H3L_lSig_ndfscan[i]->Write();
      h3H3L_l_ndfscan[i]->Write();
      h3H3L_ldlSig_ndfscan[i]->Write();
      h3H3L_ldl_ndfscan[i]->Write();
    }

    for (int i=0;i<7;i++) {
      hH3LMassPtYTopoCut[i]->Write(); 
      h3H3L_chi2ndfSig_toposcan[i]->Write();
      h3H3L_chi2ndf_toposcan[i]->Write();
      h3H3L_lSig_toposcan[i]->Write();
      h3H3L_l_toposcan[i]->Write();
      h3H3L_ldlSig_toposcan[i]->Write();
      h3H3L_ldl_toposcan[i]->Write();
    }

    for (int i=0;i<nsys;i++) {
      hH3LMassPtYSysCut[i]->Write(); 
      h3H3L_chi2ndf_sysscan[i]->Write();
      h3H3L_chi2ndfSig_sysscan[i]->Write();
      h3H3L_chi2topo_sysscan[i]->Write();
      h3H3L_chi2topoSig_sysscan[i]->Write();
      h3H3L_l_sysscan[i]->Write();
      h3H3L_lSig_sysscan[i]->Write();
      h3H3L_ldl_sysscan[i]->Write();
      h3H3L_ldlSig_sysscan[i]->Write();

    }
    hH3LMassPtCent->Write();
    hH3LMassPtY_5_40->Write();
  }

  // hrefmult->Write();
  hcentwt->Write();
  hcent->Write();

  fout->Close();
  // htriton3_tree->Close();
  time.Stop();
  time.Print();
}
