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

void readMc(TString mInputlist="Lambda_tree_mc.root", int const mode = 1,   TString outfile="fout_Lambda.root", int const mcState=1, int centL=0, int centH=8)
{
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
  double const mass_p = 0.93827;
  double const mass_pi = 0.13957;

  TF1* levyfit4;
  TF1* levyfit5;
  TF1* levyfit6;
  TF1* levyfit7;
  TF1* levyfit8;
  TF1* t_quadr;
  TF1* t_quadr0;
  TF1* t_quadr1;
  TF1* bolt0;
  TF1* bolt1;
  TF1* bolt2;

  TF1* fH3Ldydpt[3];
  double par[3][3]={{0.267422, 43782.6, 2.99339}, {0.200836, 92736.2,2.99339}, {0.136789, 136882, 2.99339} };
  for (int irap=0;irap<3;irap++) {
      fH3Ldydpt[irap] = new TF1(Form("fH3Ldydpt[%d]",irap), "x*[1]*exp(-(sqrt([2]*[2]+x*x)-[2])/[0])", 0,4);
      fH3Ldydpt[irap]->SetParameters(par[irap]);
  }



  //weighting function for MC
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
    // ptmptexp[] = new TF1(Form("ptmptexp[%d]",irap), "x*[1]*exp(-(sqrt([2]*[2]+x*x)-[2])/[0])", 0,4);

    // TF1* ptmptexp[3];
    // double par[3][3]={{0.267422, 43782.6, 2.99339}, {0.200836, 92736.2,2.99339}, {0.136789, 136882, 2.99339} };
    // for (int irap=0;irap<3;irap++) {
    //     ptmptexp[irap] = new TF1(Form("ptmptexp[%d]",irap), "x*[1]*exp(-(sqrt([2]*[2]+x*x)-[2])/[0])", 0,4);
    //     ptmptexp[irap]->SetParameters(par[irap]);
    // }
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
    // if (mode==0) fgpt_0= new TFile("ht_input_mc_fine_v20.root","READ");
    // if (mode==0 && mcState==1) fgpt_0= new TFile("fH3L_phase_quasi_wt.root","READ");
    TFile* fgpt_0;
    fgpt_0= new TFile("fMC_H3L_wt.root");
    TH2F* g_pt_fine_in = (TH2F*)fgpt_0->Get("hPhase")->Clone("g_pt_fine_in");
    g_pt_fine_in->Draw("colz");
    cout <<"finish book mc weighting functions." << endl;
  ///////////////////////end of MC weighting////////////////

  cout <<"start read tree!" << endl;
  TString treename;
  if (mode==1) treename = "lambda_mc_tree";
  if (mode==0) treename = "htriton_mc_tree";
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

  htriton3_tree.SetBranchAddress("bmcparticleid",&bparticleid);

  htriton3_tree.SetBranchAddress("reweight", &reweight);
  htriton3_tree.SetBranchAddress("cent9", &cent9);

  htriton3_tree.SetBranchAddress("bmcrawpx",&bmcpx);
  htriton3_tree.SetBranchAddress("bmcrawpy",&bmcpy);                                               
  htriton3_tree.SetBranchAddress("bmcrawpz",&bmcpz);
  htriton3_tree.SetBranchAddress("bmcrawl",&bmcl);
  htriton3_tree.SetBranchAddress("bmcrawpl",&bmcpl);

  TH2F* hptH3Lmass  = new TH2F("hptH3Lmass","hptH3Lmass;p_{T};H3L mass",100,0,5,200,2.95,3.05);
  hptH3Lmass->Sumw2();
  TH3F* hH3LMassPtY= new TH3F("hH3LMassPtY","hH3LMassPtY;p_{T};H3L mass;Rapidity",100,0,5,200,2.95,3.05, 100, -1.,0);
  hH3LMassPtY->Sumw2();
  TH3F* hH3LMassPtCent= new TH3F("hH3LMassPtCent","hH3LMassPtCent;p_{T};H3L mass;Centrality",100,0,5,200,2.95,3.05, 9, -0.5,8.5);
  hH3LMassPtCent->Sumw2();
  TH3F* hH3LMassPtY_5_40= new TH3F("hH3LMassPtY_5_40","hH3LMassPtY;p_{T};H3L mass;Rapidity",100,0,5,200,2.95,3.05, 100, -1.,0);
  hH3LMassPtY_5_40->Sumw2();

  TH2F* hPhase = new TH2F("hPhase", "hPhase;y;pt", 1200,-1.2,1.2, 1000,0,5 );
  TH2F* hPhase_wt = new TH2F("hPhase_wt", "hPhase_wt;y;pt", 120,-1.2,1.2, 100,0,5 );
  ////////////////////////////////////////////////////////////////

  Long64_t n_Entries = htriton3_tree.GetEntries();
  cout <<"start process "<< n_Entries<<" events" << endl;

  for (int i=0;i<n_Entries;i++)
  {
    htriton3_tree.GetEntry(i); 
    if (i%100000000==0) cout <<"read "<<i<<" events!" << endl;
    // if (cent9<=2 || cent9>8) continue; //remove 60-80%
    if (cent9<centL || cent9>centH) continue; //remove 60-80%

    double ptweight = 1; //reserved
    double rapweight = 1;

    double mcweight = 1;
    if (mcState!=0) {
      TLorentzVector mcptc;
      if (mode==0 && mcState==1 ) 
        mcptc.SetXYZM(bmcpx,bmcpy,bmcpz,mass_ht);
      else if (mode==1 && mcState==1) 
        mcptc.SetXYZM(bmcpx,bmcpy,bmcpz,mass_ld);
      double bmcrap = -1* (mcptc.Rapidity() - ycm);
      double mcMotherPt = mcptc.Pt();
      // g_pt_fine_in->Draw("colz");
      double mccounts = g_pt_fine_in->GetBinContent(g_pt_fine_in->GetXaxis()->FindBin(bmcrap),g_pt_fine_in->GetYaxis()->FindBin(mcMotherPt));
      if (mccounts>0) mcweight = 1./mccounts; 
      // rapweight = t_quadr->Eval(bmcrap);
      // if (rapweight<0) rapweight=0; 
      rapweight = 1;
      // ptweight = bolt1->Eval(mcMotherPt);
       // ptweight = fH3Ldydpt[1]->Eval(mcMotherPt);
      if (bmcrap<0 && bmcrap>=-0.25) ptweight = fH3Ldydpt[0]->Eval(mcMotherPt);
      if (bmcrap<-0.25 && bmcrap>=-0.5) ptweight = fH3Ldydpt[1]->Eval(mcMotherPt);
      if (bmcrap<-0.5 && bmcrap>=-0.75) ptweight = fH3Ldydpt[2]->Eval(mcMotherPt);
      if (bmcrap<-0.75) ptweight = fH3Ldydpt[2]->Eval(mcMotherPt);
      // if ( mcptc.Rapidity()<0 && mcptc.Rapidity()>-0.25) ptweight = ptmptexp[0]->Eval(mcMotherPt);
      // if ( mcptc.Rapidity()<-0.25 && mcptc.Rapidity()>-0.5) ptweight = ptmptexp[1]->Eval(mcMotherPt);
      // if ( mcptc.Rapidity()<-0.5 && mcptc.Rapidity()>-0.8) ptweight = ptmptexp[2]->Eval(mcMotherPt);
    }
    
    // gweight=reweight; //centrality weight
    gweight=1; //centrality weight

    double weight = ptweight*rapweight*mcweight*gweight;
    // double weight = mcweight*gweight;
    // double weight = 1;
    if (mcState==0) weight = gweight; 
    // if (mcState==0) weight = 1; 
    // double ppi_pt = sqrt(bpx*bpx+bpy*bpy); // lambda case 
    // if (mode==0) ppi_pt= sqrt(bpionpx*bpionpx+bpionpy*bpionpy+bprotonpy*bprotonpy+bprotonpx*bprotonpx);

    //compare H3L or  background
    if (mode==0  ) {
      TLorentzVector H3L;
      H3L.SetXYZM(bmcpx,bmcpy,bmcpz,mass_ht);
      double H3LpT = sqrt(bmcpx*bmcpx+bmcpy*bmcpy);
      double H3Ly = -1*(H3L.Rapidity()-ycm);
      hptH3Lmass->Fill(H3LpT, bparticlemass, weight); 
      hH3LMassPtY->Fill(H3LpT, bparticlemass, H3Ly, weight);
      if (fabs(H3Ly)<1) hH3LMassPtCent->Fill(H3LpT, bparticlemass, cent9, weight);
      hPhase->Fill( H3Ly, H3LpT);
      hPhase_wt->Fill( H3Ly, H3LpT, mcweight);

      if ( cent9<=7 && cent9 >=4  ) 
        hH3LMassPtY_5_40->Fill(H3LpT, bparticlemass, H3Ly, weight);
        // hH3LptProtonPt->Fill(H3LpT, p_pt , weight );
        // hH3LptPionPt->Fill(H3LpT, pi_pt, weight);       
        // hH3LPPiMassPt->Fill(H3LpT, bparticlemass, mass_01, weight);
    }
  }

  TFile* fout = new TFile(outfile.Data(),"recreate");

  ///////////////////
  if (mode==0)
  {
    // hH3LptPionPt->Write();
    // hH3LptProtonPt->Write();
    // hH3LPPiMassPt->Write();
    hH3LMassPtY->Write();
    hH3LMassPtCent->Write();
    hH3LMassPtY_5_40->Write();
    hPhase->Write();
    hPhase_wt->Write();
  }

  // hrefmult->Write();
  hcentwt->Write();
  hcent->Write();

  fout->Close();
  // htriton3_tree->Close();
  time.Stop();
  time.Print();
}
