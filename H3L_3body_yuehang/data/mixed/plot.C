#include "sPhenixStyle.h"
void plot()
{
  SetsPhenixStyle();
  // TFile* fLd = new TFile("fout_Lambda_MC_Cuts_0050_015pt_sys_mixpid3_ncorr.root");
  TFile* fLd = new TFile("fout_Lambda_MC_Cuts_0050_015pt_sys_mixpid3.root");
  TH3F* h3Ld = (TH3F*)fLd->Get("hH3LMassPtY");
  TH1F* hLd = (TH1F*)h3Ld->ProjectionY();
  hLd->SetDirectory(0);
  hLd->Draw();

  // TFile* fLdcor = new TFile("fout_Lambda_MC_Cuts_0050_015pt_sys_mixpid3_corr.root");
  TFile* fLdcor = new TFile("fout_Lambda_MC_Cuts_0050_015pt_sys_pid3.root");
  TH3F* h3Ldcor = (TH3F*)fLdcor->Get("hH3LMassPtY")->Clone("aaa");
  TH1F* hLdcor = (TH1F*)h3Ldcor->ProjectionY();
  hLdcor->SetDirectory(0);

  // TFile* fH3L = new TFile("rootfile/fout_H3L_data_SE_Sep02_0050_015pt_sys.root");
  // TH3F* h3H3L = (TH3F*)fH3L->Get("hH3LMassPtY")->Clone("bbb");
  // TH1F* hH3L = (TH1F*)h3H3L->ProjectionY();
  // hH3L->SetDirectory(0);

  hLdcor->Rebin(4);
  hLd->Rebin(4);


  hLdcor->Rebin(4);
  hLd->Rebin(4);

  
  hLd->Scale(hLdcor->Integral()/hLd->Integral());
  // hLdcor->Add(hLd, -1);
  // hLdcor->Divide(hLd);
  // hLd->Draw();
  // TF1* fun = new TF1( "fun", "gaus(0)", 0, 5);
  // hLdcor->Fit(fun);
  hLdcor->Draw(" p");
  hLd->Draw(" p same");
  hLd->SetMarkerColor(kRed);

}
