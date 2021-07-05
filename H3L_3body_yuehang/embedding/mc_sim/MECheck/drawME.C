#include "sPhenixStyle.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TLine.h"

void drawME()
{

  SetsPhenixStyle();
  TCanvas* c = new TCanvas("c","c");
  TPDF* pdf = new TPDF("ME_dLd.pdf");
  pdf->Off();
  const Int_t mode = 1;
  const int rebin=40;

  TH2D *hInvMassPt;
  TH2D *hInvMassPtRotate;
  TH2D *hInvMassPtMix;
  TH2D *hInvMassPtMixRotate;
  TH1D *hMass;
  TH1D *hMass_rotate;
  TH1D *hMass_ME;
  TH1D *hMass_MErt;
  TH1D *hRatio;
  TH1D* hRatioRt;
  // TFile *fin = new TFile("dLdInvMass_ME_5bin.root");
  // TFile *fin = new TFile("dLdInvMass_ME.root");
  // TFile *fin = new TFile("dLdInvMass_ME_20bin.root");
  TFile *fin = new TFile("dLdInvMass_ME_10bin.root");
  hInvMassPt = (TH2D *)fin->Get(Form("invmass_pt_%d",mode));
  hInvMassPtRotate = (TH2D *)fin->Get(Form("invmass_rt_pt_%d",mode));
  hMass = (TH1D *)fin->Get(Form("mass_%d",mode));
  // hMass->Scale(1e-6);
  hMass->Rebin(rebin);
  hMass_rotate = (TH1D *)fin->Get(Form("mass_rotate_%d",mode));
  hMass_rotate->Rebin(rebin);
  // hMass_rotate->Scale(1e-6);
  hMass_ME = (TH1D *)fin->Get(Form("mass_ME_%d",mode));
  hMass_ME->Rebin(rebin);
  hMass_MErt = (TH1D *)fin->Get(Form("mass_ME_rotate_%d",mode));
  hMass_MErt->Rebin(rebin);
  // double scale =hMass_rotate->Integral(hMass_rotate->GetXaxis()->FindBin(3.1),hMass_rotate->GetXaxis()->FindBin(3.25)) / hMass_MErt->Integral(hMass_rotate->GetXaxis()->FindBin(3.1),hMass_rotate->GetXaxis()->FindBin(3.25));
  double scale =hMass->Integral(hMass_rotate->GetXaxis()->FindBin(3.15),hMass_rotate->GetXaxis()->FindBin(3.25)) / hMass_ME->Integral(hMass_rotate->GetXaxis()->FindBin(3.15),hMass_rotate->GetXaxis()->FindBin(3.25));
  cout <<scale <<endl;
  hMass_MErt->Scale(scale);
  hMass_ME->Scale(scale);
  hRatio = (TH1D*)hMass_ME->Clone(Form("ratio_%d",mode));    
  hRatio->Divide(hMass);    
  hRatioRt = (TH1D*)hMass_MErt->Clone(Form("ratioRt_%d",mode));    
  hRatioRt->Divide(hMass_rotate);    


  int color[2]={kRed,kBlue};
  hMass->GetXaxis()->SetTitle("M_{d#Lambda} (GeV/c^{2})");
  hMass->Draw("hist");

  hMass->SetMarkerSize(1.0);
  hMass->SetMarkerColor(color[mode]);
  hMass->SetMarkerStyle(kFullCircle);
  
  hMass_rotate->SetMarkerSize(1.0);
  hMass_rotate->SetMarkerColor(color[mode]);
  hMass_rotate->SetMarkerStyle(kOpenCircle);

  hMass_ME->SetMarkerSize(1.0);
  hMass_ME->SetMarkerColor(kRed);
  hMass_ME->SetMarkerStyle(kFullSquare);

  hMass_MErt->SetMarkerSize(1.0);
  hMass_MErt->SetMarkerColor(kRed);
  hMass_MErt->SetMarkerStyle(kOpenSquare);

  hMass->Draw("same");
  hMass_rotate->Draw("same");
  hMass_ME->Draw("same");
  hMass_MErt->Draw("same");

  TLegend *leg = new TLegend(0.55, 0.22, 0.94, 0.44);
  leg->AddEntry(hMass, "v1#neq0","pl");  
  leg->AddEntry(hMass_rotate, "v1#neq0 rotate","pl");  
  leg->AddEntry(hMass_ME, "v1#neq0 ME scale","pl");  
  leg->AddEntry(hMass_MErt, "v1#neq0 ME rotate scale","pl");  
  leg->Draw();
  addpdf(pdf);

  TH1D *h0 = new TH1D("h0","",1,2.98, 3.28);
  h0->SetMinimum(0.9);
  h0->SetMaximum(1.1);
  h0->GetXaxis()->SetTitle("M_{d#Lambda} (GeV/c^{2})");
  h0->GetYaxis()->SetTitle("ME/Original");
  h0->Draw();

  hRatio->SetMarkerSize(1.0);
  hRatio->SetMarkerColor(color[mode]);
  hRatio->SetMarkerStyle(kFullSquare);
  hRatio->SetLineColor(color[mode]);
  hRatio->Draw("same");
  hRatio->GetYaxis()->SetTitle("ME/Original");

  hRatioRt->SetMarkerSize(1.0);
  hRatioRt->SetMarkerColor(color[mode]);
  hRatioRt->SetMarkerStyle(kOpenSquare);
  hRatioRt->SetLineColor(color[mode]);
  hRatioRt->Draw("same");
  hRatioRt->GetYaxis()->SetTitle("ME-Rotate/Original");


  drawLine(2.98,1,3.18,1,1.5,2,kBlack);

  TLegend *leg1 = new TLegend(0.6, 0.22, 0.94, 0.34);
  leg1->AddEntry(hRatio, "ME w/o rotation","pl");
  leg1->AddEntry(hRatioRt, "ME w/ rotation","pl");
  leg1->Draw();
  addpdf(pdf);

  pdf->On();
  pdf->Close();
}
