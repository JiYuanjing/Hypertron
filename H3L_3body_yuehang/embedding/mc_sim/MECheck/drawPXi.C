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

void drawPXi()
{
  
  SetsPhenixStyle();
  TCanvas* c = new TCanvas("c","c");
  TPDF* pdf = new TPDF("invmass_dLd.pdf");
  pdf->Off();
  const Int_t NV = 2;
  // const Double_t v1[NV] = {0.0, 0.2, 0.4};

  TH2D *hInvMassPt[NV];
  TH2D *hInvMassPtRotate[NV];
  TH1D *hMass[NV];
  TH1D *hMass_rotate[NV];
  TH1D *hRatio[NV];
  // TFile *fin = new TFile("invmass_dLd_0_5M.root");
  TFile *fin = new TFile("invmass_dLd_3M.root");
  for(int i=0;i<NV;i++) {
    hInvMassPt[i] = (TH2D *)fin->Get(Form("invmass_pt_%d",i));
    hInvMassPtRotate[i] = (TH2D *)fin->Get(Form("invmass_rt_pt_%d",i));
    hMass[i] = (TH1D *)fin->Get(Form("mass_%d",i));
    hMass[i]->Scale(1e-6);
    hMass_rotate[i] = (TH1D *)fin->Get(Form("mass_rotate_%d",i));
    hMass_rotate[i]->Scale(1e-6);
    hMass_rotate[i]->Rebin();
    hMass[i]->Rebin();
    hRatio[i] = (TH1D*)hMass_rotate[i]->Clone(Form("ratio_%d",i));    
    hRatio[i]->Divide(hMass[i]);    
  }
  
  int color[2]={kRed,kBlue};
  hMass[0]->GetYaxis()->SetTitle("Counts (weighted)");
  hMass[0]->GetXaxis()->SetTitle("M_{d#Lambda} (GeV/c^{2})");
  // hMass[0]->SetMaximum(hMass[NV-1]->GetMaximum()/0.8);
  hMass[0]->Draw("hist");
  for(int i=0;i<NV;i++) {
    hMass[i]->SetMarkerSize(1.0);
    hMass[i]->SetMarkerColor(color[i]);
    hMass[i]->SetMarkerStyle(kFullCircle);
    hMass_rotate[i]->SetMarkerSize(1.0);
    hMass_rotate[i]->SetMarkerColor(color[i]);
    hMass_rotate[i]->SetMarkerStyle(kOpenCircle);
    hMass[i]->Draw("same");
    hMass_rotate[i]->Draw("same");

  }

  TLegend *leg = new TLegend(0.7, 0.22, 0.94, 0.44);
  leg->AddEntry(hMass[0], "v1=0","lp");  
  leg->AddEntry(hMass[1], "v1#neq0","lp");  
  leg->AddEntry(hMass_rotate[0], "v1=0 rotate","lp");  
  leg->AddEntry(hMass_rotate[1], "v1#neq0 rotate","lp");  
  leg->Draw();
  addpdf(pdf);

  TH1D *h0 = new TH1D("h0","",1,2.98, 3.18);
  h0->SetMinimum(0.5);
  h0->SetMaximum(1.2);
  h0->GetXaxis()->SetTitle("M_{d#Lambda} (GeV/c^{2})");
  h0->GetYaxis()->SetTitle("Rotation/No-rotation");
  h0->Draw();

  for(int i=0;i<NV;i++) {
    hRatio[i]->SetMarkerSize(1.0);
    hRatio[i]->SetMarkerColor(color[i]);
    hRatio[i]->SetMarkerStyle(kOpenCircle);
    hRatio[i]->SetLineColor(color[i]);
    hRatio[i]->Draw("same");
    hRatio[i]->GetYaxis()->SetTitle("Rotation/No-rotation");
  }
    
  drawLine(2.98,1,3.18,1,1.5,2,kBlack);

  TLegend *leg1 = new TLegend(0.6, 0.22, 0.94, 0.34);
  leg1->AddEntry(hRatio[0], "v1=0","pl");
  leg1->AddEntry(hRatio[1], "v1#neq0","pl");
  leg1->Draw();
  addpdf(pdf);

  TH1D* hflowratio = (TH1D*)hMass[1]->Clone("hflowratio");
  hflowratio->Divide(hMass[0]);
  hflowratio->GetXaxis()->SetTitle("M_{d#Lambda} (GeV/c^{2})");
  hflowratio->GetYaxis()->SetTitle("v1#neq0/v1=0");
  hflowratio->Draw();
  hflowratio->GetYaxis()->SetRangeUser(0.9,1.2);
  drawLine(2.98,1,3.18,1,1.5,2,kBlack);
  addpdf(pdf);

  pdf->On();
  pdf->Close();
}
