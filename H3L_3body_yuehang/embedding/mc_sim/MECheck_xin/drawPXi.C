#include "style.C+"
#include "draw.C+"
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
  
  style();

  const Int_t NV = 3;
  const Double_t v1[NV] = {0.0, 0.2, 0.4};

  TH2D *hInvMassPt[NV];
  TH1D *hMass[NV];
  TH1D *hRatio[NV];
  TFile *fin = new TFile("invmass_pXi.root");
  for(int i=0;i<NV;i++) {
    hInvMassPt[i] = (TH2D *)fin->Get(Form("invmass_pt_%d",i));
    hMass[i] = (TH1D *)fin->Get(Form("mass_%d",i));
    hRatio[i] = (TH1D *)fin->Get(Form("ratio_%d",i));    
  }
  
  TCanvas *c1 = new TCanvas("c1", "c1", 0, 0, 600, 800);
  c1->Divide(1,2);
  c1->Draw();

  c1->cd(1);
  hMass[0]->GetXaxis()->SetTitle("M_{p#Xi} (GeV/c^{2})");
  hMass[0]->SetMaximum(hMass[NV-1]->GetMaximum()/0.8);
  hMass[0]->GetXaxis()->SetRange(9,40);
  hMass[0]->Draw("hist");
  for(int i=1;i<NV;i++) {
    hMass[i]->SetMarkerSize(1.0);
    hMass[i]->SetMarkerColor(i);
    hMass[i]->Draw("same");
  }

  TLegend *leg = new TLegend(0.7, 0.22, 0.94, 0.44);
  for(int i=NV-1;i>=1;i--) {
    leg->AddEntry(hMass[i], Form("v_{1} = %3.1f", v1[i]),"p");
  }
  leg->AddEntry(hMass[0], Form("v_{1} = %3.1f", v1[0]),"l");  
  leg->Draw();

  c1->cd(2);
  TH1D *h0 = new TH1D("h0","",1,2.24, 2.4);
  h0->SetMinimum(0.5);
  h0->SetMaximum(1.6);
  h0->GetXaxis()->SetTitle("M_{p#Xi} (GeV/c^{2})");
  h0->GetYaxis()->SetTitle("Ratio");
  h0->Draw();

  for(int i=1;i<NV;i++) {
    hRatio[i]->SetMarkerSize(1.0);
    hRatio[i]->SetMarkerColor(i);
    hRatio[i]->SetLineColor(i);
    hRatio[i]->Draw("same");
  }
    
  TLegend *leg1 = new TLegend(0.6, 0.22, 0.94, 0.34);
  for(int i=NV-1;i>=1;i--) {
    leg1->AddEntry(hRatio[i], Form("v_{1} = %3.1f / v_{1} = 0", v1[i]),"p");
  }
  leg1->Draw();

  c1->Update();
  
  c1->SaveAs("invMass_pXi.pdf");
  c1->SaveAs("invMass_pXi.png");

  
}
