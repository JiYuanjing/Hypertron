#include "style.h"
#include "draw.h"
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

void pXi()
{
  /* gSystem->Load("style_C.so"); */
  /* gSystem->Load("draw_C.so"); */
  
  style();
  const Double_t M_H = 2.230;
  const Double_t M_P = 0.938272;
  const Double_t M_Xi = 1.32171;
  const Double_t yMax = 1.05; // beam rapidity
  const Double_t Teff = 0.2;
  const Int_t NM = 100000;
  const Int_t NRebin = 5;
  const Int_t NV = 3;
  const Double_t v1[NV] = {0.0, 0.2, 0.4};
  
  TRandom3 *gRandom = new TRandom3();

  TF1 *func_phi = new TF1("func_phi","1+2.*[0]*cos(x)",0.,TMath::Pi()*2.);

  TF1 *func_pT = new TF1("func_pT","x*exp(-(sqrt(x*x+[1]*[1])-[1])/[0])",0.,5.0);
  func_pT->SetParameters(0.1, M_P);

  TH2D *hInvMassPt[NV];
  for(int i=0;i<NV;i++) {
    hInvMassPt[i] = new TH2D(Form("invmass_pt_%d",i),"",1000, 2.2, 3.2, 100, 0., 5.0);
    hInvMassPt[i]->Sumw2();
  }

  for(int iv=0;iv<NV;iv++) {
    cout << " Calculate invmass for v1 = " << v1[iv] << endl;
    func_phi->SetParameter(0, v1[iv]);  
    for(int i=0;i<NM;i++) {
      if(i%100==0)
	cout << " Processing " << i << "-th event " << endl;
      
      TLorentzVector proton, xi;
      
      double psi_EP = gRandom->Rndm()*2*TMath::Pi();
      double y_p = gRandom->Rndm()*yMax;
      //    double y_p = 0;
      func_pT->SetParameters(Teff, M_P);
      double pT_p = 0;
      do {
	pT_p = func_pT->GetRandom();
      } while (pT_p<0.15);
      double phi_p = func_phi->GetRandom() + psi_EP;
      double mT_p = TMath::Sqrt(pT_p*pT_p+M_P*M_P);
      double pz_p = mT_p*TMath::SinH(y_p);
      double p_p = TMath::Sqrt(pT_p*pT_p+pz_p*pz_p);
      double eta_p = 0.5*TMath::Log((p_p+pz_p)/(p_p-pz_p));
      proton.SetPtEtaPhiM(pT_p, eta_p, phi_p, M_P);
      
      double y_xi = gRandom->Rndm()*yMax;
      //    double y_xi = 0;
      func_pT->SetParameters(0.1, M_Xi);
      double pT_xi = 0;
      do {
	pT_xi = func_pT->GetRandom();
      } while (pT_xi<0.15);    
      double phi_xi = func_phi->GetRandom() + psi_EP;
      double mT_xi = TMath::Sqrt(pT_xi*pT_xi+M_Xi*M_Xi);
      double pz_xi = mT_xi*TMath::SinH(y_xi);
      double p_xi = TMath::Sqrt(pT_xi*pT_xi+pz_xi*pz_xi);
      double eta_xi = 0.5*TMath::Log((p_xi+pz_xi)/(p_xi-pz_xi));
      xi.SetPtEtaPhiM(pT_xi, eta_xi, phi_xi, M_Xi);
      
      TLorentzVector H = proton + xi;
      hInvMassPt[iv]->Fill(H.M(), H.Pt());
      
      
    } // end int i - number of events
  }

  TH1D *hMass[NV];
  TH1D *hRatio[NV];
  for(int iv=0;iv<NV;iv++) {
    hInvMassPt[iv]->RebinX(NRebin);
    hMass[iv] = (TH1D *)hInvMassPt[iv]->ProjectionX(Form("mass_%d",iv));
    hRatio[iv] = (TH1D *)hMass[iv]->Clone(Form("ratio_%d",iv));
    hRatio[iv]->Divide(hMass[0]);
  }
  
  TCanvas *c1 = new TCanvas("c1", "c1", 0, 0, 600, 800);
  c1->Divide(1,2);
  c1->Draw();

  c1->cd(1);
  hMass[0]->SetMaximum(hMass[NV-1]->GetMaximum()/0.8);
  hMass[0]->GetXaxis()->SetRange(9,40);
  hMass[0]->Draw("hist");
  for(int i=1;i<NV;i++) {
    hMass[i]->SetMarkerSize(1.0);
    hMass[i]->SetMarkerColor(i);
    hMass[i]->Draw("same");
  }

  c1->cd(2);
  TH1D *h0 = new TH1D("h0","",1,2.24, 2.4);
  h0->SetMinimum(0);
  h0->SetMaximum(2.0);
  h0->Draw();

  for(int i=1;i<NV;i++) {
    hRatio[i]->SetMarkerSize(1.0);
    hRatio[i]->SetMarkerColor(i);
    hRatio[i]->SetLineColor(i);
    hRatio[i]->Draw("same");
  }
    
  c1->Update();
  
  c1->SaveAs("invMass_pXi.pdf");
  c1->SaveAs("invMass_pXi.png");

  TFile *fout = new TFile("invmass_pXi.root","recreate");
  for(int i=0;i<NV;i++) {
    hInvMassPt[i]->Write();
    hMass[i]->Write();
    hRatio[i]->Write();
  }
  fout->Close();
  
  
  
}
