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
  
  TGraph* halandrad = new TGraph("corr.txt"); 
  halandrad->Draw();

  style();
  // const Double_t M_H = 2.230;
  const Double_t M_H = 2.9913;
  const Double_t M_P = 1.8756;
  // const Double_t M_Xi = 1.32171;
  const Double_t M_Xi = 1.1158;
  const Double_t yMax = 1.05; // beam rapidity
  const Double_t Teff = 0.2;
  const Int_t NM = 100000;
  // const Int_t NM = 1000;
  const Int_t NRebin = 5;
  const Int_t NRebin2 = 10;
  const Int_t NV = 3;
//  const Double_t v1[NV] = {0.0, 0.2, 0.4};
//  v1;

  TRandom3 *gRandom = new TRandom3();

  TF1 *func_phi = new TF1("func_phi","fabs(1+2.*[0]*cos(x))",0.,TMath::Pi()*2.);
  TF1 *func_pT = new TF1("func_pT","x*exp(-(sqrt(x*x+[1]*[1])-[1])/[0])",0.,5.0);
  func_pT->SetParameters(0.1, M_P);

  TH2D *hInvMassPt[NV];
  for(int i=0;i<NV;i++) {
    hInvMassPt[i] = new TH2D(Form("invmass_pt_%d",i),"",1000, 2.95, 3.2, 100, 0., 5.0);
    hInvMassPt[i]->Sumw2();
  }

  TH1D *hQ[NV];
  for(int i=0;i<NV;i++) {
    hQ[i] = new TH1D(Form("rawcorr_%d",i),"",400, 0,0.2);
    hQ[i]->Sumw2();
  }

  for(int iv=0;iv<NV;iv++) {
    //cout << " Calculate invmass for v1 = " << v1[iv] << endl;
    //func_phi->SetParameter(0, v1[iv]);
    for(int i=0;i<NM;i++) {
      if (i%100==0)
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

      if(iv==0){
	       func_phi->SetParameter(0, 0.0);
       } else if(iv==1 || iv==2){
         func_phi->SetParameter(0, (pT_p-0.15)*0.3);
      }

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

      // func_phi->SetParameter(0, (pT_xi-0.15)*0.3);

	     if(iv==0){
        func_phi->SetParameter(0, 0.0);
       } 
       else if(iv==1 || iv==2){
	      func_phi->SetParameter(0, (pT_xi-0.15)*0.3);
       }

      double phi_xi = func_phi->GetRandom() + psi_EP;
      double mT_xi = TMath::Sqrt(pT_xi*pT_xi+M_Xi*M_Xi);
      double pz_xi = mT_xi*TMath::SinH(y_xi);
      double p_xi = TMath::Sqrt(pT_xi*pT_xi+pz_xi*pz_xi);
      double eta_xi = 0.5*TMath::Log((p_xi+pz_xi)/(p_xi-pz_xi));
      xi.SetPtEtaPhiM(pT_xi, eta_xi, phi_xi, M_Xi);

      TLorentzVector H = proton + xi;
      TLorentzVector Qvect = (xi-proton);

      double corr_weight;
      double Pinv = H.Mag();
      double Q1 = (M_Xi*M_Xi-M_P*M_P)/Pinv;
      double Q=sqrt(Q1*Q1-Qvect.Mag2());
      double kstar = Q/2.0;

      if (kstar<0.1) corr_weight=halandrad->Eval(kstar*1e3);
      else if (kstar<0.003) corr_weight = halandrad->Eval(3);
      else corr_weight=1.;

      if(iv==0||iv==1){
        hInvMassPt[iv]->Fill(H.M(), H.Pt());
        hQ[iv]->Fill(kstar);
      }else{
        hInvMassPt[iv]->Fill(H.M(), H.Pt(),corr_weight);
        hQ[iv]->Fill(kstar,corr_weight);
      }

    } // end int i - number of events
  }

  TH1D *hMass[NV];
  TH1D *hRatio[NV];
  TH1D *hRatioQ[NV];

  for(int iv=0;iv<NV;iv++) {
    hInvMassPt[iv]->RebinX(NRebin);
    hMass[iv] = (TH1D *)hInvMassPt[iv]->ProjectionX(Form("mass_%d",iv));
    hRatio[iv] = (TH1D *)hMass[iv]->Clone(Form("ratio_%d",iv));
    hRatio[iv]->Divide(hMass[0]);

    hQ[iv]->Rebin(NRebin2);
    hRatioQ[iv] = (TH1D *)hQ[iv]->Clone(Form("ratio_q_%d",iv));
    hRatioQ[iv]->Divide(hQ[0]);
  }

  TCanvas *c1 = new TCanvas("c1", "c1", 0, 0, 600, 800);
  c1->Divide(1,2);
  c1->Draw();

  c1->cd(1);
  hMass[0]->SetMaximum(hMass[NV-1]->GetMaximum()/0.8);
  hMass[0]->GetXaxis()->SetRangeUser(2.24, 2.4);
  hMass[0]->GetXaxis()->SetTitle("m_{p#Xi}[GeV/c^{2}]");
  hMass[0]->Draw("hist");
  for(int i=1;i<NV;i++) {
    hMass[i]->SetMarkerSize(1.0);
    hMass[i]->SetMarkerColor(i);
    hMass[i]->Draw("same");
  }

  c1->cd(2);
  TH1D *h0 = new TH1D("h0","",1,2.95, 3.2);
  h0->SetMinimum(0);
  h0->SetMaximum(2.0);
  h0->GetXaxis()->SetTitle("m_{p#Xi}[GeV/c^{2}]");
  h0->GetYaxis()->SetTitle("Ratio");
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


  TCanvas *c2 = new TCanvas("c2", "c2", 0, 0, 600, 800);
  c2->cd();

  c2->Divide(1,2);
  c2->Draw();

  c2->cd(1);
  hQ[0]->SetMaximum(hQ[NV-1]->GetMaximum()/0.8);
  hQ[0]->GetXaxis()->SetRangeUser(0,1);
  hQ[0]->GetXaxis()->SetTitle("k^{*}[GeV/c]");
  hQ[0]->Draw("hist");
  for(int i=1;i<NV;i++) {
    hQ[i]->SetMarkerSize(1.0);
    hQ[i]->SetMarkerColor(i);
    hQ[i]->Draw("same");
  }

  c2->cd(2);
  TH1D *h1 = new TH1D("h1","",1000,0,0.2);
  h1->SetMinimum(0);
  h1->SetMaximum(3.0);
  h1->GetXaxis()->SetTitle("k^{*}[GeV/c]");
  h1->GetYaxis()->SetTitle("Ratio");
  h1->Draw();

  for(int i=1;i<NV;i++) {
    hRatioQ[i]->SetMarkerSize(1.0);
    hRatioQ[i]->SetMarkerColor(i);
    hRatioQ[i]->SetLineColor(i);
    hRatioQ[i]->Draw("same");
  }

  c2->Update();

  c2->SaveAs("corrQ_pXi.pdf");
  c2->SaveAs("corrQ_pXi.png");

  // hQ[0]->Draw("hist");
  // for(int i=1;i<NV;i++) {
  //   hQ[i]->SetMarkerSize(1.0);
  //   hQ[i]->SetMarkerColor(i);
  // hQ[i]->Draw("same");
  // }
  // c2->SaveAs("rawcorr_pXi.pdf");

  TFile *fout = new TFile("invmass_pXi.root","recreate");
  for(int i=0;i<NV;i++) {
    hInvMassPt[i]->Write();
    hMass[i]->Write();
    hRatio[i]->Write();
    hQ[i]->Write();
    hRatioQ[i]->Write();
  }
  fout->Close();



}
