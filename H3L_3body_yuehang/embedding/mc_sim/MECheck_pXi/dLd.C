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
#include "flow.h"

TF1 *func_phi = new TF1("func_phi","1+2.*[0]*cos(x)",0.,TMath::Pi()*2.);
TF1 *func_pT = new TF1("func_pT","x*exp(-(sqrt(x*x+[1]*[1])-[1])/[0])",0.,5.0);
const Double_t M_H = 2.23;
const Double_t M_d = 0.938;
const Double_t M_Ld = 1.32;
const Double_t yMax = 1.05; // beam rapidity
const Double_t Teff = 0.2;
void sampleKenimatics(double v1, double mass, double Teff, double psi_EP, TLorentzVector& Ld, double pT_min=0, double pT_max=5)
{
  double y_Ld = gRandom->Rndm()*yMax;
  func_pT->SetParameters(Teff, mass);
  double pT_Ld = 0;
  do {
    pT_Ld = func_pT->GetRandom();
  } while (pT_Ld<pT_max && pT_Ld>pT_min);    

  func_phi->SetParameter(0, v1);
  double phi_Ld = func_phi->GetRandom() + psi_EP;
  double mT_Ld = TMath::Sqrt(pT_Ld*pT_Ld+M_Ld*M_Ld);
  double pz_Ld = mT_Ld*TMath::SinH(y_Ld);
  double p_Ld = TMath::Sqrt(pT_Ld*pT_Ld+pz_Ld*pz_Ld);
  double eta_Ld = 0.5*TMath::Log((p_Ld+pz_Ld)/(p_Ld-pz_Ld));
  Ld.SetPtEtaPhiM(pT_Ld, eta_Ld, phi_Ld, M_Ld); 
}
void dLd()
{
  style();

  const Int_t NM = 100000;
  const Int_t NRebin = 10;
  const Int_t NV = 2;
  /* const Double_t v1[NV] = {0.0, 0.2, 0.4}; */

  TRandom3 *gRandom = new TRandom3();
  gRandom->SetSeed(0);

  TH2D *hInvMassPt[NV];
  TH2D *hInvMassPtRotate[NV];
  TH2D* hPhiRap[NV];
  for(int i=0;i<NV;i++) {
    hInvMassPt[i] = new TH2D(Form("invmass_pt_%d",i),"",1000, 2.2, 3.2,  100, 0., 5.0);
    hInvMassPt[i]->Sumw2();
    hInvMassPtRotate[i] = new TH2D(Form("invmass_rt_pt_%d",i),"",1000, 2.2, 3.2, 100, 0., 5.0);
    hInvMassPtRotate[i]->Sumw2();
    hPhiRap[i] = new TH2D(Form("Phi_%d",i),"",100,0,TMath::Pi()*2.0, 20, 0., yMax);
  }

  for(int iv=0;iv<NV;iv++) {
    cout << " Calculate v1 mode " << iv << endl;
    for(int i=0;i<NM;i++) {
      if ((iv==0 && i%1000==0) || (iv==1 && i%10000==0))
        cout << " Processing " << i << "-th event " << endl;

      TLorentzVector deutetron, Ld;
      double psi_EP = gRandom->Rndm()*2*TMath::Pi();
      /* double psi_EP = 0; */

      double y_d = gRandom->Rndm()*yMax;
      double pT_d = 0;
      func_pT->SetParameters(0.2, M_d);
      do {
        /* pT_d = func_pT->GetRandom(); */
           pT_d = gRandom->Rndm()*5;
      } 
      while (pT_d<0.15);
      double w_d = func_pT->Eval(pT_d)*(5-0.15); 

      double mT_d = TMath::Sqrt(pT_d*pT_d+M_d*M_d);
      double pz_d = mT_d*TMath::SinH(y_d);
      double p_d = TMath::Sqrt(pT_d*pT_d+pz_d*pz_d);
      double eta_d = 0.5*TMath::Log((p_d+pz_d)/(p_d-pz_d));
      /* v1_d = gLdv1->Eval(y_d); */
      double v1_d = 0;
      if (iv==1){
        /* if (y_d>0) v1_d = gdv1->Eval(-1*y_d)*-1; */
        /* if (y_d<0) v1_d = gdv1->Eval(y_d); */
        /* if (v1_d>=0.5) v1_d = 0.5; */
        v1_d = 0.5;
        func_phi->SetParameter(0, v1_d);  
        if (func_phi->Eval(3.1415)<0) {func_phi->Draw();return;}
      }
      double phi_d=0;
      if (iv==1) phi_d = func_phi->GetRandom() + psi_EP;
      if (iv==0) phi_d = gRandom->Rndm()*2*TMath::Pi() + psi_EP; 
      deutetron.SetPtEtaPhiM(pT_d, eta_d, phi_d, M_d);
      TLorentzVector deutetron_rotate;
      deutetron_rotate.SetPtEtaPhiM(pT_d, eta_d, phi_d+TMath::Pi(), M_d);

      double y_Ld = gRandom->Rndm()*yMax;
      double pT_Ld = 0;
      func_pT->SetParameters(0.2, M_Ld);
      do {
        /* pT_Ld = func_pT->GetRandom(); */
       pT_Ld = gRandom->Rndm()*5; 
      } while (pT_Ld<0.15);    

      double w_Ld = func_pT->Eval(pT_Ld)*(5-0.15); 

      double v1_Ld=0;
      if (iv==1){
        /* if (y_Ld<0) v1_Ld = gLdv1->Eval(y_Ld); */
        /* if (y_Ld>0) v1_Ld = gLdv1->Eval(-1*y_Ld)*-1; */
        /* if (v1_Ld>=0.5) v1_Ld = 0.4999; */
        v1_Ld = 0.5;
        func_phi->SetParameter(0, v1_Ld);  
        if (func_phi->Eval(3.14159265)<0) {func_phi->Draw();cout<<v1_Ld<<endl;return;}
      }

      /* func_phi->Draw(); */
      /* cout<<y_Ld <<v1_Ld << endl; */
      /* return; */
      double phi_Ld=0;
      if (iv==1) phi_Ld = func_phi->GetRandom() + psi_EP;
      if (iv==0) phi_Ld = gRandom->Rndm()*2*TMath::Pi() + psi_EP; 
      double mT_Ld = TMath::Sqrt(pT_Ld*pT_Ld+M_Ld*M_Ld);
      double pz_Ld = mT_Ld*TMath::SinH(y_Ld);
      double p_Ld = TMath::Sqrt(pT_Ld*pT_Ld+pz_Ld*pz_Ld);
      double eta_Ld = 0.5*TMath::Log((p_Ld+pz_Ld)/(p_Ld-pz_Ld));
      Ld.SetPtEtaPhiM(pT_Ld, eta_Ld, phi_Ld, M_Ld); 

      TLorentzVector H3L = deutetron + Ld;
      TLorentzVector H3L_rotate = deutetron_rotate + Ld;
      double weight = w_d*w_Ld;
      /* double weight = 1; */
      
      hInvMassPt[iv]->Fill(H3L.M(), H3L.Pt(), weight);
      /* hInvMassPt[0]->Fill(H3L.M(), H3L.Pt(), weight); */
      /* hInvMassPt[1]->Fill(H3L.M(), H3L.Pt(), weight*(1+2*cos(phi_Ld)*(1+2*cos(phi_d)))); */
      /* hInvMassPt[1]->Fill(H3L.M(), H3L.Pt(), weight); */
      /* hInvMassPtNoFlowRotate[iv]->Fill(H3L.M(), H3L.Pt()); */
      hInvMassPtRotate[iv]->Fill(H3L_rotate.M(), H3L_rotate.Pt(),weight);
      /* hInvMassPtRotate[0]->Fill(H3L_rotate.M(), H3L_rotate.Pt(),weight); */
      /* hInvMassPtRotate[1]->Fill(H3L_rotate.M(), H3L_rotate.Pt(),weight); */

    } // end int i - number of events
  }

  TH1D *hMass[NV];
  TH1D *hMass_rotate[NV];
  TH1D *hRatio[NV];
  int color[2]={kRed,kBlue};
  for(int iv=0;iv<NV;iv++) {
    hInvMassPt[iv]->RebinX(NRebin);
    hInvMassPtRotate[iv]->RebinX(NRebin);
    hMass[iv] = (TH1D *)hInvMassPt[iv]->ProjectionX(Form("mass_%d",iv));
    hMass_rotate[iv] = (TH1D *)hInvMassPtRotate[iv]->ProjectionX(Form("mass_rotate_%d",iv));
    hRatio[iv] = (TH1D *)hMass_rotate[iv]->Clone(Form("ratio_%d",iv));
    hRatio[iv]->Divide(hMass[iv]);
  }

  TCanvas *c1 = new TCanvas("c1", "c1", 0, 0, 600, 800);
  c1->Divide(1,2);
  c1->Draw();

  c1->cd(1);
  /* hMass[0]->SetMaximum(hMass[NV-1]->GetMaximum()/0.8); */
  hMass[0]->GetXaxis()->SetRangeUser( 2.24, 2.4);
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

  c1->cd(2);
  TH1D *h0 = new TH1D("h0","",1,2.24, 2.4);
  h0->SetMinimum(0);
  h0->SetMaximum(2.0);
  h0->Draw();

  for(int i=1;i<NV;i++) {
    hRatio[i]->SetMarkerSize(1.0);
    hRatio[i]->SetMarkerColor(color[i]);
    hRatio[i]->SetMarkerStyle(kOpenCircle);
    hRatio[i]->SetLineColor(color[i]);
    hRatio[i]->Draw("same");
  }

  c1->Update();

  c1->SaveAs("invMass_dLd.pdf");
  c1->SaveAs("invMass_dLd.png");

  TFile *fout = new TFile("invmass_dLd.root","recreate");
  for(int i=0;i<NV;i++) {
    hInvMassPt[i]->Write();
    hInvMassPtRotate[i]->Write();
    hMass[i]->Write();
    hMass_rotate[i]->Write();
    hRatio[i]->Write();
  }
  fout->Close();
}
