#include "style.h"
/* #include "sPhenixStyle.h" */
#include "draw.h"
#include "flow.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TStopwatch.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TLine.h"
#include "TNtuple.h"

TF1 *func_phi;
TF1 *func_pT;
const Double_t M_H = 2.99;
const Double_t M_d = 1.875;
const Double_t M_Ld = 1.115;
const Double_t yMax = 1.05; // beam rapidity
const Double_t Teff = 0.2;

TF1 *levyfit4;
TF1 *levyfit5;
TF1 *levyfit6;
TF1 *levyfit7;
TF1 *levyfit8;
TF1 *t_quadr;
TF1* flevyd[10];

void book()
{
  func_phi = new TF1("func_phi","1+2.*[0]*cos(x)",0.,TMath::Pi()*2.);
  func_pT = new TF1("func_pT","x*exp(-(sqrt(x*x+[1]*[1])-[1])/[0])",0.,5.0);

  double dpara[10][3]=
  {{0.241832,-22.8805,4.10711},
    {0.231218,-27.529,4.41142},
    {0.243687,-19.255,3.74636},
    {0.232387,-20.6127,4.07431},
    {0.226336,-19.886,4.26132},
    {0.206177,-24.6636,4.82047},
    {0.192037,-24.6565,5.39291},
    {0.166661,-36.4815,6.59758},
    {0.140038,-74.0464,8.24021},
    {0.120619,-226.527,8.46702}
  };
  for (int i=0;i<10;i++){
    flevyd[i] = new TF1(Form("flevyd_%d",i+3),"2*3.1415*x*[2]*pow((1+(sqrt(1.875*1.875+x*x)-1.875)/[1]/[0]),(-1)*[1])", 0.2,4);
    flevyd[i]->SetParameters(dpara[i]);
  }

  levyfit4 = new TF1("levyfit4","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.4,2.4);
  levyfit5 = new TF1("levyfit5","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.2,2.4);
  levyfit6 = new TF1("levyfit6","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.2,2.4);
  levyfit7 = new TF1("levyfit7","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.2,2.4);
  levyfit8 = new TF1("levyfit8","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.2,2.4);
  t_quadr = new TF1("t_quadr","[0]+[1]*x+[2]*x*x",-2,2);
  levyfit4->SetParameters(0.141849, -1.21226e+08, 1.879e+07);
  levyfit5->SetParameters(0.145031, -1.19088e+08, 1.43524e+07);
  levyfit6->SetParameters(0.141403, -2.20994e+08, 1.2987e+07);
  levyfit7->SetParameters(0.128239, -1.60612e+08, 1.17238e+07);
  levyfit8->SetParameters(0.114106, -1.35632e+08, 8.78861e+06);
  t_quadr->SetParameters(1.15426e+00, 0 ,-8.55891e-01 );
}

void getWeight(double pt, double y, int PDG, double &ptwt, double &rapwt)
{
  ptwt=1;rapwt=1;
  if (PDG==3122)
  {
    if (y<-0.7){ ptwt = levyfit8->Eval(pt); }
    else if (y<-0.5){ ptwt = levyfit7->Eval(pt); }
    else if (y<-0.3){ ptwt = levyfit6->Eval(pt); }
    else if (y<-0.1){ ptwt = levyfit5->Eval(pt); }
    else if (y<0.1) { ptwt = levyfit4->Eval(pt); }
    else if (y<0.3) { ptwt = levyfit5->Eval(pt); }
    else if (y<0.5) { ptwt = levyfit6->Eval(pt); }
    else if (y<0.7) { ptwt = levyfit7->Eval(pt); }
    else { ptwt = levyfit8->Eval(pt); }
    rapwt= t_quadr->Eval(y);   
  } 
  else if (PDG==1000010020)
  {
    rapwt=1;
    double rapbins[11]={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
    int idx=-99;
    for (idx=0;idx<10;idx++)
    {
      if (fabs(y)<rapbins[idx+1]&&fabs(y)>rapbins[idx])
        break;
    }
    if (fabs(y)>1) idx=9;
    ptwt = flevyd[idx]->Eval(pt);
  }
  else {cout <<"input correct PDG id" <<endl; }
}

void sampleKenimatics(double v1, double mass, double Teff, double psi_EP, TLorentzVector& Ld, double pT_min=0, double pT_max=5)
{
  double y_Ld = gRandom->Rndm()*yMax;
  func_pT->SetParameters(Teff, mass);
  double pT_Ld = 0;
  do {
    /* pT_Ld = func_pT->GetRandom(); */
    pT_Ld = gRandom->Rndm()*(pT_max);
  } while (pT_Ld<pT_min || pT_Ld>pT_max);    

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
  TStopwatch time;
  time.Start();
  style();

  const Int_t NM = 20000000;
  const Int_t NRebin = 10;
  const Int_t NV = 2;
  /* const Double_t v1[NV] = {0.0, 0.2, 0.4}; */

  TRandom3 *gRandom = new TRandom3();
  gRandom->SetSeed(0);

  book();

  TH2D *hInvMassPt[NV];
  TH2D *hInvMassPtRotate[NV];
  for(int i=0;i<NV;i++) {
    hInvMassPt[i] = new TH2D(Form("invmass_pt_%d",i),"",500, 2.98, 3.18, 100, 0., 5.0);
    hInvMassPt[i]->Sumw2();
    hInvMassPtRotate[i] = new TH2D(Form("invmass_rt_pt_%d",i),"",500, 2.98, 3.18, 100, 0., 5.0);
    hInvMassPtRotate[i]->Sumw2();
  }
  TNtuple* tree[2]; 
  tree[0] = new TNtuple("dLbTreeNoFlow","dLbTreeNoFlow","ep:dpt:deta:dphi:dv1:dwt:Lbpt:Lbeta:Lbphi:Lbv1:Lbptwt:Lbywt");
  tree[1] = new TNtuple("dLbTreeFlow","dLbTreeFlow","ep:dpt:deta:dphi:dv1:dwt:Lbpt:Lbeta:Lbphi:Lbv1:Lbptwt:Lbywt");

  for(int iv=0;iv<NV;iv++) {
    cout << " Calculate v1 mode " << iv << endl;
    for(int i=0;i<NM;i++) {
      if ((iv==0 && i%100000==0) || (iv==1 && i%10000==0))
        cout << " Processing " << i << "-th event " << endl;

      TLorentzVector deutetron, Ld;
      double psi_EP = gRandom->Rndm()*2*TMath::Pi();
      /* double psi_EP = 0; */

      double y_d = gRandom->Rndm()*yMax;
      double pT_d = 0;
      do {
        /* pT_d = func_pT->GetRandom(); */
        pT_d = gRandom->Rndm()*5;
      } 
      while (pT_d<0.15);
      // func_pT->SetParameters(Teff, M_d);
      // double w_d = func_pT->Eval(pT_d)*(5.0-0.15)*gd_10_40->Eval(y_d); 
      double w_d=1,ptw_d=1,rapw_d=1;
      getWeight(pT_d, y_d, 1000010020, ptw_d, rapw_d);
      w_d = ptw_d*rapw_d * (5-0.15); 
      double mT_d = TMath::Sqrt(pT_d*pT_d+M_d*M_d);
      double pz_d = mT_d*TMath::SinH(y_d);
      double p_d = TMath::Sqrt(pT_d*pT_d+pz_d*pz_d);
      double eta_d = 0.5*TMath::Log((p_d+pz_d)/(p_d-pz_d));
      /* v1_d = gLdv1->Eval(y_d); */
      double v1_d = 0;
      if (iv==1){
        if (y_d>0) v1_d = gdv1->Eval(-1*y_d)*-1;
        if (y_d<0) v1_d = gdv1->Eval(y_d);
        if (v1_d>=0.5) v1_d = 0.5;
        /* v1_d=0.4; */
        func_phi->SetParameter(0, v1_d);  
        if (func_phi->Eval(3.1415)<0) {func_phi->Draw();return;}
      }
      double phi_d=0;
      if (iv==0) phi_d = gRandom->Rndm()*2*TMath::Pi() + psi_EP; 
      if (iv==1) phi_d = func_phi->GetRandom() + psi_EP;

      deutetron.SetPtEtaPhiM(pT_d, eta_d, phi_d, M_d);
      TLorentzVector deutetron_rotate;
      deutetron_rotate.SetPtEtaPhiM(pT_d, eta_d, phi_d+TMath::Pi(), M_d);
      // cout << "d "<< pT_d<<" "<<y_d<<endl; 

      double y_Ld = gRandom->Rndm()*yMax;
      // func_pT->SetParameters(0.1, M_Ld);
      double pT_Ld = 0;
      do {
        /* pT_Ld = func_pT->GetRandom(); */
        pT_Ld = gRandom->Rndm()*5;
      } while (pT_Ld<0.15);    
      // double w_Ld = func_pT->Eval(pT_Ld)*(5.0-0.15); 
      double w_Ld=1,ptwt_Ld=1,ywt_Ld=1;
      getWeight(pT_Ld,y_Ld,3122,ptwt_Ld,ywt_Ld);
      // if (y_Ld<-0.7){ ptwt_Ld = levyfit8->Eval(pT_Ld); }
      // else if (y_Ld<-0.5){ ptwt_Ld = levyfit7->Eval(pT_Ld); }
      // else if (y_Ld<-0.3){ ptwt_Ld = levyfit6->Eval(pT_Ld); }
      // else if (y_Ld<-0.1){ ptwt_Ld = levyfit5->Eval(pT_Ld); }
      // else if (y_Ld<0.1) { ptwt_Ld = levyfit4->Eval(pT_Ld); }
      // else if (y_Ld<0.3) { ptwt_Ld = levyfit5->Eval(pT_Ld); }
      // else if (y_Ld<0.5) { ptwt_Ld = levyfit6->Eval(pT_Ld); }
      // else if (y_Ld<0.7) { ptwt_Ld = levyfit7->Eval(pT_Ld); }
      // else { ptwt_Ld = levyfit8->Eval(pT_Ld); }
      // ywt_Ld=t_quadr->Eval(y_Ld);
      w_Ld=ptwt_Ld*ywt_Ld*(5-0.15);

      double v1_Ld=0;
      if (iv==1){
        if (y_Ld<0) v1_Ld = gLdv1->Eval(y_Ld);
        if (y_Ld>0) v1_Ld = gLdv1->Eval(-1*y_Ld)*-1;
        if (v1_Ld>=0.5) v1_Ld = 0.5;
        /* v1_Ld = 0.4; */
        func_phi->SetParameter(0, v1_Ld);  
        if (func_phi->Eval(3.14159265)<0) {func_phi->Draw();cout<<v1_Ld<<endl;return;}
      }

      /* func_phi->Draw(); */
      /* cout<<y_Ld <<v1_Ld << endl; */
      /* return; */
      double phi_Ld=0;
      if (iv==0) phi_Ld = gRandom->Rndm()*2*TMath::Pi() + psi_EP; 
      if (iv==1) phi_Ld = func_phi->GetRandom() + psi_EP;
      double mT_Ld = TMath::Sqrt(pT_Ld*pT_Ld+M_Ld*M_Ld);
      double pz_Ld = mT_Ld*TMath::SinH(y_Ld);
      double p_Ld = TMath::Sqrt(pT_Ld*pT_Ld+pz_Ld*pz_Ld);
      double eta_Ld = 0.5*TMath::Log((p_Ld+pz_Ld)/(p_Ld-pz_Ld));
      Ld.SetPtEtaPhiM(pT_Ld, eta_Ld, phi_Ld, M_Ld); 

      TLorentzVector H3L = deutetron + Ld;
      TLorentzVector H3L_rotate = deutetron_rotate + Ld;
      double weight = w_d*w_Ld;
      hInvMassPt[iv]->Fill(H3L.M(), H3L.Pt(), weight);
      /* hInvMassPt[0]->Fill(H3L.M(), H3L.Pt(), weight); */
      /* hInvMassPt[1]->Fill(H3L.M(), H3L.Pt(), weight*(1+2*cos(phi_Ld)*(1+2*cos(phi_d)))); */
      /* hInvMassPt[1]->Fill(H3L.M(), H3L.Pt(), weight); */
      /* hInvMassPtNoFlowRotate[iv]->Fill(H3L.M(), H3L.Pt()); */
      hInvMassPtRotate[iv]->Fill(H3L_rotate.M(), H3L_rotate.Pt(),weight);
      /* hInvMassPtRotate[0]->Fill(H3L_rotate.M(), H3L_rotate.Pt(),weight); */
      /* hInvMassPtRotate[1]->Fill(H3L_rotate.M(), H3L_rotate.Pt(),weight); */

      //fill tree     
      tree[iv]->Fill(psi_EP,pT_d,eta_d,phi_d,v1_d,w_d,pT_Ld,eta_Ld,phi_Ld,v1_Ld,ptwt_Ld,ywt_Ld);
      // cout << psi_EP<<" "<<pT_d<<" "<<eta_d<<" "<<phi_d<<" "<<v1_d<<" "<<w_d<<" "<<pT_Ld<<" "<<eta_Ld<<" "<<phi_Ld<<" "<<v1_Ld<<" "<<ptwt_Ld<<" "<<ywt_Ld<<endl;
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
  /* hMass[0]->GetXaxis()->SetRangeUser( 2.95, 3.05); */
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
  TH1D *h0 = new TH1D("h0","",1,2.98, 3.18);
  h0->SetMinimum(0);
  h0->SetMaximum(2.0);
  h0->Draw();

  for(int i=0;i<NV;i++) {
    hRatio[i]->SetMarkerSize(1.0);
    hRatio[i]->SetMarkerColor(color[i]);
    hRatio[i]->SetLineColor(color[i]);
    hRatio[i]->Draw("same");
  }
  drawLine(2.98,1,3.18,1,1.5,2);

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
    tree[i]->Write();
  }
  fout->Close();
  time.Stop();
  time.Print();
}
