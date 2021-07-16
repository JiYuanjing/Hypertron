#include "sPhenixStyle.h"
void Norm(TH1F* h)
{
  h->Scale(1./h->Integral()*h->GetBinWidth(1));
  h->Scale(1./h->GetMaximum());
  h->SetDirectory(0);
  h->GetYaxis()->SetTitle("Arb. Unit");
}

/* void projAndComp(TH2F* h2H3L, TH2F* h2La,TCanvas* c,TPDF* pdf ) */
void projAndComp(TString name, TFile* fH3L, TFile* fLa, TCanvas* c,TPDF* pdf,TString xTitle , TString drawstyle="", int rebin=1,  TString legtitle1="H3L quasi", TString legtitle2="#Lambda+d", TString text="dp#pi")
{
  fH3L->cd();
  TH2F* h2H3L = (TH2F*)fH3L->Get(name.Data());
  h2H3L->Sumw2();
  h2H3L->SetDirectory(0);
  fLa->cd();
  TH2F* h2La = (TH2F*)fLa->Get(name.Data());
  h2La->Sumw2();
  h2La->SetDirectory(0);
  c->Clear();
  /* c->Divide(3,2); */
  /* int  const bins=6; */
  /* double ptedge[bins+1]={0,0.5,1,2,3,4,5}; */
  c->Divide(2,2);
  int  const bins=4;
  double ptedge[bins+1]={0,0.5,1,2,3};

  for (int i=0;i<bins;i++)
  {
    c->cd(i+1);
    TH1F* h1H3= (TH1F*)h2H3L->ProjectionY(Form("h1H3%d",i), h2H3L->GetXaxis()->FindBin(ptedge[i]),h2H3L->GetXaxis()->FindBin(ptedge[i+1]));
    TH1F* h1La= (TH1F*)h2La->ProjectionY(Form("h1La%d",i), h2La->GetXaxis()->FindBin(ptedge[i]),h2La->GetXaxis()->FindBin(ptedge[i+1]));
    h1La->SetLineColor(kRed);
    h1La->SetMarkerColor(kRed);
    h1La->SetMarkerStyle(kOpenCircle);
    h1La->SetMarkerSize(1);
    h1H3->SetLineColor(kBlue);
    h1H3->SetMarkerColor(kBlue);
    h1H3->SetMarkerSize(1);
    if (rebin>1) {
      h1H3->Rebin(rebin);
      h1La->Rebin(rebin);
    }
    Norm(h1La);
    Norm(h1H3);
    h1H3->GetXaxis()->SetTitle(xTitle.Data());
    h1H3->GetYaxis()->SetRangeUser(-0.1 , 1.15);
    /* h1H3->DrawCopy(drawstyle.Data()); */
    h1H3->DrawCopy("");
    h1La->DrawCopy((drawstyle+"same").Data()); 
    // TRatioPlot* rp = new TRatioPlot(h1La, h1H3);
    // gPad->SetTicks(0, 1);    rp->Draw();
    drawLatex( 0.18, 0.88, Form("%0.1f<(%s) p_{T}<%0.1f GeV/c", ptedge[i],text.Data(),ptedge[i+1]),0.055);
    TLegend* l = new TLegend(0.65,0.68,0.9,0.87);
    l->AddEntry( h1La, legtitle2.Data(), "pl");
    l->AddEntry(h1H3,legtitle1.Data(),"pl");
    l->Draw();
  }
  c->cd();
  addpdf(pdf);
}

void projAndCompH3L(TString name, TFile* fH3L, TFile* fquasi, TCanvas* c,TPDF* pdf,TString xTitle , TString drawstyle="", int rebin=1)
{
  fH3L->cd();
  TH2F* h2H3L = (TH2F*)fH3L->Get(name.Data());
  h2H3L->Sumw2();
  h2H3L->SetDirectory(0);
  fquasi->cd();
  TH2F* h2La = (TH2F*)fquasi->Get(name.Data());
  h2La->Sumw2();
  h2La->SetDirectory(0);
  c->Clear();
  c->Divide(2,2);
  /* int  const bins=6; */
  /* double ptedge[bins+1]={0,0.2,0.4,1,2,3,4}; */
  int  const bins=4;
  double ptedge[bins+1]={0,0.5,1,2,3};

  for (int i=0;i<bins;i++)
  {
    c->cd(i+1);
    TH1F* h1H3= (TH1F*)h2H3L->ProjectionY("h1H3", h2H3L->GetXaxis()->FindBin(ptedge[i]),h2H3L->GetXaxis()->FindBin(ptedge[i+1]));
    TH1F* h1La= (TH1F*)h2La->ProjectionY("h1La", h2La->GetXaxis()->FindBin(ptedge[i]),h2La->GetXaxis()->FindBin(ptedge[i+1]));
    h1La->SetLineColor(kRed);
    h1La->SetMarkerColor(kRed);
    h1La->SetMarkerStyle(kOpenCircle);
    h1La->SetMarkerSize(1);
    h1H3->SetLineColor(kBlue);
    h1H3->SetMarkerColor(kBlue);
    h1H3->SetMarkerSize(1);
    if (rebin>1) {
      h1H3->Rebin(rebin);
      h1La->Rebin(rebin);
    }
    Norm(h1La);
    Norm(h1H3);
    h1H3->GetXaxis()->SetTitle(xTitle.Data());
    h1H3->GetYaxis()->SetRangeUser(0 , h1H3->GetMaximum()*1.3);
    /* h1H3->DrawCopy(drawstyle.Data()); */
    h1H3->DrawCopy("");
    h1La->DrawCopy((drawstyle+"same").Data()); 
    drawLatex( 0.18, 0.88, Form("%0.1f<(p#pi) p_{T}<%0.1f GeV/c", ptedge[i],ptedge[i+1]),0.055);
    TLegend* l = new TLegend(0.65,0.68,0.9,0.87);
    l->AddEntry( h1H3, "Phase", "pl");
    l->AddEntry(h1La,"Quasi","pl");
    l->Draw();
  }
  c->cd();
  addpdf(pdf);
}

void projAndCompLaDeu(TString name, TFile* fH3L, TFile* fLa, TCanvas* c,TPDF* pdf,TString xTitle , TString drawstyle="", int rebin=1,  TString legtitle1="H3L", TString legtitle2="#Lambda")
{
  fH3L->cd();
  TH2F* h2H3L = (TH2F*)fH3L->Get(name.Data());
  h2H3L->Sumw2();
  h2H3L->SetDirectory(0);
  fLa->cd();
  TH2F* h2La = (TH2F*)fLa->Get(name.Data());
  h2La->Sumw2();
  h2La->SetDirectory(0);
  c->Clear();
  c->Divide(1,1);
  int  const bins=1;
  /* double ptedge[bins+1]={0,0.5,1,2,3,4,5}; */
  double ptedge[bins+1]={0,5};
  /* c->Divide(2,2); */
  /* int  const bins=4; */
  /* double ptedge[bins+1]={0,0.5,1,2,3}; */

  for (int i=0;i<bins;i++)
  {
    c->cd(i+1);
    TH1F* h1H3= (TH1F*)h2H3L->ProjectionY("h1H3", h2H3L->GetXaxis()->FindBin(ptedge[i]),h2H3L->GetXaxis()->FindBin(ptedge[i+1]));
    TH1F* h1La= (TH1F*)h2La->ProjectionY("h1La", h2La->GetXaxis()->FindBin(ptedge[i]),h2La->GetXaxis()->FindBin(ptedge[i+1]));
    h1La->SetLineColor(kRed);
    h1La->SetMarkerColor(kRed);
    h1La->SetMarkerStyle(kOpenCircle);
    h1La->SetMarkerSize(1);
    h1H3->SetLineColor(kBlue);
    h1H3->SetMarkerColor(kBlue);
    h1H3->SetMarkerSize(1);
    if (rebin>1) {
      h1H3->Rebin(rebin);
      h1La->Rebin(rebin);
    }
    /* Norm(h1La); */
    /* Norm(h1H3); */
    h1H3->GetXaxis()->SetTitle(xTitle.Data());
    h1H3->GetYaxis()->SetRangeUser(0 , h1H3->GetMaximum()*1.3);
    h1H3->GetXaxis()->SetRangeUser(2.97, 3.02);
    h1H3->GetXaxis()->SetNdivisions(206);
    /* h1H3->DrawCopy(drawstyle.Data()); */
    h1H3->DrawCopy("");
    h1La->DrawCopy((drawstyle+"same").Data()); 

    drawLatex( 0.18, 0.88, Form("%0.1f< p_{T}<%0.1f GeV/c", ptedge[i],ptedge[i+1]),0.055);
    TLegend* l = new TLegend(0.65,0.68,0.9,0.87);
    l->AddEntry( h1La, legtitle2.Data(), "pl");
    l->AddEntry(h1H3,legtitle1.Data(),"pl");
    l->Draw();
    drawLine(2.9913,0,2.9913,  h1H3->GetMaximum()*0.8,1.2,2,kBlack );
    /* drawLine(2.99,0,2.99,  h1H3->GetMaximum()*0.8,1.2,2,kBlack ); */
    drawLine(3,0,3,  h1H3->GetMaximum()*0.8,1.2,2,kBlack );
    /* h1La->GetYaxis()->SetRangeUser(0. , 2); */
    /* h1La->GetXaxis()->SetRangeUser(2.97 , 3.02); */
    /* h1La->Divide(h1H3); */
    /* h1La->Draw(); */
    /* #<{(| drawLine(2.97,1,3.02,1,1.2,2,kBlack ); |)}># */
    /* h1La->Fit("pol0"); */
  }
  c->cd();
  addpdf(pdf);
}
void projAndScale(TString name, TFile* fSig, TFile* fBk, double scale, TFile* fMc, TCanvas* c,TPDF* pdf,TString xTitle , TString drawstyle="", int rebin=1,  TString legtitle1="H3L quasi", TString legtitle2="#Lambda+d", TString text="dp#pi")
{
  fSig->cd();
  TH2F* h2H3L = (TH2F*)fSig->Get(name.Data());
  h2H3L->Sumw2();
  TH2F* h2Bk = (TH2F*)fBk->Get(name.Data());
  h2Bk->Sumw2();
  h2Bk->Scale(scale);
  h2H3L->Add(h2Bk);;
  h2H3L->Sumw2();
  h2H3L->SetDirectory(0);
  fMc->cd();
  TH2F* h2La = (TH2F*)fMc->Get(name.Data());
  h2La->Sumw2();
  h2La->SetDirectory(0);
  c->Clear();
  /* c->Divide(3,2); */
  /* int  const bins=6; */
  /* double ptedge[bins+1]={0,0.5,1,2,3,4,5}; */
  c->Divide(2,2);
  int  const bins=4;
  double ptedge[bins+1]={0,0.5,1,2,3};

  for (int i=0;i<bins;i++)
  {
    c->cd(i+1);
    TH1F* h1H3= (TH1F*)h2H3L->ProjectionY(Form("h1H3%d",i), h2H3L->GetXaxis()->FindBin(ptedge[i]),h2H3L->GetXaxis()->FindBin(ptedge[i+1]));
    TH1F* h1La= (TH1F*)h2La->ProjectionY(Form("h1La%d",i), h2La->GetXaxis()->FindBin(ptedge[i]),h2La->GetXaxis()->FindBin(ptedge[i+1]));
    h1La->SetLineColor(kRed);
    h1La->SetMarkerColor(kRed);
    h1La->SetMarkerStyle(kOpenCircle);
    h1La->SetMarkerSize(1);
    h1H3->SetLineColor(kBlue);
    h1H3->SetMarkerColor(kBlue);
    h1H3->SetMarkerSize(1);
    if (rebin>1) {
      h1H3->Rebin(rebin);
      h1La->Rebin(rebin);
    }
    Norm(h1La);
    Norm(h1H3);
    h1H3->GetXaxis()->SetTitle(xTitle.Data());
    /* h1H3->GetYaxis()->SetRangeUser(0 , h1H3->GetMaximum()*1.3); */
    /* h1H3->DrawCopy(drawstyle.Data()); */
    h1H3->DrawCopy("");
    h1La->DrawCopy((drawstyle+"same").Data()); 
    drawLatex( 0.18, 0.88, Form("%0.1f<(%s) p_{T}<%0.1f GeV/c", ptedge[i],text.Data(),ptedge[i+1]),0.055);
    TLegend* l = new TLegend(0.65,0.68,0.9,0.87);
    l->AddEntry( h1La, legtitle2.Data(), "pl");
    l->AddEntry(h1H3,legtitle1.Data(),"pl");
    l->Draw();
  }
  c->cd();
  addpdf(pdf);
}
void projSig(TString name, TFile* fH3L, TFile* fBk, TCanvas* c,TPDF* pdf,TString xTitle , TString drawstyle="", int rebin=1,  TString legtitle1="H3L", TString legtitle2="#Bkmbda")
{
  fH3L->cd();
  TH2F* h2H3L = (TH2F*)fH3L->Get(name.Data());
  h2H3L->Sumw2();
  h2H3L->SetDirectory(0);
  fBk->cd();
  TH2F* h2Bk = (TH2F*)fBk->Get(name.Data());
  h2Bk->Sumw2();
  h2Bk->SetDirectory(0);
  c->Clear();
  // c->Divide(1,1);
  // int  const bins=1;
  /* double ptedge[bins+1]={0,0.5,1,2,3,4,5}; */
  // double ptedge[bins+1]={0,5};
  c->Divide(2,2);
  int  const bins=4;
  double ptedge[bins+1]={0,0.5,1,2,3};

  for (int i=0;i<bins;i++)
  {
    c->cd(i+1);
    TH1F* h1H3= (TH1F*)h2H3L->ProjectionY("h1H3", h2H3L->GetXaxis()->FindBin(ptedge[i]),h2H3L->GetXaxis()->FindBin(ptedge[i+1]));
    TH1F* h1Bk= (TH1F*)h2Bk->ProjectionY("h1Bk", h2Bk->GetXaxis()->FindBin(ptedge[i]),h2Bk->GetXaxis()->FindBin(ptedge[i+1]));
    h1Bk->SetLineColor(kRed);
    h1Bk->SetMarkerColor(kRed);
    h1Bk->SetMarkerStyle(kOpenCircle);
    h1Bk->SetMarkerSize(1);
    h1H3->SetLineColor(kBlue);
    h1H3->SetMarkerColor(kBlue);
    h1H3->SetMarkerSize(1);
    if (rebin>1) {
      h1H3->Rebin(rebin);
      h1Bk->Rebin(rebin);
    }
    h1Bk->Scale(h1H3->Integral(h1H3->GetXaxis()->FindBin(3.04), h1H3->GetXaxis()->FindBin(3.15)) / h1Bk->Integral(h1Bk->GetXaxis()->FindBin(3.04), h1Bk->GetXaxis()->FindBin(3.15)));
    h1H3->GetXaxis()->SetTitle(xTitle.Data());
    h1H3->GetYaxis()->SetRangeUser(0 , h1H3->GetMaximum()*1.3);
    h1H3->GetXaxis()->SetRangeUser(2.97, 3.02);
    h1H3->GetXaxis()->SetNdivisions(206);
    /* h1H3->DrawCopy(drawstyle.Data()); */
    h1H3->DrawCopy("");
    h1Bk->DrawCopy((drawstyle+"same").Data()); 

    drawLatex( 0.18, 0.88, Form("%0.1f< p_{T}<%0.1f GeV/c", ptedge[i],ptedge[i+1]),0.055);
    TLegend* l = new TLegend(0.65,0.68,0.9,0.87);
    l->AddEntry( h1Bk, legtitle2.Data(), "pl");
    l->AddEntry(h1H3,legtitle1.Data(),"pl");
    l->Draw();
    drawLine(2.9913,0,2.9913,  h1H3->GetMaximum()*0.8,1.2,2,kBlack );
    /* drawLine(2.99,0,2.99,  h1H3->GetMaximum()*0.8,1.2,2,kBlack ); */
    drawLine(3,0,3,  h1H3->GetMaximum()*0.8,1.2,2,kBlack );
  }
  c->cd();
  addpdf(pdf);
}
void drawCompLaH3L()
{
  SetsPhenixStyle();
  /* TFile* fH3L= new TFile("fout_H3L3b_phase.root"); */
  TFile* fH3L= new TFile("fout_H3L3b_quasi.root");
  /* TFile* fLa = new TFile("fout_Lambda.root"); */
  TFile* fLa = new TFile("fout_Lambda_wLaD.root");

  TCanvas* c = new TCanvas("c","c");
  /* TPDF* pdf = new TPDF("Lambda_H3L_phase.pdf"); */
  /* TPDF* pdf = new TPDF("Lambda_H3L_quasi.pdf"); */
  TPDF* pdf = new TPDF("Lambda_H3L_LambdaDeu.pdf");
  pdf->Off();

  projAndComp("hptppimass", fH3L, fLa, c,pdf,"p#pi Mass (GeV/c^{2})", "plhist",1 ,"H3L quasi" ,"#Lambda+d","p#pi");
  projAndComp("hptppichi2ndf", fH3L, fLa, c,pdf,"p#pi #chi^{2}_{NDF}","" ,1 , "H3L quasi" ,"#Lambda+d","p#pi");
  projAndComp("hptppichi2prim", fH3L, fLa, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "H3L quasi" ,"#Lambda+d","p#pi");
  projAndComp("hptppil", fH3L, fLa, c,pdf,"(p#pi) l" ,"",1 , "H3L quasi" ,"#Lambda+d","p#pi");
  projAndComp("hptppildl", fH3L, fLa, c,pdf,"(p#pi) l/#deltal" ,"",1 , "H3L quasi" ,"#Lambda+d","p#pi");
  projAndComp("hptpichi2prim", fH3L, fLa, c,pdf,"#pi #chi^{2}_{prim}" ,"",1 , "H3L quasi" ,"#Lambda+d","p#pi");
  projAndComp("hptpchi2prim", fH3L, fLa, c,pdf,"p #chi^{2}_{prim}" ,"",1 , "H3L quasi" ,"#Lambda+d","p#pi");
  projAndComp("hptpidca", fH3L, fLa, c,pdf,"#pi DCA","" ,2, "H3L quasi" ,"#Lambda+d","p#pi");
  projAndComp("hptpdca", fH3L, fLa, c,pdf , "p DCA","" ,2, "H3L quasi" ,"#Lambda+d", "p#pi");
  projAndComp("hptsumdca", fH3L, fLa, c,pdf , "p+pi DCA","" ,2, "H3L quasi" ,"#Lambda+d", "p#pi");

  projAndComp("hptH3L_l", fH3L, fLa, c,pdf,"l", "p",4 );
  projAndComp("hptH3L_ldl", fH3L, fLa, c,pdf,"l/dl", "p",4 );
  projAndComp("hptH3L_dchi2prim", fH3L, fLa, c,pdf,"d chi2primary", "p",4 );
  projAndComp("hptH3L_pchi2prim", fH3L, fLa, c,pdf,"p chi2primary", "p",4 );
  projAndComp("hptH3L_pichi2prim", fH3L, fLa, c,pdf,"#pi chi2primary", "p",4);
  projAndComp("hptH3L_dDca", fH3L, fLa, c,pdf,"d DCA", "p",1);
  projAndComp("hptH3L_dpDca", fH3L, fLa, c,pdf,"dp pair DCA", "p",4);
  projAndComp("hptH3L_chi2topo", fH3L, fLa, c,pdf,"chi2topo", "p",4);
  projAndComp("hptH3L_chi2ndf", fH3L, fLa, c,pdf,"chi2NDF", "p",4);

  pdf->On();
  pdf->Close();
}

void drawCompQuasi()
{
  SetsPhenixStyle();
  TFile* fH3L= new TFile("fout_H3L3b_phase.root");
  TFile* fH3L3bquasi = new TFile("fout_H3L3b_quasi.root");

  TCanvas* c = new TCanvas("c","c");
  TPDF* pdf = new TPDF("H3L_phase_vs_quasi.pdf");
  pdf->Off();

  projAndComp("hptH3Lmass", fH3L, fH3L3bquasi, c,pdf,"H3L Mass (GeV/c^{2})", "plhist" );
  projAndCompH3L("hptppimass", fH3L, fH3L3bquasi, c,pdf,"p#pi Mass (GeV/c^{2})", "plhist" );
  projAndCompH3L("hptppichi2ndf", fH3L, fH3L3bquasi, c,pdf,"p#pi #chi^{2}_{NDF}" );
  projAndCompH3L("hptppichi2prim", fH3L, fH3L3bquasi, c,pdf,"(p#pi) #chi^{2}_{prim}" );
  projAndCompH3L("hptppil", fH3L, fH3L3bquasi, c,pdf,"(p#pi) l" );
  projAndCompH3L("hptppildl", fH3L, fH3L3bquasi, c,pdf,"(p#pi) l/#deltal" );
  projAndCompH3L("hptpichi2prim", fH3L, fH3L3bquasi, c,pdf,"#pi #chi^{2}_{prim}" );
  projAndCompH3L("hptpchi2prim", fH3L, fH3L3bquasi, c,pdf,"p #chi^{2}_{prim}" );
  projAndComp("hptpidca", fH3L, fH3L3bquasi, c,pdf,"#pi DCA","" ,2);
  projAndComp("hptpdca", fH3L, fH3L3bquasi, c,pdf , "p DCA","" ,2);
  projAndComp("hptsumdca", fH3L, fH3L3bquasi, c,pdf , "p+pi DCA","" ,2);

  pdf->On();
  pdf->Close();
}

void drawCompLaDeu()
{
  SetsPhenixStyle();
  TFile* fH3L= new TFile("fout_Lambda_wLaD.root");
  TFile* fH3L3bquasi = new TFile("fout_Lambda_wLaD_rotate.root");

  TCanvas* c = new TCanvas("c","c");
  TPDF* pdf = new TPDF("LambdaDeu_rotate.pdf");
  pdf->Off();

  projAndComp("hptH3Lmass", fH3L, fH3L3bquasi, c,pdf,"Mass(dp#pi) (GeV/c^{2})", "p",1,"#Lambda+d","#Lambda+Rotate d", "dp#pi");
  projAndComp("hptH3L_l", fH3L, fH3L3bquasi, c,pdf,"l", "p",4,"#Lambda+d","#Lambda+Rotate d", "dp#pi");
  projAndComp("hptH3L_ldl", fH3L, fH3L3bquasi, c,pdf,"l/dl", "p",4,"#Lambda+d","#Lambda+Rotate d", "dp#pi");
  projAndComp("hptH3L_dchi2prim", fH3L, fH3L3bquasi, c,pdf,"d chi2primary", "p",4,"#Lambda+d","#Lambda+Rotate d", "dp#pi");
  projAndComp("hptH3L_pchi2prim", fH3L, fH3L3bquasi, c,pdf,"p chi2primary", "p",4,"#Lambda+d","#Lambda+Rotate d", "dp#pi");
  projAndComp("hptH3L_pichi2prim", fH3L, fH3L3bquasi, c,pdf,"#pi chi2primary", "p",4,"#Lambda+d","#Lambda+Rotate d", "dp#pi");
  projAndComp("hptH3L_dDca", fH3L, fH3L3bquasi, c,pdf,"d DCA", "p",1,"#Lambda+d","#Lambda+Rotate d","dp#pi"  );
  projAndComp("hptH3L_dpDca", fH3L, fH3L3bquasi, c,pdf,"dp pair DCA", "p",4,"#Lambda+d","#Lambda+Rotate d", "dp#pi" );
  projAndComp("hptH3L_chi2topo", fH3L, fH3L3bquasi, c,pdf,"chi2topo", "p",4,"#Lambda+d","#Lambda+Rotate d", "dp#pi" );

  pdf->On();
  pdf->Close();
}

void drawMixData()
{
  SetsPhenixStyle();
  TCanvas* c = new TCanvas("c","c");
  TPDF* pdf = new TPDF("MixEventQA_check.pdf");
  pdf->Off();

  TFile *f1 = TFile::Open("fout_H3L_data_large.root"); 
  TH2F* h2sig = (TH2F*)f1->Get("hptH3Lmass")->Clone("hptH3Lmass_sig");
  TH1F* hsig = (TH1F*)h2sig->ProjectionY("hsig");
  
  TFile *f2 = TFile::Open("fout_H3L_data_ME_large.root"); 
  // TFile *f2 = TFile::Open("fout_H3L_data_SE.root"); 
  TH2F* h2bk = (TH2F*)f2->Get("hptH3Lmass")->Clone("hptH3Lmass_ME");
  TH1F* hbk = (TH1F*)h2bk->ProjectionY("hbk");
  
  double sig_sb = hsig->Integral(hsig->GetXaxis()->FindBin(2.97), hsig->GetXaxis()->FindBin(2.98)) + hsig->Integral(hsig->GetXaxis()->FindBin(3.005),  hsig->GetXaxis()->FindBin(3.03));
  // double sig_sb = hsig->Integral(hsig->GetXaxis()->FindBin(2.95), hsig->GetXaxis()->FindBin(2.98)) ;
  double bk_sb =  hbk->Integral(hbk->GetXaxis()->FindBin(2.97), hbk->GetXaxis()->FindBin(2.98)) + hbk->Integral(hbk->GetXaxis()->FindBin(3.005),  hbk->GetXaxis()->FindBin(3.03));
  // double bk_sb =  hbk->Integral(hbk->GetXaxis()->FindBin(2.95), hbk->GetXaxis()->FindBin(2.98)) ;
  double scale = sig_sb/bk_sb;
  cout <<1./scale << endl;
  hbk->Scale(scale);
  // hsig->Rebin();
  // hbk->Rebin();

  hsig->Draw();
  hsig->GetXaxis()->SetTitle("Mass(p#pid) (GeV/c^{2})");
  hsig->GetYaxis()->SetTitle("Counts");
  // hsig->GetXaxis()->SetRangeUser(2.98, 3.02);
  hsig->GetYaxis()->SetRangeUser(0, hsig->GetMaximum()*1.1);
  hbk->SetMarkerStyle(kOpenCircle);
  hbk->SetMarkerColor(kRed);
  hbk->SetLineColor(kRed);
  hbk->Draw("same");

  TH1F* hsig_bk = (TH1F*)hsig->Clone("hsig_bk");
  hsig_bk->Add(hbk,-1);
  hsig_bk->SetMarkerColor(kBlue);
  hsig_bk->SetMarkerStyle(kFullCircle);
  hsig_bk->Draw("same");

  TLegend* leg = new TLegend( 0.2 , 0.7 ,0.5,0.9  );
  leg->AddEntry(hbk, "mix-event(ME)","pl");
  leg->AddEntry(hsig, "same-event(SE)","pl");
  leg->Draw();
  drawLatex( 0.2,0.6,Form("ME/SE=%0.2f", 1./scale), 0.055);
 
  addpdf(pdf);
  
  projAndComp("hptH3L_l", f1, f2, c,pdf,"l", "p",4,"Sig","BK", "dp#pi");
  projAndComp("hptH3L_ldl", f1, f2, c,pdf,"l/dl", "p",4,"Sig","BK", "dp#pi");
  projAndComp("hptH3L_chi2ndf", f1, f2, c,pdf,"chi2NDF", "p",4,"Sig","BK", "dp#pi");
  projAndComp("hptH3L_chi2topo", f1, f2, c,pdf,"chi2topo", "p",4,"Sig","BK", "dp#pi" );
  projAndComp("hptH3L_dchi2prim", f1, f2, c,pdf,"d chi2primary", "p",4,"Sig","BK", "dp#pi");
  projAndComp("hptH3L_pchi2prim", f1, f2, c,pdf,"p chi2primary", "p",4,"Sig","BK", "dp#pi");
  projAndComp("hptH3L_pichi2prim", f1, f2, c,pdf,"#pi chi2primary", "p",4,"Sig","BK", "dp#pi");
  projAndComp("hptH3L_dDca", f1, f2, c,pdf,"d DCA", "p",1,"Sig","BK","dp#pi"  );
  projAndComp("hptH3L_piDca", f1, f2, c,pdf,"#pi DCA", "p",4,"Sig","BK", "dp#pi" );
  projAndComp("hptH3L_pDca", f1, f2, c,pdf,"#p DCA", "p",4,"Sig","BK", "dp#pi" );
  projAndComp("hptH3L_dpDca", f1, f2, c,pdf,"dp pair DCA", "p",4,"Sig","BK", "dp#pi" );

  projAndComp("hptppimass", f1, f2, c,pdf,"p#pi Mass (GeV/c^{2})", "plhist",1 ,"Sig" ,"BK","p#pi");
  projAndComp("hptppichi2ndf", f1, f2, c,pdf,"p#pi #chi^{2}_{NDF}","" ,1 , "Sig" ,"BK","p#pi");
  projAndComp("hptppichi2prim", f1, f2, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "Sig" ,"BK","p#pi");
  projAndComp("hptppil", f1, f2, c,pdf,"(p#pi) l" ,"",1 , "Sig" ,"BK","p#pi");
  projAndComp("hptppildl", f1, f2, c,pdf,"(p#pi) l/#deltal" ,"",1 , "Sig" ,"BK","p#pi");
  // projAndComp("hptpichi2prim", f1, f2, c,pdf,"#pi #chi^{2}_{prim}" ,"",1 , "Sig" ,"BK","p#pi");
  // projAndComp("hptpchi2prim", f1, f2, c,pdf,"p #chi^{2}_{prim}" ,"",1 , "Sig" ,"BK","p#pi");
  // projAndComp("hptpidca", f1, f2, c,pdf,"#pi DCA","" ,2, "Sig" ,"BK","p#pi");
  // projAndComp("hptpdca", f1, f2, c,pdf , "p DCA","" ,2, "Sig" ,"BK", "p#pi");
  // projAndComp("hptsumdca", f1, f2, c,pdf , "p+pi DCA","" ,2, "Sig" ,"BK", "p#pi");

  pdf->On();
  pdf->Close();

}
void drawSEData()
{
  SetsPhenixStyle();
  TCanvas* c = new TCanvas("c","c");
  TPDF* pdf = new TPDF("SameEventQA_check.pdf");
  pdf->Off();

  TFile *f1 = TFile::Open("fout_H3L_data_large.root"); 
  TH2F* h2sig = (TH2F*)f1->Get("hptH3Lmass")->Clone("hptH3Lmass_sig");
  TH1F* hsig = (TH1F*)h2sig->ProjectionY("hsig");
  
  // TFile *f2 = TFile::Open("fout_H3L_data_ME.root"); 
  TFile *f2 = TFile::Open("fout_H3L_data_SE_large.root"); 
  TH2F* h2bk = (TH2F*)f2->Get("hptH3Lmass")->Clone("hptH3Lmass_ME");
  TH1F* hbk = (TH1F*)h2bk->ProjectionY("hbk");
  
  double sig_sb = hsig->Integral(hsig->GetXaxis()->FindBin(2.95), hsig->GetXaxis()->FindBin(2.98)) + hsig->Integral(hsig->GetXaxis()->FindBin(3.0),  hsig->GetXaxis()->FindBin(3.04));
  // double sig_sb = hsig->Integral(hsig->GetXaxis()->FindBin(2.95), hsig->GetXaxis()->FindBin(2.98)) ;
  double bk_sb =  hbk->Integral(hbk->GetXaxis()->FindBin(2.95), hbk->GetXaxis()->FindBin(2.98)) + hbk->Integral(hbk->GetXaxis()->FindBin(3.0),  hbk->GetXaxis()->FindBin(3.04));
  // double bk_sb =  hbk->Integral(hbk->GetXaxis()->FindBin(2.95), hbk->GetXaxis()->FindBin(2.98)) ;
  double scale = sig_sb/bk_sb;
  cout <<1./scale << endl;
  hbk->Scale(scale);

  hsig->Draw();
  hsig->GetXaxis()->SetTitle("Mass(p#pid) (GeV/c^{2})");
  hsig->GetYaxis()->SetTitle("Counts");
  hbk->SetMarkerStyle(kOpenCircle);
  hbk->SetMarkerColor(kRed);
  hbk->SetLineColor(kRed);
  hbk->Draw("same");

  TLegend* leg = new TLegend( 0.2 , 0.7 ,0.5,0.9  );
  leg->AddEntry(hsig, "Official KF","pl");
  leg->AddEntry(hbk, "My Contruction","pl");
  leg->Draw();
  drawLatex( 0.2,0.6,Form("My/KF=%0.4f", 1./scale), 0.055);
 
  addpdf(pdf);
  
  projAndComp("hptH3L_l", f1, f2, c,pdf,"l", "p",4,"KF","My", "dp#pi");
  projAndComp("hptH3L_ldl", f1, f2, c,pdf,"l/dl", "p",4,"KF","My", "dp#pi");
  projAndComp("hptH3L_chi2ndf", f1, f2, c,pdf,"chi2NDF", "p",4,"KF","My", "dp#pi");
  projAndComp("hptH3L_chi2topo", f1, f2, c,pdf,"chi2topo", "p",4,"KF","My", "dp#pi" );
  projAndComp("hptH3L_dchi2prim", f1, f2, c,pdf,"d chi2primary", "p",4,"KF","My", "dp#pi");
  projAndComp("hptH3L_pchi2prim", f1, f2, c,pdf,"p chi2primary", "p",4,"KF","My", "dp#pi");
  projAndComp("hptH3L_pichi2prim", f1, f2, c,pdf,"#pi chi2primary", "p",4,"KF","My", "dp#pi");
  projAndComp("hptH3L_dDca", f1, f2, c,pdf,"d DCA", "p",1,"KF","My","dp#pi"  );
  projAndComp("hptH3L_piDca", f1, f2, c,pdf,"#pi DCA", "p",4,"KF","My", "dp#pi" );
  projAndComp("hptH3L_pDca", f1, f2, c,pdf,"#p DCA", "p",4,"KF","My", "dp#pi" );
  projAndComp("hptH3L_dpDca", f1, f2, c,pdf,"dp pair DCA", "p",4,"KF","My", "dp#pi" );

  projAndComp("hptppimass", f1, f2, c,pdf,"p#pi Mass (GeV/c^{2})", "plhist",1 ,"KF" ,"My","p#pi");
  projAndComp("hptppichi2ndf", f1, f2, c,pdf,"p#pi #chi^{2}_{NDF}","" ,1 , "KF" ,"My","p#pi");
  projAndComp("hptppichi2prim", f1, f2, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "KF" ,"My","p#pi");
  projAndComp("hptppil", f1, f2, c,pdf,"(p#pi) l" ,"",1 , "KF" ,"My","p#pi");
  projAndComp("hptppildl", f1, f2, c,pdf,"(p#pi) l/#deltal" ,"",1 , "KF" ,"My","p#pi");
  // projAndComp("hptpichi2prim", f1, f2, c,pdf,"#pi #chi^{2}_{prim}" ,"",1 , "KF" ,"My","p#pi");
  // projAndComp("hptpchi2prim", f1, f2, c,pdf,"p #chi^{2}_{prim}" ,"",1 , "KF" ,"My","p#pi");
  // projAndComp("hptpidca", f1, f2, c,pdf,"#pi DCA","" ,2, "KF" ,"My","p#pi");
  // projAndComp("hptpdca", f1, f2, c,pdf , "p DCA","" ,2, "KF" ,"My", "p#pi");
  // projAndComp("hptsumdca", f1, f2, c,pdf , "p+pi DCA","" ,2, "KF" ,"My", "p#pi");

  pdf->On();
  pdf->Close();

}
void drawData()
{
  drawMixData();
  drawSEData();
}
