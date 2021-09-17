#include "sPhenixStyle.h"
void Norm(TH1F* h)
{
  h->Scale(1./h->Integral()*h->GetBinWidth(1));
  h->Scale(1./h->GetMaximum());
  h->SetDirectory(0);
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
 
  projAndComp("hptH3L_l", fH3L, fH3L3bquasi, c,pdf,"l", "p",4 );
  projAndComp("hptH3L_ldl", fH3L, fH3L3bquasi, c,pdf,"l/dl", "p",4 );
  projAndComp("hptH3L_dchi2prim", fH3L, fH3L3bquasi, c,pdf,"d chi2primary", "p",4 );
  projAndComp("hptH3L_pchi2prim", fH3L, fH3L3bquasi, c,pdf,"p chi2primary", "p",4 );
  projAndComp("hptH3L_pichi2prim", fH3L, fH3L3bquasi, c,pdf,"#pi chi2primary", "p",4);
  projAndComp("hptH3L_dDca", fH3L, fH3L3bquasi, c,pdf,"d DCA", "p",1);
  projAndComp("hptH3L_dpDca", fH3L, fH3L3bquasi, c,pdf,"dp pair DCA", "p",4);
  projAndComp("hptH3L_chi2topo", fH3L, fH3L3bquasi, c,pdf,"chi2topo", "p",4);
  projAndComp("hptH3L_chi2ndf", fH3L, fH3L3bquasi, c,pdf,"chi2NDF", "p",4);

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

void drawComp()
{
  // drawCompLaH3L();
  drawCompQuasi();
}
