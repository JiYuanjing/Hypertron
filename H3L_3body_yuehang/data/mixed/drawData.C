#include "sPhenixStyle.h"

TH1F* hchi2ndf_Ld;
TH1F* hchi2ndf_H3L;
Double_t MyFitFun(double* x, double* p)
{
  double a = hchi2ndf_H3L->GetBinContent(hchi2ndf_H3L->GetXaxis()->FindBin(x[0]));
  double b = hchi2ndf_Ld->GetBinContent(hchi2ndf_Ld->GetXaxis()->FindBin(x[0]));
  return p[0]*a+p[1]*b;
}
bool Norm(TH1F* h)
{
  double scale = 1./h->Integral();
  h->Scale(fabs(scale));
  h->Scale(fabs(1./h->GetMaximum()));
  h->SetDirectory(0);
  h->GetYaxis()->SetTitle("Arb. Unit");
  if (scale>0) return true;
  else return false;
}

template <typename Hist>
void setHistStyle(Hist* h, int color, int markerstyle, double size , int mode = 0)
{
  h->SetLineColor(color);
  if (mode==0) //hist and graph
  {
    h->SetMarkerColor(color);
    h->SetMarkerStyle(markerstyle);
    h->SetMarkerSize(size);
  }
  else if (mode==1) //TF1
  {
    h->SetLineStyle(markerstyle);
    h->SetLineWidth(size);
  }
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
  // c->Divide(2,2);
  // int  const bins=4;
  // double ptedge[bins+1]={0,0.5,1,2,3};
  c->Divide(1,1);
  int  const bins=1;
  double ptedge[bins+1]={0,4};
  // double ptedge[bins+1]={1,2.5};

  for (int i=0;i<bins;i++)
  {
    c->cd(i+1);
    // TPad* p1 = (TPad*)c->cd(i+1);
    // p1->Divide(2,1,0,0);
    // TPad* p11 = (TPad*)p1->cd(1);
    // TPad* p12 = (TPad*)p1->cd(2);

    TPad*    p11 = new TPad("upperPad", "upperPad",.005, .4, .995, .995);
    TPad*    p12 = new TPad("lowerPad", "lowerPad", .005, .005, .995, .4);
    p11->Draw(); 			       
    p12->Draw();  
    p11->SetBottomMargin(0);
    p12->SetTopMargin(0);
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
    // h1H3->GetYaxis()->SetRangeUser(-0.1 , 1.15);
    h1H3->GetYaxis()->SetRangeUser(-0.1*h1H3->GetMaximum() , 1.3*h1H3->GetMaximum());
    /* h1H3->DrawCopy(drawstyle.Data()); */
    p11->cd();
    h1H3->DrawCopy("");
    h1La->DrawCopy((drawstyle+"same").Data()); 
    // TRatioPlot* rp = new TRatioPlot(h1La, h1H3);
    // gPad->SetTicks(0, 1);    rp->Draw();
    drawLatex( 0.18, 0.88, Form("%0.1f<(%s) p_{T}<%0.1f GeV/c", ptedge[i],text.Data(),ptedge[i+1]),0.055);
    TLegend* l = new TLegend(0.65,0.48,0.9,0.67);
    l->AddEntry( h1La, legtitle2.Data(), "pl");
    l->AddEntry(h1H3,legtitle1.Data(),"pl");
    l->Draw();
    p12->cd();
    TH1F* hratio = (TH1F*)h1H3->Clone("hratio");
    hratio->Divide(h1La);
    hratio->Draw();
    // hratio->GetYaxis()->SetRangeUser(hratio->GetMinimum()<0.95?hratio->GetMinimum():0.95,hratio->GetMaximum()>1.05?hratio->GetMaximum():1.05);
    hratio->GetYaxis()->SetRangeUser(0.5,1.2);
    hratio->GetYaxis()->SetNdivisions(206);
    hratio->GetYaxis()->SetTitleOffset(0.8);
    hratio->GetYaxis()->SetTitleSize(0.08);
    hratio->GetYaxis()->SetTitle(Form("%s/%s  ", legtitle1.Data(),legtitle2.Data()));
    hratio->GetXaxis()->SetTitleSize(0.08);
    hratio->GetXaxis()->SetTitleOffset(0.8);
    drawLine( hratio->GetXaxis()->GetXmin(), 1., hratio->GetXaxis()->GetXmax(), 1, 1.5, 9, 1);
  }
  c->cd();
  addpdf(pdf);
}
void projAndCompNoRatio(TString name, TFile* fH3L, TFile* fLa, TCanvas* c,TPDF* pdf,TString xTitle , TString drawstyle="", int rebin=1,  TString legtitle1="H3L quasi", TString legtitle2="#Lambda+d", TString text="dp#pi")
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
  // c->Divide(2,2);
  // int  const bins=4;
  // double ptedge[bins+1]={0,0.5,1,2,3};
  c->Divide(1,1);
  int  const bins=1;
  double ptedge[bins+1]={0,4};
  // double ptedge[bins+1]={1,2.5};

  for (int i=0;i<bins;i++)
  {
    c->cd(i+1);
    // TPad* p1 = (TPad*)c->cd(i+1);
    // p1->Divide(2,1,0,0);
    // TPad* p11 = (TPad*)p1->cd(1);
    // TPad* p12 = (TPad*)p1->cd(2);

    TPad*    p11 = new TPad("upperPad", "upperPad",.005, .4, .995, .995);
    TPad*    p12 = new TPad("lowerPad", "lowerPad", .005, .005, .995, .4);
    p11->Draw(); 			       
    p12->Draw();  
    p11->SetBottomMargin(0);
    p12->SetTopMargin(0);
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
    // h1H3->GetYaxis()->SetRangeUser(-0.1 , 1.15);
    h1H3->GetYaxis()->SetRangeUser(-0.1*h1H3->GetMaximum() , 1.3*h1H3->GetMaximum());
    /* h1H3->DrawCopy(drawstyle.Data()); */
    p11->cd();
    h1H3->DrawCopy("");
    h1La->DrawCopy((drawstyle+"same").Data()); 
    // TRatioPlot* rp = new TRatioPlot(h1La, h1H3);
    // gPad->SetTicks(0, 1);    rp->Draw();
    drawLatex( 0.18, 0.88, Form("%0.1f<(%s) p_{T}<%0.1f GeV/c", ptedge[i],text.Data(),ptedge[i+1]),0.055);
    TLegend* l = new TLegend(0.65,0.48,0.9,0.67);
    l->AddEntry( h1La, legtitle2.Data(), "pl");
    l->AddEntry(h1H3,legtitle1.Data(),"pl");
    l->Draw();
    p12->cd();
    TH1F* hratio = (TH1F*)h1H3->Clone("hratio");
    hratio->Divide(h1La);
    hratio->Draw();
    hratio->GetYaxis()->SetRangeUser(hratio->GetMinimum()<0.95?hratio->GetMinimum():0.95,hratio->GetMaximum()>1.05?hratio->GetMaximum():1.05);
    hratio->GetYaxis()->SetNdivisions(206);
    hratio->GetYaxis()->SetTitleOffset(0.8);
    hratio->GetYaxis()->SetTitleSize(0.08);
    hratio->GetYaxis()->SetTitle(Form("%s/%s  ", legtitle1.Data(),legtitle2.Data()));
    hratio->GetXaxis()->SetTitleSize(0.08);
    hratio->GetXaxis()->SetTitleOffset(0.8);
    drawLine( hratio->GetXaxis()->GetXmin(), 1., hratio->GetXaxis()->GetXmax(), 1, 1.5, 9, 1);
  }
  c->cd();
  addpdf(pdf);
}
void projAndScaleComp(TString name, double scale, TFile* fH3L, TFile* fH3Lbk, TFile* fLa, TCanvas* c,TPDF* pdf,TString xTitle , TString drawstyle="", int rebin=1,  TString legtitle1="H3L quasi", TString legtitle2="#Lambda+d", TString text="dp#pi", int scalemode=0, int signalmode=0)
{
  fH3L->cd();
  TH2F* h2H3Lsb = (TH2F*)fH3L->Get((name+"SB").Data())->Clone("hsb_s");
  TH2F* h2H3Lbksb = (TH2F*)fH3Lbk->Get((name+"SB").Data())->Clone("hsb_bk");
  double scale_auto = h2H3Lsb->Integral()/h2H3Lbksb->Integral();
  TH2F* h2H3L = (TH2F*)fH3L->Get((name+"Sig").Data())->Clone("hsig");
  h2H3L->Sumw2();
  h2H3L->SetDirectory(0);
  fH3Lbk->cd();
  TH2F* h2H3Lbk = (TH2F*)fH3Lbk->Get((name+"Sig").Data())->Clone("hbk");
  h2H3Lbk->Sumw2();
  h2H3Lbk->SetDirectory(0);
  // double sc = scale_auto;
  // double sc = scale_auto>scale? scale : scale_auto;
  double sc = scale;
  if (scalemode==1) sc = scale_auto;
  // h2H3Lbk->Scale(sc);
  h2H3Lbk->Scale(sc);
  h2H3L->Add(h2H3Lbk, -1);
  cout <<"ME scale: "<<scale<<" h2H3Lbk"<< scale_auto<< endl;

  fLa->cd();
  TH2F* h2La;
  if (signalmode==0) h2La= (TH2F*)fLa->Get((name+"Sig").Data());
  else if (signalmode==1) h2La =  (TH2F*)fLa->Get((name).Data());
  h2La->Sumw2();
  h2La->SetDirectory(0);
  h2La->Scale(sc);
  c->Clear();
  /* c->Divide(3,2); */
  /* int  const bins=6; */
  /* double ptedge[bins+1]={0,0.5,1,2,3,4,5}; */
  // c->Divide(2,2);
  // int  const bins=4;
  // double ptedge[bins+1]={0,0.5,1,2,3};
  c->Divide(1,1);
  int  const bins=1;
  // double ptedge[bins+1]={1,2.5};
  double ptedge[bins+1]={0,4};


  for (int i=0;i<bins;i++)
  {
    c->cd(i+1);
    TPad*    p11 = new TPad("upperPad", "upperPad",.005, .4, .995, .995);
    TPad*    p12 = new TPad("lowerPad", "lowerPad", .005, .005, .995, .4);
    p11->Draw(); 			       
    p12->Draw();  
    p11->SetBottomMargin(0);
    p12->SetTopMargin(0);

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
    if (!Norm(h1La)) cout<<"warning:"<< fLa->GetName()<< " "<<xTitle.Data()<<" is nagetive!" <<endl;
    if (!Norm(h1H3)) cout<<"warning:"<< fH3L->GetName()<< " "<<xTitle.Data()<<" is nagetive!" <<endl;
    h1H3->GetXaxis()->SetTitle(xTitle.Data());
    // h1H3->GetYaxis()->SetRangeUser(-0.1 , 1.15);
    h1H3->GetYaxis()->SetRangeUser(-0.1*h1H3->GetMaximum() , 1.2*h1H3->GetMaximum());
    /* h1H3->DrawCopy(drawstyle.Data()); */
    p11->cd();
    h1H3->DrawCopy("");
    h1La->DrawCopy((drawstyle+"same").Data()); 
    // TRatioPlot* rp = new TRatioPlot(h1La, h1H3);
    // gPad->SetTicks(0, 1);    rp->Draw();
    drawLatex( 0.18, 0.88, Form("%0.1f<(%s) p_{T}<%0.1f GeV/c", ptedge[i],text.Data(),ptedge[i+1]),0.055);
    TLegend* l = new TLegend(0.8,0.72,0.9,0.87);
    l->AddEntry( h1La, legtitle2.Data(), "pl");
    l->AddEntry(h1H3,legtitle1.Data(),"pl");
    l->Draw();
    p12->cd();
    TH1F* hratio = (TH1F*)h1H3->Clone("hratio");
    hratio->Divide(h1La);
    hratio->Draw();
    hratio->GetYaxis()->SetRangeUser(hratio->GetMinimum()*0.5,hratio->GetMaximum()*2);
    hratio->GetYaxis()->SetNdivisions(206);
    hratio->GetYaxis()->SetTitleOffset(0.8);
    hratio->GetYaxis()->SetTitleSize(0.08);
    hratio->GetYaxis()->SetTitle(Form("%s/%s  ", legtitle1.Data(),legtitle2.Data()));
    hratio->GetXaxis()->SetTitleSize(0.08);
    hratio->GetXaxis()->SetTitleOffset(0.8);
    drawLine( hratio->GetXaxis()->GetXmin(), 1., hratio->GetXaxis()->GetXmax(), 1, 1.5, 9, 1);
  }
  c->cd();
  addpdf(pdf);
}
void projAndScaleCompNoRatio(TString name, double scale, TFile* fH3L, TFile* fH3Lbk, TFile* fLa, TCanvas* c,TPDF* pdf,TString xTitle , TString drawstyle="", int rebin=1,  TString legtitle1="H3L quasi", TString legtitle2="#Lambda+d", TString text="dp#pi", int scalemode=0)
{
  fH3L->cd();
  TH2F* h2H3Lsb = (TH2F*)fH3L->Get((name+"SB").Data())->Clone("hsb_s");
  TH2F* h2H3Lbksb = (TH2F*)fH3Lbk->Get((name+"SB").Data())->Clone("hsb_bk");
  double scale_auto = h2H3Lsb->Integral()/h2H3Lbksb->Integral();
  TH2F* h2H3L = (TH2F*)fH3L->Get((name+"Sig").Data())->Clone("hsig");
  h2H3L->Sumw2();
  h2H3L->SetDirectory(0);
  fH3Lbk->cd();
  TH2F* h2H3Lbk = (TH2F*)fH3Lbk->Get((name+"Sig").Data())->Clone("hbk");
  h2H3Lbk->Sumw2();
  h2H3Lbk->SetDirectory(0);
  // double sc = scale_auto;
  // double sc = scale_auto>scale? scale : scale_auto;
  double sc = scale;
  if (scalemode==1) sc = scale_auto;
  // h2H3Lbk->Scale(sc);
  h2H3Lbk->Scale(sc);
  h2H3L->Add(h2H3Lbk, -1);
  cout <<"ME scale: "<<scale<<" h2H3Lbk"<< scale_auto<< endl;

  fLa->cd();
  TH2F* h2La = (TH2F*)fLa->Get((name+"Sig").Data());
  h2La->Sumw2();
  h2La->SetDirectory(0);
  h2La->Scale(sc);
  c->Clear();
  /* c->Divide(3,2); */
  /* int  const bins=6; */
  /* double ptedge[bins+1]={0,0.5,1,2,3,4,5}; */
  // c->Divide(2,2);
  // int  const bins=4;
  // double ptedge[bins+1]={0,0.5,1,2,3};
  c->Divide(1,1);
  int  const bins=1;
  // double ptedge[bins+1]={1,2.5};
  double ptedge[bins+1]={0,4};


  for (int i=0;i<bins;i++)
  {
    c->cd(i+1);
    // TPad*    p11 = new TPad("upperPad", "upperPad",.005, .4, .995, .995);
    // TPad*    p12 = new TPad("lowerPad", "lowerPad", .005, .005, .995, .4);
    // p11->Draw(); 			       
    // p12->Draw();  
    // p11->SetBottomMargin(0);
    // p12->SetTopMargin(0);

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
    if (!Norm(h1La)) cout<<"warning:"<< fLa->GetName()<< " "<<xTitle.Data()<<" is nagetive!" <<endl;
    if (!Norm(h1H3)) cout<<"warning:"<< fH3L->GetName()<< " "<<xTitle.Data()<<" is nagetive!" <<endl;
    h1H3->GetXaxis()->SetTitle(xTitle.Data());
    // h1H3->GetYaxis()->SetRangeUser(-0.1 , 1.15);
    h1H3->GetYaxis()->SetRangeUser(-0.1*h1H3->GetMaximum() , 1.2*h1H3->GetMaximum());
    /* h1H3->DrawCopy(drawstyle.Data()); */
    // p11->cd();
    h1H3->DrawCopy("");
    h1La->DrawCopy((drawstyle+"same").Data()); 
    // TRatioPlot* rp = new TRatioPlot(h1La, h1H3);
    // gPad->SetTicks(0, 1);    rp->Draw();
    drawLatex( 0.18, 0.88, Form("%0.1f<(%s) p_{T}<%0.1f GeV/c", ptedge[i],text.Data(),ptedge[i+1]),0.055);
    TLegend* l = new TLegend(0.8,0.72,0.9,0.87);
    l->AddEntry( h1La, legtitle2.Data(), "pl");
    l->AddEntry(h1H3,legtitle1.Data(),"pl");
    l->Draw();
    // p12->cd();
    // TH1F* hratio = (TH1F*)h1H3->Clone("hratio");
    // hratio->Divide(h1La);
    // hratio->Draw();
    // hratio->GetYaxis()->SetRangeUser(hratio->GetMinimum()*0.5,hratio->GetMaximum()*2);
    // hratio->GetYaxis()->SetNdivisions(206);
    // hratio->GetYaxis()->SetTitleOffset(0.8);
    // hratio->GetYaxis()->SetTitleSize(0.08);
    // hratio->GetYaxis()->SetTitle(Form("%s/%s  ", legtitle1.Data(),legtitle2.Data()));
    // hratio->GetXaxis()->SetTitleSize(0.08);
    // hratio->GetXaxis()->SetTitleOffset(0.8);
    // drawLine( hratio->GetXaxis()->GetXmin(), 1., hratio->GetXaxis()->GetXmax(), 1, 1.5, 9, 1);
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
double calpurity(TString fithistname, double cutslow, double cutshigh, double lowpt, double highpt,  TFile* f1, TFile* f2, double scale, TFile* fMc, TFile* fMc_ld, TCanvas* c, TPDF* pdf, TString xTitle)
{
  c->cd();
  TPad*    p11 = new TPad("upperPad", "upperPad",.005, .4, .995, .995);
  TPad*    p12 = new TPad("lowerPad", "lowerPad", .005, .005, .995, .4);
  p11->Draw(); 			       
  p12->Draw();  
  p11->SetBottomMargin(0);
  p12->SetTopMargin(0);
  p12->SetBottomMargin(0.18);

  p11->cd();
  //fiting test
  TH2F* h2chi2ndf_Ld = (TH2F*)fMc_ld->Get(fithistname.Data())->Clone("h2chi2ndf_Ld");
  h2chi2ndf_Ld->SetDirectory(0);
  hchi2ndf_Ld = (TH1F*)h2chi2ndf_Ld->ProjectionY("hchi2ndf_Ld", h2chi2ndf_Ld->GetXaxis()->FindBin(lowpt), h2chi2ndf_Ld->GetXaxis()->FindBin(highpt));
  hchi2ndf_Ld->Rebin(2);
  
  TH2F* h2chi2ndf_H3L = (TH2F*)fMc->Get(fithistname.Data())->Clone("h2chi2ndf_H3L");
  h2chi2ndf_H3L->SetDirectory(0);
  hchi2ndf_H3L = (TH1F*)h2chi2ndf_H3L->ProjectionY("hchi2ndf_H3L", h2chi2ndf_H3L->GetXaxis()->FindBin(lowpt), h2chi2ndf_H3L->GetXaxis()->FindBin(highpt));
  hchi2ndf_H3L->Rebin(2);
 
  Norm(hchi2ndf_H3L);
  Norm(hchi2ndf_Ld);
  TH2F* h2chi2ndf_Sig = (TH2F*)f1->Get((fithistname+"Sig").Data())->Clone("h2chi2ndf_Sig");
  h2chi2ndf_Sig->SetDirectory(0);
  TH2F* h2chi2ndf_Bk = (TH2F*)f2->Get((fithistname+"Sig").Data())->Clone("h2chi2ndf_Bk");
  h2chi2ndf_Bk->SetDirectory(0);
  h2chi2ndf_Sig->Add( h2chi2ndf_Bk , -1*scale ); 
  TH1F* hchi2ndf_Sig = (TH1F*)h2chi2ndf_Sig->ProjectionY("hchi2ndf_Sig", h2chi2ndf_Sig->GetXaxis()->FindBin(lowpt), h2chi2ndf_Sig->GetXaxis()->FindBin(highpt));
  hchi2ndf_Sig->Rebin(2);

  setHistStyle(hchi2ndf_Sig, kBlue, kFullCircle, 1.5);
  hchi2ndf_Sig->Draw();
  TH1F* hRatio = (TH1F*)hchi2ndf_Sig->Clone("hratio");
  TF1* fMyFit = new TF1("fMyFit", MyFitFun, 0,15,2);
  fMyFit->SetLineColor(kRed);
  fMyFit->SetParLimits(0, 0, 100000);
  fMyFit->SetParLimits(1, 0, 100000);
  hchi2ndf_Sig->Fit(fMyFit,"B");
  TH1F* hLd = (TH1F*)hchi2ndf_Ld->Clone("hLd");
  hLd->Scale(fMyFit->GetParameter(1));
  hLd->Draw("same");
  TH1F* hH3L = (TH1F*)hchi2ndf_H3L->Clone("hH3L");
  hH3L->Scale(fMyFit->GetParameter(0));
  setHistStyle(hLd, kBlack, kOpenStar, 2);
  setHistStyle(hH3L, kRed, kOpenCircle, 1.5);
  hH3L->Draw("same");
  gPad->SetLogy();
  double purity = hH3L->Integral(hH3L->GetXaxis()->FindBin(cutslow+1e-6),hH3L->GetXaxis()->FindBin(cutshigh))/fMyFit->Integral(cutslow,cutshigh)*hH3L->GetBinWidth(1);
  drawLatex(0.65,0.6,Form("purity = %0.2f",purity),0.055);
  TLegend* legF = new TLegend(0.65,0.65,0.9,0.9);
  legF->AddEntry(hH3L ,"H3L MC", "pe");
  legF->AddEntry(hLd ,"#Lambda+d MC", "pe");
  legF->AddEntry(hchi2ndf_Sig ,"Data", "pe");
  legF->AddEntry(fMyFit,"fit", "l");
  legF->Draw();
  p12->cd();
  for (int i=1;i<hRatio->GetNbinsX();i++)
  {
    if (hRatio->GetBinContent(i)>0) {
      hRatio->SetBinContent(i, fMyFit->Eval(hRatio->GetBinCenter(i))/hRatio->GetBinContent(i) );
      hRatio->SetBinError(i, fMyFit->Eval(hRatio->GetBinCenter(i))/pow(hRatio->GetBinError(i),2 ));
    }
  }
  hRatio->GetYaxis()->SetTitle("Fit/Data");
  hRatio->GetYaxis()->SetTitleOffset(0.8);
  hRatio->GetYaxis()->SetRangeUser(0.5,1.5);
  hRatio->GetYaxis()->SetTitleSize(0.075);
  hRatio->GetYaxis()->SetLabelSize(0.075);
  hRatio->GetXaxis()->SetTitleSize(0.075);
  hRatio->GetXaxis()->SetLabelSize(0.075);
  hRatio->GetXaxis()->SetTitleOffset(0.9);
  hRatio->Draw();
  drawLine(hRatio->GetXaxis()->GetXmin(), 1, hRatio->GetXaxis()->GetXmax(),1,1.5,2,1);
  addpdf(pdf);
  // delete hchi2ndf_Ld;
  // delete hchi2ndf_H3L;
  return purity;
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
  double highpt = 4, lowpt = 0., lowy=-1.5, highy = 0;
  // double highpt = 2.5, lowpt = 1, lowy=-0.9, highy = -0.1;
  // double highpt = 4, lowpt = 0, lowy=-1.5, highy = 0;
  TString histname="hH3LMassPtY";
  // TString histname="hH3LMassPtY_5_40";
  SetsPhenixStyle();
  TCanvas* c = new TCanvas("c1","c1");
  // TPDF* pdf = new TPDF("MixEventQA_check.pdf");
  TPDF* pdf = new TPDF("MixEventQA_Jul25.pdf");
  pdf->Off();

  // TFile *f1 = TFile::Open("fout_H3L_data_SE_large.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_large.root"); 

  // TFile *f1 = TFile::Open("out_KF_test/fout_H3L_KF_test.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_KF_Jul22_part.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_KF_Jul25.root"); 
  TFile *f1 = TFile::Open("fout_H3L_data_SE_Jul25.root"); 
  // TH2F* h2sig = (TH2F*)f1->Get("hptH3Lmass")->Clone("hptH3Lmass_sig");
  // h2sig->SetDirectory(0);
  // TH1F* hsig = (TH1F*)h2sig->ProjectionY("hsig");
  // TH3F* h2sig = (TH3F*)f1->Get("hH3LMassPtY_5_40")->Clone("hptH3Lmass_sig");
  TH3F* h2sig = (TH3F*)f1->Get(histname.Data())->Clone("hptH3Lmass_sig");
  h2sig->SetDirectory(0);
  // TH3F* h2sig->Project3D("xz");
  TH1F* hsig = (TH1F*)h2sig->ProjectionY("hsig", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2sig->GetZaxis()->FindBin(lowy), h2sig->GetZaxis()->FindBin(highy));
  hsig->SetDirectory(0);
  TH1F* hcent_se = (TH1F*)f1->Get("hcent")->Clone("hcent_se");
  double nEvents_se = hcent_se->Integral();

  // TFile *f2 = TFile::Open("fout_H3L_data_ME_large.root"); 
  // TFile *f2 = TFile::Open("fout_H3L_ME_Jul22_part.root"); 
  TFile *f2 = TFile::Open("fout_H3L_ME_Jul25.root"); 
  // TFile *f2 = TFile::Open("fout_H3L_data_SE.root"); 
  // TH2F* h2bk = (TH2F*)f2->Get("hptH3Lmass")->Clone("hptH3Lmass_ME");
  // TH3F* h2bk = (TH3F*)f2->Get("hH3LMassPtY_5_40")->Clone("hptH3Lmass_bk");
  TH3F* h2bk = (TH3F*)f2->Get(histname.Data())->Clone("hptH3Lmass_bk");
  h2bk->SetDirectory(0);
  TH1F* hbk = (TH1F*)h2bk->ProjectionY("hbk", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2bk->GetZaxis()->FindBin(lowy), h2bk->GetZaxis()->FindBin(highy));
  hbk->SetDirectory(0);

  // TFile* f3 = TFile::Open("fout_H3L_data_RT_Jul22_full.root");
  TFile* f3 = TFile::Open("fout_H3L_data_RT_Jul25.root");
  // TFile* f3 = TFile::Open("out_RT_test/fout_H3L_RT_test.root");
  // TH1F* hcent_rt= (TH1F*)f3->Get("hcent");
  // double nEvents_rt = hcent_rt->Integral(); 
  // TH2F* h2rt = (TH2F*)f3->Get("hptH3Lmass")->Clone("hptH3Lmass_RT");
  // h2rt->SetDirectory(0);
  // TH1F* hrt = (TH1F*)h2rt->ProjectionY("hrt");
  TH3F* h2rt = (TH3F*)f3->Get(histname.Data())->Clone("hptH3Lmass_rt");
  h2rt->SetDirectory(0);
  TH1F* hrt = (TH1F*)h2rt->ProjectionY("hrt", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2rt->GetZaxis()->FindBin(lowy), h2rt->GetZaxis()->FindBin(highy));
  hrt->SetDirectory(0);
  // hrt->Scale(nEvents_se/(1.*nEvents_rt) );
  // cout << nEvents_se/(1.*nEvents_rt)<<endl;
  //
  //scale
  // double sig_sb = hsig->Integral(hsig->GetXaxis()->FindBin(2.96), hsig->GetXaxis()->FindBin(2.98)) ;
  // double sig_sb = hsig->Integral(hsig->GetXaxis()->FindBin(2.97), hsig->GetXaxis()->FindBin(2.98)) + hsig->Integral(hsig->GetXaxis()->FindBin(3.0),  hsig->GetXaxis()->FindBin(3.02));
  double sig_sb =  hsig->Integral(hsig->GetXaxis()->FindBin(3.0),  hsig->GetXaxis()->FindBin(3.02));
  // double rt_sb =  hrt->Integral(hrt->GetXaxis()->FindBin(2.97), hrt->GetXaxis()->FindBin(2.98)) + hrt->Integral(hrt->GetXaxis()->FindBin(3.0),  hrt->GetXaxis()->FindBin(3.02));
  double rt_sb = hrt->Integral(hrt->GetXaxis()->FindBin(3.0),  hrt->GetXaxis()->FindBin(3.02));
  // double rt_sb =  hrt->Integral(hrt->GetXaxis()->FindBin(2.96), hrt->GetXaxis()->FindBin(2.98)) ;
  double scale_rt = sig_sb/rt_sb;
  // sig_sb = hsig->Integral(hsig->GetXaxis()->FindBin(2.97), hsig->GetXaxis()->FindBin(2.98)) + hsig->Integral(hsig->GetXaxis()->FindBin(3.0),  hsig->GetXaxis()->FindBin(3.02));
  double bk_sb =   hbk->Integral(hbk->GetXaxis()->FindBin(3.0),  hbk->GetXaxis()->FindBin(3.02));
  // double bk_sb =  hbk->Integral(hbk->GetXaxis()->FindBin(2.97), hbk->GetXaxis()->FindBin(2.98)) + hbk->Integral(hbk->GetXaxis()->FindBin(3.0),  hbk->GetXaxis()->FindBin(3.02));
  // double rt_sb =  hrt->Integral(hrt->GetXaxis()->FindBin(2.97), hrt->GetXaxis()->FindBin(2.98)) + hrt->Integral(hrt->GetXaxis()->FindBin(3.01),  hrt->GetXaxis()->FindBin(3.02));
  // double bk_sb =  hbk->Integral(hbk->GetXaxis()->FindBin(2.96), hbk->GetXaxis()->FindBin(2.98)) ;
  double scale = sig_sb/bk_sb;
  cout<<"ME scale: " <<1./scale << endl;
  hbk->Scale(scale);
  hrt->Scale(scale_rt);
  cout <<"rotation: " <<sig_sb/rt_sb<< endl;

  hsig->Draw();
  hsig->GetXaxis()->SetTitle("Mass(p#pid) (GeV/c^{2})");
  hsig->GetYaxis()->SetTitle("Counts");
  hsig->GetYaxis()->SetRangeUser(-0.1*hsig->GetMaximum(), hsig->GetMaximum()*1.1);
  
  setHistStyle(hbk, kRed, kOpenCircle, 1.5);
  hbk->Draw("same");
  setHistStyle(hrt, kGreen+2, kDiamond, 1.5);
  hrt->Draw("same");
  TLegend* leg_sig = new TLegend(0.65,0.25,0.88,0.45);
  leg_sig->AddEntry(hbk, "ME", "pl");
  leg_sig->AddEntry(hrt, "RT", "pl");
  leg_sig->AddEntry(hsig, "SE", "pl");
  leg_sig->Draw();

  drawLatex( 0.65,0.61,Form("%0.1f<y<%0.1f",lowy, highy ), 0.055);
  drawLatex( 0.65,0.54,Form("%0.1f<p_{T}<%0.1f GeV/c^{2}", lowpt, highpt), 0.055);
  // drawLatex( 0.65,0.47,Form("5-40%s", "%"), 0.055);
  drawLatex( 0.65,0.47,Form("0-50%s", "%"), 0.055);
  drawBox( 2.97, hsig->GetMinimum(),2.98, hsig->GetMaximum()*0.5 , kBlue-9, 1001 , 1,0.3); 
  // drawBox( 2.989, hsig_bk->GetMinimum(),2.995, hsig_bk->GetMaximum()*0.5 , kOrange-4, 1001,1,0.3); 
  drawBox( 3., hsig->GetMinimum(),3.02, hsig->GetMaximum()*0.5 , kBlue-9, 1001,1,0.3 ); 

  addpdf(pdf);

  c->cd();
  // calculate the significance
  TH1F* hsig_bk = (TH1F*)hsig->Clone("hsig_bk");
  hsig_bk->Add(hbk,-1);
  setHistStyle(hsig_bk, kBlue, kFullCircle, 1.5);
  TH1F* hsig_rt = (TH1F*)hsig->Clone("hsig_rt");
  hsig_rt->Add(hrt,-1);
  setHistStyle(hsig_rt, kGreen+2, kDiamond, 1.5);
  hsig_rt->Rebin();
  hsig_bk->Rebin();
  hsig_bk->Draw();

  // TF1* fit = new TF1("fit" ,"gaus(0)+pol1(3)", 2.97,3.02 );
  TF1* fit = new TF1("fit" ,"gausn(0)+pol1(3)", 2.97,3.02 );
  TF1* resfit = new TF1("resfit" ,"pol1", 2.95,3.05 );
  hsig_bk->GetXaxis()->SetRangeUser(2.97,2.985);
  hsig_bk->Fit(resfit,"R");
  fit->SetLineColor(kRed);
  double yield_bc = hsig_bk->Integral(hsig_bk->GetXaxis()->FindBin(2.98), hsig_bk->GetXaxis()->FindBin(3.));
  // double para[5]={yield_bc*hsig_bk->GetBinWidth(1)/sqrt(2*3.1415), 2.991, 0.0015,  resfit->GetParameter(0), resfit->GetParameter(1)};
  double para[5]={yield_bc*hsig_bk->GetBinWidth(1), 2.991, 0.0014,  resfit->GetParameter(0), resfit->GetParameter(1)};
  fit->SetParameters(para);
  double lowx=2.97 ,highx =3.02;
  hsig_bk->GetXaxis()->SetRangeUser(lowx,highx);
  hsig_rt->GetXaxis()->SetRangeUser(lowx,highy);
  hsig_bk->Draw("same");
  hsig_rt->Draw("same");
  hsig_bk->Fit(fit,"R");
  resfit->SetParameter(0, fit->GetParameter(3));
  resfit->SetParameter(1, fit->GetParameter(4));
  resfit->Draw("same");
  setHistStyle(resfit, kRed-2, 9, 2.5 ,1);
  drawLine(lowx, 0, highx, 0, 1.5, 2, 1 );
  // cout<<"binwidth: "<< hsig_bk->GetBinWidth(1)<< endl;
  TF1* fit_rt = (TF1*)fit->Clone("fit_rt");
  fit_rt->SetParameters(para);
  TF1* resfit_rt = (TF1*)resfit->Clone("resfit_rt");
  setHistStyle(resfit_rt, kGreen+2, 9, 2.5 ,1);
  setHistStyle(fit_rt, kGreen+2, 9, 2.5 ,1);
  hsig_rt->GetXaxis()->SetRangeUser(lowx,highx);
  hsig_rt->Fit(resfit_rt,"R");
  fit_rt->SetParameter(3, resfit_rt->GetParameter(0) );
  fit_rt->SetParameter(4, resfit_rt->GetParameter(1) );
  hsig_rt->Fit(fit_rt,"R");
  resfit_rt->SetParameter(0, fit_rt->GetParameter(3));
  resfit_rt->SetParameter(1, fit_rt->GetParameter(4));
  resfit_rt->Draw("same");

  double sigma = fit->GetParameter(2);
  double mean = fit->GetParameter(1);
  // double yield_me = fit->GetParameter(0)/hsig_bk->GetBinWidth(1)*fit->GetParameter(2)*sqrt(2*3.1415);
  double yield_me = fit->GetParameter(0)/hsig_bk->GetBinWidth(1);
  // double yield_rt = fit_rt->GetParameter(0)/hsig_rt->GetBinWidth(1)*fit_rt->GetParameter(2)*sqrt(2*3.1415);
  double yield_rt = fit_rt->GetParameter(0)/hsig_rt->GetBinWidth(1);
  double yield_counts = hsig_bk->Integral(hsig_bk->GetXaxis()->FindBin(mean-2.5*sigma), hsig_bk->GetXaxis()->FindBin(mean+2.5*sigma));
  cout<<"rotate: " << yield_rt<<" ME: "<<yield_me<<" bin counting: "<<yield_counts << endl;
  double bk_counts = hbk->Integral(hbk->GetXaxis()->FindBin(mean-2.5*sigma), hbk->GetXaxis()->FindBin(mean+2.5*sigma));
  double sp_counts = hsig->Integral(hbk->GetXaxis()->FindBin(mean-2.5*sigma), hbk->GetXaxis()->FindBin(mean+2.5*sigma));
  // double significance = yield_counts/sqrt(yield_counts+bk_counts);
  // double significance = yield_me/sqrt(sp_counts);
  double significance = yield_rt/sqrt(sp_counts);
  double s_me = yield_me/(fit->GetParError(0)/hsig_bk->GetBinWidth(1));
  double s_rt = yield_rt/(fit_rt->GetParError(0)/hsig_rt->GetBinWidth(1));
  cout<<"significance: " <<significance << endl;

  TLegend* leg = new TLegend( 0.72, 0.68 ,0.9,0.9 );
  // leg->AddEntry(hbk, "mix-event(ME)","pl");
  // leg->AddEntry(hrt, "rotate d(RT) (scale)","pl");
  leg->AddEntry(hsig, "SE","pl");
  leg->AddEntry(hsig_bk, "SE-ME","pl");
  leg->AddEntry(hsig_rt, "SE-RT","pl");
  leg->Draw();
  drawLatex( 0.2,0.82,Form("ME/SE=%0.2f", 1./scale), 0.055);
  drawLatex( 0.2,0.75,Form("Yield=%0.2f", yield_me), 0.055);
  drawLatex( 0.2,0.68,Form("#sigma=%0.2f MeV", sigma*1000.), 0.055);
  drawLatex( 0.2,0.61,Form("nEvents=%0.0f M", nEvents_se/1e6), 0.055);
  drawLatex( 0.2,0.54,Form("S/#sqrt{S+B}=%0.0f", significance), 0.055);
  drawLatex( 0.2,0.47,Form("S/#DeltaS=%0.0f (ME)", s_me), 0.055);
  // drawLatex( 0.2,0.4,Form("S/#DeltaS=%0.0f (RT)", s_rt), 0.055);
  drawLatex( 0.2,0.4,Form("Mean=%0.3f", mean), 0.055);
  drawLatex( 0.62,0.61,Form("%0.1f<y<%0.1f",lowy, highy ), 0.055);
  drawLatex( 0.62,0.54,Form("%0.1f<p_{T}<%0.1f GeV/c^{2}", lowpt, highpt), 0.055);
  // drawLatex( 0.62,0.47,Form("5-40%s", "%"), 0.055);
  drawLatex( 0.62,0.47,Form("0-50%s", "%"), 0.055);
  drawBox( 2.97, hsig_bk->GetMinimum(),2.98, hsig_bk->GetMaximum()*0.5 , kBlue-9, 1001 , 1,0.3); 
  drawBox( 2.989, hsig_bk->GetMinimum(),2.995, hsig_bk->GetMaximum()*0.5 , kOrange-4, 1001,1,0.3); 
  drawBox( 3., hsig_bk->GetMinimum(),3.02, hsig_bk->GetMaximum()*0.5 , kBlue-9, 1001,1,0.3 ); 

  addpdf(pdf);

  hsig_bk->Draw();
  // hsig_rt->Draw("same");
  TLegend* leg2 = new TLegend( 0.2 , 0.7 ,0.4,0.9  );
  leg2->AddEntry(hsig_rt, "SE-RT","pl");
  leg2->AddEntry(hsig_bk, "SE-ME","pl");
  leg2->Draw();

  addpdf(pdf);

  // check the difference between rotation and mix event
  c->Clear(); 
  c->Divide(1,2); 
  c->cd(1);
  TH1F* hrtcn = (TH1F*)hrt->Clone("hrtcn");
  TH1F* hbkcn = (TH1F*)hbk->Clone("hbkcn");
  TH1F* hsigcn = (TH1F*)hsig->Clone("hsigcn");
  hsigcn->Rebin(4);
  hrtcn->Rebin(4);
  hbkcn->Rebin(4);
  hrtcn->Divide(hsigcn);
  hbkcn->Divide(hsigcn);
  hrtcn->Draw();
  hbkcn->Draw("same");
  hrtcn->GetYaxis()->SetTitleOffset(0.8);
  hrtcn->GetYaxis()->SetRangeUser(0.8,1.2);
  hrtcn->GetYaxis()->SetTitle("BK/Sig");
  hrtcn->GetXaxis()->SetTitle("Mass(dp#pi) GeV/c^{2}");
  TLegend* leg_r = new TLegend( 0.7, 0.5, 0.9, 0.7); 
  leg_r->AddEntry( hbkcn, "MixEvent", "pl");
  leg_r->AddEntry(hrtcn, "Rotation", "pl");
  leg_r->Draw();
  drawLine(2.987,0.8, 2.987, 1.2, 1.5, 2, 1  );
  drawLine(2.997,0.8, 2.997, 1.2, 1.5, 2, 1  );
  drawLine(2.95,1, 3.05, 1, 1.5, 9, 1  );
  // addpdf(pdf);

  c->cd(2);
  TH1F* hrtcn2 = (TH1F*)hrt->Clone("hrtcn2");
  TH1F* hbkcn2 = (TH1F*)hbk->Clone("hbkcn2");
  hrtcn2->Rebin(4);
  hbkcn2->Rebin(4);
  hrtcn2->Divide(hbkcn2);
  hrtcn2->Draw();
  hrtcn2->GetYaxis()->SetRangeUser(0.8, 1.2);
  hrtcn2->GetYaxis()->SetTitleOffset(0.8);
  hrtcn2->GetYaxis()->SetTitle("RT/ME");
  hrtcn2->GetXaxis()->SetTitle("Mass(dp#pi) GeV/c^{2}");
  drawLine(2.987,0.8, 2.987, 1.2, 1.5, 2, 1  );
  drawLine(2.997,0.8, 2.997, 1.2, 1.5, 2, 1  );
  drawLine(2.95,1, 3.05, 1, 1.5, 9, 1  );
  c->cd();
  addpdf(pdf);

  c = new TCanvas( "c","c", 800,800);
  projAndComp("hptH3L_lSB", f1, f2, c,pdf,"l", "p",4,"Sig_SB","ME_SB", "dp#pi");
  projAndComp("hptH3L_ldlSB", f1, f2, c,pdf,"l/dl", "p",4,"Sig_SB","ME_SB", "dp#pi");
  projAndComp("hptH3L_chi2ndfSB", f1, f2, c,pdf,"chi2NDF", "p",4,"Sig_SB","ME_SB", "dp#pi");
  projAndComp("hptH3L_chi2topoSB", f1, f2, c,pdf,"chi2topo", "p",4,"Sig_SB","ME_SB", "dp#pi" );
  projAndComp("hptH3L_dchi2primSB", f1, f2, c,pdf,"d chi2primary", "p",4,"Sig_SB","ME_SB", "dp#pi");
  projAndComp("hptH3L_pchi2primSB", f1, f2, c,pdf,"p chi2primary", "p",4,"Sig_SB","ME_SB", "dp#pi");
  projAndComp("hptH3L_pichi2primSB", f1, f2, c,pdf,"#pi chi2primary", "p",4,"Sig_SB","ME_SB", "dp#pi");
  projAndComp("hptH3L_dDcaSB", f1, f2, c,pdf,"d DCA", "p",1,"Sig_SB","ME_SB","dp#pi"  );
  projAndComp("hptH3L_piDcaSB", f1, f2, c,pdf,"#pi DCA", "p",4,"Sig_SB","ME_SB", "dp#pi" );
  projAndComp("hptH3L_pDcaSB", f1, f2, c,pdf,"#p DCA", "p",4,"Sig_SB","ME_SB", "dp#pi" );
  projAndComp("hptH3L_dpDcaSB", f1, f2, c,pdf,"dp pair DCA", "p",4,"Sig_SB","ME_SB", "dp#pi" );
  projAndComp("hptppimassSB", f1, f2, c,pdf,"p#pi Mass", "p",4,"Sig_SB","ME_SB", "dp#pi" );
  //
  projAndComp("hptH3L_lSB", f1, f3, c,pdf,"l", "p",4,"Sig_SB","RT_SB", "dp#pi");
  projAndComp("hptH3L_ldlSB", f1, f3, c,pdf,"l/dl", "p",4,"Sig_SB","RT_SB", "dp#pi");
  projAndComp("hptH3L_chi2ndfSB", f1, f3, c,pdf,"chi2NDF", "p",4,"Sig_SB","RT_SB", "dp#pi");
  projAndComp("hptH3L_chi2topoSB", f1, f3, c,pdf,"chi2topo", "p",4,"Sig_SB","RT_SB", "dp#pi" );
  projAndComp("hptH3L_dchi2primSB", f1, f3, c,pdf,"d chi2primary", "p",4,"Sig_SB","RT_SB", "dp#pi");
  projAndComp("hptH3L_pchi2primSB", f1, f3, c,pdf,"p chi2primary", "p",4,"Sig_SB","RT_SB", "dp#pi");
  projAndComp("hptH3L_pichi2primSB", f1, f3, c,pdf,"#pi chi2primary", "p",4,"Sig_SB","RT_SB", "dp#pi");
  projAndComp("hptH3L_dDcaSB", f1, f3, c,pdf,"d DCA", "p",1,"Sig_SB","RT_SB","dp#pi"  );
  projAndComp("hptH3L_piDcaSB", f1, f3, c,pdf,"#pi DCA", "p",4,"Sig_SB","RT_SB", "dp#pi" );
  projAndComp("hptH3L_pDcaSB", f1, f3, c,pdf,"#p DCA", "p",4,"Sig_SB","RT_SB", "dp#pi" );
  projAndComp("hptH3L_dpDcaSB", f1, f3, c,pdf,"dp pair DCA", "p",4,"Sig_SB","RT_SB", "dp#pi" );
  projAndComp("hptppimassSB", f1, f3, c,pdf,"p#pi Mass (GeV/c^{2})", "plhist",1 ,"Sig" ,"BK","p#pi");
  // projAndComp("hptppichi2ndfSB", f1, f2, c,pdf,"p#pi #chi^{2}_{NDF}","" ,1 , "Sig" ,"BK","p#pi");
  // projAndComp("hptppichi2primSB", f1, f2, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "Sig" ,"BK","p#pi");
  // projAndComp("hptppilSB", f1, f2, c,pdf,"(p#pi) l" ,"",1 , "Sig" ,"BK","p#pi");
  // projAndComp("hptppildlSB", f1, f2, c,pdf,"(p#pi) l/#deltal" ,"",1 , "Sig" ,"BK","p#pi");

  projAndScaleComp("hptH3L_l", scale,f1, f2, f2, c,pdf,"l", "p",1,"Sig-ME","ME", "dp#pi");
  projAndScaleComp("hptppimass", scale, f1, f2, f2, c,pdf,"M(p#pi)", "p",1,"Sig-ME","ME", "p#pi");
  // projAndScaleCompNoRatio("hptppimass", scale, f1, f2, f2, c,pdf,"M(p#pi)", "p",1,"Sig-ME","ME", "p#pi");
  projAndScaleComp("hptH3L_ldl", scale,f1, f2, f2, c,pdf,"l/dl", "p",1,"Sig-ME","ME", "dp#pi");
  projAndScaleComp("hptH3L_chi2ndf", scale,f1, f2,f2, c,pdf,"chi2NDF", "p",1,"Sig-ME","ME", "dp#pi");
  projAndScaleComp("hptH3L_chi2topo", scale,f1, f2, f2,c,pdf,"chi2topo", "p",1,"Sig-ME","ME", "dp#pi" );
  projAndScaleComp("hptH3L_dchi2prim", scale,f1, f2, f2,c,pdf,"d chi2primary", "p",1,"Sig-ME","ME", "dp#pi");
  projAndScaleComp("hptH3L_pchi2prim", scale,f1, f2, f2,c,pdf,"p chi2primary", "p",1,"Sig-ME","ME", "dp#pi");
  projAndScaleComp("hptH3L_pichi2prim", scale,f1, f2, f2,c,pdf,"#pi chi2primary", "p",1,"Sig-ME","ME", "dp#pi");
  projAndScaleComp("hptH3L_dDca", scale,f1, f2, f2,c,pdf,"d DCA", "p",1,"Sig-ME","ME","dp#pi"  );
  projAndScaleComp("hptH3L_piDca", scale,f1, f2, f2,c,pdf,"#pi DCA", "p",1,"Sig-ME","ME", "dp#pi" );
  projAndScaleComp("hptH3L_pDca", scale,f1, f2, f2,c,pdf,"#p DCA", "p",1,"Sig-ME","ME", "dp#pi" );
  projAndScaleComp("hptH3L_dpDca", scale,f1, f2, f2,c,pdf,"dp pair DCA", "p",1,"Sig-ME","ME", "dp#pi" );
  projAndScaleComp("hptppichi2prim",scale, f1, f2,f2, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "Sig-ME" ,"ME","p#pi");

  // projAndScaleComp("hptH3L_l", scale_rt,f1, f3, f3, c,pdf,"l", "p",1,"Sig-RT","RT", "dp#pi");
  // projAndScaleComp("hptppimass", scale_rt, f1, f3, f3, c,pdf,"M(p#pi)", "p",1,"Sig-RT","RT", "p#pi");
  // projAndScaleComp("hptH3L_ldl", scale_rt,f1, f3, f3, c,pdf,"l/dl", "p",1,"Sig-RT","RT", "dp#pi");
  // projAndScaleComp("hptH3L_chi2ndf", scale_rt,f1, f3,f3, c,pdf,"chi2NDF", "p",1,"Sig-RT","RT", "dp#pi");
  // projAndScaleComp("hptH3L_chi2topo", scale_rt,f1, f3, f3,c,pdf,"chi2topo", "p",1,"Sig-RT","RT", "dp#pi" );
  // projAndScaleComp("hptH3L_dchi2prim", scale_rt,f1, f3, f3,c,pdf,"d chi2primary", "p",1,"Sig-RT","RT", "dp#pi");
  // projAndScaleComp("hptH3L_pchi2prim", scale_rt,f1, f3, f3,c,pdf,"p chi2primary", "p",1,"Sig-RT","RT", "dp#pi");
  // projAndScaleComp("hptH3L_pichi2prim", scale_rt,f1, f3, f3,c,pdf,"#pi chi2primary", "p",1,"Sig-RT","RT", "dp#pi");
  // projAndScaleComp("hptH3L_dDca", scale_rt,f1, f3, f3,c,pdf,"d DCA", "p",1,"Sig-RT","RT","dp#pi"  );
  // projAndScaleComp("hptH3L_piDca", scale_rt,f1, f3, f3,c,pdf,"#pi DCA", "p",1,"Sig-RT","RT", "dp#pi" );
  // projAndScaleComp("hptH3L_pDca", scale_rt,f1, f3, f3,c,pdf,"#p DCA", "p",1,"Sig-RT","RT", "dp#pi" );
  // projAndScaleComp("hptH3L_dpDca", scale_rt,f1, f3, f3,c,pdf,"dp pair DCA", "p",1,"Sig-RT","RT", "dp#pi" );
  // projAndScaleComp("hptppichi2prim",scale_rt, f1, f3,f3, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "Sig-RT" ,"RT","p#pi");

  TFile* fMc = TFile::Open("fout_H3L_MC.root");
  projAndScaleComp("hptH3L_l", scale,f1, f2, fMc, c,pdf,"l", "p",4,"Sig-ME","MC", "dp#pi");
  projAndScaleComp("hptppimass", scale, f1, f2, fMc, c,pdf,"M(p#pi)", "p",4,"Sig-ME","MC", "p#pi");
  projAndScaleCompNoRatio("hptppimass", scale, f1, f2, fMc, c,pdf,"M(p#pi)", "p",4,"Sig-ME","MC", "p#pi");
  projAndScaleComp("hptH3L_ldl", scale,f1, f2, fMc, c,pdf,"l/dl", "p",4,"Sig-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_chi2ndf", scale,f1, f2,fMc, c,pdf,"chi2NDF","p",2,"Sig-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_chi2topo", scale,f1, f2, fMc,c,pdf,"chi2topo","p",2,"Sig-ME","MC", "dp#pi" );
  // projAndScaleComp("hptH3L_chi2", scale,f1, f2, fMc,c,pdf,"chi2topo","p",2,"Sig-ME","MC", "dp#pi" );
  projAndScaleComp("hptH3L_dchi2prim", scale,f1, f2, fMc,c,pdf,"d chi2primary", "p",5,"Sig-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_pchi2prim", scale,f1, f2, fMc,c,pdf,"p chi2primary", "p",2,"Sig-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_pichi2prim", scale,f1, f2, fMc,c,pdf,"#pi chi2primary", "p",2,"Sig-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_dDca", scale,f1, f2, fMc,c,pdf,"d DCA", "p",2,"Sig-ME","MC","dp#pi"  );
  projAndScaleComp("hptH3L_piDca", scale,f1, f2, fMc,c,pdf,"#pi DCA", "p",2,"Sig-ME","MC", "dp#pi" );
  projAndScaleComp("hptH3L_pDca", scale,f1, f2, fMc,c,pdf,"#p DCA", "p",2,"Sig-ME","MC", "dp#pi" );
  projAndScaleComp("hptH3L_dpDca", scale,f1, f2, fMc,c,pdf,"dp pair DCA", "p",2,"Sig-ME","MC", "dp#pi" );
  projAndScaleComp("hptppichi2prim",scale, f1, f2,fMc, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",2 , "Sig-ME" ,"MC","p#pi");


  // TFile* fMc_ph = TFile::Open("fout_H3L_phase_MC.root");
  // //  projAndScaleComp("hptH3L_l", scale,f1, f2, fMc_ph, c,pdf,"l", "p",4,"Sig-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptppimass", scale, f1, f2, fMc_ph, c,pdf,"M(p#pi)", "p",4,"Sig-ME","phase_MC", "p#pi");
  // projAndScaleComp("hptH3L_ldl", scale,f1, f2, fMc_ph, c,pdf,"l/dl", "p",4,"Sig-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptH3L_chi2ndf", scale,f1, f2,fMc_ph, c,pdf,"chi2NDF","p",2,"Sig-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptH3L_chi2topo", scale,f1, f2, fMc_ph,c,pdf,"chi2topo","p",2,"Sig-ME","phase_MC", "dp#pi" );
  // // projAndScaleComp("hptH3L_chi2", scale,f1, f2, fMc_ph,c,pdf,"chi2topo","p",2,"Sig-ME","phase_MC", "dp#pi" );
  // projAndScaleComp("hptH3L_dchi2prim", scale,f1, f2, fMc_ph,c,pdf,"d chi2primary", "p",2,"Sig-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptH3L_pchi2prim", scale,f1, f2, fMc_ph,c,pdf,"p chi2primary", "p",2,"Sig-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptH3L_pichi2prim", scale,f1, f2, fMc_ph,c,pdf,"#pi chi2primary", "p",2,"Sig-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptH3L_dDca", scale,f1, f2, fMc_ph,c,pdf,"d DCA", "p",2,"Sig-ME","phase_MC","dp#pi"  );
  // projAndScaleComp("hptH3L_piDca", scale,f1, f2, fMc_ph,c,pdf,"#pi DCA", "p",2,"Sig-ME","phase_MC", "dp#pi" );
  // projAndScaleComp("hptH3L_pDca", scale,f1, f2, fMc_ph,c,pdf,"#p DCA", "p",2,"Sig-ME","phase_MC", "dp#pi" );
  // projAndScaleComp("hptH3L_dpDca", scale,f1, f2, fMc_ph,c,pdf,"dp pair DCA", "p",2,"Sig-ME","phase_MC", "dp#pi" );
  // projAndScaleComp("hptppichi2prim",scale, f1, f2,fMc_ph, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",2 , "Sig-ME" ,"phase_MC","p#pi");

  TFile* fMc_ld = TFile::Open("fout_Lambda_MC_Cuts.root");
  projAndScaleComp("hptH3L_l", scale,f1, f2, fMc_ld, c,pdf,"l", "p",4,"Sig-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptppimass", scale, f1, f2, fMc_ld, c,pdf,"M(p#pi)", "p",4,"Sig-ME","MC_#Lambda", "p#pi",0,1);
  projAndScaleCompNoRatio("hptppimass", scale, f1, f2, fMc_ld, c,pdf,"M(p#pi)", "p",4,"Sig-ME","MC_#Lambda", "p#pi",0);
  projAndScaleComp("hptH3L_ldl", scale,f1, f2, fMc_ld, c,pdf,"l/dl", "p",4,"Sig-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_chi2ndf", scale,f1, f2,fMc_ld, c,pdf,"chi2NDF","p",2,"Sig-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_chi2topo", scale,f1, f2, fMc_ld,c,pdf,"chi2topo","p",2,"Sig-ME","MC_#Lambda", "dp#pi" ,0,1);
  // projAndScaleComp("hptH3L_chi2", scale,f1, f2, fMc_ld,c,pdf,"chi2topo","p",2,"Sig-ME","MC_#Lambda", "dp#pi" ,0,1);
  projAndScaleComp("hptH3L_dchi2prim", scale,f1, f2, fMc_ld,c,pdf,"d chi2primary", "p",5,"Sig-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_pchi2prim", scale,f1, f2, fMc_ld,c,pdf,"p chi2primary", "p",2,"Sig-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_pichi2prim", scale,f1, f2, fMc_ld,c,pdf,"#pi chi2primary", "p",2,"Sig-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_dDca", scale,f1, f2, fMc_ld,c,pdf,"d DCA", "p",2,"Sig-ME","MC_#Lambda","dp#pi"  ,0,1);
  projAndScaleComp("hptH3L_piDca", scale,f1, f2, fMc_ld,c,pdf,"#pi DCA", "p",2,"Sig-ME","MC_#Lambda", "dp#pi" ,0,1);
  projAndScaleComp("hptH3L_pDca", scale,f1, f2, fMc_ld,c,pdf,"#p DCA", "p",2,"Sig-ME","MC_#Lambda", "dp#pi" ,0,1);
  projAndScaleComp("hptH3L_dpDca", scale,f1, f2, fMc_ld,c,pdf,"dp pair DCA", "p",2,"Sig-ME","MC_#Lambda", "dp#pi",0,1 );
  projAndScaleComp("hptppichi2prim",scale, f1, f2,fMc_ld, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",2 , "Sig-ME" ,"MC_#Lambda","p#pi", 0, 1);

  cout<<"chi2ndf" << calpurity("hptH3L_chi2ndf",0, 4, lowpt, highpt, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}" )<<endl;
  cout<<"chi2topo" << calpurity("hptH3L_chi2topo",0, 4, lowpt, highpt, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}" )<<endl;
  // cout<<"decaylength" << calpurity("hptH3L_l",8, 19,  lowpt, highpt,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l" )<<endl;
  // cout<<"l/dl" << calpurity("hptH3L_ldl",5, 19,  lowpt, highpt,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l/dl" )<<endl;
  // cout<<"dchi2prim" << calpurity("hptH3L_dchi2prim",0, 19,  lowpt, highpt,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L d #chi^{2}" )<<endl;
  // cout<<"pichi2prim" << calpurity("hptH3L_pichi2prim",0, 19,  lowpt, highpt,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L pi #chi^{2}" )<<endl;
  // cout<<"pip" << calpurity("hptppil",0, 39,  lowpt, highpt,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L pi #chi^{2}" )<<endl;
  pdf->On();
  pdf->Close();

}
void drawSEData()
{
  SetsPhenixStyle();
  TCanvas* c = new TCanvas("c","c", 800,800);
  TPDF* pdf = new TPDF("SameEventQA_check.pdf");
  pdf->Off();

  TFile *f1 = TFile::Open("fout_H3L_data_KF_Jul25.root"); 
  // TH2F* h2sig = (TH2F*)f1->Get("hptH3Lmass")->Clone("hptH3Lmass_sig");
  TH2F* h2sig = (TH2F*)f1->Get("hptH3Lmass")->Clone("hptH3Lmass_sig");
  TH1F* hsig = (TH1F*)h2sig->ProjectionY("hsig");

  // TFile *f2 = TFile::Open("fout_H3L_data_ME.root"); 
  TFile *f2 = TFile::Open("fout_H3L_data_SE_Jul25.root"); 
  TH2F* h2bk = (TH2F*)f2->Get("hptH3Lmass")->Clone("hptH3Lmass_ME");
  TH1F* hbk = (TH1F*)h2bk->ProjectionY("hbk");

  double sig_sb = hsig->Integral(hsig->GetXaxis()->FindBin(2.97), hsig->GetXaxis()->FindBin(2.98)) + hsig->Integral(hsig->GetXaxis()->FindBin(3.0),  hsig->GetXaxis()->FindBin(3.02));
  // double sig_sb = hsig->Integral(hsig->GetXaxis()->FindBin(2.95), hsig->GetXaxis()->FindBin(2.98)) ;
  double bk_sb =  hbk->Integral(hbk->GetXaxis()->FindBin(2.97), hbk->GetXaxis()->FindBin(2.98)) + hbk->Integral(hbk->GetXaxis()->FindBin(3.0),  hbk->GetXaxis()->FindBin(3.02));
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
  // TRatioPlot* rp = new TRatioPlot(hsig,hbk);
  // rp->Draw();
  TLegend* leg = new TLegend( 0.2 , 0.7 ,0.5,0.9  );
  leg->AddEntry(hsig, "Official KF","pl");
  leg->AddEntry(hbk, "My Contruction","pl");
  leg->Draw();
  drawLatex( 0.2,0.6,Form("My/KF=%0.4f", 1./scale), 0.055);
  addpdf(pdf);

  TH1F* hratio = (TH1F*)hbk->Clone("ratio");
  hratio->Divide(hsig);
  hratio->Draw();
  hratio->GetYaxis()->SetTitle("My/KF");
  drawLine(2.95,1,3.05,1,1.5, 2, 1);
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
void drawRT_check()
{
  // double highpt = 4, lowpt = 0, lowy=-1.5, highy = 0;
  double highpt = 2.5, lowpt = 1, lowy=-0.9, highy = -0.1;
  // TString histname="hH3LMassPtY";
  TString histname="hH3LMassPtY_5_40";
  SetsPhenixStyle();
  TCanvas* c = new TCanvas("c","c");
  TPDF* pdf = new TPDF("Rotation_QA_check.pdf");
  pdf->Off();

  // TFile *f1 = TFile::Open("fout_H3L_data_SE_large.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_large.root"); 
  // TFile *f1 = TFile::Open("out_KF_test/fout_H3L_KF_test.root"); 
  TFile *f1 = TFile::Open("fout_H3L_data_KF_Jul25.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_KF_Jul22_part.root"); 
  // TH2F* h2sig = (TH2F*)f1->Get("hptH3Lmass")->Clone("hptH3Lmass_sig");
  // h2sig->SetDirectory(0);
  // TH1F* hsig = (TH1F*)h2sig->ProjectionY("hsig");
  // TH3F* h2sig = (TH3F*)f1->Get("hH3LMassPtY_5_40")->Clone("hptH3Lmass_sig");
  TH3F* h2sig = (TH3F*)f1->Get(histname.Data())->Clone("hptH3Lmass_sig");
  h2sig->SetDirectory(0);
  // TH3F* h2sig->Project3D("xz");
  TH1F* hsig = (TH1F*)h2sig->ProjectionY("hsig", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2sig->GetZaxis()->FindBin(lowy), h2sig->GetZaxis()->FindBin(highy));
  hsig->SetDirectory(0);
  TH1F* hcent_se = (TH1F*)f1->Get("hcent")->Clone("hcent_se");
  double nEvents_se = hcent_se->Integral();

  // TFile *f2 = TFile::Open("fout_H3L_data_ME_large.root"); 
  TFile *f2 = TFile::Open("fout_H3L_ME_Jul25.root"); 
  // TFile *f2 = TFile::Open("fout_H3L_data_SE.root"); 
  // TH2F* h2bk = (TH2F*)f2->Get("hptH3Lmass")->Clone("hptH3Lmass_ME");
  // TH3F* h2bk = (TH3F*)f2->Get("hH3LMassPtY_5_40")->Clone("hptH3Lmass_bk");
  TH3F* h2bk = (TH3F*)f2->Get(histname.Data())->Clone("hptH3Lmass_bk");
  h2bk->SetDirectory(0);
  TH1F* hbk = (TH1F*)h2bk->ProjectionY("hbk", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2bk->GetZaxis()->FindBin(lowy), h2bk->GetZaxis()->FindBin(highy));
  hbk->SetDirectory(0);

  TFile* f3 = TFile::Open("fout_H3L_data_RT_Jul25.root");
  // TFile* f3 = TFile::Open("fout_H3L_data_RT_Jul22_full.root");
  // TFile* f3 = TFile::Open("out_RT_test/fout_H3L_RT_test.root");
  // TH1F* hcent_rt= (TH1F*)f3->Get("hcent");
  // double nEvents_rt = hcent_rt->Integral(); 
  // TH2F* h2rt = (TH2F*)f3->Get("hptH3Lmass")->Clone("hptH3Lmass_RT");
  // h2rt->SetDirectory(0);
  // TH1F* hrt = (TH1F*)h2rt->ProjectionY("hrt");
  TH3F* h2rt = (TH3F*)f3->Get(histname.Data())->Clone("hptH3Lmass_rt");
  h2rt->SetDirectory(0);
  TH1F* hrt = (TH1F*)h2rt->ProjectionY("hrt", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2rt->GetZaxis()->FindBin(lowy), h2rt->GetZaxis()->FindBin(highy));
  hrt->SetDirectory(0);
  // hrt->Scale(nEvents_se/(1.*nEvents_rt) );
  // cout << nEvents_se/(1.*nEvents_rt)<<endl;
  //
  //scale
  // double sig_sb = hsig->Integral(hsig->GetXaxis()->FindBin(2.96), hsig->GetXaxis()->FindBin(2.98)) ;
  // double sig_sb = hsig->Integral(hsig->GetXaxis()->FindBin(2.97), hsig->GetXaxis()->FindBin(2.98)) + hsig->Integral(hsig->GetXaxis()->FindBin(3.0),  hsig->GetXaxis()->FindBin(3.02));
  double sig_sb = hsig->Integral(hsig->GetXaxis()->FindBin(3.0),  hsig->GetXaxis()->FindBin(3.02));
  // double rt_sb =  hrt->Integral(hrt->GetXaxis()->FindBin(2.97), hrt->GetXaxis()->FindBin(2.98)) + hrt->Integral(hrt->GetXaxis()->FindBin(3.0),  hrt->GetXaxis()->FindBin(3.02));
  double rt_sb =   hrt->Integral(hrt->GetXaxis()->FindBin(3.0),  hrt->GetXaxis()->FindBin(3.02));
  // double rt_sb =  hrt->Integral(hrt->GetXaxis()->FindBin(2.96), hrt->GetXaxis()->FindBin(2.98)) ;
  double scale_rt = sig_sb/rt_sb;
  // sig_sb = hsig->Integral(hsig->GetXaxis()->FindBin(2.97), hsig->GetXaxis()->FindBin(2.98)) + hsig->Integral(hsig->GetXaxis()->FindBin(3.0),  hsig->GetXaxis()->FindBin(3.02));
  sig_sb = hsig->Integral(hsig->GetXaxis()->FindBin(3.0),  hsig->GetXaxis()->FindBin(3.02));
  // double bk_sb =  hbk->Integral(hbk->GetXaxis()->FindBin(2.97), hbk->GetXaxis()->FindBin(2.98)) + hbk->Integral(hbk->GetXaxis()->FindBin(3.0),  hbk->GetXaxis()->FindBin(3.02));
  double bk_sb =  hbk->Integral(hbk->GetXaxis()->FindBin(3.0),  hbk->GetXaxis()->FindBin(3.02));
  // double rt_sb =  hrt->Integral(hrt->GetXaxis()->FindBin(2.97), hrt->GetXaxis()->FindBin(2.98)) + hrt->Integral(hrt->GetXaxis()->FindBin(3.01),  hrt->GetXaxis()->FindBin(3.02));
  // double bk_sb =  hbk->Integral(hbk->GetXaxis()->FindBin(2.96), hbk->GetXaxis()->FindBin(2.98)) ;
  double scale = sig_sb/bk_sb;
  cout<<"ME scale: " <<1./scale << endl;
  hbk->Scale(scale);
  hrt->Scale(scale_rt);
  cout <<"rotation: " <<sig_sb/rt_sb<< endl;

  hsig->Draw();
  hsig->GetXaxis()->SetTitle("Mass(p#pid) (GeV/c^{2})");
  hsig->GetYaxis()->SetTitle("Counts");
  hsig->GetYaxis()->SetRangeUser(-0.1*hsig->GetMaximum(), hsig->GetMaximum()*1.1);
  
  setHistStyle(hbk, kRed, kOpenCircle, 1.5);
  hbk->Draw("same");
  setHistStyle(hrt, kGreen+2, kDiamond, 1.5);
  hrt->Draw("same");
  addpdf(pdf);
  cout<<"binwidth: "<< hrt->GetBinWidth(1)<< endl;

  c->cd();
  // calculate the significance
  TH1F* hsig_bk = (TH1F*)hsig->Clone("hsig_bk");
  hsig_bk->Add(hbk,-1);
  setHistStyle(hsig_bk, kBlue, kFullCircle, 1.5);
  TH1F* hsig_rt = (TH1F*)hsig->Clone("hsig_rt");
  hsig_rt->Add(hrt,-1);
  setHistStyle(hsig_rt, kGreen+2, kDiamond, 1.5);
  hsig_rt->Rebin();
  // hsig_bk->Rebin();
  // hsig_bk->Draw();
  hsig_rt->Draw();

  // TF1* fit = new TF1("fit" ,"[0]*TMath::Gaus(x, [1], [2], 1)+[3]*x+[4]", 2.97,3.02 );
  TF1* fit = new TF1("fit" ,"gaus(0)+pol1(3)", 2.97,3.02 );
  TF1* resfit = new TF1("resfit" ,"pol1", 2.95,3.05 );
  hsig_bk->GetXaxis()->SetRangeUser(2.97,2.985);
  hsig_bk->Fit(resfit,"R");
  fit->SetLineColor(kRed);
  double yield_bc = hsig_bk->Integral(hsig_bk->GetXaxis()->FindBin(2.98), hsig_bk->GetXaxis()->FindBin(3.));
  double para[5]={yield_bc*hsig_bk->GetBinWidth(1)/sqrt(2*3.1415), 2.991, 0.0015,  resfit->GetParameter(0), resfit->GetParameter(1)};
  fit->SetParameters(para);
  // double lowx=2.98 ,highx =3.01;
  double lowx=2.97 ,highx =3.02;
  hsig_bk->GetXaxis()->SetRangeUser(lowx,highx);
  hsig_rt->GetXaxis()->SetRangeUser(lowx,highy);
  // hsig_bk->Draw("same");
  hsig_rt->Draw();
  // hsig_bk->Fit(fit,"R");
  resfit->SetParameter(0, fit->GetParameter(3));
  resfit->SetParameter(1, fit->GetParameter(4));
  // resfit->Draw("same");
  setHistStyle(resfit, kRed-2, 9, 2.5 ,1);
  // hsig_bk->Rebin(4);
  drawLine(lowx, 0, highx, 0, 1.5, 2, 1 );
  // cout<<"binwidth: "<< hsig_bk->GetBinWidth(1)<< endl;
  TF1* fit_rt = (TF1*)fit->Clone("fit_rt");
  fit_rt->SetParameters(para);
  TF1* resfit_rt = (TF1*)resfit->Clone("resfit_rt");
  setHistStyle(resfit_rt, kGreen+2, 9, 2.5 ,1);
  setHistStyle(fit_rt, kGreen+2, 9, 2.5 ,1);
  hsig_rt->GetXaxis()->SetRangeUser(lowx,highx);
  hsig_rt->Fit(resfit_rt,"R");
  fit_rt->SetParameter(3, resfit_rt->GetParameter(0) );
  fit_rt->SetParameter(4, resfit_rt->GetParameter(1) );
  hsig_rt->Fit(fit_rt,"R");
  resfit_rt->SetParameter(0, fit_rt->GetParameter(3));
  resfit_rt->SetParameter(1, fit_rt->GetParameter(4));
  resfit_rt->Draw("same");

  double sigma = fit_rt->GetParameter(2);
  double mean = fit_rt->GetParameter(1);
  double yield_me = fit->GetParameter(0)/hsig_bk->GetBinWidth(1)*fit->GetParameter(2)*sqrt(2*3.1415);
  double yield_rt = fit_rt->GetParameter(0)/hsig_rt->GetBinWidth(1)*fit_rt->GetParameter(2)*sqrt(2*3.1415);
  double yield_counts = hsig_bk->Integral(hsig_bk->GetXaxis()->FindBin(mean-2.5*sigma), hsig_bk->GetXaxis()->FindBin(mean+2.5*sigma));
  cout<<"rotate: " << yield_rt<<" ME: "<<yield_me<<" bin counting: "<<yield_counts << endl;
  double bk_counts = hbk->Integral(hbk->GetXaxis()->FindBin(mean-2.5*sigma), hbk->GetXaxis()->FindBin(mean+2.5*sigma));
  double sp_counts = hsig->Integral(hbk->GetXaxis()->FindBin(mean-2.5*sigma), hbk->GetXaxis()->FindBin(mean+2.5*sigma));
  // double significance = yield_counts/sqrt(yield_counts+bk_counts);
  // double significance = yield_me/sqrt(sp_counts);
  double significance = yield_rt/sqrt(sp_counts);
  cout<<"significance: " <<significance << endl;

  TLegend* leg = new TLegend( 0.55, 0.35 ,0.9,0.7  );
  // leg->AddEntry(hbk, "mix-event(ME)","pl");
  // leg->AddEntry(hrt, "rotate d(RT) (scale)","pl");
  leg->AddEntry(hsig, "same-event(SE)","pl");
  // leg->AddEntry(hsig_bk, "SE-ME","pl");
  leg->AddEntry(hsig_rt, "SE-RT","pl");
  leg->Draw();
  // drawLatex( 0.2,0.82,Form("ME/SE=%0.2f", 1./scale_rt), 0.055);
  drawLatex( 0.2,0.75,Form("Yield=%0.2f", yield_rt), 0.055);
  drawLatex( 0.2,0.68,Form("sigma=%0.2f MeV", sigma*1000.), 0.055);
  drawLatex( 0.2,0.61,Form("nEvents=%0.0f M", nEvents_se/1e6), 0.055);

  addpdf(pdf);
  return;

  hsig_bk->Draw();
  // hsig_rt->Draw("same");
  TLegend* leg2 = new TLegend( 0.2 , 0.7 ,0.4,0.9  );
  // leg2->AddEntry(hsig_rt, "SE-RT","pl");
  leg2->AddEntry(hsig_bk, "SE-ME","pl");
  leg2->Draw();

  addpdf(pdf);

  // check the difference between rotation and mix event
  c->Clear(); 
  c->Divide(1,2); 
  c->cd(1);
  TH1F* hrtcn = (TH1F*)hrt->Clone("hrtcn");
  TH1F* hbkcn = (TH1F*)hbk->Clone("hbkcn");
  TH1F* hsigcn = (TH1F*)hsig->Clone("hsigcn");
  hsigcn->Rebin(4);
  hrtcn->Rebin(4);
  hbkcn->Rebin(4);
  hrtcn->Divide(hsigcn);
  hbkcn->Divide(hsigcn);
  hrtcn->Draw();
  hbkcn->Draw("same");
  hrtcn->GetYaxis()->SetTitleOffset(0.8);
  hrtcn->GetYaxis()->SetTitle("BK/Sig");
  hrtcn->GetXaxis()->SetTitle("Mass(dp#pi) GeV/c^{2}");
  hrtcn->GetYaxis()->SetRangeUser(0.8,1.2);
  TLegend* leg_r = new TLegend( 0.7, 0.3, 0.9, 0.5); 
  leg_r->AddEntry( hbkcn, "MixEvent", "pl");
  leg_r->AddEntry(hrtcn, "Rotation", "pl");
  leg_r->Draw();
  // drawLine(2.987,0.8, 2.987, 1.2, 1.5, 2, 1  );
  // drawLine(2.997,0.8, 2.997, 1.2, 1.5, 2, 1  );
  drawLine(mean+2.5*sigma,0.8, mean+2.5*sigma, 1.2, 1.5, 2, 1  );
  drawLine(mean-2.5*sigma,0.8, mean-2.5*sigma, 1.2, 1.5, 2, 1  );
  drawLine(2.95,1, 3.05, 1, 1.5, 9, 1  );
  // addpdf(pdf);

  c->cd(2);
  TH1F* hrtcn2 = (TH1F*)hrt->Clone("hrtcn2");
  TH1F* hbkcn2 = (TH1F*)hbk->Clone("hbkcn2");
  hrtcn2->Rebin(4);
  hbkcn2->Rebin(4);
  hrtcn2->Divide(hbkcn2);
  hrtcn2->Draw();
  hrtcn2->GetYaxis()->SetRangeUser(0.8, 1.2);
  hrtcn2->GetYaxis()->SetTitleOffset(0.8);
  hrtcn2->GetYaxis()->SetTitle("RT/ME");
  hrtcn2->GetXaxis()->SetTitle("Mass(dp#pi) GeV/c^{2}");
  drawLine(2.989,0.8, 2.989, 1.2, 1.5, 2, 1  );
  drawLine(2.996,0.8, 2.996, 1.2, 1.5, 2, 1  );
  drawLine(2.95,1, 3.05, 1, 1.5, 9, 1  );
  c->cd();
  addpdf(pdf);
  projAndComp("hptH3L_lSB", f1, f3, c,pdf,"l", "p",4,"Sig_SB","RT_SB", "dp#pi");
  projAndComp("hptH3L_ldlSB", f1, f3, c,pdf,"l/dl", "p",4,"Sig_SB","RT_SB", "dp#pi");
  projAndComp("hptH3L_chi2ndfSB", f1, f3, c,pdf,"chi2NDF", "p",4,"Sig_SB","RT_SB", "dp#pi");
  projAndComp("hptH3L_chi2topoSB", f1, f3, c,pdf,"chi2topo", "p",4,"Sig_SB","RT_SB", "dp#pi" );
  projAndComp("hptH3L_dchi2primSB", f1, f3, c,pdf,"d chi2primary", "p",4,"Sig_SB","RT_SB", "dp#pi");
  projAndComp("hptH3L_pchi2primSB", f1, f3, c,pdf,"p chi2primary", "p",4,"Sig_SB","RT_SB", "dp#pi");
  projAndComp("hptH3L_pichi2primSB", f1, f3, c,pdf,"#pi chi2primary", "p",4,"Sig_SB","RT_SB", "dp#pi");
  projAndComp("hptH3L_dDcaSB", f1, f3, c,pdf,"d DCA", "p",1,"Sig_SB","RT_SB","dp#pi"  );
  projAndComp("hptH3L_piDcaSB", f1, f3, c,pdf,"#pi DCA", "p",4,"Sig_SB","RT_SB", "dp#pi" );
  projAndComp("hptH3L_pDcaSB", f1, f3, c,pdf,"#p DCA", "p",4,"Sig_SB","RT_SB", "dp#pi" );
  projAndComp("hptH3L_dpDcaSB", f1, f3, c,pdf,"dp pair DCA", "p",4,"Sig_SB","RT_SB", "dp#pi" );
  //
  // projAndComp("hptppimassSB", f1, f2, c,pdf,"p#pi Mass (GeV/c^{2})", "plhist",1 ,"Sig" ,"BK","p#pi");
  // projAndComp("hptppichi2ndfSB", f1, f2, c,pdf,"p#pi #chi^{2}_{NDF}","" ,1 , "Sig" ,"BK","p#pi");
  // projAndComp("hptppichi2primSB", f1, f2, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "Sig" ,"BK","p#pi");
  // projAndComp("hptppilSB", f1, f2, c,pdf,"(p#pi) l" ,"",1 , "Sig" ,"BK","p#pi");
  // projAndComp("hptppildlSB", f1, f2, c,pdf,"(p#pi) l/#deltal" ,"",1 , "Sig" ,"BK","p#pi");

  // projAndScaleComp("hptH3L_l", scale,f1, f2, f2, c,pdf,"l", "p",1,"Sig-ME","ME", "dp#pi");
  // projAndScaleComp("hptppimass", scale, f1, f2, f2, c,pdf,"M(p#pi)", "p",1,"Sig-ME","ME", "p#pi");
  // projAndScaleComp("hptH3L_ldl", scale,f1, f2, f2, c,pdf,"l/dl", "p",1,"Sig-ME","ME", "dp#pi");
  // projAndScaleComp("hptH3L_chi2ndf", scale,f1, f2,f2, c,pdf,"chi2NDF", "p",1,"Sig-ME","ME", "dp#pi");
  // projAndScaleComp("hptH3L_chi2topo", scale,f1, f2, f2,c,pdf,"chi2topo", "p",1,"Sig-ME","ME", "dp#pi" );
  // projAndScaleComp("hptH3L_dchi2prim", scale,f1, f2, f2,c,pdf,"d chi2primary", "p",1,"Sig-ME","ME", "dp#pi");
  // projAndScaleComp("hptH3L_pchi2prim", scale,f1, f2, f2,c,pdf,"p chi2primary", "p",1,"Sig-ME","ME", "dp#pi");
  // projAndScaleComp("hptH3L_pichi2prim", scale,f1, f2, f2,c,pdf,"#pi chi2primary", "p",1,"Sig-ME","ME", "dp#pi");
  // projAndScaleComp("hptH3L_dDca", scale,f1, f2, f2,c,pdf,"d DCA", "p",1,"Sig-ME","ME","dp#pi"  );
  // projAndScaleComp("hptH3L_piDca", scale,f1, f2, f2,c,pdf,"#pi DCA", "p",1,"Sig-ME","ME", "dp#pi" );
  // projAndScaleComp("hptH3L_pDca", scale,f1, f2, f2,c,pdf,"#p DCA", "p",1,"Sig-ME","ME", "dp#pi" );
  // projAndScaleComp("hptH3L_dpDca", scale,f1, f2, f2,c,pdf,"dp pair DCA", "p",1,"Sig-ME","ME", "dp#pi" );
  // projAndScaleComp("hptppichi2prim",scale, f1, f2,f2, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "Sig-ME" ,"ME","p#pi");

  projAndScaleComp("hptH3L_l", scale_rt,f1, f3, f3, c,pdf,"l", "p",1,"Sig-RT","RT", "dp#pi");
  projAndScaleComp("hptppimass", scale_rt, f1, f3, f3, c,pdf,"M(p#pi)", "p",1,"Sig-RT","RT", "p#pi");
  projAndScaleComp("hptH3L_ldl", scale_rt,f1, f3, f3, c,pdf,"l/dl", "p",1,"Sig-RT","RT", "dp#pi");
  projAndScaleComp("hptH3L_chi2ndf", scale_rt,f1, f3,f3, c,pdf,"chi2NDF", "p",1,"Sig-RT","RT", "dp#pi");
  projAndScaleComp("hptH3L_chi2topo", scale_rt,f1, f3, f3,c,pdf,"chi2topo", "p",1,"Sig-RT","RT", "dp#pi" );
  projAndScaleComp("hptH3L_dchi2prim", scale_rt,f1, f3, f3,c,pdf,"d chi2primary", "p",1,"Sig-RT","RT", "dp#pi");
  projAndScaleComp("hptH3L_pchi2prim", scale_rt,f1, f3, f3,c,pdf,"p chi2primary", "p",1,"Sig-RT","RT", "dp#pi");
  projAndScaleComp("hptH3L_pichi2prim", scale_rt,f1, f3, f3,c,pdf,"#pi chi2primary", "p",1,"Sig-RT","RT", "dp#pi");
  projAndScaleComp("hptH3L_dDca", scale_rt,f1, f3, f3,c,pdf,"d DCA", "p",1,"Sig-RT","RT","dp#pi"  );
  projAndScaleComp("hptH3L_piDca", scale_rt,f1, f3, f3,c,pdf,"#pi DCA", "p",1,"Sig-RT","RT", "dp#pi" );
  projAndScaleComp("hptH3L_pDca", scale_rt,f1, f3, f3,c,pdf,"#p DCA", "p",1,"Sig-RT","RT", "dp#pi" );
  projAndScaleComp("hptH3L_dpDca", scale_rt,f1, f3, f3,c,pdf,"dp pair DCA", "p",1,"Sig-RT","RT", "dp#pi" );
  projAndScaleComp("hptppichi2prim",scale_rt, f1, f3,f3, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "Sig-RT" ,"RT","p#pi");

  pdf->On();
  pdf->Close();

}
void drawData()
{
  // drawRT_check();
  drawMixData();
  // drawSEData();
}
