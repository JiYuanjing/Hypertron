#include "sPhenixStyle.h"
// #include "style.h"

TH1F* hchi2ndf_Ld;
TH1F* hchi2ndf_H3L;
Double_t MyFitFun(double* x, double* p)
{
  double a = hchi2ndf_H3L->GetBinContent(hchi2ndf_H3L->GetXaxis()->FindBin(x[0]));
  double b = hchi2ndf_Ld->GetBinContent(hchi2ndf_Ld->GetXaxis()->FindBin(x[0]));
  return p[0]*(a+p[1]*b);
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
double BrErr(double y2b, double y2berr, double y3b, double y3berr)
{
   return sqrt( pow((y2b/(y3b+y2b)/(y3b+y2b))*y3berr, 2) +  pow((y3b/(y3b+y2b)/(y3b+y2b))*y2berr, 2));
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
void projAndComp(TString name, TFile* fH3L, TFile* fLa, TCanvas* c,TPDF* pdf,TString xTitle , TString drawstyle="", int rebin=1,  TString legtitle1="H3L quasi", TString legtitle2="#Lambda+d", TString text="dp#pi", double scale = 0 )
{
  cout <<fH3L->GetName()<<" "<<xTitle.Data() << endl;
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
  double ptedge[bins+1]={1, 1.5,2,2.5, 3.5};
  // c->Divide(1,1);
  // int  const bins=1;
  // double ptedge[bins+1]={0,3.};
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
    // TH1F* h1H3= (TH1F*)h2H3L->ProjectionX(Form("h1H3%d",i), h2H3L->GetYaxis()->FindBin(ptedge[i]),h2H3L->GetYaxis()->FindBin(ptedge[i+1]));
    // TH1F* h1La= (TH1F*)h2La->ProjectionX(Form("h1La%d",i), h2La->GetYaxis()->FindBin(ptedge[i]),h2La->GetYaxis()->FindBin(ptedge[i+1]));

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
    if (scale>0) 
    {
      h1La->Scale(scale);
      Norm(h1La);
      Norm(h1H3);

    }
    else  {
      Norm(h1La);
      Norm(h1H3);
    }
    h1H3->GetXaxis()->SetTitle(xTitle.Data());
    // h1H3->GetYaxis()->SetRangeUser(-0.1 , 1.15);
    // h1H3->GetYaxis()->SetRangeUser(-0.1*h1H3->GetMaximum() , 1.3*h1H3->GetMaximum());
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
    hratio->GetYaxis()->SetRangeUser(0.5,1.5);
    // hratio->GetYaxis()->SetRangeUser(0.9,1.1);
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
void projAndComp2(TString name1, TString name2, TFile* fH3L, TFile* fLa, TCanvas* c,TPDF* pdf,TString xTitle , TString drawstyle="", int rebin=1,  TString legtitle1="H3L quasi", TString legtitle2="#Lambda+d", TString text="dp#pi", double scale = 0 )
{
  cout <<fH3L->GetName()<<" "<<xTitle.Data() << endl;
  fH3L->cd();
  TH2F* h2H3L = (TH2F*)fH3L->Get(name1.Data());
  h2H3L->Sumw2();
  h2H3L->SetDirectory(0);
  fLa->cd();
  TH2F* h2La = (TH2F*)fLa->Get(name2.Data());
  h2La->Sumw2();
  h2La->SetDirectory(0);
  c->Clear();

  c->Divide(1,1);
  int  const bins=1;
  // double ptedge[bins+1]={0,4};
  double ptedge[bins+1]={1,2.5};

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
    if (scale>0) 
    {
      h1La->Scale(scale);
      Norm(h1La);
      Norm(h1H3);

    }
    else  {
      Norm(h1La);
      Norm(h1H3);
    }
    h1H3->GetXaxis()->SetTitle(xTitle.Data());
    // h1H3->GetYaxis()->SetRangeUser(-0.1 , 1.15);
    // h1H3->GetYaxis()->SetRangeUser(-0.1*h1H3->GetMaximum() , 1.3*h1H3->GetMaximum());
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
    // hratio->GetYaxis()->SetRangeUser(0.5,1.5);
    hratio->GetYaxis()->SetRangeUser(0.9,1.1);
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
void projAndComp3(TString name1, TString name2, TString name3, TFile* fH3L, TFile* fLa, TFile* fMc, TCanvas* c,TPDF* pdf,TString xTitle , TString drawstyle="", int rebin=1,  TString legtitle1="H3L quasi", TString legtitle2="#Lambda+d", TString legtitle3="MC",TString text="dp#pi", double scale = 0 )
{
  cout <<fH3L->GetName()<<" "<<xTitle.Data() << endl;
  fH3L->cd();
  TH2F* h2H3L = (TH2F*)fH3L->Get(name1.Data());
  h2H3L->Sumw2();
  h2H3L->SetDirectory(0);
  fLa->cd();
  TH2F* h2La = (TH2F*)fLa->Get(name2.Data());
  h2La->Sumw2();
  h2La->SetDirectory(0);
  TH2F* h2Mc = (TH2F*)fMc->Get(name3.Data());
  h2Mc->Sumw2();
  h2Mc->SetDirectory(0);
  c->Clear();

  c->Divide(1,1);
  int  const bins=1;
  // double ptedge[bins+1]={0,4};
  double ptedge[bins+1]={1,2.5};

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
    TH1F* h1Mc= (TH1F*)h2Mc->ProjectionY(Form("h1Mc%d",i), h2La->GetXaxis()->FindBin(ptedge[i]),h2La->GetXaxis()->FindBin(ptedge[i+1]));
    h1La->SetLineColor(kRed);
    h1La->SetMarkerColor(kRed);
    h1La->SetMarkerStyle(kOpenCircle);
    h1La->SetMarkerSize(1);
    h1Mc->SetLineColor(kMagenta);
    h1Mc->SetMarkerColor(kMagenta);
    h1Mc->SetMarkerStyle(kOpenStar);
    h1Mc->SetMarkerSize(1.4);
    h1H3->SetLineColor(kBlue);
    h1H3->SetMarkerColor(kBlue);
    h1H3->SetMarkerSize(1);
    if (rebin>1) {
      h1H3->Rebin(rebin);
      h1La->Rebin(rebin);
      h1Mc->Rebin(rebin);
    }
    if (scale>0) 
    {
      h1La->Scale(scale);
      Norm(h1La);
      Norm(h1H3);
      Norm(h1Mc);

    }
    else  {
      Norm(h1La);
      Norm(h1H3);
      Norm(h1Mc);
    }
    h1H3->GetXaxis()->SetTitle(xTitle.Data());
    // h1H3->GetYaxis()->SetRangeUser(-0.1 , 1.15);
    // h1H3->GetYaxis()->SetRangeUser(-0.1*h1H3->GetMaximum() , 1.3*h1H3->GetMaximum());
    /* h1H3->DrawCopy(drawstyle.Data()); */
    p11->cd();
    h1H3->DrawCopy("");
    h1La->DrawCopy((drawstyle+"same").Data()); 
    h1Mc->DrawCopy((drawstyle+"same").Data()); 
    // TRatioPlot* rp = new TRatioPlot(h1La, h1H3);
    // gPad->SetTicks(0, 1);    rp->Draw();
    drawLatex( 0.18, 0.88, Form("%0.1f<(%s) p_{T}<%0.1f GeV/c", ptedge[i],text.Data(),ptedge[i+1]),0.055);
    TLegend* l = new TLegend(0.65,0.48,0.9,0.67);
    l->AddEntry(h1H3,legtitle1.Data(),"pl");
    l->AddEntry( h1La, legtitle2.Data(), "pl");
    l->AddEntry(h1Mc,legtitle3.Data(),"pl");
    l->Draw();
    p12->cd();
    TH1F* hratio = (TH1F*)h1H3->Clone("hratio");
    hratio->Divide(h1La);
    hratio->Draw();
    // hratio->GetYaxis()->SetRangeUser(hratio->GetMinimum()<0.95?hratio->GetMinimum():0.95,hratio->GetMaximum()>1.05?hratio->GetMaximum():1.05);
    // hratio->GetYaxis()->SetRangeUser(0.5,1.5);
    hratio->GetYaxis()->SetRangeUser(0.9,1.1);
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
  cout <<fH3L->GetName()<<" "<<xTitle.Data() << endl;
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
  // double ptedge[bins+1]={0,4};
  double ptedge[bins+1]={1,2.5};

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
  cout <<fH3L->GetName()<<" "<<xTitle.Data() << endl;
  fH3L->cd();
  // TH2F* h2H3Lsb = (TH2F*)fH3L->Get((name+"SBR").Data())->Clone("hsb_s");
  // TH2F* h2H3Lbksb = (TH2F*)fH3Lbk->Get((name+"SBR").Data())->Clone("hsb_bk");
  // double scale_auto = h2H3Lsb->Integral()/h2H3Lbksb->Integral();
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
  // if (scalemode==1) sc = scale_auto;
  // h2H3Lbk->Scale(sc);
  h2H3Lbk->Scale(sc);
  h2H3L->Add(h2H3Lbk, -1);
  // cout <<"ME scale: "<<scale<<" h2H3Lbk"<< scale_auto<< endl;

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
  // double ptedge[bins+1]={0,1,1.5,2,2.5};
  c->Divide(1,1);
  int  const bins=1;
  double ptedge[bins+1]={1,2.5};
  // double ptedge[bins+1]={0,4};

  for (int i=0;i<bins;i++)
  {
    c->cd(i+1);
    TPad* p11 = new TPad("upperPad", "upperPad",.005, .4, .995, .995);
    TPad* p12 = new TPad("lowerPad", "lowerPad", .005, .005, .995, .4);
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
  cout <<fH3L->GetName()<<" "<<xTitle.Data() << endl;
  fH3L->cd();
  TH2F* h2H3Lsb = (TH2F*)fH3L->Get((name+"SBR").Data())->Clone("hsb_s");
  TH2F* h2H3Lbksb = (TH2F*)fH3Lbk->Get((name+"SBR").Data())->Clone("hsb_bk");
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
TH1F* reBinHist(double xmid,double xlow, double xhigh, TH1F* h, int n)
{
  int const nMax = 2000;
  double xedge[nMax];
  int nbin=0;
  xedge[0]=xlow;
  double rebin = n;
  while (xedge[nbin]<=xmid)
  {
    nbin++;
    xedge[nbin]=xedge[nbin-1]+h->GetBinWidth(1)*rebin;
  }
  int turningpoint = nbin+1;
  rebin = n*4;
  while (xedge[nbin]<xhigh)
  {
    nbin++;
    xedge[nbin]=xedge[nbin-1]+h->GetBinWidth(1)*rebin;
    // h->SetBinContent(h->FindBin(xedge[nbin]-h->GetBinWidth(1)*1.5),h->GetBinContent(h->FindBin(xedge[nbin]-h->GetBinWidth(1)*1.5))*0.5);
    // h->SetBinContent(h->FindBin(xedge[nbin]-h->GetBinWidth(1)*0.5),h->GetBinContent(h->FindBin(xedge[nbin]-h->GetBinWidth(1)*0.5))*0.5);
  }
  // return (TH1F*)h->Rebin(nbin,h->GetName(),xedge);
  h = (TH1F*)h->Rebin(nbin,h->GetName(),xedge);
  for (int ibin=0;ibin<turningpoint;ibin++)
  {
    h->SetBinContent(ibin,h->GetBinContent(ibin)/(1.0*n));
    h->SetBinError(ibin,h->GetBinError(ibin)/sqrt(1.0*n));
  }
  for (int ibin=turningpoint;ibin<h->GetNbinsX();ibin++)
  {
    h->SetBinContent(ibin,h->GetBinContent(ibin)/(1.0*rebin));
    h->SetBinError(ibin,h->GetBinError(ibin)/sqrt(1.0*rebin));
  }
  return h;
}
double calpurity(TString fithistname, double cutslow, double cutshigh, double lowpt, double highpt,  TFile* f1, TFile* f2, double scale, TFile* fMc, TFile* fMc_ld, TCanvas* c, TPDF* pdf, TString xTitle, int rebin, double Xrange1=-1, double Xrange2=-1)
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
  // TH2F* h2chi2ndf_Ld = (TH2F*)fMc_ld->Get((fithistname+"Sig").Data())->Clone("h2chi2ndf_Ld");
  h2chi2ndf_Ld->SetDirectory(0);
  hchi2ndf_Ld = (TH1F*)h2chi2ndf_Ld->ProjectionY("hchi2ndf_Ld", h2chi2ndf_Ld->GetXaxis()->FindBin(lowpt), h2chi2ndf_Ld->GetXaxis()->FindBin(highpt));
  hchi2ndf_Ld->Rebin(rebin);
  
  TH2F* h2chi2ndf_H3L = (TH2F*)fMc->Get(fithistname.Data())->Clone("h2chi2ndf_H3L");
  h2chi2ndf_H3L->SetDirectory(0);
  hchi2ndf_H3L = (TH1F*)h2chi2ndf_H3L->ProjectionY("hchi2ndf_H3L", h2chi2ndf_H3L->GetXaxis()->FindBin(lowpt), h2chi2ndf_H3L->GetXaxis()->FindBin(highpt));
  if (rebin>0) hchi2ndf_H3L->Rebin(rebin);
 
  Norm(hchi2ndf_H3L);
  Norm(hchi2ndf_Ld);
   
  TH2F* h2chi2ndf_Sig = (TH2F*)f1->Get((fithistname+"Sig").Data())->Clone("h2chi2ndf_Sig");
  h2chi2ndf_Sig->SetDirectory(0);
  // TH2F* h2chi2ndf_SBL = (TH2F*)f1->Get((fithistname+"SBL").Data())->Clone("h2chi2ndf_SBL");
  // h2chi2ndf_SBL->SetDirectory(0);
  // TH2F* h2chi2ndf_SBR = (TH2F*)f1->Get((fithistname+"SBR").Data())->Clone("h2chi2ndf_SBR");
  // h2chi2ndf_SBR->SetDirectory(0);

  TH2F* h2chi2ndf_Bk = (TH2F*)f2->Get((fithistname+"Sig").Data())->Clone("h2chi2ndf_Bk");
  h2chi2ndf_Bk->SetDirectory(0);
  h2chi2ndf_Sig->Add( h2chi2ndf_Bk , -1*scale ); 
  // h2chi2ndf_SBL->Add(h2chi2ndf_SBR);
  // h2chi2ndf_SBL->Scale(scale*h2chi2ndf_Bk->Integral()/h2chi2ndf_SBL->Integral());
  // h2chi2ndf_Sig->Add(h2chi2ndf_SBL, -1);

  TH1F* hchi2ndf_Sig = (TH1F*)h2chi2ndf_Sig->ProjectionY("hchi2ndf_Sig", h2chi2ndf_Sig->GetXaxis()->FindBin(lowpt), h2chi2ndf_Sig->GetXaxis()->FindBin(highpt));
  if (rebin>0) hchi2ndf_Sig->Rebin(rebin);

  setHistStyle(hchi2ndf_Sig, kBlue, kFullCircle, 1.5);
  hchi2ndf_Sig->Draw();
  TH1F* hRatio = (TH1F*)hchi2ndf_Sig->Clone("hratio");
  TF1* fMyFit = new TF1("fMyFit", MyFitFun, 0,15,2);
  fMyFit->SetLineColor(kRed);
  fMyFit->SetParLimits(0, 0, 10000000);
  fMyFit->SetParLimits(1, 0, 10000000);
  if (Xrange1>=0 && Xrange2>=0) hchi2ndf_Sig->GetXaxis()->SetRangeUser( Xrange1, Xrange2);
  hchi2ndf_Sig->Fit(fMyFit,"RB");
  // fMyFit->SetParLimits(1, fMyFit->GetParameter(0)*0.2, fMyFit->GetParameter(0)*1);
  // hchi2ndf_Sig->Fit(fMyFit,"RB");

  TH1F* hLd = (TH1F*)hchi2ndf_Ld->Clone("hLd");
  hLd->Scale(fMyFit->GetParameter(1));
  hLd->Draw("same");
  TH1F* hH3L = (TH1F*)hchi2ndf_H3L->Clone("hH3L");
  hH3L->Scale(fMyFit->GetParameter(0));
  setHistStyle(hLd, kBlack, kOpenStar, 2);
  setHistStyle(hH3L, kRed, kOpenCircle, 1.5);
  hH3L->Draw("same");
  gPad->SetLogy();
  double purity = hH3L->Integral(hH3L->GetXaxis()->FindBin(cutslow+1e-6),hH3L->GetXaxis()->FindBin(cutshigh-1e-6))/fMyFit->Integral(cutslow,cutshigh)*hH3L->GetBinWidth(1);
  drawLatex(0.65,0.6,Form("purity = %0.2f",purity),0.055);
  drawLatex(0.65,0.5,Form("%0.1f<pT<%0.1f", lowpt, highpt),0.055);
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
      hRatio->SetBinError(i, fMyFit->Eval(hRatio->GetBinCenter(i))/pow(hRatio->GetBinContent(i),2 )*hRatio->GetBinError(i));
      hRatio->SetBinContent(i, fMyFit->Eval(hRatio->GetBinCenter(i))/hRatio->GetBinContent(i) );
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
  c->cd();
  addpdf(pdf);
  // delete hchi2ndf_Ld;
  // delete hchi2ndf_H3L;
  return purity;
}
double calpurity(TString fithistname, double cutslow, double cutshigh, double highpt, double lowpt, double highy, double lowy, TFile* f1, TFile* f2, double scale_i, TFile* fMc, TFile* fMc_ld, TCanvas* c, TPDF* pdf, TString xTitle, int rebin, double Xrange1, double Xrange2,double &error, double  rebinedge)
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
  TH3F* h2chi2ndf_Ld = (TH3F*)fMc_ld->Get(fithistname.Data())->Clone("h2chi2ndf_Ld");
  // TH3F* h2chi2ndf_Ld = (TH3F*)fMc_ld->Get((fithistname+"Sig").Data())->Clone("h2chi2ndf_Ld");
  h2chi2ndf_Ld->SetDirectory(0);
  hchi2ndf_Ld = (TH1F*)h2chi2ndf_Ld->ProjectionY("hchi2ndf_Ld", h2chi2ndf_Ld->GetXaxis()->FindBin(lowpt), h2chi2ndf_Ld->GetXaxis()->FindBin(highpt),
                h2chi2ndf_Ld->GetZaxis()->FindBin(lowy), h2chi2ndf_Ld->GetZaxis()->FindBin(highy)  );
  if (rebin>0) {
    hchi2ndf_Ld->Rebin(rebin);
    if (rebinedge>0) hchi2ndf_Ld = (TH1F*)reBinHist( (Xrange1+Xrange2)*rebinedge, Xrange1, Xrange2, hchi2ndf_Ld, 1 ); 
  }
  
  TH3F* h2chi2ndf_H3L = (TH3F*)fMc->Get(fithistname.Data())->Clone("h2chi2ndf_H3L");
  // TH3F* h2chi2ndf_H3L = (TH3F*)fMc->Get((fithistname+"Sig").Data())->Clone("h2chi2ndf_H3L");
  h2chi2ndf_H3L->SetDirectory(0);
  hchi2ndf_H3L = (TH1F*)h2chi2ndf_H3L->ProjectionY("hchi2ndf_H3L", h2chi2ndf_H3L->GetXaxis()->FindBin(lowpt), h2chi2ndf_H3L->GetXaxis()->FindBin(highpt),
                  h2chi2ndf_H3L->GetZaxis()->FindBin(lowy), h2chi2ndf_H3L->GetZaxis()->FindBin(highy));
  if (rebin>0) {
    hchi2ndf_H3L->Rebin(rebin);
    if (rebinedge>0) hchi2ndf_H3L = (TH1F*)reBinHist( (Xrange1+Xrange2)*rebinedge, Xrange1, Xrange2, hchi2ndf_H3L, 1 ); 
  }
 
  Norm(hchi2ndf_H3L);
  Norm(hchi2ndf_Ld);

  TH1F* hcent_se = (TH1F*)f1->Get("hcent")->Clone("hcent_se");
  double nEvents_se = hcent_se->Integral(4, 9);

  TH3F* h2sig = (TH3F*)f1->Get("hH3LMassPtY")->Clone("hptH3Lmass_sig");
  h2sig->SetDirectory(0);
  TH1F* hsig = (TH1F*)h2sig->ProjectionY("hsig", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2sig->GetZaxis()->FindBin(lowy), h2sig->GetZaxis()->FindBin(highy));
  hsig->SetDirectory(0);

  TH3F* h2bk = (TH3F*)f2->Get("hH3LMassPtY")->Clone("hptH3Lmass_bk");
  h2bk->SetDirectory(0);
  TH1F* hbk = (TH1F*)h2bk->ProjectionY("hbk", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2bk->GetZaxis()->FindBin(lowy), h2bk->GetZaxis()->FindBin(highy));
  hbk->SetDirectory(0);

  double sig_sb =  hsig->Integral(hsig->GetXaxis()->FindBin(3.01),  hsig->GetXaxis()->FindBin(3.02));
  double bk_sb =   hbk->Integral(hbk->GetXaxis()->FindBin(3.01),  hbk->GetXaxis()->FindBin(3.02));
  double scale = sig_sb/bk_sb;
  hbk->Scale(scale);
  TH1F* hsig_bk = (TH1F*)hsig->Clone("hsig_bk");
  hsig_bk->Add(hbk,-1);
  setHistStyle(hsig_bk, kBlue, kFullCircle, 1.5);
  hsig_bk->Rebin();
  hsig_bk->Draw();

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
  hsig_bk->Draw("same");
  hsig_bk->Fit(fit,"R");
  resfit->SetParameter(0, fit->GetParameter(3));
  resfit->SetParameter(1, fit->GetParameter(4));
  resfit->Draw("same");
  setHistStyle(resfit, kRed-2, 9, 2.5 ,1);
  drawLine(lowx, 0, highx, 0, 1.5, 2, 1 );
  // cout<<"binwidth: "<< hsig_bk->GetBinWidth(1)<< endl;

  double sigma = fit->GetParameter(2);
  double mean = fit->GetParameter(1);
  double yield_me = fit->GetParameter(0)/hsig_bk->GetBinWidth(1);
  double yield_counts = hsig_bk->Integral(hsig_bk->GetXaxis()->FindBin(mean-2.5*sigma), hsig_bk->GetXaxis()->FindBin(mean+2.5*sigma));
  double bk_counts = hbk->Integral(hbk->GetXaxis()->FindBin(mean-2.5*sigma), hbk->GetXaxis()->FindBin(mean+2.5*sigma));
  double sp_counts = hsig->Integral(hbk->GetXaxis()->FindBin(mean-2.5*sigma), hbk->GetXaxis()->FindBin(mean+2.5*sigma));
  // double significance = yield_counts/sqrt(yield_counts+bk_counts);
  double significance = yield_me/sqrt(sp_counts);
  double s_me = yield_me/(fit->GetParError(0)/hsig_bk->GetBinWidth(1));

  TLegend* leg = new TLegend( 0.72, 0.68 ,0.9,0.9 );
  // leg->AddEntry(hbk, "mix-event(ME)","pl");
  // leg->AddEntry(hrt, "rotate d(RT) (scale)","pl");
  leg->AddEntry(hsig, "SE","pl");
  leg->AddEntry(hsig_bk, "SE-ME","pl");
  leg->Draw();
  drawLatex( 0.2,0.82,Form("ME/SE=%0.2f", 1./scale), 0.055);
  drawLatex( 0.2,0.75,Form("Yield=%0.2f", yield_me), 0.055);
  drawLatex( 0.2,0.68,Form("#sigma=%0.2f MeV", sigma*1000.), 0.055);
  drawLatex( 0.2,0.61,Form("nEvents=%0.0f M", nEvents_se/1e6), 0.055);
  drawLatex( 0.2,0.54,Form("S/#sqrt{S+B}=%0.0f", significance), 0.055);
  // drawLatex( 0.2,0.47,Form("S/#DeltaS=%0.0f (ME)", s_me), 0.055);
  // drawLatex( 0.2,0.4,Form("S/#DeltaS=%0.0f (RT)", s_rt), 0.055);
  drawLatex( 0.2,0.4,Form("Mean=%0.3f", mean), 0.055);
  drawLatex( 0.62,0.61,Form("%0.1f<y<%0.1f",lowy, highy ), 0.055);
  drawLatex( 0.62,0.54,Form("%0.1f<p_{T}<%0.1f GeV/c^{2}", lowpt, highpt), 0.055);
  // drawLatex( 0.62,0.47,Form("5-40%s", "%"), 0.055);
  drawLatex( 0.62,0.47,Form("0-80%s", "%"), 0.055);
  // drawBox( 2.97, hsig_bk->GetMinimum(),2.98, hsig_bk->GetMaximum()*0.5 , kBlue-9, 1001 , 1,0.3); 
  // drawBox( 2.989, hsig_bk->GetMinimum(),2.995, hsig_bk->GetMaximum()*0.5 , kOrange-4, 1001,1,0.3); 
  // drawBox( 3., hsig_bk->GetMinimum(),3.02, hsig_bk->GetMaximum()*0.5 , kBlue-9, 1001,1,0.3 ); 

  // addpdf(pdf);

  TH3F* h2chi2ndf_Sig = (TH3F*)f1->Get((fithistname+"Sig").Data())->Clone("h2chi2ndf_Sig");
  h2chi2ndf_Sig->SetDirectory(0);
  // TH3F* h2chi2ndf_SBL = (TH3F*)f1->Get((fithistname+"SBL").Data())->Clone("h2chi2ndf_SBL");
  // h2chi2ndf_SBL->SetDirectory(0);
  // TH3F* h2chi2ndf_SBR = (TH3F*)f1->Get((fithistname+"SBR").Data())->Clone("h2chi2ndf_SBR");
  // h2chi2ndf_SBR->SetDirectory(0);

  TH3F* h2chi2ndf_Bk = (TH3F*)f2->Get((fithistname+"Sig").Data())->Clone("h2chi2ndf_Bk");
  h2chi2ndf_Bk->SetDirectory(0);
  h2chi2ndf_Sig->Add( h2chi2ndf_Bk , -1*scale ); 
  h2chi2ndf_Sig->Draw();
  // h2chi2ndf_SBL->Add(h2chi2ndf_SBR);
  // h2chi2ndf_SBL->Scale(scale*h2chi2ndf_Bk->Integral()/h2chi2ndf_SBL->Integral());
  // h2chi2ndf_Sig->Add(h2chi2ndf_SBL, -1);

  TH1F* hchi2ndf_Sig = (TH1F*)h2chi2ndf_Sig->ProjectionY("hchi2ndf_Sig", h2chi2ndf_Sig->GetXaxis()->FindBin(lowpt), h2chi2ndf_Sig->GetXaxis()->FindBin(highpt),
                        h2chi2ndf_Sig->GetZaxis()->FindBin(lowy), h2chi2ndf_Sig->GetZaxis()->FindBin(highy));
  if (rebin>0) {
    hchi2ndf_Sig->Rebin(rebin);
    if (rebinedge>0) hchi2ndf_Sig = (TH1F*)reBinHist( (Xrange1+Xrange2)*rebinedge, Xrange1, Xrange2, hchi2ndf_Sig, 1 ); 
  }

  setHistStyle(hchi2ndf_Sig, kBlue, kFullCircle, 1.5);
  hchi2ndf_Sig->Draw();
  double realeff = (hchi2ndf_Sig->Integral(hchi2ndf_Sig->GetXaxis()->FindBin(cutslow), hchi2ndf_Sig->GetXaxis()->FindBin(cutshigh))/hchi2ndf_Sig->Integral(hchi2ndf_Sig->GetXaxis()->FindBin(cutslow), hchi2ndf_Sig->GetXaxis()->FindBin(10)));
  double fakeeff = (hchi2ndf_H3L->Integral(hchi2ndf_H3L->GetXaxis()->FindBin(cutslow), hchi2ndf_H3L->GetXaxis()->FindBin(cutshigh))/hchi2ndf_H3L->Integral(hchi2ndf_H3L->GetXaxis()->FindBin(cutslow), hchi2ndf_H3L->GetXaxis()->FindBin(10)));
  error = 0;
  return fakeeff/realeff;

  TH1F* hRatio = (TH1F*)hchi2ndf_Sig->Clone("hratio");
  TF1* fMyFit = new TF1("fMyFit", MyFitFun, 0,100,2);
  TF1* fL = new TF1("fL", MyFitFun, 0,100,2);
  TF1* fH = new TF1("fH", MyFitFun, 0,100,2);
  fMyFit->SetLineColor(kRed);
  fL->SetLineColor(kRed);
  fH->SetLineColor(kRed);
  fMyFit->SetParLimits(0, 0, 1e10);
  fMyFit->SetParLimits(1, 0, 1e10);
  if (Xrange1>=0 && Xrange2>=0) hchi2ndf_Sig->GetXaxis()->SetRangeUser( Xrange1, Xrange2);
  hchi2ndf_Sig->Fit(fMyFit,"RBL");
  double par[2];
  fMyFit->GetParameters( par );
  fMyFit->SetParameters( par );
  if (lowy<-0.7 && highy<-0.4) fMyFit->SetParLimits(0, fMyFit->GetParameter(0)*0.2, fMyFit->GetParameter(0)*0.9);
  hchi2ndf_Sig->Fit(fMyFit,"RBL");
  fMyFit->GetParameters( par );

  TH1F* hLd = (TH1F*)hchi2ndf_Ld->Clone("hLd");
  hLd->Scale(fMyFit->GetParameter(0)*fMyFit->GetParameter(1));
  hLd->Draw("same");
  TH1F* hH3L = (TH1F*)hchi2ndf_H3L->Clone("hH3L");
  TH1F* hH3L_L = (TH1F*)hchi2ndf_H3L->Clone("hH3L_L");
  TH1F* hH3L_H = (TH1F*)hchi2ndf_H3L->Clone("hH3L_H");
  hH3L->Scale(fMyFit->GetParameter(0));

  setHistStyle(hLd, kBlack, kOpenStar, 2);
  setHistStyle(hH3L, kRed, kOpenCircle, 1.5);
  hH3L->Draw("same");
  gPad->SetLogy();
  double purity = hH3L->Integral(hH3L->GetXaxis()->FindBin(cutslow+1e-6),hH3L->GetXaxis()->FindBin(cutshigh-1e-6))/fMyFit->Integral(cutslow,cutshigh)*hH3L->GetBinWidth(1);
  hH3L_L->Scale(fMyFit->GetParameter(0));
  fL->SetParameter(1, fMyFit->GetParameter(1)-fMyFit->GetParError(1));
  fL->SetParameter(0, fMyFit->GetParameter(0));
  fH->SetParameter(1, fMyFit->GetParameter(1)+fMyFit->GetParError(1));
  fH->SetParameter(0, fMyFit->GetParameter(0));
  fH->Draw("same");
  fL->Draw("same");
  double purity_low = hH3L_L->Integral(hH3L_L->GetXaxis()->FindBin(cutslow+1e-6),hH3L_L->GetXaxis()->FindBin(cutshigh-1e-6))/fL->Integral(cutslow,cutshigh)*hH3L_L->GetBinWidth(1);
  hH3L_H->Scale(fMyFit->GetParameter(0));
  double purity_high = hH3L_H->Integral(hH3L_H->GetXaxis()->FindBin(cutslow+1e-6),hH3L_H->GetXaxis()->FindBin(cutshigh-1e-6))/fH->Integral(cutslow,cutshigh)*hH3L_H->GetBinWidth(1);
  error = abs(purity_high-purity)>abs(purity_low-purity)? abs(purity_high-purity) : abs(purity_low-purity);
  // drawLatex(0.5,0.6,Form("purity = %0.2f  + %0.4f -%0.4f",purity, purity_low-purity, purity_high-purity),0.055);
  drawLatex(0.5,0.6,Form("purity = %0.2f#pm %0.4f",purity, purity_low-purity),0.055);
  drawLatex(0.65,0.5,Form("%0.1f<pT<%0.1f", lowpt, highpt),0.055);
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
      hRatio->SetBinError(i, fMyFit->Eval(hRatio->GetBinCenter(i))/pow(hRatio->GetBinContent(i),2 )*hRatio->GetBinError(i));
      hRatio->SetBinContent(i, fMyFit->Eval(hRatio->GetBinCenter(i))/hRatio->GetBinContent(i) );
    }
  }
  hRatio->GetYaxis()->SetTitle("Data/Fit");
  hRatio->GetYaxis()->SetTitleOffset(0.8);
  hRatio->GetYaxis()->SetRangeUser(0.5,1.5);
  hRatio->GetYaxis()->SetTitleSize(0.075);
  hRatio->GetYaxis()->SetLabelSize(0.075);
  hRatio->GetXaxis()->SetTitleSize(0.075);
  hRatio->GetXaxis()->SetLabelSize(0.075);
  hRatio->GetXaxis()->SetTitleOffset(0.9);
  hRatio->GetXaxis()->SetRangeUser(Xrange1,Xrange2);
  hRatio->Draw();
  drawLine(hRatio->GetXaxis()->GetXmin(), 1, hRatio->GetXaxis()->GetXmax(),1,1.5,2,1);
  c->cd();
  addpdf(pdf);
  // delete hchi2ndf_Ld;
  // delete hchi2ndf_H3L;
  return purity;
}
double calpurityCut(TString fithistname, TString scalehistname, double cutslow, double cutshigh, double highpt, double lowpt, double highy, double lowy, TFile* f1, TFile* f2, double scale_i, TFile* fMc, TFile* fMc_ld, TCanvas* c, TPDF* pdf, TString xTitle, int rebin, double Xrange1, double Xrange2,double &error, double  rebinedge)
{
  //fiting test
  f1->cd();
  TH1F* hcent_se = (TH1F*)f1->Get("hcent")->Clone("hcent_se");
  double nEvents_se = hcent_se->Integral(4, 9);

  TH3F* h2sig = (TH3F*)f1->Get(scalehistname.Data())->Clone("hptH3Lmass_sig");
  h2sig->SetDirectory(0);
  TH1F* hsig = (TH1F*)h2sig->ProjectionY("hsig", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2sig->GetZaxis()->FindBin(lowy), h2sig->GetZaxis()->FindBin(highy));
  hsig->SetDirectory(0);

  TH3F* h2bk = (TH3F*)f2->Get(scalehistname.Data())->Clone("hptH3Lmass_bk");
  h2bk->SetDirectory(0);
  TH1F* hbk = (TH1F*)h2bk->ProjectionY("hbk", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2bk->GetZaxis()->FindBin(lowy), h2bk->GetZaxis()->FindBin(highy));
  hbk->SetDirectory(0);

  double sig_sb =  hsig->Integral(hsig->GetXaxis()->FindBin(3.01),  hsig->GetXaxis()->FindBin(3.02));
  double bk_sb =   hbk->Integral(hbk->GetXaxis()->FindBin(3.01),  hbk->GetXaxis()->FindBin(3.02));
  double scale = sig_sb/bk_sb;
  hbk->Scale(scale);
  TH1F* hsig_bk = (TH1F*)hsig->Clone("hsig_bk");
  hsig_bk->Add(hbk,-1);
  setHistStyle(hsig_bk, kBlue, kFullCircle, 1.5);
  hsig_bk->Rebin();
  hsig_bk->Draw();

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
  hsig_bk->Draw("same");
  hsig_bk->Fit(fit,"R");
  resfit->SetParameter(0, fit->GetParameter(3));
  resfit->SetParameter(1, fit->GetParameter(4));
  resfit->Draw("same");
  setHistStyle(resfit, kRed-2, 9, 2.5 ,1);
  drawLine(lowx, 0, highx, 0, 1.5, 2, 1 );
  // cout<<"binwidth: "<< hsig_bk->GetBinWidth(1)<< endl;

  double sigma = fit->GetParameter(2);
  double mean = fit->GetParameter(1);
  double yield_me = fit->GetParameter(0)/hsig_bk->GetBinWidth(1);
  double yield_counts = hsig_bk->Integral(hsig_bk->GetXaxis()->FindBin(mean-2.5*sigma), hsig_bk->GetXaxis()->FindBin(mean+2.5*sigma));
  double bk_counts = hbk->Integral(hbk->GetXaxis()->FindBin(mean-2.5*sigma), hbk->GetXaxis()->FindBin(mean+2.5*sigma));
  double sp_counts = hsig->Integral(hbk->GetXaxis()->FindBin(mean-2.5*sigma), hbk->GetXaxis()->FindBin(mean+2.5*sigma));
  // double significance = yield_counts/sqrt(yield_counts+bk_counts);
  double significance = yield_me/sqrt(sp_counts);
  double s_me = yield_me/(fit->GetParError(0)/hsig_bk->GetBinWidth(1));

  TLegend* leg = new TLegend( 0.72, 0.68 ,0.9,0.9 );
  // leg->AddEntry(hbk, "mix-event(ME)","pl");
  // leg->AddEntry(hrt, "rotate d(RT) (scale)","pl");
  leg->AddEntry(hsig, "SE","pl");
  leg->AddEntry(hsig_bk, "SE-ME","pl");
  leg->Draw();
  drawLatex( 0.2,0.82,Form("ME/SE=%0.2f", 1./scale), 0.055);
  drawLatex( 0.2,0.75,Form("Yield=%0.2f", yield_me), 0.055);
  drawLatex( 0.2,0.68,Form("#sigma=%0.2f MeV", sigma*1000.), 0.055);
  drawLatex( 0.2,0.61,Form("nEvents=%0.0f M", nEvents_se/1e6), 0.055);
  drawLatex( 0.2,0.54,Form("S/#sqrt{S+B}=%0.0f", significance), 0.055);
  // drawLatex( 0.2,0.47,Form("S/#DeltaS=%0.0f (ME)", s_me), 0.055);
  // drawLatex( 0.2,0.4,Form("S/#DeltaS=%0.0f (RT)", s_rt), 0.055);
  drawLatex( 0.2,0.4,Form("Mean=%0.3f", mean), 0.055);
  drawLatex( 0.62,0.61,Form("%0.1f<y<%0.1f",lowy, highy ), 0.055);
  drawLatex( 0.62,0.54,Form("%0.1f<p_{T}<%0.1f GeV/c^{2}", lowpt, highpt), 0.055);
  // drawLatex( 0.62,0.47,Form("5-40%s", "%"), 0.055);
  drawLatex( 0.62,0.47,Form("0-80%s", "%"), 0.055);
  // drawBox( 2.97, hsig_bk->GetMinimum(),2.98, hsig_bk->GetMaximum()*0.5 , kBlue-9, 1001 , 1,0.3); 
  // drawBox( 2.989, hsig_bk->GetMinimum(),2.995, hsig_bk->GetMaximum()*0.5 , kOrange-4, 1001,1,0.3); 
  // drawBox( 3., hsig_bk->GetMinimum(),3.02, hsig_bk->GetMaximum()*0.5 , kBlue-9, 1001,1,0.3 ); 
  // addpdf(pdf);

  // TH3F* h2chi2ndf_Ld = (TH3F*)fMc_ld->Get(fithistname.Data())->Clone("h2chi2ndf_Ld");
  TH3F* h2chi2ndf_Ld = (TH3F*)fMc_ld->Get((fithistname+"Sig").Data())->Clone("h2chi2ndf_Ld");
  h2chi2ndf_Ld->SetDirectory(0);
  hchi2ndf_Ld = (TH1F*)h2chi2ndf_Ld->ProjectionY("hchi2ndf_Ld", h2chi2ndf_Ld->GetXaxis()->FindBin(lowpt), h2chi2ndf_Ld->GetXaxis()->FindBin(highpt),
                h2chi2ndf_Ld->GetZaxis()->FindBin(lowy), h2chi2ndf_Ld->GetZaxis()->FindBin(highy)  );
  if (rebin>0) {
    hchi2ndf_Ld->Rebin(rebin);
    if (rebinedge>0) hchi2ndf_Ld = (TH1F*)reBinHist( (Xrange1+Xrange2)*rebinedge, Xrange1, Xrange2, hchi2ndf_Ld, 1 ); 
  }
  
  // TH3F* h2chi2ndf_H3L = (TH3F*)fMc->Get(fithistname.Data())->Clone("h2chi2ndf_H3L");
  TH3F* h2chi2ndf_H3L = (TH3F*)fMc->Get((fithistname+"Sig").Data())->Clone("h2chi2ndf_H3L");
  h2chi2ndf_H3L->SetDirectory(0);
  hchi2ndf_H3L = (TH1F*)h2chi2ndf_H3L->ProjectionY("hchi2ndf_H3L", h2chi2ndf_H3L->GetXaxis()->FindBin(lowpt), h2chi2ndf_H3L->GetXaxis()->FindBin(highpt),
                  h2chi2ndf_H3L->GetZaxis()->FindBin(lowy), h2chi2ndf_H3L->GetZaxis()->FindBin(highy));
  if (rebin>0) {
    hchi2ndf_H3L->Rebin(rebin);
    if (rebinedge>0) hchi2ndf_H3L = (TH1F*)reBinHist( (Xrange1+Xrange2)*rebinedge, Xrange1, Xrange2, hchi2ndf_H3L, 1 ); 
  }
 
  Norm(hchi2ndf_H3L);
  Norm(hchi2ndf_Ld);
  hchi2ndf_Ld->Draw();
  hchi2ndf_Ld->GetXaxis()->SetTitle(xTitle.Data());
  hchi2ndf_H3L->Draw("same");
  hchi2ndf_H3L->SetMarkerColor(kRed);
  hchi2ndf_H3L->SetMarkerStyle(kOpenCircle);

  f1->cd();
  TH3F* h2chi2ndf_Sig = (TH3F*)f1->Get((fithistname+"Sig").Data())->Clone("h2chi2ndf_Sig");
  h2chi2ndf_Sig->SetDirectory(0);
  // TH3F* h2chi2ndf_SBL = (TH3F*)f1->Get((fithistname+"SBL").Data())->Clone("h2chi2ndf_SBL");
  // h2chi2ndf_SBL->SetDirectory(0);
  // TH3F* h2chi2ndf_SBR = (TH3F*)f1->Get((fithistname+"SBR").Data())->Clone("h2chi2ndf_SBR");
  // h2chi2ndf_SBR->SetDirectory(0);

  TH3F* h2chi2ndf_Bk = (TH3F*)f2->Get((fithistname+"Sig").Data())->Clone("h2chi2ndf_Bk");
  h2chi2ndf_Bk->SetDirectory(0);
  h2chi2ndf_Sig->Add( h2chi2ndf_Bk , -1*scale ); 

  // h2chi2ndf_Sig->Draw();
  // h2chi2ndf_SBL->Add(h2chi2ndf_SBR);
  // h2chi2ndf_SBL->Scale(scale*h2chi2ndf_Bk->Integral()/h2chi2ndf_SBL->Integral());
  // h2chi2ndf_Sig->Add(h2chi2ndf_SBL, -1);

  TH1F* hchi2ndf_Sig = (TH1F*)h2chi2ndf_Sig->ProjectionY("hchi2ndf_Sig", h2chi2ndf_Sig->GetXaxis()->FindBin(lowpt), h2chi2ndf_Sig->GetXaxis()->FindBin(highpt),
                        h2chi2ndf_Sig->GetZaxis()->FindBin(lowy), h2chi2ndf_Sig->GetZaxis()->FindBin(highy));
  if (rebin>0) {
    hchi2ndf_Sig->Rebin(rebin);
    if (rebinedge>0) hchi2ndf_Sig = (TH1F*)reBinHist( (Xrange1+Xrange2)*rebinedge, Xrange1, Xrange2, hchi2ndf_Sig, 1 ); 
  }

  setHistStyle(hchi2ndf_Sig, kBlue, kFullCircle, 1.5);

  TH1F* hSig = (TH1F*)hchi2ndf_Sig->Clone("hSig");
  Norm(hSig);
  hSig->SetMarkerColor(kBlue);
  hSig->Draw("same");
  addpdf(pdf);
  // double realeff = (hchi2ndf_Sig->Integral(hchi2ndf_Sig->GetXaxis()->FindBin(cutslow), hchi2ndf_Sig->GetXaxis()->FindBin(cutshigh))/hchi2ndf_Sig->Integral(hchi2ndf_Sig->GetXaxis()->FindBin(cutslow), hchi2ndf_Sig->GetXaxis()->FindBin(10)));
  // double fakeeff = (hchi2ndf_H3L->Integral(hchi2ndf_H3L->GetXaxis()->FindBin(cutslow), hchi2ndf_H3L->GetXaxis()->FindBin(cutshigh))/hchi2ndf_H3L->Integral(hchi2ndf_H3L->GetXaxis()->FindBin(cutslow), hchi2ndf_H3L->GetXaxis()->FindBin(10)));
  // error = 0;
  // return fakeeff/realeff;

  c->cd();
  TPad*    p11 = new TPad("upperPad", "upperPad",.005, .4, .995, .995);
  TPad*    p12 = new TPad("lowerPad", "lowerPad", .005, .005, .995, .4);
  p11->Draw(); 			       
  p12->Draw();  
  p11->SetBottomMargin(0);
  p12->SetTopMargin(0);
  p12->SetBottomMargin(0.18);

  p11->cd();

  hchi2ndf_Sig->Draw();

  TH1F* hRatio = (TH1F*)hchi2ndf_Sig->Clone("hratio");
  TF1* fMyFit = new TF1("fMyFit", MyFitFun, 0,100,2);
  TF1* fL = new TF1("fL", MyFitFun, 0,100,2);
  TF1* fH = new TF1("fH", MyFitFun, 0,100,2);
  fMyFit->SetLineColor(kRed);
  fL->SetLineColor(kRed);
  fH->SetLineColor(kRed);
  fMyFit->SetParLimits(0, 0, 1e10);
  fMyFit->SetParLimits(1, 0, 4);
  if (Xrange1>=0 && Xrange2>=0) hchi2ndf_Sig->GetXaxis()->SetRangeUser( Xrange1, Xrange2);
  hchi2ndf_Sig->Fit(fMyFit,"RBL");
  double par[2];
  fMyFit->GetParameters( par );
  // par[1]=0.5*par[1];
  // par[0]=1.5*par[0];
  if (par[1]>3)
    fMyFit->SetParLimits(1, 0.3, 3);
  if (par[1]<0.3)
    fMyFit->SetParLimits(1, 0.3, 3);

  // fMyFit->SetParameters( par );
  // fMyFit->SetParLimits(1, fMyFit->GetParameter(0)*0.2, fMyFit->GetParameter(0)*1);
  hchi2ndf_Sig->Fit(fMyFit,"RBL");

  TH1F* hLd = (TH1F*)hchi2ndf_Ld->Clone("hLd");
  hLd->Scale(fMyFit->GetParameter(0)*fMyFit->GetParameter(1));
  hLd->Draw("same");
  TH1F* hH3L = (TH1F*)hchi2ndf_H3L->Clone("hH3L");
  TH1F* hH3L_L = (TH1F*)hchi2ndf_H3L->Clone("hH3L_L");
  TH1F* hH3L_H = (TH1F*)hchi2ndf_H3L->Clone("hH3L_H");
  hH3L->Scale(fMyFit->GetParameter(0));

  setHistStyle(hLd, kBlack, kOpenStar, 2);
  setHistStyle(hH3L, kRed, kOpenCircle, 1.5);
  hH3L->Draw("same");
  gPad->SetLogy();
  double purity = hH3L->Integral(hH3L->GetXaxis()->FindBin(cutslow+1e-6),hH3L->GetXaxis()->FindBin(cutshigh-1e-6))/fMyFit->Integral(cutslow,cutshigh)*hH3L->GetBinWidth(1);
  hH3L_L->Scale(fMyFit->GetParameter(0));
  fL->SetParameter(1, fMyFit->GetParameter(1)-fMyFit->GetParError(1));
  fL->SetParameter(0, fMyFit->GetParameter(0));
  fH->SetParameter(1, fMyFit->GetParameter(1)+fMyFit->GetParError(1));
  fH->SetParameter(0, fMyFit->GetParameter(0));
  fH->Draw("same");
  fL->Draw("same");
  double purity_low = hH3L_L->Integral(hH3L_L->GetXaxis()->FindBin(cutslow+1e-6),hH3L_L->GetXaxis()->FindBin(cutshigh-1e-6))/fL->Integral(cutslow,cutshigh)*hH3L_L->GetBinWidth(1);
  hH3L_H->Scale(fMyFit->GetParameter(0));
  double purity_high = hH3L_H->Integral(hH3L_H->GetXaxis()->FindBin(cutslow+1e-6),hH3L_H->GetXaxis()->FindBin(cutshigh-1e-6))/fH->Integral(cutslow,cutshigh)*hH3L_H->GetBinWidth(1);
  error = abs(purity_high-purity)>abs(purity_low-purity)? abs(purity_high-purity) : abs(purity_low-purity);
  // drawLatex(0.5,0.6,Form("purity = %0.2f  + %0.4f -%0.4f",purity, purity_low-purity, purity_high-purity),0.055);
  drawLatex(0.5,0.6,Form("purity = %0.2f#pm %0.4f",purity, purity_low-purity),0.055);
  drawLatex(0.65,0.5,Form("%0.1f<pT<%0.1f", lowpt, highpt),0.055);
  drawLatex(0.85,0.5,Form("f=%0.2f",fMyFit->GetParameter(1)),0.055);
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
      hRatio->SetBinError(i, fMyFit->Eval(hRatio->GetBinCenter(i))/pow(hRatio->GetBinContent(i),2 )*hRatio->GetBinError(i));
      hRatio->SetBinContent(i, fMyFit->Eval(hRatio->GetBinCenter(i))/hRatio->GetBinContent(i) );
    }
  }
  hRatio->GetYaxis()->SetTitle("Data/Fit");
  hRatio->GetYaxis()->SetTitleOffset(0.8);
  hRatio->GetYaxis()->SetRangeUser(0.5,1.5);
  hRatio->GetYaxis()->SetTitleSize(0.075);
  hRatio->GetYaxis()->SetLabelSize(0.075);
  hRatio->GetXaxis()->SetTitleSize(0.075);
  hRatio->GetXaxis()->SetLabelSize(0.075);
  hRatio->GetXaxis()->SetTitleOffset(0.9);
  hRatio->GetXaxis()->SetRangeUser(Xrange1,Xrange2);
  hRatio->Draw();
  drawLine(hRatio->GetXaxis()->GetXmin(), 1, hRatio->GetXaxis()->GetXmax(),1,1.5,2,1);
  c->cd();
  addpdf(pdf);
  // delete hchi2ndf_Ld;
  // delete hchi2ndf_H3L;
  return purity;
}
double fityield( double lowpt, double highpt, double lowy, double highy, double& err, TFile* f1, TFile* f2, TCanvas* c, TPDF* pdf, int centL=4, int centH=9)
{ 
  c->cd();
  // TPad*    p11 = new TPad("upperPad", "upperPad",.005, .4, .995, .995);
  // TPad*    p12 = new TPad("lowerPad", "lowerPad", .005, .005, .995, .4);
  // p11->Draw(); 			       
  // p12->Draw();  
  // p11->SetBottomMargin(0);
  // p12->SetTopMargin(0);
  // p12->SetBottomMargin(0.18);
  //
  // p11->cd();
  // c->Divide(2,1);
  // c->cd(1);
  TH1F* hcent_se = (TH1F*)f1->Get("hcent")->Clone("hcent_se");
  // double nEvents_se = hcent_se->Integral(4,9);
  double nEvents_se = hcent_se->Integral(centL,centH);

  TH3F* h2sig = (TH3F*)f1->Get("hH3LMassPtY")->Clone("hptH3Lmass_sig");
  h2sig->SetDirectory(0);
  TH1F* hsig = (TH1F*)h2sig->ProjectionY("hsig", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2sig->GetZaxis()->FindBin(lowy+1e-6), h2sig->GetZaxis()->FindBin(highy-1e-6));
  hsig->SetDirectory(0);

  TH3F* h2bk = (TH3F*)f2->Get("hH3LMassPtY")->Clone("hptH3Lmass_bk");
  h2bk->SetDirectory(0);
  TH1F* hbk = (TH1F*)h2bk->ProjectionY("hbk", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2bk->GetZaxis()->FindBin(lowy+1e-6), h2bk->GetZaxis()->FindBin(highy-1e-6));
  hbk->SetDirectory(0);

  double sig_sb =  hsig->Integral(hsig->GetXaxis()->FindBin(3.00),  hsig->GetXaxis()->FindBin(3.02));
  double bk_sb =   hbk->Integral(hbk->GetXaxis()->FindBin(3.00),  hbk->GetXaxis()->FindBin(3.02));
  double scale = sig_sb/bk_sb;
  hbk->Scale(scale);
  // hsig->Draw();
  // hbk->Draw("same");

  TH1F* hsig_bk = (TH1F*)hsig->Clone("hsig_bk");
  hsig_bk->Add(hbk,-1);
  setHistStyle(hsig_bk, kBlue, kFullCircle, 1.5);
  hsig_bk->Rebin(4);
  hsig_bk->Draw();

  TF1* fit = new TF1("fit" ,"gausn(0)+pol1(3)", 2.97,3.02 );
  TF1* resfit = new TF1("resfit" ,"pol1", 2.95,3.05 );
  hsig_bk->GetXaxis()->SetRangeUser(2.97,2.985);
  hsig_bk->Fit(resfit,"R");
  fit->SetLineColor(kRed);
  double yield_bc = hsig_bk->Integral(hsig_bk->GetXaxis()->FindBin(2.98), hsig_bk->GetXaxis()->FindBin(3.));
  // double para[5]={yield_bc*hsig_bk->GetBinWidth(1)/sqrt(2*3.1415), 2.991, 0.0015,  resfit->GetParameter(0), resfit->GetParameter(1)};
  double para[5]={yield_bc*hsig_bk->GetBinWidth(1), 2.992, 0.0014,  resfit->GetParameter(0), resfit->GetParameter(1)};
  fit->SetParameters(para);
  double lowx=2.97 ,highx =3.02;
  hsig_bk->GetXaxis()->SetRangeUser(lowx,highx);
  hsig_bk->Draw("same");
  hsig_bk->Fit(fit,"R");
  resfit->SetParameter(0, fit->GetParameter(3));
  resfit->SetParameter(1, fit->GetParameter(4));
  resfit->Draw("same");
  setHistStyle(resfit, kRed-2, 9, 2.5 ,1);
  drawLine(lowx, 0, highx, 0, 1.5, 2, 1 );
  // cout<<"binwidth: "<< hsig_bk->GetBinWidth(1)<< endl;

  double sigma = fit->GetParameter(2);
  double mean = fit->GetParameter(1);
  double yield_me = fit->GetParameter(0)/hsig_bk->GetBinWidth(1);
  double yield_counts = hsig_bk->Integral(hsig_bk->GetXaxis()->FindBin(mean-2.5*sigma), hsig_bk->GetXaxis()->FindBin(mean+2.5*sigma));
  double bk_counts = hbk->Integral(hbk->GetXaxis()->FindBin(mean-2.5*sigma), hbk->GetXaxis()->FindBin(mean+2.5*sigma));
  double sp_counts = hsig->Integral(hbk->GetXaxis()->FindBin(mean-2.5*sigma), hbk->GetXaxis()->FindBin(mean+2.5*sigma));
  // double significance = yield_counts/sqrt(yield_counts+bk_counts);
  double significance = yield_me/sqrt(sp_counts);
  double s_me = yield_me/(fit->GetParError(0)/hsig_bk->GetBinWidth(1));
  err = (fit->GetParError(0)/hsig_bk->GetBinWidth(1));

  TLegend* leg = new TLegend( 0.72, 0.68 ,0.9,0.9 );
  // leg->AddEntry(hbk, "mix-event(ME)","pl");
  // leg->AddEntry(hrt, "rotate d(RT) (scale)","pl");
  leg->AddEntry(hsig, "SE","pl");
  leg->AddEntry(hsig_bk, "SE-ME","pl");
  leg->Draw();
  drawLatex( 0.2,0.82,Form("ME/SE=%0.2f", 1./scale), 0.055);
  drawLatex( 0.2,0.75,Form("Yield=%0.2f", yield_me), 0.055);
  drawLatex( 0.2,0.68,Form("#sigma=%0.2f MeV", sigma*1000.), 0.055);
  drawLatex( 0.2,0.61,Form("nEvents=%0.0f M", nEvents_se/1e6), 0.055);
  drawLatex( 0.2,0.54,Form("S/#sqrt{S+B}=%0.0f", significance), 0.055);
  // drawLatex( 0.2,0.47,Form("S/#DeltaS=%0.0f (ME)", s_me), 0.055);
  // drawLatex( 0.2,0.4,Form("S/#DeltaS=%0.0f (RT)", s_rt), 0.055);
  drawLatex( 0.2,0.4,Form("Mean=%0.3f", mean), 0.055);
  drawLatex( 0.62,0.61,Form("%0.2f<y<%0.2f",lowy, highy ), 0.055);
  drawLatex( 0.62,0.54,Form("%0.1f<p_{T}<%0.1f GeV/c^{2}", lowpt, highpt), 0.055);
  // drawLatex( 0.62,0.47,Form("5-40%s", "%"), 0.055);
  drawLatex( 0.62,0.47,Form("0-80%s", "%"), 0.055);
  // drawBox( 2.97, hsig_bk->GetMinimum(),2.98, hsig_bk->GetMaximum()*0.5 , kBlue-9, 1001 , 1,0.3); 
  // drawBox( 2.989, hsig_bk->GetMinimum(),2.995, hsig_bk->GetMaximum()*0.5 , kOrange-4, 1001,1,0.3); 
  // drawBox( 3., hsig_bk->GetMinimum(),3.02, hsig_bk->GetMaximum()*0.5 , kBlue-9, 1001,1,0.3 ); 

  addpdf(pdf);

  return (yield_counts>20 && sigma<0.005 && significance>3.5)?yield_me:0;
}
double fityield2(  TString histname, double lowpt, double highpt, double lowy, double highy, double& err, TFile* f1, TFile* f2, TCanvas* c, TPDF* pdf, int centL=4, int centH=9)
{ 
  c->cd();
  // TPad*    p11 = new TPad("upperPad", "upperPad",.005, .4, .995, .995);
  // TPad*    p12 = new TPad("lowerPad", "lowerPad", .005, .005, .995, .4);
  // p11->Draw(); 			       
  // p12->Draw();  
  // p11->SetBottomMargin(0);
  // p12->SetTopMargin(0);
  // p12->SetBottomMargin(0.18);
  //
  // p11->cd();
  // c->Divide(2,1);
  // c->cd(1);
  TH1F* hcent_se = (TH1F*)f1->Get("hcent")->Clone("hcent_se");
  // double nEvents_se = hcent_se->Integral(4,9);
  double nEvents_se = hcent_se->Integral(centL,centH);

  TH3F* h2sig = (TH3F*)f1->Get(Form("%s", histname.Data()))->Clone("hptH3Lmass_sig");
  h2sig->SetDirectory(0);
  TH1F* hsig = (TH1F*)h2sig->ProjectionY("hsig", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2sig->GetZaxis()->FindBin(lowy+1e-6), h2sig->GetZaxis()->FindBin(highy-1e-6));
  hsig->SetDirectory(0);

  TH3F* h2bk = (TH3F*)f2->Get(Form("%s", histname.Data()))->Clone("hptH3Lmass_bk");
  h2bk->SetDirectory(0);
  TH1F* hbk = (TH1F*)h2bk->ProjectionY("hbk", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2bk->GetZaxis()->FindBin(lowy+1e-6), h2bk->GetZaxis()->FindBin(highy-1e-6));
  hbk->SetDirectory(0);

  double sig_sb =  hsig->Integral(hsig->GetXaxis()->FindBin(3.00),  hsig->GetXaxis()->FindBin(3.02));
  double bk_sb =   hbk->Integral(hbk->GetXaxis()->FindBin(3.00),  hbk->GetXaxis()->FindBin(3.02));
  double scale = sig_sb/bk_sb;
  hbk->Scale(scale);
  // hsig->Draw();
  // hbk->Draw("same");

  TH1F* hsig_bk = (TH1F*)hsig->Clone("hsig_bk");
  hsig_bk->Add(hbk,-1);
  setHistStyle(hsig_bk, kBlue, kFullCircle, 1.5);
  hsig_bk->Rebin();
  hsig_bk->Draw();

  TF1* fit = new TF1("fit" ,"gausn(0)+pol1(3)", 2.97,3.02 );
  TF1* resfit = new TF1("resfit" ,"pol1", 2.95,3.05 );
  hsig_bk->GetXaxis()->SetRangeUser(2.97,2.985);
  hsig_bk->Fit(resfit,"R");
  fit->SetLineColor(kRed);
  double yield_bc = hsig_bk->Integral(hsig_bk->GetXaxis()->FindBin(2.98), hsig_bk->GetXaxis()->FindBin(3.));
  // double para[5]={yield_bc*hsig_bk->GetBinWidth(1)/sqrt(2*3.1415), 2.991, 0.0015,  resfit->GetParameter(0), resfit->GetParameter(1)};
  double para[5]={yield_bc*hsig_bk->GetBinWidth(1), 2.992, 0.0014,  resfit->GetParameter(0), resfit->GetParameter(1)};
  fit->SetParameters(para);
  double lowx=2.97 ,highx =3.02;
  hsig_bk->GetXaxis()->SetRangeUser(lowx,highx);
  hsig_bk->Draw("same");
  hsig_bk->Fit(fit,"R");
  resfit->SetParameter(0, fit->GetParameter(3));
  resfit->SetParameter(1, fit->GetParameter(4));
  resfit->Draw("same");
  setHistStyle(resfit, kRed-2, 9, 2.5 ,1);
  drawLine(lowx, 0, highx, 0, 1.5, 2, 1 );
  // cout<<"binwidth: "<< hsig_bk->GetBinWidth(1)<< endl;

  double sigma = fit->GetParameter(2);
  double mean = fit->GetParameter(1);
  double yield_me = fit->GetParameter(0)/hsig_bk->GetBinWidth(1);
  double yield_counts = hsig_bk->Integral(hsig_bk->GetXaxis()->FindBin(mean-2.5*sigma), hsig_bk->GetXaxis()->FindBin(mean+2.5*sigma));
  double bk_counts = hbk->Integral(hbk->GetXaxis()->FindBin(mean-2.5*sigma), hbk->GetXaxis()->FindBin(mean+2.5*sigma));
  double sp_counts = hsig->Integral(hbk->GetXaxis()->FindBin(mean-2.5*sigma), hbk->GetXaxis()->FindBin(mean+2.5*sigma));
  // double significance = yield_counts/sqrt(yield_counts+bk_counts);
  double significance = yield_me/sqrt(sp_counts);
  double s_me = yield_me/(fit->GetParError(0)/hsig_bk->GetBinWidth(1));
  err = (fit->GetParError(0)/hsig_bk->GetBinWidth(1));

  TLegend* leg = new TLegend( 0.72, 0.68 ,0.9,0.9 );
  // leg->AddEntry(hbk, "mix-event(ME)","pl");
  // leg->AddEntry(hrt, "rotate d(RT) (scale)","pl");
  leg->AddEntry(hsig, "SE","pl");
  leg->AddEntry(hsig_bk, "SE-ME","pl");
  leg->Draw();
  drawLatex( 0.2,0.82,Form("ME/SE=%0.2f", 1./scale), 0.055);
  drawLatex( 0.2,0.75,Form("Yield=%0.2f", yield_me), 0.055);
  drawLatex( 0.2,0.68,Form("#sigma=%0.2f MeV", sigma*1000.), 0.055);
  drawLatex( 0.2,0.61,Form("nEvents=%0.0f M", nEvents_se/1e6), 0.055);
  drawLatex( 0.2,0.54,Form("S/#sqrt{S+B}=%0.0f", significance), 0.055);
  // drawLatex( 0.2,0.47,Form("S/#DeltaS=%0.0f (ME)", s_me), 0.055);
  // drawLatex( 0.2,0.4,Form("S/#DeltaS=%0.0f (RT)", s_rt), 0.055);
  drawLatex( 0.2,0.4,Form("Mean=%0.3f", mean), 0.055);
  drawLatex( 0.62,0.61,Form("%0.2f<y<%0.2f",lowy, highy ), 0.055);
  drawLatex( 0.62,0.54,Form("%0.1f<p_{T}<%0.1f GeV/c^{2}", lowpt, highpt), 0.055);
  // drawLatex( 0.62,0.47,Form("5-40%s", "%"), 0.055);
  drawLatex( 0.62,0.47,Form("0-80%s", "%"), 0.055);
  // drawBox( 2.97, hsig_bk->GetMinimum(),2.98, hsig_bk->GetMaximum()*0.5 , kBlue-9, 1001 , 1,0.3); 
  // drawBox( 2.989, hsig_bk->GetMinimum(),2.995, hsig_bk->GetMaximum()*0.5 , kOrange-4, 1001,1,0.3); 
  // drawBox( 3., hsig_bk->GetMinimum(),3.02, hsig_bk->GetMaximum()*0.5 , kBlue-9, 1001,1,0.3 ); 

  addpdf(pdf);

  return (yield_counts>20 && sigma<0.005 && significance>3.5)?yield_me:0;
}
double fityield2(  TString histname, int icut, double lowpt, double highpt, double lowy, double highy, double& err, TFile* f1, TFile* f2, TCanvas* c, TPDF* pdf, int centL=4, int centH=9)
{ 
  c->cd();

  TH1F* hcent_se = (TH1F*)f1->Get("hcent")->Clone("hcent_se");
  // double nEvents_se = hcent_se->Integral(4,9);
  double nEvents_se = hcent_se->Integral(centL,centH);

  TH3F* h2sig = (TH3F*)f1->Get(Form("%s%d", histname.Data(), icut))->Clone("hptH3Lmass_sig");
  h2sig->SetDirectory(0);
  TH1F* hsig = (TH1F*)h2sig->ProjectionY("hsig", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2sig->GetZaxis()->FindBin(lowy+1e-6), h2sig->GetZaxis()->FindBin(highy-1e-6));
  hsig->SetDirectory(0);

  TH3F* h2bk = (TH3F*)f2->Get(Form("%s%d", histname.Data(),icut))->Clone("hptH3Lmass_bk");
  h2bk->SetDirectory(0);
  TH1F* hbk = (TH1F*)h2bk->ProjectionY("hbk", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2bk->GetZaxis()->FindBin(lowy+1e-6), h2bk->GetZaxis()->FindBin(highy-1e-6));
  hbk->SetDirectory(0);

  double sig_sb =  hsig->Integral(hsig->GetXaxis()->FindBin(3.00),  hsig->GetXaxis()->FindBin(3.02));
  double bk_sb =   hbk->Integral(hbk->GetXaxis()->FindBin(3.00),  hbk->GetXaxis()->FindBin(3.02));
  double scale = sig_sb/bk_sb;
  hbk->Scale(scale);
  // hsig->Draw();
  // hbk->Draw("same");

  TH1F* hsig_bk = (TH1F*)hsig->Clone("hsig_bk");
  hsig_bk->Add(hbk,-1);
  setHistStyle(hsig_bk, kBlue, kFullCircle, 1.5);
  hsig_bk->Rebin();
  hsig_bk->Draw();

  TF1* fit = new TF1("fit" ,"gausn(0)+pol1(3)", 2.97,3.02 );
  TF1* resfit = new TF1("resfit" ,"pol1", 2.95,3.05 );
  hsig_bk->GetXaxis()->SetRangeUser(2.97,2.985);
  hsig_bk->Fit(resfit,"R");
  fit->SetLineColor(kRed);
  double yield_bc = hsig_bk->Integral(hsig_bk->GetXaxis()->FindBin(2.98), hsig_bk->GetXaxis()->FindBin(3.));
  // double para[5]={yield_bc*hsig_bk->GetBinWidth(1)/sqrt(2*3.1415), 2.991, 0.0015,  resfit->GetParameter(0), resfit->GetParameter(1)};
  double para[5]={yield_bc*hsig_bk->GetBinWidth(1), 2.992, 0.0014,  resfit->GetParameter(0), resfit->GetParameter(1)};
  fit->SetParameters(para);
  double lowx=2.97 ,highx =3.02;
  hsig_bk->GetXaxis()->SetRangeUser(lowx,highx);
  hsig_bk->Draw("same");
  hsig_bk->Fit(fit,"R");
  resfit->SetParameter(0, fit->GetParameter(3));
  resfit->SetParameter(1, fit->GetParameter(4));
  resfit->Draw("same");
  setHistStyle(resfit, kRed-2, 9, 2.5 ,1);
  drawLine(lowx, 0, highx, 0, 1.5, 2, 1 );
  // cout<<"binwidth: "<< hsig_bk->GetBinWidth(1)<< endl;

  double sigma = fit->GetParameter(2);
  double mean = fit->GetParameter(1);
  double yield_me = fit->GetParameter(0)/hsig_bk->GetBinWidth(1);
  double yield_counts = hsig_bk->Integral(hsig_bk->GetXaxis()->FindBin(mean-2.5*sigma), hsig_bk->GetXaxis()->FindBin(mean+2.5*sigma));
  double bk_counts = hbk->Integral(hbk->GetXaxis()->FindBin(mean-2.5*sigma), hbk->GetXaxis()->FindBin(mean+2.5*sigma));
  double sp_counts = hsig->Integral(hbk->GetXaxis()->FindBin(mean-2.5*sigma), hbk->GetXaxis()->FindBin(mean+2.5*sigma));
  // double significance = yield_counts/sqrt(yield_counts+bk_counts);
  double significance = yield_me/sqrt(sp_counts);
  double s_me = yield_me/(fit->GetParError(0)/hsig_bk->GetBinWidth(1));
  err = (fit->GetParError(0)/hsig_bk->GetBinWidth(1));

  TLegend* leg = new TLegend( 0.72, 0.68 ,0.9,0.9 );
  // leg->AddEntry(hbk, "mix-event(ME)","pl");
  // leg->AddEntry(hrt, "rotate d(RT) (scale)","pl");
  leg->AddEntry(hsig, "SE","pl");
  leg->AddEntry(hsig_bk, "SE-ME","pl");
  leg->Draw();
  drawLatex( 0.2,0.82,Form("ME/SE=%0.2f", 1./scale), 0.055);
  drawLatex( 0.2,0.75,Form("Yield=%0.2f", yield_me), 0.055);
  drawLatex( 0.2,0.68,Form("#sigma=%0.2f MeV", sigma*1000.), 0.055);
  drawLatex( 0.2,0.61,Form("nEvents=%0.0f M", nEvents_se/1e6), 0.055);
  drawLatex( 0.2,0.54,Form("S/#sqrt{S+B}=%0.0f", significance), 0.055);
  // drawLatex( 0.2,0.47,Form("S/#DeltaS=%0.0f (ME)", s_me), 0.055);
  // drawLatex( 0.2,0.4,Form("S/#DeltaS=%0.0f (RT)", s_rt), 0.055);
  drawLatex( 0.2,0.4,Form("Mean=%0.3f", mean), 0.055);
  drawLatex( 0.62,0.61,Form("%0.2f<y<%0.2f",lowy, highy ), 0.055);
  drawLatex( 0.62,0.54,Form("%0.1f<p_{T}<%0.1f GeV/c^{2}", lowpt, highpt), 0.055);
  // drawLatex( 0.62,0.47,Form("5-40%s", "%"), 0.055);
  drawLatex( 0.62,0.47,Form("0-80%s", "%"), 0.055);
  // drawBox( 2.97, hsig_bk->GetMinimum(),2.98, hsig_bk->GetMaximum()*0.5 , kBlue-9, 1001 , 1,0.3); 
  // drawBox( 2.989, hsig_bk->GetMinimum(),2.995, hsig_bk->GetMaximum()*0.5 , kOrange-4, 1001,1,0.3); 
  // drawBox( 3., hsig_bk->GetMinimum(),3.02, hsig_bk->GetMaximum()*0.5 , kBlue-9, 1001,1,0.3 ); 

  addpdf(pdf);

  return (yield_counts>20 && sigma<0.005 && significance>3.5)?yield_me:0;
}
void drawCompLaH3L()
{
  SetsPhenixStyle();
  /* TFile* fH3L= new TFile("fout_H3L3b_phase.root"); */
  // TFile* fH3L= new TFile("fout_H3L3b_quasi.root");
  // TFile* fH3L= new TFile("fout_H3L_MC_0050_015pt.root");
  // TFile* fH3L= new TFile("samecut/fout_Lambda_MC_Cuts_0050_015pt_sys_mixpid3_corr.root");
  // TFile* fH3L= new TFile("fout_Lambda_MC_Cuts_0050_015pt_sys_mixpid3.root");
  TFile* fH3L= new TFile("fout_Lambda_MC_Cuts_0050_015pt_sys_sepid3.root");
  /* TFile* fLa = new TFile("fout_Lambda.root"); */
  // TFile* fLa = new TFile("fout_Lambda_wLaD.root");
  // TFile* fLa = new TFile("fout_Lambda_MC_Cuts_0050_015pt.root");
  // TFile* fLa = new TFile("fout_Lambda_MC_Cuts_0050_015pt_sys.root");
  TFile* fLa = new TFile("fout_Lambda_MC_Cuts_0050_015pt_sys_pid3.root");

  TCanvas* c = new TCanvas("c","c");
  /* TPDF* pdf = new TPDF("Lambda_H3L_phase.pdf"); */
  /* TPDF* pdf = new TPDF("Lambda_H3L_quasi.pdf"); */
  TPDF* pdf = new TPDF("Lambda_ME.pdf");
  pdf->Off();

  // projAndComp("hptppimass", fH3L, fLa, c,pdf,"p#pi Mass (GeV/c^{2})", "plhist",1 ,"H3L quasi" ,"#Lambda+d","p#pi");
  // projAndComp("hptppichi2ndf", fH3L, fLa, c,pdf,"p#pi #chi^{2}_{NDF}","" ,1 , "H3L quasi" ,"#Lambda+d","p#pi");
  // projAndComp("hptppichi2prim", fH3L, fLa, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "H3L quasi" ,"#Lambda+d","p#pi");
  // projAndComp("hptppil", fH3L, fLa, c,pdf,"(p#pi) l" ,"",1 , "H3L quasi" ,"#Lambda+d","p#pi");
  // projAndComp("hptppildl", fH3L, fLa, c,pdf,"(p#pi) l/#deltal" ,"",1 , "H3L quasi" ,"#Lambda+d","p#pi");
  // projAndComp("hptpichi2prim", fH3L, fLa, c,pdf,"#pi #chi^{2}_{prim}" ,"",1 , "H3L quasi" ,"#Lambda+d","p#pi");
  // projAndComp("hptpchi2prim", fH3L, fLa, c,pdf,"p #chi^{2}_{prim}" ,"",1 , "H3L quasi" ,"#Lambda+d","p#pi");
  // projAndComp("hptpidca", fH3L, fLa, c,pdf,"#pi DCA","" ,2, "H3L quasi" ,"#Lambda+d","p#pi");
  // projAndComp("hptpdca", fH3L, fLa, c,pdf , "p DCA","" ,2, "H3L quasi" ,"#Lambda+d", "p#pi");
  // projAndComp("hptsumdca", fH3L, fLa, c,pdf , "p+pi DCA","" ,2, "H3L quasi" ,"#Lambda+d", "p#pi");

  // projAndComp("hptH3L_l", fH3L, fLa, c,pdf,"l", "p",4 , "#Lambda+d ME", "#Lambda+d SE" );
  projAndComp("hptH3L_ldl", fH3L, fLa, c,pdf,"l/dl", "p",4 , "#Lambda+d ME", "#Lambda+d SE");
  projAndComp("hptH3L_dchi2prim", fH3L, fLa, c,pdf,"d chi2primary", "p",4 );
  // // projAndComp("hptH3L_pchi2prim", fH3L, fLa, c,pdf,"p chi2primary", "p",4 );
  // // projAndComp("hptH3L_pichi2prim", fH3L, fLa, c,pdf,"#pi chi2primary", "p",4);
  projAndComp("hptH3L_dDca", fH3L, fLa, c,pdf,"d DCA", "p",1);
  // // projAndComp("hptH3L_dpDca", fH3L, fLa, c,pdf,"dp pair DCA", "p",4);
  // projAndComp("hptH3L_chi2topo", fH3L, fLa, c,pdf,"chi2topo", "p",4, "#Lambda+d C(K*)-1", "#Lambda+d no wt");
  // projAndComp("hptH3L_chi2ndf", fH3L, fLa, c,pdf,"chi2NDF", "p",4, "#Lambda+d C(K*)-1", "#Lambda+d no wt");
  // // projAndComp("hptH3L_ppichi2prim", fH3L, fLa, c,pdf,"p#pi chi2primary", "p",1);
  //
  // projAndComp("hptH3L_lSig", fH3L, fLa, c,pdf,"l (2.989<M<2.995)", "p",8 , "#Lambda+d ME", "#Lambda+d SE" );
  // projAndComp("hptH3L_ldlSig", fH3L, fLa, c,pdf,"l/dl (2.989<M<2.995)", "p",4 , "#Lambda+d ME", "#Lambda+d SE");
  // projAndComp("hptH3L_chi2topoSig", fH3L, fLa, c,pdf,"chi2topo (2.989<M<2.995)", "p",8, "#Lambda+d ME", "#Lambda+d SE");
  // projAndComp("hH3LptProtonPt", fH3L, fLa, c,pdf,"chi2topo (2.989<M<2.995)", "p",1, "#Lambda+d ME", "#Lambda+d SE");
  // projAndComp("hH3LptPionPt", fH3L, fLa, c,pdf,"chi2topo (2.989<M<2.995)", "p",1, "#Lambda+d ME", "#Lambda+d SE");
  projAndComp("hH3LptDeuPt", fH3L, fLa, c,pdf,"chi2topo (2.989<M<2.995)", "p",1, "#Lambda+d ME", "#Lambda+d SE");
  //
  // projAndComp("hptH3L_chi2ndfSig", fH3L, fLa, c,pdf,"chi2NDF (2.989<M<2.995)", "p",4, "#Lambda+d ME", "#Lambda+d SE");
  // projAndScaleComp("hptH3L_chi2topo", 1, fH3L, fLa, fLa,c,pdf,"chi2topo", "p",4,"#Lambda+d C(k*)-1","#Lambda+d no wt", "dp#pi" );
  // projAndScaleComp("hptH3L_chi2ndf", 1, fH3L, fLa, fLa,c,pdf,"chi2NDF", "p",8,"#Lambda+d C(k*)-1","#Lambda+d no wt", "dp#pi" );
  // projAndScaleComp("hptH3L_l", 1, fH3L, fLa, fLa,c,pdf,"l", "p",8,"#Lambda+d C(k*)-1","#Lambda+d no wt", "dp#pi" );

  pdf->On();
  pdf->Close();
}

void drawCompQuasi()
{
  SetsPhenixStyle();
  // TFile* fH3L= new TFile("fout_H3L3b_phase.root");
  TFile* fH3L= new TFile("fout_H3L_phaseMC_0050_015pt_sys.root");
  TFile* fH3L3bquasi = new TFile("fout_H3L_MC_0050_015pt_sys.root");
  // TFile* fH3L3bquasi = new TFile("fout_H3L3b_quasi.root");

  TCanvas* c = new TCanvas("c","c");
  TPDF* pdf = new TPDF("H3L_phase_vs_quasi.pdf");
  pdf->Off();

  projAndComp("hptH3L_chi2topo", fH3L, fH3L3bquasi, c,pdf,"chi2topo", "p",4,"phase","quasi", "dp#pi" );
  // projAndComp("hptH3Lmass", fH3L, fH3L3bquasi, c,pdf,"H3L Mass (GeV/c^{2})", "plhist" );
  // projAndCompH3L("hptppimass", fH3L, fH3L3bquasi, c,pdf,"p#pi Mass (GeV/c^{2})", "plhist" );
  // projAndCompH3L("hptppichi2ndf", fH3L, fH3L3bquasi, c,pdf,"p#pi #chi^{2}_{NDF}" );
  // projAndCompH3L("hptppichi2prim", fH3L, fH3L3bquasi, c,pdf,"(p#pi) #chi^{2}_{prim}" );
  // projAndCompH3L("hptppil", fH3L, fH3L3bquasi, c,pdf,"(p#pi) l" );
  // projAndCompH3L("hptppildl", fH3L, fH3L3bquasi, c,pdf,"(p#pi) l/#deltal" );
  // projAndCompH3L("hptpichi2prim", fH3L, fH3L3bquasi, c,pdf,"#pi #chi^{2}_{prim}" );
  // projAndCompH3L("hptpchi2prim", fH3L, fH3L3bquasi, c,pdf,"p #chi^{2}_{prim}" );
  // projAndComp("hptpidca", fH3L, fH3L3bquasi, c,pdf,"#pi DCA","" ,2);
  // projAndComp("hptpdca", fH3L, fH3L3bquasi, c,pdf , "p DCA","" ,2);
  // projAndComp("hptsumdca", fH3L, fH3L3bquasi, c,pdf , "p+pi DCA","" ,2);

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
void drawMixDataMB()
{
  double highpt = 2.5, lowpt = 1., lowy=-0.8, highy = -0.2;
  // double highpt = 2.5, lowpt = 1, lowy=-0.9, highy = -0.1;
  // double highpt = 4, lowpt = 0, lowy=-1.5, highy = 0;
  TString histname="hH3LMassPtY";
  // TString histname="hH3LMassPtY_5_40";
  SetsPhenixStyle();
  TCanvas* c = new TCanvas("c1","c1");
  // TPDF* pdf = new TPDF("MixEventQA_check.pdf");
  // TPDF* pdf = new TPDF("MixEventQA_Jul27.pdf");
  // TPDF* pdf = new TPDF("MixEventQA_beforeDcacut.pdf");
  TPDF* pdf = new TPDF("MixEventQA_MB.pdf");
  pdf->Off();
  gStyle->SetPalette(1);

  // TFile *f1 = TFile::Open("fout_H3L_data_SE_large.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_large.root"); 

  // TFile *f1 = TFile::Open("out_KF_test/fout_H3L_KF_test.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_KF_Jul22_part.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_KF_Jul25.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_SE_Jul25.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_SE_Jul27.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Aug09.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_SE_Aug16.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_SE_Aug20.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_SE_Aug23.root"); 
  TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02_nodca.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02_nodca_nopt.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_SE_Aug22.root"); 
  // TFile *f1 = TFile::Open("fout_KF_test.root"); 
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
  // double nEvents_se = hcent_se->Integral( 4, 9);
  double nEvents_se = hcent_se->Integral( 1, 9);

  // TFile *f2 = TFile::Open("fout_H3L_data_ME_large.root"); 
  // TFile *f2 = TFile::Open("fout_H3L_ME_Jul22_part.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Aug16.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Aug20.root"); 
  // TFile *f2 = TFile::Open("output2/fout_H3L_data_ME_Aug24.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Sep02.root"); 
  TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Sep02_nodca.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Sep02_nodca_nopt.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Aug23.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Aug22.root"); 
  // TFile *f2 = TFile::Open("fout_ME_test.root"); 
  // TFile *f2 = TFile::Open("fout_H3L_ME_Jul25.root"); 
  // TFile *f2 = TFile::Open("fout_H3L_data_SE.root"); 
  // TH2F* h2bk = (TH2F*)f2->Get("hptH3Lmass")->Clone("hptH3Lmass_ME");
  // TH3F* h2bk = (TH3F*)f2->Get("hH3LMassPtY_5_40")->Clone("hptH3Lmass_bk");
  TH3F* h2bk = (TH3F*)f2->Get(histname.Data())->Clone("hptH3Lmass_bk");
  h2bk->SetDirectory(0);
  TH1F* hbk = (TH1F*)h2bk->ProjectionY("hbk", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2bk->GetZaxis()->FindBin(lowy), h2bk->GetZaxis()->FindBin(highy));
  hbk->SetDirectory(0);

  // TFile* f3 = TFile::Open("fout_H3L_data_RT_Jul22_full.root");
  // TFile* f3 = TFile::Open("fout_H3L_data_RT_Jul25.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_RT_Jul27.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_KF_Aug20.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_SE_Aug23.root");
  // TFile* f3 = TFile::Open("output2/fout_H3L_data_KF_Aug24.root");
  TFile* f3 = TFile::Open("rootfile/fout_H3L_data_SE_Sep02_nodca.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_RT_Aug29.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_KF_Aug22.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_KF_Aug16.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_RT_Aug13.root");
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
  double sig_sb =  hsig->Integral(hsig->GetXaxis()->FindBin(3.001),  hsig->GetXaxis()->FindBin(3.02));
  // double rt_sb =  hrt->Integral(hrt->GetXaxis()->FindBin(2.97), hrt->GetXaxis()->FindBin(2.98)) + hrt->Integral(hrt->GetXaxis()->FindBin(3.0),  hrt->GetXaxis()->FindBin(3.02));
  double rt_sb = hrt->Integral(hrt->GetXaxis()->FindBin(3.001),  hrt->GetXaxis()->FindBin(3.02));
  // double rt_sb =  hrt->Integral(hrt->GetXaxis()->FindBin(2.96), hrt->GetXaxis()->FindBin(2.98)) ;
  double scale_rt = sig_sb/rt_sb;
  // sig_sb = hsig->Integral(hsig->GetXaxis()->FindBin(2.97), hsig->GetXaxis()->FindBin(2.98)) + hsig->Integral(hsig->GetXaxis()->FindBin(3.0),  hsig->GetXaxis()->FindBin(3.02));
  double bk_sb =   hbk->Integral(hbk->GetXaxis()->FindBin(3.001),  hbk->GetXaxis()->FindBin(3.02));
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
  // hrt->Draw("same");
  TLegend* leg_sig = new TLegend(0.65,0.25,0.88,0.45);
  leg_sig->AddEntry(hbk, "ME", "pl");
  // leg_sig->AddEntry(hrt, "RT", "pl");
  // leg_sig->AddEntry(hrt, "KF", "pl");
  leg_sig->AddEntry(hsig, "SE", "pl");
  leg_sig->Draw();

  drawLatex( 0.65,0.61,Form("%0.1f<y<%0.1f",lowy, highy ), 0.055);
  drawLatex( 0.65,0.54,Form("%0.1f<p_{T}<%0.1f GeV/c^{2}", lowpt, highpt), 0.055);
  // drawLatex( 0.65,0.47,Form("5-40%s", "%"), 0.055);
  drawLatex( 0.65,0.47,Form("0-50%s", "%"), 0.055);
  drawBox( 2.97, hsig->GetMinimum(),2.98, hsig->GetMaximum()*0.5 , kBlue-9, 1001 , 1,0.3); 
  // drawBox( 2.989, hsig_bk->GetMinimum(),2.995, hsig_bk->GetMaximum()*0.5 , kOrange-4, 1001,1,0.3); 
  drawBox( 3.01, hsig->GetMinimum(),3.02, hsig->GetMaximum()*0.5 , kBlue-9, 1001,1,0.3 ); 

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
  // hsig_rt->Draw("same");
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
  // hsig_rt->Fit(resfit_rt,"R");
  fit_rt->SetParameter(3, resfit_rt->GetParameter(0) );
  fit_rt->SetParameter(4, resfit_rt->GetParameter(1) );
  // hsig_rt->Fit(fit_rt,"R");
  resfit_rt->SetParameter(0, fit_rt->GetParameter(3));
  resfit_rt->SetParameter(1, fit_rt->GetParameter(4));
  // resfit_rt->Draw("same");

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
  double significance = yield_me/sqrt(sp_counts);
  // double significance = yield_rt/sqrt(sp_counts);
  double s_me = yield_me/(fit->GetParError(0)/hsig_bk->GetBinWidth(1));
  double s_rt = yield_rt/(fit_rt->GetParError(0)/hsig_rt->GetBinWidth(1));
  cout<<"significance: " <<significance << endl;

  TLegend* leg = new TLegend( 0.72, 0.68 ,0.9,0.9 );
  // leg->AddEntry(hbk, "mix-event(ME)","pl");
  // leg->AddEntry(hrt, "rotate d(RT) (scale)","pl");
  leg->AddEntry(hsig, "SE","pl");
  leg->AddEntry(hsig_bk, "SE-ME","pl");
  // leg->AddEntry(hsig_rt, "SE-KF","pl");
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
  drawBox( 3.01, hsig_bk->GetMinimum(),3.02, hsig_bk->GetMaximum()*0.5 , kBlue-9, 1001,1,0.3 ); 

  addpdf(pdf);

  hsig_bk->Draw();
  // hsig_rt->Draw("same");
  TLegend* leg2 = new TLegend( 0.2 , 0.7 ,0.4,0.9  );
  leg2->AddEntry(hsig_rt, "SE-RT","pl");
  leg2->AddEntry(hsig_bk, "SE-ME","pl");
  leg2->Draw();

  addpdf(pdf);
  //
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

  TH2F* hYield = new TH2F( "hYield", "hYield;y,pt", 5, -1., 0., 6, 0, 3);
  double xbw = hYield->GetXaxis()->GetBinWidth(1);
  double ybw = hYield->GetYaxis()->GetBinWidth(1);
  for (int ix=1; ix<=hYield->GetNbinsX(); ix++) {
    for (int iy=1; iy<=hYield->GetNbinsY();iy++){
      double xlow = hYield->GetXaxis()->GetBinLowEdge(ix); // y
      double ylow = hYield->GetYaxis()->GetBinLowEdge(iy); // pt
      double err;
      hYield->SetBinContent( ix , iy, fityield( ylow, ylow+ybw, xlow, xlow+xbw, err, f1, f2, c, pdf, 1,9) );
      hYield->SetBinError( ix , iy, err);
    }
  }
  hYield->Draw("colz text");
  hYield->SaveAs("fH3L_yield_MB.root");
  addpdf(pdf);

  TFile* fMc_ld = TFile::Open("fout_Lambda_MC_Cuts.root");
  // TFile* fMc = TFile::Open("fout_H3L_MC_0050_dcacut.root");
  TFile* fMc = TFile::Open("fout_H3L_MC_0080_010pt.root");
  TH1F* hPhase[3];
  TH1F* hPhaseCor[3];
  TH1F* hBr[3];
  TGraphErrors* g[3];
  double edge[3][4]={
    { 1.7, 2.1, 2.7, 3.1},
    { 1.2, 1.6, 2.0, 2.4},
    { 1.2, 1.8, 2.2, 2.6}
  };

  double ybin[4]={ 0, -0.25, -0.50, -0.75};
  hPhase[0] = new TH1F("hPhase0", "hPhase0;pt", 25, 1, 3.5);
  hPhase[0]=(TH1F*)hPhase[0]->Rebin(3, "hPhase0", edge[0]);
  hPhase[1] = new TH1F("hPhase1", "hPhase1;pt", 25, 1, 3.5);
  hPhase[1]=(TH1F*)hPhase[1]->Rebin(3, "hPhase1", edge[1]);
  hPhase[2] = new TH1F("hPhase2", "hPhase2;pt", 25, 1, 3.5);
  hPhase[2]=(TH1F*)hPhase[2]->Rebin(3, "hPhase2", edge[2]);
  // TFile* fMcH3L = new TFile("fMC_H3L_0050.root");
  // TFile* fRcH3L = new TFile("fout_H3L_MC_0050_dcacut.root");
  TFile* fMcH3L = new TFile("fMC_H3L_0080.root");
  TFile* fRcH3L = new TFile("fout_H3L_MC_0080_010pt.root");
  TH3F* h3Mc = (TH3F*)fMcH3L->Get("hH3LMassPtY")->Clone("h3Mc");
  h3Mc->SetDirectory(0);
  TH3F* h3Rc = (TH3F*)fRcH3L->Get("hH3LMassPtY")->Clone("h3Rc");
  h3Rc->SetDirectory(0);

  h3Mc->Sumw2();
  h3Rc->Sumw2();
  TH2F* h2MC = (TH2F*)h3Mc->Project3D("xz");
  TH2F* h2Rc = (TH2F*)h3Rc->Project3D("xz");

  TH2F* h2Eff = (TH2F*)h2Rc->Clone("h2Eff");
  TH2F* h2temp = (TH2F*)h2MC->Clone("h2temp");
  h2temp->RebinY(10);
  h2temp->RebinX(10);
  h2Eff->RebinY(10);
  h2Eff->RebinX(10);
  h2Eff->Divide(h2temp);
  h2Eff->Draw("colz text");
  h2Eff->GetYaxis()->SetRangeUser(0,4.5);
  addpdf(pdf);

  TH1F* heff[3];
  for (int ij=0;ij<3;ij++){ 
    TH1F* h1Mc = (TH1F*)h2MC->ProjectionY(Form("h1Mc%d", ij), h2MC->GetXaxis()->FindBin( ybin[ij+1]+1e-6), h2MC->GetXaxis()->FindBin(ybin[ij]-1e-6) );
    h1Mc = (TH1F*)h1Mc->Rebin( 3, Form("hMc%d",ij), edge[ij]);
    TH1F* h1Rc = (TH1F*)h2Rc->ProjectionY(Form("h1Rc%d", ij), h2Rc->GetXaxis()->FindBin(ybin[ij+1]+1e-6), h2Rc->GetXaxis()->FindBin(ybin[ij]-1e-6) );
    heff[ij] = (TH1F*)h1Rc->Rebin( 3, Form("heff%d",ij), edge[ij]);
    heff[ij]->Divide(h1Mc);
    heff[ij]->GetYaxis()->SetTitle("Eff.");
    heff[ij]->Draw();
    addpdf(pdf);
  }

  double dy=0.25;
  double dpt,pt;
  
  TFile* f2b = TFile::Open("h_h3l_corr_yield.root_cut00002_cent0_yuehang.root");

  double ny = 3;
  for (int ij=0; ij<ny;ij++){
    double npt =3;
    if (ij==0) npt =2; 
    for (int ipt=0; ipt<npt;ipt++){
      double err;
      dpt = edge[ij][ipt+1]-edge[ij][ipt];
      pt = 0.5*(edge[ij][ipt+1]+edge[ij][ipt]);
      double yield = fityield( edge[ij][ipt], edge[ij][ipt+1], ybin[ij+1], ybin[ij], err, f1, f2, c, pdf);
      // cout << edge[ij][ipt]<<" "<<edge[ij][ipt+1]<<" "<<ybin[ij+1]<<" "<<ybin[ij]<<" "<< yield<<" "<<err <<endl;
      if (ij==0 && ipt==2) { yield=0; err=0;}
      hPhase[ij]->SetBinContent( ipt+1, yield/dy/dpt/pt/2./3.1415926/nEvents_se);
      hPhase[ij]->SetBinError( ipt+1, err/dy/dpt/pt/2./3.1415926/nEvents_se );
    }
    hPhase[ij]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPhase[ij]->GetYaxis()->SetTitle("Raw d^{2}N/(N_{evt}2#pip_{T}dydp_{T})");
    hPhase[ij]->Draw();
    hPhase[ij]->SetDirectory(0);
    addpdf(pdf); 

    hPhaseCor[ij]=(TH1F*)hPhase[ij]->Clone(Form("hYieldCor_%d",ij));
    hPhaseCor[ij]->SetDirectory(0);
    for (int ipt=0;ipt<npt;ipt++)
    {
      double y3b = hPhaseCor[ij]->GetBinContent(ipt+1);
      double y3berr = hPhaseCor[ij]->GetBinError(ipt+1)/y3b;
      double eff = heff[ij]->GetBinContent(ipt+1);
      double efferr = heff[ij]->GetBinError(ipt+1)/eff;
      double perr;
      double purity1 = calpurity("h3H3L_chi2ndf",0, 3.5, edge[ij][ipt+1], edge[ij][ipt], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 5, perr, 0);
      // double purity2 = calpurity("h3H3L_chi2topo",0, 3., edge[ij][ipt+1], edge[ij][ipt], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
      // double purity = (purity1+purity2)*0.5;
      // perr = fabs(purity2-purity1);
      // perr = perr/purity;
      double purity = purity1;
      purity=1;
      perr = 0;
      double yield_cor = y3b/eff*purity;
      double err = sqrt(y3berr*y3berr + efferr*efferr + perr*perr )*yield_cor;
      hPhaseCor[ij]->SetBinContent(ipt+1,yield_cor );
      hPhaseCor[ij]->SetBinError(ipt+1, err);
      cout <<"pt "<<hPhaseCor[ij]->GetBinCenter(ipt+1) << "purity "<<purity <<" eff " <<eff <<" yield " << yield_cor << " err "<< err<< " s="<<yield_cor/err<<endl;
    }
    hPhaseCor[ij]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPhaseCor[ij]->GetYaxis()->SetTitle("d^{2}N/(N_{evt}2#pip_{T}dydp_{T})");
    hPhaseCor[ij]->Draw();
    hPhaseCor[ij]->SetDirectory(0);
    addpdf(pdf);

    hBr[ij]=(TH1F*)hPhaseCor[ij]->Clone(Form("hBr_%d",ij));
    hBr[ij]->SetDirectory(0);
    g[ij] = (TGraphErrors*)f2b->Get(Form("t_sgct1_corr_yield[0][%d]", ij));
    g[ij]->Draw();
    for (int ipt=0;ipt<npt;ipt++)
    {
      double y3b = hBr[ij]->GetBinContent(ipt+1);
      double y3berr = hBr[ij]->GetBinError(ipt+1);
      double y2b = g[ij]->GetPointY(ipt);
      double y2bpt = g[ij]->GetPointX(ipt);
      double y2berr = g[ij]->GetErrorY(ipt);
      if (ij==1) {y2berr=y2berr*10; y2b=y2b*10.;}
      if (ij==2) {y2berr=y2berr*100; y2b=y2b*100.;}
      // cout <<"y2b significance:"<<y2b/y2berr<<" yield "<<y2b << endl;
      hBr[ij]->SetBinContent(ipt+1, y2b/(y2b+y3b) );
      hBr[ij]->SetBinError(ipt+1, sqrt((y2b*y2b)/pow((y2b+y3b), 4)*(y3berr*y3berr) + (y3b*y3b)/pow((y2b+y3b), 4)*(y2berr*y2berr)));
      if (y2b/y2berr<2) {
         hBr[ij]->SetBinContent(ipt+1, 0);
         hBr[ij]->SetBinError(ipt+1, 0);
      } 
      cout <<"bincenter: " << hBr[ij]->GetBinCenter(ipt+1) << " content:"<< hBr[ij]->GetBinContent(ipt+1)<<" err:"<<hBr[ij]->GetBinError(ipt+1)<< "y2b significance:"<<y2b/y2berr<<" yield "<<y2b << " y2bpt "<< y2bpt<< " y3b yield "<< y3b<<endl;
    } 
    hBr[ij]->Draw();
    hBr[ij]->GetYaxis()->SetTitle("2body/(2body + 3body)");
    drawLatex( 0.65,0.88 , Form("%0.2f<y<%0.2f", ybin[ij+1],ybin[ij]) , 0.055);
    addpdf(pdf);
  }

  TLegend* legb = new TLegend(0.2, 0.2, 0.4, 0.4);
  TH1F* htemp = new TH1F("htemp","htemp;p_{T} (GeV/c); 2-body/(2-body+3-body)",  5, 1, 3);
  htemp->GetYaxis()->SetRangeUser( 0, 0.7);
  htemp->Draw();
  int color[]={kRed, kGreen, kBlue };
  for (int i=0;i<3;i++)
  {
    hBr[i]->Draw("same");
    hBr[i]->SetLineColor(color[i]);
    hBr[i]->SetMarkerColor(color[i]);
    legb->AddEntry( hBr[i], Form("%0.2f<y<%0.2f", ybin[i+1], ybin[i]), "pl");
  }
  legb->Draw();
  addpdf(pdf);

  TFile * fout = new TFile("fout.root","recreate");
  fout->cd();
  for (int i=0;i<3;i++){
    hPhase[i]->Write();
    hPhaseCor[i]->Write();
    hBr[i]->Write();
    heff[i]->Write();
  }
  fout->Close();

  pdf->On();
  pdf->Close();
  return;

  c = new TCanvas( "c","c", 800,800);
  projAndComp("hptH3L_lSBL", f1, f2, c,pdf,"l", "p",4,"SE_SBL","ME_SBL", "dp#pi", scale);
  projAndComp("hptH3L_lSBR", f1, f2, c,pdf,"l", "p",4,"SE_SBR","ME_SBR", "dp#pi", scale);
  projAndComp("hptH3L_ldlSBL", f1, f2, c,pdf,"l/dl", "p",4,"SE_SBL","ME_SBL", "dp#pi", scale);
  projAndComp("hptH3L_ldlSBR", f1, f2, c,pdf,"l/dl", "p",4,"SE_SBR","ME_SBR", "dp#pi", scale);
  projAndComp("hptH3L_chi2ndfSBL", f1, f2, c,pdf,"chi2NDF", "p",4,"SE_SBL","ME_SBL", "dp#pi", scale);
  projAndComp("hptH3L_chi2ndfSBR", f1, f2, c,pdf,"chi2NDF", "p",4,"SE_SBR","ME_SBR", "dp#pi", scale);
  projAndComp("hptH3L_chi2topoSBL", f1, f2, c,pdf,"chi2topo", "p",4,"SE_SBL","ME_SBL", "dp#pi" , scale);
  projAndComp("hptH3L_chi2topoSBR", f1, f2, c,pdf,"chi2topo", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  projAndComp("hptH3L_dchi2primSBL", f1, f2, c,pdf,"d chi2primary", "p",4,"SE_SBL","ME_SBL", "dp#pi", scale);
  projAndComp("hptH3L_dchi2primSBR", f1, f2, c,pdf,"d chi2primary", "p",4,"SE_SBR","ME_SBR", "dp#pi", scale);
  projAndComp("hptH3L_pchi2primSBL", f1, f2, c,pdf,"p chi2primary", "p",4,"SE_SBL","ME_SBL", "dp#pi", scale);
  projAndComp("hptH3L_pchi2primSBR", f1, f2, c,pdf,"p chi2primary", "p",4,"SE_SBR","ME_SBR", "dp#pi", scale);
  projAndComp("hptH3L_pichi2primSBL", f1, f2, c,pdf,"#pi chi2primary", "p",4,"SE_SBL","ME_SBL", "dp#pi", scale);
  projAndComp("hptH3L_pichi2primSBR", f1, f2, c,pdf,"#pi chi2primary", "p",4,"SE_SBR","ME_SBR", "dp#pi", scale);
  projAndComp("hptH3L_dDcaSBL", f1, f2, c,pdf,"d DCA", "p",1,"SE_SBL","ME_SBL","dp#pi"  , scale);
  projAndComp("hptH3L_dDcaSBR", f1, f2, c,pdf,"d DCA", "p",1,"SE_SBR","ME_SBR","dp#pi"  , scale);
  projAndComp("hptH3L_piDcaSBL", f1, f2, c,pdf,"#pi DCA", "p",4,"SE_SBL","ME_SBL", "dp#pi" , scale);
  projAndComp("hptH3L_piDcaSBR", f1, f2, c,pdf,"#pi DCA", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  projAndComp("hptH3L_pDcaSBL", f1, f2, c,pdf,"#p DCA", "p",4,"SE_SBL","ME_SBL", "dp#pi" , scale);
  projAndComp("hptH3L_pDcaSBR", f1, f2, c,pdf,"#p DCA", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  projAndComp("hptH3L_dpDcaSBL", f1, f2, c,pdf,"dp pair DCA", "p",4,"SE_SBL","ME_SBL", "dp#pi" , scale);
  projAndComp("hptH3L_dpDcaSBR", f1, f2, c,pdf,"dp pair DCA", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  // projAndComp("hptppimassSBL", f1, f2, c,pdf,"p#pi Mass", "p",4,"SE_SBL","ME_SBL", "dp#pi" );
  // projAndComp("hptH3L_ppimassSBR", f1, f2, c,pdf,"p#pi Mass", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  // projAndComp("hptH3L_ppichi2ndfSBR", f1, f2, c,pdf,"p#pi Mass", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  // projAndComp("hptH3L_ppichi2ndfSBL", f1, f2, c,pdf,"p#pi chi2ndf", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  // projAndComp("hptH3L_ppichi2primSBL", f1, f2, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "SE" ,"BK","p#pi",scale);

  // projAndComp("hptH3L_lSBR", f1, f3, c,pdf,"l", "p",4,"SE_SBR","RT_SBR", "dp#pi");
  // projAndComp("hptH3L_ldlSBR", f1, f3, c,pdf,"l/dl", "p",4,"SE_SBR","RT_SBR", "dp#pi");
  // projAndComp("hptH3L_chi2ndfSBR", f1, f3, c,pdf,"chi2NDF", "p",4,"SE_SBR","RT_SBR", "dp#pi");
  // projAndComp("hptH3L_chi2topoSBR", f1, f3, c,pdf,"chi2topo", "p",4,"SE_SBR","RT_SBR", "dp#pi" );
  // projAndComp("hptH3L_dchi2primSBR", f1, f3, c,pdf,"d chi2primary", "p",4,"SE_SBR","RT_SBR", "dp#pi");
  // projAndComp("hptH3L_pchi2primSBR", f1, f3, c,pdf,"p chi2primary", "p",4,"SE_SBR","RT_SBR", "dp#pi");
  // projAndComp("hptH3L_pichi2primSBR", f1, f3, c,pdf,"#pi chi2primary", "p",4,"SE_SBR","RT_SBR", "dp#pi");
  // projAndComp("hptH3L_dDcaSBR", f1, f3, c,pdf,"d DCA", "p",1,"SE_SBR","RT_SBR","dp#pi"  );
  // projAndComp("hptH3L_piDcaSBR", f1, f3, c,pdf,"#pi DCA", "p",4,"SE_SBR","RT_SBR", "dp#pi" );
  // projAndComp("hptH3L_pDcaSBR", f1, f3, c,pdf,"#p DCA", "p",4,"SE_SBR","RT_SBR", "dp#pi" );
  // projAndComp("hptH3L_dpDcaSBR", f1, f3, c,pdf,"dp pair DCA", "p",4,"SE_SBR","RT_SBR", "dp#pi" );
  // projAndComp("hptH3L_ppimassSBR", f1, f3, c,pdf,"p#pi Mass (GeV/c^{2})", "plhist",1 ,"SE" ,"BK","p#pi");
  // // projAndComp("hptppichi2ndfSBR", f1, f2, c,pdf,"p#pi #chi^{2}_{NDF}","" ,1 , "SE" ,"BK","p#pi");
  // // projAndComp("hptppichi2primSBR", f1, f2, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "SE" ,"BK","p#pi");
  // // projAndComp("hptppilSBR", f1, f2, c,pdf,"(p#pi) l" ,"",1 , "SE" ,"BK","p#pi");
  // // projAndComp("hptppildlSBR", f1, f2, c,pdf,"(p#pi) l/#deltal" ,"",1 , "SE" ,"BK","p#pi");

  // projAndComp("hptH3L_lSBL", f1, f3, c,pdf,"l", "p",4,"SE_SBL","RT_SBL", "dp#pi");
  // projAndComp("hptH3L_ldlSBL", f1, f3, c,pdf,"l/dl", "p",4,"SE_SBL","RT_SBL", "dp#pi");
  // projAndComp("hptH3L_chi2ndfSBL", f1, f3, c,pdf,"chi2NDF", "p",4,"SE_SBL","RT_SBL", "dp#pi");
  // projAndComp("hptH3L_chi2topoSBL", f1, f3, c,pdf,"chi2topo", "p",4,"SE_SBL","RT_SBL", "dp#pi" );
  // projAndComp("hptH3L_dchi2primSBL", f1, f3, c,pdf,"d chi2primary", "p",4,"SE_SBL","RT_SBL", "dp#pi");
  // projAndComp("hptH3L_pchi2primSBL", f1, f3, c,pdf,"p chi2primary", "p",4,"SE_SBL","RT_SBL", "dp#pi");
  // projAndComp("hptH3L_pichi2primSBL", f1, f3, c,pdf,"#pi chi2primary", "p",4,"SE_SBL","RT_SBL", "dp#pi");
  // projAndComp("hptH3L_dDcaSBL", f1, f3, c,pdf,"d DCA", "p",1,"SE_SBL","RT_SBL","dp#pi"  );
  // projAndComp("hptH3L_piDcaSBL", f1, f3, c,pdf,"#pi DCA", "p",4,"SE_SBL","RT_SBL", "dp#pi" );
  // projAndComp("hptH3L_pDcaSBL", f1, f3, c,pdf,"#p DCA", "p",4,"SE_SBL","RT_SBL", "dp#pi" );
  // projAndComp("hptH3L_dpDcaSBL", f1, f3, c,pdf,"dp pair DCA", "p",4,"SE_SBL","RT_SBL", "dp#pi" );
  // // projAndComp("hptppimassSBL", f1, f3, c,pdf,"p#pi Mass (GeV/c^{2})", "plhist",1 ,"SE" ,"BK","p#pi");
  // projAndComp("hptH3L_ppimassSBL", f1, f3, c,pdf,"p#pi Mass (GeV/c^{2})", "plhist",1 ,"SE" ,"BK","dp#pi");
  // // projAndComp("hptppichi2ndfSBL", f1, f2, c,pdf,"p#pi #chi^{2}_{NDF}","" ,1 , "SE" ,"BK","p#pi");
  // // projAndComp("hptppichi2primSBL", f1, f2, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "SE" ,"BK","p#pi");
  // // projAndComp("hptppilSBL", f1, f2, c,pdf,"(p#pi) l" ,"",1 , "SE" ,"BK","p#pi");
  // // projAndComp("hptppildlSBL", f1, f2, c,pdf,"(p#pi) l/#deltal" ,"",1 , "SE" ,"BK","p#pi");

  // projAndScaleComp("hptH3L_l", scale,f1, f2, f2, c,pdf,"l", "p",1,"SE-ME","ME", "dp#pi");
  // projAndScaleComp("hptH3L_ppimass", scale, f1, f2, f2, c,pdf,"M(p#pi)", "p",1,"SE-ME","ME", "p#pi");
  // // projAndScaleCompNoRatio("hptppimass", scale, f1, f2, f2, c,pdf,"M(p#pi)", "p",1,"SE-ME","ME", "p#pi");
  // projAndScaleComp("hptH3L_ldl", scale,f1, f2, f2, c,pdf,"l/dl", "p",1,"SE-ME","ME", "dp#pi");
  // projAndScaleComp("hptH3L_chi2ndf", scale,f1, f2,f2, c,pdf,"chi2NDF", "p",1,"SE-ME","ME", "dp#pi");
  projAndScaleComp("hptH3L_chi2topo", scale,f1, f2, f2,c,pdf,"chi2topo", "p",1,"SE-ME","ME", "dp#pi" );
  // projAndScaleComp("hptH3L_dchi2prim", scale,f1, f2, f2,c,pdf,"d chi2primary", "p",1,"SE-ME","ME", "dp#pi");
  // projAndScaleComp("hptH3L_pchi2prim", scale,f1, f2, f2,c,pdf,"p chi2primary", "p",1,"SE-ME","ME", "dp#pi");
  // projAndScaleComp("hptH3L_pichi2prim", scale,f1, f2, f2,c,pdf,"#pi chi2primary", "p",1,"SE-ME","ME", "dp#pi");
  // projAndScaleComp("hptH3L_dDca", scale,f1, f2, f2,c,pdf,"d DCA", "p",1,"SE-ME","ME","dp#pi"  );
  // projAndScaleComp("hptH3L_piDca", scale,f1, f2, f2,c,pdf,"#pi DCA", "p",1,"SE-ME","ME", "dp#pi" );
  // projAndScaleComp("hptH3L_pDca", scale,f1, f2, f2,c,pdf,"#p DCA", "p",1,"SE-ME","ME", "dp#pi" );
  // projAndScaleComp("hptH3L_dpDca", scale,f1, f2, f2,c,pdf,"dp pair DCA", "p",1,"SE-ME","ME", "dp#pi" );
  // // projAndScaleComp("hptppichi2prim",scale, f1, f2,f2, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "SE-ME" ,"ME","p#pi");

  // projAndScaleComp("hptH3L_l", scale_rt,f1, f3, f3, c,pdf,"l", "p",1,"SE-RT","RT", "dp#pi");
  // projAndScaleComp("hptppimass", scale_rt, f1, f3, f3, c,pdf,"M(p#pi)", "p",1,"SE-RT","RT", "p#pi");
  // projAndScaleComp("hptH3L_ldl", scale_rt,f1, f3, f3, c,pdf,"l/dl", "p",1,"SE-RT","RT", "dp#pi");
  // projAndScaleComp("hptH3L_chi2ndf", scale_rt,f1, f3,f3, c,pdf,"chi2NDF", "p",1,"SE-RT","RT", "dp#pi");
  // projAndScaleComp("hptH3L_chi2topo", scale_rt,f1, f3, f3,c,pdf,"chi2topo", "p",1,"SE-RT","RT", "dp#pi" );
  // projAndScaleComp("hptH3L_dchi2prim", scale_rt,f1, f3, f3,c,pdf,"d chi2primary", "p",1,"SE-RT","RT", "dp#pi");
  // projAndScaleComp("hptH3L_pchi2prim", scale_rt,f1, f3, f3,c,pdf,"p chi2primary", "p",1,"SE-RT","RT", "dp#pi");
  // projAndScaleComp("hptH3L_pichi2prim", scale_rt,f1, f3, f3,c,pdf,"#pi chi2primary", "p",1,"SE-RT","RT", "dp#pi");
  // projAndScaleComp("hptH3L_dDca", scale_rt,f1, f3, f3,c,pdf,"d DCA", "p",1,"SE-RT","RT","dp#pi"  );
  // projAndScaleComp("hptH3L_piDca", scale_rt,f1, f3, f3,c,pdf,"#pi DCA", "p",1,"SE-RT","RT", "dp#pi" );
  // projAndScaleComp("hptH3L_pDca", scale_rt,f1, f3, f3,c,pdf,"#p DCA", "p",1,"SE-RT","RT", "dp#pi" );
  // projAndScaleComp("hptH3L_dpDca", scale_rt,f1, f3, f3,c,pdf,"dp pair DCA", "p",1,"SE-RT","RT", "dp#pi" );
  // projAndScaleComp("hptppichi2prim",scale_rt, f1, f3,f3, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "SE-RT" ,"RT","p#pi");

  cout <<"compare with MC" <<endl;
  // projAndScaleComp("hptH3L_l", scale,f1, f2, fMc, c,pdf,"l", "p",4,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_ppimass", scale, f1, f2, fMc, c,pdf,"M(p#pi)", "p",4,"SE-ME","MC", "p#pi");
  // // projAndScaleCompNoRatio("hptH3L_ppimass", scale, f1, f2, fMc, c,pdf,"M(p#pi)", "p",4,"SE-ME","MC", "p#pi");
  // projAndScaleComp("hptH3L_ldl", scale,f1, f2, fMc, c,pdf,"l/dl", "p",4,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_chi2ndf", scale,f1, f2,fMc, c,pdf,"chi2NDF","p",2,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_chi2topo", scale,f1, f2, fMc,c,pdf,"chi2topo","p",2,"SE-ME","MC", "dp#pi" );
  // // projAndScaleComp("hptH3L_chi2", scale,f1, f2, fMc,c,pdf,"chi2topo","p",2,"SE-ME","MC", "dp#pi" );
  // projAndScaleComp("hptH3L_dchi2prim", scale,f1, f2, fMc,c,pdf,"d chi2primary", "p",5,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_pchi2prim", scale,f1, f2, fMc,c,pdf,"p chi2primary", "p",2,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_pichi2prim", scale,f1, f2, fMc,c,pdf,"#pi chi2primary", "p",2,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_dDca", scale,f1, f2, fMc,c,pdf,"d DCA", "p",2,"SE-ME","MC","dp#pi"  );
  // projAndScaleComp("hptH3L_piDca", scale,f1, f2, fMc,c,pdf,"#pi DCA", "p",2,"SE-ME","MC", "dp#pi" );
  // projAndScaleComp("hptH3L_pDca", scale,f1, f2, fMc,c,pdf,"#p DCA", "p",2,"SE-ME","MC", "dp#pi" );
  // projAndScaleComp("hptH3L_dpDca", scale,f1, f2, fMc,c,pdf,"dp pair DCA", "p",2,"SE-ME","MC", "dp#pi" );
  // // projAndScaleComp("hptppichi2prim",scale, f1, f2,fMc, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",2 , "SE-ME" ,"MC","p#pi");
  //

  // TFile* fMc_ph = TFile::Open("fout_H3L_phase_MC.root");
  // //  projAndScaleComp("hptH3L_l", scale,f1, f2, fMc_ph, c,pdf,"l", "p",4,"SE-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptppimass", scale, f1, f2, fMc_ph, c,pdf,"M(p#pi)", "p",4,"SE-ME","phase_MC", "p#pi");
  // projAndScaleComp("hptH3L_ldl", scale,f1, f2, fMc_ph, c,pdf,"l/dl", "p",4,"SE-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptH3L_chi2ndf", scale,f1, f2,fMc_ph, c,pdf,"chi2NDF","p",2,"SE-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptH3L_chi2topo", scale,f1, f2, fMc_ph,c,pdf,"chi2topo","p",2,"SE-ME","phase_MC", "dp#pi" );
  // // projAndScaleComp("hptH3L_chi2", scale,f1, f2, fMc_ph,c,pdf,"chi2topo","p",2,"SE-ME","phase_MC", "dp#pi" );
  // projAndScaleComp("hptH3L_dchi2prim", scale,f1, f2, fMc_ph,c,pdf,"d chi2primary", "p",2,"SE-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptH3L_pchi2prim", scale,f1, f2, fMc_ph,c,pdf,"p chi2primary", "p",2,"SE-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptH3L_pichi2prim", scale,f1, f2, fMc_ph,c,pdf,"#pi chi2primary", "p",2,"SE-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptH3L_dDca", scale,f1, f2, fMc_ph,c,pdf,"d DCA", "p",2,"SE-ME","phase_MC","dp#pi"  );
  // projAndScaleComp("hptH3L_piDca", scale,f1, f2, fMc_ph,c,pdf,"#pi DCA", "p",2,"SE-ME","phase_MC", "dp#pi" );
  // projAndScaleComp("hptH3L_pDca", scale,f1, f2, fMc_ph,c,pdf,"#p DCA", "p",2,"SE-ME","phase_MC", "dp#pi" );
  // projAndScaleComp("hptH3L_dpDca", scale,f1, f2, fMc_ph,c,pdf,"dp pair DCA", "p",2,"SE-ME","phase_MC", "dp#pi" );
  // projAndScaleComp("hptppichi2prim",scale, f1, f2,fMc_ph, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",2 , "SE-ME" ,"phase_MC","p#pi");

  cout <<"compare with Lambda MC" <<endl;
  projAndScaleComp("hptH3L_l", scale,f1, f2, fMc_ld, c,pdf,"l", "p",4,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_l", scale,f1, f2, fMc, c,pdf,"l", "p",4,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_ppimass", scale, f1, f2, fMc_ld, c,pdf,"M(p#pi)", "p",4,"SE-ME","MC_#Lambda", "p#pi",0,1);
  // projAndScaleCompNoRatio("hptppimass", scale, f1, f2, fMc_ld, c,pdf,"M(p#pi)", "p",4,"SE-ME","MC_#Lambda", "p#pi",0);
  // projAndScaleCompNoRatio("hptH3L_ppimass", scale, f1, f2, fMc_ld, c,pdf,"M(p#pi)", "p",4,"SE-ME","MC_#Lambda", "p#pi",0);
  projAndScaleComp("hptH3L_ldl", scale,f1, f2, fMc_ld, c,pdf,"l/dl", "p",4,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_ldl", scale,f1, f2, fMc, c,pdf,"l/dl", "p",4,"SE-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_chi2ndf", scale,f1, f2,fMc_ld, c,pdf,"chi2NDF","p",2,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_chi2ndf", scale,f1, f2,fMc, c,pdf,"chi2NDF","p",2,"SE-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_chi2topo", scale,f1, f2, fMc_ld,c,pdf,"chi2topo","p",2,"SE-ME","MC_#Lambda", "dp#pi" ,0,1);
  projAndScaleComp("hptH3L_chi2topo", scale,f1, f2, fMc,c,pdf,"chi2topo","p",2,"SE-ME","MC", "dp#pi" );
  // projAndScaleComp("hptH3L_chi2", scale,f1, f2, fMc_ld,c,pdf,"chi2topo","p",2,"SE-ME","MC_#Lambda", "dp#pi" ,0,1);
  projAndScaleComp("hptH3L_dchi2prim", scale,f1, f2, fMc_ld,c,pdf,"d chi2primary", "p",5,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_dchi2prim", scale,f1, f2, fMc,c,pdf,"d chi2primary", "p",5,"SE-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_pchi2prim", scale,f1, f2, fMc_ld,c,pdf,"p chi2primary", "p",2,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_pchi2prim", scale,f1, f2, fMc,c,pdf,"p chi2primary", "p",2,"SE-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_pichi2prim", scale,f1, f2, fMc_ld,c,pdf,"#pi chi2primary", "p",2,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_pichi2prim", scale,f1, f2, fMc,c,pdf,"#pi chi2primary", "p",2,"SE-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_dDca", scale,f1, f2, fMc_ld,c,pdf,"d DCA", "p",2,"SE-ME","MC_#Lambda","dp#pi"  ,0,1);
  projAndScaleComp("hptH3L_dDca", scale,f1, f2, fMc,c,pdf,"d DCA", "p",2,"SE-ME","MC","dp#pi"  );
  projAndScaleComp("hptH3L_piDca", scale,f1, f2, fMc_ld,c,pdf,"#pi DCA", "p",2,"SE-ME","MC_#Lambda", "dp#pi" ,0,1);
  projAndScaleComp("hptH3L_piDca", scale,f1, f2, fMc,c,pdf,"#pi DCA", "p",2,"SE-ME","MC", "dp#pi" );
  projAndScaleComp("hptH3L_pDca", scale,f1, f2, fMc_ld,c,pdf,"#p DCA", "p",2,"SE-ME","MC_#Lambda", "dp#pi" ,0,1);
  projAndScaleComp("hptH3L_pDca", scale,f1, f2, fMc,c,pdf,"#p DCA", "p",2,"SE-ME","MC", "dp#pi" );
  projAndScaleComp("hptH3L_dpDca", scale,f1, f2, fMc_ld,c,pdf,"dp pair DCA", "p",2,"SE-ME","MC_#Lambda", "dp#pi",0,1 );
  // projAndScaleComp("hptppichi2prim",scale, f1, f2,fMc_ld, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",2 , "SE-ME" ,"MC_#Lambda","p#pi", 0, 1);

  // projAndComp2("hptH3L_chi2topoSig", "hptH3L_chi2topo", fMc_ld, fMc_ld, c,pdf,"(p#pid) #chi^{2}_{topo}" ,"",8 , "Sig" ,"Wide","p#pi");
  // projAndComp2("hptH3L_ldlSig", "hptH3L_ldl", fMc_ld, fMc_ld, c,pdf,"ldl" ,"",8 , "Sig" ,"Wide","p#pi");
  // projAndComp2("hptH3L_lSig", "hptH3L_l", fMc_ld, fMc_ld, c,pdf,"ldl" ,"",8 , "Sig" ,"Wide","p#pi");
  // projAndComp2("hptH3L_lSBR", "hptH3L_l", fMc_ld, fMc_ld, c,pdf,"ldl" ,"",8 , "Sig" ,"Wide","p#pi");
  projAndComp2("hptH3L_ldlSig",  "hptH3L_ldl", fMc, fMc_ld, c,pdf,"ldl" ,"",4 , "Sig" ,"Ld","p#pi");
  projAndComp2("hptH3L_lSig",  "hptH3L_l", fMc, fMc_ld, c,pdf,"l" ,"",4 , "Sig" ,"Ld","p#pi");
  projAndComp2("hptH3L_chi2topoSig",  "hptH3L_chi2topo", fMc, fMc_ld, c,pdf,"chi2topo" ,"",4 , "Sig" ,"Ld","p#pi");
  projAndComp2("hptH3L_chi2ndfSig",  "hptH3L_chi2ndf", fMc, fMc_ld, c,pdf,"chi2ndf" ,"",2 , "Sig" ,"Ld","p#pi");
         
  //add 3 D case
  // cout<<"chi2ndf" << calpurity("h3H3L_chi2ndf",0, 3.5, 1.5, 1, -0.4, -0.6, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 1, 0, 10)<<endl;
  // cout<<"chi2topo" << calpurity("h3H3L_chi2topo",0, 3, 1.5, 1., -0.4, -0.6, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 1, 0, 10)<<endl;
  // cout<<"decaylength" << calpurity("h3H3L_l",8, 99, 1.5, 1, -0.4, -0.6,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 10, 0 ,100 )<<endl;

  // highy = -0.4; lowy = -0.6;
  double error;
  cout<<"chi2ndf" << calpurity("h3H3L_chi2ndf",0, 3.5, highpt, lowpt, highy , lowy, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 1, 0, 5, error, 0.25 )<<endl;
  // // cout<<"chi2ndf" << calpurity("hptH3L_chi2ndf",0, 3.5, 0., 1, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 5, error )<<endl;
  cout<<"chi2ndf" << calpurity("h3H3L_chi2ndf",0, 3.5, 1.5, 1., highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 1, 0, 5, error , 0.25)<<endl;
  cout<<"chi2ndf" << calpurity("h3H3L_chi2ndf",0, 3.5, 2, 1.5, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 1, 0, 5, error , 0.25)<<endl;
  cout<<"chi2ndf" << calpurity("h3H3L_chi2ndf",0, 3.5, 2.5, 2., highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 1, 0, 5, error , 0.25)<<endl;
  cout<<"chi2topo" << calpurity("h3H3L_chi2topo",0, 3, highpt, lowpt, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 1, 0, 10, error , 0.25)<<endl;
  cout<<"chi2topo" << calpurity("h3H3L_chi2topo",0, 3, 1.5, 1., highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 1 , 0, 10, error , 0.25)<<endl;
  cout<<"chi2topo" << calpurity("h3H3L_chi2topo",0, 3, 2, 1.5, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 1 , 0, 10, error , 0.25)<<endl;
  cout<<"chi2topo" << calpurity("h3H3L_chi2topo",0, 3, 2.5, 2., highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 1 , 0, 10, error, 0.25 )<<endl;
  // cout<<"chi2topo" << calpurity("hptH3L_chi2topo",0, 3, lowpt, highpt, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 4 , 0, 50, error )<<endl;
  cout<<"decaylength" << calpurity("h3H3L_l",8, 99,  highpt, lowpt, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 2, 0 ,100 , error , 0.25)<<endl;
  cout<<"decaylength" << calpurity("h3H3L_l",8, 99,  1.5, 1, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 5, 0 ,100 , error , 0.25)<<endl;
  cout<<"decaylength" << calpurity("h3H3L_l",8, 99,  2, 1.5, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 5, 0 ,100 , error , 0.25)<<endl;
  cout<<"decaylength" << calpurity("h3H3L_l",8, 99,  2.5, 2., highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 5, 0 ,100, error , 0.25)<<endl;
  // cout<<"decaylength" << calpurity("hptH3L_l",8, 100,  0.2, 2,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 4 , 0 ,50 )<<endl;
  // cout<<"l/dl" << calpurity("h3H3L_ldl",5, 50, highpt,  lowpt, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l/dl", 2, 3, 50 )<<endl;

  //2D
  // cout<<"chi2ndf" << calpurity("hptH3L_chi2ndf",0, 3.5, 1, 2.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 5)<<endl;
  // cout<<"chi2ndf" << calpurity("hptH3L_chi2ndf",0, 3.5, 0., 1, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 5)<<endl;
  // cout<<"chi2ndf" << calpurity("hptH3L_chi2ndf",0, 3.5, 1., 1.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 5)<<endl;
  // cout<<"chi2ndf" << calpurity("hptH3L_chi2ndf",0, 3.5, 1.5, 2, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 5)<<endl;
  // cout<<"chi2ndf" << calpurity("hptH3L_chi2ndf",0, 3.5, 2., 2.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 5)<<endl;
  // cout<<"chi2topo" << calpurity("hptH3L_chi2topo",0, 3, 1, 2.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 4 , 0, 50)<<endl;
  // cout<<"chi2topo" << calpurity("hptH3L_chi2topo",0, 3, 1, 1.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 4 , 0, 50)<<endl;
  // cout<<"chi2topo" << calpurity("hptH3L_chi2topo",0, 3, 1.5, 2., f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 4 , 0, 50)<<endl;
  // cout<<"chi2topo" << calpurity("hptH3L_chi2topo",0, 3, 2., 2.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 4 , 0, 50)<<endl;
  // // cout<<"chi2topo" << calpurity("hptH3L_chi2topo",0, 3, lowpt, highpt, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 4 , 0, 50)<<endl;
  // cout<<"decaylength" << calpurity("hptH3L_l",8, 100,  lowpt, highpt,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 4 , 0 ,50 )<<endl;
  // // cout<<"decaylength" << calpurity("hptH3L_l",8, 100,  0.2, 2,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 4 , 0 ,50 )<<endl;
  // cout<<"l/dl" << calpurity("hptH3L_ldl",5, 50,  lowpt, highpt,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l/dl", 2, 3, 50 )<<endl;
  // cout<<"dchi2prim" << calpurity("hptH3L_dchi2prim",0, 19,  lowpt, highpt,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L d #chi^{2}",0 )<<endl;
  // cout<<"pichi2prim" << calpurity("hptH3L_pichi2prim",0, 19,  lowpt, highpt,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L pi #chi^{2}",0 )<<endl;
  // cout<<"pip" << calpurity("hptH3L_ppil",0, 39,  lowpt, highpt,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L ppi" , 0)<<endl;
  pdf->On();
  pdf->Close();

}
void drawMixData_015pt()
{
  double highpt = 2.5, lowpt = 1., lowy=-0.8, highy = -0.2;
  // double highpt = 2.5, lowpt = 1, lowy=-0.9, highy = -0.1;
  // double highpt = 4, lowpt = 0, lowy=-1.5, highy = 0;
  TString histname="hH3LMassPtY";
  // TString histname="hH3LMassPtY_5_40";
  SetsPhenixStyle();
  TCanvas* c = new TCanvas("c1","c1");
  // TPDF* pdf = new TPDF("MixEventQA_check.pdf");
  // TPDF* pdf = new TPDF("MixEventQA_Jul27.pdf");
  // TPDF* pdf = new TPDF("MixEventQA_beforeDcacut.pdf");
  TPDF* pdf = new TPDF("MixEventQA_0050_015pt.pdf");
  pdf->Off();
  gStyle->SetPalette(1);

  // TFile *f1 = TFile::Open("fout_H3L_data_SE_large.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_large.root"); 

  // TFile *f1 = TFile::Open("out_KF_test/fout_H3L_KF_test.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_KF_Jul22_part.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_KF_Jul25.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_SE_Jul25.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_SE_Jul27.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Aug09.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_SE_Aug16.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_SE_Aug20.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_SE_Aug23.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02_nodca.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02.root"); 
  TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02_0050_015pt.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02_nodca_nopt.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_SE_Aug22.root"); 
  // TFile *f1 = TFile::Open("fout_KF_test.root"); 
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
  double nEvents_se = hcent_se->Integral( 4, 9);
  // double nEvents_se = hcent_se->Integral( 1, 9);

  // TFile *f2 = TFile::Open("fout_H3L_data_ME_large.root"); 
  // TFile *f2 = TFile::Open("fout_H3L_ME_Jul22_part.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Aug16.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Aug20.root"); 
  // TFile *f2 = TFile::Open("output2/fout_H3L_data_ME_Aug24.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Sep02.root"); 
  TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Sep02_0050_015pt.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Sep02_nodca.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Sep02_nodca_nopt.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Aug23.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Aug22.root"); 
  // TFile *f2 = TFile::Open("fout_ME_test.root"); 
  // TFile *f2 = TFile::Open("fout_H3L_ME_Jul25.root"); 
  // TFile *f2 = TFile::Open("fout_H3L_data_SE.root"); 
  // TH2F* h2bk = (TH2F*)f2->Get("hptH3Lmass")->Clone("hptH3Lmass_ME");
  // TH3F* h2bk = (TH3F*)f2->Get("hH3LMassPtY_5_40")->Clone("hptH3Lmass_bk");
  TH3F* h2bk = (TH3F*)f2->Get(histname.Data())->Clone("hptH3Lmass_bk");
  h2bk->SetDirectory(0);
  TH1F* hbk = (TH1F*)h2bk->ProjectionY("hbk", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2bk->GetZaxis()->FindBin(lowy), h2bk->GetZaxis()->FindBin(highy));
  hbk->SetDirectory(0);

  // TFile* f3 = TFile::Open("fout_H3L_data_RT_Jul22_full.root");
  // TFile* f3 = TFile::Open("fout_H3L_data_RT_Jul25.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_RT_Jul27.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_KF_Aug20.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_SE_Aug23.root");
  // TFile* f3 = TFile::Open("output2/fout_H3L_data_KF_Aug24.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_SE_Sep02_nodca.root");
  TFile* f3 = TFile::Open("rootfile/fout_H3L_data_SE_Sep02_0050_015pt.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_RT_Aug29.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_KF_Aug22.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_KF_Aug16.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_RT_Aug13.root");
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
  double sig_sb =  hsig->Integral(hsig->GetXaxis()->FindBin(3.01),  hsig->GetXaxis()->FindBin(3.02));
  // double rt_sb =  hrt->Integral(hrt->GetXaxis()->FindBin(2.97), hrt->GetXaxis()->FindBin(2.98)) + hrt->Integral(hrt->GetXaxis()->FindBin(3.0),  hrt->GetXaxis()->FindBin(3.02));
  double rt_sb = hrt->Integral(hrt->GetXaxis()->FindBin(3.01),  hrt->GetXaxis()->FindBin(3.02));
  // double rt_sb =  hrt->Integral(hrt->GetXaxis()->FindBin(2.96), hrt->GetXaxis()->FindBin(2.98)) ;
  double scale_rt = sig_sb/rt_sb;
  // sig_sb = hsig->Integral(hsig->GetXaxis()->FindBin(2.97), hsig->GetXaxis()->FindBin(2.98)) + hsig->Integral(hsig->GetXaxis()->FindBin(3.0),  hsig->GetXaxis()->FindBin(3.02));
  double bk_sb =   hbk->Integral(hbk->GetXaxis()->FindBin(3.01),  hbk->GetXaxis()->FindBin(3.02));
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
  // hrt->Draw("same");
  TLegend* leg_sig = new TLegend(0.65,0.25,0.88,0.45);
  leg_sig->AddEntry(hbk, "ME", "pl");
  // leg_sig->AddEntry(hrt, "RT", "pl");
  // leg_sig->AddEntry(hrt, "KF", "pl");
  leg_sig->AddEntry(hsig, "SE", "pl");
  leg_sig->Draw();

  drawLatex( 0.65,0.61,Form("%0.1f<y<%0.1f",lowy, highy ), 0.055);
  drawLatex( 0.65,0.54,Form("%0.1f<p_{T}<%0.1f GeV/c^{2}", lowpt, highpt), 0.055);
  // drawLatex( 0.65,0.47,Form("5-40%s", "%"), 0.055);
  drawLatex( 0.65,0.47,Form("0-50%s", "%"), 0.055);
  // drawBox( 2.97, hsig->GetMinimum(),2.98, hsig->GetMaximum()*0.5 , kBlue-9, 1001 , 1,0.3); 
  // drawBox( 2.989, hsig_bk->GetMinimum(),2.995, hsig_bk->GetMaximum()*0.5 , kOrange-4, 1001,1,0.3); 
  // drawBox( 3.01, hsig->GetMinimum(),3.02, hsig->GetMaximum()*0.5 , kBlue-9, 1001,1,0.3 ); 

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
  // hsig_rt->Draw("same");
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
  // hsig_rt->Fit(resfit_rt,"R");
  fit_rt->SetParameter(3, resfit_rt->GetParameter(0) );
  fit_rt->SetParameter(4, resfit_rt->GetParameter(1) );
  // hsig_rt->Fit(fit_rt,"R");
  resfit_rt->SetParameter(0, fit_rt->GetParameter(3));
  resfit_rt->SetParameter(1, fit_rt->GetParameter(4));
  // resfit_rt->Draw("same");

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
  double significance = yield_me/sqrt(sp_counts);
  // double significance = yield_rt/sqrt(sp_counts);
  double s_me = yield_me/(fit->GetParError(0)/hsig_bk->GetBinWidth(1));
  double s_rt = yield_rt/(fit_rt->GetParError(0)/hsig_rt->GetBinWidth(1));
  cout<<"significance: " <<significance << endl;

  TLegend* leg = new TLegend( 0.72, 0.68 ,0.9,0.9 );
  // leg->AddEntry(hbk, "mix-event(ME)","pl");
  // leg->AddEntry(hrt, "rotate d(RT) (scale)","pl");
  leg->AddEntry(hsig, "SE","pl");
  leg->AddEntry(hsig_bk, "SE-ME","pl");
  // leg->AddEntry(hsig_rt, "SE-KF","pl");
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
  // drawBox( 2.97, hsig_bk->GetMinimum(),2.98, hsig_bk->GetMaximum()*0.5 , kBlue-9, 1001 , 1,0.3); 
  // drawBox( 2.989, hsig_bk->GetMinimum(),2.995, hsig_bk->GetMaximum()*0.5 , kOrange-4, 1001,1,0.3); 
  drawBox( 3.01, hsig_bk->GetMinimum(),3.02, hsig_bk->GetMaximum()*0.5 , kBlue-9, 1001,1,0.3 ); 

  addpdf(pdf);

  hsig_bk->Draw();
  // hsig_rt->Draw("same");
  TLegend* leg2 = new TLegend( 0.2 , 0.7 ,0.4,0.9  );
  leg2->AddEntry(hsig_rt, "SE-RT","pl");
  leg2->AddEntry(hsig_bk, "SE-ME","pl");
  leg2->Draw();

  addpdf(pdf);
  //
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

  TH2F* hYield = new TH2F( "hYield", "hYield;y,pt", 5, -1., 0., 6, 0, 3);
  double xbw = hYield->GetXaxis()->GetBinWidth(1);
  double ybw = hYield->GetYaxis()->GetBinWidth(1);
  for (int ix=1; ix<=hYield->GetNbinsX(); ix++) {
    for (int iy=1; iy<=hYield->GetNbinsY();iy++){
      double xlow = hYield->GetXaxis()->GetBinLowEdge(ix); // y
      double ylow = hYield->GetYaxis()->GetBinLowEdge(iy); // pt
      double err;
      hYield->SetBinContent( ix , iy, fityield( ylow, ylow+ybw, xlow, xlow+xbw, err, f1, f2, c, pdf) );
      hYield->SetBinError( ix , iy, err);
    }
  }
  hYield->Draw("colz text");
  hYield->SaveAs("fH3L_yield_0050.root");
  addpdf(pdf);

  // TFile* fMc_ld = TFile::Open("fout_Lambda_MC_Cuts_0050_015pt_sys_pid3.root");
  TFile* fMc_ld = TFile::Open("fout_Lambda_MC_Cuts_0050_015pt_sys_mixpid3.root");
  TFile* fMc = TFile::Open("fout_H3L_MC_0050_015pt_sys.root");
  // TFile* fMc = TFile::Open("fout_H3L_MC_0080_010pt.root");
  TH1F* hPhase[3];
  TH1F* hPhaseCor[3];
  TH1F* hBr[3];
  TGraphErrors* g[3];
  double edge[3][4]={
    { 1.7, 2.1, 2.7, 3.1},
    { 1.2, 1.6, 2.0, 2.4},
    { 1.2, 1.8, 2.2, 2.6}
  };

  double ybin[4]={ 0, -0.25, -0.50, -0.75};
  hPhase[0] = new TH1F("hPhase0", "hPhase0;pt", 25, 1, 3.5);
  hPhase[0]=(TH1F*)hPhase[0]->Rebin(3, "hPhase0", edge[0]);
  hPhase[1] = new TH1F("hPhase1", "hPhase1;pt", 25, 1, 3.5);
  hPhase[1]=(TH1F*)hPhase[1]->Rebin(3, "hPhase1", edge[1]);
  hPhase[2] = new TH1F("hPhase2", "hPhase2;pt", 25, 1, 3.5);
  hPhase[2]=(TH1F*)hPhase[2]->Rebin(3, "hPhase2", edge[2]);
  TFile* fMcH3L = new TFile("fMC_H3L_0050.root");
  // TFile* fRcH3L = new TFile("fout_H3L_MC_0050_015pt.root");
  TFile* fRcH3L = new TFile("fout_H3L_MC_0050_015pt_sys.root");
  // TFile* fMcH3L = new TFile("fMC_H3L_0080.root");
  // TFile* fRcH3L = new TFile("fout_H3L_MC_0080_010pt.root");
  TH3F* h3Mc = (TH3F*)fMcH3L->Get("hH3LMassPtY")->Clone("h3Mc");
  h3Mc->SetDirectory(0);
  TH3F* h3Rc = (TH3F*)fRcH3L->Get("hH3LMassPtY")->Clone("h3Rc");
  h3Rc->SetDirectory(0);

  h3Mc->Sumw2();
  h3Rc->Sumw2();
  TH2F* h2MC = (TH2F*)h3Mc->Project3D("xz");
  TH2F* h2Rc = (TH2F*)h3Rc->Project3D("xz");

  TH2F* h2Eff = (TH2F*)h2Rc->Clone("h2Eff");
  TH2F* h2temp = (TH2F*)h2MC->Clone("h2temp");
  h2temp->RebinY(10);
  h2temp->RebinX(10);
  h2Eff->RebinY(10);
  h2Eff->RebinX(10);
  h2Eff->Divide(h2temp);
  h2Eff->Draw("colz text");
  h2Eff->GetYaxis()->SetRangeUser(0,4.5);
  addpdf(pdf);

  TH1F* heff[3];
  for (int ij=0;ij<3;ij++){ 
    TH1F* h1Mc = (TH1F*)h2MC->ProjectionY(Form("h1Mc%d", ij), h2MC->GetXaxis()->FindBin(ybin[ij+1] +1e-6), h2MC->GetXaxis()->FindBin(ybin[ij]-1e-6) );
    h1Mc = (TH1F*)h1Mc->Rebin( 3, Form("hMc%d",ij), edge[ij]);
    TH1F* h1Rc = (TH1F*)h2Rc->ProjectionY(Form("h1Rc%d", ij), h2Rc->GetXaxis()->FindBin(ybin[ij+1]+1e-6), h2Rc->GetXaxis()->FindBin(ybin[ij]-1e-6) );
    heff[ij] = (TH1F*)h1Rc->Rebin( 3, Form("heff%d",ij), edge[ij]);
    heff[ij]->Divide(h1Mc);
    heff[ij]->GetYaxis()->SetTitle("Eff.");
    heff[ij]->Draw();
    addpdf(pdf);
  }

  double dy=0.25;
  double dpt,pt;
  
  TFile* f2b = TFile::Open("h_h3l_corr_yield.root_cut00002_cent0_yuehang.root");

  double ny = 3;
  for (int ij=0; ij<ny;ij++){
    double npt =3;
    if (ij==0) npt =2; 
    for (int ipt=0; ipt<npt;ipt++){
      double err;
      dpt = edge[ij][ipt+1]-edge[ij][ipt];
      pt = 0.5*(edge[ij][ipt+1]+edge[ij][ipt]);
      double yield = fityield( edge[ij][ipt], edge[ij][ipt+1], ybin[ij+1], ybin[ij], err, f1, f2, c, pdf);
      // cout << edge[ij][ipt]<<" "<<edge[ij][ipt+1]<<" "<<ybin[ij+1]<<" "<<ybin[ij]<<" "<< yield<<" "<<err <<endl;
      if (ij==0 && ipt==2) { yield=0; err=0;}
      hPhase[ij]->SetBinContent( ipt+1, yield/dy/dpt/pt/2./3.1415926/nEvents_se);
      hPhase[ij]->SetBinError( ipt+1, err/dy/dpt/pt/2./3.1415926/nEvents_se );
    }
    hPhase[ij]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPhase[ij]->GetYaxis()->SetTitle("Raw d^{2}N/(N_{evt}2#pip_{T}dydp_{T})");
    hPhase[ij]->Draw();
    hPhase[ij]->SetDirectory(0);
    addpdf(pdf); 

    hPhaseCor[ij]=(TH1F*)hPhase[ij]->Clone(Form("hYieldCor_%d",ij));
    hPhaseCor[ij]->SetDirectory(0);

    double perr;
    double purity1_s = calpurity("h3H3L_chi2ndf",0, 3.5, edge[ij][2], edge[ij][0], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 10, perr, 0.25);

    // calpurityCut("h3H3L_chi2topo", "hH3LMassPtY", 0, 3.,  edge[ij][2],  edge[ij][0], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);

    for (int ipt=0;ipt<npt;ipt++)
    {
      double y3b = hPhaseCor[ij]->GetBinContent(ipt+1);
      double y3berr = hPhaseCor[ij]->GetBinError(ipt+1)/y3b;
      double eff = heff[ij]->GetBinContent(ipt+1);
      double efferr = heff[ij]->GetBinError(ipt+1)/eff;
      // double purity1 = calpurity("h3H3L_chi2ndf",0, 3.5, edge[ij][ipt+1], edge[ij][ipt], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 1, 0, 10, perr, 0.25);
      // double purity2 = calpurity("h3H3L_chi2topo",0, 3., edge[ij][ipt+1], edge[ij][ipt], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 1, 0, 10, perr, 0.25);
      // double purity = (purity1+purity2)*0.5;
      // perr = fabs(purity2-purity1);
      // perr = perr/purity;
    double purity2_s = calpurity("h3H3L_chi2topo",0, 3.,  edge[ij][ipt],  edge[ij][ipt], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
      double purity = purity2_s;
      // double purity;
      // purity=1;
      // perr = 0;
      double yield_cor = y3b/eff*purity;
      double err = sqrt(y3berr*y3berr + efferr*efferr + perr*perr )*yield_cor;
      hPhaseCor[ij]->SetBinContent(ipt+1,yield_cor );
      hPhaseCor[ij]->SetBinError(ipt+1, err);
      cout <<"pt "<<hPhaseCor[ij]->GetBinCenter(ipt+1) << "purity "<<purity <<" eff " <<eff <<" yield " << yield_cor << " err "<< err<< " s="<<yield_cor/err<<endl;
    }
    hPhaseCor[ij]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPhaseCor[ij]->GetYaxis()->SetTitle("d^{2}N/(N_{evt}2#pip_{T}dydp_{T})");
    hPhaseCor[ij]->Draw();
    hPhaseCor[ij]->SetDirectory(0);
    addpdf(pdf);

    hBr[ij]=(TH1F*)hPhaseCor[ij]->Clone(Form("hBr_%d",ij));
    hBr[ij]->SetDirectory(0);
    g[ij] = (TGraphErrors*)f2b->Get(Form("t_sgct1_corr_yield[0][%d]", ij));
    g[ij]->Draw();
    for (int ipt=0;ipt<npt;ipt++)
    {
      double y3b = hBr[ij]->GetBinContent(ipt+1);
      double y3berr = hBr[ij]->GetBinError(ipt+1);
      double y2b = g[ij]->GetPointY(ipt);
      double y2bpt = g[ij]->GetPointX(ipt);
      double y2berr = g[ij]->GetErrorY(ipt);
      if (ij==1) {y2berr=y2berr*10; y2b=y2b*10.;}
      if (ij==2) {y2berr=y2berr*100; y2b=y2b*100.;}
      // cout <<"y2b significance:"<<y2b/y2berr<<" yield "<<y2b << endl;
      hBr[ij]->SetBinContent(ipt+1, y2b/(y2b+y3b) );
      hBr[ij]->SetBinError(ipt+1, sqrt((y2b*y2b)/pow((y2b+y3b), 4)*(y3berr*y3berr) + (y3b*y3b)/pow((y2b+y3b), 4)*(y2berr*y2berr)));
      if (y2b/y2berr<2) {
         hBr[ij]->SetBinContent(ipt+1, 0);
         hBr[ij]->SetBinError(ipt+1, 0);
      } 
      cout <<"bincenter: " << hBr[ij]->GetBinCenter(ipt+1) << " content:"<< hBr[ij]->GetBinContent(ipt+1)<<" err:"<<hBr[ij]->GetBinError(ipt+1)<< "y2b significance:"<<y2b/y2berr<<" yield "<<y2b << " y2bpt "<< y2bpt<< " y3b yield "<< y3b<<endl;
    } 
    hBr[ij]->Draw();
    hBr[ij]->GetYaxis()->SetTitle("2body/(2body + 3body)");
    drawLatex( 0.65,0.88 , Form("%0.2f<y<%0.2f", ybin[ij+1],ybin[ij]) , 0.055);
    addpdf(pdf);
  }

  TLegend* legb = new TLegend(0.2, 0.2, 0.4, 0.4);
  TH1F* htemp = new TH1F("htemp","htemp;p_{T} (GeV/c); 2-body/(2-body+3-body)",  5, 1, 3);
  htemp->GetYaxis()->SetRangeUser( 0, 0.7);
  htemp->Draw();
  int color[]={kRed-4, kGreen+2, kBlue-4 };
  for (int i=0;i<3;i++)
  {
    hBr[i]->Draw("same");
    hBr[i]->SetLineColor(color[i]);
    hBr[i]->SetMarkerColor(color[i]);
    legb->AddEntry( hBr[i], Form("%0.2f<y<%0.2f", ybin[i+1], ybin[i]), "pl");
  }
  legb->Draw();
  addpdf(pdf);

  TFile * fout = new TFile("fout_0050.root","recreate");
  fout->cd();
  for (int i=0;i<3;i++){
    hPhase[i]->Write();
    hPhaseCor[i]->Write();
    hBr[i]->Write();
    heff[i]->Write();
  }
  fout->Close();

  pdf->On();
  pdf->Close();
  return;

  c = new TCanvas( "c","c", 800,800);
  projAndComp("hptH3L_lSBL", f1, f2, c,pdf,"l", "p",4,"SE_SBL","ME_SBL", "dp#pi", scale);
  projAndComp("hptH3L_lSBR", f1, f2, c,pdf,"l", "p",4,"SE_SBR","ME_SBR", "dp#pi", scale);
  projAndComp("hptH3L_ldlSBL", f1, f2, c,pdf,"l/dl", "p",4,"SE_SBL","ME_SBL", "dp#pi", scale);
  projAndComp("hptH3L_ldlSBR", f1, f2, c,pdf,"l/dl", "p",4,"SE_SBR","ME_SBR", "dp#pi", scale);
  projAndComp("hptH3L_chi2ndfSBL", f1, f2, c,pdf,"chi2NDF", "p",4,"SE_SBL","ME_SBL", "dp#pi", scale);
  projAndComp("hptH3L_chi2ndfSBR", f1, f2, c,pdf,"chi2NDF", "p",4,"SE_SBR","ME_SBR", "dp#pi", scale);
  projAndComp("hptH3L_chi2topoSBL", f1, f2, c,pdf,"chi2topo", "p",4,"SE_SBL","ME_SBL", "dp#pi" , scale);
  projAndComp("hptH3L_chi2topoSBR", f1, f2, c,pdf,"chi2topo", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  projAndComp("hptH3L_dchi2primSBL", f1, f2, c,pdf,"d chi2primary", "p",4,"SE_SBL","ME_SBL", "dp#pi", scale);
  projAndComp("hptH3L_dchi2primSBR", f1, f2, c,pdf,"d chi2primary", "p",4,"SE_SBR","ME_SBR", "dp#pi", scale);
  projAndComp("hptH3L_pchi2primSBL", f1, f2, c,pdf,"p chi2primary", "p",4,"SE_SBL","ME_SBL", "dp#pi", scale);
  projAndComp("hptH3L_pchi2primSBR", f1, f2, c,pdf,"p chi2primary", "p",4,"SE_SBR","ME_SBR", "dp#pi", scale);
  projAndComp("hptH3L_pichi2primSBL", f1, f2, c,pdf,"#pi chi2primary", "p",4,"SE_SBL","ME_SBL", "dp#pi", scale);
  projAndComp("hptH3L_pichi2primSBR", f1, f2, c,pdf,"#pi chi2primary", "p",4,"SE_SBR","ME_SBR", "dp#pi", scale);
  projAndComp("hptH3L_dDcaSBL", f1, f2, c,pdf,"d DCA", "p",1,"SE_SBL","ME_SBL","dp#pi"  , scale);
  projAndComp("hptH3L_dDcaSBR", f1, f2, c,pdf,"d DCA", "p",1,"SE_SBR","ME_SBR","dp#pi"  , scale);
  projAndComp("hptH3L_piDcaSBL", f1, f2, c,pdf,"#pi DCA", "p",4,"SE_SBL","ME_SBL", "dp#pi" , scale);
  projAndComp("hptH3L_piDcaSBR", f1, f2, c,pdf,"#pi DCA", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  projAndComp("hptH3L_pDcaSBL", f1, f2, c,pdf,"#p DCA", "p",4,"SE_SBL","ME_SBL", "dp#pi" , scale);
  projAndComp("hptH3L_pDcaSBR", f1, f2, c,pdf,"#p DCA", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  projAndComp("hptH3L_dpDcaSBL", f1, f2, c,pdf,"dp pair DCA", "p",4,"SE_SBL","ME_SBL", "dp#pi" , scale);
  projAndComp("hptH3L_dpDcaSBR", f1, f2, c,pdf,"dp pair DCA", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  // projAndComp("hptppimassSBL", f1, f2, c,pdf,"p#pi Mass", "p",4,"SE_SBL","ME_SBL", "dp#pi" );
  // projAndComp("hptH3L_ppimassSBR", f1, f2, c,pdf,"p#pi Mass", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  // projAndComp("hptH3L_ppichi2ndfSBR", f1, f2, c,pdf,"p#pi Mass", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  // projAndComp("hptH3L_ppichi2ndfSBL", f1, f2, c,pdf,"p#pi chi2ndf", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  // projAndComp("hptH3L_ppichi2primSBL", f1, f2, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "SE" ,"BK","p#pi",scale);

  // projAndComp("hptH3L_lSBR", f1, f3, c,pdf,"l", "p",4,"SE_SBR","RT_SBR", "dp#pi");
  // projAndComp("hptH3L_ldlSBR", f1, f3, c,pdf,"l/dl", "p",4,"SE_SBR","RT_SBR", "dp#pi");
  // projAndComp("hptH3L_chi2ndfSBR", f1, f3, c,pdf,"chi2NDF", "p",4,"SE_SBR","RT_SBR", "dp#pi");
  // projAndComp("hptH3L_chi2topoSBR", f1, f3, c,pdf,"chi2topo", "p",4,"SE_SBR","RT_SBR", "dp#pi" );
  // projAndComp("hptH3L_dchi2primSBR", f1, f3, c,pdf,"d chi2primary", "p",4,"SE_SBR","RT_SBR", "dp#pi");
  // projAndComp("hptH3L_pchi2primSBR", f1, f3, c,pdf,"p chi2primary", "p",4,"SE_SBR","RT_SBR", "dp#pi");
  // projAndComp("hptH3L_pichi2primSBR", f1, f3, c,pdf,"#pi chi2primary", "p",4,"SE_SBR","RT_SBR", "dp#pi");
  // projAndComp("hptH3L_dDcaSBR", f1, f3, c,pdf,"d DCA", "p",1,"SE_SBR","RT_SBR","dp#pi"  );
  // projAndComp("hptH3L_piDcaSBR", f1, f3, c,pdf,"#pi DCA", "p",4,"SE_SBR","RT_SBR", "dp#pi" );
  // projAndComp("hptH3L_pDcaSBR", f1, f3, c,pdf,"#p DCA", "p",4,"SE_SBR","RT_SBR", "dp#pi" );
  // projAndComp("hptH3L_dpDcaSBR", f1, f3, c,pdf,"dp pair DCA", "p",4,"SE_SBR","RT_SBR", "dp#pi" );
  // projAndComp("hptH3L_ppimassSBR", f1, f3, c,pdf,"p#pi Mass (GeV/c^{2})", "plhist",1 ,"SE" ,"BK","p#pi");
  // // projAndComp("hptppichi2ndfSBR", f1, f2, c,pdf,"p#pi #chi^{2}_{NDF}","" ,1 , "SE" ,"BK","p#pi");
  // // projAndComp("hptppichi2primSBR", f1, f2, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "SE" ,"BK","p#pi");
  // // projAndComp("hptppilSBR", f1, f2, c,pdf,"(p#pi) l" ,"",1 , "SE" ,"BK","p#pi");
  // // projAndComp("hptppildlSBR", f1, f2, c,pdf,"(p#pi) l/#deltal" ,"",1 , "SE" ,"BK","p#pi");

  // projAndComp("hptH3L_lSBL", f1, f3, c,pdf,"l", "p",4,"SE_SBL","RT_SBL", "dp#pi");
  // projAndComp("hptH3L_ldlSBL", f1, f3, c,pdf,"l/dl", "p",4,"SE_SBL","RT_SBL", "dp#pi");
  // projAndComp("hptH3L_chi2ndfSBL", f1, f3, c,pdf,"chi2NDF", "p",4,"SE_SBL","RT_SBL", "dp#pi");
  // projAndComp("hptH3L_chi2topoSBL", f1, f3, c,pdf,"chi2topo", "p",4,"SE_SBL","RT_SBL", "dp#pi" );
  // projAndComp("hptH3L_dchi2primSBL", f1, f3, c,pdf,"d chi2primary", "p",4,"SE_SBL","RT_SBL", "dp#pi");
  // projAndComp("hptH3L_pchi2primSBL", f1, f3, c,pdf,"p chi2primary", "p",4,"SE_SBL","RT_SBL", "dp#pi");
  // projAndComp("hptH3L_pichi2primSBL", f1, f3, c,pdf,"#pi chi2primary", "p",4,"SE_SBL","RT_SBL", "dp#pi");
  // projAndComp("hptH3L_dDcaSBL", f1, f3, c,pdf,"d DCA", "p",1,"SE_SBL","RT_SBL","dp#pi"  );
  // projAndComp("hptH3L_piDcaSBL", f1, f3, c,pdf,"#pi DCA", "p",4,"SE_SBL","RT_SBL", "dp#pi" );
  // projAndComp("hptH3L_pDcaSBL", f1, f3, c,pdf,"#p DCA", "p",4,"SE_SBL","RT_SBL", "dp#pi" );
  // projAndComp("hptH3L_dpDcaSBL", f1, f3, c,pdf,"dp pair DCA", "p",4,"SE_SBL","RT_SBL", "dp#pi" );
  // // projAndComp("hptppimassSBL", f1, f3, c,pdf,"p#pi Mass (GeV/c^{2})", "plhist",1 ,"SE" ,"BK","p#pi");
  // projAndComp("hptH3L_ppimassSBL", f1, f3, c,pdf,"p#pi Mass (GeV/c^{2})", "plhist",1 ,"SE" ,"BK","dp#pi");
  // // projAndComp("hptppichi2ndfSBL", f1, f2, c,pdf,"p#pi #chi^{2}_{NDF}","" ,1 , "SE" ,"BK","p#pi");
  // // projAndComp("hptppichi2primSBL", f1, f2, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "SE" ,"BK","p#pi");
  // // projAndComp("hptppilSBL", f1, f2, c,pdf,"(p#pi) l" ,"",1 , "SE" ,"BK","p#pi");
  // // projAndComp("hptppildlSBL", f1, f2, c,pdf,"(p#pi) l/#deltal" ,"",1 , "SE" ,"BK","p#pi");

  // projAndScaleComp("hptH3L_l", scale,f1, f2, f2, c,pdf,"l", "p",1,"SE-ME","ME", "dp#pi");
  // projAndScaleComp("hptH3L_ppimass", scale, f1, f2, f2, c,pdf,"M(p#pi)", "p",1,"SE-ME","ME", "p#pi");
  // // projAndScaleCompNoRatio("hptppimass", scale, f1, f2, f2, c,pdf,"M(p#pi)", "p",1,"SE-ME","ME", "p#pi");
  // projAndScaleComp("hptH3L_ldl", scale,f1, f2, f2, c,pdf,"l/dl", "p",1,"SE-ME","ME", "dp#pi");
  // projAndScaleComp("hptH3L_chi2ndf", scale,f1, f2,f2, c,pdf,"chi2NDF", "p",1,"SE-ME","ME", "dp#pi");
  projAndScaleComp("hptH3L_chi2topo", scale,f1, f2, f2,c,pdf,"chi2topo", "p",1,"SE-ME","ME", "dp#pi" );
  // projAndScaleComp("hptH3L_dchi2prim", scale,f1, f2, f2,c,pdf,"d chi2primary", "p",1,"SE-ME","ME", "dp#pi");
  // projAndScaleComp("hptH3L_pchi2prim", scale,f1, f2, f2,c,pdf,"p chi2primary", "p",1,"SE-ME","ME", "dp#pi");
  // projAndScaleComp("hptH3L_pichi2prim", scale,f1, f2, f2,c,pdf,"#pi chi2primary", "p",1,"SE-ME","ME", "dp#pi");
  // projAndScaleComp("hptH3L_dDca", scale,f1, f2, f2,c,pdf,"d DCA", "p",1,"SE-ME","ME","dp#pi"  );
  // projAndScaleComp("hptH3L_piDca", scale,f1, f2, f2,c,pdf,"#pi DCA", "p",1,"SE-ME","ME", "dp#pi" );
  // projAndScaleComp("hptH3L_pDca", scale,f1, f2, f2,c,pdf,"#p DCA", "p",1,"SE-ME","ME", "dp#pi" );
  // projAndScaleComp("hptH3L_dpDca", scale,f1, f2, f2,c,pdf,"dp pair DCA", "p",1,"SE-ME","ME", "dp#pi" );
  // // projAndScaleComp("hptppichi2prim",scale, f1, f2,f2, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "SE-ME" ,"ME","p#pi");

  // projAndScaleComp("hptH3L_l", scale_rt,f1, f3, f3, c,pdf,"l", "p",1,"SE-RT","RT", "dp#pi");
  // projAndScaleComp("hptppimass", scale_rt, f1, f3, f3, c,pdf,"M(p#pi)", "p",1,"SE-RT","RT", "p#pi");
  // projAndScaleComp("hptH3L_ldl", scale_rt,f1, f3, f3, c,pdf,"l/dl", "p",1,"SE-RT","RT", "dp#pi");
  // projAndScaleComp("hptH3L_chi2ndf", scale_rt,f1, f3,f3, c,pdf,"chi2NDF", "p",1,"SE-RT","RT", "dp#pi");
  // projAndScaleComp("hptH3L_chi2topo", scale_rt,f1, f3, f3,c,pdf,"chi2topo", "p",1,"SE-RT","RT", "dp#pi" );
  // projAndScaleComp("hptH3L_dchi2prim", scale_rt,f1, f3, f3,c,pdf,"d chi2primary", "p",1,"SE-RT","RT", "dp#pi");
  // projAndScaleComp("hptH3L_pchi2prim", scale_rt,f1, f3, f3,c,pdf,"p chi2primary", "p",1,"SE-RT","RT", "dp#pi");
  // projAndScaleComp("hptH3L_pichi2prim", scale_rt,f1, f3, f3,c,pdf,"#pi chi2primary", "p",1,"SE-RT","RT", "dp#pi");
  // projAndScaleComp("hptH3L_dDca", scale_rt,f1, f3, f3,c,pdf,"d DCA", "p",1,"SE-RT","RT","dp#pi"  );
  // projAndScaleComp("hptH3L_piDca", scale_rt,f1, f3, f3,c,pdf,"#pi DCA", "p",1,"SE-RT","RT", "dp#pi" );
  // projAndScaleComp("hptH3L_pDca", scale_rt,f1, f3, f3,c,pdf,"#p DCA", "p",1,"SE-RT","RT", "dp#pi" );
  // projAndScaleComp("hptH3L_dpDca", scale_rt,f1, f3, f3,c,pdf,"dp pair DCA", "p",1,"SE-RT","RT", "dp#pi" );
  // projAndScaleComp("hptppichi2prim",scale_rt, f1, f3,f3, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "SE-RT" ,"RT","p#pi");

  cout <<"compare with MC" <<endl;
  // projAndScaleComp("hptH3L_l", scale,f1, f2, fMc, c,pdf,"l", "p",4,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_ppimass", scale, f1, f2, fMc, c,pdf,"M(p#pi)", "p",4,"SE-ME","MC", "p#pi");
  // // projAndScaleCompNoRatio("hptH3L_ppimass", scale, f1, f2, fMc, c,pdf,"M(p#pi)", "p",4,"SE-ME","MC", "p#pi");
  // projAndScaleComp("hptH3L_ldl", scale,f1, f2, fMc, c,pdf,"l/dl", "p",4,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_chi2ndf", scale,f1, f2,fMc, c,pdf,"chi2NDF","p",2,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_chi2topo", scale,f1, f2, fMc,c,pdf,"chi2topo","p",2,"SE-ME","MC", "dp#pi" );
  // // projAndScaleComp("hptH3L_chi2", scale,f1, f2, fMc,c,pdf,"chi2topo","p",2,"SE-ME","MC", "dp#pi" );
  // projAndScaleComp("hptH3L_dchi2prim", scale,f1, f2, fMc,c,pdf,"d chi2primary", "p",5,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_pchi2prim", scale,f1, f2, fMc,c,pdf,"p chi2primary", "p",2,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_pichi2prim", scale,f1, f2, fMc,c,pdf,"#pi chi2primary", "p",2,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_dDca", scale,f1, f2, fMc,c,pdf,"d DCA", "p",2,"SE-ME","MC","dp#pi"  );
  // projAndScaleComp("hptH3L_piDca", scale,f1, f2, fMc,c,pdf,"#pi DCA", "p",2,"SE-ME","MC", "dp#pi" );
  // projAndScaleComp("hptH3L_pDca", scale,f1, f2, fMc,c,pdf,"#p DCA", "p",2,"SE-ME","MC", "dp#pi" );
  // projAndScaleComp("hptH3L_dpDca", scale,f1, f2, fMc,c,pdf,"dp pair DCA", "p",2,"SE-ME","MC", "dp#pi" );
  // // projAndScaleComp("hptppichi2prim",scale, f1, f2,fMc, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",2 , "SE-ME" ,"MC","p#pi");
  //

  // TFile* fMc_ph = TFile::Open("fout_H3L_phase_MC.root");
  // //  projAndScaleComp("hptH3L_l", scale,f1, f2, fMc_ph, c,pdf,"l", "p",4,"SE-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptppimass", scale, f1, f2, fMc_ph, c,pdf,"M(p#pi)", "p",4,"SE-ME","phase_MC", "p#pi");
  // projAndScaleComp("hptH3L_ldl", scale,f1, f2, fMc_ph, c,pdf,"l/dl", "p",4,"SE-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptH3L_chi2ndf", scale,f1, f2,fMc_ph, c,pdf,"chi2NDF","p",2,"SE-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptH3L_chi2topo", scale,f1, f2, fMc_ph,c,pdf,"chi2topo","p",2,"SE-ME","phase_MC", "dp#pi" );
  // // projAndScaleComp("hptH3L_chi2", scale,f1, f2, fMc_ph,c,pdf,"chi2topo","p",2,"SE-ME","phase_MC", "dp#pi" );
  // projAndScaleComp("hptH3L_dchi2prim", scale,f1, f2, fMc_ph,c,pdf,"d chi2primary", "p",2,"SE-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptH3L_pchi2prim", scale,f1, f2, fMc_ph,c,pdf,"p chi2primary", "p",2,"SE-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptH3L_pichi2prim", scale,f1, f2, fMc_ph,c,pdf,"#pi chi2primary", "p",2,"SE-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptH3L_dDca", scale,f1, f2, fMc_ph,c,pdf,"d DCA", "p",2,"SE-ME","phase_MC","dp#pi"  );
  // projAndScaleComp("hptH3L_piDca", scale,f1, f2, fMc_ph,c,pdf,"#pi DCA", "p",2,"SE-ME","phase_MC", "dp#pi" );
  // projAndScaleComp("hptH3L_pDca", scale,f1, f2, fMc_ph,c,pdf,"#p DCA", "p",2,"SE-ME","phase_MC", "dp#pi" );
  // projAndScaleComp("hptH3L_dpDca", scale,f1, f2, fMc_ph,c,pdf,"dp pair DCA", "p",2,"SE-ME","phase_MC", "dp#pi" );
  // projAndScaleComp("hptppichi2prim",scale, f1, f2,fMc_ph, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",2 , "SE-ME" ,"phase_MC","p#pi");

  cout <<"compare with Lambda MC" <<endl;
  projAndScaleComp("hptH3L_l", scale,f1, f2, fMc_ld, c,pdf,"l", "p",4,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_l", scale,f1, f2, fMc, c,pdf,"l", "p",4,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_ppimass", scale, f1, f2, fMc_ld, c,pdf,"M(p#pi)", "p",4,"SE-ME","MC_#Lambda", "p#pi",0,1);
  // projAndScaleCompNoRatio("hptppimass", scale, f1, f2, fMc_ld, c,pdf,"M(p#pi)", "p",4,"SE-ME","MC_#Lambda", "p#pi",0);
  // projAndScaleCompNoRatio("hptH3L_ppimass", scale, f1, f2, fMc_ld, c,pdf,"M(p#pi)", "p",4,"SE-ME","MC_#Lambda", "p#pi",0);
  projAndScaleComp("hptH3L_ldl", scale,f1, f2, fMc_ld, c,pdf,"l/dl", "p",4,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_ldl", scale,f1, f2, fMc, c,pdf,"l/dl", "p",4,"SE-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_chi2ndf", scale,f1, f2,fMc_ld, c,pdf,"chi2NDF","p",2,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_chi2ndf", scale,f1, f2,fMc, c,pdf,"chi2NDF","p",2,"SE-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_chi2topo", scale,f1, f2, fMc_ld,c,pdf,"chi2topo","p",2,"SE-ME","MC_#Lambda", "dp#pi" ,0,1);
  projAndScaleComp("hptH3L_chi2topo", scale,f1, f2, fMc,c,pdf,"chi2topo","p",2,"SE-ME","MC", "dp#pi" );
  // projAndScaleComp("hptH3L_chi2", scale,f1, f2, fMc_ld,c,pdf,"chi2topo","p",2,"SE-ME","MC_#Lambda", "dp#pi" ,0,1);
  projAndScaleComp("hptH3L_dchi2prim", scale,f1, f2, fMc_ld,c,pdf,"d chi2primary", "p",5,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_dchi2prim", scale,f1, f2, fMc,c,pdf,"d chi2primary", "p",5,"SE-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_pchi2prim", scale,f1, f2, fMc_ld,c,pdf,"p chi2primary", "p",2,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_pchi2prim", scale,f1, f2, fMc,c,pdf,"p chi2primary", "p",2,"SE-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_pichi2prim", scale,f1, f2, fMc_ld,c,pdf,"#pi chi2primary", "p",2,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_pichi2prim", scale,f1, f2, fMc,c,pdf,"#pi chi2primary", "p",2,"SE-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_dDca", scale,f1, f2, fMc_ld,c,pdf,"d DCA", "p",2,"SE-ME","MC_#Lambda","dp#pi"  ,0,1);
  projAndScaleComp("hptH3L_dDca", scale,f1, f2, fMc,c,pdf,"d DCA", "p",2,"SE-ME","MC","dp#pi"  );
  projAndScaleComp("hptH3L_piDca", scale,f1, f2, fMc_ld,c,pdf,"#pi DCA", "p",2,"SE-ME","MC_#Lambda", "dp#pi" ,0,1);
  projAndScaleComp("hptH3L_piDca", scale,f1, f2, fMc,c,pdf,"#pi DCA", "p",2,"SE-ME","MC", "dp#pi" );
  projAndScaleComp("hptH3L_pDca", scale,f1, f2, fMc_ld,c,pdf,"#p DCA", "p",2,"SE-ME","MC_#Lambda", "dp#pi" ,0,1);
  projAndScaleComp("hptH3L_pDca", scale,f1, f2, fMc,c,pdf,"#p DCA", "p",2,"SE-ME","MC", "dp#pi" );
  projAndScaleComp("hptH3L_dpDca", scale,f1, f2, fMc_ld,c,pdf,"dp pair DCA", "p",2,"SE-ME","MC_#Lambda", "dp#pi",0,1 );
  // projAndScaleComp("hptppichi2prim",scale, f1, f2,fMc_ld, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",2 , "SE-ME" ,"MC_#Lambda","p#pi", 0, 1);

  // projAndComp2("hptH3L_chi2topoSig", "hptH3L_chi2topo", fMc_ld, fMc_ld, c,pdf,"(p#pid) #chi^{2}_{topo}" ,"",8 , "Sig" ,"Wide","p#pi");
  // projAndComp2("hptH3L_ldlSig", "hptH3L_ldl", fMc_ld, fMc_ld, c,pdf,"ldl" ,"",8 , "Sig" ,"Wide","p#pi");
  // projAndComp2("hptH3L_lSig", "hptH3L_l", fMc_ld, fMc_ld, c,pdf,"ldl" ,"",8 , "Sig" ,"Wide","p#pi");
  // projAndComp2("hptH3L_lSBR", "hptH3L_l", fMc_ld, fMc_ld, c,pdf,"ldl" ,"",8 , "Sig" ,"Wide","p#pi");
  projAndComp2("hptH3L_ldlSig",  "hptH3L_ldl", fMc, fMc_ld, c,pdf,"ldl" ,"",4 , "Sig" ,"Ld","p#pi");
  projAndComp2("hptH3L_lSig",  "hptH3L_l", fMc, fMc_ld, c,pdf,"l" ,"",4 , "Sig" ,"Ld","p#pi");
  projAndComp2("hptH3L_chi2topoSig",  "hptH3L_chi2topo", fMc, fMc_ld, c,pdf,"chi2topo" ,"",4 , "Sig" ,"Ld","p#pi");
  projAndComp2("hptH3L_chi2ndfSig",  "hptH3L_chi2ndf", fMc, fMc_ld, c,pdf,"chi2ndf" ,"",2 , "Sig" ,"Ld","p#pi");
         
  //add 3 D case
  // cout<<"chi2ndf" << calpurity("h3H3L_chi2ndf",0, 3.5, 1.5, 1, -0.4, -0.6, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 1, 0, 10)<<endl;
  // cout<<"chi2topo" << calpurity("h3H3L_chi2topo",0, 3, 1.5, 1., -0.4, -0.6, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 1, 0, 10)<<endl;
  // cout<<"decaylength" << calpurity("h3H3L_l",8, 99, 1.5, 1, -0.4, -0.6,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 10, 0 ,100 )<<endl;

  // highy = -0.4; lowy = -0.6;
  double error;
  cout<<"chi2ndf" << calpurity("h3H3L_chi2ndf",0, 3.5, highpt, lowpt, highy , lowy, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 1, 0, 5, error, 0.25 )<<endl;
  // // cout<<"chi2ndf" << calpurity("hptH3L_chi2ndf",0, 3.5, 0., 1, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 5, error )<<endl;
  cout<<"chi2ndf" << calpurity("h3H3L_chi2ndf",0, 3.5, 1.5, 1., highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 1, 0, 5, error , 0.25)<<endl;
  cout<<"chi2ndf" << calpurity("h3H3L_chi2ndf",0, 3.5, 2, 1.5, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 1, 0, 5, error , 0.25)<<endl;
  cout<<"chi2ndf" << calpurity("h3H3L_chi2ndf",0, 3.5, 2.5, 2., highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 1, 0, 5, error , 0.25)<<endl;
  cout<<"chi2topo" << calpurity("h3H3L_chi2topo",0, 3, highpt, lowpt, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 1, 0, 10, error , 0.25)<<endl;
  cout<<"chi2topo" << calpurity("h3H3L_chi2topo",0, 3, 1.5, 1., highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 1 , 0, 10, error , 0.25)<<endl;
  cout<<"chi2topo" << calpurity("h3H3L_chi2topo",0, 3, 2, 1.5, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 1 , 0, 10, error , 0.25)<<endl;
  cout<<"chi2topo" << calpurity("h3H3L_chi2topo",0, 3, 2.5, 2., highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 1 , 0, 10, error, 0.25 )<<endl;
  // cout<<"chi2topo" << calpurity("hptH3L_chi2topo",0, 3, lowpt, highpt, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 4 , 0, 50, error )<<endl;
  cout<<"decaylength" << calpurity("h3H3L_l",8, 99,  highpt, lowpt, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 2, 0 ,100 , error , 0.25)<<endl;
  cout<<"decaylength" << calpurity("h3H3L_l",8, 99,  1.5, 1, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 5, 0 ,100 , error , 0.25)<<endl;
  cout<<"decaylength" << calpurity("h3H3L_l",8, 99,  2, 1.5, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 5, 0 ,100 , error , 0.25)<<endl;
  cout<<"decaylength" << calpurity("h3H3L_l",8, 99,  2.5, 2., highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 5, 0 ,100, error , 0.25)<<endl;
  // cout<<"decaylength" << calpurity("hptH3L_l",8, 100,  0.2, 2,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 4 , 0 ,50 )<<endl;
  // cout<<"l/dl" << calpurity("h3H3L_ldl",5, 50, highpt,  lowpt, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l/dl", 2, 3, 50 )<<endl;

  //2D
  // cout<<"chi2ndf" << calpurity("hptH3L_chi2ndf",0, 3.5, 1, 2.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 5)<<endl;
  // cout<<"chi2ndf" << calpurity("hptH3L_chi2ndf",0, 3.5, 0., 1, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 5)<<endl;
  // cout<<"chi2ndf" << calpurity("hptH3L_chi2ndf",0, 3.5, 1., 1.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 5)<<endl;
  // cout<<"chi2ndf" << calpurity("hptH3L_chi2ndf",0, 3.5, 1.5, 2, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 5)<<endl;
  // cout<<"chi2ndf" << calpurity("hptH3L_chi2ndf",0, 3.5, 2., 2.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 5)<<endl;
  // cout<<"chi2topo" << calpurity("hptH3L_chi2topo",0, 3, 1, 2.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 4 , 0, 50)<<endl;
  // cout<<"chi2topo" << calpurity("hptH3L_chi2topo",0, 3, 1, 1.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 4 , 0, 50)<<endl;
  // cout<<"chi2topo" << calpurity("hptH3L_chi2topo",0, 3, 1.5, 2., f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 4 , 0, 50)<<endl;
  // cout<<"chi2topo" << calpurity("hptH3L_chi2topo",0, 3, 2., 2.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 4 , 0, 50)<<endl;
  // // cout<<"chi2topo" << calpurity("hptH3L_chi2topo",0, 3, lowpt, highpt, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 4 , 0, 50)<<endl;
  // cout<<"decaylength" << calpurity("hptH3L_l",8, 100,  lowpt, highpt,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 4 , 0 ,50 )<<endl;
  // // cout<<"decaylength" << calpurity("hptH3L_l",8, 100,  0.2, 2,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 4 , 0 ,50 )<<endl;
  // cout<<"l/dl" << calpurity("hptH3L_ldl",5, 50,  lowpt, highpt,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l/dl", 2, 3, 50 )<<endl;
  // cout<<"dchi2prim" << calpurity("hptH3L_dchi2prim",0, 19,  lowpt, highpt,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L d #chi^{2}",0 )<<endl;
  // cout<<"pichi2prim" << calpurity("hptH3L_pichi2prim",0, 19,  lowpt, highpt,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L pi #chi^{2}",0 )<<endl;
  // cout<<"pip" << calpurity("hptH3L_ppil",0, 39,  lowpt, highpt,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L ppi" , 0)<<endl;
  pdf->On();
  pdf->Close();

}
void drawMixData()
{
  // double highpt = 2.5, lowpt = 1., lowy=-0.8, highy = -0.2;
  // double highpt = 2.5, lowpt = 1, lowy=-0.9, highy = -0.1;
  double highpt = 4, lowpt = 0, lowy=-1., highy = 0;
  TString histname="hH3LMassPtY";
  // TString histname="hH3LMassPtY_5_40";
  SetsPhenixStyle();
  TCanvas* c = new TCanvas("c1","c1");
  // TPDF* pdf = new TPDF("MixEventQA_check.pdf");
  // TPDF* pdf = new TPDF("MixEventQA_Jul27.pdf");
  // TPDF* pdf = new TPDF("MixEventQA_beforeDcacut.pdf");
  TPDF* pdf = new TPDF("MixEventQA_0050.pdf");
  pdf->Off();
  gStyle->SetPalette(1);

  // TFile *f1 = TFile::Open("fout_H3L_data_SE_large.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_large.root"); 

  // TFile *f1 = TFile::Open("out_KF_test/fout_H3L_KF_test.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_KF_Jul22_part.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_KF_Jul25.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_SE_Jul25.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_SE_Jul27.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Aug09.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_SE_Aug16.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_SE_Aug20.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_SE_Aug23.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02_nodca.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02.root"); 
  TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02_0050_015pt.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02_0050_015pt_sys.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02_0050_010pt.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02_nodca_nopt.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_SE_Aug22.root"); 
  // TFile *f1 = TFile::Open("fout_KF_test.root"); 
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
  double nEvents_se = hcent_se->Integral( 4, 9);
  // double nEvents_se = hcent_se->Integral( 1, 9);

  // TFile *f2 = TFile::Open("fout_H3L_data_ME_large.root"); 
  // TFile *f2 = TFile::Open("fout_H3L_ME_Jul22_part.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Aug16.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Aug20.root"); 
  // TFile *f2 = TFile::Open("output2/fout_H3L_data_ME_Aug24.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Sep02.root"); 
  TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Sep02_0050_015pt.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Sep02_0050_015pt_sys.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Sep02_0050_010pt.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Sep02_nodca.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Sep02_nodca_nopt.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Aug23.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Aug22.root"); 
  // TFile *f2 = TFile::Open("fout_ME_test.root"); 
  // TFile *f2 = TFile::Open("fout_H3L_ME_Jul25.root"); 
  // TFile *f2 = TFile::Open("fout_H3L_data_SE.root"); 
  // TH2F* h2bk = (TH2F*)f2->Get("hptH3Lmass")->Clone("hptH3Lmass_ME");
  // TH3F* h2bk = (TH3F*)f2->Get("hH3LMassPtY_5_40")->Clone("hptH3Lmass_bk");
  TH3F* h2bk = (TH3F*)f2->Get(histname.Data())->Clone("hptH3Lmass_bk");
  h2bk->SetDirectory(0);
  TH1F* hbk = (TH1F*)h2bk->ProjectionY("hbk", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2bk->GetZaxis()->FindBin(lowy), h2bk->GetZaxis()->FindBin(highy));
  hbk->SetDirectory(0);

  // TFile* f3 = TFile::Open("fout_H3L_data_RT_Jul22_full.root");
  // TFile* f3 = TFile::Open("fout_H3L_data_RT_Jul25.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_RT_Jul27.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_KF_Aug20.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_SE_Aug23.root");
  // TFile* f3 = TFile::Open("output2/fout_H3L_data_KF_Aug24.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_SE_Sep02_nodca.root");
  TFile* f3 = TFile::Open("rootfile/fout_H3L_data_SE_Sep02_0050_015pt.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_RT_Aug29.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_KF_Aug22.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_KF_Aug16.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_RT_Aug13.root");
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
  double sig_sb =  hsig->Integral(hsig->GetXaxis()->FindBin(3.01),  hsig->GetXaxis()->FindBin(3.02));
  // double rt_sb =  hrt->Integral(hrt->GetXaxis()->FindBin(2.97), hrt->GetXaxis()->FindBin(2.98)) + hrt->Integral(hrt->GetXaxis()->FindBin(3.0),  hrt->GetXaxis()->FindBin(3.02));
  double rt_sb = hrt->Integral(hrt->GetXaxis()->FindBin(3.01),  hrt->GetXaxis()->FindBin(3.02));
  // double rt_sb =  hrt->Integral(hrt->GetXaxis()->FindBin(2.96), hrt->GetXaxis()->FindBin(2.98)) ;
  double scale_rt = sig_sb/rt_sb;
  // sig_sb = hsig->Integral(hsig->GetXaxis()->FindBin(2.97), hsig->GetXaxis()->FindBin(2.98)) + hsig->Integral(hsig->GetXaxis()->FindBin(3.0),  hsig->GetXaxis()->FindBin(3.02));
  double bk_sb =   hbk->Integral(hbk->GetXaxis()->FindBin(3.01),  hbk->GetXaxis()->FindBin(3.02));
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
  hsig->Rebin();
  hbk->Rebin();
  
  setHistStyle(hbk, kRed, kOpenCircle, 1.5);
  hbk->Draw("same");
  setHistStyle(hrt, kGreen+2, kDiamond, 1.5);
  // hrt->Draw("same");
  TLegend* leg_sig = new TLegend(0.65,0.25,0.88,0.45);
  leg_sig->AddEntry(hbk, "ME", "pl");
  // leg_sig->AddEntry(hrt, "RT", "pl");
  // leg_sig->AddEntry(hrt, "KF", "pl");
  leg_sig->AddEntry(hsig, "SE", "pl");
  leg_sig->Draw();

  drawLatex( 0.65,0.61,Form("%0.1f<y<%0.1f",lowy, highy ), 0.055);
  drawLatex( 0.65,0.54,Form("%0.1f<p_{T}<%0.1f GeV/c^{2}", lowpt, highpt), 0.055);
  // drawLatex( 0.65,0.47,Form("5-40%s", "%"), 0.055);
  drawLatex( 0.65,0.47,Form("0-50%s", "%"), 0.055);
  drawBox( 2.97, hsig->GetMinimum(),2.98, hsig->GetMaximum()*0.5 , kBlue-9, 1001 , 1,0.3); 
  // drawBox( 2.989, hsig_bk->GetMinimum(),2.995, hsig_bk->GetMaximum()*0.5 , kOrange-4, 1001,1,0.3); 
  drawBox( 3.01, hsig->GetMinimum(),3.02, hsig->GetMaximum()*0.5 , kBlue-9, 1001,1,0.3 ); 

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
  // hsig_rt->Draw("same");
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
  // hsig_rt->Fit(resfit_rt,"R");
  fit_rt->SetParameter(3, resfit_rt->GetParameter(0) );
  fit_rt->SetParameter(4, resfit_rt->GetParameter(1) );
  // hsig_rt->Fit(fit_rt,"R");
  resfit_rt->SetParameter(0, fit_rt->GetParameter(3));
  resfit_rt->SetParameter(1, fit_rt->GetParameter(4));
  // resfit_rt->Draw("same");

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
  double significance = yield_me/sqrt(sp_counts);
  // double significance = yield_rt/sqrt(sp_counts);
  double s_me = yield_me/(fit->GetParError(0)/hsig_bk->GetBinWidth(1));
  double s_rt = yield_rt/(fit_rt->GetParError(0)/hsig_rt->GetBinWidth(1));
  cout<<"significance: " <<significance << endl;

  TLegend* leg = new TLegend( 0.72, 0.68 ,0.9,0.9 );
  // leg->AddEntry(hbk, "mix-event(ME)","pl");
  // leg->AddEntry(hrt, "rotate d(RT) (scale)","pl");
  leg->AddEntry(hsig, "SE","pl");
  leg->AddEntry(hsig_bk, "SE-ME","pl");
  // leg->AddEntry(hsig_rt, "SE-KF","pl");
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
  drawBox( 3.01, hsig_bk->GetMinimum(),3.02, hsig_bk->GetMaximum()*0.5 , kBlue-9, 1001,1,0.3 ); 

  addpdf(pdf);

  hsig_bk->Draw();
  // hsig_rt->Draw("same");
  TLegend* leg2 = new TLegend( 0.2 , 0.7 ,0.4,0.9  );
  leg2->AddEntry(hsig_rt, "SE-RT","pl");
  leg2->AddEntry(hsig_bk, "SE-ME","pl");
  leg2->Draw();

  addpdf(pdf);
  //
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

  TH2F* hYield = new TH2F( "hYield", "hYield;y,pt", 5, -1., 0., 6, 0, 3);
  double xbw = hYield->GetXaxis()->GetBinWidth(1);
  double ybw = hYield->GetYaxis()->GetBinWidth(1);
  for (int ix=1; ix<=hYield->GetNbinsX(); ix++) {
    for (int iy=1; iy<=hYield->GetNbinsY();iy++){
      double xlow = hYield->GetXaxis()->GetBinLowEdge(ix); // y
      double ylow = hYield->GetYaxis()->GetBinLowEdge(iy); // pt
      double err;
      hYield->SetBinContent( ix , iy, fityield( ylow, ylow+ybw, xlow, xlow+xbw, err, f1, f2, c, pdf) );
      hYield->SetBinError( ix , iy, err);
    }
  }
  hYield->Draw("colz text");
  hYield->SaveAs("fH3L_yield_0050.root");
  addpdf(pdf);

  TFile* fMc_ld = TFile::Open("fout_Lambda_MC_Cuts_0050_015pt.root");
  // TFile* fMc_ld = TFile::Open("fout_Lambda_MC_Cuts_0050_015pt_sys_pid3.root");
  TFile* fMc = TFile::Open("fout_H3L_MC_0050_015pt.root");
  // TFile* fMc = TFile::Open("fout_H3L_MC_0050_015pt_sys.root");
  // TFile* fMc = TFile::Open("fout_H3L_MC_0080_010pt.root");
  TH1F* hPhase[3];
  TH1F* hPhaseCor[3];
  TH1F* hPurityCor[3];
  TH1F* hBr[3];
  TH1F* hBrRaw[3];
  TGraphErrors* g[3];
  double edge[3][4]={
    { 1.7, 2.1, 2.7, 3.1},
    { 1.2, 1.6, 2.0, 2.4},
    { 1.2, 1.8, 2.2, 2.6}
  };

  double ybin[4]={ 0, -0.25, -0.50, -0.75};
  hPhase[0] = new TH1F("hPhase0", "hPhase0;pt", 25, 1, 3.5);
  hPhase[0]=(TH1F*)hPhase[0]->Rebin(3, "hPhase0", edge[0]);
  hPhase[1] = new TH1F("hPhase1", "hPhase1;pt", 25, 1, 3.5);
  hPhase[1]=(TH1F*)hPhase[1]->Rebin(3, "hPhase1", edge[1]);
  hPhase[2] = new TH1F("hPhase2", "hPhase2;pt", 25, 1, 3.5);
  hPhase[2]=(TH1F*)hPhase[2]->Rebin(3, "hPhase2", edge[2]);
  TFile* fMcH3L = new TFile("fMC_H3L_0050.root");
  TFile* fRcH3L = new TFile("fout_H3L_MC_0050_015pt.root");
  // TFile* fRcH3L = new TFile("fout_H3L_MC_0050_015pt_sys.root");
  // TFile* fMcH3L = new TFile("fMC_H3L_0080.root");
  // TFile* fRcH3L = new TFile("fout_H3L_MC_0080_010pt.root");
  TH3F* h3Mc = (TH3F*)fMcH3L->Get("hH3LMassPtY")->Clone("h3Mc");
  h3Mc->SetDirectory(0);
  TH3F* h3Rc = (TH3F*)fRcH3L->Get("hH3LMassPtY")->Clone("h3Rc");
  h3Rc->SetDirectory(0);

  h3Mc->Sumw2();
  h3Rc->Sumw2();
  TH2F* h2MC = (TH2F*)h3Mc->Project3D("xz");
  TH2F* h2Rc = (TH2F*)h3Rc->Project3D("xz");

  TH2F* h2Eff = (TH2F*)h2Rc->Clone("h2Eff");
  TH2F* h2temp = (TH2F*)h2MC->Clone("h2temp");
  h2temp->RebinY(10);
  h2temp->RebinX(10);
  h2Eff->RebinY(10);
  h2Eff->RebinX(10);
  h2Eff->Divide(h2temp);
  h2Eff->Draw("col text");
  h2Eff->GetYaxis()->SetRangeUser(0,4.5);
  addpdf(pdf);

  TH1F* heff[3];
  for (int ij=0;ij<3;ij++){ 
    TH1F* h1Mc = (TH1F*)h2MC->ProjectionY(Form("h1Mc%d", ij), h2MC->GetXaxis()->FindBin(ybin[ij+1] +1e-6), h2MC->GetXaxis()->FindBin(ybin[ij]-1e-6) );
    h1Mc = (TH1F*)h1Mc->Rebin( 3, Form("hMc%d",ij), edge[ij]);
    TH1F* h1Rc = (TH1F*)h2Rc->ProjectionY(Form("h1Rc%d", ij), h2Rc->GetXaxis()->FindBin(ybin[ij+1]+1e-6), h2Rc->GetXaxis()->FindBin(ybin[ij]-1e-6) );
    heff[ij] = (TH1F*)h1Rc->Rebin( 3, Form("heff%d",ij), edge[ij]);
    heff[ij]->Divide(h1Mc);
    heff[ij]->GetYaxis()->SetTitle("Eff.");
    heff[ij]->Draw();
    addpdf(pdf);
  }

  double ptesterr[7];
  double chi2topo[7]={0.2,0.6,1,1.5,2,2.5,3 };
  double puritytest[7];
  for (int i=0;i<7;i++) {
    puritytest[i] = calpurity("h3H3L_chi2topo",0,chi2topo[i],  3., 1.5, -0.25, -0.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, ptesterr[i], 0.25);
    // puritytest[i] = calpurity("h3H3L_chi2topo",0, 3,  3., 1., -0., -0.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, ptesterr[i], 0.25);
  }
  for (int i=0;i<7;i++) 
    cout <<chi2topo[i]<<" "<<puritytest[i]<<" "<<ptesterr[i] << endl;
  TGraphErrors* gpurity = new TGraphErrors(7, chi2topo , puritytest,0, ptesterr  );
  gpurity->Draw("pA");
  gpurity->GetXaxis()->SetTitle("Chi2Topo");
  gpurity->GetYaxis()->SetTitle("Purity");
  addpdf(pdf);

  double dy=0.25;
  double dpt,pt;
  
  TFile* f2b = TFile::Open("h_h3l_corr_yield.root_cut00002_cent0_yuehang.root");

  double ny = 3;
  for (int ij=0; ij<ny;ij++){
    double npt =3;
    if (ij==0) npt =2; 
    for (int ipt=0; ipt<npt;ipt++){
      double err;
      dpt = edge[ij][ipt+1]-edge[ij][ipt];
      pt = 0.5*(edge[ij][ipt+1]+edge[ij][ipt]);
      double yield = fityield( edge[ij][ipt], edge[ij][ipt+1], ybin[ij+1], ybin[ij], err, f1, f2, c, pdf);
      // cout << edge[ij][ipt]<<" "<<edge[ij][ipt+1]<<" "<<ybin[ij+1]<<" "<<ybin[ij]<<" "<< yield<<" "<<err <<endl;
      if (ij==0 && ipt==2) { yield=0; err=0;}
      hPhase[ij]->SetBinContent( ipt+1, yield/dy/dpt/pt/2./3.1415926/nEvents_se);
      hPhase[ij]->SetBinError( ipt+1, err/dy/dpt/pt/2./3.1415926/nEvents_se );
    }
    hPhase[ij]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPhase[ij]->GetYaxis()->SetTitle("Raw d^{2}N/(N_{evt}2#pip_{T}dydp_{T})");
    hPhase[ij]->Draw();
    hPhase[ij]->SetDirectory(0);
    addpdf(pdf); 

    hPhaseCor[ij]=(TH1F*)hPhase[ij]->Clone(Form("hYieldCor_%d",ij));
    hPhaseCor[ij]->SetDirectory(0);
    hPurityCor[ij]=(TH1F*)hPhase[ij]->Clone(Form("hPurityCor_%d",ij));
    hPurityCor[ij]->SetDirectory(0);

    double perr;
    double purity1_s = calpurity("h3H3L_chi2ndf",0, 3.5, edge[ij][3], edge[ij][0], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 10, perr, 0.25);
    // double purity2_s = calpurity("h3H3L_chi2topo",0, 3.,  edge[ij][3],  edge[ij][0], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
    double purity2_s = calpurity("h3H3L_chi2topo",0, 3.,  3.,  1.5, ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
    // double puritytest = calpurity("h3H3L_chi2topo",0, 3.,  3.,  1.,ybin[0], ybin[2], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
    // double purity3_s = calpurity("h3H3L_l",8, 99,  3.,  1., ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L decaylength", 2, 0, 100, perr, 0.3);

    for (int ipt=0;ipt<npt;ipt++)
    {
      double y3b = hPhaseCor[ij]->GetBinContent(ipt+1);
      double y3berr = hPhaseCor[ij]->GetBinError(ipt+1)/y3b;
      double eff = heff[ij]->GetBinContent(ipt+1);
      double efferr = heff[ij]->GetBinError(ipt+1)/eff;
      // double purity1 = calpurity("h3H3L_chi2ndf",0, 3.5, edge[ij][ipt+1], edge[ij][ipt], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 10, perr, 0.25);
      // double purity2 = calpurity("h3H3L_chi2topo",0, 3., edge[ij][ipt+1], edge[ij][ipt], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
      // double purity = (purity1+purity2)*0.5;
      // perr = fabs(purity2-purity1);
      // perr = perr/purity;
      double purity = purity2_s;
      // double purity;
      // purity=1;
      // perr = 0;
      double yield_cor = y3b/eff;
      double err = sqrt(y3berr*y3berr + efferr*efferr)*yield_cor;
      hPhaseCor[ij]->SetBinContent(ipt+1,yield_cor );
      hPhaseCor[ij]->SetBinError(ipt+1, err);
      yield_cor = y3b/eff*purity;
      err = sqrt(y3berr*y3berr + efferr*efferr + perr*perr )*yield_cor;
      hPurityCor[ij]->SetBinContent(ipt+1,yield_cor );
      hPurityCor[ij]->SetBinError(ipt+1, err);

      cout <<"pt "<<hPhaseCor[ij]->GetBinCenter(ipt+1) << "purity "<<purity <<" eff " <<eff <<" yield " << yield_cor << " err "<< err<< " s="<<yield_cor/err<<endl;
    }
    hPhaseCor[ij]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPhaseCor[ij]->GetYaxis()->SetTitle("d^{2}N/(N_{evt}2#pip_{T}dydp_{T}(no purity corr.))");
    hPurityCor[ij]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPurityCor[ij]->GetYaxis()->SetTitle("d^{2}N/(N_{evt}2#pip_{T}dydp_{T})");

    hPhaseCor[ij]->Draw();
    hPhaseCor[ij]->SetDirectory(0);
    addpdf(pdf);

    hBrRaw[ij]=(TH1F*)hPhaseCor[ij]->Clone(Form("hBrRaw_%d",ij));
    hBrRaw[ij]->SetDirectory(0);
    g[ij] = (TGraphErrors*)f2b->Get(Form("t_sgct1_corr_yield[0][%d]", ij));
    // g[ij]->Draw();
    for (int ipt=0;ipt<npt;ipt++)
    {
      double y3b = hBrRaw[ij]->GetBinContent(ipt+1);
      double y3berr = hBrRaw[ij]->GetBinError(ipt+1);
      double y2b = g[ij]->GetPointY(ipt);
      double y2bpt = g[ij]->GetPointX(ipt);
      double y2berr = g[ij]->GetErrorY(ipt);
      if (ij==1) {y2berr=y2berr*10; y2b=y2b*10.;}
      if (ij==2) {y2berr=y2berr*100; y2b=y2b*100.;}
      // cout <<"y2b significance:"<<y2b/y2berr<<" yield "<<y2b << endl;
      hBrRaw[ij]->SetBinContent(ipt+1, y2b/(y2b+y3b) );
      hBrRaw[ij]->SetBinError(ipt+1, sqrt((y2b*y2b)/pow((y2b+y3b), 4)*(y3berr*y3berr) + (y3b*y3b)/pow((y2b+y3b), 4)*(y2berr*y2berr)));
      if (y2b/y2berr<2) {
         hBrRaw[ij]->SetBinContent(ipt+1, 0);
         hBrRaw[ij]->SetBinError(ipt+1, 0);
      } 
      // cout <<"bincenter: " << hBr[ij]->GetBinCenter(ipt+1) << " content:"<< hBr[ij]->GetBinContent(ipt+1)<<" err:"<<hBr[ij]->GetBinError(ipt+1)<< "y2b significance:"<<y2b/y2berr<<" yield "<<y2b << " y2bpt "<< y2bpt<< " y3b yield "<< y3b<<endl;
    } 
    hBrRaw[ij]->Draw();
    hBrRaw[ij]->GetYaxis()->SetTitle("R_{3} Raw");
    drawLatex( 0.65,0.88 , Form("%0.2f<y<%0.2f", ybin[ij+1],ybin[ij]) , 0.055);
    addpdf(pdf);

    hBr[ij]=(TH1F*)hPurityCor[ij]->Clone(Form("hBr_%d",ij));
    hBr[ij]->SetDirectory(0);
    g[ij] = (TGraphErrors*)f2b->Get(Form("t_sgct1_corr_yield[0][%d]", ij));
    g[ij]->Draw();
    for (int ipt=0;ipt<npt;ipt++)
    {
      double y3b = hBr[ij]->GetBinContent(ipt+1);
      double y3berr = hBr[ij]->GetBinError(ipt+1);
      double y2b = g[ij]->GetPointY(ipt);
      double y2bpt = g[ij]->GetPointX(ipt);
      double y2berr = g[ij]->GetErrorY(ipt);
      if (ij==1) {y2berr=y2berr*10; y2b=y2b*10.;}
      if (ij==2) {y2berr=y2berr*100; y2b=y2b*100.;}
      // cout <<"y2b significance:"<<y2b/y2berr<<" yield "<<y2b << endl;
      hBr[ij]->SetBinContent(ipt+1, y2b/(y2b+y3b) );
      hBr[ij]->SetBinError(ipt+1, sqrt((y2b*y2b)/pow((y2b+y3b), 4)*(y3berr*y3berr) + (y3b*y3b)/pow((y2b+y3b), 4)*(y2berr*y2berr)));
      if (y2b/y2berr<2) {
         hBr[ij]->SetBinContent(ipt+1, 0);
         hBr[ij]->SetBinError(ipt+1, 0);
      } 
      cout <<"bincenter: " << hBr[ij]->GetBinCenter(ipt+1) << " content:"<< hBr[ij]->GetBinContent(ipt+1)<<" err:"<<hBr[ij]->GetBinError(ipt+1)<< "y2b significance:"<<y2b/y2berr<<" yield "<<y2b << " y2bpt "<< y2bpt<< " y3b yield "<< y3b<<endl;
    } 
    hBr[ij]->Draw();
    hBr[ij]->GetYaxis()->SetTitle("R_{3}");
    drawLatex( 0.65,0.88 , Form("%0.2f<y<%0.2f", ybin[ij+1],ybin[ij]) , 0.055);
    addpdf(pdf);
  }

  TLegend* legb = new TLegend(0.7, 0.7, 0.9, 0.9);
  TH1F* htemp = new TH1F("htemp","htemp;p_{T} (GeV/c); 2-body/(2-body+3-body)",  5, 1, 3);
  htemp->GetYaxis()->SetRangeUser( 0, 0.5);
  htemp->Draw();
  int color[]={kRed-4, kGreen+2, kBlue-4 };
  for (int i=0;i<3;i++)
  {
    hBr[i]->Draw("same");
    hBr[i]->SetLineColor(color[i]);
    hBr[i]->SetMarkerColor(color[i]);
    legb->AddEntry( hBr[i], Form("%0.2f<y<%0.2f", ybin[i+1], ybin[i]), "pl");
  }
  legb->Draw();
  addpdf(pdf);

  double tmp=0, tmperr = 0;
  double sum=0, wtsum=0, sumerr=0;
  for (int i=0;i<2;i++)
  {
    for (int ip=1;ip<=hBr[i]->GetNbinsX();ip++)
    {
      tmp = hBr[i]->GetBinContent(ip);
      if (tmp<1e-31) continue; 
      tmperr = hBr[i]->GetBinError(ip);
      cout <<tmperr <<" "<<tmp << endl;
      sum += tmp/tmperr/tmperr;
      sumerr += 1./(tmperr*tmperr);
      wtsum++;
    }
  } 
  cout<<"after correct:" << sum/sumerr<<" "<<sqrt(1./sumerr)<<endl;

  tmp=0; tmperr = 0;
  sum=0; wtsum=0; sumerr=0;
  for (int i=0;i<2;i++)
  {
    for (int ip=1;ip<=hBr[i]->GetNbinsX();ip++)
    {
      tmp = hBrRaw[i]->GetBinContent(ip);
      if (tmp<1e-31) continue; 
      tmperr = hBrRaw[i]->GetBinError(ip);
      cout <<tmperr <<" "<<tmp << endl;
      sum += tmp/tmperr/tmperr;
      sumerr += 1./(tmperr*tmperr);
      wtsum++;
    }
  } 
  cout<<"Before correct:" << sum/sumerr<<" "<<sqrt(1./sumerr)<<endl;

  TFile * fout = new TFile("fout_0050.root","recreate");
  fout->cd();
  for (int i=0;i<3;i++){
    hPhase[i]->Write();
    hPhaseCor[i]->Write();
    hPurityCor[i]->Write();
    hBr[i]->Write();
    heff[i]->Write();
  }
  fout->Close();

  // pdf->On();
  // pdf->Close();
  //
  // return;
  c = new TCanvas( "c","c", 800,800);
  projAndScaleComp("hptH3L_ppichi2prim", scale,f1, f2, fMc_ld,c,pdf,"p#pi chi2primary", "p",2,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_ppichi2prim", scale,f1, f2, fMc,c,pdf,"p#pi chi2primary", "p",2,"SE-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_ppimass", scale,f1, f2, fMc,c,pdf,"p#pi chi2primary", "p",2,"SE-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_ppimass", scale,f1, f2, fMc_ld,c,pdf,"p#pi chi2primary", "p",2,"SE-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_l", scale, f1, f2, fMc_ld, c,pdf,"l", "p",8,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_l", scale, f1, f2, fMc, c,pdf,"l", "p",8,"SE-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_ldl", scale, f1, f2, fMc_ld, c,pdf,"l/dl", "p",4,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_ldl", scale, f1, f2, fMc, c,pdf,"l/dl", "p",4,"SE-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_chi2ndf", scale, f1, f2,fMc_ld, c,pdf,"chi2NDF","p",2,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_chi2ndf", scale, f1, f2,fMc, c,pdf,"chi2NDF","p",2,"SE-ME","MC", "dp#pi", 0, 1);
  projAndScaleComp("hptH3L_chi2topo", scale, f1, f2, fMc_ld,c,pdf,"chi2topo","p",2,"SE-ME","MC_#Lambda", "dp#pi" ,0,1);
  projAndScaleComp("hptH3L_chi2topo", scale, f1, f2, fMc,c,pdf,"chi2topo","p",2,"SE-ME","MC", "dp#pi" );
  // projAndComp("hptH3L_lSBL", f1, f2, c,pdf,"l", "p",4,"SE_SBL","ME_SBL", "dp#pi", scale);
  // projAndComp("hptH3L_lSBR", f1, f2, c,pdf,"l", "p",4,"SE_SBR","ME_SBR", "dp#pi", scale);
  // projAndComp("hptH3L_ldlSBL", f1, f2, c,pdf,"l/dl", "p",4,"SE_SBL","ME_SBL", "dp#pi", scale);
  // projAndComp("hptH3L_ldlSBR", f1, f2, c,pdf,"l/dl", "p",4,"SE_SBR","ME_SBR", "dp#pi", scale);
  // projAndComp("hptH3L_chi2ndfSBL", f1, f2, c,pdf,"chi2NDF", "p",4,"SE_SBL","ME_SBL", "dp#pi", scale);
  // projAndComp("hptH3L_chi2ndfSBR", f1, f2, c,pdf,"chi2NDF", "p",4,"SE_SBR","ME_SBR", "dp#pi", scale);
  // projAndComp("hptH3L_chi2topoSBL", f1, f2, c,pdf,"chi2topo", "p",4,"SE_SBL","ME_SBL", "dp#pi" , scale);
  // projAndComp("hptH3L_chi2topoSBR", f1, f2, c,pdf,"chi2topo", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  projAndComp2("hptH3L_ldlSig",  "hptH3L_ldl", fMc, fMc_ld, c,pdf,"ldl" ,"",4 , "Sig" ,"Ld","p#pi");
  projAndComp2("hptH3L_lSig",  "hptH3L_l", fMc, fMc_ld, c,pdf,"l" ,"",4 , "Sig" ,"Ld","p#pi");
  projAndComp2("hptH3L_chi2topoSig",  "hptH3L_chi2topo", fMc, fMc_ld, c,pdf,"chi2topo" ,"",4 , "Sig" ,"Ld","p#pi");
  projAndComp2("hptH3L_chi2ndfSig",  "hptH3L_chi2ndf", fMc, fMc_ld, c,pdf,"chi2ndf" ,"",2 , "Sig" ,"Ld","p#pi");
  TFile* fMc_ph = TFile::Open("fout_H3L_phase_MC.root");
   projAndScaleComp("hptH3L_l", scale,f1, f2, fMc_ph, c,pdf,"l", "p",4,"SE-ME","phase_MC", "dp#pi");
  projAndScaleComp("hptH3L_ldl", scale,f1, f2, fMc_ph, c,pdf,"l/dl", "p",4,"SE-ME","phase_MC", "dp#pi");
  projAndScaleComp("hptH3L_chi2ndf", scale,f1, f2,fMc_ph, c,pdf,"chi2NDF","p",2,"SE-ME","phase_MC", "dp#pi");
  projAndScaleComp("hptH3L_chi2topo", scale,f1, f2, fMc_ph,c,pdf,"chi2topo","p",2,"SE-ME","phase_MC", "dp#pi" );
  projAndComp2("hptH3L_ldlSig",  "hptH3L_ldlSig", fMc_ph, fMc, c,pdf,"ldl" ,"",4 , "phase" ,"quasi","p#pi");
  projAndComp2("hptH3L_lSig",  "hptH3L_lSig", fMc_ph, fMc, c,pdf,"l" ,"",4 , "phase" ,"quasi","p#pi");
  projAndComp2("hptH3L_chi2topoSig",  "hptH3L_chi2topoSig", fMc_ph, fMc, c,pdf,"chi2topo" ,"",4 , "phase" ,"quasi","p#pi");
  projAndComp2("hptH3L_chi2ndfSig",  "hptH3L_chi2ndfSig", fMc_ph, fMc, c,pdf,"chi2ndf" ,"",0 , "phase" ,"quasi","p#pi");
  pdf->On();
  pdf->Close();

  return;


  projAndComp("hptH3L_lSBL", f1, f2, c,pdf,"l", "p",4,"SE_SBL","ME_SBL", "dp#pi", scale);
  projAndComp("hptH3L_lSBR", f1, f2, c,pdf,"l", "p",4,"SE_SBR","ME_SBR", "dp#pi", scale);
  projAndComp("hptH3L_ldlSBL", f1, f2, c,pdf,"l/dl", "p",4,"SE_SBL","ME_SBL", "dp#pi", scale);
  projAndComp("hptH3L_ldlSBR", f1, f2, c,pdf,"l/dl", "p",4,"SE_SBR","ME_SBR", "dp#pi", scale);
  projAndComp("hptH3L_chi2ndfSBL", f1, f2, c,pdf,"chi2NDF", "p",4,"SE_SBL","ME_SBL", "dp#pi", scale);
  projAndComp("hptH3L_chi2ndfSBR", f1, f2, c,pdf,"chi2NDF", "p",4,"SE_SBR","ME_SBR", "dp#pi", scale);
  projAndComp("hptH3L_chi2topoSBL", f1, f2, c,pdf,"chi2topo", "p",4,"SE_SBL","ME_SBL", "dp#pi" , scale);
  projAndComp("hptH3L_chi2topoSBR", f1, f2, c,pdf,"chi2topo", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  projAndComp("hptH3L_dchi2primSBL", f1, f2, c,pdf,"d chi2primary", "p",4,"SE_SBL","ME_SBL", "dp#pi", scale);
  projAndComp("hptH3L_dchi2primSBR", f1, f2, c,pdf,"d chi2primary", "p",4,"SE_SBR","ME_SBR", "dp#pi", scale);
  projAndComp("hptH3L_pchi2primSBL", f1, f2, c,pdf,"p chi2primary", "p",4,"SE_SBL","ME_SBL", "dp#pi", scale);
  projAndComp("hptH3L_pchi2primSBR", f1, f2, c,pdf,"p chi2primary", "p",4,"SE_SBR","ME_SBR", "dp#pi", scale);
  projAndComp("hptH3L_pichi2primSBL", f1, f2, c,pdf,"#pi chi2primary", "p",4,"SE_SBL","ME_SBL", "dp#pi", scale);
  projAndComp("hptH3L_pichi2primSBR", f1, f2, c,pdf,"#pi chi2primary", "p",4,"SE_SBR","ME_SBR", "dp#pi", scale);
  projAndComp("hptH3L_dDcaSBL", f1, f2, c,pdf,"d DCA", "p",1,"SE_SBL","ME_SBL","dp#pi"  , scale);
  projAndComp("hptH3L_dDcaSBR", f1, f2, c,pdf,"d DCA", "p",1,"SE_SBR","ME_SBR","dp#pi"  , scale);
  projAndComp("hptH3L_piDcaSBL", f1, f2, c,pdf,"#pi DCA", "p",4,"SE_SBL","ME_SBL", "dp#pi" , scale);
  projAndComp("hptH3L_piDcaSBR", f1, f2, c,pdf,"#pi DCA", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  projAndComp("hptH3L_pDcaSBL", f1, f2, c,pdf,"#p DCA", "p",4,"SE_SBL","ME_SBL", "dp#pi" , scale);
  projAndComp("hptH3L_pDcaSBR", f1, f2, c,pdf,"#p DCA", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  projAndComp("hptH3L_dpDcaSBL", f1, f2, c,pdf,"dp pair DCA", "p",4,"SE_SBL","ME_SBL", "dp#pi" , scale);
  projAndComp("hptH3L_dpDcaSBR", f1, f2, c,pdf,"dp pair DCA", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  // projAndComp("hptppimassSBL", f1, f2, c,pdf,"p#pi Mass", "p",4,"SE_SBL","ME_SBL", "dp#pi" );
  // projAndComp("hptH3L_ppimassSBR", f1, f2, c,pdf,"p#pi Mass", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  // projAndComp("hptH3L_ppichi2ndfSBR", f1, f2, c,pdf,"p#pi Mass", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  // projAndComp("hptH3L_ppichi2ndfSBL", f1, f2, c,pdf,"p#pi chi2ndf", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  // projAndComp("hptH3L_ppichi2primSBL", f1, f2, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "SE" ,"BK","p#pi",scale);

  // projAndComp("hptH3L_lSBR", f1, f3, c,pdf,"l", "p",4,"SE_SBR","RT_SBR", "dp#pi");
  // projAndComp("hptH3L_ldlSBR", f1, f3, c,pdf,"l/dl", "p",4,"SE_SBR","RT_SBR", "dp#pi");
  // projAndComp("hptH3L_chi2ndfSBR", f1, f3, c,pdf,"chi2NDF", "p",4,"SE_SBR","RT_SBR", "dp#pi");
  // projAndComp("hptH3L_chi2topoSBR", f1, f3, c,pdf,"chi2topo", "p",4,"SE_SBR","RT_SBR", "dp#pi" );
  // projAndComp("hptH3L_dchi2primSBR", f1, f3, c,pdf,"d chi2primary", "p",4,"SE_SBR","RT_SBR", "dp#pi");
  // projAndComp("hptH3L_pchi2primSBR", f1, f3, c,pdf,"p chi2primary", "p",4,"SE_SBR","RT_SBR", "dp#pi");
  // projAndComp("hptH3L_pichi2primSBR", f1, f3, c,pdf,"#pi chi2primary", "p",4,"SE_SBR","RT_SBR", "dp#pi");
  // projAndComp("hptH3L_dDcaSBR", f1, f3, c,pdf,"d DCA", "p",1,"SE_SBR","RT_SBR","dp#pi"  );
  // projAndComp("hptH3L_piDcaSBR", f1, f3, c,pdf,"#pi DCA", "p",4,"SE_SBR","RT_SBR", "dp#pi" );
  // projAndComp("hptH3L_pDcaSBR", f1, f3, c,pdf,"#p DCA", "p",4,"SE_SBR","RT_SBR", "dp#pi" );
  // projAndComp("hptH3L_dpDcaSBR", f1, f3, c,pdf,"dp pair DCA", "p",4,"SE_SBR","RT_SBR", "dp#pi" );
  // projAndComp("hptH3L_ppimassSBR", f1, f3, c,pdf,"p#pi Mass (GeV/c^{2})", "plhist",1 ,"SE" ,"BK","p#pi");
  // // projAndComp("hptppichi2ndfSBR", f1, f2, c,pdf,"p#pi #chi^{2}_{NDF}","" ,1 , "SE" ,"BK","p#pi");
  // // projAndComp("hptppichi2primSBR", f1, f2, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "SE" ,"BK","p#pi");
  // // projAndComp("hptppilSBR", f1, f2, c,pdf,"(p#pi) l" ,"",1 , "SE" ,"BK","p#pi");
  // // projAndComp("hptppildlSBR", f1, f2, c,pdf,"(p#pi) l/#deltal" ,"",1 , "SE" ,"BK","p#pi");

  // projAndComp("hptH3L_lSBL", f1, f3, c,pdf,"l", "p",4,"SE_SBL","RT_SBL", "dp#pi");
  // projAndComp("hptH3L_ldlSBL", f1, f3, c,pdf,"l/dl", "p",4,"SE_SBL","RT_SBL", "dp#pi");
  // projAndComp("hptH3L_chi2ndfSBL", f1, f3, c,pdf,"chi2NDF", "p",4,"SE_SBL","RT_SBL", "dp#pi");
  // projAndComp("hptH3L_chi2topoSBL", f1, f3, c,pdf,"chi2topo", "p",4,"SE_SBL","RT_SBL", "dp#pi" );
  // projAndComp("hptH3L_dchi2primSBL", f1, f3, c,pdf,"d chi2primary", "p",4,"SE_SBL","RT_SBL", "dp#pi");
  // projAndComp("hptH3L_pchi2primSBL", f1, f3, c,pdf,"p chi2primary", "p",4,"SE_SBL","RT_SBL", "dp#pi");
  // projAndComp("hptH3L_pichi2primSBL", f1, f3, c,pdf,"#pi chi2primary", "p",4,"SE_SBL","RT_SBL", "dp#pi");
  // projAndComp("hptH3L_dDcaSBL", f1, f3, c,pdf,"d DCA", "p",1,"SE_SBL","RT_SBL","dp#pi"  );
  // projAndComp("hptH3L_piDcaSBL", f1, f3, c,pdf,"#pi DCA", "p",4,"SE_SBL","RT_SBL", "dp#pi" );
  // projAndComp("hptH3L_pDcaSBL", f1, f3, c,pdf,"#p DCA", "p",4,"SE_SBL","RT_SBL", "dp#pi" );
  // projAndComp("hptH3L_dpDcaSBL", f1, f3, c,pdf,"dp pair DCA", "p",4,"SE_SBL","RT_SBL", "dp#pi" );
  // // projAndComp("hptppimassSBL", f1, f3, c,pdf,"p#pi Mass (GeV/c^{2})", "plhist",1 ,"SE" ,"BK","p#pi");
  // projAndComp("hptH3L_ppimassSBL", f1, f3, c,pdf,"p#pi Mass (GeV/c^{2})", "plhist",1 ,"SE" ,"BK","dp#pi");
  // // projAndComp("hptppichi2ndfSBL", f1, f2, c,pdf,"p#pi #chi^{2}_{NDF}","" ,1 , "SE" ,"BK","p#pi");
  // // projAndComp("hptppichi2primSBL", f1, f2, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "SE" ,"BK","p#pi");
  // // projAndComp("hptppilSBL", f1, f2, c,pdf,"(p#pi) l" ,"",1 , "SE" ,"BK","p#pi");
  // // projAndComp("hptppildlSBL", f1, f2, c,pdf,"(p#pi) l/#deltal" ,"",1 , "SE" ,"BK","p#pi");

  // projAndScaleComp("hptH3L_l", scale,f1, f2, f2, c,pdf,"l", "p",1,"SE-ME","ME", "dp#pi");
  // projAndScaleComp("hptH3L_ppimass", scale, f1, f2, f2, c,pdf,"M(p#pi)", "p",1,"SE-ME","ME", "p#pi");
  // // projAndScaleCompNoRatio("hptppimass", scale, f1, f2, f2, c,pdf,"M(p#pi)", "p",1,"SE-ME","ME", "p#pi");
  // projAndScaleComp("hptH3L_ldl", scale,f1, f2, f2, c,pdf,"l/dl", "p",1,"SE-ME","ME", "dp#pi");
  // projAndScaleComp("hptH3L_chi2ndf", scale,f1, f2,f2, c,pdf,"chi2NDF", "p",1,"SE-ME","ME", "dp#pi");
  projAndScaleComp("hptH3L_chi2topo", scale,f1, f2, f2,c,pdf,"chi2topo", "p",1,"SE-ME","ME", "dp#pi" );
  // projAndScaleComp("hptH3L_dchi2prim", scale,f1, f2, f2,c,pdf,"d chi2primary", "p",1,"SE-ME","ME", "dp#pi");
  // projAndScaleComp("hptH3L_pchi2prim", scale,f1, f2, f2,c,pdf,"p chi2primary", "p",1,"SE-ME","ME", "dp#pi");
  // projAndScaleComp("hptH3L_pichi2prim", scale,f1, f2, f2,c,pdf,"#pi chi2primary", "p",1,"SE-ME","ME", "dp#pi");
  // projAndScaleComp("hptH3L_dDca", scale,f1, f2, f2,c,pdf,"d DCA", "p",1,"SE-ME","ME","dp#pi"  );
  // projAndScaleComp("hptH3L_piDca", scale,f1, f2, f2,c,pdf,"#pi DCA", "p",1,"SE-ME","ME", "dp#pi" );
  // projAndScaleComp("hptH3L_pDca", scale,f1, f2, f2,c,pdf,"#p DCA", "p",1,"SE-ME","ME", "dp#pi" );
  // projAndScaleComp("hptH3L_dpDca", scale,f1, f2, f2,c,pdf,"dp pair DCA", "p",1,"SE-ME","ME", "dp#pi" );
  // // projAndScaleComp("hptppichi2prim",scale, f1, f2,f2, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "SE-ME" ,"ME","p#pi");

  // projAndScaleComp("hptH3L_l", scale_rt,f1, f3, f3, c,pdf,"l", "p",1,"SE-RT","RT", "dp#pi");
  // projAndScaleComp("hptppimass", scale_rt, f1, f3, f3, c,pdf,"M(p#pi)", "p",1,"SE-RT","RT", "p#pi");
  // projAndScaleComp("hptH3L_ldl", scale_rt,f1, f3, f3, c,pdf,"l/dl", "p",1,"SE-RT","RT", "dp#pi");
  // projAndScaleComp("hptH3L_chi2ndf", scale_rt,f1, f3,f3, c,pdf,"chi2NDF", "p",1,"SE-RT","RT", "dp#pi");
  // projAndScaleComp("hptH3L_chi2topo", scale_rt,f1, f3, f3,c,pdf,"chi2topo", "p",1,"SE-RT","RT", "dp#pi" );
  // projAndScaleComp("hptH3L_dchi2prim", scale_rt,f1, f3, f3,c,pdf,"d chi2primary", "p",1,"SE-RT","RT", "dp#pi");
  // projAndScaleComp("hptH3L_pchi2prim", scale_rt,f1, f3, f3,c,pdf,"p chi2primary", "p",1,"SE-RT","RT", "dp#pi");
  // projAndScaleComp("hptH3L_pichi2prim", scale_rt,f1, f3, f3,c,pdf,"#pi chi2primary", "p",1,"SE-RT","RT", "dp#pi");
  // projAndScaleComp("hptH3L_dDca", scale_rt,f1, f3, f3,c,pdf,"d DCA", "p",1,"SE-RT","RT","dp#pi"  );
  // projAndScaleComp("hptH3L_piDca", scale_rt,f1, f3, f3,c,pdf,"#pi DCA", "p",1,"SE-RT","RT", "dp#pi" );
  // projAndScaleComp("hptH3L_pDca", scale_rt,f1, f3, f3,c,pdf,"#p DCA", "p",1,"SE-RT","RT", "dp#pi" );
  // projAndScaleComp("hptH3L_dpDca", scale_rt,f1, f3, f3,c,pdf,"dp pair DCA", "p",1,"SE-RT","RT", "dp#pi" );
  // projAndScaleComp("hptppichi2prim",scale_rt, f1, f3,f3, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "SE-RT" ,"RT","p#pi");

  cout <<"compare with MC" <<endl;
  // projAndScaleComp("hptH3L_l", scale,f1, f2, fMc, c,pdf,"l", "p",4,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_ppimass", scale, f1, f2, fMc, c,pdf,"M(p#pi)", "p",4,"SE-ME","MC", "p#pi");
  // // projAndScaleCompNoRatio("hptH3L_ppimass", scale, f1, f2, fMc, c,pdf,"M(p#pi)", "p",4,"SE-ME","MC", "p#pi");
  // projAndScaleComp("hptH3L_ldl", scale,f1, f2, fMc, c,pdf,"l/dl", "p",4,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_chi2ndf", scale,f1, f2,fMc, c,pdf,"chi2NDF","p",2,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_chi2topo", scale,f1, f2, fMc,c,pdf,"chi2topo","p",2,"SE-ME","MC", "dp#pi" );
  // // projAndScaleComp("hptH3L_chi2", scale,f1, f2, fMc,c,pdf,"chi2topo","p",2,"SE-ME","MC", "dp#pi" );
  // projAndScaleComp("hptH3L_dchi2prim", scale,f1, f2, fMc,c,pdf,"d chi2primary", "p",5,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_pchi2prim", scale,f1, f2, fMc,c,pdf,"p chi2primary", "p",2,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_pichi2prim", scale,f1, f2, fMc,c,pdf,"#pi chi2primary", "p",2,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_dDca", scale,f1, f2, fMc,c,pdf,"d DCA", "p",2,"SE-ME","MC","dp#pi"  );
  // projAndScaleComp("hptH3L_piDca", scale,f1, f2, fMc,c,pdf,"#pi DCA", "p",2,"SE-ME","MC", "dp#pi" );
  // projAndScaleComp("hptH3L_pDca", scale,f1, f2, fMc,c,pdf,"#p DCA", "p",2,"SE-ME","MC", "dp#pi" );
  // projAndScaleComp("hptH3L_dpDca", scale,f1, f2, fMc,c,pdf,"dp pair DCA", "p",2,"SE-ME","MC", "dp#pi" );
  // // projAndScaleComp("hptppichi2prim",scale, f1, f2,fMc, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",2 , "SE-ME" ,"MC","p#pi");
  //

  // TFile* fMc_ph = TFile::Open("fout_H3L_phase_MC.root");
  // //  projAndScaleComp("hptH3L_l", scale,f1, f2, fMc_ph, c,pdf,"l", "p",4,"SE-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptppimass", scale, f1, f2, fMc_ph, c,pdf,"M(p#pi)", "p",4,"SE-ME","phase_MC", "p#pi");
  // projAndScaleComp("hptH3L_ldl", scale,f1, f2, fMc_ph, c,pdf,"l/dl", "p",4,"SE-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptH3L_chi2ndf", scale,f1, f2,fMc_ph, c,pdf,"chi2NDF","p",2,"SE-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptH3L_chi2topo", scale,f1, f2, fMc_ph,c,pdf,"chi2topo","p",2,"SE-ME","phase_MC", "dp#pi" );
  // // projAndScaleComp("hptH3L_chi2", scale,f1, f2, fMc_ph,c,pdf,"chi2topo","p",2,"SE-ME","phase_MC", "dp#pi" );
  // projAndScaleComp("hptH3L_dchi2prim", scale,f1, f2, fMc_ph,c,pdf,"d chi2primary", "p",2,"SE-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptH3L_pchi2prim", scale,f1, f2, fMc_ph,c,pdf,"p chi2primary", "p",2,"SE-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptH3L_pichi2prim", scale,f1, f2, fMc_ph,c,pdf,"#pi chi2primary", "p",2,"SE-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptH3L_dDca", scale,f1, f2, fMc_ph,c,pdf,"d DCA", "p",2,"SE-ME","phase_MC","dp#pi"  );
  // projAndScaleComp("hptH3L_piDca", scale,f1, f2, fMc_ph,c,pdf,"#pi DCA", "p",2,"SE-ME","phase_MC", "dp#pi" );
  // projAndScaleComp("hptH3L_pDca", scale,f1, f2, fMc_ph,c,pdf,"#p DCA", "p",2,"SE-ME","phase_MC", "dp#pi" );
  // projAndScaleComp("hptH3L_dpDca", scale,f1, f2, fMc_ph,c,pdf,"dp pair DCA", "p",2,"SE-ME","phase_MC", "dp#pi" );
  // projAndScaleComp("hptppichi2prim",scale, f1, f2,fMc_ph, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",2 , "SE-ME" ,"phase_MC","p#pi");

  cout <<"compare with Lambda MC" <<endl;
  projAndScaleComp("hptH3L_l", scale,f1, f2, fMc_ld, c,pdf,"l", "p",4,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_l", scale,f1, f2, fMc, c,pdf,"l", "p",4,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_ppimass", scale, f1, f2, fMc_ld, c,pdf,"M(p#pi)", "p",4,"SE-ME","MC_#Lambda", "p#pi",0,1);
  // projAndScaleCompNoRatio("hptppimass", scale, f1, f2, fMc_ld, c,pdf,"M(p#pi)", "p",4,"SE-ME","MC_#Lambda", "p#pi",0);
  // projAndScaleCompNoRatio("hptH3L_ppimass", scale, f1, f2, fMc_ld, c,pdf,"M(p#pi)", "p",4,"SE-ME","MC_#Lambda", "p#pi",0);
  projAndScaleComp("hptH3L_ldl", scale,f1, f2, fMc_ld, c,pdf,"l/dl", "p",4,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_ldl", scale,f1, f2, fMc, c,pdf,"l/dl", "p",4,"SE-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_chi2ndf", scale,f1, f2,fMc_ld, c,pdf,"chi2NDF","p",2,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_chi2ndf", scale,f1, f2,fMc, c,pdf,"chi2NDF","p",2,"SE-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_chi2topo", scale,f1, f2, fMc_ld,c,pdf,"chi2topo","p",2,"SE-ME","MC_#Lambda", "dp#pi" ,0,1);
  projAndScaleComp("hptH3L_chi2topo", scale,f1, f2, fMc,c,pdf,"chi2topo","p",2,"SE-ME","MC", "dp#pi" );
  // projAndScaleComp("hptH3L_chi2", scale,f1, f2, fMc_ld,c,pdf,"chi2topo","p",2,"SE-ME","MC_#Lambda", "dp#pi" ,0,1);
  projAndScaleComp("hptH3L_dchi2prim", scale,f1, f2, fMc_ld,c,pdf,"d chi2primary", "p",5,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_dchi2prim", scale,f1, f2, fMc,c,pdf,"d chi2primary", "p",5,"SE-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_pchi2prim", scale,f1, f2, fMc_ld,c,pdf,"p chi2primary", "p",2,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_pchi2prim", scale,f1, f2, fMc,c,pdf,"p chi2primary", "p",2,"SE-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_pichi2prim", scale,f1, f2, fMc_ld,c,pdf,"#pi chi2primary", "p",2,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_pichi2prim", scale,f1, f2, fMc,c,pdf,"#pi chi2primary", "p",2,"SE-ME","MC", "dp#pi");

  projAndScaleComp("hptH3L_dDca", scale,f1, f2, fMc_ld,c,pdf,"d DCA", "p",2,"SE-ME","MC_#Lambda","dp#pi"  ,0,1);
  projAndScaleComp("hptH3L_dDca", scale,f1, f2, fMc,c,pdf,"d DCA", "p",2,"SE-ME","MC","dp#pi"  );
  projAndScaleComp("hptH3L_piDca", scale,f1, f2, fMc_ld,c,pdf,"#pi DCA", "p",2,"SE-ME","MC_#Lambda", "dp#pi" ,0,1);
  projAndScaleComp("hptH3L_piDca", scale,f1, f2, fMc,c,pdf,"#pi DCA", "p",2,"SE-ME","MC", "dp#pi" );
  projAndScaleComp("hptH3L_pDca", scale,f1, f2, fMc_ld,c,pdf,"#p DCA", "p",2,"SE-ME","MC_#Lambda", "dp#pi" ,0,1);
  projAndScaleComp("hptH3L_pDca", scale,f1, f2, fMc,c,pdf,"#p DCA", "p",2,"SE-ME","MC", "dp#pi" );
  projAndScaleComp("hptH3L_dpDca", scale,f1, f2, fMc_ld,c,pdf,"dp pair DCA", "p",2,"SE-ME","MC_#Lambda", "dp#pi",0,1 );
  // projAndScaleComp("hptppichi2prim",scale, f1, f2,fMc_ld, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",2 , "SE-ME" ,"MC_#Lambda","p#pi", 0, 1);

  // projAndComp2("hptH3L_chi2topoSig", "hptH3L_chi2topo", fMc_ld, fMc_ld, c,pdf,"(p#pid) #chi^{2}_{topo}" ,"",8 , "Sig" ,"Wide","p#pi");
  // projAndComp2("hptH3L_ldlSig", "hptH3L_ldl", fMc_ld, fMc_ld, c,pdf,"ldl" ,"",8 , "Sig" ,"Wide","p#pi");
  // projAndComp2("hptH3L_lSig", "hptH3L_l", fMc_ld, fMc_ld, c,pdf,"ldl" ,"",8 , "Sig" ,"Wide","p#pi");
  // projAndComp2("hptH3L_lSBR", "hptH3L_l", fMc_ld, fMc_ld, c,pdf,"ldl" ,"",8 , "Sig" ,"Wide","p#pi");
  projAndComp2("hptH3L_ldlSig",  "hptH3L_ldl", fMc, fMc_ld, c,pdf,"ldl" ,"",4 , "Sig" ,"Ld","p#pi");
  projAndComp2("hptH3L_lSig",  "hptH3L_l", fMc, fMc_ld, c,pdf,"l" ,"",4 , "Sig" ,"Ld","p#pi");
  projAndComp2("hptH3L_chi2topoSig",  "hptH3L_chi2topo", fMc, fMc_ld, c,pdf,"chi2topo" ,"",4 , "Sig" ,"Ld","p#pi");
  projAndComp2("hptH3L_chi2ndfSig",  "hptH3L_chi2ndf", fMc, fMc_ld, c,pdf,"chi2ndf" ,"",2 , "Sig" ,"Ld","p#pi");
         
  //add 3 D case
  // cout<<"chi2ndf" << calpurity("h3H3L_chi2ndf",0, 3.5, 1.5, 1, -0.4, -0.6, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 1, 0, 10)<<endl;
  // cout<<"chi2topo" << calpurity("h3H3L_chi2topo",0, 3, 1.5, 1., -0.4, -0.6, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 1, 0, 10)<<endl;
  // cout<<"decaylength" << calpurity("h3H3L_l",8, 99, 1.5, 1, -0.4, -0.6,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 10, 0 ,100 )<<endl;

  // // highy = -0.4; lowy = -0.6;
  // double error;
  // cout<<"chi2ndf" << calpurity("h3H3L_chi2ndf",0, 3.5, highpt, lowpt, highy , lowy, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 1, 0, 5, error, 0.25 )<<endl;
  // // // cout<<"chi2ndf" << calpurity("hptH3L_chi2ndf",0, 3.5, 0., 1, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 5, error )<<endl;
  // cout<<"chi2ndf" << calpurity("h3H3L_chi2ndf",0, 3.5, 1.5, 1., highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 1, 0, 5, error , 0.25)<<endl;
  // cout<<"chi2ndf" << calpurity("h3H3L_chi2ndf",0, 3.5, 2, 1.5, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 1, 0, 5, error , 0.25)<<endl;
  // cout<<"chi2ndf" << calpurity("h3H3L_chi2ndf",0, 3.5, 2.5, 2., highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 1, 0, 5, error , 0.25)<<endl;
  // cout<<"chi2topo" << calpurity("h3H3L_chi2topo",0, 3, highpt, lowpt, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 1, 0, 10, error , 0.25)<<endl;
  // cout<<"chi2topo" << calpurity("h3H3L_chi2topo",0, 3, 1.5, 1., highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 1 , 0, 10, error , 0.25)<<endl;
  // cout<<"chi2topo" << calpurity("h3H3L_chi2topo",0, 3, 2, 1.5, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 1 , 0, 10, error , 0.25)<<endl;
  // cout<<"chi2topo" << calpurity("h3H3L_chi2topo",0, 3, 2.5, 2., highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 1 , 0, 10, error, 0.25 )<<endl;
  // // cout<<"chi2topo" << calpurity("hptH3L_chi2topo",0, 3, lowpt, highpt, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 4 , 0, 50, error )<<endl;
  // cout<<"decaylength" << calpurity("h3H3L_l",8, 99,  highpt, lowpt, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 2, 0 ,100 , error , 0.25)<<endl;
  // cout<<"decaylength" << calpurity("h3H3L_l",8, 99,  1.5, 1, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 5, 0 ,100 , error , 0.25)<<endl;
  // cout<<"decaylength" << calpurity("h3H3L_l",8, 99,  2, 1.5, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 5, 0 ,100 , error , 0.25)<<endl;
  // cout<<"decaylength" << calpurity("h3H3L_l",8, 99,  2.5, 2., highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 5, 0 ,100, error , 0.25)<<endl;
  // // cout<<"decaylength" << calpurity("hptH3L_l",8, 100,  0.2, 2,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 4 , 0 ,50 )<<endl;
  // // cout<<"l/dl" << calpurity("h3H3L_ldl",5, 50, highpt,  lowpt, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l/dl", 2, 3, 50 )<<endl;

  //2D
  // cout<<"chi2ndf" << calpurity("hptH3L_chi2ndf",0, 3.5, 1, 2.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 5)<<endl;
  // cout<<"chi2ndf" << calpurity("hptH3L_chi2ndf",0, 3.5, 0., 1, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 5)<<endl;
  // cout<<"chi2ndf" << calpurity("hptH3L_chi2ndf",0, 3.5, 1., 1.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 5)<<endl;
  // cout<<"chi2ndf" << calpurity("hptH3L_chi2ndf",0, 3.5, 1.5, 2, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 5)<<endl;
  // cout<<"chi2ndf" << calpurity("hptH3L_chi2ndf",0, 3.5, 2., 2.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 5)<<endl;
  // cout<<"chi2topo" << calpurity("hptH3L_chi2topo",0, 3, 1, 2.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 4 , 0, 50)<<endl;
  // cout<<"chi2topo" << calpurity("hptH3L_chi2topo",0, 3, 1, 1.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 4 , 0, 50)<<endl;
  // cout<<"chi2topo" << calpurity("hptH3L_chi2topo",0, 3, 1.5, 2., f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 4 , 0, 50)<<endl;
  // cout<<"chi2topo" << calpurity("hptH3L_chi2topo",0, 3, 2., 2.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 4 , 0, 50)<<endl;
  // // cout<<"chi2topo" << calpurity("hptH3L_chi2topo",0, 3, lowpt, highpt, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 4 , 0, 50)<<endl;
  // cout<<"decaylength" << calpurity("hptH3L_l",8, 100,  lowpt, highpt,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 4 , 0 ,50 )<<endl;
  // // cout<<"decaylength" << calpurity("hptH3L_l",8, 100,  0.2, 2,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 4 , 0 ,50 )<<endl;
  // cout<<"l/dl" << calpurity("hptH3L_ldl",5, 50,  lowpt, highpt,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l/dl", 2, 3, 50 )<<endl;
  // cout<<"dchi2prim" << calpurity("hptH3L_dchi2prim",0, 19,  lowpt, highpt,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L d #chi^{2}",0 )<<endl;
  // cout<<"pichi2prim" << calpurity("hptH3L_pichi2prim",0, 19,  lowpt, highpt,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L pi #chi^{2}",0 )<<endl;
  // cout<<"pip" << calpurity("hptH3L_ppil",0, 39,  lowpt, highpt,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L ppi" , 0)<<endl;
  pdf->On();
  pdf->Close();

}
void drawComp()
{
  // double highpt = 2.5, lowpt = 1., lowy=-0.8, highy = -0.2;
  // double highpt = 2.5, lowpt = 1, lowy=-0.9, highy = -0.1;
  double highpt = 2.5, lowpt = 1, lowy=-0.5, highy = 0;
  TString histname="hH3LMassPtY";
  // TString histname="hH3LMassPtY_5_40";
  SetsPhenixStyle();
  TCanvas* c = new TCanvas("c1","c1");
  // TPDF* pdf = new TPDF("MixEventQA_check.pdf");
  // TPDF* pdf = new TPDF("MixEventQA_Jul27.pdf");
  // TPDF* pdf = new TPDF("MixEventQA_beforeDcacut.pdf");
  // TPDF* pdf = new TPDF("MixEventQA_0050_tight.pdf");
  TPDF* pdf = new TPDF("compare_corr.pdf");
  pdf->Off();
  gStyle->SetPalette(1);

  // TFile *f1 = TFile::Open("fout_H3L_data_SE_large.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_large.root"); 

  // TFile *f1 = TFile::Open("out_KF_test/fout_H3L_KF_test.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_KF_Jul22_part.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_KF_Jul25.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_SE_Jul25.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_SE_Jul27.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Aug09.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_SE_Aug16.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_SE_Aug20.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_SE_Aug23.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02_nodca.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02.root"); 
  TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_tight.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02_0050_010pt.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02_nodca_nopt.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_SE_Aug22.root"); 
  // TFile *f1 = TFile::Open("fout_KF_test.root"); 
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
  double nEvents_se = hcent_se->Integral( 4, 9);
  // double nEvents_se = hcent_se->Integral( 1, 9);

  // TFile *f2 = TFile::Open("fout_H3L_data_ME_large.root"); 
  // TFile *f2 = TFile::Open("fout_H3L_ME_Jul22_part.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Aug16.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Aug20.root"); 
  // TFile *f2 = TFile::Open("output2/fout_H3L_data_ME_Aug24.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Sep02.root"); 
  TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_tight.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Sep02_0050_010pt.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Sep02_nodca.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Sep02_nodca_nopt.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Aug23.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Aug22.root"); 
  // TFile *f2 = TFile::Open("fout_ME_test.root"); 
  // TFile *f2 = TFile::Open("fout_H3L_ME_Jul25.root"); 
  // TFile *f2 = TFile::Open("fout_H3L_data_SE.root"); 
  // TH2F* h2bk = (TH2F*)f2->Get("hptH3Lmass")->Clone("hptH3Lmass_ME");
  // TH3F* h2bk = (TH3F*)f2->Get("hH3LMassPtY_5_40")->Clone("hptH3Lmass_bk");
  TH3F* h2bk = (TH3F*)f2->Get(histname.Data())->Clone("hptH3Lmass_bk");
  h2bk->SetDirectory(0);
  TH1F* hbk = (TH1F*)h2bk->ProjectionY("hbk", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2bk->GetZaxis()->FindBin(lowy), h2bk->GetZaxis()->FindBin(highy));
  hbk->SetDirectory(0);

  // TFile* f3 = TFile::Open("fout_H3L_data_RT_Jul22_full.root");
  // TFile* f3 = TFile::Open("fout_H3L_data_RT_Jul25.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_RT_Jul27.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_KF_Aug20.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_SE_Aug23.root");
  // TFile* f3 = TFile::Open("output2/fout_H3L_data_KF_Aug24.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_SE_Sep02_nodca.root");
  TFile* f3 = TFile::Open("rootfile/fout_H3L_data_SE_tight.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_RT_Aug29.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_KF_Aug22.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_KF_Aug16.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_RT_Aug13.root");
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
  double sig_sb =  hsig->Integral(hsig->GetXaxis()->FindBin(3.01),  hsig->GetXaxis()->FindBin(3.02));
  // double rt_sb =  hrt->Integral(hrt->GetXaxis()->FindBin(2.97), hrt->GetXaxis()->FindBin(2.98)) + hrt->Integral(hrt->GetXaxis()->FindBin(3.0),  hrt->GetXaxis()->FindBin(3.02));
  double rt_sb = hrt->Integral(hrt->GetXaxis()->FindBin(3.01),  hrt->GetXaxis()->FindBin(3.02));
  // double rt_sb =  hrt->Integral(hrt->GetXaxis()->FindBin(2.96), hrt->GetXaxis()->FindBin(2.98)) ;
  double scale_rt = sig_sb/rt_sb;
  // sig_sb = hsig->Integral(hsig->GetXaxis()->FindBin(2.97), hsig->GetXaxis()->FindBin(2.98)) + hsig->Integral(hsig->GetXaxis()->FindBin(3.0),  hsig->GetXaxis()->FindBin(3.02));
  double bk_sb =   hbk->Integral(hbk->GetXaxis()->FindBin(3.01),  hbk->GetXaxis()->FindBin(3.02));
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
  // hrt->Draw("same");
  TLegend* leg_sig = new TLegend(0.65,0.25,0.88,0.45);
  leg_sig->AddEntry(hbk, "ME", "pl");
  // leg_sig->AddEntry(hrt, "RT", "pl");
  // leg_sig->AddEntry(hrt, "KF", "pl");
  leg_sig->AddEntry(hsig, "SE", "pl");
  leg_sig->Draw();

  drawLatex( 0.65,0.61,Form("%0.1f<y<%0.1f",lowy, highy ), 0.055);
  drawLatex( 0.65,0.54,Form("%0.1f<p_{T}<%0.1f GeV/c^{2}", lowpt, highpt), 0.055);
  // drawLatex( 0.65,0.47,Form("5-40%s", "%"), 0.055);
  drawLatex( 0.65,0.47,Form("0-50%s", "%"), 0.055);
  drawBox( 2.97, hsig->GetMinimum(),2.98, hsig->GetMaximum()*0.5 , kBlue-9, 1001 , 1,0.3); 
  // drawBox( 2.989, hsig_bk->GetMinimum(),2.995, hsig_bk->GetMaximum()*0.5 , kOrange-4, 1001,1,0.3); 
  drawBox( 3.01, hsig->GetMinimum(),3.02, hsig->GetMaximum()*0.5 , kBlue-9, 1001,1,0.3 ); 

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
  // hsig_rt->Draw("same");
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
  // hsig_rt->Fit(resfit_rt,"R");
  fit_rt->SetParameter(3, resfit_rt->GetParameter(0) );
  fit_rt->SetParameter(4, resfit_rt->GetParameter(1) );
  // hsig_rt->Fit(fit_rt,"R");
  resfit_rt->SetParameter(0, fit_rt->GetParameter(3));
  resfit_rt->SetParameter(1, fit_rt->GetParameter(4));
  // resfit_rt->Draw("same");

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
  double significance = yield_me/sqrt(sp_counts);
  // double significance = yield_rt/sqrt(sp_counts);
  double s_me = yield_me/(fit->GetParError(0)/hsig_bk->GetBinWidth(1));
  double s_rt = yield_rt/(fit_rt->GetParError(0)/hsig_rt->GetBinWidth(1));
  cout<<"significance: " <<significance << endl;

  TLegend* leg = new TLegend( 0.72, 0.68 ,0.9,0.9 );
  // leg->AddEntry(hbk, "mix-event(ME)","pl");
  // leg->AddEntry(hrt, "rotate d(RT) (scale)","pl");
  leg->AddEntry(hsig, "SE","pl");
  leg->AddEntry(hsig_bk, "SE-ME","pl");
  // leg->AddEntry(hsig_rt, "SE-KF","pl");
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
  drawBox( 3.01, hsig_bk->GetMinimum(),3.02, hsig_bk->GetMaximum()*0.5 , kBlue-9, 1001,1,0.3 ); 

  addpdf(pdf);

  hsig_bk->Draw();
  // hsig_rt->Draw("same");
  TLegend* leg2 = new TLegend( 0.2 , 0.7 ,0.4,0.9  );
  leg2->AddEntry(hsig_rt, "SE-RT","pl");
  leg2->AddEntry(hsig_bk, "SE-ME","pl");
  leg2->Draw();

  addpdf(pdf);
  //
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

  TH2F* hYield = new TH2F( "hYield", "hYield;y,pt", 5, -1., 0., 6, 0, 3);
  double xbw = hYield->GetXaxis()->GetBinWidth(1);
  double ybw = hYield->GetYaxis()->GetBinWidth(1);
  for (int ix=1; ix<=hYield->GetNbinsX(); ix++) {
    for (int iy=1; iy<=hYield->GetNbinsY();iy++){
      double xlow = hYield->GetXaxis()->GetBinLowEdge(ix); // y
      double ylow = hYield->GetYaxis()->GetBinLowEdge(iy); // pt
      double err;
      hYield->SetBinContent( ix , iy, fityield( ylow, ylow+ybw, xlow, xlow+xbw, err, f1, f2, c, pdf) );
      hYield->SetBinError( ix , iy, err);
    }
  }
  hYield->Draw("colz text");
  hYield->SaveAs("fH3L_yield_0050.root");
  addpdf(pdf);

  TFile* fMc_ld = TFile::Open("fout_Lambda_MC_Cuts_0050_tight.root");
  TFile* fMc = TFile::Open("fout_H3L_MC_0050_tight.root");
  // TFile* fMc = TFile::Open("fout_H3L_MC_0080_010pt.root");
  TH1F* hPhase[3];
  TH1F* hPhaseCor[3];
  TH1F* hPurityCor[3];
  TH1F* hBr[3];
  TH1F* hBrRaw[3];
  TGraphErrors* g[3];
  double edge[3][4]={
    { 1.7, 2.1, 2.7, 3.1},
    { 1.2, 1.6, 2.0, 2.4},
    { 1.2, 1.8, 2.2, 2.6}
  };

  double ybin[4]={ 0, -0.25, -0.50, -0.75};
  hPhase[0] = new TH1F("hPhase0", "hPhase0;pt", 25, 1, 3.5);
  hPhase[0]=(TH1F*)hPhase[0]->Rebin(3, "hPhase0", edge[0]);
  hPhase[1] = new TH1F("hPhase1", "hPhase1;pt", 25, 1, 3.5);
  hPhase[1]=(TH1F*)hPhase[1]->Rebin(3, "hPhase1", edge[1]);
  hPhase[2] = new TH1F("hPhase2", "hPhase2;pt", 25, 1, 3.5);
  hPhase[2]=(TH1F*)hPhase[2]->Rebin(3, "hPhase2", edge[2]);
  TFile* fMcH3L = new TFile("fMC_H3L_0050.root");
  TFile* fRcH3L = new TFile("fout_H3L_MC_0050_015pt.root");
  // TFile* fMcH3L = new TFile("fMC_H3L_0080.root");
  // TFile* fRcH3L = new TFile("fout_H3L_MC_0080_010pt.root");
  TH3F* h3Mc = (TH3F*)fMcH3L->Get("hH3LMassPtY")->Clone("h3Mc");
  h3Mc->SetDirectory(0);
  TH3F* h3Rc = (TH3F*)fRcH3L->Get("hH3LMassPtY")->Clone("h3Rc");
  h3Rc->SetDirectory(0);

  h3Mc->Sumw2();
  h3Rc->Sumw2();
  TH2F* h2MC = (TH2F*)h3Mc->Project3D("xz");
  TH2F* h2Rc = (TH2F*)h3Rc->Project3D("xz");

  TH2F* h2Eff = (TH2F*)h2Rc->Clone("h2Eff");
  TH2F* h2temp = (TH2F*)h2MC->Clone("h2temp");
  h2temp->RebinY(10);
  h2temp->RebinX(10);
  h2Eff->RebinY(10);
  h2Eff->RebinX(10);
  h2Eff->Divide(h2temp);
  h2Eff->Draw("col text");
  h2Eff->GetYaxis()->SetRangeUser(0,4.5);
  addpdf(pdf);

  TH1F* heff[3];
  for (int ij=0;ij<3;ij++){ 
    TH1F* h1Mc = (TH1F*)h2MC->ProjectionY(Form("h1Mc%d", ij), h2MC->GetXaxis()->FindBin(ybin[ij+1] +1e-6), h2MC->GetXaxis()->FindBin(ybin[ij]-1e-6) );
    h1Mc = (TH1F*)h1Mc->Rebin( 3, Form("hMc%d",ij), edge[ij]);
    TH1F* h1Rc = (TH1F*)h2Rc->ProjectionY(Form("h1Rc%d", ij), h2Rc->GetXaxis()->FindBin(ybin[ij+1]+1e-6), h2Rc->GetXaxis()->FindBin(ybin[ij]-1e-6) );
    heff[ij] = (TH1F*)h1Rc->Rebin( 3, Form("heff%d",ij), edge[ij]);
    heff[ij]->Divide(h1Mc);
    heff[ij]->GetYaxis()->SetTitle("Eff.");
    heff[ij]->Draw();
    addpdf(pdf);
  }

  double ptesterr[7];
  double chi2topo[7]={0.2,0.6,1,1.5,2,2.5,3 };
  double puritytest[7];
  for (int i=0;i<7;i++) {
    puritytest[i] = calpurity("h3H3L_chi2topo",0,chi2topo[i],  3., 1.5, -0.25, -0.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, ptesterr[i], 0.25);
    // puritytest[i] = calpurity("h3H3L_chi2topo",0, 3,  3., 1., -0., -0.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, ptesterr[i], 0.25);
  }
  for (int i=0;i<7;i++) 
    cout <<chi2topo[i]<<" "<<puritytest[i]<<" "<<ptesterr[i] << endl;
  TGraphErrors* gpurity = new TGraphErrors(7, chi2topo , puritytest,0, ptesterr  );
  gpurity->Draw("pA");
  gpurity->GetXaxis()->SetTitle("Chi2Topo");
  gpurity->GetYaxis()->SetTitle("Purity");
  addpdf(pdf);

  double dy=0.25;
  double dpt,pt;
  
  TFile* f2b = TFile::Open("h_h3l_corr_yield.root_cut00002_cent0_yuehang.root");

  double ny = 3;
  for (int ij=0; ij<ny;ij++){
    double npt =3;
    if (ij==0) npt =2; 
    for (int ipt=0; ipt<npt;ipt++){
      double err;
      dpt = edge[ij][ipt+1]-edge[ij][ipt];
      pt = 0.5*(edge[ij][ipt+1]+edge[ij][ipt]);
      double yield = fityield( edge[ij][ipt], edge[ij][ipt+1], ybin[ij+1], ybin[ij], err, f1, f2, c, pdf);
      // cout << edge[ij][ipt]<<" "<<edge[ij][ipt+1]<<" "<<ybin[ij+1]<<" "<<ybin[ij]<<" "<< yield<<" "<<err <<endl;
      if (ij==0 && ipt==2) { yield=0; err=0;}
      hPhase[ij]->SetBinContent( ipt+1, yield/dy/dpt/pt/2./3.1415926/nEvents_se);
      hPhase[ij]->SetBinError( ipt+1, err/dy/dpt/pt/2./3.1415926/nEvents_se );
    }
    hPhase[ij]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPhase[ij]->GetYaxis()->SetTitle("Raw d^{2}N/(N_{evt}2#pip_{T}dydp_{T})");
    hPhase[ij]->Draw();
    hPhase[ij]->SetDirectory(0);
    addpdf(pdf); 

    hPhaseCor[ij]=(TH1F*)hPhase[ij]->Clone(Form("hYieldCor_%d",ij));
    hPhaseCor[ij]->SetDirectory(0);
    hPurityCor[ij]=(TH1F*)hPhase[ij]->Clone(Form("hPurityCor_%d",ij));
    hPurityCor[ij]->SetDirectory(0);

    double perr;
    double purity1_s = calpurity("h3H3L_chi2ndf",0, 3.5, edge[ij][3], edge[ij][0], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 10, perr, 0.25);
    // double purity2_s = calpurity("h3H3L_chi2topo",0, 3.,  edge[ij][3],  edge[ij][0], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
    double purity2_s = calpurity("h3H3L_chi2topo",0, 3.,  3.,  1.5, ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
    // double puritytest = calpurity("h3H3L_chi2topo",0, 3.,  3.,  1.,ybin[0], ybin[2], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
    // double purity3_s = calpurity("h3H3L_l",8, 99,  3.,  1., ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L decaylength", 2, 0, 100, perr, 0.3);

    for (int ipt=0;ipt<npt;ipt++)
    {
      double y3b = hPhaseCor[ij]->GetBinContent(ipt+1);
      double y3berr = hPhaseCor[ij]->GetBinError(ipt+1)/y3b;
      double eff = heff[ij]->GetBinContent(ipt+1);
      double efferr = heff[ij]->GetBinError(ipt+1)/eff;
      // double purity1 = calpurity("h3H3L_chi2ndf",0, 3.5, edge[ij][ipt+1], edge[ij][ipt], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 10, perr, 0.25);
      // double purity2 = calpurity("h3H3L_chi2topo",0, 3., edge[ij][ipt+1], edge[ij][ipt], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
      // double purity = (purity1+purity2)*0.5;
      // perr = fabs(purity2-purity1);
      // perr = perr/purity;
      double purity = purity2_s;
      // double purity;
      // purity=1;
      // perr = 0;
      double yield_cor = y3b/eff;
      double err = sqrt(y3berr*y3berr + efferr*efferr)*yield_cor;
      hPhaseCor[ij]->SetBinContent(ipt+1,yield_cor );
      hPhaseCor[ij]->SetBinError(ipt+1, err);
      yield_cor = y3b/eff*purity;
      err = sqrt(y3berr*y3berr + efferr*efferr + perr*perr )*yield_cor;
      hPurityCor[ij]->SetBinContent(ipt+1,yield_cor );
      hPurityCor[ij]->SetBinError(ipt+1, err);

      cout <<"pt "<<hPhaseCor[ij]->GetBinCenter(ipt+1) << "purity "<<purity <<" eff " <<eff <<" yield " << yield_cor << " err "<< err<< " s="<<yield_cor/err<<endl;
    }
    hPhaseCor[ij]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPhaseCor[ij]->GetYaxis()->SetTitle("d^{2}N/(N_{evt}2#pip_{T}dydp_{T}(no purity corr.))");
    hPurityCor[ij]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPurityCor[ij]->GetYaxis()->SetTitle("d^{2}N/(N_{evt}2#pip_{T}dydp_{T})");

    hPhaseCor[ij]->Draw();
    hPhaseCor[ij]->SetDirectory(0);
    addpdf(pdf);

    hBrRaw[ij]=(TH1F*)hPhaseCor[ij]->Clone(Form("hBrRaw_%d",ij));
    hBrRaw[ij]->SetDirectory(0);
    g[ij] = (TGraphErrors*)f2b->Get(Form("t_sgct1_corr_yield[0][%d]", ij));
    // g[ij]->Draw();
    for (int ipt=0;ipt<npt;ipt++)
    {
      double y3b = hBrRaw[ij]->GetBinContent(ipt+1);
      double y3berr = hBrRaw[ij]->GetBinError(ipt+1);
      double y2b = g[ij]->GetPointY(ipt);
      double y2bpt = g[ij]->GetPointX(ipt);
      double y2berr = g[ij]->GetErrorY(ipt);
      if (ij==1) {y2berr=y2berr*10; y2b=y2b*10.;}
      if (ij==2) {y2berr=y2berr*100; y2b=y2b*100.;}
      // cout <<"y2b significance:"<<y2b/y2berr<<" yield "<<y2b << endl;
      hBrRaw[ij]->SetBinContent(ipt+1, y2b/(y2b+y3b) );
      hBrRaw[ij]->SetBinError(ipt+1, sqrt((y2b*y2b)/pow((y2b+y3b), 4)*(y3berr*y3berr) + (y3b*y3b)/pow((y2b+y3b), 4)*(y2berr*y2berr)));
      if (y2b/y2berr<2) {
         hBrRaw[ij]->SetBinContent(ipt+1, 0);
         hBrRaw[ij]->SetBinError(ipt+1, 0);
      } 
      // cout <<"bincenter: " << hBr[ij]->GetBinCenter(ipt+1) << " content:"<< hBr[ij]->GetBinContent(ipt+1)<<" err:"<<hBr[ij]->GetBinError(ipt+1)<< "y2b significance:"<<y2b/y2berr<<" yield "<<y2b << " y2bpt "<< y2bpt<< " y3b yield "<< y3b<<endl;
    } 
    hBrRaw[ij]->Draw();
    hBrRaw[ij]->GetYaxis()->SetTitle("R_{3} Raw");
    drawLatex( 0.65,0.88 , Form("%0.2f<y<%0.2f", ybin[ij+1],ybin[ij]) , 0.055);
    addpdf(pdf);

    hBr[ij]=(TH1F*)hPurityCor[ij]->Clone(Form("hBr_%d",ij));
    hBr[ij]->SetDirectory(0);
    g[ij] = (TGraphErrors*)f2b->Get(Form("t_sgct1_corr_yield[0][%d]", ij));
    g[ij]->Draw();
    for (int ipt=0;ipt<npt;ipt++)
    {
      double y3b = hBr[ij]->GetBinContent(ipt+1);
      double y3berr = hBr[ij]->GetBinError(ipt+1);
      double y2b = g[ij]->GetPointY(ipt);
      double y2bpt = g[ij]->GetPointX(ipt);
      double y2berr = g[ij]->GetErrorY(ipt);
      if (ij==1) {y2berr=y2berr*10; y2b=y2b*10.;}
      if (ij==2) {y2berr=y2berr*100; y2b=y2b*100.;}
      // cout <<"y2b significance:"<<y2b/y2berr<<" yield "<<y2b << endl;
      hBr[ij]->SetBinContent(ipt+1, y2b/(y2b+y3b) );
      hBr[ij]->SetBinError(ipt+1, sqrt((y2b*y2b)/pow((y2b+y3b), 4)*(y3berr*y3berr) + (y3b*y3b)/pow((y2b+y3b), 4)*(y2berr*y2berr)));
      if (y2b/y2berr<2) {
         hBr[ij]->SetBinContent(ipt+1, 0);
         hBr[ij]->SetBinError(ipt+1, 0);
      } 
      cout <<"bincenter: " << hBr[ij]->GetBinCenter(ipt+1) << " content:"<< hBr[ij]->GetBinContent(ipt+1)<<" err:"<<hBr[ij]->GetBinError(ipt+1)<< "y2b significance:"<<y2b/y2berr<<" yield "<<y2b << " y2bpt "<< y2bpt<< " y3b yield "<< y3b<<endl;
    } 
    hBr[ij]->Draw();
    hBr[ij]->GetYaxis()->SetTitle("R_{3}");
    drawLatex( 0.65,0.88 , Form("%0.2f<y<%0.2f", ybin[ij+1],ybin[ij]) , 0.055);
    addpdf(pdf);
  }

  TLegend* legb = new TLegend(0.7, 0.7, 0.9, 0.9);
  TH1F* htemp = new TH1F("htemp","htemp;p_{T} (GeV/c); 2-body/(2-body+3-body)",  5, 1, 3);
  htemp->GetYaxis()->SetRangeUser( 0, 0.5);
  htemp->Draw();
  int color[]={kRed, kGreen, kBlue };
  for (int i=0;i<3;i++)
  {
    hBr[i]->Draw("same");
    hBr[i]->SetLineColor(color[i]);
    hBr[i]->SetMarkerColor(color[i]);
    legb->AddEntry( hBr[i], Form("%0.2f<y<%0.2f", ybin[i+1], ybin[i]), "pl");
  }
  legb->Draw();
  addpdf(pdf);

  double tmp=0, tmperr = 0;
  double sum=0, wtsum=0, sumerr=0;
  for (int i=0;i<2;i++)
  {
    for (int ip=1;ip<=hBr[i]->GetNbinsX();ip++)
    {
      tmp = hBr[i]->GetBinContent(ip);
      if (tmp<1e-31) continue; 
      tmperr = hBr[i]->GetBinError(ip);
      cout <<tmperr <<" "<<tmp << endl;
      sum += tmp/tmperr/tmperr;
      sumerr += 1./(tmperr*tmperr);
      wtsum++;
    }
  } 
  cout<<"after correct:" << sum/sumerr<<" "<<sqrt(1./sumerr)<<endl;

  tmp=0; tmperr = 0;
  sum=0; wtsum=0; sumerr=0;
  for (int i=0;i<2;i++)
  {
    for (int ip=1;ip<=hBr[i]->GetNbinsX();ip++)
    {
      tmp = hBrRaw[i]->GetBinContent(ip);
      if (tmp<1e-31) continue; 
      tmperr = hBrRaw[i]->GetBinError(ip);
      cout <<tmperr <<" "<<tmp << endl;
      sum += tmp/tmperr/tmperr;
      sumerr += 1./(tmperr*tmperr);
      wtsum++;
    }
  } 
  cout<<"Before correct:" << sum/sumerr<<" "<<sqrt(1./sumerr)<<endl;

  TFile * fout = new TFile("fout_0050.root","recreate");
  fout->cd();
  for (int i=0;i<3;i++){
    hPhase[i]->Write();
    hPhaseCor[i]->Write();
    hPurityCor[i]->Write();
    hBr[i]->Write();
    heff[i]->Write();
  }
  fout->Close();

  // pdf->On();
  // pdf->Close();
  //
  // return;
  c = new TCanvas( "c","c", 800,800);
  projAndScaleComp("hptH3L_ppichi2prim", scale,f1, f2, fMc_ld,c,pdf,"p#pi chi2primary", "p",2,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_ppichi2prim", scale,f1, f2, fMc,c,pdf,"p#pi chi2primary", "p",2,"SE-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_ppimass", scale,f1, f2, fMc,c,pdf,"p#pi mass", "p",2,"SE-ME","MC", "dp#pi");
  return;
  projAndScaleComp("hptH3L_ppimass", scale,f1, f2, fMc_ld,c,pdf,"p#pi mass", "p",2,"SE-ME","MC", "dp#pi");



  projAndComp("hptH3L_lSBL", f1, f2, c,pdf,"l", "p",4,"SE_SBL","ME_SBL", "dp#pi", scale);
  projAndComp("hptH3L_lSBR", f1, f2, c,pdf,"l", "p",4,"SE_SBR","ME_SBR", "dp#pi", scale);
  projAndComp("hptH3L_ldlSBL", f1, f2, c,pdf,"l/dl", "p",4,"SE_SBL","ME_SBL", "dp#pi", scale);
  projAndComp("hptH3L_ldlSBR", f1, f2, c,pdf,"l/dl", "p",4,"SE_SBR","ME_SBR", "dp#pi", scale);
  projAndComp("hptH3L_chi2ndfSBL", f1, f2, c,pdf,"chi2NDF", "p",4,"SE_SBL","ME_SBL", "dp#pi", scale);
  projAndComp("hptH3L_chi2ndfSBR", f1, f2, c,pdf,"chi2NDF", "p",4,"SE_SBR","ME_SBR", "dp#pi", scale);
  projAndComp("hptH3L_chi2topoSBL", f1, f2, c,pdf,"chi2topo", "p",4,"SE_SBL","ME_SBL", "dp#pi" , scale);
  projAndComp("hptH3L_chi2topoSBR", f1, f2, c,pdf,"chi2topo", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  projAndComp("hptH3L_dchi2primSBL", f1, f2, c,pdf,"d chi2primary", "p",4,"SE_SBL","ME_SBL", "dp#pi", scale);
  projAndComp("hptH3L_dchi2primSBR", f1, f2, c,pdf,"d chi2primary", "p",4,"SE_SBR","ME_SBR", "dp#pi", scale);
  projAndComp("hptH3L_pchi2primSBL", f1, f2, c,pdf,"p chi2primary", "p",4,"SE_SBL","ME_SBL", "dp#pi", scale);
  projAndComp("hptH3L_pchi2primSBR", f1, f2, c,pdf,"p chi2primary", "p",4,"SE_SBR","ME_SBR", "dp#pi", scale);
  projAndComp("hptH3L_pichi2primSBL", f1, f2, c,pdf,"#pi chi2primary", "p",4,"SE_SBL","ME_SBL", "dp#pi", scale);
  projAndComp("hptH3L_pichi2primSBR", f1, f2, c,pdf,"#pi chi2primary", "p",4,"SE_SBR","ME_SBR", "dp#pi", scale);
  projAndComp("hptH3L_dDcaSBL", f1, f2, c,pdf,"d DCA", "p",1,"SE_SBL","ME_SBL","dp#pi"  , scale);
  projAndComp("hptH3L_dDcaSBR", f1, f2, c,pdf,"d DCA", "p",1,"SE_SBR","ME_SBR","dp#pi"  , scale);
  projAndComp("hptH3L_piDcaSBL", f1, f2, c,pdf,"#pi DCA", "p",4,"SE_SBL","ME_SBL", "dp#pi" , scale);
  projAndComp("hptH3L_piDcaSBR", f1, f2, c,pdf,"#pi DCA", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  projAndComp("hptH3L_pDcaSBL", f1, f2, c,pdf,"#p DCA", "p",4,"SE_SBL","ME_SBL", "dp#pi" , scale);
  projAndComp("hptH3L_pDcaSBR", f1, f2, c,pdf,"#p DCA", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  projAndComp("hptH3L_dpDcaSBL", f1, f2, c,pdf,"dp pair DCA", "p",4,"SE_SBL","ME_SBL", "dp#pi" , scale);
  projAndComp("hptH3L_dpDcaSBR", f1, f2, c,pdf,"dp pair DCA", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  // projAndComp("hptppimassSBL", f1, f2, c,pdf,"p#pi Mass", "p",4,"SE_SBL","ME_SBL", "dp#pi" );
  // projAndComp("hptH3L_ppimassSBR", f1, f2, c,pdf,"p#pi Mass", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  // projAndComp("hptH3L_ppichi2ndfSBR", f1, f2, c,pdf,"p#pi Mass", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  // projAndComp("hptH3L_ppichi2ndfSBL", f1, f2, c,pdf,"p#pi chi2ndf", "p",4,"SE_SBR","ME_SBR", "dp#pi" , scale);
  // projAndComp("hptH3L_ppichi2primSBL", f1, f2, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "SE" ,"BK","p#pi",scale);

  // projAndComp("hptH3L_lSBR", f1, f3, c,pdf,"l", "p",4,"SE_SBR","RT_SBR", "dp#pi");
  // projAndComp("hptH3L_ldlSBR", f1, f3, c,pdf,"l/dl", "p",4,"SE_SBR","RT_SBR", "dp#pi");
  // projAndComp("hptH3L_chi2ndfSBR", f1, f3, c,pdf,"chi2NDF", "p",4,"SE_SBR","RT_SBR", "dp#pi");
  // projAndComp("hptH3L_chi2topoSBR", f1, f3, c,pdf,"chi2topo", "p",4,"SE_SBR","RT_SBR", "dp#pi" );
  // projAndComp("hptH3L_dchi2primSBR", f1, f3, c,pdf,"d chi2primary", "p",4,"SE_SBR","RT_SBR", "dp#pi");
  // projAndComp("hptH3L_pchi2primSBR", f1, f3, c,pdf,"p chi2primary", "p",4,"SE_SBR","RT_SBR", "dp#pi");
  // projAndComp("hptH3L_pichi2primSBR", f1, f3, c,pdf,"#pi chi2primary", "p",4,"SE_SBR","RT_SBR", "dp#pi");
  // projAndComp("hptH3L_dDcaSBR", f1, f3, c,pdf,"d DCA", "p",1,"SE_SBR","RT_SBR","dp#pi"  );
  // projAndComp("hptH3L_piDcaSBR", f1, f3, c,pdf,"#pi DCA", "p",4,"SE_SBR","RT_SBR", "dp#pi" );
  // projAndComp("hptH3L_pDcaSBR", f1, f3, c,pdf,"#p DCA", "p",4,"SE_SBR","RT_SBR", "dp#pi" );
  // projAndComp("hptH3L_dpDcaSBR", f1, f3, c,pdf,"dp pair DCA", "p",4,"SE_SBR","RT_SBR", "dp#pi" );
  // projAndComp("hptH3L_ppimassSBR", f1, f3, c,pdf,"p#pi Mass (GeV/c^{2})", "plhist",1 ,"SE" ,"BK","p#pi");
  // // projAndComp("hptppichi2ndfSBR", f1, f2, c,pdf,"p#pi #chi^{2}_{NDF}","" ,1 , "SE" ,"BK","p#pi");
  // // projAndComp("hptppichi2primSBR", f1, f2, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "SE" ,"BK","p#pi");
  // // projAndComp("hptppilSBR", f1, f2, c,pdf,"(p#pi) l" ,"",1 , "SE" ,"BK","p#pi");
  // // projAndComp("hptppildlSBR", f1, f2, c,pdf,"(p#pi) l/#deltal" ,"",1 , "SE" ,"BK","p#pi");

  // projAndComp("hptH3L_lSBL", f1, f3, c,pdf,"l", "p",4,"SE_SBL","RT_SBL", "dp#pi");
  // projAndComp("hptH3L_ldlSBL", f1, f3, c,pdf,"l/dl", "p",4,"SE_SBL","RT_SBL", "dp#pi");
  // projAndComp("hptH3L_chi2ndfSBL", f1, f3, c,pdf,"chi2NDF", "p",4,"SE_SBL","RT_SBL", "dp#pi");
  // projAndComp("hptH3L_chi2topoSBL", f1, f3, c,pdf,"chi2topo", "p",4,"SE_SBL","RT_SBL", "dp#pi" );
  // projAndComp("hptH3L_dchi2primSBL", f1, f3, c,pdf,"d chi2primary", "p",4,"SE_SBL","RT_SBL", "dp#pi");
  // projAndComp("hptH3L_pchi2primSBL", f1, f3, c,pdf,"p chi2primary", "p",4,"SE_SBL","RT_SBL", "dp#pi");
  // projAndComp("hptH3L_pichi2primSBL", f1, f3, c,pdf,"#pi chi2primary", "p",4,"SE_SBL","RT_SBL", "dp#pi");
  // projAndComp("hptH3L_dDcaSBL", f1, f3, c,pdf,"d DCA", "p",1,"SE_SBL","RT_SBL","dp#pi"  );
  // projAndComp("hptH3L_piDcaSBL", f1, f3, c,pdf,"#pi DCA", "p",4,"SE_SBL","RT_SBL", "dp#pi" );
  // projAndComp("hptH3L_pDcaSBL", f1, f3, c,pdf,"#p DCA", "p",4,"SE_SBL","RT_SBL", "dp#pi" );
  // projAndComp("hptH3L_dpDcaSBL", f1, f3, c,pdf,"dp pair DCA", "p",4,"SE_SBL","RT_SBL", "dp#pi" );
  // // projAndComp("hptppimassSBL", f1, f3, c,pdf,"p#pi Mass (GeV/c^{2})", "plhist",1 ,"SE" ,"BK","p#pi");
  // projAndComp("hptH3L_ppimassSBL", f1, f3, c,pdf,"p#pi Mass (GeV/c^{2})", "plhist",1 ,"SE" ,"BK","dp#pi");
  // // projAndComp("hptppichi2ndfSBL", f1, f2, c,pdf,"p#pi #chi^{2}_{NDF}","" ,1 , "SE" ,"BK","p#pi");
  // // projAndComp("hptppichi2primSBL", f1, f2, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "SE" ,"BK","p#pi");
  // // projAndComp("hptppilSBL", f1, f2, c,pdf,"(p#pi) l" ,"",1 , "SE" ,"BK","p#pi");
  // // projAndComp("hptppildlSBL", f1, f2, c,pdf,"(p#pi) l/#deltal" ,"",1 , "SE" ,"BK","p#pi");

  // projAndScaleComp("hptH3L_l", scale,f1, f2, f2, c,pdf,"l", "p",1,"SE-ME","ME", "dp#pi");
  // projAndScaleComp("hptH3L_ppimass", scale, f1, f2, f2, c,pdf,"M(p#pi)", "p",1,"SE-ME","ME", "p#pi");
  // // projAndScaleCompNoRatio("hptppimass", scale, f1, f2, f2, c,pdf,"M(p#pi)", "p",1,"SE-ME","ME", "p#pi");
  // projAndScaleComp("hptH3L_ldl", scale,f1, f2, f2, c,pdf,"l/dl", "p",1,"SE-ME","ME", "dp#pi");
  // projAndScaleComp("hptH3L_chi2ndf", scale,f1, f2,f2, c,pdf,"chi2NDF", "p",1,"SE-ME","ME", "dp#pi");
  projAndScaleComp("hptH3L_chi2topo", scale,f1, f2, f2,c,pdf,"chi2topo", "p",1,"SE-ME","ME", "dp#pi" );
  // projAndScaleComp("hptH3L_dchi2prim", scale,f1, f2, f2,c,pdf,"d chi2primary", "p",1,"SE-ME","ME", "dp#pi");
  // projAndScaleComp("hptH3L_pchi2prim", scale,f1, f2, f2,c,pdf,"p chi2primary", "p",1,"SE-ME","ME", "dp#pi");
  // projAndScaleComp("hptH3L_pichi2prim", scale,f1, f2, f2,c,pdf,"#pi chi2primary", "p",1,"SE-ME","ME", "dp#pi");
  // projAndScaleComp("hptH3L_dDca", scale,f1, f2, f2,c,pdf,"d DCA", "p",1,"SE-ME","ME","dp#pi"  );
  // projAndScaleComp("hptH3L_piDca", scale,f1, f2, f2,c,pdf,"#pi DCA", "p",1,"SE-ME","ME", "dp#pi" );
  // projAndScaleComp("hptH3L_pDca", scale,f1, f2, f2,c,pdf,"#p DCA", "p",1,"SE-ME","ME", "dp#pi" );
  // projAndScaleComp("hptH3L_dpDca", scale,f1, f2, f2,c,pdf,"dp pair DCA", "p",1,"SE-ME","ME", "dp#pi" );
  // // projAndScaleComp("hptppichi2prim",scale, f1, f2,f2, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "SE-ME" ,"ME","p#pi");

  // projAndScaleComp("hptH3L_l", scale_rt,f1, f3, f3, c,pdf,"l", "p",1,"SE-RT","RT", "dp#pi");
  // projAndScaleComp("hptppimass", scale_rt, f1, f3, f3, c,pdf,"M(p#pi)", "p",1,"SE-RT","RT", "p#pi");
  // projAndScaleComp("hptH3L_ldl", scale_rt,f1, f3, f3, c,pdf,"l/dl", "p",1,"SE-RT","RT", "dp#pi");
  // projAndScaleComp("hptH3L_chi2ndf", scale_rt,f1, f3,f3, c,pdf,"chi2NDF", "p",1,"SE-RT","RT", "dp#pi");
  // projAndScaleComp("hptH3L_chi2topo", scale_rt,f1, f3, f3,c,pdf,"chi2topo", "p",1,"SE-RT","RT", "dp#pi" );
  // projAndScaleComp("hptH3L_dchi2prim", scale_rt,f1, f3, f3,c,pdf,"d chi2primary", "p",1,"SE-RT","RT", "dp#pi");
  // projAndScaleComp("hptH3L_pchi2prim", scale_rt,f1, f3, f3,c,pdf,"p chi2primary", "p",1,"SE-RT","RT", "dp#pi");
  // projAndScaleComp("hptH3L_pichi2prim", scale_rt,f1, f3, f3,c,pdf,"#pi chi2primary", "p",1,"SE-RT","RT", "dp#pi");
  // projAndScaleComp("hptH3L_dDca", scale_rt,f1, f3, f3,c,pdf,"d DCA", "p",1,"SE-RT","RT","dp#pi"  );
  // projAndScaleComp("hptH3L_piDca", scale_rt,f1, f3, f3,c,pdf,"#pi DCA", "p",1,"SE-RT","RT", "dp#pi" );
  // projAndScaleComp("hptH3L_pDca", scale_rt,f1, f3, f3,c,pdf,"#p DCA", "p",1,"SE-RT","RT", "dp#pi" );
  // projAndScaleComp("hptH3L_dpDca", scale_rt,f1, f3, f3,c,pdf,"dp pair DCA", "p",1,"SE-RT","RT", "dp#pi" );
  // projAndScaleComp("hptppichi2prim",scale_rt, f1, f3,f3, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "SE-RT" ,"RT","p#pi");

  cout <<"compare with MC" <<endl;
  // projAndScaleComp("hptH3L_l", scale,f1, f2, fMc, c,pdf,"l", "p",4,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_ppimass", scale, f1, f2, fMc, c,pdf,"M(p#pi)", "p",4,"SE-ME","MC", "p#pi");
  // // projAndScaleCompNoRatio("hptH3L_ppimass", scale, f1, f2, fMc, c,pdf,"M(p#pi)", "p",4,"SE-ME","MC", "p#pi");
  // projAndScaleComp("hptH3L_ldl", scale,f1, f2, fMc, c,pdf,"l/dl", "p",4,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_chi2ndf", scale,f1, f2,fMc, c,pdf,"chi2NDF","p",2,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_chi2topo", scale,f1, f2, fMc,c,pdf,"chi2topo","p",2,"SE-ME","MC", "dp#pi" );
  // // projAndScaleComp("hptH3L_chi2", scale,f1, f2, fMc,c,pdf,"chi2topo","p",2,"SE-ME","MC", "dp#pi" );
  // projAndScaleComp("hptH3L_dchi2prim", scale,f1, f2, fMc,c,pdf,"d chi2primary", "p",5,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_pchi2prim", scale,f1, f2, fMc,c,pdf,"p chi2primary", "p",2,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_pichi2prim", scale,f1, f2, fMc,c,pdf,"#pi chi2primary", "p",2,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_dDca", scale,f1, f2, fMc,c,pdf,"d DCA", "p",2,"SE-ME","MC","dp#pi"  );
  // projAndScaleComp("hptH3L_piDca", scale,f1, f2, fMc,c,pdf,"#pi DCA", "p",2,"SE-ME","MC", "dp#pi" );
  // projAndScaleComp("hptH3L_pDca", scale,f1, f2, fMc,c,pdf,"#p DCA", "p",2,"SE-ME","MC", "dp#pi" );
  // projAndScaleComp("hptH3L_dpDca", scale,f1, f2, fMc,c,pdf,"dp pair DCA", "p",2,"SE-ME","MC", "dp#pi" );
  // // projAndScaleComp("hptppichi2prim",scale, f1, f2,fMc, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",2 , "SE-ME" ,"MC","p#pi");
  //

  // TFile* fMc_ph = TFile::Open("fout_H3L_phase_MC.root");
  // //  projAndScaleComp("hptH3L_l", scale,f1, f2, fMc_ph, c,pdf,"l", "p",4,"SE-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptppimass", scale, f1, f2, fMc_ph, c,pdf,"M(p#pi)", "p",4,"SE-ME","phase_MC", "p#pi");
  // projAndScaleComp("hptH3L_ldl", scale,f1, f2, fMc_ph, c,pdf,"l/dl", "p",4,"SE-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptH3L_chi2ndf", scale,f1, f2,fMc_ph, c,pdf,"chi2NDF","p",2,"SE-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptH3L_chi2topo", scale,f1, f2, fMc_ph,c,pdf,"chi2topo","p",2,"SE-ME","phase_MC", "dp#pi" );
  // // projAndScaleComp("hptH3L_chi2", scale,f1, f2, fMc_ph,c,pdf,"chi2topo","p",2,"SE-ME","phase_MC", "dp#pi" );
  // projAndScaleComp("hptH3L_dchi2prim", scale,f1, f2, fMc_ph,c,pdf,"d chi2primary", "p",2,"SE-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptH3L_pchi2prim", scale,f1, f2, fMc_ph,c,pdf,"p chi2primary", "p",2,"SE-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptH3L_pichi2prim", scale,f1, f2, fMc_ph,c,pdf,"#pi chi2primary", "p",2,"SE-ME","phase_MC", "dp#pi");
  // projAndScaleComp("hptH3L_dDca", scale,f1, f2, fMc_ph,c,pdf,"d DCA", "p",2,"SE-ME","phase_MC","dp#pi"  );
  // projAndScaleComp("hptH3L_piDca", scale,f1, f2, fMc_ph,c,pdf,"#pi DCA", "p",2,"SE-ME","phase_MC", "dp#pi" );
  // projAndScaleComp("hptH3L_pDca", scale,f1, f2, fMc_ph,c,pdf,"#p DCA", "p",2,"SE-ME","phase_MC", "dp#pi" );
  // projAndScaleComp("hptH3L_dpDca", scale,f1, f2, fMc_ph,c,pdf,"dp pair DCA", "p",2,"SE-ME","phase_MC", "dp#pi" );
  // projAndScaleComp("hptppichi2prim",scale, f1, f2,fMc_ph, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",2 , "SE-ME" ,"phase_MC","p#pi");

  cout <<"compare with Lambda MC" <<endl;
  projAndScaleComp("hptH3L_l", scale,f1, f2, fMc_ld, c,pdf,"l", "p",4,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_l", scale,f1, f2, fMc, c,pdf,"l", "p",4,"SE-ME","MC", "dp#pi");
  // projAndScaleComp("hptH3L_ppimass", scale, f1, f2, fMc_ld, c,pdf,"M(p#pi)", "p",4,"SE-ME","MC_#Lambda", "p#pi",0,1);
  // projAndScaleCompNoRatio("hptppimass", scale, f1, f2, fMc_ld, c,pdf,"M(p#pi)", "p",4,"SE-ME","MC_#Lambda", "p#pi",0);
  // projAndScaleCompNoRatio("hptH3L_ppimass", scale, f1, f2, fMc_ld, c,pdf,"M(p#pi)", "p",4,"SE-ME","MC_#Lambda", "p#pi",0);
  projAndScaleComp("hptH3L_ldl", scale,f1, f2, fMc_ld, c,pdf,"l/dl", "p",4,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_ldl", scale,f1, f2, fMc, c,pdf,"l/dl", "p",4,"SE-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_chi2ndf", scale,f1, f2,fMc_ld, c,pdf,"chi2NDF","p",2,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_chi2ndf", scale,f1, f2,fMc, c,pdf,"chi2NDF","p",2,"SE-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_chi2topo", scale,f1, f2, fMc_ld,c,pdf,"chi2topo","p",2,"SE-ME","MC_#Lambda", "dp#pi" ,0,1);
  projAndScaleComp("hptH3L_chi2topo", scale,f1, f2, fMc,c,pdf,"chi2topo","p",2,"SE-ME","MC", "dp#pi" );
  // projAndScaleComp("hptH3L_chi2", scale,f1, f2, fMc_ld,c,pdf,"chi2topo","p",2,"SE-ME","MC_#Lambda", "dp#pi" ,0,1);
  projAndScaleComp("hptH3L_dchi2prim", scale,f1, f2, fMc_ld,c,pdf,"d chi2primary", "p",5,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_dchi2prim", scale,f1, f2, fMc,c,pdf,"d chi2primary", "p",5,"SE-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_pchi2prim", scale,f1, f2, fMc_ld,c,pdf,"p chi2primary", "p",2,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_pchi2prim", scale,f1, f2, fMc,c,pdf,"p chi2primary", "p",2,"SE-ME","MC", "dp#pi");
  projAndScaleComp("hptH3L_pichi2prim", scale,f1, f2, fMc_ld,c,pdf,"#pi chi2primary", "p",2,"SE-ME","MC_#Lambda", "dp#pi",0,1);
  projAndScaleComp("hptH3L_pichi2prim", scale,f1, f2, fMc,c,pdf,"#pi chi2primary", "p",2,"SE-ME","MC", "dp#pi");

  projAndScaleComp("hptH3L_dDca", scale,f1, f2, fMc_ld,c,pdf,"d DCA", "p",2,"SE-ME","MC_#Lambda","dp#pi"  ,0,1);
  projAndScaleComp("hptH3L_dDca", scale,f1, f2, fMc,c,pdf,"d DCA", "p",2,"SE-ME","MC","dp#pi"  );
  projAndScaleComp("hptH3L_piDca", scale,f1, f2, fMc_ld,c,pdf,"#pi DCA", "p",2,"SE-ME","MC_#Lambda", "dp#pi" ,0,1);
  projAndScaleComp("hptH3L_piDca", scale,f1, f2, fMc,c,pdf,"#pi DCA", "p",2,"SE-ME","MC", "dp#pi" );
  projAndScaleComp("hptH3L_pDca", scale,f1, f2, fMc_ld,c,pdf,"#p DCA", "p",2,"SE-ME","MC_#Lambda", "dp#pi" ,0,1);
  projAndScaleComp("hptH3L_pDca", scale,f1, f2, fMc,c,pdf,"#p DCA", "p",2,"SE-ME","MC", "dp#pi" );
  projAndScaleComp("hptH3L_dpDca", scale,f1, f2, fMc_ld,c,pdf,"dp pair DCA", "p",2,"SE-ME","MC_#Lambda", "dp#pi",0,1 );
  // projAndScaleComp("hptppichi2prim",scale, f1, f2,fMc_ld, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",2 , "SE-ME" ,"MC_#Lambda","p#pi", 0, 1);

  // projAndComp2("hptH3L_chi2topoSig", "hptH3L_chi2topo", fMc_ld, fMc_ld, c,pdf,"(p#pid) #chi^{2}_{topo}" ,"",8 , "Sig" ,"Wide","p#pi");
  // projAndComp2("hptH3L_ldlSig", "hptH3L_ldl", fMc_ld, fMc_ld, c,pdf,"ldl" ,"",8 , "Sig" ,"Wide","p#pi");
  // projAndComp2("hptH3L_lSig", "hptH3L_l", fMc_ld, fMc_ld, c,pdf,"ldl" ,"",8 , "Sig" ,"Wide","p#pi");
  // projAndComp2("hptH3L_lSBR", "hptH3L_l", fMc_ld, fMc_ld, c,pdf,"ldl" ,"",8 , "Sig" ,"Wide","p#pi");
  projAndComp2("hptH3L_ldlSig",  "hptH3L_ldl", fMc, fMc_ld, c,pdf,"ldl" ,"",4 , "Sig" ,"Ld","p#pi");
  projAndComp2("hptH3L_lSig",  "hptH3L_l", fMc, fMc_ld, c,pdf,"l" ,"",4 , "Sig" ,"Ld","p#pi");
  projAndComp2("hptH3L_chi2topoSig",  "hptH3L_chi2topo", fMc, fMc_ld, c,pdf,"chi2topo" ,"",4 , "Sig" ,"Ld","p#pi");
  projAndComp2("hptH3L_chi2ndfSig",  "hptH3L_chi2ndf", fMc, fMc_ld, c,pdf,"chi2ndf" ,"",2 , "Sig" ,"Ld","p#pi");
         
  //add 3 D case
  // cout<<"chi2ndf" << calpurity("h3H3L_chi2ndf",0, 3.5, 1.5, 1, -0.4, -0.6, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 1, 0, 10)<<endl;
  // cout<<"chi2topo" << calpurity("h3H3L_chi2topo",0, 3, 1.5, 1., -0.4, -0.6, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 1, 0, 10)<<endl;
  // cout<<"decaylength" << calpurity("h3H3L_l",8, 99, 1.5, 1, -0.4, -0.6,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 10, 0 ,100 )<<endl;

  // // highy = -0.4; lowy = -0.6;
  // double error;
  // cout<<"chi2ndf" << calpurity("h3H3L_chi2ndf",0, 3.5, highpt, lowpt, highy , lowy, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 1, 0, 5, error, 0.25 )<<endl;
  // // // cout<<"chi2ndf" << calpurity("hptH3L_chi2ndf",0, 3.5, 0., 1, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 5, error )<<endl;
  // cout<<"chi2ndf" << calpurity("h3H3L_chi2ndf",0, 3.5, 1.5, 1., highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 1, 0, 5, error , 0.25)<<endl;
  // cout<<"chi2ndf" << calpurity("h3H3L_chi2ndf",0, 3.5, 2, 1.5, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 1, 0, 5, error , 0.25)<<endl;
  // cout<<"chi2ndf" << calpurity("h3H3L_chi2ndf",0, 3.5, 2.5, 2., highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 1, 0, 5, error , 0.25)<<endl;
  // cout<<"chi2topo" << calpurity("h3H3L_chi2topo",0, 3, highpt, lowpt, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 1, 0, 10, error , 0.25)<<endl;
  // cout<<"chi2topo" << calpurity("h3H3L_chi2topo",0, 3, 1.5, 1., highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 1 , 0, 10, error , 0.25)<<endl;
  // cout<<"chi2topo" << calpurity("h3H3L_chi2topo",0, 3, 2, 1.5, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 1 , 0, 10, error , 0.25)<<endl;
  // cout<<"chi2topo" << calpurity("h3H3L_chi2topo",0, 3, 2.5, 2., highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 1 , 0, 10, error, 0.25 )<<endl;
  // // cout<<"chi2topo" << calpurity("hptH3L_chi2topo",0, 3, lowpt, highpt, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 4 , 0, 50, error )<<endl;
  // cout<<"decaylength" << calpurity("h3H3L_l",8, 99,  highpt, lowpt, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 2, 0 ,100 , error , 0.25)<<endl;
  // cout<<"decaylength" << calpurity("h3H3L_l",8, 99,  1.5, 1, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 5, 0 ,100 , error , 0.25)<<endl;
  // cout<<"decaylength" << calpurity("h3H3L_l",8, 99,  2, 1.5, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 5, 0 ,100 , error , 0.25)<<endl;
  // cout<<"decaylength" << calpurity("h3H3L_l",8, 99,  2.5, 2., highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 5, 0 ,100, error , 0.25)<<endl;
  // // cout<<"decaylength" << calpurity("hptH3L_l",8, 100,  0.2, 2,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 4 , 0 ,50 )<<endl;
  // // cout<<"l/dl" << calpurity("h3H3L_ldl",5, 50, highpt,  lowpt, highy , lowy,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l/dl", 2, 3, 50 )<<endl;

  //2D
  // cout<<"chi2ndf" << calpurity("hptH3L_chi2ndf",0, 3.5, 1, 2.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 5)<<endl;
  // cout<<"chi2ndf" << calpurity("hptH3L_chi2ndf",0, 3.5, 0., 1, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 5)<<endl;
  // cout<<"chi2ndf" << calpurity("hptH3L_chi2ndf",0, 3.5, 1., 1.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 5)<<endl;
  // cout<<"chi2ndf" << calpurity("hptH3L_chi2ndf",0, 3.5, 1.5, 2, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 5)<<endl;
  // cout<<"chi2ndf" << calpurity("hptH3L_chi2ndf",0, 3.5, 2., 2.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 5)<<endl;
  // cout<<"chi2topo" << calpurity("hptH3L_chi2topo",0, 3, 1, 2.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 4 , 0, 50)<<endl;
  // cout<<"chi2topo" << calpurity("hptH3L_chi2topo",0, 3, 1, 1.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 4 , 0, 50)<<endl;
  // cout<<"chi2topo" << calpurity("hptH3L_chi2topo",0, 3, 1.5, 2., f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 4 , 0, 50)<<endl;
  // cout<<"chi2topo" << calpurity("hptH3L_chi2topo",0, 3, 2., 2.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 4 , 0, 50)<<endl;
  // // cout<<"chi2topo" << calpurity("hptH3L_chi2topo",0, 3, lowpt, highpt, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 4 , 0, 50)<<endl;
  // cout<<"decaylength" << calpurity("hptH3L_l",8, 100,  lowpt, highpt,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 4 , 0 ,50 )<<endl;
  // // cout<<"decaylength" << calpurity("hptH3L_l",8, 100,  0.2, 2,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l", 4 , 0 ,50 )<<endl;
  // cout<<"l/dl" << calpurity("hptH3L_ldl",5, 50,  lowpt, highpt,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L l/dl", 2, 3, 50 )<<endl;
  // cout<<"dchi2prim" << calpurity("hptH3L_dchi2prim",0, 19,  lowpt, highpt,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L d #chi^{2}",0 )<<endl;
  // cout<<"pichi2prim" << calpurity("hptH3L_pichi2prim",0, 19,  lowpt, highpt,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L pi #chi^{2}",0 )<<endl;
  // cout<<"pip" << calpurity("hptH3L_ppil",0, 39,  lowpt, highpt,f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L ppi" , 0)<<endl;
  pdf->On();
  pdf->Close();

}
void drawMixDataScanNDF(int icut, TString histname, TString topohistname, double topocut,  double& Br, double & error)
{
  // double highpt = 2.5, lowpt = 1., lowy=-0.8, highy = -0.2;
  // double highpt = 2.5, lowpt = 1, lowy=-0.9, highy = -0.1;
  double highpt = 4, lowpt = 0, lowy=-1.5, highy = 0;
  // TString histname=Form("hH3LMassPtYTopoCut%d",icut);
  cout << histname.Data()<<endl;
  // TString histname="hH3LMassPtY_5_40";
  SetsPhenixStyle();
  TCanvas* c = new TCanvas("c1","c1");
  // TPDF* pdf = new TPDF("MixEventQA_check.pdf");
  // TPDF* pdf = new TPDF("MixEventQA_Jul27.pdf");
  // TPDF* pdf = new TPDF("MixEventQA_beforeDcacut.pdf");
  TPDF* pdf = new TPDF(Form("MixEventQA_0050_ndfscan_%d.pdf",icut ));
  pdf->Off();
  gStyle->SetPalette(1);

  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02_0050_015pt_scan.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02_0050_015pt_sys2.root"); 
  TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02_0050_015pt_sys.root"); 
  TH3F* h2sig = (TH3F*)f1->Get(histname.Data())->Clone("hptH3Lmass_sig");
  h2sig->SetDirectory(0);
  // TH3F* h2sig->Project3D("xz");
  TH1F* hsig = (TH1F*)h2sig->ProjectionY("hsig", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2sig->GetZaxis()->FindBin(lowy), h2sig->GetZaxis()->FindBin(highy));
  hsig->SetDirectory(0);
  TH1F* hcent_se = (TH1F*)f1->Get("hcent")->Clone("hcent_se");
  double nEvents_se = hcent_se->Integral( 4, 9);
  // double nEvents_se = hcent_se->Integral( 1, 9);

  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Sep02_0050_015pt_sys2.root"); 
  TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Sep02_0050_015pt_sys.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Sep02_0050_015pt_scan.root"); 
  TH3F* h2bk = (TH3F*)f2->Get(histname.Data())->Clone("hptH3Lmass_bk");
  h2bk->SetDirectory(0);
  TH1F* hbk = (TH1F*)h2bk->ProjectionY("hbk", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2bk->GetZaxis()->FindBin(lowy), h2bk->GetZaxis()->FindBin(highy));
  hbk->SetDirectory(0);

  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02_0050_015pt_scan.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_SE_Sep02_0050_015pt_sys2.root");
  TFile* f3 = TFile::Open("rootfile/fout_H3L_data_SE_Sep02_0050_015pt_sys.root");
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
  double sig_sb =  hsig->Integral(hsig->GetXaxis()->FindBin(3.01),  hsig->GetXaxis()->FindBin(3.02));
  // double rt_sb =  hrt->Integral(hrt->GetXaxis()->FindBin(2.97), hrt->GetXaxis()->FindBin(2.98)) + hrt->Integral(hrt->GetXaxis()->FindBin(3.0),  hrt->GetXaxis()->FindBin(3.02));
  double rt_sb = hrt->Integral(hrt->GetXaxis()->FindBin(3.01),  hrt->GetXaxis()->FindBin(3.02));
  // double rt_sb =  hrt->Integral(hrt->GetXaxis()->FindBin(2.96), hrt->GetXaxis()->FindBin(2.98)) ;
  double scale_rt = sig_sb/rt_sb;
  // sig_sb = hsig->Integral(hsig->GetXaxis()->FindBin(2.97), hsig->GetXaxis()->FindBin(2.98)) + hsig->Integral(hsig->GetXaxis()->FindBin(3.0),  hsig->GetXaxis()->FindBin(3.02));
  double bk_sb =   hbk->Integral(hbk->GetXaxis()->FindBin(3.01),  hbk->GetXaxis()->FindBin(3.02));
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
  // hrt->Draw("same");
  TLegend* leg_sig = new TLegend(0.65,0.25,0.88,0.45);
  leg_sig->AddEntry(hbk, "ME", "pl");
  // leg_sig->AddEntry(hrt, "RT", "pl");
  // leg_sig->AddEntry(hrt, "KF", "pl");
  leg_sig->AddEntry(hsig, "SE", "pl");
  leg_sig->Draw();


  drawLatex( 0.65,0.61,Form("%0.1f<y<%0.1f",lowy, highy ), 0.055);
  drawLatex( 0.65,0.54,Form("%0.1f<p_{T}<%0.1f GeV/c^{2}", lowpt, highpt), 0.055);
  // drawLatex( 0.65,0.47,Form("5-40%s", "%"), 0.055);
  drawLatex( 0.65,0.47,Form("0-50%s", "%"), 0.055);
  drawBox( 2.97, hsig->GetMinimum(),2.98, hsig->GetMaximum()*0.5 , kBlue-9, 1001 , 1,0.3); 
  // drawBox( 2.989, hsig_bk->GetMinimum(),2.995, hsig_bk->GetMaximum()*0.5 , kOrange-4, 1001,1,0.3); 
  drawBox( 3.01, hsig->GetMinimum(),3.02, hsig->GetMaximum()*0.5 , kBlue-9, 1001,1,0.3 ); 

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
  // hsig_rt->Draw("same");
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
  // hsig_rt->Fit(resfit_rt,"R");
  fit_rt->SetParameter(3, resfit_rt->GetParameter(0) );
  fit_rt->SetParameter(4, resfit_rt->GetParameter(1) );
  // hsig_rt->Fit(fit_rt,"R");
  resfit_rt->SetParameter(0, fit_rt->GetParameter(3));
  resfit_rt->SetParameter(1, fit_rt->GetParameter(4));
  // resfit_rt->Draw("same");

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
  double significance = yield_me/sqrt(sp_counts);
  // double significance = yield_rt/sqrt(sp_counts);
  double s_me = yield_me/(fit->GetParError(0)/hsig_bk->GetBinWidth(1));
  double s_rt = yield_rt/(fit_rt->GetParError(0)/hsig_rt->GetBinWidth(1));
  cout<<"significance: " <<significance << endl;

  TLegend* leg = new TLegend( 0.72, 0.68 ,0.9,0.9 );
  // leg->AddEntry(hbk, "mix-event(ME)","pl");
  // leg->AddEntry(hrt, "rotate d(RT) (scale)","pl");
  leg->AddEntry(hsig, "SE","pl");
  leg->AddEntry(hsig_bk, "SE-ME","pl");
  // leg->AddEntry(hsig_rt, "SE-KF","pl");
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
  drawBox( 3.01, hsig_bk->GetMinimum(),3.02, hsig_bk->GetMaximum()*0.5 , kBlue-9, 1001,1,0.3 ); 

  addpdf(pdf);

  hsig_bk->Draw();
  // hsig_rt->Draw("same");
  TLegend* leg2 = new TLegend( 0.2 , 0.7 ,0.4,0.9  );
  leg2->AddEntry(hsig_rt, "SE-RT","pl");
  leg2->AddEntry(hsig_bk, "SE-ME","pl");
  leg2->Draw();

  addpdf(pdf);

  //
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

  TH2F* hYield = new TH2F( "hYield", "hYield;y,pt", 5, -1., 0., 6, 0, 3);
  double xbw = hYield->GetXaxis()->GetBinWidth(1);
  double ybw = hYield->GetYaxis()->GetBinWidth(1);
  for (int ix=1; ix<=hYield->GetNbinsX(); ix++) {
    for (int iy=1; iy<=hYield->GetNbinsY();iy++){
      double xlow = hYield->GetXaxis()->GetBinLowEdge(ix); // y
      double ylow = hYield->GetYaxis()->GetBinLowEdge(iy); // pt
      double err;
      hYield->SetBinContent( ix , iy, fityield2( histname,   ylow, ylow+ybw, xlow, xlow+xbw, err, f1, f2, c, pdf) );
      hYield->SetBinError( ix , iy, err);
    }
  }
  hYield->Draw("colz text");
  hYield->SaveAs(Form("fH3L_yield_0050_ndf%d.root",icut));
  addpdf(pdf);

  // TFile* fMc_ld = TFile::Open("fout_Lambda_MC_Cuts_0050_015pt_sys_pid3.root");
  TFile* fMc_ld = TFile::Open("fout_Lambda_MC_Cuts_0050_015pt_sys_mixpid3.root");
  TFile* fMc = TFile::Open("fout_H3L_MC_0050_015pt_sys.root");
  // TFile* fMc = TFile::Open("fout_H3L_MC_0080_010pt.root");
  TH1F* hPhase[3];
  TH1F* hPhaseCor[3];
  TH1F* hPurityCor[3];
  TH1F* hBr[3];
  TH1F* hBrRaw[3];
  TGraphErrors* g[3];
  double edge[3][4]={
    { 1.7, 2.1, 2.7, 3.1},
    { 1.2, 1.6, 2.0, 2.4},
    { 1.2, 1.8, 2.2, 2.6}
  };

  double ybin[4]={ 0, -0.25, -0.50, -0.75};
  hPhase[0] = new TH1F("hPhase0", "hPhase0;pt", 25, 1, 3.5);
  hPhase[0]=(TH1F*)hPhase[0]->Rebin(3, "hPhase0", edge[0]);
  hPhase[1] = new TH1F("hPhase1", "hPhase1;pt", 25, 1, 3.5);
  hPhase[1]=(TH1F*)hPhase[1]->Rebin(3, "hPhase1", edge[1]);
  hPhase[2] = new TH1F("hPhase2", "hPhase2;pt", 25, 1, 3.5);
  hPhase[2]=(TH1F*)hPhase[2]->Rebin(3, "hPhase2", edge[2]);
  TFile* fMcH3L = new TFile("fMC_H3L_0050.root");
  TFile* fRcH3L = new TFile("fout_H3L_MC_0050_015pt_sys.root");
  // TFile* fMcH3L = new TFile("fMC_H3L_0080.root");
  // TFile* fRcH3L = new TFile("fout_H3L_MC_0080_010pt.root");
  TH3F* h3Mc = (TH3F*)fMcH3L->Get("hH3LMassPtY")->Clone("h3Mc");
  h3Mc->SetDirectory(0);
  TH3F* h3Rc = (TH3F*)fRcH3L->Get(histname.Data())->Clone("h3Rc");
  h3Rc->SetDirectory(0);

  h3Mc->Sumw2();
  h3Rc->Sumw2();
  TH2F* h2MC = (TH2F*)h3Mc->Project3D("xz");
  TH2F* h2Rc = (TH2F*)h3Rc->Project3D("xz");

  TH2F* h2Eff = (TH2F*)h2Rc->Clone("h2Eff");
  TH2F* h2temp = (TH2F*)h2MC->Clone("h2temp");
  h2temp->RebinY(10);
  h2temp->RebinX(10);
  h2Eff->RebinY(10);
  h2Eff->RebinX(10);
  h2Eff->Divide(h2temp);
  h2Eff->Draw("col text");
  h2Eff->GetYaxis()->SetRangeUser(0,4.5);
  addpdf(pdf);

  TH1F* heff[3];
  for (int ij=0;ij<3;ij++){ 
    TH1F* h1Mc = (TH1F*)h2MC->ProjectionY(Form("h1Mc%d", ij), h2MC->GetXaxis()->FindBin(ybin[ij+1] +1e-6), h2MC->GetXaxis()->FindBin(ybin[ij]-1e-6) );
    h1Mc = (TH1F*)h1Mc->Rebin( 3, Form("hMc%d",ij), edge[ij]);
    TH1F* h1Rc = (TH1F*)h2Rc->ProjectionY(Form("h1Rc%d", ij), h2Rc->GetXaxis()->FindBin(ybin[ij+1]+1e-6), h2Rc->GetXaxis()->FindBin(ybin[ij]-1e-6) );
    heff[ij] = (TH1F*)h1Rc->Rebin( 3, Form("heff%d",ij), edge[ij]);
    heff[ij]->Divide(h1Mc);
    heff[ij]->GetYaxis()->SetTitle("Eff.");
    heff[ij]->Draw();
    addpdf(pdf);
  }
  double ptesterr[10];
  double chi2topo[10]={1,1.5,1.75,2,2.25,2.5,2.75,3,3.25, 3.5 };
  double puritytest[10];
  for (int i=0;i<10;i++) {
    double perr;
  // puritytest[i] = calpurityCut( Form("hH3LMassPtYNDFCut%d", i),  Form("h3H3L_chi2topoNDFCut%d", i), 0, 2, 3., 1., -0.75, 0, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 0, 0, 10, perr, 0.25);
  //   // puritytest[i] = calpurity("h3H3L_chi2topo",0, 3,  3., 1., -0., -0.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, ptesterr[i], 0.25);
    calpurityCut( Form("h3H3L_chi2topoNDFCut%d", i), Form("hH3LMassPtYNDFCut%d",i ), 0, 2., 3, 1.5, 0, -0.75, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
    calpurityCut( "h3H3L_chi2ndf", Form("hH3LMassPtYNDFCut%d", i), 0, chi2topo[i], 3, 1.5, 0, -0.75, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{ndf}", 1, 0, 10, perr, 0.25);
  }
  // for (int i=0;i<7;i++) 
  //   cout <<chi2topo[i]<<" "<<puritytest[i]<<" "<<ptesterr[i] << endl;
  TGraphErrors* gpurity = new TGraphErrors(10, chi2topo , puritytest,0, ptesterr  );
  gpurity->Draw("pA");
  gpurity->GetXaxis()->SetTitle("Chi2Topo");
  gpurity->GetYaxis()->SetTitle("Purity");
  addpdf(pdf);

  double dy=0.25;
  double dpt,pt;
  
  TFile* f2b = TFile::Open("h_h3l_corr_yield.root_cut00002_cent0_yuehang.root");

  double ny = 3;
  for (int ij=0; ij<ny;ij++){
    double npt =3;
    if (ij==0) npt =2; 
    for (int ipt=0; ipt<npt;ipt++){
      double err;
      dpt = edge[ij][ipt+1]-edge[ij][ipt];
      pt = 0.5*(edge[ij][ipt+1]+edge[ij][ipt]);
      double yield = fityield2(histname,  edge[ij][ipt], edge[ij][ipt+1], ybin[ij+1], ybin[ij], err, f1, f2, c, pdf);
      // cout << edge[ij][ipt]<<" "<<edge[ij][ipt+1]<<" "<<ybin[ij+1]<<" "<<ybin[ij]<<" "<< yield<<" "<<err <<endl;
      if (ij==0 && ipt==2) { yield=0; err=0;}
      hPhase[ij]->SetBinContent( ipt+1, yield/dy/dpt/pt/2./3.1415926/nEvents_se);
      hPhase[ij]->SetBinError( ipt+1, err/dy/dpt/pt/2./3.1415926/nEvents_se );
    }
    hPhase[ij]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPhase[ij]->GetYaxis()->SetTitle("Raw d^{2}N/(N_{evt}2#pip_{T}dydp_{T})");
    hPhase[ij]->Draw();
    hPhase[ij]->SetDirectory(0);
    addpdf(pdf); 

    hPhaseCor[ij]=(TH1F*)hPhase[ij]->Clone(Form("hYieldCor_%d",ij));
    hPhaseCor[ij]->SetDirectory(0);
    hPurityCor[ij]=(TH1F*)hPhase[ij]->Clone(Form("hPurityCor_%d",ij));
    hPurityCor[ij]->SetDirectory(0);

    cout <<"ok?" <<endl;
    double perr;
    // double purity1_s = calpurityCut("h3H3L_chi2ndf", "hH3LMassPtY",0, 3.5, edge[ij][3], edge[ij][0], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 10, perr, 0.25);
    // double purity2_s = calpurity("h3H3L_chi2topo",0, 3.,  edge[ij][3],  edge[ij][0], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
    // cout << "ok?"<<endl;
    // double purity2_s = calpurityCut("h3H3L_chi2topo", Form("hH3LMassPtYTopoCut%d", icut) ,0, chi2topo[icut], edge[ij][3], edge[ij][0], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
    // double purity2_s = calpurityCut( topohistname, histname, 0, topocut, edge[ij][3], edge[ij][0], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
    // double purity2_s = calpurityCut("h3H3L_chi2topo", "hH3LMassPtY" ,0, chi2topo[icut], edge[ij][3], edge[ij][0], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
    // double purity2_s = calpurity("h3H3L_chi2topo",0, chi2topo[icut], edge[ij][3], edge[ij][0], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
    // double puritytest = calpurity("h3H3L_chi2topo",0, 3.,  3.,  1.,ybin[0], ybin[2], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
    // double purity3_s = calpurity("h3H3L_l",8, 99,  3.,  1., ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L decaylength", 2, 0, 100, perr, 0.3);
    // double purity2_s = calpurityCut( topohistname, histname, 0, topocut, edge[ij][ipt+1], edge[ij][ipt], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 4, 0, 10, perr, 0.25);
    double purity2_s = calpurityCut( topohistname, histname, 0, topocut, edge[ij][3], edge[ij][0], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);

    for (int ipt=0;ipt<npt;ipt++)
    {
      double y3b = hPhaseCor[ij]->GetBinContent(ipt+1);
      double y3berr = hPhaseCor[ij]->GetBinError(ipt+1)/y3b;
      double eff = heff[ij]->GetBinContent(ipt+1);
      double efferr = heff[ij]->GetBinError(ipt+1)/eff;
      // double purity1 = calpurity("h3H3L_chi2ndf",0, 3.5, edge[ij][ipt+1], edge[ij][ipt], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 10, perr, 0.25);
      // double purity2 = calpurity("h3H3L_chi2topo", 0, chi2topo[icut], edge[ij][ipt+1], edge[ij][ipt], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
      // double purity = (purity1+purity2)*0.5;
      // perr = fabs(purity2-purity1);
      double purity = purity2_s;
      perr = perr/purity;
      // double purity = purity2;
      // double purity = purity1_s;
      // double purity;
      // purity=1;
      // perr = 0;
      double yield_cor = y3b/eff;
      double err = sqrt(y3berr*y3berr + efferr*efferr)*yield_cor;
      hPhaseCor[ij]->SetBinContent(ipt+1,yield_cor );
      hPhaseCor[ij]->SetBinError(ipt+1, err);
      yield_cor = y3b/eff*purity;
      err = sqrt(y3berr*y3berr + efferr*efferr + perr*perr )*yield_cor;
      hPurityCor[ij]->SetBinContent(ipt+1,yield_cor );
      hPurityCor[ij]->SetBinError(ipt+1, err);

      cout <<"pt "<<hPhaseCor[ij]->GetBinCenter(ipt+1) << "purity "<<purity <<" eff " <<eff <<" yield " << yield_cor << " err "<< err<< " s="<<yield_cor/err<<endl;
    }
    hPhaseCor[ij]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPhaseCor[ij]->GetYaxis()->SetTitle("d^{2}N/(N_{evt}2#pip_{T}dydp_{T}(no purity corr.))");
    hPurityCor[ij]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPurityCor[ij]->GetYaxis()->SetTitle("d^{2}N/(N_{evt}2#pip_{T}dydp_{T})");

    hPhaseCor[ij]->Draw();
    hPhaseCor[ij]->SetDirectory(0);
    addpdf(pdf);

    hBrRaw[ij]=(TH1F*)hPhaseCor[ij]->Clone(Form("hBrRaw_%d",ij));
    hBrRaw[ij]->SetDirectory(0);
    g[ij] = (TGraphErrors*)f2b->Get(Form("t_sgct1_corr_yield[0][%d]", ij));
    // g[ij]->Draw();
    for (int ipt=0;ipt<npt;ipt++)
    {
      double y3b = hBrRaw[ij]->GetBinContent(ipt+1);
      double y3berr = hBrRaw[ij]->GetBinError(ipt+1);
      double y2b = g[ij]->GetPointY(ipt);
      double y2bpt = g[ij]->GetPointX(ipt);
      double y2berr = g[ij]->GetErrorY(ipt);
      if (ij==1) {y2berr=y2berr*10; y2b=y2b*10.;}
      if (ij==2) {y2berr=y2berr*100; y2b=y2b*100.;}
      // cout <<"y2b significance:"<<y2b/y2berr<<" yield "<<y2b << endl;
      hBrRaw[ij]->SetBinContent(ipt+1, y2b/(y2b+y3b) );
      hBrRaw[ij]->SetBinError(ipt+1, sqrt((y2b*y2b)/pow((y2b+y3b), 4)*(y3berr*y3berr) + (y3b*y3b)/pow((y2b+y3b), 4)*(y2berr*y2berr)));
      if (y2b/y2berr<2) {
         hBrRaw[ij]->SetBinContent(ipt+1, 0);
         hBrRaw[ij]->SetBinError(ipt+1, 0);
      } 
      // cout <<"bincenter: " << hBr[ij]->GetBinCenter(ipt+1) << " content:"<< hBr[ij]->GetBinContent(ipt+1)<<" err:"<<hBr[ij]->GetBinError(ipt+1)<< "y2b significance:"<<y2b/y2berr<<" yield "<<y2b << " y2bpt "<< y2bpt<< " y3b yield "<< y3b<<endl;
    } 
    hBrRaw[ij]->Draw();
    hBrRaw[ij]->GetYaxis()->SetTitle("R_{3} Raw");
    drawLatex( 0.65,0.88 , Form("%0.2f<y<%0.2f", ybin[ij+1],ybin[ij]) , 0.055);
    addpdf(pdf);

    hBr[ij]=(TH1F*)hPurityCor[ij]->Clone(Form("hBr_%d",ij));
    hBr[ij]->SetDirectory(0);
    g[ij] = (TGraphErrors*)f2b->Get(Form("t_sgct1_corr_yield[0][%d]", ij));
    g[ij]->Draw();
    for (int ipt=0;ipt<npt;ipt++)
    {
      double y3b = hBr[ij]->GetBinContent(ipt+1);
      double y3berr = hBr[ij]->GetBinError(ipt+1);
      double y2b = g[ij]->GetPointY(ipt);
      double y2bpt = g[ij]->GetPointX(ipt);
      double y2berr = g[ij]->GetErrorY(ipt);
      if (ij==1) {y2berr=y2berr*10; y2b=y2b*10.;}
      if (ij==2) {y2berr=y2berr*100; y2b=y2b*100.;}
      // cout <<"y2b significance:"<<y2b/y2berr<<" yield "<<y2b << endl;
      hBr[ij]->SetBinContent(ipt+1, y2b/(y2b+y3b) );
      hBr[ij]->SetBinError(ipt+1, sqrt((y2b*y2b)/pow((y2b+y3b), 4)*(y3berr*y3berr) + (y3b*y3b)/pow((y2b+y3b), 4)*(y2berr*y2berr)));
      if (y2b/y2berr<2) {
         hBr[ij]->SetBinContent(ipt+1, 0);
         hBr[ij]->SetBinError(ipt+1, 0);
      } 
      cout <<"bincenter: " << hBr[ij]->GetBinCenter(ipt+1) << " content:"<< hBr[ij]->GetBinContent(ipt+1)<<" err:"<<hBr[ij]->GetBinError(ipt+1)<< "y2b significance:"<<y2b/y2berr<<" yield "<<y2b << " y2bpt "<< y2bpt<< " y3b yield "<< y3b<<endl;
    } 
    hBr[ij]->Draw();
    hBr[ij]->GetYaxis()->SetTitle("R_{3}");
    drawLatex( 0.65,0.88 , Form("%0.2f<y<%0.2f", ybin[ij+1],ybin[ij]) , 0.055);
    addpdf(pdf);
  }

  TLegend* legb = new TLegend(0.7, 0.7, 0.9, 0.9);
  TH1F* htemp = new TH1F("htemp","htemp;p_{T} (GeV/c); 2-body/(2-body+3-body)",  5, 1, 3);
  htemp->GetYaxis()->SetRangeUser( 0, 0.5);
  htemp->Draw();
  int color[]={kRed-4, kGreen+2, kBlue-4 };
  int style[]={kOpenCircle, kOpenSquare, kOpenStar };
  for (int i=0;i<3;i++)
  {
    hBr[i]->Draw("same");
    hBr[i]->SetLineColor(color[i]);
    hBr[i]->SetMarkerColor(color[i]);
    legb->AddEntry( hBr[i], Form("%0.2f<y<%0.2f", ybin[i+1], ybin[i]), "pl");
  }
  legb->Draw();
  addpdf(pdf);

  double tmp=0, tmperr = 0;
  double sum=0, wtsum=0, sumerr=0;
  for (int i=0;i<2;i++)
  {
    for (int ip=1;ip<=hBr[i]->GetNbinsX();ip++)
    {
      tmp = hBr[i]->GetBinContent(ip);
      if (tmp<1e-31) continue; 
      tmperr = hBr[i]->GetBinError(ip);
      cout <<tmperr <<" "<<tmp << endl;
      sum += tmp/tmperr/tmperr;
      sumerr += 1./(tmperr*tmperr);
      wtsum++;
    }
  } 
  Br = sum/sumerr;
  error = sqrt(1./sumerr);
  cout<<"after correct:" << sum/sumerr<<" "<<sqrt(1./sumerr)<<endl;

  tmp=0; tmperr = 0;
  sum=0; wtsum=0; sumerr=0;
  for (int i=0;i<2;i++)
  {
    for (int ip=1;ip<=hBr[i]->GetNbinsX();ip++)
    {
      tmp = hBrRaw[i]->GetBinContent(ip);
      if (tmp<1e-31) continue; 
      tmperr = hBrRaw[i]->GetBinError(ip);
      cout <<tmperr <<" "<<tmp << endl;
      sum += tmp/tmperr/tmperr;
      sumerr += 1./(tmperr*tmperr);
      wtsum++;
    }
  } 
  cout<<"Before correct:" << sum/sumerr<<" "<<sqrt(1./sumerr)<<endl;

  TFile * fout = new TFile(Form("fout_0050_ndfscan_%d.root", icut),"recreate");
  fout->cd();
  for (int i=0;i<3;i++){
    hPhase[i]->Write();
    hPhaseCor[i]->Write();
    hPurityCor[i]->Write();
    hBr[i]->Write();
    heff[i]->Write();
  }
  fout->Close();

  pdf->On();
  pdf->Close();

}
void drawMixDataScanTopo(int icut, TString histname, TString topohistname, double topocut,  double& Br, double & error)
{
  // double highpt = 2.5, lowpt = 1., lowy=-0.8, highy = -0.2;
  // double highpt = 2.5, lowpt = 1, lowy=-0.9, highy = -0.1;
  double highpt = 4, lowpt = 0, lowy=-1.5, highy = 0;
  // TString histname=Form("hH3LMassPtYTopoCut%d",icut);
  cout << histname.Data()<<endl;
  // TString histname="hH3LMassPtY_5_40";
  SetsPhenixStyle();
  TCanvas* c = new TCanvas("c1","c1");
  // TPDF* pdf = new TPDF("MixEventQA_check.pdf");
  // TPDF* pdf = new TPDF("MixEventQA_Jul27.pdf");
  // TPDF* pdf = new TPDF("MixEventQA_beforeDcacut.pdf");
  TPDF* pdf = new TPDF(Form("MixEventQA_0050_scan_%d.pdf",icut ));
  pdf->Off();
  gStyle->SetPalette(1);

  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02_0050_015pt_scan.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02_0050_015pt_sys2.root"); 
  TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02_0050_015pt_sys.root"); 
  TH3F* h2sig = (TH3F*)f1->Get(histname.Data())->Clone("hptH3Lmass_sig");
  h2sig->SetDirectory(0);
  // TH3F* h2sig->Project3D("xz");
  TH1F* hsig = (TH1F*)h2sig->ProjectionY("hsig", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2sig->GetZaxis()->FindBin(lowy), h2sig->GetZaxis()->FindBin(highy));
  hsig->SetDirectory(0);
  TH1F* hcent_se = (TH1F*)f1->Get("hcent")->Clone("hcent_se");
  double nEvents_se = hcent_se->Integral( 4, 9);
  // double nEvents_se = hcent_se->Integral( 1, 9);

  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Sep02_0050_015pt_sys2.root"); 
  TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Sep02_0050_015pt_sys.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Sep02_0050_015pt_scan.root"); 
  TH3F* h2bk = (TH3F*)f2->Get(histname.Data())->Clone("hptH3Lmass_bk");
  h2bk->SetDirectory(0);
  TH1F* hbk = (TH1F*)h2bk->ProjectionY("hbk", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2bk->GetZaxis()->FindBin(lowy), h2bk->GetZaxis()->FindBin(highy));
  hbk->SetDirectory(0);

  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02_0050_015pt_scan.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_SE_Sep02_0050_015pt_sys2.root");
  TFile* f3 = TFile::Open("rootfile/fout_H3L_data_SE_Sep02_0050_015pt_sys.root");
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
  double sig_sb =  hsig->Integral(hsig->GetXaxis()->FindBin(3.01),  hsig->GetXaxis()->FindBin(3.02));
  // double rt_sb =  hrt->Integral(hrt->GetXaxis()->FindBin(2.97), hrt->GetXaxis()->FindBin(2.98)) + hrt->Integral(hrt->GetXaxis()->FindBin(3.0),  hrt->GetXaxis()->FindBin(3.02));
  double rt_sb = hrt->Integral(hrt->GetXaxis()->FindBin(3.01),  hrt->GetXaxis()->FindBin(3.02));
  // double rt_sb =  hrt->Integral(hrt->GetXaxis()->FindBin(2.96), hrt->GetXaxis()->FindBin(2.98)) ;
  double scale_rt = sig_sb/rt_sb;
  // sig_sb = hsig->Integral(hsig->GetXaxis()->FindBin(2.97), hsig->GetXaxis()->FindBin(2.98)) + hsig->Integral(hsig->GetXaxis()->FindBin(3.0),  hsig->GetXaxis()->FindBin(3.02));
  double bk_sb =   hbk->Integral(hbk->GetXaxis()->FindBin(3.01),  hbk->GetXaxis()->FindBin(3.02));
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
  // hrt->Draw("same");
  TLegend* leg_sig = new TLegend(0.65,0.25,0.88,0.45);
  leg_sig->AddEntry(hbk, "ME", "pl");
  // leg_sig->AddEntry(hrt, "RT", "pl");
  // leg_sig->AddEntry(hrt, "KF", "pl");
  leg_sig->AddEntry(hsig, "SE", "pl");
  leg_sig->Draw();


  drawLatex( 0.65,0.61,Form("%0.1f<y<%0.1f",lowy, highy ), 0.055);
  drawLatex( 0.65,0.54,Form("%0.1f<p_{T}<%0.1f GeV/c^{2}", lowpt, highpt), 0.055);
  // drawLatex( 0.65,0.47,Form("5-40%s", "%"), 0.055);
  drawLatex( 0.65,0.47,Form("0-50%s", "%"), 0.055);
  drawBox( 2.97, hsig->GetMinimum(),2.98, hsig->GetMaximum()*0.5 , kBlue-9, 1001 , 1,0.3); 
  // drawBox( 2.989, hsig_bk->GetMinimum(),2.995, hsig_bk->GetMaximum()*0.5 , kOrange-4, 1001,1,0.3); 
  drawBox( 3.01, hsig->GetMinimum(),3.02, hsig->GetMaximum()*0.5 , kBlue-9, 1001,1,0.3 ); 

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
  // hsig_rt->Draw("same");
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
  // hsig_rt->Fit(resfit_rt,"R");
  fit_rt->SetParameter(3, resfit_rt->GetParameter(0) );
  fit_rt->SetParameter(4, resfit_rt->GetParameter(1) );
  // hsig_rt->Fit(fit_rt,"R");
  resfit_rt->SetParameter(0, fit_rt->GetParameter(3));
  resfit_rt->SetParameter(1, fit_rt->GetParameter(4));
  // resfit_rt->Draw("same");

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
  double significance = yield_me/sqrt(sp_counts);
  // double significance = yield_rt/sqrt(sp_counts);
  double s_me = yield_me/(fit->GetParError(0)/hsig_bk->GetBinWidth(1));
  double s_rt = yield_rt/(fit_rt->GetParError(0)/hsig_rt->GetBinWidth(1));
  cout<<"significance: " <<significance << endl;

  TLegend* leg = new TLegend( 0.72, 0.68 ,0.9,0.9 );
  // leg->AddEntry(hbk, "mix-event(ME)","pl");
  // leg->AddEntry(hrt, "rotate d(RT) (scale)","pl");
  leg->AddEntry(hsig, "SE","pl");
  leg->AddEntry(hsig_bk, "SE-ME","pl");
  // leg->AddEntry(hsig_rt, "SE-KF","pl");
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
  drawBox( 3.01, hsig_bk->GetMinimum(),3.02, hsig_bk->GetMaximum()*0.5 , kBlue-9, 1001,1,0.3 ); 

  addpdf(pdf);

  hsig_bk->Draw();
  // hsig_rt->Draw("same");
  TLegend* leg2 = new TLegend( 0.2 , 0.7 ,0.4,0.9  );
  leg2->AddEntry(hsig_rt, "SE-RT","pl");
  leg2->AddEntry(hsig_bk, "SE-ME","pl");
  leg2->Draw();

  addpdf(pdf);

  //
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

  cout << "????????"<< endl;
  TH2F* hYield = new TH2F( "hYield", "hYield;y,pt", 5, -1., 0., 6, 0, 3);
  double xbw = hYield->GetXaxis()->GetBinWidth(1);
  double ybw = hYield->GetYaxis()->GetBinWidth(1);
  for (int ix=1; ix<=hYield->GetNbinsX(); ix++) {
    for (int iy=1; iy<=hYield->GetNbinsY();iy++){
      double xlow = hYield->GetXaxis()->GetBinLowEdge(ix); // y
      double ylow = hYield->GetYaxis()->GetBinLowEdge(iy); // pt
      double err;
      hYield->SetBinContent( ix , iy, fityield2( histname,   ylow, ylow+ybw, xlow, xlow+xbw, err, f1, f2, c, pdf) );
      hYield->SetBinError( ix , iy, err);
    }
  }
  hYield->Draw("colz text");
  hYield->SaveAs(Form("fH3L_yield_0050_%d.root",icut));
  addpdf(pdf);

  // TFile* fMc_ld = TFile::Open("fout_Lambda_MC_Cuts_0050_015pt_sys_pid3.root");
  TFile* fMc_ld = TFile::Open("fout_Lambda_MC_Cuts_0050_015pt_sys_mixpid3.root");
  TFile* fMc = TFile::Open("fout_H3L_MC_0050_015pt_sys.root");
  // TFile* fMc = TFile::Open("fout_H3L_MC_0080_010pt.root");
  TH1F* hPhase[3];
  TH1F* hPhaseCor[3];
  TH1F* hPurityCor[3];
  TH1F* hBr[3];
  TH1F* hBrRaw[3];
  TGraphErrors* g[3];
  double edge[3][4]={
    { 1.7, 2.1, 2.7, 3.1},
    { 1.2, 1.6, 2.0, 2.4},
    { 1.2, 1.8, 2.2, 2.6}
  };

  double ybin[4]={ 0, -0.25, -0.50, -0.75};
  hPhase[0] = new TH1F("hPhase0", "hPhase0;pt", 25, 1, 3.5);
  hPhase[0]=(TH1F*)hPhase[0]->Rebin(3, "hPhase0", edge[0]);
  hPhase[1] = new TH1F("hPhase1", "hPhase1;pt", 25, 1, 3.5);
  hPhase[1]=(TH1F*)hPhase[1]->Rebin(3, "hPhase1", edge[1]);
  hPhase[2] = new TH1F("hPhase2", "hPhase2;pt", 25, 1, 3.5);
  hPhase[2]=(TH1F*)hPhase[2]->Rebin(3, "hPhase2", edge[2]);
  TFile* fMcH3L = new TFile("fMC_H3L_0050.root");
  TFile* fRcH3L = new TFile("fout_H3L_MC_0050_015pt_sys.root");
  // TFile* fMcH3L = new TFile("fMC_H3L_0080.root");
  // TFile* fRcH3L = new TFile("fout_H3L_MC_0080_010pt.root");
  TH3F* h3Mc = (TH3F*)fMcH3L->Get("hH3LMassPtY")->Clone("h3Mc");
  h3Mc->SetDirectory(0);
  TH3F* h3Rc = (TH3F*)fRcH3L->Get(histname.Data())->Clone("h3Rc");
  h3Rc->SetDirectory(0);

  h3Mc->Sumw2();
  h3Rc->Sumw2();
  TH2F* h2MC = (TH2F*)h3Mc->Project3D("xz");
  TH2F* h2Rc = (TH2F*)h3Rc->Project3D("xz");

  TH2F* h2Eff = (TH2F*)h2Rc->Clone("h2Eff");
  TH2F* h2temp = (TH2F*)h2MC->Clone("h2temp");
  h2temp->RebinY(10);
  h2temp->RebinX(10);
  h2Eff->RebinY(10);
  h2Eff->RebinX(10);
  h2Eff->Divide(h2temp);
  h2Eff->Draw("col text");
  h2Eff->GetYaxis()->SetRangeUser(0,4.5);
  addpdf(pdf);

  TH1F* heff[3];
  for (int ij=0;ij<3;ij++){ 
    TH1F* h1Mc = (TH1F*)h2MC->ProjectionY(Form("h1Mc%d", ij), h2MC->GetXaxis()->FindBin(ybin[ij+1] +1e-6), h2MC->GetXaxis()->FindBin(ybin[ij]-1e-6) );
    h1Mc = (TH1F*)h1Mc->Rebin( 3, Form("hMc%d",ij), edge[ij]);
    TH1F* h1Rc = (TH1F*)h2Rc->ProjectionY(Form("h1Rc%d", ij), h2Rc->GetXaxis()->FindBin(ybin[ij+1]+1e-6), h2Rc->GetXaxis()->FindBin(ybin[ij]-1e-6) );
    heff[ij] = (TH1F*)h1Rc->Rebin( 3, Form("heff%d",ij), edge[ij]);
    heff[ij]->Divide(h1Mc);
    heff[ij]->GetYaxis()->SetTitle("Eff.");
    heff[ij]->Draw();
    addpdf(pdf);
  }

  double ptesterr[7];
  double chi2topo[7]={0.2,0.6,1,1.5,2,2.5,3 };
  double puritytest[7];
  for (int i=0;i<7;i++) {
    puritytest[i] = calpurity("h3H3L_chi2topo",0,chi2topo[i],  3., 1., -0.25, -0.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, ptesterr[i], 0.25);
    // puritytest[i] = calpurity("h3H3L_chi2topo",0, 3,  3., 1., -0., -0.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, ptesterr[i], 0.25);
  }
  for (int i=0;i<7;i++) 
    cout <<chi2topo[i]<<" "<<puritytest[i]<<" "<<ptesterr[i] << endl;
  TGraphErrors* gpurity = new TGraphErrors(7, chi2topo , puritytest,0, ptesterr  );
  gpurity->Draw("pA");
  gpurity->GetXaxis()->SetTitle("Chi2Topo");
  gpurity->GetYaxis()->SetTitle("Purity");
  addpdf(pdf);

  double dy=0.25;
  double dpt,pt;
  
  TFile* f2b = TFile::Open("h_h3l_corr_yield.root_cut00002_cent0_yuehang.root");

  double ny = 3;
  for (int ij=0; ij<ny;ij++){
    double npt =3;
    if (ij==0) npt =2; 
    for (int ipt=0; ipt<npt;ipt++){
      double err;
      dpt = edge[ij][ipt+1]-edge[ij][ipt];
      pt = 0.5*(edge[ij][ipt+1]+edge[ij][ipt]);
      double yield = fityield2(histname,  edge[ij][ipt], edge[ij][ipt+1], ybin[ij+1], ybin[ij], err, f1, f2, c, pdf);
      // cout << edge[ij][ipt]<<" "<<edge[ij][ipt+1]<<" "<<ybin[ij+1]<<" "<<ybin[ij]<<" "<< yield<<" "<<err <<endl;
      if (ij==0 && ipt==2) { yield=0; err=0;}
      hPhase[ij]->SetBinContent( ipt+1, yield/dy/dpt/pt/2./3.1415926/nEvents_se);
      hPhase[ij]->SetBinError( ipt+1, err/dy/dpt/pt/2./3.1415926/nEvents_se );
    }
    hPhase[ij]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPhase[ij]->GetYaxis()->SetTitle("Raw d^{2}N/(N_{evt}2#pip_{T}dydp_{T})");
    hPhase[ij]->Draw();
    hPhase[ij]->SetDirectory(0);
    addpdf(pdf); 

    hPhaseCor[ij]=(TH1F*)hPhase[ij]->Clone(Form("hYieldCor_%d",ij));
    hPhaseCor[ij]->SetDirectory(0);
    hPurityCor[ij]=(TH1F*)hPhase[ij]->Clone(Form("hPurityCor_%d",ij));
    hPurityCor[ij]->SetDirectory(0);

    cout <<"ok?" <<endl;
    double perr;
    // double purity1_s = calpurityCut("h3H3L_chi2ndf", "hH3LMassPtY",0, 3.5, edge[ij][3], edge[ij][0], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 10, perr, 0.25);
    // double purity2_s = calpurity("h3H3L_chi2topo",0, 3.,  edge[ij][3],  edge[ij][0], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
    // cout << "ok?"<<endl;
    // double purity2_s = calpurityCut("h3H3L_chi2topo", Form("hH3LMassPtYTopoCut%d", icut) ,0, chi2topo[icut], edge[ij][3], edge[ij][0], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
    double purity2_s = calpurityCut( topohistname, histname, 0, topocut, edge[ij][3], edge[ij][0], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
    // double purity2_s = calpurityCut("h3H3L_chi2topo", "hH3LMassPtY" ,0, chi2topo[icut], edge[ij][3], edge[ij][0], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
    // double purity2_s = calpurity("h3H3L_chi2topo",0, chi2topo[icut], edge[ij][3], edge[ij][0], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
    // double puritytest = calpurity("h3H3L_chi2topo",0, 3.,  3.,  1.,ybin[0], ybin[2], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
    // double purity3_s = calpurity("h3H3L_l",8, 99,  3.,  1., ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L decaylength", 2, 0, 100, perr, 0.3);

    for (int ipt=0;ipt<npt;ipt++)
    {
      double y3b = hPhaseCor[ij]->GetBinContent(ipt+1);
      double y3berr = hPhaseCor[ij]->GetBinError(ipt+1)/y3b;
      double eff = heff[ij]->GetBinContent(ipt+1);
      double efferr = heff[ij]->GetBinError(ipt+1)/eff;
      // double purity1 = calpurity("h3H3L_chi2ndf",0, 3.5, edge[ij][ipt+1], edge[ij][ipt], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 10, perr, 0.25);
      // double purity2 = calpurity("h3H3L_chi2topo", 0, chi2topo[icut], edge[ij][ipt+1], edge[ij][ipt], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
      // double purity = (purity1+purity2)*0.5;
      // perr = fabs(purity2-purity1);
      double purity = purity2_s;
      perr = perr/purity;
      // double purity = purity2;
      // double purity = purity1_s;
      // double purity;
      // purity=1;
      // perr = 0;
      double yield_cor = y3b/eff;
      double err = sqrt(y3berr*y3berr + efferr*efferr)*yield_cor;
      hPhaseCor[ij]->SetBinContent(ipt+1,yield_cor );
      hPhaseCor[ij]->SetBinError(ipt+1, err);
      yield_cor = y3b/eff*purity;
      err = sqrt(y3berr*y3berr + efferr*efferr + perr*perr )*yield_cor;
      hPurityCor[ij]->SetBinContent(ipt+1,yield_cor );
      hPurityCor[ij]->SetBinError(ipt+1, err);

      cout <<"pt "<<hPhaseCor[ij]->GetBinCenter(ipt+1) << "purity "<<purity <<" eff " <<eff <<" yield " << yield_cor << " err "<< err<< " s="<<yield_cor/err<<endl;
    }
    hPhaseCor[ij]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPhaseCor[ij]->GetYaxis()->SetTitle("d^{2}N/(N_{evt}2#pip_{T}dydp_{T}(no purity corr.))");
    hPurityCor[ij]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPurityCor[ij]->GetYaxis()->SetTitle("d^{2}N/(N_{evt}2#pip_{T}dydp_{T})");

    hPhaseCor[ij]->Draw();
    hPhaseCor[ij]->SetDirectory(0);
    addpdf(pdf);

    hBrRaw[ij]=(TH1F*)hPhaseCor[ij]->Clone(Form("hBrRaw_%d",ij));
    hBrRaw[ij]->SetDirectory(0);
    g[ij] = (TGraphErrors*)f2b->Get(Form("t_sgct1_corr_yield[0][%d]", ij));
    // g[ij]->Draw();
    for (int ipt=0;ipt<npt;ipt++)
    {
      double y3b = hBrRaw[ij]->GetBinContent(ipt+1);
      double y3berr = hBrRaw[ij]->GetBinError(ipt+1);
      double y2b = g[ij]->GetPointY(ipt);
      double y2bpt = g[ij]->GetPointX(ipt);
      double y2berr = g[ij]->GetErrorY(ipt);
      if (ij==1) {y2berr=y2berr*10; y2b=y2b*10.;}
      if (ij==2) {y2berr=y2berr*100; y2b=y2b*100.;}
      // cout <<"y2b significance:"<<y2b/y2berr<<" yield "<<y2b << endl;
      hBrRaw[ij]->SetBinContent(ipt+1, y2b/(y2b+y3b) );
      hBrRaw[ij]->SetBinError(ipt+1, sqrt((y2b*y2b)/pow((y2b+y3b), 4)*(y3berr*y3berr) + (y3b*y3b)/pow((y2b+y3b), 4)*(y2berr*y2berr)));
      if (y2b/y2berr<2) {
         hBrRaw[ij]->SetBinContent(ipt+1, 0);
         hBrRaw[ij]->SetBinError(ipt+1, 0);
      } 
      // cout <<"bincenter: " << hBr[ij]->GetBinCenter(ipt+1) << " content:"<< hBr[ij]->GetBinContent(ipt+1)<<" err:"<<hBr[ij]->GetBinError(ipt+1)<< "y2b significance:"<<y2b/y2berr<<" yield "<<y2b << " y2bpt "<< y2bpt<< " y3b yield "<< y3b<<endl;
    } 
    hBrRaw[ij]->Draw();
    hBrRaw[ij]->GetYaxis()->SetTitle("R_{3} Raw");
    drawLatex( 0.65,0.88 , Form("%0.2f<y<%0.2f", ybin[ij+1],ybin[ij]) , 0.055);
    addpdf(pdf);

    hBr[ij]=(TH1F*)hPurityCor[ij]->Clone(Form("hBr_%d",ij));
    hBr[ij]->SetDirectory(0);
    g[ij] = (TGraphErrors*)f2b->Get(Form("t_sgct1_corr_yield[0][%d]", ij));
    g[ij]->Draw();
    for (int ipt=0;ipt<npt;ipt++)
    {
      double y3b = hBr[ij]->GetBinContent(ipt+1);
      double y3berr = hBr[ij]->GetBinError(ipt+1);
      double y2b = g[ij]->GetPointY(ipt);
      double y2bpt = g[ij]->GetPointX(ipt);
      double y2berr = g[ij]->GetErrorY(ipt);
      if (ij==1) {y2berr=y2berr*10; y2b=y2b*10.;}
      if (ij==2) {y2berr=y2berr*100; y2b=y2b*100.;}
      // cout <<"y2b significance:"<<y2b/y2berr<<" yield "<<y2b << endl;
      hBr[ij]->SetBinContent(ipt+1, y2b/(y2b+y3b) );
      hBr[ij]->SetBinError(ipt+1, sqrt((y2b*y2b)/pow((y2b+y3b), 4)*(y3berr*y3berr) + (y3b*y3b)/pow((y2b+y3b), 4)*(y2berr*y2berr)));
      if (y2b/y2berr<2) {
         hBr[ij]->SetBinContent(ipt+1, 0);
         hBr[ij]->SetBinError(ipt+1, 0);
      } 
      cout <<"bincenter: " << hBr[ij]->GetBinCenter(ipt+1) << " content:"<< hBr[ij]->GetBinContent(ipt+1)<<" err:"<<hBr[ij]->GetBinError(ipt+1)<< "y2b significance:"<<y2b/y2berr<<" yield "<<y2b << " y2bpt "<< y2bpt<< " y3b yield "<< y3b<<endl;
    } 
    hBr[ij]->Draw();
    hBr[ij]->GetYaxis()->SetTitle("R_{3}");
    drawLatex( 0.65,0.88 , Form("%0.2f<y<%0.2f", ybin[ij+1],ybin[ij]) , 0.055);
    addpdf(pdf);
  }

  TLegend* legb = new TLegend(0.7, 0.7, 0.9, 0.9);
  TH1F* htemp = new TH1F("htemp","htemp;p_{T} (GeV/c); 2-body/(2-body+3-body)",  5, 1, 3);
  htemp->GetYaxis()->SetRangeUser( 0, 0.5);
  htemp->Draw();
  int color[]={kRed-4, kGreen+2, kBlue-4 };
  int style[]={kOpenCircle, kOpenSquare, kOpenStar };
  for (int i=0;i<3;i++)
  {
    hBr[i]->Draw("same");
    hBr[i]->SetLineColor(color[i]);
    hBr[i]->SetMarkerColor(color[i]);
    legb->AddEntry( hBr[i], Form("%0.2f<y<%0.2f", ybin[i+1], ybin[i]), "pl");
  }
  legb->Draw();
  addpdf(pdf);

  double tmp=0, tmperr = 0;
  double sum=0, wtsum=0, sumerr=0;
  for (int i=0;i<2;i++)
  {
    for (int ip=1;ip<=hBr[i]->GetNbinsX();ip++)
    {
      tmp = hBr[i]->GetBinContent(ip);
      if (tmp<1e-31) continue; 
      tmperr = hBr[i]->GetBinError(ip);
      cout <<tmperr <<" "<<tmp << endl;
      sum += tmp/tmperr/tmperr;
      sumerr += 1./(tmperr*tmperr);
      wtsum++;
    }
  } 
  Br = sum/sumerr;
  error = sqrt(1./sumerr);
  cout<<"after correct:" << sum/sumerr<<" "<<sqrt(1./sumerr)<<endl;

  tmp=0; tmperr = 0;
  sum=0; wtsum=0; sumerr=0;
  for (int i=0;i<2;i++)
  {
    for (int ip=1;ip<=hBr[i]->GetNbinsX();ip++)
    {
      tmp = hBrRaw[i]->GetBinContent(ip);
      if (tmp<1e-31) continue; 
      tmperr = hBrRaw[i]->GetBinError(ip);
      cout <<tmperr <<" "<<tmp << endl;
      sum += tmp/tmperr/tmperr;
      sumerr += 1./(tmperr*tmperr);
      wtsum++;
    }
  } 
  cout<<"Before correct:" << sum/sumerr<<" "<<sqrt(1./sumerr)<<endl;

  TFile * fout = new TFile(Form("fout_0050_scan_%d.root", icut),"recreate");
  fout->cd();
  for (int i=0;i<3;i++){
    hPhase[i]->Write();
    hPhaseCor[i]->Write();
    hPurityCor[i]->Write();
    hBr[i]->Write();
    heff[i]->Write();
  }
  fout->Close();

  pdf->On();
  pdf->Close();

}
void drawMixDataScanSys(int icut, TString histname, TString topohistname, TString pdfname, double topocut, double& Br, double & error)
{
  // double highpt = 2.5, lowpt = 1., lowy=-0.8, highy = -0.2;
  // double highpt = 2.5, lowpt = 1, lowy=-0.9, highy = -0.1;
  double highpt = 4, lowpt = 0, lowy=-1.5, highy = 0;
  // TString histname=Form("hH3LMassPtYTopoCut%d",icut);
  cout << histname.Data()<<endl;
  // TString histname="hH3LMassPtY_5_40";
  SetsPhenixStyle();
  TCanvas* c = new TCanvas("c1","c1");
  // TPDF* pdf = new TPDF("MixEventQA_check.pdf");
  // TPDF* pdf = new TPDF("MixEventQA_Jul27.pdf");
  // TPDF* pdf = new TPDF("MixEventQA_beforeDcacut.pdf");
  // TPDF* pdf = new TPDF(Form("MixEventQA_0050_scan_%d.pdf",icut ));
  TPDF* pdf = new TPDF(pdfname);
  pdf->Off();
  gStyle->SetPalette(1);

  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02_0050_015pt_scan.root"); 
  TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02_0050_015pt_sys.root"); 
  TH3F* h2sig = (TH3F*)f1->Get(histname.Data())->Clone("hptH3Lmass_sig");
  h2sig->SetDirectory(0);
  // TH3F* h2sig->Project3D("xz");
  TH1F* hsig = (TH1F*)h2sig->ProjectionY("hsig", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2sig->GetZaxis()->FindBin(lowy), h2sig->GetZaxis()->FindBin(highy));
  hsig->SetDirectory(0);
  TH1F* hcent_se = (TH1F*)f1->Get("hcent")->Clone("hcent_se");
  double nEvents_se = hcent_se->Integral( 4, 9);
  // double nEvents_se = hcent_se->Integral( 1, 9);

  TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Sep02_0050_015pt_sys.root"); 
  TH3F* h2bk = (TH3F*)f2->Get(histname.Data())->Clone("hptH3Lmass_bk");
  h2bk->SetDirectory(0);
  TH1F* hbk = (TH1F*)h2bk->ProjectionY("hbk", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2bk->GetZaxis()->FindBin(lowy), h2bk->GetZaxis()->FindBin(highy));
  hbk->SetDirectory(0);

  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02_0050_015pt_scan.root");
  TFile* f3 = TFile::Open("rootfile/fout_H3L_data_SE_Sep02_0050_015pt_sys.root");
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
  double sig_sb =  hsig->Integral(hsig->GetXaxis()->FindBin(3.01),  hsig->GetXaxis()->FindBin(3.02));
  // double rt_sb =  hrt->Integral(hrt->GetXaxis()->FindBin(2.97), hrt->GetXaxis()->FindBin(2.98)) + hrt->Integral(hrt->GetXaxis()->FindBin(3.0),  hrt->GetXaxis()->FindBin(3.02));
  double rt_sb = hrt->Integral(hrt->GetXaxis()->FindBin(3.01),  hrt->GetXaxis()->FindBin(3.02));
  // double rt_sb =  hrt->Integral(hrt->GetXaxis()->FindBin(2.96), hrt->GetXaxis()->FindBin(2.98)) ;
  double scale_rt = sig_sb/rt_sb;
  // sig_sb = hsig->Integral(hsig->GetXaxis()->FindBin(2.97), hsig->GetXaxis()->FindBin(2.98)) + hsig->Integral(hsig->GetXaxis()->FindBin(3.0),  hsig->GetXaxis()->FindBin(3.02));
  double bk_sb =   hbk->Integral(hbk->GetXaxis()->FindBin(3.01),  hbk->GetXaxis()->FindBin(3.02));
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
  // hrt->Draw("same");
  TLegend* leg_sig = new TLegend(0.65,0.25,0.88,0.45);
  leg_sig->AddEntry(hbk, "ME", "pl");
  // leg_sig->AddEntry(hrt, "RT", "pl");
  // leg_sig->AddEntry(hrt, "KF", "pl");
  leg_sig->AddEntry(hsig, "SE", "pl");
  leg_sig->Draw();


  drawLatex( 0.65,0.61,Form("%0.1f<y<%0.1f",lowy, highy ), 0.055);
  drawLatex( 0.65,0.54,Form("%0.1f<p_{T}<%0.1f GeV/c^{2}", lowpt, highpt), 0.055);
  // drawLatex( 0.65,0.47,Form("5-40%s", "%"), 0.055);
  drawLatex( 0.65,0.47,Form("0-50%s", "%"), 0.055);
  drawBox( 2.97, hsig->GetMinimum(),2.98, hsig->GetMaximum()*0.5 , kBlue-9, 1001 , 1,0.3); 
  // drawBox( 2.989, hsig_bk->GetMinimum(),2.995, hsig_bk->GetMaximum()*0.5 , kOrange-4, 1001,1,0.3); 
  drawBox( 3.01, hsig->GetMinimum(),3.02, hsig->GetMaximum()*0.5 , kBlue-9, 1001,1,0.3 ); 

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
  // hsig_rt->Draw("same");
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
  // hsig_rt->Fit(resfit_rt,"R");
  fit_rt->SetParameter(3, resfit_rt->GetParameter(0) );
  fit_rt->SetParameter(4, resfit_rt->GetParameter(1) );
  // hsig_rt->Fit(fit_rt,"R");
  resfit_rt->SetParameter(0, fit_rt->GetParameter(3));
  resfit_rt->SetParameter(1, fit_rt->GetParameter(4));
  // resfit_rt->Draw("same");

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
  double significance = yield_me/sqrt(sp_counts);
  // double significance = yield_rt/sqrt(sp_counts);
  double s_me = yield_me/(fit->GetParError(0)/hsig_bk->GetBinWidth(1));
  double s_rt = yield_rt/(fit_rt->GetParError(0)/hsig_rt->GetBinWidth(1));
  cout<<"significance: " <<significance << endl;

  TLegend* leg = new TLegend( 0.72, 0.68 ,0.9,0.9 );
  // leg->AddEntry(hbk, "mix-event(ME)","pl");
  // leg->AddEntry(hrt, "rotate d(RT) (scale)","pl");
  leg->AddEntry(hsig, "SE","pl");
  leg->AddEntry(hsig_bk, "SE-ME","pl");
  // leg->AddEntry(hsig_rt, "SE-KF","pl");
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
  drawBox( 3.01, hsig_bk->GetMinimum(),3.02, hsig_bk->GetMaximum()*0.5 , kBlue-9, 1001,1,0.3 ); 

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

  TH2F* hYield = new TH2F( "hYield", "hYield;y,pt", 5, -1., 0., 6, 0, 3);
  double xbw = hYield->GetXaxis()->GetBinWidth(1);
  double ybw = hYield->GetYaxis()->GetBinWidth(1);
  for (int ix=1; ix<=hYield->GetNbinsX(); ix++) {
    for (int iy=1; iy<=hYield->GetNbinsY();iy++){
      double xlow = hYield->GetXaxis()->GetBinLowEdge(ix); // y
      double ylow = hYield->GetYaxis()->GetBinLowEdge(iy); // pt
      double err;
      hYield->SetBinContent( ix , iy, fityield2( histname,  ylow, ylow+ybw, xlow, xlow+xbw, err, f1, f2, c, pdf) );
      hYield->SetBinError( ix , iy, err);
    }
  }
  hYield->Draw("colz text");
  hYield->SaveAs(Form("fH3L_yield_0050_sys_%d.root",icut));
  addpdf(pdf);

  // TFile* fMc_ld = TFile::Open("fout_Lambda_MC_Cuts_0050_015pt_sys.root");
  TFile* fMc_ld = TFile::Open("fout_Lambda_MC_Cuts_0050_015pt_sys_mixpid3.root");
  // TFile* fMc_ld = TFile::Open("fout_Lambda_MC_Cuts_0050_015pt_sys_cut2.root");
  // TFile* fMc_ld = TFile::Open("fout_Lambda_MC_Cuts_0050_015pt_sys.root");
  TFile* fMc = TFile::Open("fout_H3L_MC_0050_015pt_sys.root");
  // TFile* fMc = TFile::Open("fout_H3L_MC_0050_015pt_sys.root");
  // TFile* fMc = TFile::Open("fout_H3L_MC_0080_010pt.root");
  TH1F* hPhase[3];
  TH1F* hPhaseCor[3];
  TH1F* hPurityCor[3];
  TH1F* hBr[3];
  TH1F* hBrRaw[3];
  TGraphErrors* g[3];
  double edge[3][4]={
    { 1.7, 2.1, 2.7, 3.1}, //1.9, 2.4, 2.9
    { 1.2, 1.6, 2.0, 2.4},
    { 1.2, 1.8, 2.2, 2.6}
  };

  double ybin[4]={ 0, -0.25, -0.50, -0.75};
  hPhase[0] = new TH1F("hPhase0", "hPhase0;pt", 25, 1, 3.5);
  hPhase[0]=(TH1F*)hPhase[0]->Rebin(3, "hPhase0", edge[0]);
  hPhase[1] = new TH1F("hPhase1", "hPhase1;pt", 25, 1, 3.5);
  hPhase[1]=(TH1F*)hPhase[1]->Rebin(3, "hPhase1", edge[1]);
  hPhase[2] = new TH1F("hPhase2", "hPhase2;pt", 25, 1, 3.5);
  hPhase[2]=(TH1F*)hPhase[2]->Rebin(3, "hPhase2", edge[2]);
  TFile* fMcH3L = new TFile("fMC_H3L_0050.root");
  // TFile* fRcH3L = new TFile("fout_H3L_MC_0050_015pt_sys.root");
  TFile* fRcH3L = new TFile("fout_H3L_MC_0050_015pt_sys.root");
  // TFile* fMcH3L = new TFile("fMC_H3L_0080.root");
  // TFile* fRcH3L = new TFile("fout_H3L_MC_0080_010pt.root");
  TH3F* h3Mc = (TH3F*)fMcH3L->Get("hH3LMassPtY")->Clone("h3Mc");
  h3Mc->SetDirectory(0);
  TH3F* h3Rc = (TH3F*)fRcH3L->Get(histname.Data())->Clone("h3Rc");
  h3Rc->SetDirectory(0);

  h3Mc->Sumw2();
  h3Rc->Sumw2();
  TH2F* h2MC = (TH2F*)h3Mc->Project3D("xz");
  TH2F* h2Rc = (TH2F*)h3Rc->Project3D("xz");

  TH2F* h2Eff = (TH2F*)h2Rc->Clone("h2Eff");
  TH2F* h2temp = (TH2F*)h2MC->Clone("h2temp");
  h2temp->RebinY(10);
  h2temp->RebinX(10);
  h2Eff->RebinY(10);
  h2Eff->RebinX(10);
  h2Eff->Divide(h2temp);
  h2Eff->Draw("col text");
  h2Eff->GetYaxis()->SetRangeUser(0,4.5);
  addpdf(pdf);

  TH1F* heff[3];
  for (int ij=0;ij<3;ij++){ 
    TH1F* h1Mc = (TH1F*)h2MC->ProjectionY(Form("h1Mc%d", ij), h2MC->GetXaxis()->FindBin(ybin[ij+1] +1e-6), h2MC->GetXaxis()->FindBin(ybin[ij]-1e-6) );
    h1Mc = (TH1F*)h1Mc->Rebin( 3, Form("hMc%d",ij), edge[ij]);
    TH1F* h1Rc = (TH1F*)h2Rc->ProjectionY(Form("h1Rc%d", ij), h2Rc->GetXaxis()->FindBin(ybin[ij+1]+1e-6), h2Rc->GetXaxis()->FindBin(ybin[ij]-1e-6) );
    heff[ij] = (TH1F*)h1Rc->Rebin( 3, Form("heff%d",ij), edge[ij]);
    heff[ij]->Divide(h1Mc);
    heff[ij]->GetYaxis()->SetTitle("Eff.");
    heff[ij]->Draw();
    addpdf(pdf);
  }

  double ptesterr[7];
  double chi2topo[7]={0.2,0.6,1,1.5,2,2.5,3 };
  double puritytest[7];
  for (int i=0;i<7;i++) {
    puritytest[i] = calpurity("h3H3L_chi2topo",0,chi2topo[i],  3., 1., -0.25, -0.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 0, 0, 10, ptesterr[i], 0.25);
    // puritytest[i] = calpurity("h3H3L_chi2topo",0, 3,  3., 1., -0., -0.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, ptesterr[i], 0.25);
  }
  for (int i=0;i<7;i++) 
    cout <<chi2topo[i]<<" "<<puritytest[i]<<" "<<ptesterr[i] << endl;
  TGraphErrors* gpurity = new TGraphErrors(7, chi2topo , puritytest,0, ptesterr  );
  gpurity->Draw("pA");
  gpurity->GetXaxis()->SetTitle("Chi2Topo");
  gpurity->GetYaxis()->SetTitle("Purity");
  addpdf(pdf);

  double dy=0.25;
  double dpt,pt;
  
  TFile* f2b = TFile::Open("h_h3l_corr_yield.root_cut00002_cent0_yuehang.root");

  double ny = 3;
  for (int ij=0; ij<ny;ij++){
    double npt =3;
    if (ij==0) npt =2; 
    for (int ipt=0; ipt<npt;ipt++){
      double err;
      dpt = edge[ij][ipt+1]-edge[ij][ipt];
      pt = 0.5*(edge[ij][ipt+1]+edge[ij][ipt]);
      double yield = fityield2(histname, edge[ij][ipt], edge[ij][ipt+1], ybin[ij+1], ybin[ij], err, f1, f2, c, pdf);
      // cout << edge[ij][ipt]<<" "<<edge[ij][ipt+1]<<" "<<ybin[ij+1]<<" "<<ybin[ij]<<" "<< yield<<" "<<err <<endl;
      if (ij==0 && ipt==2) { yield=0; err=0;}
      hPhase[ij]->SetBinContent( ipt+1, yield/dy/dpt/pt/2./3.1415926/nEvents_se);
      hPhase[ij]->SetBinError( ipt+1, err/dy/dpt/pt/2./3.1415926/nEvents_se );
    }
    hPhase[ij]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPhase[ij]->GetYaxis()->SetTitle("Raw d^{2}N/(N_{evt}2#pip_{T}dydp_{T})");
    hPhase[ij]->Draw();
    hPhase[ij]->SetDirectory(0);
    addpdf(pdf); 

    hPhaseCor[ij]=(TH1F*)hPhase[ij]->Clone(Form("hYieldCor_%d",ij));
    hPhaseCor[ij]->SetDirectory(0);
    hPurityCor[ij]=(TH1F*)hPhase[ij]->Clone(Form("hPurityCor_%d",ij));
    hPurityCor[ij]->SetDirectory(0);

    double perr;
    // double purity1_s = calpurityCut("h3H3L_chi2ndf", "hH3LMassPtY",0, 3.5, edge[ij][3], edge[ij][0], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 10, perr, 0.25);
    // double purity2_s = calpurityCut(Form("h3H3L_chi2ndfSysCut%d", icut), histname,0, 3.5,  edge[ij][3],  edge[ij][0], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 10, perr, 0.25);
    double purity2_s = calpurityCut( topohistname, histname , 0, topocut, edge[ij][3], edge[ij][0], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 0, 0, 10, perr, 0.25);
    // double purity2_s = calpurityCut( topohistname, "hH3LMassPtY", 0, topocut, edge[ij][3], edge[ij][0], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
    // double purity3_s = calpurity("h3H3L_l",8, 99,  3.,  1., ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L decaylength", 2, 0, 100, perr, 0.3);

    for (int ipt=0;ipt<npt;ipt++)
    {
      double y3b = hPhaseCor[ij]->GetBinContent(ipt+1);
      double y3berr = hPhaseCor[ij]->GetBinError(ipt+1)/y3b;
      double eff = heff[ij]->GetBinContent(ipt+1);
      double efferr = heff[ij]->GetBinError(ipt+1)/eff;
      // double purity1 = calpurity("h3H3L_chi2ndf",0, 3.5, edge[ij][ipt+1], edge[ij][ipt], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 10, perr, 0.25);
      // double purity2 = calpurity("h3H3L_chi2topo", 0, chi2topo[icut], edge[ij][ipt+1], edge[ij][ipt], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
    // double purity2_s = calpurityCut( topohistname, histname , 0, topocut, edge[ij][ipt+1], edge[ij][ipt], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.);
      // double purity = (purity1+purity2)*0.5;
      // perr = fabs(purity2-purity1);
      double purity = purity2_s;
      perr = perr/purity;
      // double purity = purity2;
      // double purity = purity1_s;
      // double purity=1;
      // perr = 0;
      double yield_cor = y3b/eff;
      double err = sqrt(y3berr*y3berr + efferr*efferr)*yield_cor;
      hPhaseCor[ij]->SetBinContent(ipt+1,yield_cor );
      hPhaseCor[ij]->SetBinError(ipt+1, err);
      yield_cor = y3b/eff*purity;
      err = sqrt(y3berr*y3berr + efferr*efferr + perr*perr )*yield_cor;
      hPurityCor[ij]->SetBinContent(ipt+1,yield_cor );
      hPurityCor[ij]->SetBinError(ipt+1, err);

      cout <<"pt "<<hPhaseCor[ij]->GetBinCenter(ipt+1) << "purity "<<purity <<" eff " <<eff <<" yield " << yield_cor << " err "<< err<< " s="<<yield_cor/err<<endl;
    }
    hPhaseCor[ij]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPhaseCor[ij]->GetYaxis()->SetTitle("d^{2}N/(N_{evt}2#pip_{T}dydp_{T}(no purity corr.))");
    hPurityCor[ij]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPurityCor[ij]->GetYaxis()->SetTitle("d^{2}N/(N_{evt}2#pip_{T}dydp_{T})");

    hPhaseCor[ij]->Draw();
    hPhaseCor[ij]->SetDirectory(0);
    addpdf(pdf);

    hBrRaw[ij]=(TH1F*)hPhaseCor[ij]->Clone(Form("hBrRaw_%d",ij));
    hBrRaw[ij]->SetDirectory(0);
    g[ij] = (TGraphErrors*)f2b->Get(Form("t_sgct1_corr_yield[0][%d]", ij));
    // g[ij]->Draw();
    for (int ipt=0;ipt<npt;ipt++)
    {
      double y3b = hBrRaw[ij]->GetBinContent(ipt+1);
      double y3berr = hBrRaw[ij]->GetBinError(ipt+1);
      double y2b = g[ij]->GetPointY(ipt);
      double y2bpt = g[ij]->GetPointX(ipt);
      double y2berr = g[ij]->GetErrorY(ipt);
      if (ij==1) {y2berr=y2berr*10; y2b=y2b*10.;}
      if (ij==2) {y2berr=y2berr*100; y2b=y2b*100.;}
      // cout <<"y2b significance:"<<y2b/y2berr<<" yield "<<y2b << endl;
      hBrRaw[ij]->SetBinContent(ipt+1, y2b/(y2b+y3b) );
      hBrRaw[ij]->SetBinError(ipt+1, sqrt((y2b*y2b)/pow((y2b+y3b), 4)*(y3berr*y3berr) + (y3b*y3b)/pow((y2b+y3b), 4)*(y2berr*y2berr)));
      if (y2b/y2berr<2) {
         hBrRaw[ij]->SetBinContent(ipt+1, 0);
         hBrRaw[ij]->SetBinError(ipt+1, 0);
      } 
      // cout <<"bincenter: " << hBr[ij]->GetBinCenter(ipt+1) << " content:"<< hBr[ij]->GetBinContent(ipt+1)<<" err:"<<hBr[ij]->GetBinError(ipt+1)<< "y2b significance:"<<y2b/y2berr<<" yield "<<y2b << " y2bpt "<< y2bpt<< " y3b yield "<< y3b<<endl;
    } 
    hBrRaw[ij]->Draw();
    hBrRaw[ij]->GetYaxis()->SetTitle("R_{3} Raw");
    drawLatex( 0.65,0.88 , Form("%0.2f<y<%0.2f", ybin[ij+1],ybin[ij]) , 0.055);
    addpdf(pdf);

    hBr[ij]=(TH1F*)hPurityCor[ij]->Clone(Form("hBr_%d",ij));
    hBr[ij]->SetDirectory(0);
    g[ij] = (TGraphErrors*)f2b->Get(Form("t_sgct1_corr_yield[0][%d]", ij));
    g[ij]->Draw();
    for (int ipt=0;ipt<npt;ipt++)
    {
      double y3b = hBr[ij]->GetBinContent(ipt+1);
      double y3berr = hBr[ij]->GetBinError(ipt+1);
      double y2b = g[ij]->GetPointY(ipt);
      double y2bpt = g[ij]->GetPointX(ipt);
      double y2berr = g[ij]->GetErrorY(ipt);
      if (ij==1) {y2berr=y2berr*10; y2b=y2b*10.;}
      if (ij==2) {y2berr=y2berr*100; y2b=y2b*100.;}
      // cout <<"y2b significance:"<<y2b/y2berr<<" yield "<<y2b << endl;
      hBr[ij]->SetBinContent(ipt+1, y2b/(y2b+y3b) );
      hBr[ij]->SetBinError(ipt+1, sqrt((y2b*y2b)/pow((y2b+y3b), 4)*(y3berr*y3berr) + (y3b*y3b)/pow((y2b+y3b), 4)*(y2berr*y2berr)));
      if (y2b/y2berr<2) {
         hBr[ij]->SetBinContent(ipt+1, 0);
         hBr[ij]->SetBinError(ipt+1, 0);
      } 
      cout <<"bincenter: " << hBr[ij]->GetBinCenter(ipt+1) << " content:"<< hBr[ij]->GetBinContent(ipt+1)<<" err:"<<hBr[ij]->GetBinError(ipt+1)<< "y2b significance:"<<y2b/y2berr<<" yield "<<y2b << " y2bpt "<< y2bpt<< " y3b yield "<< y3b<<endl;
    } 
    hBr[ij]->Draw();
    hBr[ij]->GetYaxis()->SetTitle("R_{3}");
    drawLatex( 0.65,0.88 , Form("%0.2f<y<%0.2f", ybin[ij+1],ybin[ij]) , 0.055);
    addpdf(pdf);
  }

  TLegend* legb = new TLegend(0.7, 0.7, 0.9, 0.9);
  TH1F* htemp = new TH1F("htemp","htemp;p_{T} (GeV/c); 2-body/(2-body+3-body)",  5, 1, 3);
  htemp->GetYaxis()->SetRangeUser( 0, 0.5);
  htemp->Draw();
  int color[]={kRed-4, kGreen+2, kBlue-4 };
  int style[]={kOpenCircle, kOpenSquare, kOpenStar };
  for (int i=0;i<3;i++)
  {
    hBr[i]->Draw("same");
    hBr[i]->SetLineColor(color[i]);
    hBr[i]->SetMarkerColor(color[i]);
    legb->AddEntry( hBr[i], Form("%0.2f<y<%0.2f", ybin[i+1], ybin[i]), "pl");
  }
  legb->Draw();
  addpdf(pdf);

  //treat different measurements as seperate 
  // double tmp=0, tmperr = 0;
  // double sum=0, wtsum=0, sumerr=0;
  // for (int i=1;i<2;i++)
  // {
  //   for (int ip=1;ip<=hBr[i]->GetNbinsX();ip++)
  //   {
  //     tmp = hBr[i]->GetBinContent(ip);
  //     if (tmp<1e-31) continue; 
  //     tmperr = hBr[i]->GetBinError(ip);
  //     cout <<tmperr <<" "<<tmp << endl;
  //     sum += tmp/tmperr/tmperr;
  //     sumerr += 1./(tmperr*tmperr);
  //     wtsum++;
  //   }
  // } 
  // Br = sum/sumerr;
  // error = sqrt(1./sumerr);
  // cout<<"after correct:" << sum/sumerr<<" "<<sqrt(1./sumerr)<<endl;
  //
  // tmp=0; tmperr = 0;
  // sum=0; wtsum=0; sumerr=0;
  // for (int i=1;i<2;i++)
  // {
  //   for (int ip=1;ip<=hBr[i]->GetNbinsX();ip++)
  //   {
  //     tmp = hBrRaw[i]->GetBinContent(ip);
  //     if (tmp<1e-31) continue; 
  //     tmperr = hBrRaw[i]->GetBinError(ip);
  //     cout <<tmperr <<" "<<tmp << endl;
  //     sum += tmp/tmperr/tmperr;
  //     sumerr += 1./(tmperr*tmperr);
  //     wtsum++;
  //   }
  // } 
  // cout<<"Before correct:" << sum/sumerr<<" "<<sqrt(1./sumerr)<<endl;

  //Direct combine yields
  double sum=0, err=0, y2b=0,y2berr=0, tmp=0, tmperr=0, tmp2=0,tmperr2=0;
  for (int i=0;i<2;i++)
  {
    double sf;
    if (i==0) sf=1;
    else if (i==1) sf=10;
    else if (i==2) sf=100;
    for (int ip=1;ip<=hBr[i]->GetNbinsX();ip++)
    {
      tmp = hPurityCor[i]->GetBinContent(ip);
      if (tmp<1e-31) continue; 
      tmp*=hPurityCor[i]->GetBinCenter(ip)*hPurityCor[i]->GetBinWidth(ip)*2*3.1415*0.25;
      tmperr = hPurityCor[i]->GetBinError(ip)*hPurityCor[i]->GetBinCenter(ip)*hPurityCor[i]->GetBinWidth(ip)*2*3.1415*0.25;
      tmp2 = g[i]->GetPointY(ip-1)*hPurityCor[i]->GetBinCenter(ip)*hPurityCor[i]->GetBinWidth(ip)*2*3.1415*0.25*sf;
      tmperr2 = g[i]->GetErrorY(ip-1)*hPurityCor[i]->GetBinCenter(ip)*hPurityCor[i]->GetBinWidth(ip)*2*3.1415*0.25*sf;
      
      y2b+=tmp2;
      y2berr=tmperr2*tmperr2;
      sum+=tmp;
      err+=tmperr*tmperr;
    }
  } 

  err = sqrt(err);
  y2berr = sqrt(y2berr);
  Br = y2b/(y2b+sum);
  error =BrErr(y2b, y2berr, sum, err); 
  cout << "Br: "<<Br<<" y2berr: "<<y2berr/y2b<<" y3b: "<<err/sum<<" total Br:"<<error<<" "<<error/Br <<endl;

  TFile * fout = new TFile(Form("fout_0050_sys_%d.root", icut),"recreate");
  fout->cd();
  for (int i=0;i<3;i++){
    hPhase[i]->Write();
    hPhaseCor[i]->Write();
    hPurityCor[i]->Write();
    hBr[i]->Write();
    heff[i]->Write();
  }
  fout->Close();

  pdf->On();
  pdf->Close();

}
void drawMixDataScanOld(int icut, double& Br, double & error)
{
  // double highpt = 2.5, lowpt = 1., lowy=-0.8, highy = -0.2;
  // double highpt = 2.5, lowpt = 1, lowy=-0.9, highy = -0.1;
  double highpt = 4, lowpt = 0, lowy=-1.5, highy = 0;
  TString histname=Form("hH3LMassPtYCut%d",icut);
  cout << histname.Data()<<endl;
  // TString histname="hH3LMassPtY_5_40";
  SetsPhenixStyle();
  TCanvas* c = new TCanvas("c1","c1");
  // TPDF* pdf = new TPDF("MixEventQA_check.pdf");
  // TPDF* pdf = new TPDF("MixEventQA_Jul27.pdf");
  // TPDF* pdf = new TPDF("MixEventQA_beforeDcacut.pdf");
  TPDF* pdf = new TPDF(Form("MixEventQA_0050_scan_%d.pdf",icut ));
  pdf->Off();
  gStyle->SetPalette(1);

  TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02_0050_015pt_scan.root"); 
  TH3F* h2sig = (TH3F*)f1->Get(histname.Data())->Clone("hptH3Lmass_sig");
  h2sig->SetDirectory(0);
  // TH3F* h2sig->Project3D("xz");
  TH1F* hsig = (TH1F*)h2sig->ProjectionY("hsig", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2sig->GetZaxis()->FindBin(lowy), h2sig->GetZaxis()->FindBin(highy));
  hsig->SetDirectory(0);
  TH1F* hcent_se = (TH1F*)f1->Get("hcent")->Clone("hcent_se");
  double nEvents_se = hcent_se->Integral( 4, 9);
  // double nEvents_se = hcent_se->Integral( 1, 9);

  TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Sep02_0050_015pt_scan.root"); 
  TH3F* h2bk = (TH3F*)f2->Get(histname.Data())->Clone("hptH3Lmass_bk");
  h2bk->SetDirectory(0);
  TH1F* hbk = (TH1F*)h2bk->ProjectionY("hbk", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2bk->GetZaxis()->FindBin(lowy), h2bk->GetZaxis()->FindBin(highy));
  hbk->SetDirectory(0);

  TFile* f3 = TFile::Open("rootfile/fout_H3L_data_KF_Sep02_0050_015pt_scan.root");
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
  double sig_sb =  hsig->Integral(hsig->GetXaxis()->FindBin(3.01),  hsig->GetXaxis()->FindBin(3.02));
  // double rt_sb =  hrt->Integral(hrt->GetXaxis()->FindBin(2.97), hrt->GetXaxis()->FindBin(2.98)) + hrt->Integral(hrt->GetXaxis()->FindBin(3.0),  hrt->GetXaxis()->FindBin(3.02));
  double rt_sb = hrt->Integral(hrt->GetXaxis()->FindBin(3.01),  hrt->GetXaxis()->FindBin(3.02));
  // double rt_sb =  hrt->Integral(hrt->GetXaxis()->FindBin(2.96), hrt->GetXaxis()->FindBin(2.98)) ;
  double scale_rt = sig_sb/rt_sb;
  // sig_sb = hsig->Integral(hsig->GetXaxis()->FindBin(2.97), hsig->GetXaxis()->FindBin(2.98)) + hsig->Integral(hsig->GetXaxis()->FindBin(3.0),  hsig->GetXaxis()->FindBin(3.02));
  double bk_sb =   hbk->Integral(hbk->GetXaxis()->FindBin(3.01),  hbk->GetXaxis()->FindBin(3.02));
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
  // hrt->Draw("same");
  TLegend* leg_sig = new TLegend(0.65,0.25,0.88,0.45);
  leg_sig->AddEntry(hbk, "ME", "pl");
  // leg_sig->AddEntry(hrt, "RT", "pl");
  // leg_sig->AddEntry(hrt, "KF", "pl");
  leg_sig->AddEntry(hsig, "SE", "pl");
  leg_sig->Draw();


  drawLatex( 0.65,0.61,Form("%0.1f<y<%0.1f",lowy, highy ), 0.055);
  drawLatex( 0.65,0.54,Form("%0.1f<p_{T}<%0.1f GeV/c^{2}", lowpt, highpt), 0.055);
  // drawLatex( 0.65,0.47,Form("5-40%s", "%"), 0.055);
  drawLatex( 0.65,0.47,Form("0-50%s", "%"), 0.055);
  drawBox( 2.97, hsig->GetMinimum(),2.98, hsig->GetMaximum()*0.5 , kBlue-9, 1001 , 1,0.3); 
  // drawBox( 2.989, hsig_bk->GetMinimum(),2.995, hsig_bk->GetMaximum()*0.5 , kOrange-4, 1001,1,0.3); 
  drawBox( 3.01, hsig->GetMinimum(),3.02, hsig->GetMaximum()*0.5 , kBlue-9, 1001,1,0.3 ); 

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
  // hsig_rt->Draw("same");
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
  // hsig_rt->Fit(resfit_rt,"R");
  fit_rt->SetParameter(3, resfit_rt->GetParameter(0) );
  fit_rt->SetParameter(4, resfit_rt->GetParameter(1) );
  // hsig_rt->Fit(fit_rt,"R");
  resfit_rt->SetParameter(0, fit_rt->GetParameter(3));
  resfit_rt->SetParameter(1, fit_rt->GetParameter(4));
  // resfit_rt->Draw("same");

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
  double significance = yield_me/sqrt(sp_counts);
  // double significance = yield_rt/sqrt(sp_counts);
  double s_me = yield_me/(fit->GetParError(0)/hsig_bk->GetBinWidth(1));
  double s_rt = yield_rt/(fit_rt->GetParError(0)/hsig_rt->GetBinWidth(1));
  cout<<"significance: " <<significance << endl;

  TLegend* leg = new TLegend( 0.72, 0.68 ,0.9,0.9 );
  // leg->AddEntry(hbk, "mix-event(ME)","pl");
  // leg->AddEntry(hrt, "rotate d(RT) (scale)","pl");
  leg->AddEntry(hsig, "SE","pl");
  leg->AddEntry(hsig_bk, "SE-ME","pl");
  // leg->AddEntry(hsig_rt, "SE-KF","pl");
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
  drawBox( 3.01, hsig_bk->GetMinimum(),3.02, hsig_bk->GetMaximum()*0.5 , kBlue-9, 1001,1,0.3 ); 

  addpdf(pdf);

  hsig_bk->Draw();
  // hsig_rt->Draw("same");
  TLegend* leg2 = new TLegend( 0.2 , 0.7 ,0.4,0.9  );
  leg2->AddEntry(hsig_rt, "SE-RT","pl");
  leg2->AddEntry(hsig_bk, "SE-ME","pl");
  leg2->Draw();

  addpdf(pdf);
  //
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

  TH2F* hYield = new TH2F( "hYield", "hYield;y,pt", 5, -1., 0., 6, 0, 3);
  double xbw = hYield->GetXaxis()->GetBinWidth(1);
  double ybw = hYield->GetYaxis()->GetBinWidth(1);
  for (int ix=1; ix<=hYield->GetNbinsX(); ix++) {
    for (int iy=1; iy<=hYield->GetNbinsY();iy++){
      double xlow = hYield->GetXaxis()->GetBinLowEdge(ix); // y
      double ylow = hYield->GetYaxis()->GetBinLowEdge(iy); // pt
      double err;
      hYield->SetBinContent( ix , iy, fityield2("hH3LMassPtYTopoCut", icut,  ylow, ylow+ybw, xlow, xlow+xbw, err, f1, f2, c, pdf) );
      hYield->SetBinError( ix , iy, err);
    }
  }
  hYield->Draw("colz text");
  hYield->SaveAs(Form("fH3L_yield_0050_%d.root",icut));
  addpdf(pdf);

  TFile* fMc_ld = TFile::Open("fout_Lambda_MC_Cuts_0050_015pt.root");
  TFile* fMc = TFile::Open("fout_H3L_MC_0050_015pt.root");
  // TFile* fMc = TFile::Open("fout_H3L_MC_0080_010pt.root");
  TH1F* hPhase[3];
  TH1F* hPhaseCor[3];
  TH1F* hPurityCor[3];
  TH1F* hBr[3];
  TH1F* hBrRaw[3];
  TGraphErrors* g[3];
  double edge[3][4]={
    { 1.7, 2.1, 2.7, 3.1},
    { 1.2, 1.6, 2.0, 2.4},
    { 1.2, 1.8, 2.2, 2.6}
  };

  double ybin[4]={ 0, -0.25, -0.50, -0.75};
  hPhase[0] = new TH1F("hPhase0", "hPhase0;pt", 25, 1, 3.5);
  hPhase[0]=(TH1F*)hPhase[0]->Rebin(3, "hPhase0", edge[0]);
  hPhase[1] = new TH1F("hPhase1", "hPhase1;pt", 25, 1, 3.5);
  hPhase[1]=(TH1F*)hPhase[1]->Rebin(3, "hPhase1", edge[1]);
  hPhase[2] = new TH1F("hPhase2", "hPhase2;pt", 25, 1, 3.5);
  hPhase[2]=(TH1F*)hPhase[2]->Rebin(3, "hPhase2", edge[2]);
  TFile* fMcH3L = new TFile("fMC_H3L_0050.root");
  TFile* fRcH3L = new TFile("fout_H3L_MC_0050_015pt.root");
  // TFile* fMcH3L = new TFile("fMC_H3L_0080.root");
  // TFile* fRcH3L = new TFile("fout_H3L_MC_0080_010pt.root");
  TH3F* h3Mc = (TH3F*)fMcH3L->Get("hH3LMassPtY")->Clone("h3Mc");
  h3Mc->SetDirectory(0);
  TH3F* h3Rc = (TH3F*)fRcH3L->Get(histname.Data())->Clone("h3Rc");
  h3Rc->SetDirectory(0);

  h3Mc->Sumw2();
  h3Rc->Sumw2();
  TH2F* h2MC = (TH2F*)h3Mc->Project3D("xz");
  TH2F* h2Rc = (TH2F*)h3Rc->Project3D("xz");

  TH2F* h2Eff = (TH2F*)h2Rc->Clone("h2Eff");
  TH2F* h2temp = (TH2F*)h2MC->Clone("h2temp");
  h2temp->RebinY(10);
  h2temp->RebinX(10);
  h2Eff->RebinY(10);
  h2Eff->RebinX(10);
  h2Eff->Divide(h2temp);
  h2Eff->Draw("col text");
  h2Eff->GetYaxis()->SetRangeUser(0,4.5);
  addpdf(pdf);

  TH1F* heff[3];
  for (int ij=0;ij<3;ij++){ 
    TH1F* h1Mc = (TH1F*)h2MC->ProjectionY(Form("h1Mc%d", ij), h2MC->GetXaxis()->FindBin(ybin[ij+1] +1e-6), h2MC->GetXaxis()->FindBin(ybin[ij]-1e-6) );
    h1Mc = (TH1F*)h1Mc->Rebin( 3, Form("hMc%d",ij), edge[ij]);
    TH1F* h1Rc = (TH1F*)h2Rc->ProjectionY(Form("h1Rc%d", ij), h2Rc->GetXaxis()->FindBin(ybin[ij+1]+1e-6), h2Rc->GetXaxis()->FindBin(ybin[ij]-1e-6) );
    heff[ij] = (TH1F*)h1Rc->Rebin( 3, Form("heff%d",ij), edge[ij]);
    heff[ij]->Divide(h1Mc);
    heff[ij]->GetYaxis()->SetTitle("Eff.");
    heff[ij]->Draw();
    addpdf(pdf);
  }

  double ptesterr[7];
  double chi2topo[7]={0.2,0.6,1,1.5,2,2.5,3 };
  double puritytest[7];
  for (int i=0;i<7;i++) {
    puritytest[i] = calpurity("h3H3L_chi2topo",0,chi2topo[i],  3., 1., -0.25, -0.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, ptesterr[i], 0.25);
    // puritytest[i] = calpurity("h3H3L_chi2topo",0, 3,  3., 1., -0., -0.5, f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, ptesterr[i], 0.25);
  }
  for (int i=0;i<7;i++) 
    cout <<chi2topo[i]<<" "<<puritytest[i]<<" "<<ptesterr[i] << endl;
  TGraphErrors* gpurity = new TGraphErrors(7, chi2topo , puritytest,0, ptesterr  );
  gpurity->Draw("pA");
  gpurity->GetXaxis()->SetTitle("Chi2Topo");
  gpurity->GetYaxis()->SetTitle("Purity");
  addpdf(pdf);

  double dy=0.25;
  double dpt,pt;
  
  TFile* f2b = TFile::Open("h_h3l_corr_yield.root_cut00002_cent0_yuehang.root");

  double ny = 3;
  for (int ij=0; ij<ny;ij++){
    double npt =3;
    if (ij==0) npt =2; 
    for (int ipt=0; ipt<npt;ipt++){
      double err;
      dpt = edge[ij][ipt+1]-edge[ij][ipt];
      pt = 0.5*(edge[ij][ipt+1]+edge[ij][ipt]);
      double yield = fityield2("hH3LMassPtYTopoCut", icut, edge[ij][ipt], edge[ij][ipt+1], ybin[ij+1], ybin[ij], err, f1, f2, c, pdf);
      // cout << edge[ij][ipt]<<" "<<edge[ij][ipt+1]<<" "<<ybin[ij+1]<<" "<<ybin[ij]<<" "<< yield<<" "<<err <<endl;
      if (ij==0 && ipt==2) { yield=0; err=0;}
      hPhase[ij]->SetBinContent( ipt+1, yield/dy/dpt/pt/2./3.1415926/nEvents_se);
      hPhase[ij]->SetBinError( ipt+1, err/dy/dpt/pt/2./3.1415926/nEvents_se );
    }
    hPhase[ij]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPhase[ij]->GetYaxis()->SetTitle("Raw d^{2}N/(N_{evt}2#pip_{T}dydp_{T})");
    hPhase[ij]->Draw();
    hPhase[ij]->SetDirectory(0);
    addpdf(pdf); 

    hPhaseCor[ij]=(TH1F*)hPhase[ij]->Clone(Form("hYieldCor_%d",ij));
    hPhaseCor[ij]->SetDirectory(0);
    hPurityCor[ij]=(TH1F*)hPhase[ij]->Clone(Form("hPurityCor_%d",ij));
    hPurityCor[ij]->SetDirectory(0);

    double perr;
    double purity1_s = calpurity("h3H3L_chi2ndf",0, 3.5, edge[ij][3], edge[ij][0], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 10, perr, 0.25);
    double purity2_s = calpurity(Form("h3H3L_chi2ndfTopoCut%d", icut),0, 3.5, edge[ij][3], edge[ij][0], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 10, perr, 0.25);
    // double purity2_s = calpurity("h3H3L_chi2topo",0, 3.,  edge[ij][3],  edge[ij][0], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
    // double purity2_s = calpurity("h3H3L_chi2topo",0, chi2topo[icut], edge[ij][2], edge[ij][0], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
    // double puritytest = calpurity("h3H3L_chi2topo",0, 3.,  3.,  1.,ybin[0], ybin[2], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
    // double purity3_s = calpurity("h3H3L_l",8, 99,  3.,  1., ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L decaylength", 2, 0, 100, perr, 0.3);

    for (int ipt=0;ipt<npt;ipt++)
    {
      double y3b = hPhaseCor[ij]->GetBinContent(ipt+1);
      double y3berr = hPhaseCor[ij]->GetBinError(ipt+1)/y3b;
      double eff = heff[ij]->GetBinContent(ipt+1);
      double efferr = heff[ij]->GetBinError(ipt+1)/eff;
      // double purity1 = calpurity("h3H3L_chi2ndf",0, 3.5, edge[ij][ipt+1], edge[ij][ipt], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{NDF}", 2, 0, 10, perr, 0.25);
      // double purity2 = calpurity("h3H3L_chi2topo", 0, chi2topo[icut], edge[ij][ipt+1], edge[ij][ipt], ybin[ij], ybin[ij+1], f1, f2, scale, fMc, fMc_ld, c, pdf, "H3L #chi^{2}_{topo}", 2, 0, 10, perr, 0.25);
      // double purity = (purity1+purity2)*0.5;
      // perr = fabs(purity2-purity1);
      // perr = perr/purity;
      double purity = purity2_s;
      // double purity = purity2;
      // double purity = purity1_s;
      // double purity;
      // purity=1;
      // perr = 0;
      double yield_cor = y3b/eff;
      double err = sqrt(y3berr*y3berr + efferr*efferr)*yield_cor;
      hPhaseCor[ij]->SetBinContent(ipt+1,yield_cor );
      hPhaseCor[ij]->SetBinError(ipt+1, err);
      yield_cor = y3b/eff*purity;
      err = sqrt(y3berr*y3berr + efferr*efferr + perr*perr )*yield_cor;
      hPurityCor[ij]->SetBinContent(ipt+1,yield_cor );
      hPurityCor[ij]->SetBinError(ipt+1, err);

      cout <<"pt "<<hPhaseCor[ij]->GetBinCenter(ipt+1) << "purity "<<purity <<" eff " <<eff <<" yield " << yield_cor << " err "<< err<< " s="<<yield_cor/err<<endl;
    }
    hPhaseCor[ij]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPhaseCor[ij]->GetYaxis()->SetTitle("d^{2}N/(N_{evt}2#pip_{T}dydp_{T}(no purity corr.))");
    hPurityCor[ij]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hPurityCor[ij]->GetYaxis()->SetTitle("d^{2}N/(N_{evt}2#pip_{T}dydp_{T})");

    hPhaseCor[ij]->Draw();
    hPhaseCor[ij]->SetDirectory(0);
    addpdf(pdf);

    hBrRaw[ij]=(TH1F*)hPhaseCor[ij]->Clone(Form("hBrRaw_%d",ij));
    hBrRaw[ij]->SetDirectory(0);
    g[ij] = (TGraphErrors*)f2b->Get(Form("t_sgct1_corr_yield[0][%d]", ij));
    // g[ij]->Draw();
    for (int ipt=0;ipt<npt;ipt++)
    {
      double y3b = hBrRaw[ij]->GetBinContent(ipt+1);
      double y3berr = hBrRaw[ij]->GetBinError(ipt+1);
      double y2b = g[ij]->GetPointY(ipt);
      double y2bpt = g[ij]->GetPointX(ipt);
      double y2berr = g[ij]->GetErrorY(ipt);
      if (ij==1) {y2berr=y2berr*10; y2b=y2b*10.;}
      if (ij==2) {y2berr=y2berr*100; y2b=y2b*100.;}
      // cout <<"y2b significance:"<<y2b/y2berr<<" yield "<<y2b << endl;
      hBrRaw[ij]->SetBinContent(ipt+1, y2b/(y2b+y3b) );
      hBrRaw[ij]->SetBinError(ipt+1, sqrt((y2b*y2b)/pow((y2b+y3b), 4)*(y3berr*y3berr) + (y3b*y3b)/pow((y2b+y3b), 4)*(y2berr*y2berr)));
      if (y2b/y2berr<2) {
         hBrRaw[ij]->SetBinContent(ipt+1, 0);
         hBrRaw[ij]->SetBinError(ipt+1, 0);
      } 
      // cout <<"bincenter: " << hBr[ij]->GetBinCenter(ipt+1) << " content:"<< hBr[ij]->GetBinContent(ipt+1)<<" err:"<<hBr[ij]->GetBinError(ipt+1)<< "y2b significance:"<<y2b/y2berr<<" yield "<<y2b << " y2bpt "<< y2bpt<< " y3b yield "<< y3b<<endl;
    } 
    hBrRaw[ij]->Draw();
    hBrRaw[ij]->GetYaxis()->SetTitle("R_{3} Raw");
    drawLatex( 0.65,0.88 , Form("%0.2f<y<%0.2f", ybin[ij+1],ybin[ij]) , 0.055);
    addpdf(pdf);

    hBr[ij]=(TH1F*)hPurityCor[ij]->Clone(Form("hBr_%d",ij));
    hBr[ij]->SetDirectory(0);
    g[ij] = (TGraphErrors*)f2b->Get(Form("t_sgct1_corr_yield[0][%d]", ij));
    g[ij]->Draw();
    for (int ipt=0;ipt<npt;ipt++)
    {
      double y3b = hBr[ij]->GetBinContent(ipt+1);
      double y3berr = hBr[ij]->GetBinError(ipt+1);
      double y2b = g[ij]->GetPointY(ipt);
      double y2bpt = g[ij]->GetPointX(ipt);
      double y2berr = g[ij]->GetErrorY(ipt);
      if (ij==1) {y2berr=y2berr*10; y2b=y2b*10.;}
      if (ij==2) {y2berr=y2berr*100; y2b=y2b*100.;}
      // cout <<"y2b significance:"<<y2b/y2berr<<" yield "<<y2b << endl;
      hBr[ij]->SetBinContent(ipt+1, y2b/(y2b+y3b) );
      hBr[ij]->SetBinError(ipt+1, sqrt((y2b*y2b)/pow((y2b+y3b), 4)*(y3berr*y3berr) + (y3b*y3b)/pow((y2b+y3b), 4)*(y2berr*y2berr)));
      if (y2b/y2berr<2) {
         hBr[ij]->SetBinContent(ipt+1, 0);
         hBr[ij]->SetBinError(ipt+1, 0);
      } 
      cout <<"bincenter: " << hBr[ij]->GetBinCenter(ipt+1) << " content:"<< hBr[ij]->GetBinContent(ipt+1)<<" err:"<<hBr[ij]->GetBinError(ipt+1)<< "y2b significance:"<<y2b/y2berr<<" yield "<<y2b << " y2bpt "<< y2bpt<< " y3b yield "<< y3b<<endl;
    } 
    hBr[ij]->Draw();
    hBr[ij]->GetYaxis()->SetTitle("R_{3}");
    drawLatex( 0.65,0.88 , Form("%0.2f<y<%0.2f", ybin[ij+1],ybin[ij]) , 0.055);
    addpdf(pdf);
  }

  TLegend* legb = new TLegend(0.7, 0.7, 0.9, 0.9);
  TH1F* htemp = new TH1F("htemp","htemp;p_{T} (GeV/c); 2-body/(2-body+3-body)",  5, 1, 3);
  htemp->GetYaxis()->SetRangeUser( 0, 0.5);
  htemp->Draw();
  int color[]={kRed, kGreen, kBlue };
  for (int i=0;i<3;i++)
  {
    hBr[i]->Draw("same");
    hBr[i]->SetLineColor(color[i]);
    hBr[i]->SetMarkerColor(color[i]);
    legb->AddEntry( hBr[i], Form("%0.2f<y<%0.2f", ybin[i+1], ybin[i]), "pl");
  }
  legb->Draw();
  addpdf(pdf);

  double tmp=0, tmperr = 0;
  double sum=0, wtsum=0, sumerr=0;
  for (int i=0;i<2;i++)
  {
    for (int ip=1;ip<=hBr[i]->GetNbinsX();ip++)
    {
      tmp = hBr[i]->GetBinContent(ip);
      if (tmp<1e-31) continue; 
      tmperr = hBr[i]->GetBinError(ip);
      cout <<tmperr <<" "<<tmp << endl;
      sum += tmp/tmperr/tmperr;
      sumerr += 1./(tmperr*tmperr);
      wtsum++;
    }
  } 
  Br = sum/sumerr;
  error = sqrt(1./sumerr);
  cout<<"after correct:" << sum/sumerr<<" "<<sqrt(1./sumerr)<<endl;

  tmp=0; tmperr = 0;
  sum=0; wtsum=0; sumerr=0;
  for (int i=0;i<2;i++)
  {
    for (int ip=1;ip<=hBr[i]->GetNbinsX();ip++)
    {
      tmp = hBrRaw[i]->GetBinContent(ip);
      if (tmp<1e-31) continue; 
      tmperr = hBrRaw[i]->GetBinError(ip);
      cout <<tmperr <<" "<<tmp << endl;
      sum += tmp/tmperr/tmperr;
      sumerr += 1./(tmperr*tmperr);
      wtsum++;
    }
  } 
  cout<<"Before correct:" << sum/sumerr<<" "<<sqrt(1./sumerr)<<endl;

  TFile * fout = new TFile(Form("fout_0050_scan_%d.root", icut),"recreate");
  fout->cd();
  for (int i=0;i<3;i++){
    hPhase[i]->Write();
    hPhaseCor[i]->Write();
    hPurityCor[i]->Write();
    hBr[i]->Write();
    heff[i]->Write();
  }
  fout->Close();

  pdf->On();
  pdf->Close();

}
void drawMixDataTest()
{
  // double highpt = 2.5, lowpt = 1., lowy=-0.8, highy = -0.2;
  // double highpt = 2.5, lowpt = 1, lowy=-0.9, highy = -0.1;
  double highpt = 4, lowpt = 0, lowy=-1.5, highy = 0;
  TString histname="hH3LMassPtY";
  // TString histname="hH3LMassPtY_5_40";
  SetsPhenixStyle();
  TCanvas* c = new TCanvas("c1","c1");
  // TPDF* pdf = new TPDF("MixEventQA_check.pdf");
  // TPDF* pdf = new TPDF("MixEventQA_Jul27.pdf");
  // TPDF* pdf = new TPDF("MixEventQA_beforeDcacut.pdf");
  TPDF* pdf = new TPDF("MixEventQA_vertex.pdf");
  pdf->Off();
  gStyle->SetPalette(1);

  // TFile *f1 = TFile::Open("fout_H3L_data_SE_large.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_large.root"); 

  // TFile *f1 = TFile::Open("out_KF_test/fout_H3L_KF_test.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_KF_Jul22_part.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_KF_Jul25.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_SE_Jul25.root"); 
  // TFile *f1 = TFile::Open("fout_H3L_data_SE_Jul27.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Aug09.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_SE_Aug16.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_SE_Aug20.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_SE_Aug24.root"); 
  TFile *f1 = TFile::Open("rootfile/fout_H3L_data_SE_Sep02_0050_015pt.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_SE_Aug22.root"); 
  // TFile *f1 = TFile::Open("fout_KF_test.root"); 
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
  double nEvents_se = hcent_se->Integral( 3, 9);

  // TFile *f2 = TFile::Open("fout_H3L_data_ME_large.root"); 
  // TFile *f2 = TFile::Open("fout_H3L_ME_Jul22_part.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Aug16.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Aug20.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Aug24.root"); 
  TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Sep02_0050_015pt.root"); 
  // TFile *f2 = TFile::Open("rootfile/fout_H3L_data_ME_Aug22.root"); 
  // TFile *f2 = TFile::Open("fout_ME_test.root"); 
  // TFile *f2 = TFile::Open("fout_H3L_ME_Jul25.root"); 
  // TFile *f2 = TFile::Open("fout_H3L_data_SE.root"); 
  // TH2F* h2bk = (TH2F*)f2->Get("hptH3Lmass")->Clone("hptH3Lmass_ME");
  // TH3F* h2bk = (TH3F*)f2->Get("hH3LMassPtY_5_40")->Clone("hptH3Lmass_bk");
  TH3F* h2bk = (TH3F*)f2->Get(histname.Data())->Clone("hptH3Lmass_bk");
  h2bk->SetDirectory(0);
  TH1F* hbk = (TH1F*)h2bk->ProjectionY("hbk", h2sig->GetXaxis()->FindBin(lowpt+1e-6),  h2sig->GetXaxis()->FindBin(highpt-1e-6), h2bk->GetZaxis()->FindBin(lowy), h2bk->GetZaxis()->FindBin(highy));
  hbk->SetDirectory(0);

  // TFile* f3 = TFile::Open("fout_H3L_data_RT_Jul22_full.root");
  // TFile* f3 = TFile::Open("fout_H3L_data_RT_Jul25.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_RT_Jul27.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_KF_Aug20.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_SE_Aug24.root");
  TFile* f3 = TFile::Open("rootfile/fout_H3L_data_SE_Sep02_0050_015pt.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_KF_Aug22.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_KF_Aug16.root");
  // TFile* f3 = TFile::Open("rootfile/fout_H3L_data_RT_Aug13.root");
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
  // double sig_sb = hsig->Integral(hsig->GetXaxis()->FindBin(2.976), hsig->GetXaxis()->FindBin(2.986)) ;
  // double sig_sb = hsig->Integral(hsig->GetXaxis()->FindBin(2.97), hsig->GetXaxis()->FindBin(2.98)) + hsig->Integral(hsig->GetXaxis()->FindBin(3.0),  hsig->GetXaxis()->FindBin(3.02));
  double sig_sb =  hsig->Integral(hsig->GetXaxis()->FindBin(3.00),  hsig->GetXaxis()->FindBin(3.02));
  // double rt_sb =  hrt->Integral(hrt->GetXaxis()->FindBin(2.97), hrt->GetXaxis()->FindBin(2.98)) + hrt->Integral(hrt->GetXaxis()->FindBin(3.0),  hrt->GetXaxis()->FindBin(3.02));
  double rt_sb = hrt->Integral(hrt->GetXaxis()->FindBin(3.00),  hrt->GetXaxis()->FindBin(3.02));
  // double rt_sb =  hrt->Integral(hrt->GetXaxis()->FindBin(2.96), hrt->GetXaxis()->FindBin(2.98)) ;
  double scale_rt = sig_sb/rt_sb;
  // sig_sb = hsig->Integral(hsig->GetXaxis()->FindBin(2.97), hsig->GetXaxis()->FindBin(2.98)) + hsig->Integral(hsig->GetXaxis()->FindBin(3.0),  hsig->GetXaxis()->FindBin(3.02));
  double bk_sb =   hbk->Integral(hbk->GetXaxis()->FindBin(3.00),  hbk->GetXaxis()->FindBin(3.02));
  // double bk_sb =   hbk->Integral(hbk->GetXaxis()->FindBin(2.976),  hbk->GetXaxis()->FindBin(2.986));
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
  //
  TH3F* h2sigv12[7];
  TH3F* h2bkv12[7];
  TH2F* h2v12[7];
  TH2F* h2v12_xy[7];
  TH2F* h2v12_yz[7];
  TH2F* h2v12_xz[7];
  gPad->SetLogz();
  for (int i=2;i<7;i++){
    if (i>=2) h2sigv12[i] = (TH3F*)f1->Get(Form("h3H3L_v12xySig%d", i));
    if (i==1) h2sigv12[i] = (TH3F*)f1->Get("h3H3L_v12xySig");
    h2sigv12[i]->RebinX(10);
    h2sigv12[i]->RebinY(20);
    h2sigv12[i]->RebinZ(10);
    h2sigv12[i]->SetDirectory(0);
    if (i==5) h2sigv12[i]->GetZaxis()->SetTitle("d Dca");
    if (i==6) h2sigv12[i]->GetZaxis()->SetTitle("p Dca");
    h2sigv12[i]->GetYaxis()->SetTitle("dp vertex in XY");
    if (i>=2) h2bkv12[i] = (TH3F*)f2->Get(Form("h3H3L_v12xySig%d",i))->Clone(Form("MEsbr%d",i));
    if (i==1) h2bkv12[i] = (TH3F*)f2->Get("h3H3L_v12xySig")->Clone(Form("MEsbr%d",i));
    
    h2bkv12[i]->RebinX(10);
    h2bkv12[i]->RebinY(20);
    h2bkv12[i]->RebinZ(10);
    h2bkv12[i]->SetDirectory(0);
    h2sigv12[i]->Add( h2bkv12[i], -1*scale);
    if (i==2 || i == 5 || i==6)  h2sigv12[i]->GetZaxis()->SetRangeUser(0, 3.5);
    if (i==3 || i == 4)  h2sigv12[i]->GetZaxis()->SetRangeUser(0, 3.);
    h2v12_xy[i] = (TH2F*)h2sigv12[i]->Project3D("xy");
    h2v12_xy[i]->Draw("colz");
    drawLatex( 0.5, 0.8, "SE-ME 2.989<M<2.995", 0.055);
    addpdf(pdf);
    h2sigv12[i]->GetZaxis()->SetRangeUser(0, 10);
    h2v12_yz[i] = (TH2F*)h2sigv12[i]->Project3D("zy");
    h2v12_yz[i]->Draw("colz");
    drawLatex( 0.5, 0.8, "SE-ME 2.989<M<2.995", 0.055);
    addpdf(pdf);
    h2v12_xz[i] = (TH2F*)h2sigv12[i]->Project3D("xz");
    h2v12_xz[i]->Draw("colz");
    drawLatex( 0.5, 0.8, "SE-ME 2.989<M<2.995", 0.055);
    addpdf(pdf);
    if (i==6) {
      TH1F* htest = (TH1F*)h2v12_xz[i]->ProjectionX("htest", h2v12_xz[i]->GetYaxis()->FindBin(5), h2v12_xz[i]->GetYaxis()->FindBin(10));
      htest->Draw();
    addpdf(pdf);
    }

   }


  for (int i=2;i<7;i++){
    if (i>=2) h2sigv12[i] = (TH3F*)f1->Get(Form("h3H3L_v12xySBR%d", i));
    if (i==1) h2sigv12[i] = (TH3F*)f1->Get("h3H3L_v12xySBR");
    h2sigv12[i]->RebinX(10);
    h2sigv12[i]->RebinY(20);
    h2sigv12[i]->RebinZ(10);
    h2sigv12[i]->SetDirectory(0);
    if (i==5) h2sigv12[i]->GetZaxis()->SetTitle("d Dca");
    if (i==6) h2sigv12[i]->GetZaxis()->SetTitle("p Dca");
    h2sigv12[i]->GetYaxis()->SetTitle("dp vertex in XY");
    if (i>=2) h2bkv12[i] = (TH3F*)f2->Get(Form("h3H3L_v12xySBR%d",i))->Clone(Form("MEsbr%d",i));
    if (i==1) h2bkv12[i] = (TH3F*)f2->Get("h3H3L_v12xySBR")->Clone(Form("MEsbr%d",i));
    
    h2bkv12[i]->RebinX(10);
    h2bkv12[i]->RebinY(20);
    h2bkv12[i]->RebinZ(10);
    h2bkv12[i]->SetDirectory(0);
    // h2sigv12[i]->Add( h2bkv12[i], -1*scale);
    if (i==2 || i == 5 || i==6)  h2sigv12[i]->GetZaxis()->SetRangeUser(0, 3.5);
    if (i==2 || i == 5 || i==6)  h2bkv12[i]->GetZaxis()->SetRangeUser(0, 3.5);
    if (i==3 || i == 4)  h2sigv12[i]->GetZaxis()->SetRangeUser(0, 3.);
    if (i==3 || i == 4)  h2bkv12[i]->GetZaxis()->SetRangeUser(0, 3.);
    h2v12_xy[i] = (TH2F*)h2bkv12[i]->Project3D("xy");
    h2v12_xy[i]->Draw("colz");
    drawLatex( 0.5, 0.8, "ME 3<M<3.02", 0.055);
    addpdf(pdf);
    h2v12_xy[i] = (TH2F*)h2sigv12[i]->Project3D("xy");
    h2v12_xy[i]->Draw("colz");
    drawLatex( 0.5, 0.8, "SE 3<M<3.02", 0.055);
    addpdf(pdf);

    h2bkv12[i]->GetZaxis()->SetRangeUser(0, 10);
    h2v12_yz[i] = (TH2F*)h2bkv12[i]->Project3D("zy");
    h2v12_yz[i]->Draw("colz");
    drawLatex( 0.5, 0.8, "ME 3<M<3.02", 0.055);
    addpdf(pdf);
    h2sigv12[i]->GetZaxis()->SetRangeUser(0, 10);
    h2v12_yz[i] = (TH2F*)h2sigv12[i]->Project3D("zy");
    h2v12_yz[i]->Draw("colz");
    drawLatex( 0.5, 0.8, "SE 3<M<3.02", 0.055);
    addpdf(pdf);

    h2v12_xz[i] = (TH2F*)h2bkv12[i]->Project3D("xz");
    h2v12_xz[i]->Draw("colz");
    drawLatex( 0.5, 0.8, "ME 3<M<3.02", 0.055);
    addpdf(pdf);
    h2v12_xz[i] = (TH2F*)h2sigv12[i]->Project3D("xz");
    h2v12_xz[i]->Draw("colz");
    drawLatex( 0.5, 0.8, "SE 3<M<3.02", 0.055);
    addpdf(pdf);

    if (i==6) {
      TH1F* htest = (TH1F*)h2v12_xz[i]->ProjectionX("htest", h2v12_xz[i]->GetYaxis()->FindBin(5), h2v12_xz[i]->GetYaxis()->FindBin(10));
      htest->Draw();
      addpdf(pdf);
    }

  }

  TH3F* h2sigv02[7];
  TH3F* h2bkv02[7];
  TH2F* h2v02[7];
  TH2F* h2v02_xy[7];
  TH2F* h2v02_yz[7];
  TH2F* h2v02_xz[7];
  gPad->SetLogz();
  for (int i=2;i<7;i++){
    if (i>=2) h2sigv02[i] = (TH3F*)f1->Get(Form("h3H3L_v02xySig%d", i));
    if (i==1) h2sigv02[i] = (TH3F*)f1->Get("h3H3L_v02xySig");
    h2sigv02[i]->RebinX(10);
    h2sigv02[i]->RebinY(20);
    h2sigv02[i]->RebinZ(10);
    h2sigv02[i]->SetDirectory(0);
    if (i==5) h2sigv02[i]->GetZaxis()->SetTitle("d Dca");
    if (i==6) h2sigv02[i]->GetZaxis()->SetTitle("p Dca");
    h2sigv02[i]->GetYaxis()->SetTitle("dpi vertex in XY");
    if (i>=2) h2bkv02[i] = (TH3F*)f2->Get(Form("h3H3L_v02xySig%d",i))->Clone(Form("MEsbr%d",i));
    if (i==1) h2bkv02[i] = (TH3F*)f2->Get("h3H3L_v02xySig")->Clone(Form("MEsbr%d",i));
    
    h2bkv02[i]->RebinX(10);
    h2bkv02[i]->RebinY(20);
    h2bkv02[i]->RebinZ(10);
    h2bkv02[i]->SetDirectory(0);
    h2sigv02[i]->Add( h2bkv02[i], -1*scale);
    if (i==2 || i == 5 || i==6)  h2sigv02[i]->GetZaxis()->SetRangeUser(0, 3.5);
    if (i==3 || i == 4)  h2sigv02[i]->GetZaxis()->SetRangeUser(0, 3.);
    h2v02_xy[i] = (TH2F*)h2sigv02[i]->Project3D("xy");
    h2v02_xy[i]->Draw("colz");
    drawLatex( 0.5, 0.8, "SE-ME 2.989<M<2.995", 0.055);
    addpdf(pdf);
    h2sigv02[i]->GetZaxis()->SetRangeUser(0, 10);
    h2v02_yz[i] = (TH2F*)h2sigv02[i]->Project3D("zy");
    h2v02_yz[i]->Draw("colz");
    drawLatex( 0.5, 0.8, "SE-ME 2.989<M<2.995", 0.055);
    addpdf(pdf);
    h2v02_xz[i] = (TH2F*)h2sigv02[i]->Project3D("xz");
    h2v02_xz[i]->Draw("colz");
    drawLatex( 0.5, 0.8, "SE-ME 2.989<M<2.995", 0.055);
    addpdf(pdf);
    if (i==6) {
      TH1F* htest = (TH1F*)h2v02_xz[i]->ProjectionX("htest", h2v02_xz[i]->GetYaxis()->FindBin(5), h2v02_xz[i]->GetYaxis()->FindBin(10));
      htest->Draw();
    addpdf(pdf);
    }

   }

  for (int i=2;i<7;i++){
    if (i>=2) h2sigv02[i] = (TH3F*)f1->Get(Form("h3H3L_v02xySBR%d", i));
    if (i==1) h2sigv02[i] = (TH3F*)f1->Get("h3H3L_v02xySBR");
    h2sigv02[i]->RebinX(10);
    h2sigv02[i]->RebinY(20);
    h2sigv02[i]->RebinZ(10);
    h2sigv02[i]->SetDirectory(0);
    if (i==5) h2sigv02[i]->GetZaxis()->SetTitle("d Dca");
    if (i==6) h2sigv02[i]->GetZaxis()->SetTitle("p Dca");
    h2sigv02[i]->GetYaxis()->SetTitle("dpi vertex in XY");
    if (i>=2) h2bkv02[i] = (TH3F*)f2->Get(Form("h3H3L_v02xySBR%d",i))->Clone(Form("MEsbr%d",i));
    if (i==1) h2bkv02[i] = (TH3F*)f2->Get("h3H3L_v02xySBR")->Clone(Form("MEsbr%d",i));
    
    h2bkv02[i]->RebinX(10);
    h2bkv02[i]->RebinY(20);
    h2bkv02[i]->RebinZ(10);
    h2bkv02[i]->SetDirectory(0);
    // h2sigv02[i]->Add( h2bkv02[i], -1*scale);
    if (i==2 || i == 5 || i==6)  h2sigv02[i]->GetZaxis()->SetRangeUser(0, 3.5);
    if (i==2 || i == 5 || i==6)  h2bkv02[i]->GetZaxis()->SetRangeUser(0, 3.5);
    if (i==3 || i == 4)  h2sigv02[i]->GetZaxis()->SetRangeUser(0, 3.);
    if (i==3 || i == 4)  h2bkv02[i]->GetZaxis()->SetRangeUser(0, 3.);
    h2v02_xy[i] = (TH2F*)h2bkv02[i]->Project3D("xy");
    h2v02_xy[i]->Draw("colz");
    drawLatex( 0.5, 0.8, "ME 3<M<3.02", 0.055);
    addpdf(pdf);
    h2v02_xy[i] = (TH2F*)h2sigv02[i]->Project3D("xy");
    h2v02_xy[i]->Draw("colz");
    drawLatex( 0.5, 0.8, "SE 3<M<3.02", 0.055);
    addpdf(pdf);

    h2bkv02[i]->GetZaxis()->SetRangeUser(0, 10);
    h2v02_yz[i] = (TH2F*)h2bkv02[i]->Project3D("zy");
    h2v02_yz[i]->Draw("colz");
    drawLatex( 0.5, 0.8, "ME 3<M<3.02", 0.055);
    addpdf(pdf);
    h2sigv02[i]->GetZaxis()->SetRangeUser(0, 10);
    h2v02_yz[i] = (TH2F*)h2sigv02[i]->Project3D("zy");
    h2v02_yz[i]->Draw("colz");
    drawLatex( 0.5, 0.8, "SE 3<M<3.02", 0.055);
    addpdf(pdf);

    h2v02_xz[i] = (TH2F*)h2bkv02[i]->Project3D("xz");
    h2v02_xz[i]->Draw("colz");
    drawLatex( 0.5, 0.8, "ME 3<M<3.02", 0.055);
    addpdf(pdf);
    h2v02_xz[i] = (TH2F*)h2sigv02[i]->Project3D("xz");
    h2v02_xz[i]->Draw("colz");
    drawLatex( 0.5, 0.8, "SE 3<M<3.02", 0.055);
    addpdf(pdf);

    if (i==6) {
      TH1F* htest = (TH1F*)h2v02_xz[i]->ProjectionX("htest", h2v02_xz[i]->GetYaxis()->FindBin(5), h2v02_xz[i]->GetYaxis()->FindBin(10));
      htest->Draw();
      addpdf(pdf);
    }

  }

  pdf->On();
  pdf->Close();

}
void drawSEData()
{
  SetsPhenixStyle();
  TCanvas* c = new TCanvas("c","c", 800,800);
  TPDF* pdf = new TPDF("SameEventQA_ME3bcheck.pdf");
  pdf->Off();

  // TFile *f1 = TFile::Open("fout_H3L_data_KF_Jul27.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_KF_Aug09_test.root"); 
  // TFile *f1 = TFile::Open("rootfile/fout_H3L_data_SE_Aug10_test.root"); 
  TFile *f1 = TFile::Open("fout_KF_test.root"); 
  // TH2F* h2sig = (TH2F*)f1->Get("hptH3Lmass")->Clone("hptH3Lmass_sig");
  TH2F* h2sig = (TH2F*)f1->Get("hptH3Lmass")->Clone("hptH3Lmass_sig");
  TH1F* hsig = (TH1F*)h2sig->ProjectionY("hsig");

  // TFile *f2 = TFile::Open("fout_H3L_data_ME.root"); 
  // TFile *f2 = TFile::Open("fout_H3L_data_SE_Jul27.root"); 
  TFile *f2 = TFile::Open("fout_SE_test.root"); 
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

  projAndComp("hptH3L_lSBR", f1, f2, c,pdf,"l", "p",4,"KF","SE", "dp#pi");
  projAndComp("hptH3L_ldlSBR", f1, f2, c,pdf,"l/dl", "p",4,"KF","SE", "dp#pi");
  projAndComp("hptH3L_chi2ndfSBR", f1, f2, c,pdf,"chi2NDF", "p",4,"KF","SE", "dp#pi");
  projAndComp("hptH3L_chi2topoSBR", f1, f2, c,pdf,"chi2topo", "p",4,"KF","SE", "dp#pi" );
  projAndComp("hptH3L_dchi2primSBR", f1, f2, c,pdf,"d chi2primary", "p",4,"KF","SE", "dp#pi");
  projAndComp("hptH3L_pchi2primSBR", f1, f2, c,pdf,"p chi2primary", "p",4,"KF","SE", "dp#pi");
  projAndComp("hptH3L_pichi2primSBR", f1, f2, c,pdf,"#pi chi2primary", "p",4,"KF","SE", "dp#pi");
  projAndComp("hptH3L_dDcaSBR", f1, f2, c,pdf,"d DCA", "p",1,"KF","SE","dp#pi"  );
  projAndComp("hptH3L_piDcaSBR", f1, f2, c,pdf,"#pi DCA", "p",4,"KF","SE", "dp#pi" );
  projAndComp("hptH3L_pDcaSBR", f1, f2, c,pdf,"#p DCA", "p",4,"KF","SE", "dp#pi" );
  projAndComp("hptH3L_dpDcaSBR", f1, f2, c,pdf,"dp pair DCA", "p",4,"KF","SE", "dp#pi" );

  projAndComp("hptH3L_ppimassSig", f1, f2, c,pdf,"p#pi Mass (GeV/c^{2})", "plhist",1 ,"KF" ,"SE","p#pi");
  projAndComp("hptH3L_ppichi2ndfSig", f1, f2, c,pdf,"p#pi #chi^{2}_{NDF}","" ,1 , "KF" ,"SE","p#pi");
  projAndComp("hptH3L_ppichi2primSig", f1, f2, c,pdf,"(p#pi) #chi^{2}_{prim}" ,"",1 , "KF" ,"SE","p#pi");
  projAndComp("hptH3L_ppilSig", f1, f2, c,pdf,"(p#pi) l" ,"",1 , "KF" ,"SE","p#pi");
  projAndComp("hptH3L_ppildlSig", f1, f2, c,pdf,"(p#pi) l/#deltal" ,"",1 , "KF" ,"SE","p#pi");
  // projAndComp("hptsumdcaSig", f1, f2, c,pdf , "p+pi DCA","" ,2, "KF" ,"My", "p#pi");

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
void checkbkgd()
{
  TFile* f1 = new TFile("rootfile/fout_H3L_data_ME_Aug24.root");
  TH3F* h3Bk = (TH3F*)f1->Get("hH3LMassPtCent")->Clone("hbk");
  h3Bk->SetDirectory(0);
  TFile* f2 = new TFile("rootfile/fout_H3L_data_KF_Aug24.root");
  TH3F* h3Sig = (TH3F*)f2->Get("hH3LMassPtCent");
  h3Sig->SetDirectory(0);

  TH1F* hSig = (TH1F*)h3Sig->ProjectionZ("hSig" , h3Sig->GetXaxis()->FindBin(0+1e-6), h3Sig->GetXaxis()->FindBin(4), h3Sig->GetYaxis()->FindBin(3.00),h3Sig->GetYaxis()->FindBin(3.02) );
  TH1F* hBk = (TH1F*)h3Bk->ProjectionZ("hBk" , h3Sig->GetXaxis()->FindBin(0+1e-6), h3Sig->GetXaxis()->FindBin(4), h3Sig->GetYaxis()->FindBin(3.00),h3Sig->GetYaxis()->FindBin(3.02) );
  Norm(hSig);
  Norm(hBk);
  hSig->Draw();
  hBk->SetLineColor(kRed);
  hBk->Draw("same");

}
void correctEff()
{
  SetsPhenixStyle();
  TCanvas* c = new TCanvas("c","c");
  TPDF* pdf = new TPDF("correctEff.pdf");
  pdf->Off();
  TFile* fMC = new TFile("fMC_H3L_0080.root");
  TFile* fRc = new TFile("fout_H3L_MC_0080_010pt.root");
  // TFile* fMC = new TFile("fMC_H3L_0050.root");
  // TFile* fRc = new TFile("fout_H3L_MC_0050_dcacut.root");

  TH3F* h3Mc = (TH3F*)fMC->Get("hH3LMassPtY")->Clone("h3Mc");
  h3Mc->SetDirectory(0);
  TH3F* h3Rc = (TH3F*)fRc->Get("hH3LMassPtY")->Clone("h3Rc");
  h3Rc->SetDirectory(0);

  h3Mc->GetXaxis()->SetRangeUser(0,3);
  h3Mc->GetZaxis()->SetRangeUser(-1,0);
  h3Rc->GetXaxis()->SetRangeUser(0,3);
  h3Rc->GetZaxis()->SetRangeUser(-1,0);
  TH2F* h2MC = (TH2F*)h3Mc->Project3D("xz");
  TH2F* h2Rc = (TH2F*)h3Rc->Project3D("xz");
  h2MC->RebinY(10);
  h2Rc->RebinY(10);
  h2MC->RebinX(20);
  h2Rc->RebinX(20);
  h2Rc->Draw("colz");

  h2Rc->Divide(h2MC);
  // h2Rc->Draw("colz");
  h2Rc->Draw("colz text");
  addpdf(pdf);
  h2Rc->SetDirectory(0);

  fMC->Close();
  fRc->Close();
 
  TFile* fFY_3braw = new TFile("hH3L3bRaw.root");
  TH2F* hFY_3braw = (TH2F*)fFY_3braw->Get("hH3L3bRaw");
  hFY_3braw->SetDirectory(0);
  fFY_3braw->Close();

  TFile* fFY_3bCor = new TFile("hH3L3bCor.root");
  TH2F* hFY_3bCor = (TH2F*)fFY_3bCor->Get("hH3L3bCor");
  hFY_3bCor->SetDirectory(0);
  fFY_3bCor->Close();

  // TFile* fFY_3braw = new TFile("fH3L_yield_nopt.root");
  // TH2F* hFY_3braw = (TH2F*)fFY_3braw->Get("hYield")->Clone("noptcut");
  // hFY_3braw->SetDirectory(0);
  // fFY_3braw->Close();
  
  TFile* fH3L_yield = new TFile("fH3L_yield_MB.root");
  TH2F* hPhase = (TH2F*)fH3L_yield->Get("hYield");
  hPhase->SetDirectory(0);
  hPhase->Sumw2();
  hPhase->Draw("colz text");
  addpdf(pdf);

  TH2F* hCheck = (TH2F*)hPhase->Clone("hCheck");
  hFY_3braw->Scale(272./276.);
  hFY_3braw->Sumw2();
  hCheck->Divide(hFY_3braw);
  hCheck->Draw("colz text");
  addpdf(pdf);

 
  TLegend* legy = new TLegend(0.2,0.2,0.5,0.5); 
  double yr[] = {-0.8,-0.6,-0.4,-0.2};
  for (int i=2;i<=4;i++)
  {
    // TH1F* h = (TH1F*)hFY_3braw->ProjectionY(Form("hfy%d", i), i, i); 
    TH1F* h2 = (TH1F*)hCheck->ProjectionY(Form("hcheck%d", i), i, i); 
    // h->SetMarkerColor(i);
    // h->SetLineColor(i);
    h2->SetMarkerColor(i);
    h2->SetLineColor(i);
    // h2->SetLineStyle(3);
    // h2->SetMarkerStyle(kOpenCircle);
    // h->GetXaxis()->SetRangeUser(1,2.5);
    // h->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    // if (i>2) h->Draw("same");
    // else h->Draw();
    if (i==2) h2->Draw();
    h2->Draw("same");
    h2->GetYaxis()->SetRangeUser( 0.4, 1.3);
    h2->GetXaxis()->SetRangeUser( 1, 3);
    h2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h2->GetYaxis()->SetTitle("ME/RT");
    legy->AddEntry(h2, Form("%0.1f<y<%0.1f", yr[i-2],yr[i-1]), "le");
  }
  drawLine(1, 1, 3, 1, 1.5, 2, 1);
  legy->Draw();
  addpdf(pdf);

  TH2F* hPhaseCor = (TH2F*)hPhase->Clone("hPhaseCor");
  hPhaseCor->Divide(h2Rc);
  hPhaseCor->Draw("colz text");
  addpdf(pdf);

  TH2F* hCheckCor = (TH2F*)hPhaseCor->Clone("hCheckCor");
  hFY_3bCor->Scale(272./261.);
  hCheckCor->Divide(hFY_3bCor);
  hCheckCor->Draw("colz text");
  addpdf(pdf);

  TFile* f2b = new TFile("hH3L2bCor.root"); 
  TH2F* h2b = (TH2F*)f2b->Get("hH3L2bCor");
  h2b->Draw("colz text");
  addpdf(pdf);
 
  TH1F* hyield = (TH1F*)hPhaseCor->ProjectionX("hPhaseCor", 3, 5 ); 
  // cout <<hyield->Integral(2, 4)<<endl;
  hyield->SetBinContent(1, 0);
  hyield->SetBinError(1, 0);
  hyield->SetBinContent(5, 0);
  hyield->SetBinError(5, 0);
  hyield->Rebin(5);
  hyield->Draw();
  cout <<hyield->GetBinContent(1)<<" "<< hyield->GetBinError(1)<< endl;
  addpdf(pdf);
   
  TH2F* h2Br_fy = (TH2F*)h2b->Clone("h2Br_fy");
  TH2F* hbrfy = (TH2F*)hFY_3bCor->Clone("hbrfy");
  hbrfy->Add(h2Br_fy);
  h2Br_fy->Divide(hbrfy);
  h2Br_fy->Draw("colz text");
  addpdf(pdf);
  

  TH2F* h2Br = (TH2F*)h2b->Clone("h2Br");
  for (int i=1;i<h2Br->GetNbinsX()+1;i++)
  {
    for (int j=1;j<h2Br->GetNbinsY()+1;j++)
    {
      double y2b = h2b->GetBinContent(i, j)/257.;
      // double y2b = 59139/257./1e6;
      // double y2b = 59139.;
      double y3b = hPhase->GetBinContent(i, j)/h2Rc->GetBinContent(i,j)/272;
      // double y3b = 192575./272./1e6;
      // double y3b = 219285.;
      double y2berr = h2b->GetBinError( i, j)/257.;
      // // double y2berr = 6621/257./1e6;
      // double y2berr = 6991.;
      double y3berr = hPhase->GetBinError( i, j)/h2Rc->GetBinContent( i, j)/272;
      // double y3berr = 7365.84/272/1e6;
      // double y3berr = 11850;
      // cout << "y3b " <<y3b<< " err"<<y3berr<< endl;
      // cout << "y2b " <<y2b<< " err"<<y2berr<< " significance: "<< y2b/y2berr<<endl;
      double ratio = y2b/(y3b+y2b);
      double error = sqrt( pow((y2b/(y3b+y2b)/(y3b+y2b))*y3berr, 2) +  pow((y3b/(y3b+y2b)/(y3b+y2b))*y2berr, 2));
      if (y2b/y2berr<2. || y2b<1e-6 || y3b<1e-6) {ratio =0; error = 0;}
      h2Br->SetBinContent(i, j, ratio);
      h2Br->SetBinError(i, j, error);
      cout << "y:" << h2Br->GetXaxis()->GetBinCenter(i)<<" pt:"<<h2Br->GetYaxis()->GetBinCenter(j) <<" ratio:"<< ratio<< " err:"<< error<<endl;
    }
  }

  h2Br->Draw("colz text");
  addpdf(pdf);

  
  pdf->On();
  pdf->Close();
}
void scanCuts()
{
  SetsPhenixStyle();
  double p[7],perr[7];
  double cuts[7]={0.2,0.6,1, 1.5, 2. ,2.5, 3};
  for (int i=1;i<7;i++)
  { 
    // drawMixDataScanTopo(i, Form("hH3LMassPtYTopoCut%d",i ), "h3H3L_chi2topo", cuts[i], p[i], perr[i]); 
    drawMixDataScanTopo(i, Form("hH3LMassPtYTopoCut%d",i ), "h3H3L_chi2topo", cuts[i], p[i], perr[i]); 
    // drawMixDataScanTopo(i, Form("hH3LMassPtYTopoCut%d",i ), Form("h3H3L_chi2ndfTopoCut%d",i), 3.5, p[i], perr[i]); 
  }
  // TGraphErrors* gBr = new TGraphErrors( 6, cuts+1 , p+1, 0, perr+1 );
  // gBr->Draw("pA");
  // gBr->GetXaxis()->SetTitle("Chi2Topo Cuts");
  // gBr->GetYaxis()->SetTitle("R_{3}");
  // gBr->SaveAs("scantopoBr.root");
  //
  // drawMixDataScanTopo(6, p[6], perr[6]); 
  //
  TFile* f[7];
  TH1F* h[7]; 
  TH1F* h2[7]; 
  TH1F* h3[7]; 
  f[6] = new TFile( Form("fout_0050_scan_%d.root", 6));
  h[6] = (TH1F*)f[6]->Get("hPurityCor_1")->Clone(Form("htest%d",6));
  h2[6] = (TH1F*)f[6]->Get("hPurityCor_0")->Clone(Form("htest2%d",6));
  // h[6] = (TH1F*)f[6]->Get("hYieldCor_1")->Clone(Form("htest%d",6));

  TFile* f2b = new TFile("h_h3l_corr_yield.root_cut00002_cent0_yuehang.root");
  TGraphErrors* g2b[2];
  g2b[1] = (TGraphErrors*)f2b->Get("t_sgct1_corr_yield[0][1]");
  g2b[0] = (TGraphErrors*)f2b->Get("t_sgct1_corr_yield[0][0]");
  f2b->Close();

  double y2b=0, y2berr=0;
  double binwidth1[3]={ 0.4, 0.4, 0.4};
  double binwidth0[3]={ 0.4, 0.6, 0.4};
  for (int i=0;i<3;i++){
    double x, y,err;
    g2b[1]->GetPoint( i, x, y);
    err=g2b[1]->GetErrorY(i)*x*binwidth1[i]*2*3.1415*10*0.25;
    y2b+=y*x*binwidth1[i]*2*3.1415*10*0.25;
    y2berr+=err*err;
  }
  for (int i=0;i<2;i++){
    double x, y,err;
    g2b[0]->GetPoint( i, x, y);
    err=g2b[0]->GetErrorY(i)*x*binwidth0[i]*2*3.1415*0.25;
    y2b+=y*x*binwidth0[i]*2*3.1415*0.25;
    y2berr+=err*err;
  }
  y2berr=sqrt(y2berr);
  cout<< "2 body yield: " << y2b<<" "<< y2berr<< " relative error " << y2berr/y2b<< endl;

  TH1F* hRatio[7];
  double br[7], brerr[7], y3brelerr[7];
  int color[7]={ kGray+2, kMagenta, kRed, kGreen+2, kBlue, kOrange+5, kBlack};
  TLegend* leg = new TLegend(0.2,0.2,0.5,0.5);
  double yield[7], error[7];
  for (int i=1;i<7;i++)
  {
    f[i] = new TFile( Form("fout_0050_scan_%d.root", i));
    h[i] = (TH1F*)f[i]->Get("hPurityCor_1")->Clone(Form("htest%d",i));
    h2[i] = (TH1F*)f[i]->Get("hPurityCor_0")->Clone(Form("htest2%d",i));
    h3[i] = (TH1F*)f[i]->Get("hPurityCor_2")->Clone(Form("htest2%d",i));
    // h3[i] = (TH1F*)f[i]->Get("hYieldCor_2")->Clone(Form("htest2%d",i));
    // h2[i] = (TH1F*)f[i]->Get("hYieldCor_0")->Clone(Form("htest2%d",i));
    // h[i] = (TH1F*)f[i]->Get("hYieldCor_1")->Clone(Form("htest%d",i));
    h[i]->SetDirectory(0);
    h2[i]->SetDirectory(0);
    h3[i]->SetDirectory(0);
    h[i]->SetLineColor(color[i]);
    h[i]->SetMarkerColor(color[i]);
    if (i>1) h[i]->Draw("same");
    else h[i]->Draw();
    f[i]->Close();
    leg->AddEntry(h[i], Form("chi2topo<%0.1f", cuts[i]), "pe");
    // hRatio[i] = (TH1F*)h[i]->Clone(Form("hr%d",i));
    // hRatio[i]->Divide(h[6]);
    // if (i==1) hRatio[i]->Draw();
    // hRatio[i]->Draw("same");
    // hRatio[i]->GetYaxis()->SetTitle("Yield of CurrentCuts/Yield of Chi2Topo<3");
    // leg->AddEntry(hRatio[i], Form("chi2topo<%0.1f", cuts[i]), "pe");
    yield[i]=0;error[i]=0;
    for (int ib=1;ib<=3;ib++)
    {
      yield[i]+=h[i]->GetBinContent(ib)*h[i]->GetBinWidth(ib)*h[i]->GetBinCenter(ib)*2.*3.1415*0.25;
      error[i]+=pow(h[i]->GetBinError(ib)*h[i]->GetBinWidth(ib)*h[i]->GetBinCenter(ib)*2.*3.1415*0.25, 2);
    }
    for (int ib=1;ib<=2;ib++)
    {
      if (h2[i]->GetBinContent(ib)<1e-31) continue;
      yield[i]+=h2[i]->GetBinContent(ib)*h2[i]->GetBinWidth(ib)*h2[i]->GetBinCenter(ib)*2.*3.1415*0.25;
      error[i]+=pow(h2[i]->GetBinError(ib)*h2[i]->GetBinWidth(ib)*h2[i]->GetBinCenter(ib)*2.*3.1415*0.25,2);
    }
    // for (int ib=1;ib<=3;ib++)
    // {
    //   // if (h3[i]->GetBinContent(ib)<1e-31) continue;
    //   yield[i]+=h3[i]->GetBinContent(ib)*h3[i]->GetBinWidth(ib)*h3[i]->GetBinCenter(ib)*2.*3.1415;
    //   error[i]+=pow(h3[i]->GetBinError(ib)*h3[i]->GetBinWidth(ib)*h3[i]->GetBinCenter(ib)*2.*3.1415,2);
    // }

    error[i]=sqrt(error[i]);
    y3brelerr[i]=error[i]/yield[i];
    br[i] = y2b/(yield[i]+y2b);
    brerr[i] = sqrt( pow((y2b/(yield[i]+y2b)/(yield[i]+y2b))*error[i], 2) +  pow((yield[i]/(yield[i]+y2b)/(yield[i]+y2b))*y2berr, 2));
  }
  leg->Draw();
  TGraphErrors* gyield = new TGraphErrors(6, cuts+1,yield+1, 0, error+1);
  gyield->SetMarkerStyle(kOpenCircle);
  gyield->Draw("pa");
  gyield->GetYaxis()->SetTitle("Yield");
  gyield->GetXaxis()->SetTitle("Chi2Topo Cuts");
  double sys=0, n=0, mean=0;
  for (int ip=0;ip<gyield->GetN();ip++)
  {
     double delta = fabs(yield[ip+1] - yield[4]);
     double ds = sqrt( fabs(error[ip+1]*error[ip+1]-error[4]*error[4]));
     // if (ds>delta) cout << "pass check "<< ip<< endl;
     // else {
          sys+=delta*delta; n++;
     // }
     mean+=yield[ip+1];
     // cout <<yield[ip+1] <<endl;
  }
  // cout << gyield->GetN()<< endl;
  // mean = mean*1.0/gyield->GetN();
  cout << "systematic from chi2topo cuts: "<< sqrt(sys/n)/yield[4]<< " "<< mean<<endl;
  drawLine(0.4,mean,3.1,mean, 2,1,2);
  drawLine(0.4,yield[4]-sqrt(sys/n),3.1,yield[4]-sqrt(sys/n), 2,2,2);
  drawLine(0.4,yield[4]+sqrt(sys/n),3.1,yield[4]+sqrt(sys/n), 2,2,2);
  return;

  TGraphErrors* gbr = new TGraphErrors(6, cuts+1,br+1, 0, brerr+1);
  gbr->SetMarkerStyle(kOpenCircle);
  gbr->Draw("pa");
  gbr->GetYaxis()->SetTitle("R_{3}");
  gbr->GetXaxis()->SetTitle("Chi2Topo Cuts");
  TGraphErrors* gy3berr = new TGraphErrors(6, cuts+1,y3brelerr+1, 0, 0);
  gy3berr->SetMarkerStyle(kOpenCircle);
  gy3berr->Draw("pa l");
  gy3berr->GetYaxis()->SetTitle("Relative uncertainty");
  gy3berr->GetXaxis()->SetTitle("chi2Topo Cuts");

}
void scanCutsNDF()
{
  SetsPhenixStyle();
  double p[10],perr[10];
  double cuts[10]={1, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5};
  for (int i=0;i<10;i++)
  { 
    // drawMixDataScanTopo(i, Form("hH3LMassPtYTopoCut%d",i ), "h3H3L_chi2topo", cuts[i], p[i], perr[i]); 
    // drawMixDataScanNDF(i, Form("hH3LMassPtYNDFCut%d",i ), Form("h3H3L_chi2topoNDFCut%d", i), 2, p[i], perr[i]); 
    // drawMixDataScanNDF(i, Form("hH3LMassPtYNDFCut%d",i ), "h3H3L_chi2ndf", cuts[i], p[i], perr[i]); 
  }
    // drawMixDataScanNDF(9, Form("hH3LMassPtYNDFCut%d", 9 ), Form("h3H3L_chi2topoNDFCut%d", 9), 2, p[9], perr[9]); 
  // TGraphErrors* gBr = new TGraphErrors( 9, cuts+1 , p+1, 0, perr+1 );
  // gBr->Draw("pA");
  // gBr->GetXaxis()->SetTitle("Chi2Topo Cuts");
  // gBr->GetYaxis()->SetTitle("R_{3}");
  // gBr->SaveAs("scantopoBr.root");
  //
  // drawMixDataScanTopo(6, p[6], perr[6]); 
  //
  TFile* f[10];
  TH1F* h[10]; 
  TH1F* h2[10]; 
  TH1F* h3[10]; 
  f[9] = new TFile( Form("fout_0050_ndfscan_%d.root", 9));
  h[9] = (TH1F*)f[9]->Get("hPurityCor_1")->Clone(Form("htest%d",9));
  h2[9] = (TH1F*)f[9]->Get("hPurityCor_0")->Clone(Form("htest2%d",9));
  // h[9] = (TH1F*)f[9]->Get("hYieldCor_1")->Clone(Form("htest%d",6));

  TFile* f2b = new TFile("h_h3l_corr_yield.root_cut00002_cent0_yuehang.root");
  TGraphErrors* g2b[2];
  g2b[1] = (TGraphErrors*)f2b->Get("t_sgct1_corr_yield[0][1]");
  g2b[0] = (TGraphErrors*)f2b->Get("t_sgct1_corr_yield[0][0]");
  f2b->Close();

  double y2b=0, y2berr=0;
  double binwidth1[3]={ 0.4, 0.4, 0.4};
  double binwidth0[3]={ 0.4, 0.6, 0.4};
  for (int i=0;i<3;i++){
    double x, y,err;
    g2b[1]->GetPoint( i, x, y);
    err=g2b[1]->GetErrorY(i)*x*binwidth1[i]*2*3.1415*10*0.25;
    y2b+=y*x*binwidth1[i]*2*3.1415*10*0.25;
    y2berr+=err*err;
  }
  for (int i=0;i<2;i++){
    double x, y,err;
    g2b[0]->GetPoint( i, x, y);
    err=g2b[0]->GetErrorY(i)*x*binwidth0[i]*2*3.1415*0.25;
    y2b+=y*x*binwidth0[i]*2*3.1415*0.25;
    y2berr+=err*err;
  }
  y2berr=sqrt(y2berr);
  cout<< "2 body yield: " << y2b<<" "<< y2berr<< " relative error " << y2berr/y2b<< endl;

  TH1F* hRatio[10];
  double br[10], brerr[10], y3brelerr[7];
  int color[10]={ kGray+2, kMagenta, kRed, kGreen+2, kBlue, kOrange+5, kBlack};
  TLegend* leg = new TLegend(0.2,0.2,0.5,0.5);
  double yield[10], error[10];
  for (int i=0;i<10;i++)
  {
    f[i] = new TFile( Form("fout_0050_ndfscan_%d.root", i));
    h[i] = (TH1F*)f[i]->Get("hPurityCor_1")->Clone(Form("htest%d",i));
    h2[i] = (TH1F*)f[i]->Get("hPurityCor_0")->Clone(Form("htest2%d",i));
    h3[i] = (TH1F*)f[i]->Get("hPurityCor_2")->Clone(Form("htest2%d",i));
    // h3[i] = (TH1F*)f[i]->Get("hYieldCor_2")->Clone(Form("htest2%d",i));
    // h2[i] = (TH1F*)f[i]->Get("hYieldCor_0")->Clone(Form("htest2%d",i));
    // h[i] = (TH1F*)f[i]->Get("hYieldCor_1")->Clone(Form("htest%d",i));
    h[i]->SetDirectory(0);
    h2[i]->SetDirectory(0);
    h3[i]->SetDirectory(0);
    h[i]->SetLineColor(color[i]);
    h[i]->SetMarkerColor(color[i]);
    if (i>1) h[i]->Draw("same");
    else h[i]->Draw();
    f[i]->Close();
    leg->AddEntry(h[i], Form("chi2ndf<%0.1f", cuts[i]), "pe");
    // hRatio[i] = (TH1F*)h[i]->Clone(Form("hr%d",i));
    // hRatio[i]->Divide(h[6]);
    // if (i==1) hRatio[i]->Draw();
    // hRatio[i]->Draw("same");
    // hRatio[i]->GetYaxis()->SetTitle("Yield of CurrentCuts/Yield of Chi2Topo<3");
    // leg->AddEntry(hRatio[i], Form("chi2topo<%0.1f", cuts[i]), "pe");
    yield[i]=0;error[i]=0;
    for (int ib=1;ib<=3;ib++)
    {
      yield[i]+=h[i]->GetBinContent(ib)*h[i]->GetBinWidth(ib)*h[i]->GetBinCenter(ib)*2.*3.1415*0.25;
      error[i]+=pow(h[i]->GetBinError(ib)*h[i]->GetBinWidth(ib)*h[i]->GetBinCenter(ib)*2.*3.1415*0.25, 2);
    }
    for (int ib=1;ib<=2;ib++)
    {
      if (h2[i]->GetBinContent(ib)<1e-31) continue;
      yield[i]+=h2[i]->GetBinContent(ib)*h2[i]->GetBinWidth(ib)*h2[i]->GetBinCenter(ib)*2.*3.1415*0.25;
      error[i]+=pow(h2[i]->GetBinError(ib)*h2[i]->GetBinWidth(ib)*h2[i]->GetBinCenter(ib)*2.*3.1415*0.25,2);
    }
    // for (int ib=1;ib<=3;ib++)
    // {
    //   // if (h3[i]->GetBinContent(ib)<1e-31) continue;
    //   yield[i]+=h3[i]->GetBinContent(ib)*h3[i]->GetBinWidth(ib)*h3[i]->GetBinCenter(ib)*2.*3.1415;
    //   error[i]+=pow(h3[i]->GetBinError(ib)*h3[i]->GetBinWidth(ib)*h3[i]->GetBinCenter(ib)*2.*3.1415,2);
    // }

    error[i]=sqrt(error[i]);
    y3brelerr[i]=error[i]/yield[i];
    br[i] = y2b/(yield[i]+y2b);
    brerr[i] = sqrt( pow((y2b/(yield[i]+y2b)/(yield[i]+y2b))*error[i], 2) +  pow((yield[i]/(yield[i]+y2b)/(yield[i]+y2b))*y2berr, 2));
  }
  leg->Draw();
  TGraphErrors* gyield = new TGraphErrors(10, cuts,yield, 0, error);
  gyield->SetMarkerStyle(kOpenCircle);
  gyield->Draw("pa");
  gyield->GetYaxis()->SetTitle("Yield");
  gyield->GetXaxis()->SetTitle("Chi2NDF Cuts");
  double sys=0, n=0, mean=0;
  for (int ip=0;ip<gyield->GetN();ip++)
  {
     double delta = fabs(yield[ip+1] - yield[4]);
     double ds = sqrt( fabs(error[ip+1]*error[ip+1]-error[4]*error[4]));
     // if (ds>delta) cout << "pass check "<< ip<< endl;
     // else {
          sys+=delta*delta; n++;
     // }
     mean+=yield[ip+1];
     // cout <<yield[ip+1] <<endl;
  }
  // cout << gyield->GetN()<< endl;
  // mean = mean*1.0/gyield->GetN();
  cout << "systematic from chi2topo cuts: "<< sqrt(sys/n)/yield[9]<< " "<< mean<<endl;
  drawLine(0.4,mean,3.1,mean, 2,1,2);
  drawLine(0.4,yield[9]-sqrt(sys/n),3.1,yield[9]-sqrt(sys/n), 2,2,2);
  drawLine(0.4,yield[9]+sqrt(sys/n),3.1,yield[9]+sqrt(sys/n), 2,2,2);
  return;

  TGraphErrors* gbr = new TGraphErrors(10, cuts,br, 0, brerr);
  gbr->SetMarkerStyle(kOpenCircle);
  gbr->Draw("pa");
  gbr->GetYaxis()->SetTitle("R_{3}");
  gbr->GetXaxis()->SetTitle("Chi2NDF Cuts");
  TGraphErrors* gy3berr = new TGraphErrors(6, cuts+1,y3brelerr+1, 0, 0);
  gy3berr->SetMarkerStyle(kOpenCircle);
  gy3berr->Draw("pa l");
  gy3berr->GetYaxis()->SetTitle("Relative uncertainty");
  gy3berr->GetXaxis()->SetTitle("chi2Topo Cuts");

}
void SysCuts()
{
  SetsPhenixStyle();
  int const nsys = 14;
  double p[nsys+1],perr[nsys+1];
  for (int i=0;i<nsys;i++)
  { 
    // double topocut = 2.;
    // if (i==4) topocut=1.5;
    // else if (i==5) topocut=2.5;
    // drawMixDataScanSys(i, Form("hH3LMassPtYSysCut%d",i ), Form("h3H3L_chi2topoSysCut%d",i), Form("syscuts_0050_topocut%d.pdf", i),  topocut, p[i], perr[i]); 
    double ndfcut = 3.5;
    if (i==6) ndfcut=3;
    else if (i==7) ndfcut=4;
    // drawMixDataScanSys(i, Form("hH3LMassPtYSysCut%d",i ), Form("h3H3L_chi2ndfSysCut%d",i), Form("syscuts_0050_ndfcut%d.pdf", i),  ndfcut, p[i], perr[i]); 
  }
  //  drawMixDataScanSys(6, Form("hH3LMassPtYSysCut%d",6 ), Form("h3H3L_chi2topoSysCut%d",6), Form("syscuts_0050_topocut%d.pdf", 6), 2, p[6], perr[6]); 
  drawMixDataScanSys(nsys, "hH3LMassPtY", "h3H3L_chi2topo", Form("syscuts_0050_topocut%d.pdf", 14),  2., p[nsys], perr[nsys]); 
  return;
  // drawMixDataScanSys(nsys, "hH3LMassPtY", "h3H3L_chi2ndf", Form("syscuts_0050_ndfcut%d.pdf", 14),  3.5, p[nsys], perr[nsys]); 
  double cuts[20]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
  // TGraphErrors* gBr = new TGraphErrors( nsys+1, cuts, p, 0, perr );
  // gBr->Draw("pA");
  // gBr->SaveAs("SysBr.root");
  //
  // drawMixDataScanSys(4, Form("hH3LMassPtYSysCut%d",4 ), Form("h3H3L_chi2topoSysCut%d",4), Form("syscuts_0050_topocut%d.pdf", i),  topocut, p[i], perr[i]); 
  TString cutsname[nsys]={"l>6","l>10","ldl>4.","ldl>6","chi2topo<1.5","chi2topo<2.5","chi2ndf<3","chi2ndf<4","pi chi2prim>8","pi chi2prim>12","p chi2prim>4","p chi2prim>6","nHitsFit>20","nHitsFit>17"};
  TFile* f[nsys];
  TH1F* h[nsys]; 
  TH1F* h0[nsys]; 
  TH1F* h3[nsys]; 
  // f[5] = new TFile( Form("fout_0050_sys_%d.root", 10));
  // h[5] = (TH1F*)f[5]->Get("hPurityCor_1")->Clone(Form("htest%d",5));
  //
  // TH1F* hRatio[nsys];
  // TLegend* leg = new TLegend(0.2,0.2,0.5,0.5);
  // for (int i=0;i<nsys;i++)
  // {
  //   f[i] = new TFile( Form("fout_0050_sys_%d.root", i));
  //   h[i] = (TH1F*)f[i]->Get("hPurityCor_1")->Clone(Form("htest%d",i));
  //   h[i]->SetDirectory(0);
  //   h[i]->SetLineColor(i);
  //   h[i]->SetMarkerColor(i);
  //   if (i>1) h[i]->Draw("same");
  //   else h[i]->Draw();
  //   f[i]->Close();
  //   // leg->AddEntry(h[i], Form("chi2topo<%0.1f", cuts[i]), "pe");
  //   hRatio[i] = (TH1F*)h[i]->Clone(Form("hr%d",i));
  //   hRatio[i]->Divide(h[5]);
  //   if (i==1) hRatio[i]->Draw();
  //   hRatio[i]->Draw("same");
  //   hRatio[i]->GetYaxis()->SetTitle("Yield of CurrentCuts/Yield of Chi2Topo<3");
  //   leg->AddEntry(hRatio[i], Form("%s", cutsname[i].Data()), "pe");
  // }
  // leg->Draw();

  TFile* f2b = new TFile("h_h3l_corr_yield.root_cut00002_cent0_yuehang.root");
  TGraphErrors* g2b[2];
  g2b[0] = (TGraphErrors*)f2b->Get("t_sgct1_corr_yield[0][0]");
  g2b[1] = (TGraphErrors*)f2b->Get("t_sgct1_corr_yield[0][1]");
  f2b->Close();
  double y2b=0, y2berr=0;
  double binwidth1[3]={ 0.4, 0.4, 0.4};
  double binwidth0[3]={ 0.4, 0.6, 0.4};
  for (int i=0;i<3;i++){
    double x, y,err;
    g2b[1]->GetPoint( i, x, y);
    err=g2b[1]->GetErrorY(i)*x*binwidth1[i]*2*3.1415*10*0.25;
    y2b+=y*x*binwidth1[i]*2*3.1415*10*0.25;
    y2berr+=err*err;
  }
  for (int i=0;i<2;i++){
    double x, y,err;
    g2b[0]->GetPoint( i, x, y);
    err=g2b[0]->GetErrorY(i)*x*binwidth0[i]*2*3.1415*0.25;
    y2b+=y*x*binwidth0[i]*2*3.1415*0.25;
    y2berr+=err*err;
  }
  y2berr=sqrt(y2berr);
  cout<< "2 body yield: " << y2b<<" "<< y2berr<< " relative error " << y2berr/y2b<< endl;

  double br[nsys+1], brerr[nsys+1], y3brelerr[nsys+1];
  int color[nsys+1]={ kGray+2, kMagenta, kRed, kGreen+2, kBlue, kOrange+5, kBlack};
  TLegend* leg = new TLegend(0.2,0.2,0.5,0.5);
  double yield[nsys+1], error[nsys+1];
  for (int i=0;i<nsys+1;i++)
  {
    f[i] = new TFile( Form("fout_0050_sys_%d.root", i));
    h[i] = (TH1F*)f[i]->Get("hPurityCor_1")->Clone(Form("htest%d",i));
    h0[i] = (TH1F*)f[i]->Get("hPurityCor_0")->Clone(Form("htest2%d",i));
    // h3[i] = (TH1F*)f[i]->Get("hPurityCor_2")->Clone(Form("htest2%d",i));
    h3[i] = (TH1F*)f[i]->Get("hYieldCor_2")->Clone(Form("htest2%d",i));
    // h0[i] = (TH1F*)f[i]->Get("hYieldCor_0")->Clone(Form("htest2%d",i));
    // h[i] = (TH1F*)f[i]->Get("hYieldCor_1")->Clone(Form("htest%d",i));
    h[i]->SetDirectory(0);
    h0[i]->SetDirectory(0);
    h3[i]->SetDirectory(0);
    h[i]->SetLineColor(color[i]);
    h[i]->SetMarkerColor(color[i]);
    if (i>1) h[i]->Draw("same");
    else h[i]->Draw();
    f[i]->Close();
    leg->AddEntry(h[i], Form("chi2topo<%0.1f", cuts[i]), "pe");
    // hRatio[i] = (TH1F*)h[i]->Clone(Form("hr%d",i));
    // hRatio[i]->Divide(h[6]);
    // if (i==1) hRatio[i]->Draw();
    // hRatio[i]->Draw("same");
    // hRatio[i]->GetYaxis()->SetTitle("Yield of CurrentCuts/Yield of Chi2Topo<3");
    // leg->AddEntry(hRatio[i], Form("chi2topo<%0.1f", cuts[i]), "pe");
    yield[i]=0;error[i]=0;
    for (int ib=1;ib<=3;ib++)
    {
      yield[i]+=h[i]->GetBinContent(ib)*h[i]->GetBinWidth(ib)*h[i]->GetBinCenter(ib)*2.*3.1415*0.25;
      error[i]+=pow(h[i]->GetBinError(ib)*h[i]->GetBinWidth(ib)*h[i]->GetBinCenter(ib)*2.*3.1415*0.25, 2);
    }
    for (int ib=1;ib<=2;ib++)
    {
      if (h0[i]->GetBinContent(ib)<1e-31) continue;
      yield[i]+=h0[i]->GetBinContent(ib)*h0[i]->GetBinWidth(ib)*h0[i]->GetBinCenter(ib)*2.*3.1415*0.25;
      error[i]+=pow(h0[i]->GetBinError(ib)*h0[i]->GetBinWidth(ib)*h0[i]->GetBinCenter(ib)*2.*3.1415*0.25,2);
    }
    // for (int ib=1;ib<=3;ib++)
    // {
    //   // if (h3[i]->GetBinContent(ib)<1e-31) continue;
    //   yield[i]+=h3[i]->GetBinContent(ib)*h3[i]->GetBinWidth(ib)*h3[i]->GetBinCenter(ib)*2.*3.1415;
    //   error[i]+=pow(h3[i]->GetBinError(ib)*h3[i]->GetBinWidth(ib)*h3[i]->GetBinCenter(ib)*2.*3.1415,2);
    // }

    error[i]=sqrt(error[i]);
    y3brelerr[i]=error[i]/yield[i];
    br[i] = y2b/(yield[i]+y2b);
    brerr[i] = sqrt( pow((y2b/(yield[i]+y2b)/(yield[i]+y2b))*error[i], 2) +  pow((yield[i]/(yield[i]+y2b)/(yield[i]+y2b))*y2berr, 2));
    cout<< "y3berr: " << y3brelerr[i]<<" y2berr: "<< y2berr/y2b<<" Brerr: "<< brerr[i]/br[i]<< " "<< br[i]<<" "<<brerr[i]<<endl;
  }
  
  leg->Draw();
  TGraphErrors* gyield = new TGraphErrors(nsys+1, cuts,yield, 0, error);
  TGraphErrors* gyield_d = new TGraphErrors(1, cuts+nsys,yield+nsys, 0, error+nsys);
  gyield->SetMarkerStyle(kOpenCircle);
  gyield->Draw("pa");
  gyield->GetYaxis()->SetMaxDigits(3);
  gyield->GetXaxis()->SetNdivisions(520);
  gyield_d->SetMarkerColor(kBlue);
  gyield_d->SetLineColor(kBlue);
  gyield_d->SetMarkerStyle(kFullCircle);
  gyield_d->Draw("p same");
  gyield->GetYaxis()->SetTitle("Yield");
  gyield->GetXaxis()->SetTitle("Sys. Cuts");
  
  double sys=0,n=0;
  for (int ip=0;ip<gyield->GetN()-1;ip++)
  {
     double delta = fabs(yield[ip] - yield[14]);
     double ds = sqrt( fabs(error[ip]*error[ip]-error[14]*error[14]));
     // if (ds>delta) cout << "pass check "<< ip+1<< endl;
     // else {
       sys+=delta*delta/2.;
      // else if (ip==12) sys+=delta*delta;
     // }
  }
  cout <<"systematics of yields: "<< sqrt(sys)/yield[14] <<endl;
  cout <<"systematics of R3: "<<sqrt( pow((y2b/(yield[14]+y2b)/(yield[14]+y2b))*sqrt(sys), 2))<<" "<< y2b/(y2b+yield[14])<<" "<<brerr[14]  <<endl;
  drawLine(0, yield[14], 16, yield[14], 2, 1,2);
  drawLine(0, yield[14]+sqrt(sys), 16, yield[14]+sqrt(sys), 2, 2,2);
  drawLine(0, yield[14]-sqrt(sys), 16, yield[14]-sqrt(sys), 2, 2,2);
  gyield_d->Draw("p same");

  TGraphErrors* gbr = new TGraphErrors(nsys+1, cuts,br, 0, brerr);
  gbr->SetMarkerStyle(kOpenCircle);
  gbr->Draw("pa");
  gbr->GetYaxis()->SetTitle("R_{3}");
  gbr->GetXaxis()->SetTitle("Sys. Cuts");
  return;

  TGraphErrors* gy3berr = new TGraphErrors(nsys+1, cuts,y3brelerr, 0, 0);
  gy3berr->SetMarkerStyle(kOpenCircle);
  gy3berr->Draw("pa");
  gy3berr->GetYaxis()->SetTitle("Relative error of 3 body yield in -0.5<y<-0.25");
  gy3berr->GetXaxis()->SetTitle("Sys. Cuts");

}

void drawData()
{
  // drawCompLaH3L();
  // drawRT_check();
  // drawMixData();
  // drawMixData_015pt();
  // drawMixDataMB();
  // checkbkgd();
  // drawMixDataTest();
  // drawSEData();
  // checkbkgd();
  // correctEff();

  SysCuts(); 
  // scanCuts(); 
  // scanCutsNDF(); 
  // drawComp();
  // drawCompQuasi();
  // drawR3();
}
