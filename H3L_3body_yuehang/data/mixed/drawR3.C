
#include "style.h"

void drawR3()
{
  gROOT->Reset();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0.01);

  TCanvas* c = new TCanvas("c1","c1",1200,600);
  gStyle->SetLineWidth(2);
  setPadStyle(c);
  c->cd();

  double num[7]={1,2,3,4,4,5,6};
  double val[7]={0.39, 0.36, 0.39, 0.46, 0.41, 0.3, 0.35};
  double up[7]={0.07, 0.08, 0.12, 0.04, 0.04, 0.07, 0.04};
  double low[7]={0.07, 0.06, 0.07, 0.03, 0.03, 0.07, 0.04};
  double ii[1]={6};
  double star[1]={0.32};
  double stsys[1]={0.08};
  double ststat[1]={0.05};

  double nn[1]={7};
  double mean[1]={0.296};
  double sys[1]={0.023};
  double stat[1]={0.029};

// TGraphAsymmErrors* gworld = new TGraphAsymmErrors( 7, num, val, 0, 0, low, up);
// gworld->Draw("AP");
  double x1=0, x2=10, y1=0.15, y2=0.65;
  TH1F* h = new TH1F("h", ";; ^{3}_{#Lambda}H R_{3}", 1, x1, x2);
  h->SetMinimum(y1);
  h->SetMaximum(y2);
  h->GetXaxis()->SetNdivisions(406);
  /* h->GetXaxis()->SetTitle(xtitle); */
  h->GetXaxis()->SetTitleOffset(10.);
  h->GetXaxis()->SetTitleSize(0.06);
  // h->GetXaxis()->SetLabelOffset(0.006);
  h->GetXaxis()->SetLabelSize(0.0);
  h->GetXaxis()->SetLabelFont(42);
  h->GetXaxis()->SetTitleFont(42);
  h->GetYaxis()->SetNdivisions(406);
  /* h->GetYaxis()->SetTitle(ytitle); */
  h->GetYaxis()->SetTitleOffset(1.04);
  h->GetYaxis()->SetTitleSize(0.07);
  // h->GetYaxis()->SetLabelOffset(0.01);
  h->GetYaxis()->SetLabelSize(0.055);
  h->GetYaxis()->SetLabelFont(42);
  h->GetYaxis()->SetTitleFont(42);
  h->SetLineWidth(2);
  h->DrawCopy();
  h->Draw("c");
  
  drawline(x1,y1,x2,y1,1,1,2);
  drawline(x1,y2,x2,y2,1,1,2);
  drawline(x1,y1,x1,y2,1,1,2);
  drawline(x2,y1,x2,y2,1,1,2);

  h->GetXaxis()->SetTickSize(0);

  TGraphAsymmErrors* gworld[7];
  int color[]={kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack};
  int marker[]={kOpenCircle, kOpenStar,kOpenSquare, kOpenTriangleUp, kOpenTriangleDown, kOpenDiamond, kFullCircle };
  double size[]={ 1.5, 2., 1.5,  2., 2., 2.3, 1.5 };
  for (int i=0;i<6;i++)
  {
    gworld[i] = new TGraphAsymmErrors(1, num+i, val+i, 0, 0, low+i, up+i );
    gworld[i]->Draw("p same");
    setHistStyle(gworld[i], color[i], marker[i], size[i],0 );
  }


  drawLatexAxis(0.7,0.28, "CERN Rep.", 0.035);
  drawLatexAxis(0.7,0.255, "64-1, 63(1964)", 0.035);

  drawLatexAxis(1.7, 0.5, "PRL 20,", 0.035);
  drawLatexAxis(1.7,0.475, "819(1968)", 0.035);

  drawLatexAxis(2.7, 0.57, "NC26,", 0.035);
  drawLatexAxis(2.7, 0.545, "840(1962)", 0.035);

  drawLatexAxis(3.7, 0.35, "NPB16,", 0.035);
  drawLatexAxis(3.7, 0.325, "77(1970)", 0.035);

  drawLatexAxis(4.7, 0.41, "NPB67,", 0.035);
  drawLatexAxis(4.7, 0.385, "269(1973)", 0.035);

  // drawLxtxy(0.62,0.43, "NPB67,", 0.035);
  // drawLxtxy(0.62,0.38, "269(1973)", 0.035);
  // // drxwxatey(0.62,0.65, "mean of helium bubble", 0.035);

  drawLatexAxis(5.6, 0.42, "STAR 2017", 0.035);
  drawLatexAxis(6.6, 0.23, "STAR 2021", 0.035, 42, kRed);

  TGraphAsymmErrors* gstar = new TGraphAsymmErrors(1, ii, star , 0, 0, ststat, ststat);
  gstar->Draw("p same");
  setHistStyle(gstar, kBlack, kFullStar, 2, 0);
  drawerrbar(ii[0], star[0], stsys[0], 0.06, 0.007, kBlack);


  TGraphAsymmErrors* gthis = new TGraphAsymmErrors(1, nn, mean, 0, 0, stat, stat);
  gthis->Draw("p same");
  setHistStyle(gthis, kRed, kFullStar, 2, 0);
  drawerrbar(nn[0], mean[0], sys[0], 0.06, 0.007, kRed);

  drawline( nn[0]+1 -0.2, 0.4, nn[0]+1 +0.2, 0.4, kRed, 2, 4);
  drawLatexAxis( nn[0]+1-0.35, 0.37, "PRC 57", 0.035);
  drawLatexAxis( nn[0]+1-0.35, 0.345, "1595(1998)", 0.035);

}
