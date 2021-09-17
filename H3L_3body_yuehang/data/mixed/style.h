#ifndef _style_h
#define _style_h
#include "TH1D.h"
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
void drawline(double x1, double y1, double x2, double y2, int color, int s=0, int width=1 )
{
  TLine* l = new  TLine(x1,y1,x2,y2);
  l->SetLineStyle(s);
  l->SetLineColor(color);
  l->SetLineWidth(width);
  l->Draw("same");
}
void drawerrbar(double x, double y, double ysys, double xwidth, double shortlinesize, int color){
  double xw = xwidth; //x-width or x sys
  // double yl = shortlinesize*y;
  double yl = shortlinesize;
  int s = 1;
  int w = 2;
  //top vertical line
  drawline(x-xw,y+ysys, x+xw, y+ysys, color,s,w);
  //bottom vertical line
  drawline(x-xw,y-ysys, x+xw, y-ysys, color,s,w);
  w = 1;
  //top left line
  drawline(x-xw,y+ysys, x-xw, y+ysys-yl, color,s,w);
  //top right line
  drawline(x+xw,y+ysys, x+xw, y+ysys-yl, color,s,w);
  //bottom left line
  drawline(x-xw,y-ysys, x-xw, y-ysys+yl, color,s,w);
  //bottom right line
  drawline(x+xw,y-ysys, x+xw, y-ysys+yl, color,s,w);
}
void drawerrbox(double x, double y, double ysys, double xwidth,  int color,int fillstyle,int mode=0){
  double xw = xwidth; //x-width or x sys
  TBox* box = new TBox(x-xw*1.1 ,y-ysys,x+xw*1.1 ,y+ysys );
  // box->SetFillStyle(0);
  box->SetFillStyle(fillstyle);
  // box->SetFillColor(0);
  // box->SetLineColor(color);
  // box->SetLineWidth(2);
  if (mode==0) box->SetFillColorAlpha(color,0.3);
  else if (mode==1) box->SetFillColor(color);
  else cout <<"please set error box fill mode" <<endl;
  box->Draw();
}

void drawGraphWithSys(TString fname,TString gname , TString gnamesys, int color,int style,float size)
{
  TFile* f = TFile::Open(fname.Data());
  TGraphErrors* g = (TGraphErrors*)f->Get(gname.Data());
  TGraphErrors* gsys = (TGraphErrors*)f->Get(gnamesys.Data());
  cout<<"start draw..." <<endl;
  g->SetMarkerStyle(style);
  g->SetMarkerColor(color);
  g->SetMarkerSize(size);
  g->SetLineColor(color);
  g->SetLineWidth(2);
  g->Draw("psame");
  int npoints=g->GetN();
  for (int ip=0;ip<npoints;ip++)
  {
    double x,y,yerr;
    gsys->GetPoint(ip,x,y);
    yerr = gsys->GetErrorY(ip);
    drawerrbar(x,y,yerr,0.03,0.005,color); 

  }
}
void drawGraphWithSys(TGraphErrors* g,TGraphErrors* gsys , int color,int style,float size,int mod=0,int fillcolor=0)
{
  // TFile* f = TFile::Open(fname.Data());
  // TGraphErrors* g = (TGraphErrors*)f->Get(gname.Data());
  // TGraphErrors* gsys = (TGraphErrors*)f->Get(gnamesys.Data());
  cout<<"start draw..." <<g->GetName()<<endl;
  g->SetMarkerStyle(style);
  g->SetMarkerColor(color);
  g->SetLineColor(color);
  g->SetMarkerSize(size);
  g->SetLineWidth(2.0);
  g->Draw("psame");
  int npoints=g->GetN();
  double maxX, maxY;
  g->GetPoint(npoints-1,maxX,maxY);
  for (int ip=0;ip<npoints;ip++)
  {
    double x,y,yerr;
    gsys->GetPoint(ip,x,y);
    yerr = gsys->GetErrorY(ip);
    // if (mod==0) drawerrbar(x,y,yerr,0.015*maxX,0.002,color); 
    if (mod==0) drawerrbar(x,y,yerr,0.022,0.005,color); 
    else if (mod==1) drawerrbox(x,y,yerr,0.023,color,1001);
    else if (mod==2) drawerrbox(x,y,yerr,0.022,fillcolor,1001);
    else cout<<"no error mode set!!!"<<endl;
  }
  g->Draw("psame");
}
void drawSTAR(double x,double y,double size,int font, int color)
{
  TLatex lat;
  // lat.SetTextSize(0.05);
  lat.SetTextSize(size);
  // lat.SetTextFont(72);
  lat.SetTextFont(font);
  // lat.SetTextColor(kRed);
  lat.SetTextColor(color);
  lat.DrawLatexNDC ( x, y, "STAR Preliminary");
}
void setPadStyle(TCanvas* c1)
{
  c1->SetFillColor(10);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderMode(0);

  c1->SetTopMargin(0.02);
  c1->SetBottomMargin(0.15);
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.02);
}
void setPadStyle(TPad* pad,double t,double b,double l,double r)
{
  pad->SetFillColor(10);
  pad->SetBorderMode(0);
  pad->SetBorderSize(2);
  pad->SetFrameFillColor(0);
  pad->SetFrameBorderMode(0);

  pad->SetLeftMargin(l);
  pad->SetBottomMargin(b);
  pad->SetTopMargin(t);
  pad->SetRightMargin(r);
}
void setPad(Double_t left, Double_t right, Double_t top, Double_t bottom, int color=10)
{
    gPad->SetFillColor(color);
    // gPad->SetBorderMode(0);
    // gPad->SetBorderSize(0);
    gPad->SetFrameFillColor(10);
    gPad->SetFrameBorderMode(0);
    // gPad->SetFrameBorderSize(0);
    gPad->SetLeftMargin(left);
    gPad->SetRightMargin(right);
    gPad->SetTopMargin(top);
    gPad->SetBottomMargin(bottom);
}
void setLegendStyle(TLegend* leg,int font,double size)
{
   leg->SetFillColor(10);
   leg->SetFillStyle(10);
   leg->SetLineStyle(4000);
   leg->SetLineColor(10);
   leg->SetLineWidth(0.);
   leg->SetBorderSize(0.);
   leg->SetTextFont(font);
   leg->SetTextSize(size);
}
void drawLatexAxis(double x,double y,const char* txt,double size,int font=42,int col=1)
{
  TLatex lat;
  lat.SetTextSize(size);
  lat.SetTextFont(font);
  lat.SetTextColor(col);
  lat.DrawLatex( x, y, txt);
}
void drawLatex(double x,double y,const char* txt,double size,int font=42, int col=1)
{
  TLatex lat;
  lat.SetTextSize(size);
  lat.SetTextFont(font);
  lat.SetTextColor(col);
  lat.DrawLatexNDC( x, y, txt);
}
TH1D* histo(TString histname, Double_t xlow, Double_t xup, Double_t ylow, Double_t yup, TString xTitle, TString yTitle){
  TH1D *d0 = new TH1D(histname,"",1,xlow,xup);
  d0->SetMinimum(ylow);
  d0->SetMaximum(yup);
  d0->GetXaxis()->SetNdivisions(406);
  d0->GetXaxis()->SetTitle(xTitle);
  d0->GetXaxis()->SetTitleOffset(1.0);
  d0->GetXaxis()->SetTitleSize(0.07);
  d0->GetXaxis()->SetLabelOffset(0.015);
  d0->GetXaxis()->SetLabelSize(0.055);
  d0->GetXaxis()->SetLabelFont(42);
  d0->GetXaxis()->SetTitleFont(42);
  d0->GetYaxis()->SetNdivisions(406);
  d0->GetYaxis()->SetTitle(yTitle);
  d0->GetYaxis()->SetTitleOffset(1.);
  d0->GetYaxis()->SetTitleSize(0.07);
  d0->GetYaxis()->SetLabelOffset(0.01);
  d0->GetYaxis()->SetLabelSize(0.055);
  d0->GetYaxis()->SetLabelFont(42);
  d0->GetYaxis()->SetTitleFont(42);
  d0->SetLineWidth(2);

	return d0;
}
void drawPad(double x1,double x2,double y1,double y2,TString xTitle,TString yTitle,TString histname="d0")
{
  TH1D *d0 = (TH1D*)histo(histname,x1,x2,y1,y2,xTitle,yTitle);

  d0->DrawCopy("c");

  drawline(x1,y1,x2,y1,1,1,2);
  drawline(x1,y2,x2,y2,1,1,2);
  drawline(x1,y1,x1,y2,1,1,2);
  drawline(x2,y1,x2,y2,1,1,2);
  delete d0;
}

void setTheoryStyleAs(TGraphAsymmErrors* gr, int style, int color, double alpha)
{
    gr->SetMarkerStyle(0);
    gr->SetFillStyle(style);
    gr->SetLineWidth(0);
    gr->SetLineColor(color);
    gr->SetFillColor(color);
    gr->SetFillColorAlpha(color, alpha);
}
void setTheoryStyle(TGraphErrors*gr, int style, int color, double alpha)
{
    gr->SetMarkerStyle(0);
    gr->SetFillStyle(style);
    gr->SetLineWidth(0);
    gr->SetLineColor(color);
    gr->SetFillColor(color);
    gr->SetFillColorAlpha(color, alpha);
}
//-----------------------
void drawTheorySideLineAs(TGraphAsymmErrors* gr, int color)
{
    int np = gr->GetN();
    double *x = new double[np];
    double *y = new double[np];
    double *yL = new double[np];
    double *yH = new double[np];
    for(int i=0; i<np; i++) {
        double eyH, eyL;
        gr->GetPoint(i, x[i], y[i]);
        eyH = gr->GetErrorYhigh(i);
        eyL = gr->GetErrorYlow(i);
        yL[i] = y[i] - eyL;
        yH[i] = y[i] + eyH;
    }
    TGraph *gL = new TGraph(np, x, yL);
    gL->SetLineColor(color);
    gL->Draw("sameL");
    TGraph *gH = new TGraph(np, x, yH);
    gH->SetLineColor(color);
    gH->Draw("sameL");
}
void drawTheorySideLine(TGraphErrors *gr, int color)
{
    int np = gr->GetN();
    double *x = new double[np];
    double *y = new double[np];
    double *yL = new double[np];
    double *yH = new double[np];
    for(int i=0; i<np; i++) {
        double eyH, eyL;
        gr->GetPoint(i, x[i], y[i]);
        eyH = gr->GetErrorY(i);
        eyL = gr->GetErrorY(i);
        yL[i] = y[i] - eyL;
        yH[i] = y[i] + eyH;
    }
    TGraph *gL = new TGraph(np, x, yL);
    gL->SetLineColor(color);
    gL->Draw("sameL");
    TGraph *gH = new TGraph(np, x, yH);
    gH->SetLineColor(color);
    gH->Draw("sameL");
}
//----------------------
void drawTheory(TGraphErrors *gr, int style, int color, double alpha)
{
   setTheoryStyle(gr, style, color, alpha);
   gr->Draw("same3");
   gr->SetLineWidth(1.5);
   drawTheorySideLine(gr,color);
}
void drawTheoryAs(TGraphAsymmErrors *gr, int style, int color, double alpha)
{
   setTheoryStyleAs(gr, style, color, alpha);
   gr->Draw("same E3");
   gr->SetLineWidth(1.5);
   drawTheorySideLineAs(gr,color);
}
//---------------------
void drawPoint2Band(TGraphErrors* gr, TGraphErrors* gsys, int color, int style)
{
  int const np = gr->GetN();
  double x[50], y[50], err[50], errysys[50], errtot[50];
  int ndraw=0;
  for (int i=0; i<np; i++)
  {
     gr->GetPoint(i, x[i], y[i]);
     if (x[i]>2) break; 
     err[i] = gr->GetErrorY(i);
     errysys[i] = gr->GetErrorY(i);
     errtot[i] = sqrt(errysys[i]*errysys[i]+err[i]*err[i]);
     ndraw++;
  }
  TGraphErrors* gtot = new TGraphErrors(ndraw, x, y, 0, errtot); 
  gtot->SetFillColorAlpha(color,0.6);
  gtot->SetFillStyle(style);
  gtot->Draw("same3");
  drawTheorySideLine(gtot, color);
  //for consistent
  gr->SetFillColorAlpha(color,0.6);
  gr->SetFillStyle(style);

}
void drawNonflow(TGraphErrors* gdata,TGraph* gnonflow,double wd=0.05)
{
   int np = gdata->GetN();
   double xl,xr,yh,yl,x;
   for (int ip=0;ip<np;ip++)
   {
      gdata->GetPoint(ip,x,yh);
      if (x<0.5) continue;
      yl = yh-gnonflow->Eval(x);
      TBox* box = new TBox(x-0.05,yh,x+0.05,yl);
      box->SetLineColor(kGray);
      box->SetLineWidth(kGray);
      box->SetFillStyle(1001);
      box->SetFillColorAlpha(kGray+1,0.5);
      box->Draw("same");
   }
}
void shiftGraph(TGraph* g, double shift)
{
  for (int i=0;i<g->GetN();i++){
    double x, y;
    g->GetPoint(i, x, y);
    g->SetPoint(i, x+shift, y);
  }
}
void shiftGraph(TGraph* g, double shift, int ip)
{
    double x, y;
    g->GetPoint(ip, x, y);
    g->SetPoint(ip, x+shift, y);
}
void shiftGraphLogx(TGraph* g, double shift)
{
  for (int i=0;i<g->GetN();i++){
    double x, y;
    g->GetPoint(i, x, y);
    g->SetPoint(i, x*pow(10, shift), y);
  }
}
#endif
