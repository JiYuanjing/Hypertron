// #include "style.h"
#include "sPhenixStyle.h"
// #include "draw.h"
// #include "TRandom3.h"
// #include "TCanvas.h"
// #include "TF1.h"
// #include "TH1.h"
// #include "TH2.h"
// #include "TStopwatch.h"
// #include "TGraphErrors.h"
// #include "TFile.h"
// #include "TMath.h"
// #include "TLorentzVector.h"
// #include "TSystem.h"
// #include "TLine.h"
#include "flow.h"
void plotInputs()
{
  SetsPhenixStyle();
  TCanvas* c = new TCanvas("c","c");
  TPDF* pdf = new TPDF("inputspectra.pdf");
  pdf->Off();

  TF1 *levyfit4;
  TF1 *levyfit5;
  TF1 *levyfit6;
  TF1 *levyfit7;
  TF1 *levyfit8;
  TF1 *t_quadr;
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

  levyfit8->Draw();
  levyfit8->GetYaxis()->SetTitle("dN/dp_{T}");
  levyfit8->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  levyfit8->SetLineColor(8);
  levyfit7->Draw("same");
  levyfit7->SetLineColor(7);
  levyfit6->Draw("same");
  levyfit6->SetLineColor(6);
  levyfit5->Draw("same");
  levyfit5->SetLineColor(5);
  levyfit4->Draw("same");
  levyfit4->SetLineColor(4);
  // if (y<-0.7){ ptwt = levyfit8->Eval(pt); }
  //   else if (y<-0.5){ ptwt = levyfit7->Eval(pt); }
  //   else if (y<-0.3){ ptwt = levyfit6->Eval(pt); }
  //   else if (y<-0.1){ ptwt = levyfit5->Eval(pt); }
  //   else if (y<0.1) { ptwt = levyfit4->Eval(pt); }
  //   else if (y<0.3) { ptwt = levyfit5->Eval(pt); }
  //   else if (y<0.5) { ptwt = levyfit6->Eval(pt); }
  //   else if (y<0.7) { ptwt = levyfit7->Eval(pt); }

  gPad->SetLogy();
  TLegend* leg = new TLegend( 0.7,0.7,0.88,0.88);
  leg->AddEntry(levyfit8, "0.7<|y|<1","l");
  leg->AddEntry(levyfit7, "0.5<|y|<0.7","l");
  leg->AddEntry(levyfit6, "0.3<|y|<0.5","l");
  leg->AddEntry(levyfit5, "0.1<|y|<0.3","l");
  leg->AddEntry(levyfit4, "|y|<0.1","l");
  leg->Draw();

  addpdf(pdf);
  t_quadr->Draw();
  t_quadr->GetXaxis()->SetTitle("Rapidity");
  addpdf(pdf);
  return;

  gLdv1->Draw("p e l A");
  gLdv1->SetMarkerSize(2);
  gLdv1->SetMarkerStyle(kOpenCircle);
  gLdv1->GetXaxis()->SetTitle("Rapidity");
  gLdv1->GetYaxis()->SetTitle("Direct flow");
  addpdf(pdf);
  // TGraphErrors* gd_10_40 = new TGraphErrors("xyscan_10_40.txt","%lg %lg");
  //
  //
  //deuteron
  TGraphErrors* gd[10];
  TF1* flevyd = new TF1("flevyd","[2]*pow((1+(sqrt(1.875*1.875+x*x)-1.875)/[1]/[0]),(-1)*[1])", 0.2,4);
  flevyd->SetParameters(0.17, -45, 6);

  double para[10][3];
TLegend* leg3 = new TLegend(0.2,0.15,0.4,0.6);

  for (int i=3;i<13;i++) 
  {
      gd[i-3]=new TGraphErrors(Form("deuteron_dat/d_ptt_1020_y%d.dat",i), "%lg %*lg %*lg %lg %lg");
      gd[i-3]->SetName(Form("gd_pt_1020_y%d",i)); 
      gd[i-3]->SetMarkerColor(i-2);
      if (i==12) gd[i-3]->SetMarkerColor(13);
      gd[i-3]->SetLineColor(i-2);
      if (i==12) gd[i-3]->SetLineColor(13);
      if (i==3) gd[i-3]->Draw("e p l A");
      else if (i>3) gd[i-3]->Draw("p same");
      flevyd->SetLineColor(i-2);
      if (i==12) flevyd->SetLineColor(13);
      gd[i-3]->Fit(flevyd);
      flevyd->GetParameters(para[i-3]);
      leg3->AddEntry(gd[i-3],Form("%0.1f<y<%0.1f", (i-2)*-0.1,-1*(i-3)*0.1),"pl" );
      gd[i-3]->GetYaxis()->SetTitle("dN^2/2#pip_{T}dp_{T}dy");
      gd[i-3]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  }
  leg3->Draw();

  for (int i=0;i<10;i++)
  {
    cout <<"{" <<para[i][0]<<","<<para[i][1]<<","<<para[i][2]<<"},"<< endl;
  }

}

