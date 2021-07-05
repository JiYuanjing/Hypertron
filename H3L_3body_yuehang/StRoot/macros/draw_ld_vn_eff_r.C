TH3F *h_mass_pt_fine;
TH2F *g_pt_fine;
TH3F *h_mass_pt_fine_wgt;
TH2F *g_pt_fine_wgt;
TH3F *h_mass_pt_fine2_wgt;
TH2F *g_pt_fine2_wgt;
int _cut_mode;
void draw_ld_vn_eff_r(int cut_mode){

  _cut_mode = cut_mode;

  gStyle->SetTitleFont(62,"X");
  gStyle->SetTitleFont(62,"Y");
  gStyle->SetLabelFont(62,"X");
  gStyle->SetLabelFont(62,"Y");
  gStyle->SetTextFont(62);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetStatStyle(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetFrameBorderSize(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetStatBorderSize(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat(0);

  gStyle->SetLabelSize(0.0575,"X");
  gStyle->SetTitleSize(0.0575,"X");
  gStyle->SetLabelSize(0.0575,"Y");
  gStyle->SetTitleSize(0.0575,"Y");



TFile *f0;
//f0 = new TFile("run18_3gev_ld_hist_mu_1x_vn_newreweight_centa0.root");
f0 = new TFile(Form("run18_3gev_ld_hist_mu_1x_vn_newreweight_centa0_cut%d.root",_cut_mode),"READ");

h_mass_pt_fine = (TH3F*)f0->Get("h_mass_pt_fine");
g_pt_fine = (TH2F*)f0->Get("g_pt_fine");

h_mass_pt_fine_wgt = (TH3F*)f0->Get("h_mass_pt_fine_wgt");
g_pt_fine_wgt = (TH2F*)f0->Get("g_pt_fine_wgt");

h_mass_pt_fine2_wgt = (TH3F*)f0->Get("h_mass_pt_fine2_wgt");
g_pt_fine2_wgt = (TH2F*)f0->Get("g_pt_fine2_wgt");

TH2F *g_pt_fine2_wgt_r;
g_pt_fine2_wgt_r = (TH2F*)g_pt_fine2_wgt->Clone("g_pt_fine2_wgt_r");

for(int ix=1;ix<=g_pt_fine2_wgt_r->GetNbinsX();ix++ ){//rap
    for(int iy=1;iy<=g_pt_fine2_wgt_r->GetNbinsY();iy++ ){//pt

g_pt_fine2_wgt_r->SetBinContent(g_pt_fine2_wgt_r->GetNbinsX()-ix+1,iy,g_pt_fine2_wgt->GetBinContent(ix,iy));
g_pt_fine2_wgt_r->SetBinError(g_pt_fine2_wgt_r->GetNbinsX()-ix+1,iy,g_pt_fine2_wgt->GetBinError(ix,iy));

	}
}


TH1F *h_mass_pt_wgt_px;
TH2F *h_pt_wgt_fine;
TH2F *h_pt_wgt_fine2;

TH1F *h_mass_pt_px;
TH2F *h_pt_fine;
TH2F *h_pt_fine2;

 h_pt_fine = new TH2F("h_pt_fine","",100,-1.0,1.0,50,0,5);
 h_pt_wgt_fine = new TH2F("h_pt_wgt_fine","",100,-1.0,1.0,50,0,5);
 h_pt_wgt_fine2 = new TH2F("h_pt_wgt_fine2","",20,-1.0,1.0,50,0,5);

double bincenter;
double bincontent;
double binerror;
double sgm3 = 0.0018;
double lambda_mass = 1.115683;
double counts = 0.;
double countserr = 0.;


for(int ix=1;ix<=h_pt_fine->GetNbinsX();ix++ ){//rap
    for(int iy=1;iy<=h_pt_fine->GetNbinsY();iy++ ){//pt
      h_mass_pt_px = (TH1F*)h_mass_pt_fine->ProjectionX("h_mass_pt_px",ix,ix,iy,iy);
      counts = 0;
      for(int ibin=1;ibin<=h_mass_pt_px->GetNbinsX();ibin++){
        bincenter = h_mass_pt_px->GetBinCenter(ibin);
        bincontent = h_mass_pt_px->GetBinContent(ibin);
        if( bincenter>(lambda_mass-3*sgm3) && bincenter<(lambda_mass+3*sgm3) )  {counts+=bincontent;}
      }
    h_pt_fine->SetBinContent(ix,iy,counts);
    h_pt_fine->SetBinError(ix,iy,sqrt(counts));
  }
}

TCanvas *yrrr3 = new TCanvas("yrrr3","yrrr3",800,800);
yrrr3->cd();
yrrr3->cd()->SetRightMargin(0.12);
yrrr3->cd()->SetTopMargin(0.02);
yrrr3->cd()->SetLeftMargin(0.12);
yrrr3->cd()->SetBottomMargin(0.14);

g_pt_fine->Draw("colz");
yrrr3->cd();
h_pt_fine->Draw("colz");

TH2F *h_pt_eff_fine;
h_pt_eff_fine = (TH2F*)h_pt_fine->Clone();
h_pt_eff_fine -> Divide(g_pt_fine);
h_pt_eff_fine -> Draw("colz");


counts = 0.;
countserr = 0.;

for(int ix=1;ix<=h_pt_wgt_fine->GetNbinsX();ix++ ){//rap
    for(int iy=1;iy<=h_pt_wgt_fine->GetNbinsY();iy++ ){//pt
      h_mass_pt_wgt_px = (TH1F*)h_mass_pt_fine_wgt->ProjectionX("h_mass_pt_wgt_px",ix,ix,iy,iy);
      counts = 0;
      for(int ibin=1;ibin<=h_mass_pt_wgt_px->GetNbinsX();ibin++){
        bincenter = h_mass_pt_wgt_px->GetBinCenter(ibin);
        bincontent = h_mass_pt_wgt_px->GetBinContent(ibin);
        if( bincenter>(lambda_mass-3*sgm3) && bincenter<(lambda_mass+3*sgm3) )  {counts+=bincontent;}
      }
    h_pt_wgt_fine->SetBinContent(h_pt_wgt_fine->GetNbinsX()-ix+1,iy,counts);
    h_pt_wgt_fine->SetBinError(h_pt_wgt_fine->GetNbinsX()-ix+1,iy,sqrt(counts));
  }
}


yrrr3->cd();
g_pt_fine_wgt->Draw("colz");
yrrr3->cd();
h_pt_wgt_fine->Draw("colz");

TH2F *h_pt_eff_wgt_fine;
h_pt_eff_wgt_fine = (TH2F*)h_pt_wgt_fine->Clone();
h_pt_eff_wgt_fine -> Divide(g_pt_fine_wgt);
h_pt_eff_wgt_fine -> Draw("colz");



counts = 0.;
countserr = 0.;

for(int ix=1;ix<=h_pt_wgt_fine2->GetNbinsX();ix++ ){//rap
    for(int iy=1;iy<=h_pt_wgt_fine2->GetNbinsY();iy++ ){//pt
      h_mass_pt_wgt_px = (TH1F*)h_mass_pt_fine2_wgt->ProjectionX("h_mass_pt_wgt_px",ix,ix,iy,iy);
      counts = 0;
      for(int ibin=1;ibin<=h_mass_pt_wgt_px->GetNbinsX();ibin++){
        bincenter = h_mass_pt_wgt_px->GetBinCenter(ibin);
        bincontent = h_mass_pt_wgt_px->GetBinContent(ibin);
        if( bincenter>(lambda_mass-3*sgm3) && bincenter<(lambda_mass+3*sgm3) )  {counts+=bincontent;}
      }
    h_pt_wgt_fine2->SetBinContent(h_pt_wgt_fine2->GetNbinsX()-ix+1,iy,counts);
    h_pt_wgt_fine2->SetBinError(h_pt_wgt_fine2->GetNbinsX()-ix+1,iy,sqrt(counts));
  }
}


yrrr3->cd();
g_pt_fine2_wgt->Draw("colz");
yrrr3->cd();
h_pt_wgt_fine2->Draw("colz");

TH2F *h_pt_eff_wgt_fine2;
h_pt_eff_wgt_fine2 = (TH2F*)h_pt_wgt_fine2->Clone("h_pt_eff_wgt_fine2");
//h_pt_eff_wgt_fine2 -> Divide(g_pt_fine2_wgt);
h_pt_eff_wgt_fine2 -> Divide(g_pt_fine2_wgt_r);

h_pt_eff_wgt_fine2->GetYaxis()->SetRangeUser(0,3);
h_pt_eff_wgt_fine2->GetXaxis()->SetRangeUser(-1.2,1.2);

h_pt_eff_wgt_fine2->GetXaxis()->SetTitle("y");
h_pt_eff_wgt_fine2->GetYaxis()->SetTitle("p_{T}[GeV/c]");

h_pt_eff_wgt_fine2 -> Draw("colz");

yrrr3->Print("h_pt_eff_wgt_fine2.pdf","pdf");

TFile *outhistfile = new TFile (Form("ld_vn_eff_v0_cut%d_r.root",_cut_mode), "RECREATE");
outhistfile->cd();
h_pt_eff_wgt_fine2 -> Write();
outhistfile->Close();
delete outhistfile;


}
