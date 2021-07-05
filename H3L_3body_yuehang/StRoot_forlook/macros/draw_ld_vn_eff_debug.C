TH3F *h_mass_pt_fine;
TH2F *g_pt_fine;
TH3F *h_mass_pt_fine_wgt;
TH2F *g_pt_fine_wgt;
TH3F *h_mass_pt_fine2_wgt;
TH2F *g_pt_fine2_wgt;
int _cut_mode;
void draw_ld_vn_eff_debug(int cut_mode){

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

//f0 = new TFile(Form("run18_3gev_ld_hist_mu_1x_vn_newreweight_centa0_cut%d.root",_cut_mode),"READ");
//f0 = new TFile(Form("run18_3gev_ld_hist_mu_1x_vn_newreweight_centa1_cut%d.root",_cut_mode),"READ");

//f0 = new TFile(Form("run18_3gev_ld_hist_mu_1x_vn_newreweight_centa0_cut%d_debug.root",_cut_mode),"READ");
f0 = new TFile(Form("run18_3gev_ld_hist_mu_1x_vn_newreweight_centa1_cut%d_debug.root",_cut_mode),"READ");
                     //run18_3gev_ld_hist_mu_1x_vn_newreweight_centa1_cut0_debug.root

h_mass_pt_fine = (TH3F*)f0->Get("h_mass_pt_fine");
g_pt_fine = (TH2F*)f0->Get("g_pt_fine");

h_mass_pt_fine_wgt = (TH3F*)f0->Get("h_mass_pt_fine_wgt");
g_pt_fine_wgt = (TH2F*)f0->Get("g_pt_fine_wgt");

h_mass_pt_fine2_wgt = (TH3F*)f0->Get("h_mass_pt_fine2_wgt");
g_pt_fine2_wgt = (TH2F*)f0->Get("g_pt_fine2_wgt");

TH1F *h_mass_pt_wgt_px;
TH2F *h_pt_wgt_fine;
TH2F *h_pt_wgt_fine2;

TH1F *h_mass_pt_px;
TH2F *h_pt_fine;
TH2F *h_pt_fine2;

 h_pt_wgt_fine = new TH2F("h_pt_wgt_fine",""  ,20,-0.8,1.2,15,0,3);
 h_pt_wgt_fine2 = new TH2F("h_pt_wgt_fine2","",20,-0.8,1.2,15,0,3);

double bincenter;
double bincontent;
double binerror;
double sgm3 = 0.0018;
double lambda_mass = 1.115683;
double counts = 0.;
double countserr = 0.;


TCanvas *yrrr3 = new TCanvas("yrrr3","yrrr3",800,800);
yrrr3->cd();
yrrr3->cd()->SetRightMargin(0.12);
yrrr3->cd()->SetTopMargin(0.02);
yrrr3->cd()->SetLeftMargin(0.12);
yrrr3->cd()->SetBottomMargin(0.14);

double counts = 0.;
double countserr = 0.;

for(int ix=1;ix<=h_pt_wgt_fine->GetNbinsX();ix++ ){//rap
    for(int iy=1;iy<=h_pt_wgt_fine->GetNbinsY();iy++ ){//pt
      h_mass_pt_wgt_px = (TH1F*)h_mass_pt_fine_wgt->ProjectionX("h_mass_pt_wgt_px",ix,ix,iy,iy);
      counts = 0;
      for(int ibin=1;ibin<=h_mass_pt_wgt_px->GetNbinsX();ibin++){
        bincenter = h_mass_pt_wgt_px->GetBinCenter(ibin);
        bincontent = h_mass_pt_wgt_px->GetBinContent(ibin);
        if( bincenter>(lambda_mass-3*sgm3) && bincenter<(lambda_mass+3*sgm3) )  {counts+=bincontent;}
      }
    h_pt_wgt_fine->SetBinContent(ix,iy,counts);
    h_pt_wgt_fine->SetBinError(ix,iy,sqrt(counts));
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
    h_pt_wgt_fine2->SetBinContent(ix,iy,counts);
    h_pt_wgt_fine2->SetBinError(ix,iy,sqrt(counts));
  }
}


yrrr3->cd();
g_pt_fine2_wgt->Draw("colz");
yrrr3->cd();
h_pt_wgt_fine2->Draw("colz");

TH2F *h_pt_eff_wgt_fine;
h_pt_eff_wgt_fine = (TH2F*)h_pt_wgt_fine->Clone("h_pt_eff_wgt_fine");
h_pt_eff_wgt_fine -> Divide(g_pt_fine_wgt);

TH2F *h_pt_eff_wgt_fine2;
h_pt_eff_wgt_fine2 = (TH2F*)h_pt_wgt_fine2->Clone("h_pt_eff_wgt_fine2");
h_pt_eff_wgt_fine2 -> Divide(g_pt_fine2_wgt);
h_pt_eff_wgt_fine2 -> Draw("colz");

yrrr3->Print("h_pt_eff_wgt_fine2.pdf","pdf");
//
//TFile *outhistfile = new TFile (Form("ld_vn_eff_v0_cut%d_centa0_debug.root",_cut_mode), "RECREATE");
TFile *outhistfile = new TFile (Form("ld_vn_eff_v0_cut%d_centa1_debug.root",_cut_mode), "RECREATE");

//TFile *outhistfile = new TFile (Form("ld_vn_eff_v0_cut%d_centa0_rebin.root",_cut_mode), "RECREATE");
//
outhistfile->cd();
h_pt_eff_wgt_fine2 -> Write();
h_pt_eff_wgt_fine -> Write();

h_mass_pt_fine2_wgt-> Write();
g_pt_fine2_wgt -> Write();

h_pt_wgt_fine-> Write();
h_pt_wgt_fine2-> Write();

h_mass_pt_fine_wgt-> Write();
g_pt_fine_wgt -> Write();

outhistfile->Close();
delete outhistfile;


}
