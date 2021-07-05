#include "TLorentzVector.h"
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
#include "TTree.h"

int const N = 10;
int const MixEventBufferSize=5;
int getEPidx(double ep)
{
  //split ep in 5 bins for mix
  double p=TMath::TwoPi();
  if (ep>p) ep=ep-p;
  double edge[N+1];
  edge[N]=p;

  for (int i=0;i<N;i++) {edge[i]=i*p*1.0/N; }
  int idx=-999;
  for (int i=0;i<N;i++)
  {
    if (ep>=edge[i] && ep<edge[i+1]) {idx=i;break;}
  }
  if (idx>=0 && idx<N) return idx;
  else { cout <<"event plane angle out of range! ep=" <<ep << endl; return -999;}
}
struct Events
{
  TLorentzVector lambda;
  TLorentzVector d;
  TLorentzVector d_rt;
  double Ld_wt;
  double d_wt;
};

void read(TString mInputlist="invmass_dLd_3M.root",int mode=1)
{
  const Double_t M_H = 2.99;
  const Double_t M_d = 1.875;
  const Double_t M_Ld = 1.115;

  Float_t         ep;
  Float_t         dpt;
  Float_t         deta;
  Float_t         dphi;
  Float_t         dv1;
  Float_t         dwt;
  Float_t         Lbpt;
  Float_t         Lbeta;
  Float_t         Lbphi;
  Float_t         Lbv1;
  Float_t         Lbptwt;
  Float_t         Lbywt;
  TString treename;
  if (mode==1) treename="dLbTreeFlow";
  if (mode==0) treename="dLbTreeNoFlow";
  TChain fChain(treename.Data()); 
  if (mInputlist.Contains(".root"))
  {
    fChain.Add(mInputlist.Data());
  }
  else
  {
    int nfile = 0;
    char tmp[2000];
    ifstream readlists;
    readlists.open(mInputlist.Data());
    while (readlists.good()){
      readlists.getline(tmp,2000);
      TFile *ftmp = new TFile(tmp);
      if (!ftmp||!(ftmp->IsOpen())||!(ftmp->GetNkeys())) {
        cout<<"Could not open this file: "<< tmp  <<endl;
      }
      else {
        if(nfile%50==0) cout<<"read in "<<nfile<<"th file: "<< tmp <<endl;
        fChain.Add(tmp);
        nfile++;
      }
    }
  }

  fChain.SetBranchAddress("ep", &ep);
  fChain.SetBranchAddress("dpt", &dpt);
  fChain.SetBranchAddress("deta", &deta);
  fChain.SetBranchAddress("dphi", &dphi);
  fChain.SetBranchAddress("dv1", &dv1);
  fChain.SetBranchAddress("dwt", &dwt);
  fChain.SetBranchAddress("Lbpt", &Lbpt);
  fChain.SetBranchAddress("Lbeta", &Lbeta);
  fChain.SetBranchAddress("Lbphi", &Lbphi);
  fChain.SetBranchAddress("Lbv1", &Lbv1);
  fChain.SetBranchAddress("Lbptwt", &Lbptwt);
  fChain.SetBranchAddress("Lbywt", &Lbywt);

  //init hists

  TH2D *hInvMassPt;
  TH2D *hInvMassPtRotate;
  TH2D *hInvMassPtMix;
  TH2D *hInvMassPtMixRotate;
  hInvMassPt = new TH2D(Form("invmass_pt_%d",mode),"",1000, 2.98, 3.28, 100, 0., 5.0);
  hInvMassPt->Sumw2();
  hInvMassPtRotate = new TH2D(Form("invmass_rt_pt_%d",mode),"",1000, 2.98, 3.28, 100, 0., 5.0);
  hInvMassPtRotate->Sumw2();
  hInvMassPtMix = new TH2D(Form("invmass_mx_pt_%d",mode),"",1000, 2.98, 3.28, 100, 0., 5.0);
  hInvMassPtMix->Sumw2();
  hInvMassPtMixRotate = new TH2D(Form("invmass_mxrt_pt_%d",mode),"",1000, 2.98, 3.28, 100, 0., 5.0);
  hInvMassPtMixRotate->Sumw2();

  vector<Events> MixEvent[N]; //split event plane in 10 bins
  Long_t nEntries = fChain.GetEntries();


  for (int i=0;i<nEntries;i++)
  {
    fChain.GetEntry(i); 
    if (i%500000==0) cout << "read "<<i<< endl;  
    TLorentzVector lambda, d, d_rt;
    lambda.SetPtEtaPhiM( Lbpt, Lbeta,Lbphi, M_Ld);
    // Lbywt=1;Lbptwt=1; dwt=1;
    d.SetPtEtaPhiM(dpt,deta,dphi,M_d);
    d_rt.SetPtEtaPhiM(dpt,deta,dphi+TMath::Pi(),M_d);
    TLorentzVector H3L = lambda+d; 
    TLorentzVector H3L_rt = lambda+d_rt;
    hInvMassPt->Fill(H3L.M(), H3L.Pt(), Lbywt*Lbptwt*dwt);
    hInvMassPtRotate->Fill(H3L_rt.M(), H3L_rt.Pt(), Lbywt*Lbptwt*dwt);

    int EPidx = getEPidx(ep);
    if (EPidx<0 || EPidx>N) {cout <<"error event plane idx! idx=" << EPidx << endl;break;}
    //do mixing
    if (MixEvent[EPidx].size()==MixEventBufferSize)
    {
      for (int ie=0;ie<MixEventBufferSize;ie++)
      {  
        //note that in the real analysis more than one lambda and d, will need a loop
        TLorentzVector H3Lmix1 = MixEvent[EPidx][ie].lambda + d;
        hInvMassPtMix->Fill(H3Lmix1.M(), H3Lmix1.Pt(), MixEvent[EPidx][ie].Ld_wt*dwt);
        TLorentzVector H3Lmix2 = MixEvent[EPidx][ie].d + lambda;         
        hInvMassPtMix->Fill(H3Lmix2.M(), H3Lmix2.Pt(), Lbywt*Lbptwt*MixEvent[EPidx][ie].d_wt);
        TLorentzVector H3Lmix1_rt = MixEvent[EPidx][ie].lambda + d_rt;
        hInvMassPtMixRotate->Fill(H3Lmix1_rt.M(), H3Lmix1_rt.Pt(), MixEvent[EPidx][ie].Ld_wt*dwt);
        TLorentzVector H3Lmix2_rt = MixEvent[EPidx][ie].d_rt+ lambda;
        hInvMassPtMixRotate->Fill(H3Lmix2_rt.M(), H3Lmix2_rt.Pt(), Lbywt*Lbptwt*MixEvent[EPidx][ie].d_wt);
      }
      //delete the last event 
      MixEvent[EPidx].erase(MixEvent[EPidx].begin());

    } //mix event
    //push current event into vector
    Events newEvent;
    newEvent.lambda = lambda;
    newEvent.d = d;
    newEvent.d_rt = d_rt;
    newEvent.d_wt = dwt;
    newEvent.Ld_wt=Lbywt*Lbptwt;

    MixEvent[EPidx].push_back(newEvent);
  }

  //write the files 
  TFile* fout = new TFile("dLdInvMass_ME.root","recreate");
  hInvMassPtMix->Write();
  hInvMassPtMixRotate->Write();
  hInvMassPtRotate->Write();
  hInvMassPt->Write();
  TH1D* hMass;
  TH1D* hMassRotate;
  TH1D* hMassME;
  TH1D* hMassMERotate;
  hMass = (TH1D *)hInvMassPt->ProjectionX(Form("mass_%d",mode));
  hMassRotate = (TH1D *)hInvMassPtRotate->ProjectionX(Form("mass_rotate_%d",mode));
  hMassMERotate = (TH1D *)hInvMassPtMixRotate->ProjectionX(Form("mass_ME_rotate_%d",mode));
  hMassME= (TH1D *)hInvMassPtMix->ProjectionX(Form("mass_ME_%d",mode));
  hMass->Write();
  hMassME->Write();
  hMassMERotate->Write();
  hMassRotate->Write();

  fout->Close();

} 
