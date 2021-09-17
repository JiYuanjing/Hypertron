#include "tree.h"
void readmc(TString mInputlist="data/H3L3b_tree_mc_phase.root", int mode = 0, TString outfile="fH3L_phase_phase_test.root", int const mcH3L=1)
/* void readmc(TString mInputlist="data/Lambda_tree_mc_rotate_June20.root", int const mode = 1,   TString outfile="flambda_phase_test.root", int const mcH3L=1) */
{
  double snn = 3;
  double ycm;
  if(snn==3){
    ycm = -1.045;
  }else{
    ycm = -999.;
  }
 double ht_mass = 2.99131;
 double ht_width = 0.005;
 double hl_mass = 3.9239;
 double hl_width = 0.005;
 double ld_mass = 1.11568;

 TF1 *levyfit4;
 TF1 *levyfit5;
 TF1 *levyfit6;
 TF1 *levyfit7;
 TF1 *levyfit8;
 TF1 *t_quadr;
 TF1 *t_quadr0;
 TF1 *t_quadr1;
 TF1 *bolt0;
 TF1 *bolt1;
 TF1 *bolt2;
 levyfit4 = new TF1("levyfit4","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.4,2.4);
 levyfit5 = new TF1("levyfit5","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.2,2.4);
 levyfit6 = new TF1("levyfit6","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.2,2.4);
 levyfit7 = new TF1("levyfit7","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.2,2.4);
 levyfit8 = new TF1("levyfit8","x*[2]*pow((1+(sqrt(1.11568*1.11568+x*x)-1.11568)/[1]/[0]),(-1)*[1])", 0.2,2.4);
 t_quadr = new TF1("t_quadr","[0]+[1]*x+[2]*x*x",-2,2);
 t_quadr0 = new TF1("t_quadr0","[0]+[1]*x+[2]*x*x",-2,2);
 t_quadr1 = new TF1("t_quadr1","[0]+[1]*x+[2]*x*x",-2,2);
 bolt0 = new TF1("bolt0", "1e9*x*exp(-(sqrt(2.990*2.990+x*x)-2.990)/[0])", 0.,5);
 bolt1 = new TF1("bolt1", "1e9*x*exp(-(sqrt(2.990*2.990+x*x)-2.990)/[0])", 0.,5);
 bolt2 = new TF1("bolt2", "1e9*x*exp(-(sqrt(2.990*2.990+x*x)-2.990)/[0])", 0.,5);


 if (mode==0) {
 levyfit4->SetParameters(1.53536e-01 , -1.21064e+08, 5.25287e+07);
 levyfit5->SetParameters(1.51332e-01 , -1.25927e+08, 4.86193e+07);
 levyfit6->SetParameters(1.45903e-01 , -1.60954e+08, 4.74115e+07);
 levyfit7->SetParameters(1.28382e-01 , -1.60610e+08, 4.92025e+07);
 levyfit8->SetParameters(1.11371e-01 , -2.04689e+08, 4.22255e+07);


 t_quadr->SetParameters(1.25929e+00,0,-1.80963e+00-0.5);
 t_quadr0->SetParameters(1.25929e+00,0,-1.80963e+00);
 t_quadr1->SetParameters(1.25929e+00,0,-1.80963e+00+0.5);
 

 
   bolt0->SetParameter(0,0.34);
   bolt1->SetParameter(0,0.27);
   bolt2->SetParameter(0,0.20);
 }

 if (mode==1)
 {
 levyfit4->SetParameters(0.141849, -1.21226e+08, 1.879e+07);
 levyfit5->SetParameters(0.145031, -1.19088e+08, 1.43524e+07);
 levyfit6->SetParameters(0.141403, -2.20994e+08, 1.2987e+07);
 levyfit7->SetParameters(0.128239, -1.60612e+08, 1.17238e+07);
 levyfit8->SetParameters(0.114106, -1.35632e+08, 8.78861e+06);
 t_quadr->SetParameters(1.15426e+00, 0 ,-8.55891e-01 ); 
 
 }

  TString treename;
  if (mode==1) treename = "lambda_mc_tree";
  if (mode==0) treename = "htriton_mc_tree";
  TChain htriton3_tree(treename.Data()); 

  /////////////////
 TFile*  fgpt_0;
/*  if (mode==0) fgpt_0= new TFile("ht_input_mc_fine_v20.root","READ"); */
/*  if (mode==0) fgpt_0= new TFile("fH3L_phase_wt.root","READ"); */
 if (mode==0) fgpt_0= new TFile("fH3L_phase_phase_test.root","READ");
/*  if (mode==1) fgpt_0 = new TFile("ld_input_mc_fine_v20.root","READ");//more stats */
 if (mode==1) fgpt_0 = new TFile("flambda_phase_wt.root","READ");//more stats
/*  TH2F*  g_pt_fine_in = (TH2F*)fgpt_0->Get("g_pt_fine")->Clone("g_pt_fine_in"); */
 TH2F*  g_pt_fine_in = (TH2F*)fgpt_0->Get("hPhase")->Clone("g_pt_fine_in");



  if (mInputlist.Contains(".root"))
  {
    htriton3_tree.Add(mInputlist.Data());
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
        htriton3_tree.Add(tmp);
        nfile++;
      }
    }
  }

  /* htriton3_tree.SetBranchAddress("bismc", &bismc); */
  /* htriton3_tree.SetBranchAddress("bparticlemass",&bparticlemass); */
  /* htriton3_tree.SetBranchAddress("bpx",&bpx); */
  /* htriton3_tree.SetBranchAddress("bpy",&bpy); */
  /* htriton3_tree.SetBranchAddress("bpz",&bpz); */
  /*  */
  /* htriton3_tree.SetBranchAddress("dca_proton",&dca_proton); */
  /* htriton3_tree.SetBranchAddress("chi2primary_proton", &chi2primary_proton); */
  /* htriton3_tree.SetBranchAddress("chi2primary_pi", &chi2primary_pi); */
  /* htriton3_tree.SetBranchAddress("reweight", &reweight); */
  /* htriton3_tree.SetBranchAddress("cent9", &cent9); */
  /*  */
  /* if (mode==0){ */
  /*   htriton3_tree.SetBranchAddress("chi2primary_d", &chi2primary_d); */
  /*   htriton3_tree.SetBranchAddress("bzdeuteron", &bzdeuteron); */
  /*   htriton3_tree.SetBranchAddress("bpionnsigma", & bpionnsigma); */
  /*   htriton3_tree.SetBranchAddress("bprotonsigma", & bprotonsigma); */
  /*  */
  /*   htriton3_tree.SetBranchAddress("ht_chi2topo", &ht_chi2topo); */
  /*   htriton3_tree.SetBranchAddress("ht_chi2ndf", &ht_chi2ndf); */
  /*   htriton3_tree.SetBranchAddress("ht_ldl", &ht_ldl); */
  /*   htriton3_tree.SetBranchAddress("ht_l", &ht_l); */
  /*   htriton3_tree.SetBranchAddress("ht_dl", &ht_dl); */
  /*   htriton3_tree.SetBranchAddress("dca_proton",&dca_proton); */
  /*   htriton3_tree.SetBranchAddress("dca_pion",&dca_pion); */
  /*   htriton3_tree.SetBranchAddress("dca_deuteron",&dca_deuteron); */
  /*   htriton3_tree.SetBranchAddress("nhits_pion",&nhits_pion); */
  /*   htriton3_tree.SetBranchAddress("nhits_deuteron",&nhits_deuteron); */
  /*   htriton3_tree.SetBranchAddress("nhits_proton",&nhits_proton); */
  /*   //htriton3_tree.SetBranchAddress("ht_bdfvtx",&ht_bdfvtx); */
  /*   //htriton3_tree.SetBranchAddress("ht_bdfvtx2",&ht_bdfvtx2); */
  /*   //htriton3_tree.SetBranchAddress("ht_lifetime",&ht_lifetime); */
  /*   htriton3_tree.SetBranchAddress("countrefmult",&countrefmult); */
  /*   htriton3_tree.SetBranchAddress("reweight", &reweight); */
  /*   htriton3_tree.SetBranchAddress("cent9", &cent9); */
  /*   htriton3_tree.SetBranchAddress("bismc", &bismc); */
  /*  */
  /*   htriton3_tree.SetBranchAddress("bdpx", &bdpx); */
  /*   htriton3_tree.SetBranchAddress("bdpy", &bdpy); */
  /*   htriton3_tree.SetBranchAddress("bdpz", &bdpz); */
  /*   htriton3_tree.SetBranchAddress("bpionpx", &bpionpx); */
  /*   htriton3_tree.SetBranchAddress("bpionpy", &bpionpy); */
  /*   htriton3_tree.SetBranchAddress("bpionpz", &bpionpz); */
  /*   htriton3_tree.SetBranchAddress("bprotonpx", &bprotonpx); */
  /*   htriton3_tree.SetBranchAddress("bprotonpy", &bprotonpy); */
  /*   htriton3_tree.SetBranchAddress("bprotonpz", &bprotonpz); */
  /*  */
  /*   htriton3_tree.SetBranchAddress("v_01_pvdca", &v_01_pvdca); //pair distance from PV  */
  /*   htriton3_tree.SetBranchAddress("v_01_chi2primary", &v_01_chi2primary); // */
  /*   htriton3_tree.SetBranchAddress("v_01_chi2ndf", &v_01_chi2ndf); // */
  /*   htriton3_tree.SetBranchAddress("mass_01", &mass_01); // same as below */
  /*   #<{(| htriton3_tree.SetBranchAddress("v_lambda_mass_0", &v_lambda_mass_0); |)}># */
  /*   htriton3_tree.SetBranchAddress("v_lambda_ldl_0", &v_lambda_ldl_0);  // will add later */
  /*   htriton3_tree.SetBranchAddress("v_lambda_l_0", &v_lambda_l_0);  // will add later */
  /*  */
  /*  */
  /*   #<{(| htriton3_tree.SetBranchAddress("bdpx", &bdpx); |)}># */
  /*   #<{(| htriton3_tree.SetBranchAddress("bdpy", &bdpy); |)}># */
  /*   #<{(| htriton3_tree.SetBranchAddress("bdpz", &bdpz); |)}># */
  /*   htriton3_tree.SetBranchAddress("bpionpx", &bpionpx); */
  /*   htriton3_tree.SetBranchAddress("bpionpy", &bpionpy); */
  /*   htriton3_tree.SetBranchAddress("bpionpz", &bpionpz); */
  /*   htriton3_tree.SetBranchAddress("bprotonpx", &bprotonpx); */
  /*   htriton3_tree.SetBranchAddress("bprotonpy", &bprotonpy); */
  /*   htriton3_tree.SetBranchAddress("bprotonpz", &bprotonpz); */
  /*   htriton3_tree.SetBranchAddress("b0mcpx",&b0mcpx); */
  /*   htriton3_tree.SetBranchAddress("b0mcpy",&b0mcpy); */
  /*   htriton3_tree.SetBranchAddress("b0mcpz",&b0mcpz); */
  /*   htriton3_tree.SetBranchAddress("b1mcpx",&b1mcpx); */
  /*   htriton3_tree.SetBranchAddress("b1mcpy",&b1mcpy); */
  /*   htriton3_tree.SetBranchAddress("b1mcpz",&b1mcpz); */
  /* } */
  htriton3_tree.SetBranchAddress("bmcrawpx",&bmcpx);
  htriton3_tree.SetBranchAddress("bmcrawpy",&bmcpy);                                               
  htriton3_tree.SetBranchAddress("bmcrawpz",&bmcpz);
  /*  */
  /* if (mode ==1) */
  /* { */
  /*   htriton3_tree.SetBranchAddress("ld_chi2ndf",&v_01_chi2ndf); */
  /*   htriton3_tree.SetBranchAddress("ld_chi2primary",&v_01_chi2primary); */
  /*   htriton3_tree.SetBranchAddress("ld_ldl",&v_lambda_ldl_0); */
  /*   htriton3_tree.SetBranchAddress("ld_l",&v_lambda_l_0); */
  /*   htriton3_tree.SetBranchAddress("dca_pi",&dca_pion); */
  /* } */

  /* ////for quick test */
  /* TH2F* hptppimass = new TH2F("hptppimass","hptppimass;p_{T};mass",100,0,10,50,1.06,1.16); */
  /* TH2F* hptppichi2prim= new TH2F("hptppichi2prim","hptppichi2prim;p_{T};(p#pi) #chi^{2}_{prim}",100,0,10,50,0,40); */
  /* TH2F* hptppil= new TH2F("hptppil","hptppil;p_{T};(p#pi) l",100,0,10,50,0,80); */
  /* TH2F* hptppildl= new TH2F("hptppildl","hptppildl;p_{T};(p#pi) ldl",100,0,10,50,0,40); */
  /* TH2F* hptppichi2ndf= new TH2F("hptppichi2ndf","hptppichi2ndf;p_{T};(p#pi) #chi^{2}_{ndf}",100,0,10,80,0,1); */
  /* TH2F* hptpichi2prim= new TH2F("hptpichi2prim","hptpichi2prim;p_{T};#pi #chi^{2}_{prim}",100,0,10,50,0,40); */
  /* TH2F* hptpchi2prim= new TH2F("hptpchi2prim","hptpchi2prim;p_{T};p #chi^{2}_{prim}",100,0,10,50,0,40); */
  /* TH2F* hptpidca = new TH2F("hptpidca","hptpidca;p_{T};DCA",100,0,10,100,0,15); */
  /* TH2F* hptpdca = new TH2F("hptpdca","hptpdca;p_{T};DCA",100,0,10,100,0,10); */
  /* TH2F* hptsumdca = new TH2F("hptsumdca","hptsumdca;p_{T};DCA",100,0,10,100,0,20); */
  TH2F* hPhase = new TH2F("hPhase","hPhase;y;pt",100,-1,1,250,0,5);
  TH2F* hPhase_wt = new TH2F("hPhase_wt","hPhase_wt;y;pt",100,-1,1,250,0,5);

  /* TH2F* hptH3Lmass; */
  /* if (mode==0) { */
  /*   hptH3Lmass  = new TH2F("hptH3Lmass","hptH3Lmass;p_{T};mass",100,0,10,50,2.95,3.05); */
  /* } */

  Long64_t n_lambda_Entries = htriton3_tree.GetEntries();
  for (int i=0;i<n_lambda_Entries;i++)
  {
    htriton3_tree.GetEntry(i); 
    if (i%100000==0) cout <<"read "<<i<<" events!" << endl;
    TLorentzVector mcptc;
    if (mode==0) mcptc.SetXYZT(bmcpx,bmcpy,bmcpz,sqrt(bmcpx*bmcpx+bmcpy*bmcpy+bmcpz*bmcpz+ht_mass*ht_mass));
    else if (mode==1) 
       mcptc.SetXYZT(bmcpx,bmcpy,bmcpz,sqrt(bmcpx*bmcpx+bmcpy*bmcpy+bmcpz*bmcpz+ld_mass*ld_mass));
    double bmcrap = -1*(mcptc.Rapidity() - ycm);
    hPhase->Fill(bmcrap, mcptc.Pt());
    double ptweight = 1; //reserved
    double rapweight = 1;
    if (mode==0) ptweight = bolt1->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy));
    else if (mode==1){
      if (bmcrap<-0.7){ ptweight = levyfit8->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy)); }
      else if (bmcrap<-0.5){ ptweight = levyfit7->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy)); }
      else if (bmcrap<-0.3){ ptweight = levyfit6->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy)); }
      else if (bmcrap<-0.1){ ptweight = levyfit5->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy)); }
      else if (bmcrap<0.1) { ptweight = levyfit4->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy)); }
      else if (bmcrap<0.3) { ptweight = levyfit5->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy)); }
      else if (bmcrap<0.5) { ptweight = levyfit6->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy)); }
      else if (bmcrap<0.7) { ptweight = levyfit7->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy)); }
      else { ptweight = levyfit8->Eval(sqrt(bmcpx*bmcpx+bmcpy*bmcpy)); }
    }
    rapweight = t_quadr->Eval(bmcrap);
    double mcweight = 1./g_pt_fine_in->GetBinContent(g_pt_fine_in->FindBin(bmcrap,sqrt(bmcpx*bmcpx+bmcpy*bmcpy)));;
    reweight=1;

    double weight = ptweight*rapweight*mcweight*reweight;

    hPhase_wt->Fill(bmcrap, mcptc.Pt(),weight);
   }

  TFile* fout = new TFile(outfile.Data(),"recreate");

  /* hptppimass->Write();  */
  /* hptppichi2ndf->Write(); */
  /* hptpichi2prim->Write(); */
  /* hptpchi2prim->Write(); */
  /* hptppichi2prim->Write(); */
  /* hptppil->Write(); */
  /* hptppildl->Write(); */
  /* hptpidca->Write(); */
  /* hptpdca->Write(); */
  /* hptsumdca->Write(); */
  hPhase->Write();
  hPhase_wt->Write();

  /* if (mode==0) */
  /* { */
  /*   hptH3Lmass->Write(); */
  /* } */

  fout->Close();
  // htriton3_tree->Close();
}
