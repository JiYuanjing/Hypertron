#include "tree.h"
#include "Hists.h"
// void readtree(TString mInputlist="H3L3b_tree_mc.root", int mode = 0, TString outfile="fout_H3L.root")
void readtree(TString mInputlist="Lambda_tree_mc.root", int const mode = 1,   TString outfile="fout_Lambda.root", int const mcState=1, int const isMix=0)
{
  double snn = 3;
  double ycm;
  if(snn==3) //target
  {
    ycm = -1.045;
  }else{
    ycm = -999.;
  }
  double const ht_mass = 2.99131;
  double const ht_width = 0.005;
  double const hl_mass = 3.9239;
  double const hl_width = 0.005;
  double const ld_mass = 1.11568;
  double const mass_p = 0.93827;
  double const mass_pi = 0.13957;

  //weighting function for MC
  if (mcState!=0){
    cout <<"starting book mc functions!" << endl;
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


    if (mode==0 && mcState==1) {
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

    if (mode==1 || (mode ==0 && mcState==-20)) //MC Lambda or (MC Lambda)+d
    {
      levyfit4->SetParameters(0.141849, -1.21226e+08, 1.879e+07);
      levyfit5->SetParameters(0.145031, -1.19088e+08, 1.43524e+07);
      levyfit6->SetParameters(0.141403, -2.20994e+08, 1.2987e+07);
      levyfit7->SetParameters(0.128239, -1.60612e+08, 1.17238e+07);
      levyfit8->SetParameters(0.114106, -1.35632e+08, 8.78861e+06);
      t_quadr->SetParameters(1.15426e+00, 0 ,-8.55891e-01 ); 
    }

    ////////uniform  mc pt y distribution/////////
    /*  if (mode==0) fgpt_0= new TFile("ht_input_mc_fine_v20.root","READ"); */
    if (mode==0 && mcState==1) fgpt_0= new TFile("fH3L_phase_quasi_wt.root","READ");
    if (mode==0 && mcState==1 && mInputlist.Contains("phase")) fgpt_0= new TFile("fH3L_phase_phase_wt.root","READ");  //phase space decay
    /*  if (mode==1) fgpt_0 = new TFile("ld_input_mc_fine_v20.root","READ");//more stats */
    if (mode==1 || (mcState==-20 && mode==0)) fgpt_0 = new TFile("flambda_phase_wt.root","READ");//more stats
    /*  TH2F*  g_pt_fine_in = (TH2F*)fgpt_0->Get("g_pt_fine")->Clone("g_pt_fine_in"); */
    g_pt_fine_in = (TH2F*)fgpt_0->Get("hPhase")->Clone("g_pt_fine_in");
    cout <<"finish book mc weighting functions." << endl;
  }
  ///////////////////////end of MC weighting////////////////

  cout <<"start read tree!" << endl;
  TString treename;
  if (mode==1) treename = "lambda_tree";
  if (mode==0) treename = "htriton3_tree";
  TChain htriton3_tree(treename.Data()); 

  TH1F* hrefmult  = new TH1F("hrefmult", "refmult; hrefmult; N_{evt}", 600,0,600);
  hrefmult->SetDirectory(0);

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
        TH1F* htmp = (TH1F*)ftmp->Get("hrefmult");
        hrefmult->Add(htmp);
        nfile++;
      }
    }
  }

  htriton3_tree.SetBranchAddress("bismc", &bismc);
  htriton3_tree.SetBranchAddress("bparticlemass",&bparticlemass);
  htriton3_tree.SetBranchAddress("bparticleid",&bparticleid);
  htriton3_tree.SetBranchAddress("bpx",&bpx);
  htriton3_tree.SetBranchAddress("bpy",&bpy);
  htriton3_tree.SetBranchAddress("bpz",&bpz);

  htriton3_tree.SetBranchAddress("dca_proton",&dca_proton);
  htriton3_tree.SetBranchAddress("chi2primary_proton", &chi2primary_proton);
  htriton3_tree.SetBranchAddress("chi2primary_pi", &chi2primary_pi);
  htriton3_tree.SetBranchAddress("reweight", &reweight);
  htriton3_tree.SetBranchAddress("cent9", &cent9);

  if (mode==0){
    htriton3_tree.SetBranchAddress("chi2primary_d", &chi2primary_d);
    htriton3_tree.SetBranchAddress("bzdeuteron", &bzdeuteron);
    htriton3_tree.SetBranchAddress("bpionnsigma", & bpionnsigma);
    htriton3_tree.SetBranchAddress("bprotonsigma", & bprotonsigma);

    htriton3_tree.SetBranchAddress("ht_chi2topo", &ht_chi2topo);
    htriton3_tree.SetBranchAddress("ht_chi2ndf", &ht_chi2ndf);
    htriton3_tree.SetBranchAddress("ht_NDF", &ht_NDF);
    htriton3_tree.SetBranchAddress("ht_ldl", &ht_ldl);
    htriton3_tree.SetBranchAddress("ht_l", &ht_l);
    htriton3_tree.SetBranchAddress("ht_dl", &ht_dl);
    htriton3_tree.SetBranchAddress("dca_proton",&dca_proton);
    htriton3_tree.SetBranchAddress("dca_pion",&dca_pion);
    htriton3_tree.SetBranchAddress("dca_deuteron",&dca_deuteron);
    htriton3_tree.SetBranchAddress("nhits_pion",&nhits_pion);
    htriton3_tree.SetBranchAddress("nhits_deuteron",&nhits_deuteron);
    htriton3_tree.SetBranchAddress("nhits_proton",&nhits_proton);
    //htriton3_tree.SetBranchAddress("ht_bdfvtx",&ht_bdfvtx);
    //htriton3_tree.SetBranchAddress("ht_bdfvtx2",&ht_bdfvtx2);
    //htriton3_tree.SetBranchAddress("ht_lifetime",&ht_lifetime);
    htriton3_tree.SetBranchAddress("countrefmult",&countrefmult);
    htriton3_tree.SetBranchAddress("reweight", &reweight);
    htriton3_tree.SetBranchAddress("cent9", &cent9);
    htriton3_tree.SetBranchAddress("bismc", &bismc);

    htriton3_tree.SetBranchAddress("bdpx", &bdpx);
    htriton3_tree.SetBranchAddress("bdpy", &bdpy);
    htriton3_tree.SetBranchAddress("bdpz", &bdpz);
    htriton3_tree.SetBranchAddress("bpionpx", &bpionpx);
    htriton3_tree.SetBranchAddress("bpionpy", &bpionpy);
    htriton3_tree.SetBranchAddress("bpionpz", &bpionpz);
    htriton3_tree.SetBranchAddress("bprotonpx", &bprotonpx);
    htriton3_tree.SetBranchAddress("bprotonpy", &bprotonpy);
    htriton3_tree.SetBranchAddress("bprotonpz", &bprotonpz);

    htriton3_tree.SetBranchAddress("v_01_pvdca", &v_01_pvdca); //pair distance from PV 
    htriton3_tree.SetBranchAddress("v_12_dca", &v_12_dca); //pair 
    htriton3_tree.SetBranchAddress("v_01_chi2primary", &v_01_chi2primary); //
    htriton3_tree.SetBranchAddress("v_01_chi2ndf", &v_01_chi2ndf); //
    htriton3_tree.SetBranchAddress("mass_01", &mass_01); // same as below
    /* htriton3_tree.SetBranchAddress("v_lambda_mass_0", &v_lambda_mass_0); */
    htriton3_tree.SetBranchAddress("v_lambda_ldl_0", &v_lambda_ldl_0);  // will add later
    htriton3_tree.SetBranchAddress("v_lambda_l_0", &v_lambda_l_0);  // will add later

    /* htriton3_tree.SetBranchAddress("bdpx", &bdpx); */
    /* htriton3_tree.SetBranchAddress("bdpy", &bdpy); */
    /* htriton3_tree.SetBranchAddress("bdpz", &bdpz); */
    htriton3_tree.SetBranchAddress("bpionpx", &bpionpx);
    htriton3_tree.SetBranchAddress("bpionpy", &bpionpy);
    htriton3_tree.SetBranchAddress("bpionpz", &bpionpz);
    htriton3_tree.SetBranchAddress("bprotonpx", &bprotonpx);
    htriton3_tree.SetBranchAddress("bprotonpy", &bprotonpy);
    htriton3_tree.SetBranchAddress("bprotonpz", &bprotonpz);
    htriton3_tree.SetBranchAddress("b0mcpx",&b0mcpx);
    htriton3_tree.SetBranchAddress("b0mcpy",&b0mcpy);
    htriton3_tree.SetBranchAddress("b0mcpz",&b0mcpz);
    htriton3_tree.SetBranchAddress("b1mcpx",&b1mcpx);
    htriton3_tree.SetBranchAddress("b1mcpy",&b1mcpy);
    htriton3_tree.SetBranchAddress("b1mcpz",&b1mcpz);
    htriton3_tree.SetBranchAddress("bisMix",&bisMix);
  }
  htriton3_tree.SetBranchAddress("bmcpx",&bmcpx);
  htriton3_tree.SetBranchAddress("bmcpy",&bmcpy);                                               
  htriton3_tree.SetBranchAddress("bmcpz",&bmcpz);
  htriton3_tree.SetBranchAddress("bmcl",&bmcl);
  htriton3_tree.SetBranchAddress("bmcpl",&bmcpl);

  if (mode ==1)
  {
    htriton3_tree.SetBranchAddress("ld_chi2ndf",&v_01_chi2ndf);
    htriton3_tree.SetBranchAddress("ld_chi2primary",&v_01_chi2primary);
    htriton3_tree.SetBranchAddress("ld_ldl",&v_lambda_ldl_0);
    htriton3_tree.SetBranchAddress("ld_l",&v_lambda_l_0);
    htriton3_tree.SetBranchAddress("dca_pi",&dca_pion);
  }

  ////for quick test
  TH2F* hptppimass = new TH2F("hptppimass","hptppimass;p_{T};mass",100,0,10,400,1.,1.2);
  TH2F* hptppichi2prim= new TH2F("hptppichi2prim","hptppichi2prim;p_{T};(p#pi) #chi^{2}_{prim}",100,0,10,50,0,40);
  TH2F* hptppil= new TH2F("hptppil","hptppil;p_{T};(p#pi) l",100,0,10,50,0,80);
  TH2F* hptppildl= new TH2F("hptppildl","hptppildl;p_{T};(p#pi) ldl",100,0,10,50,0,40);
  TH2F* hptppichi2ndf= new TH2F("hptppichi2ndf","hptppichi2ndf;p_{T};(p#pi) #chi^{2}_{ndf}",100,0,10,80,0,1);
  TH2F* hptpichi2prim= new TH2F("hptpichi2prim","hptpichi2prim;p_{T};#pi #chi^{2}_{prim}",100,0,10,50,0,40);
  TH2F* hptpchi2prim= new TH2F("hptpchi2prim","hptpchi2prim;p_{T};p #chi^{2}_{prim}",100,0,10,50,0,40);
  TH2F* hptpidca = new TH2F("hptpidca","hptpidca;p_{T};DCA",100,0,10,100,0,15);
  TH2F* hptpdca = new TH2F("hptpdca","hptpdca;p_{T};DCA",100,0,10,100,0,10);
  TH2F* hptsumdca = new TH2F("hptsumdca","hptsumdca;p_{T};DCA",100,0,10,100,0,20);
  /* TH2F* hPhase = new TH2F("hPhase","hPhase;y;pt",100,-1,1,250,0,5); */

  //topological variable for H3L
  TH2F* hptH3Lmass  = new TH2F("hptH3Lmass","hptH3Lmass;p_{T};mass",100,0,10,100,2.95,3.05);
  hptH3Lmass->Sumw2();
  TH2F* hptH3L_l= new TH2F("hptH3L_l","hptH3L_l;p_{T};l",100,0,10,100,0,20);
  hptH3L_l->Sumw2();
  TH2F* hptH3L_ldl= new TH2F("hptH3L_ldl","hptH3L_ldl;p_{T};ldl",100,0,10,100,0,20);
  hptH3L_ldl->Sumw2();
  TH2F* hptH3L_dchi2prim= new TH2F("hptH3L_dchi2prim","hptH3L_dchi2pri;p_{T};d chi2prim",100,0,5,100,0,20);
  hptH3L_dchi2prim->Sumw2();
  TH2F* hptH3L_pichi2prim= new TH2F("hptH3L_pichi2prim","hptH3L_pichi2pri;p_{T};pi chi2prim",100,0,5,100,0,20);
  hptH3L_pichi2prim->Sumw2();
  TH2F* hptH3L_pchi2prim= new TH2F("hptH3L_pchi2prim","hptH3L_pchi2pri;p_{T};p chi2prim",100,0,5,100,0,20);
  hptH3L_pchi2prim->Sumw2();
  TH2F* hptH3L_chi2topo= new TH2F("hptH3L_chi2topo","hptH3L_chi2topo;p_{T};chi2ptopo",100,0,5,100,0,20);
  hptH3L_chi2topo->Sumw2();
  TH2F* hptH3L_chi2ndf= new TH2F("hptH3L_chi2ndf","hptH3L_chi2topo;p_{T};chi2ndf",100,0,5,100,0,20);
  hptH3L_chi2ndf->Sumw2();
  TH2F* hptH3L_dDca= new TH2F("hptH3L_dDca","hptH3L_dDca;p_{T};dDca",100,0,5,100,0,10);
  hptH3L_dDca->Sumw2();
  TH2F* hptH3L_piDca= new TH2F("hptH3L_piDca","hptH3L_piDca;p_{T};piDca",100,0,5,100,0,20);
  hptH3L_piDca->Sumw2();
  TH2F* hptH3L_pDca= new TH2F("hptH3L_pDca","hptH3L_pDca;p_{T};pDca",100,0,5,100,0,20);
  hptH3L_pDca->Sumw2();
  TH2F* hptH3L_dpDca= new TH2F("hptH3L_dpDca","hptH3L_dpDca;p_{T};dpDca",100,0,5,100,0,20);
  hptH3L_dpDca->Sumw2();

  Long64_t n_lambda_Entries = htriton3_tree.GetEntries();
  for (int i=0;i<n_lambda_Entries;i++)
  {
    htriton3_tree.GetEntry(i); 
    if (i%100000==0) cout <<"read "<<i<<" events!" << endl;
    if ( !(bismc == mcState)  ) continue;
    if (bisMix!=isMix) continue;

    double ptweight = 1; //reserved
    double rapweight = 1;
    /* mcMotherPt =sqrt(bmcpx*bmcpx+bmcpy*bmcpy) ; */

    /* cout << bmcrap<< endl; */
    /* double mcweight = 1./g_pt_fine_in->GetBinContent(g_pt_fine_in->FindBin(bmcrap,sqrt(bmcpx*bmcpx+bmcpy*bmcpy)));; */
    double mcweight = 1;
    if (mcState!=0) {
      TLorentzVector mcptc;
      if (mode==0 && mcState==1 ) 
        mcptc.SetXYZM(bmcpx,bmcpy,bmcpz,ht_mass);
      else if (mode==0 && mcState==-20 ) {
        TLorentzVector pion;
        pion.SetXYZM(b0mcpx, b0mcpy, b0mcpz, mass_pi );
        TLorentzVector proton;
        proton.SetXYZM( b1mcpx, b1mcpy, b1mcpz, mass_p );
        mcptc = pion+proton; //using lambda pt
        /* cout <<mcptc.Pt() << endl; */
      }
      else if (mode==1 && mcState==1) 
        mcptc.SetXYZM(bmcpx,bmcpy,bmcpz,ld_mass);
      double bmcrap = mcptc.Rapidity() - ycm;
      /* hPhase->Fill(bmcrap, mcptc.Pt()); */
      double mcMotherPt = mcptc.Pt();
      rapweight = t_quadr->Eval(bmcrap);
      mcweight = 1./g_pt_fine_in->GetBinContent(g_pt_fine_in->FindBin(bmcrap,mcMotherPt));
      if (mode==0 && mcState ==1) ptweight = bolt1->Eval(mcMotherPt);
      else if (mode==1 || mcState==-20){
        if (bmcrap<-0.7){ ptweight = levyfit8->Eval(mcMotherPt); }
        else if (bmcrap<-0.5){ ptweight = levyfit7->Eval(mcMotherPt); }
        else if (bmcrap<-0.3){ ptweight = levyfit6->Eval(mcMotherPt); }
        else if (bmcrap<-0.1){ ptweight = levyfit5->Eval(mcMotherPt); }
        else if (bmcrap<0.1) { ptweight = levyfit4->Eval(mcMotherPt); }
        else if (bmcrap<0.3) { ptweight = levyfit5->Eval(mcMotherPt); }
        else if (bmcrap<0.5) { ptweight = levyfit6->Eval(mcMotherPt); }
        else if (bmcrap<0.7) { ptweight = levyfit7->Eval(mcMotherPt); }
        else { ptweight = levyfit8->Eval(mcMotherPt); }
      }
    }
    /* cout << mcMotherPt <<" "<<mcweight<<" "<<ptweight<<" "<<rapweight<<endl; */
    reweight=1; //centrality weight

    double weight = ptweight*rapweight*mcweight*reweight;
    if (mcState==0) weight = reweight; 
    double ppi_pt = sqrt(bpx*bpx+bpy*bpy); // lambda case 
    if (mode==0) ppi_pt= sqrt(bpionpx*bpionpx+bpionpy*bpionpy+bprotonpy*bprotonpy+bprotonpx*bprotonpx);
    if (mode==0) {
      if (bparticlemass<2.95 || bparticlemass > 3.05) continue;
      if (mass_01>1.15 || mass_01< 1.05) continue;
    }
    if (bparticleid<0) continue; //currently only look at particle

    //compare H3L->ppi and Lambda->ppi
    if ((mode==0 && bismc==mcState) || (mode==1 && bismc==mcState) ){ 
      // cout << reweight << " "<<ppi_pt<< " "<<dca_pion<<" "<<chi2primary_proton<<endl;
      hptpdca->Fill( ppi_pt,dca_proton, weight);
      hptpidca->Fill( ppi_pt,dca_pion, weight);
      hptsumdca->Fill( ppi_pt,dca_pion+dca_proton, weight);
      hptpchi2prim->Fill( ppi_pt,chi2primary_proton, weight);
      hptpichi2prim->Fill( ppi_pt,chi2primary_pi, weight);
      if (mode==0) 
      {
        if (bparticlemass<3.01  && bparticlemass >2.98) hptppimass->Fill( ppi_pt, mass_01, weight);
      }
      else if (mode==1) {
        hptppimass->Fill( ppi_pt, bparticlemass, weight);
      }
      hptppichi2ndf->Fill( ppi_pt, v_01_chi2ndf, weight);
      hptppichi2prim->Fill( ppi_pt, v_01_chi2primary, weight);
      hptppil->Fill( ppi_pt, v_lambda_l_0, weight);
      hptppildl->Fill( ppi_pt, v_lambda_ldl_0, weight);
    }

    //compare H3L and Lambda+d
    if (mode==0  ) {
      double H3LpT = sqrt(bpx*bpx+bpy*bpy);

      hptH3L_chi2topo->Fill( H3LpT, ht_chi2topo, weight);
      hptH3L_chi2ndf->Fill( H3LpT, ht_chi2ndf, weight);
      hptH3L_l->Fill( H3LpT, ht_l, weight);
      hptH3L_ldl->Fill( H3LpT, ht_ldl, weight);

      hptH3L_pchi2prim->Fill( H3LpT, chi2primary_proton, weight);
      hptH3L_pichi2prim->Fill( H3LpT, chi2primary_pi, weight);
      hptH3L_dchi2prim->Fill( H3LpT, chi2primary_d, weight);
      hptH3L_dDca->Fill( H3LpT, dca_deuteron, weight);
      hptH3L_pDca->Fill( H3LpT, dca_proton, weight);
      hptH3L_piDca->Fill( H3LpT, dca_pion, weight);
      hptH3L_dpDca->Fill( H3LpT, v_12_dca, weight);
      /* weight=1; */
      /* bool passTopoCuts=1; */
      // bool passTopoCuts =  ht_l >8 && ht_ldl>10 && dca_pion>1.5 && dca_proton>0.5 && dca_proton<5 && dca_deuteron<2 && ht_chi2topo<10 && 
                           // ht_chi2ndf<4.0 && mass_01 <1.15;
      bool passTopoCuts =  ht_l >10 && ht_ldl>8  && ht_chi2topo<10 && ht_chi2ndf<10;
      if ( passTopoCuts) hptH3Lmass->Fill(H3LpT, bparticlemass, weight); 
    }
  }

  TFile* fout = new TFile(outfile.Data(),"recreate");

  hptppimass->Write(); 
  hptppichi2ndf->Write();
  hptpichi2prim->Write();
  hptpchi2prim->Write();
  hptppichi2prim->Write();
  hptppil->Write();
  hptppildl->Write();
  hptpidca->Write();
  hptpdca->Write();
  hptsumdca->Write();
  /* hPhase->Write(); */

  if (mode==0)
  {
    hptH3Lmass->Write();
    hptH3L_dDca->Write();
    hptH3L_piDca->Write();
    hptH3L_pDca->Write();
    hptH3L_dpDca->Write();
    hptH3L_chi2ndf->Write();
    hptH3L_chi2topo->Write();
    hptH3L_l->Write();
    hptH3L_ldl->Write();
    hptH3L_pichi2prim->Write();
    hptH3L_pchi2prim->Write();
    hptH3L_dchi2prim->Write();
  }

  hrefmult->Write();

  fout->Close();
  // htriton3_tree->Close();
}
