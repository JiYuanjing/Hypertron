void test()
{
  TFile* f  = new TFile("ana_tree.root");

  TCanvas* c = new TCanvas("c","c");
  TTree* t = (TTree*)f->Get("htriton3_tree");
  // t->Draw("bparticlemass","bparticlemass<3.2 && bparticlemass>2.95 && mass_01 <1.15 && mass_01>1.");
  // t->Draw("bparticlemass","bparticlemass<3.05  && ht_chi2topo<10 && ht_chi2ndf<10.0 &&  bparticlemass>2.95 && mass_01 <1.15 && mass_01>1. && ht_chi2>0 && fabs(ht_chi2)<100 && bisMix==1");
  t->Draw("bparticlemass","bparticlemass<3.05 && bparticlemass>2.96&&  mass_01 <1.15 && mass_01>1. && ht_l >5 && ht_ldl>4 && chi2primary_d > 3 && ht_chi2topo<5 && ht_chi2ndf<2.5  && cent9>0&& bisMix==1");
  t->Draw("bparticlemass","bparticlemass<3.05 && bparticlemass>2.96&& mass_01 <1.15 && mass_01>1. &&  ht_l >5 && ht_ldl>4 && chi2primary_d > 3 && ht_chi2topo<5 && ht_chi2ndf<2.5 && bisMix==0  && cent9>0","same");
  // // t->Draw("bparticlemass","bparticlemass<3.05  && ht_chi2topo<10 && ht_chi2ndf<10.0 &&  bparticlemass>2.95 && mass_01 <1.15 && mass_01>1. && ht_chi2>0 && fabs(ht_chi2)<100 && bisMix==0", "same");
  // cout <<  t->GetEntries(" bparticlemass<3.05 && bparticlemass>2.96&&  mass_01 <1.15 && mass_01>1. && ht_l >5 && ht_ldl>4 && chi2primary_d > 3 && ht_chi2topo<5 && ht_chi2ndf<5 && bisMix==1  && cent9>0")<<endl;
  // cout <<  t->GetEntries(" bparticlemass<3.05 && bparticlemass>2.96&&  mass_01 <1.15 && mass_01>1. && ht_l >5 && ht_ldl>4 && chi2primary_d > 3 && ht_chi2topo<5 && ht_chi2ndf<5 && bisMix==0  && cent9>0")<<endl;
  // cout <<  t->GetEntries("bparticlemass<3.05 && ht_l >1 && ht_ldl>1 && ht_chi2topo<10 && ht_chi2ndf<4.&& bparticlemass>2.95 && mass_01 <1.15 && mass_01>1. && ht_chi2>0 && ht_l>1  && ht_ldl>3  &&fabs(ht_chi2)<100 && bisMix==1 && ht_l<200")<<endl;
  // t->Draw("mass_01","bparticlemass<3.2 && bparticlemass>2.95 && mass_01 <1.15 && mass_01>1. && bisMix==1");
  //
  // TCanvas* c1 = new TCanvas("c1","c1");
  // TTree* t1 = (TTree*)f->Get("lambda_tree");
  // t1->Draw("bparticlemass","bparticlemass<1.5 && bparticlemass>1 ");
  TFile* f2  = new TFile("ana_tree_3122_cent.root");

  // TCanvas* c = new TCanvas("c","c");
  TTree* t2 = (TTree*)f2->Get("htriton3_tree");
  // t2->Draw("bparticlemass","bparticlemass<3.2 && bparticlemass>2.95 && mass_01 <1.15 && mass_01>1.");
  // t2->Draw("mass_01","bparticlemass<3.2 && bparticlemass>2.95 && mass_01 <1.15 && mass_01>1.");
}
