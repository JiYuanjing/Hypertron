void draw()
{
  TF1* levyfit4;
  TF1* levyfit5;
  TF1* levyfit6;
  TF1* levyfit7;
  TF1* levyfit8;
  TF1* t_quadr;
  TF1* t_quadr0;
  TF1* t_quadr1;
  TF1* bolt0;
  TF1* bolt1;
  TF1* bolt2;
  //weighting function for MC
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
    // ptmptexp[] = new TF1(Form("ptmptexp[%d]",irap), "x*[1]*exp(-(sqrt([2]*[2]+x*x)-[2])/[0])", 0,4);

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

    // t_quadr0->Draw();
    // bolt1->Draw();

    // levyfit4->SetParameters(0.141849, -1.21226e+08, 1.879e+07);
    // levyfit5->SetParameters(0.145031, -1.19088e+08, 1.43524e+07);
    // levyfit6->SetParameters(0.141403, -2.20994e+08, 1.2987e+07);
    // levyfit7->SetParameters(0.128239, -1.60612e+08, 1.17238e+07);
    // levyfit8->SetParameters(0.114106, -1.35632e+08, 8.78861e+06);
    // t_quadr->SetParameters(1.15426e+00, 0 ,-8.55891e-01 ); 



    TLegend* leg = new TLegend(0.2,0.2, 0.5,0.5);
    TF1* ptmptexp[3];
    double par[3][3]={{0.267422, 43782.6, 2.99339}, {0.200836, 92736.2,2.99339}, {0.136789, 136882, 2.99339} };
    for (int irap=0;irap<3;irap++) {
        ptmptexp[irap] = new TF1(Form("ptmptexp[%d]",irap), "x*[1]*exp(-(sqrt([2]*[2]+x*x)-[2])/[0])", 0,4);
        ptmptexp[irap]->SetParameters(par[irap]);
        // cout << ptmptexp[irap]->GetParameter(0)<<endl;
        // cout << ptmptexp[irap]->GetParameter(1)<<endl;
        // cout << ptmptexp[irap]->GetParameter(2)<<endl;
        if (irap>0) ptmptexp[irap]->Draw("same");
        else ptmptexp[irap]->Draw();
        ptmptexp[irap]->GetYaxis()->SetRangeUser(1, 0.6e5);
        ptmptexp[irap]->GetYaxis()->SetTitle("dN/dpTdy");
        ptmptexp[irap]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        gPad->SetLogy();

        if (irap==0) ptmptexp[0]->SetLineColor(kRed);
        if (irap==1) ptmptexp[1]->SetLineColor(kGreen);
        if (irap==2) ptmptexp[2]->SetLineColor(kBlue);

    }
    
    leg->AddEntry(ptmptexp[0], "-0.25<y<0", "pl" );
    leg->AddEntry(ptmptexp[1], "-0.5<y<-0.25", "pl" );
    leg->AddEntry(ptmptexp[2], "-0.75<y<-0.5", "pl" );
    leg->Draw();
}
