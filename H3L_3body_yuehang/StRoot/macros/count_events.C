Char_t path_ch[10000];
TH1F   *grefmult;
TH1F   *irefmult;
int FXTMult;
int fCentrality;
TH1D *h_cent;

void count_events(){

	h_cent = new TH1D("h_cent","",10,0,10);

	int nChains = -999999;
	ifstream fin;

	fin.open("/star/u/yhleung2/pwg/33.standard_reco/taxi_ht_pico_tree_newdata_chihe/filelist.txt");
        nChains = 1225;

        for(int i=0; i<nChains; i++){
	cout<<"i:"<<i<<endl;
        fin >> path_ch;
//cout<<"ac"<<endl;
        TFile *f0 = new TFile(path_ch,"READ");
//cout<<"ac1"<<endl;
        grefmult = (TH1F*)f0->Get("hrefmult");
//cout<<"ac2"<<endl;
	if(i==0){
		irefmult = (TH1F*)grefmult->Clone("irefmult");
//		irefmult = (TH1F*)f0->Get("irefmult");
	}else{
//cout<<"ac3"<<endl;
//cout<<"int:"<<grefmult->Integral()<<endl;
		irefmult ->Add(grefmult);
//cout<<"ac4"<<endl;
	}
//cout<<"ac5"<<endl;
if(i!=0){
	f0->Close();
        delete f0;
}
	}


	fCentrality = -1;
	for(int i=1;i<=irefmult->GetNbinsX();i++){
//		hrefmult->GetBinContent(i-1);
		FXTMult = i-1;

		if(FXTMult>195) fCentrality= 9;
		else if(FXTMult>141) fCentrality=8;
		else if(FXTMult>118) fCentrality=7;
		else if(FXTMult>85) fCentrality=6;
		else if(FXTMult>59) fCentrality=5;
		else if(FXTMult>40) fCentrality=4;
		else if(FXTMult>25) fCentrality=3;
		else if(FXTMult>15) fCentrality=2;
		else if(FXTMult>8) fCentrality=1;
		else fCentrality=0;

		cout<<"fCentrality:"<<i<<" "<<fCentrality<<" "<<irefmult->GetBinContent(i-1)<<endl;
		h_cent->Fill(fCentrality,irefmult->GetBinContent(i-1));
	}

	TFile *outhistfile = new TFile ("3gev_evt_count.root", "RECREATE");
	outhistfile->cd();
	h_cent-> Write();
	irefmult-> Write();
	outhistfile->Close();
}
