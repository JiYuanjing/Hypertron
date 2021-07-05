bool _isPico;

void analysispicodst_ldxi_pol_rotate(TString picolist="test.list",  TString outFileName="test27gev1.root", bool _isPico=true)
{
    
    TStopwatch timer;
    timer.Start();
    
    bool isPico = _isPico;
    
#if !defined(__CINT__)
    std::cout << "This code cannot be compiled" << std::endl;
#else
    
    gSystem->Load("StRefMultCorr");
    gROOT->LoadMacro("lMuDst.C");
    gSystem->Load("EventClass");
    
    TString input;
    TString output;
    
    Int_t N = -999;
    //Int_t N = 1000;
    if(isPico)
    {
        input = picolist;
        output = NULL;
        lMuDst(-1,input.Data(),"ry2016,RpicoDst,mysql,kfpAna,quiet,nodefault",output);
    }
    else
    {
        input = picolist;
        output = NULL;
        lMuDst(-1,input.Data(),"ry2016,picoEvt,RMuDst,mysql,kfpAna,quiet,nodefault",output.Data());
    }
    
    StKFParticleAnalysisMaker* kfpAnalysis = (StKFParticleAnalysisMaker*) StMaker::GetTopChain()->Maker("KFParticleAnalysis");
    if(!isPico)
    {
        kfpAnalysis->AnalyseMuDst();
        kfpAnalysis->ProcessSignal();
    }
    
    //kfpAnalysis->RunCentralityAnalysis();
    
    kfpAnalysis->SetDataSet(3);
    kfpAnalysis->RunFlowAnalysis();
    kfpAnalysis->AddFlowFile("test_v2.root");
    
    chain->Init();
    
    //StKFParticleInterface::instance()->CleanLowPVTrackEvents();
    StKFParticleInterface::instance()->SetSoftKaonPIDMode();
    StKFParticleInterface::instance()->SetSoftTofPidMode();
    
    //default
    StKFParticleInterface::instance()->SetChiPrimaryCut(3);
    StKFParticleInterface::instance()->SetMaxDistanceBetweenParticlesCut(1.5); //maximum distance between secondary tracks or particles
    StKFParticleInterface::instance()->SetLCut(1.0);  //minimum distance to PV, used for K0s, Ld, Xi, Omega, hypernuclei, dibaryons.

    StKFParticleInterface::instance()->SetChiPrimaryCut2D(3); //default: > 3.0, set a cut on chi2 of the track to the primary vertex in case of 2 daughter decays
    StKFParticleInterface::instance()->SetChi2Cut2D(10); //default: < 10.0, set a cut on chi2 of the particle fit in case of 2 daughter decays
    StKFParticleInterface::instance()->SetLdLCut2D(3); //default: > 5.0, set a cut on a distance to PV normalised on the error for non-resonant 2 daughter decays
    
    StKFParticleInterface::instance()->SetLdLCutXiOmega(3);//default: < 10, set a cut on a distance to PV normalised on the error for Xi and Omega
    StKFParticleInterface::instance()->SetChi2TopoCutXiOmega(10);//default < 5, set a cut on chi2 after fitting PV to the particle for Xi and Omega
    StKFParticleInterface::instance()->SetChi2CutXiOmega(10); //default > 6, set a cut on χ2 of the particle fit for Xi and Omega


    //Add decays to the reconstruction list
    StKFParticleInterface::instance()->AddDecayToReconstructionList( 3122);//lambda
    StKFParticleInterface::instance()->AddDecayToReconstructionList(-3122);
    
    StKFParticleInterface::instance()->AddDecayToReconstructionList( 3312);//cascade
    StKFParticleInterface::instance()->AddDecayToReconstructionList(-3312);
    
    StKFParticleInterface::instance()->AddDecayToReconstructionList( 3334);//omega
    StKFParticleInterface::instance()->AddDecayToReconstructionList(-3334);
    
    StKFParticleInterface::instance()->SetRotation();
    StKFParticleInterface::instance()->SetRotationAngle(TMath::Pi());
    StKFParticleInterface::instance()->SetRotationPID(-211);//pi minus

    Long64_t nevent = N;
    if(isPico)
    {
        StPicoDstMaker* maker = (StPicoDstMaker *) StMaker::GetTopChain()->Maker("PicoDst");
        if (! maker) return;
        maker->SetStatus("*",1);
        TChain *tree = maker->chain();
        Long64_t nentries = tree->GetEntries();
        if (nentries <= 0) return;
        //nevent = TMath::Min(nevent,nentries);
        nevent = nentries;
        cout << nentries << " events in chain " << nevent << " will be read." << endl;
    }else{
        StMuDstMaker* maker = (StMuDstMaker *) StMaker::GetTopChain()->Maker("MuDst");
        if (! maker) return;
        maker->SetStatus("*",1);
        TChain *tree = maker->chain();
        Long64_t nentries = tree->GetEntries();
        if (nentries <= 0) return;
        nevent = TMath::Min(nevent,nentries);
        nevent = nentries;
        cout << nentries << " events in chain " << nevent << " will be read." << endl;
        
    }
    
    chain->EventLoop(nevent);
#endif
    
    cout<<"------"<<endl;
    timer.Stop();
    Double_t rtime = timer.RealTime();
    printf("RealTime=%f minutes", rtime/60);
    cout<<"------"<<endl;
}
