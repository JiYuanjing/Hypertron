/*
   root.exe -q -b -x 'muMc.C(1e6,"../*MuDst.root")'
   */
bool _isPico;
void analysispicodst_h3l3b_rotate_tree(TString picolist="test.list",  TString outFileName="test27gev1.root", bool _isPico=true)
{

  TStopwatch timer;
  timer.Start();

  //bool isPico = true;
  //bool isPico = false;

  bool isPico = _isPico;

#if !defined(__CINT__)
  std::cout << "This code cannot be compiled" << std::endl;
#else
  //  gSystem->SetFPEMask(kInvalid | kDivByZero | kOverflow );
  gROOT->LoadMacro("lMuDst.C");

  TString input;
  TString output;
  int year = 2014;
  //  int year = 2016;
  Int_t N = -999;

  if(isPico)
  {
    //        input = "/star/u/yhleung2/00.test/KFParticle_NewPicoFormat/datafiles/st_physics_19141051_raw_1000007.picoDst.root",      
    input = picolist;

    output = NULL;
    lMuDst(-1,input.Data(),"ry2016,RpicoDst,mysql,kfpAna,quiet,nodefault",output);
  }
  else
  {
    //    input = "/star/u/mzyzak/KFParticle_NewPicoFormat/inputData/*.MuDst.root", 
    input = picolist;    
    //    input = "/star/embed/embedding/AuAu27_production_2011/Lambda_214_2012202/P11id.SL11d_embed/2011/177/12177044/st_physics_adc_12177044_raw_3490002.MuDst.root";
    //output = "mu.root";
    output = NULL;
    lMuDst(-1,input.Data(),"ry2016,picoEvt,RMuDst,mysql,kfpAna,quiet,nodefault",output.Data());
  }

  StKFParticleAnalysisMaker* kfpAnalysis = (StKFParticleAnalysisMaker*) StMaker::GetTopChain()->Maker("KFParticleAnalysis");
  if(!isPico) 
  {
    kfpAnalysis->AnalyseMuDst();
    kfpAnalysis->ProcessSignal();
  }

  if(year == 2016)
  {
    kfpAnalysis->UseTMVA();
    // D0->Kpi
    //    kfpAnalysis->SetTMVABinsD0("0:2:3:4:5:6:7:8:9","-1:1000");
    //    kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_0_1_pt0_80_BDT.weights.xml", 0.075, 0);
  }

  //  kfpAnalysis->CollectPIDHistograms();
  //  kfpAnalysis->CollectTrackHistograms();

  kfpAnalysis->SetDataSet(3);  
  kfpAnalysis->RunCentralityAnalysis();  

  //  kfpAnalysis->StoreTMVANtuples();

  chain->Init();

  //  StKFParticleInterface::instance()->CleanLowPVTrackEvents();
  StKFParticleInterface::instance()->SetSoftKaonPIDMode();
  StKFParticleInterface::instance()->SetSoftTofPidMode();

  //default
  //  StKFParticleInterface::instance()->SetChiPrimaryCut(10);

  StKFParticleInterface::instance()->SetMaxDistanceBetweenParticlesCut(5);

  //loose
  StKFParticleInterface::instance()->SetLCut(1);
  StKFParticleInterface::instance()->SetLdLCut2D(3);

  //very loose
  // StKFParticleInterface::instance()->SetChiPrimaryCut(3);
  StKFParticleInterface::instance()->SetChiPrimaryCut(0);

  // StKFParticleInterface::instance()->SetChiPrimaryCut2D(3);
  StKFParticleInterface::instance()->SetChiPrimaryCut2D(0);


  StKFParticleInterface::instance()->SetLdLCutCharmManybodyDecays(3);
  StKFParticleInterface::instance()->SetChi2TopoCutCharmManybodyDecays(10);
  StKFParticleInterface::instance()->SetChi2CutCharmManybodyDecays(10);

  /*
     StKFParticleInterface::instance()->SetChi2Cut2D(20);
     StKFParticleInterface::instance()->SetChiPrimaryCut(3);
     StKFParticleInterface::instance()->SetChiPrimaryCut2D(3);
     StKFParticleInterface::instance()->SetChiPrimaryCut(3);
     StKFParticleInterface::instance()->SetChiPrimaryCut2D(3);
     StKFParticleInterface::instance()->SetLdLCut2D(3);
     StKFParticleInterface::instance()->SetLdLCutXiOmega(3);
     StKFParticleInterface::instance()->SetChi2TopoCutXiOmega(15);
  */  

  //Add decays to the reconstruction list
  //  StKFParticleInterface::instance()->AddDecayToReconstructionList(  310);
  //
  //  StKFParticleInterface::instance()->AddDecayToReconstructionList( 3122);//lambda0
  //  StKFParticleInterface::instance()->AddDecayToReconstructionList(-3122);//

  //StKFParticleInterface::instance()->AddDecayToReconstructionList( 3312);//cascade
  //  StKFParticleInterface::instance()->AddDecayToReconstructionList(-3312);

  //  StKFParticleInterface::instance()->AddDecayToReconstructionList( 3334);//omega baryon 
  //  StKFParticleInterface::instance()->AddDecayToReconstructionList(-3334);


  //StKFParticleInterface::instance()->AddDecayToReconstructionList(3003);

  // StKFParticleInterface::instance()->AddDecayToReconstructionList( 3122);
  StKFParticleInterface::instance()->AddDecayToReconstructionList( 3012);
  StKFParticleInterface::instance()->AddDecayToReconstructionList( 103004);


  //    StKFParticleInterface::instance()->AddDecayToReconstructionList(-3004);

  //    StKFParticleInterface::instance()->AddDecayToReconstructionList( 3005);
  //    StKFParticleInterface::instance()->AddDecayToReconstructionList(-3005);
  //StKFParticleInterface::instance()->SetMixedEventAnalysis();

  //StKFParticleInterface::instance()->SetRotation();
  //StKFParticleInterface::instance()->SetRotationAngle(TMath::Pi());
  //StKFParticleInterface::instance()->SetRotationPID(-211);//pi minus
  ////  StKFParticleInterface::instance()->SetRotationPID(-321);//K minus

  StKFParticleInterface::instance()->SetRotation();
  StKFParticleInterface::instance()->SetRotationAngle(TMath::Pi());
  StKFParticleInterface::instance()->SetRotationPID(1000010020);//deuteron
  //    StKFParticleInterface::instance()->SetRotationPID(-211);//pi minus  

  Long64_t nevent = N;
  //Long64_t nevent = 1000;
  if(isPico)
  {
    StPicoDstMaker* maker = (StPicoDstMaker *) StMaker::GetTopChain()->Maker("PicoDst");
    if (! maker) return;
    maker->SetStatus("*",1);
    TChain *tree = maker->chain();
    Long64_t nentries = tree->GetEntries();
    if (nentries <= 0) return;
    //    nevent = TMath::Min(nevent,nentries);
    nevent = nentries;    
    //nevent = 100;    
    cout << nentries << " events in chain " << nevent << " will be read." << endl;
  }
  else
  {
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

  // chain->EventLoop(nevent);
  for (int iEvent = 0; iEvent < nevent; ++iEvent)
  {
    chain->Clear();
    if(iEvent && iEvent%10000 == 0) cout<<"... finished processing "<<iEvent<<" events."<<endl;

    int iret = chain->Make();
    if (iret)
    {
      cout << "Bad return code!" << iret << endl;
      break;
    }
  }
  cout<<"Finished processing "<<nevent<<" events."<<endl;


#endif

  cout<<"------"<<endl;
  timer.Stop();
  Double_t rtime = timer.RealTime();
  printf("RealTime=%f minutes", rtime/60);
  cout<<"------"<<endl;
}
