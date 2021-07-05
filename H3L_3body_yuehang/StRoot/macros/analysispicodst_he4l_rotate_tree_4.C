/*
  root.exe -q -b -x 'muMc.C(1e6,"../*MuDst.root")'
*/
bool _isPico;
void analysispicodst_he4l_rotate_tree_4(TString picolist="test.list",  TString outFileName="test27gev1.root", bool _isPico=true)
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

  gSystem->Load("StRefMultCorr");
  gROOT->LoadMacro("lMuDst.C");
  
  TString input;
  TString output;
  int year = 2014;
//  int year = 2016;
  Int_t N = -999;
              
  if(isPico)
  {
//    input = "/star/u/mzyzak/KFParticle_NewPicoFormat/inputData/*.picoDst.root", 
//    input = "/star/u/yhleung2/00.test/KFParticle_NewPicoFormat/datafiles/*.picoDst.root",    
//      input = "root://xrdstar.rcf.bnl.gov:1095///home/starreco/reco/27GeV_production_2018/ReversedFullField/P18ih/2018/141/19141028/st_physics_19141028_raw_1500006.picoDst.root",// direct from hpss


//      input = "/star/u/yhleung2/00.test/KFParticle_NewPicoFormat/datafiles/*.root",      
//        input = "/star/u/yhleung2/00.test/KFParticle_NewPicoFormat/datafiles/st_physics_19141051_raw_1000007.picoDst.root",      
        input = picolist;
    

//    output = "pico.root";
//    output = "testpico6.root";    
//    output = "testpico8.root";    
    //output = outFileName;
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
//    kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_2_2_pt0_80_BDT.weights.xml", 0.05,  1);
//    kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_3_3_pt0_80_BDT.weights.xml", 0.05,  2);
//    kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_4_4_pt0_80_BDT.weights.xml", 0.1,   3);
//    kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_5_5_pt0_80_BDT.weights.xml", 0.1,   4);
//    kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_6_6_pt0_80_BDT.weights.xml", 0.125, 5);
//    kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_7_7_pt0_80_BDT.weights.xml", 0.125, 6);
//    kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_8_8_pt0_80_BDT.weights.xml", 0.125, 7);
    
//    kfpAnalysis->RunCentralityAnalysis();
//    kfpAnalysis->SetCentralityFile("/gpfs01/star/pwg/mzyzak/Femto/Template/Centrality/centrality_2016.txt");
  }
  if(year == 2014)
  {
/*
    kfpAnalysis->UseTMVA();
    kfpAnalysis->SetTMVAcutsD0(   "/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/D0.xml",    0.1);
    kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/DPlus.xml", 0.05);
    kfpAnalysis->SetTMVAcutsDs(   "/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/Ds.xml",    0.075);
    kfpAnalysis->SetTMVAcutsLc(   "/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/Lc.xml",    0.1);
    kfpAnalysis->SetTMVAcutsD0KK( "/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/D0KK.xml",  0.125);
    kfpAnalysis->SetTMVAcutsD04(  "/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/D04.xml",  -0.05);
    kfpAnalysis->SetTMVAcutsBPlus("/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/BPlus.xml",-0.1);
    kfpAnalysis->SetTMVAcutsB0(   "/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/B0.xml",   -0.1);
  */
//     kfpAnalysis->RunCentralityAnalysis();
//     kfpAnalysis->SetCentralityFile("/gpfs01/star/pwg/mzyzak/Femto/Template/Centrality/centrality_2014.txt");
  }


    kfpAnalysis->SetDataSet(3);

//    kfpAnalysis->RunFlowAnalysis();
//    kfpAnalysis->AddFlowFile("test_v2.root");

//  kfpAnalysis->CollectPIDHistograms();
//  kfpAnalysis->CollectTrackHistograms();
//
 
//DO NOT STORE TMVA NTUPLES: THIS IS AN IMPORTANT BOOL SETTER FOR CUTS
//  kfpAnalysis->StoreTMVANtuples();
  chain->Init();


//  StKFParticleInterface::instance()->CleanLowPVTrackEvents();
  StKFParticleInterface::instance()->SetSoftKaonPIDMode();
  StKFParticleInterface::instance()->SetSoftTofPidMode();

//default


    //loose

    StKFParticleInterface::instance()->SetMaxDistanceBetweenParticlesCut(5);
    StKFParticleInterface::instance()->SetLCut(5);
    StKFParticleInterface::instance()->SetChi2Cut2D(20);
    StKFParticleInterface::instance()->SetChiPrimaryCut(3);
    StKFParticleInterface::instance()->SetChiPrimaryCut2D(3);
    StKFParticleInterface::instance()->SetLdLCut2D(5);


  
  
//  StKFParticleInterface::instance()->SetPtCutCharm(0.2);
//  StKFParticleInterface::instance()->SetChiPrimaryCutCharm(8);
//  StKFParticleInterface::instance()->SetLdLCutCharmManybodyDecays(3);
//  StKFParticleInterface::instance()->SetChi2TopoCutCharmManybodyDecays(10);
//  StKFParticleInterface::instance()->SetChi2CutCharmManybodyDecays(3);
//  StKFParticleInterface::instance()->SetLdLCutCharm2D(3);
//  StKFParticleInterface::instance()->SetChi2TopoCutCharm2D(10);
//  StKFParticleInterface::instance()->SetChi2CutCharm2D(3);

  //Add decays to the reconstruction list
    StKFParticleInterface::instance()->AddDecayToReconstructionList(  3004);
    StKFParticleInterface::instance()->AddDecayToReconstructionList(  3006);

//  StKFParticleInterface::instance()->AddDecayToReconstructionList( 3122);//lambda0
//  StKFParticleInterface::instance()->AddDecayToReconstructionList(-3122);//

//  StKFParticleInterface::instance()->AddDecayToReconstructionList( 3312);//cascade
//  StKFParticleInterface::instance()->AddDecayToReconstructionList(-3312);

//  StKFParticleInterface::instance()->AddDecayToReconstructionList( 3334);//omega baryon 
//  StKFParticleInterface::instance()->AddDecayToReconstructionList(-3334);

    StKFParticleInterface::instance()->SetRotation();
    StKFParticleInterface::instance()->SetRotationAngle(TMath::Pi()*1.25);
    StKFParticleInterface::instance()->SetRotationPID(1000020030);//he3

  Long64_t nevent = N;
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
