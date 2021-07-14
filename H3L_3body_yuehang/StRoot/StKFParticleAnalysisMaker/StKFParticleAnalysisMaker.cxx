//*-- Author : Yuri Fisyak 02/02/2016
#include "StKFParticleAnalysisMaker.h"
#include "TDirectory.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TChain.h"
#include "TNtuple.h"
#include "TSystem.h"
//--- KF particle classes ---
#include "KFVertex.h"
#include "KFParticle.h"
#include "KFParticleSIMD.h"
#include "KFPTrack.h"
#include "KFParticleTopoReconstructor.h"
#include "StKFParticleInterface.h"
#include "StKFParticlePerformanceInterface.h"
//--- Pico classes ---
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
//--- Mu classes ---
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuMcVertex.h"
#include "StMuDSTMaker/COMMON/StMuMcTrack.h"
//--- TMVA classes ---
#include "TMVA/GeneticAlgorithm.h"
#include "TMVA/GeneticFitter.h"
#include "TMVA/IFitterTarget.h"
#include "TMVA/Factory.h"
//--- StRefMult class ---
// #include "../StRefMultCorr/StRefMultCorr.h"
// #include "../StRefMultCorr/CentralityMaker.h"
// #include "StPileupUtil/StPileupUtil.h"
//-----analysis cuts--
#include "StCuts.h"

ClassImp(StKFParticleAnalysisMaker);

//________________________________________________________________________________
StKFParticleAnalysisMaker::StKFParticleAnalysisMaker(const char *name) : StMaker(name), fNTrackTMVACuts(0), fIsPicoAnalysis(true), fsnn(3), fdEdXMode(1), 
  fStoreTmvaNTuples(false), fProcessSignal(false), fCollectTrackHistograms(false), fCollectPIDHistograms(false),fTMVAselection(false), fStoremctree(false),
  fFlowAnalysis(false), fFlowChain(NULL), fFlowRunId(-1), fFlowEventId(-1), fCentrality(-1), fFlowFiles(), fFlowMap(), 
  fRunCentralityAnalysis(1), fMixEvent(false),fRefmultCorrUtil(0), fCentralityFile(""), fAnalyseDsPhiPi(false)
{
  memset(mBeg,0,mEnd-mBeg+1);

  /*  
      fNTuplePDG[0] = 421;
      fNTuplePDG[1] = 411;
      fNTuplePDG[2] = 431;
      fNTuplePDG[3] = 4122;
      fNTuplePDG[4] = 426;
      fNTuplePDG[5] = 429;
      fNTuplePDG[6] = 521;
      fNTuplePDG[7] = 511;
      fNTuplePDG[8] = 3122;
      */
  fNTuplePDG[0] = 3122;
  //fNTuplePDG[1] = 3334;
  //TODO
  fNTuplePDG[1] = 3334;
  fNTuplePDG[2] = 3312;

  /*
     fNtupleNames[0] = "D0"; 
     fNtupleNames[1] = "DPlus"; 
     fNtupleNames[2] = "Ds"; 
     fNtupleNames[3] = "Lc";
     fNtupleNames[4] = "D0KK";
     fNtupleNames[5] = "D04";
     fNtupleNames[6] = "BPlus";
     fNtupleNames[7] = "B0";
     fNtupleNames[8] = "Ld";
     */
  fNtupleNames[0] = "Ld";
  fNtupleNames[1] = "Om";
  fNtupleNames[2] = "Xi";

  vector<TString> trackCutNames;
  trackCutNames.push_back("pt_");
  trackCutNames.push_back("chi2Primary_");
  trackCutNames.push_back("dEdXPi_");
  trackCutNames.push_back("dEdXK_");
  trackCutNames.push_back("dEdXP_");
  trackCutNames.push_back("ToFPi_");
  trackCutNames.push_back("ToFK_");
  trackCutNames.push_back("ToFP_");
  fNTrackTMVACuts = trackCutNames.size();

  /*  
      fDaughterNames[0].push_back("K");     fDaughterNames[0].push_back("Pi");                                                                              //D0 -> Kpi
      fDaughterNames[1].push_back("K");     fDaughterNames[1].push_back("Pi1");    fDaughterNames[1].push_back("Pi2");                                      //D+ -> Kpipi
      fDaughterNames[2].push_back("KPlus"); fDaughterNames[2].push_back("KMinus"); fDaughterNames[2].push_back("Pi");                                       //Ds -> KKpi
      fDaughterNames[3].push_back("K");     fDaughterNames[3].push_back("Pi");     fDaughterNames[3].push_back("P");                                        //Lc -> pKpi
      fDaughterNames[4].push_back("KPlus"); fDaughterNames[4].push_back("KMinus");                                                                          //D0 -> KK
      fDaughterNames[5].push_back("K");     fDaughterNames[5].push_back("Pi1");    fDaughterNames[5].push_back("Pi2");  fDaughterNames[5].push_back("Pi3"); //D0 -> Kpipipi
      fDaughterNames[6].push_back("PiD");   fDaughterNames[6].push_back("KD");     fDaughterNames[6].push_back("Pi");                                       //B+ -> D0_bpi
      fDaughterNames[7].push_back("Pi1D");  fDaughterNames[7].push_back("KD");     fDaughterNames[7].push_back("Pi2D"); fDaughterNames[7].push_back("Pi");  //B0 -> D-pi+
  //  fDaughterNames[8].push_back("Proton");fDaughterNames[8].push_back("Pi");
  */
  fDaughterNames[0].push_back("Proton");fDaughterNames[0].push_back("Pi");
  fDaughterNames[1].push_back("Kaon");fDaughterNames[1].push_back("Proton");fDaughterNames[1].push_back("Pion");
  fDaughterNames[2].push_back("Pion_bach");fDaughterNames[2].push_back("Proton");fDaughterNames[2].push_back("Pion");

  for(int iDecay=0; iDecay<fNNTuples; iDecay++)
  {
    for(unsigned int iDaughter=0; iDaughter<fDaughterNames[iDecay].size(); iDaughter++)
    {
      for(int iTrackTMVACut=0; iTrackTMVACut<fNTrackTMVACuts; iTrackTMVACut++)
      {
        if(iDaughter==0 && iTrackTMVACut==0)
          fNtupleCutNames[iDecay] = trackCutNames[iTrackTMVACut];  
        else
          fNtupleCutNames[iDecay] += trackCutNames[iTrackTMVACut];
        fNtupleCutNames[iDecay] += fDaughterNames[iDecay][iDaughter];
        fNtupleCutNames[iDecay] += ":";
      }
    }
    if(iDecay==0) 
      fNtupleCutNames[iDecay] += "Chi2NDF:LdL:L:Chi2Topo:refMult:refMultCor:reweight:cent9:mass:px:py:pz:vx:vy:vz:evx:evy:evz";	
    else if(iDecay==1 || iDecay==2)
      //fNtupleCutNames[iDecay] += "Chi2NDF:LdL:L:Chi2Topo:refMult:mass:pt";
      fNtupleCutNames[iDecay] += "Chi2NDF:LdL:L:Chi2Topo:refMult:refMultCor:reweight:cent9:mass:px:py:pz:vx:vy:vz:evx:evy:evz:Chi2NDF_ld:LdL_ld:L_ld:Chi2Topo_ld:pid";
    else if(iDecay>2 && iDecay<6 )
      fNtupleCutNames[iDecay] += "Chi2NDF:LdL:Chi2Topo:refMult";
    else if(iDecay>=6 && iDecay<8)
    {
      fNtupleCutNames[iDecay] += "Chi2NDF_D:LdL_D:Chi2Topo_D:Chi2NDF:LdL:Chi2Topo:refMult";
    } 

    SetTMVABins(iDecay);
  }
}
//________________________________________________________________________________
StKFParticleAnalysisMaker::~StKFParticleAnalysisMaker() 
{
  SafeDelete(fStKFParticleInterface);
  SafeDelete(fStKFParticlePerformanceInterface);
}

//_____________________________________________________________________________
Int_t StKFParticleAnalysisMaker::Init()
{
  TFile *f = GetTFile();
  if(f) 
  {
    f->cd();
    BookVertexPlots();
    if(fCollectTrackHistograms)
      fStKFParticleInterface->CollectTrackHistograms();
    if(fCollectPIDHistograms)
      fStKFParticleInterface->CollectPIDHistograms();
  }

  if(fTMVAselection || fStoreTmvaNTuples)
  {
    for(int iReader=0; iReader<fNNTuples; iReader++)
    {
      TString cutName;
      int firstSymbolOfCutName = 0;

      int nCuts = 0;
      while(fNtupleCutNames[iReader].Tokenize(cutName,firstSymbolOfCutName,":"))
        nCuts++;
      fTMVAParticleParameters[iReader].resize(nCuts);
    }
  }

  if(fTMVAselection)
  {
    for(int iReader=0; iReader<fNNTuples; iReader++)
    {
      const int nCentralityBins = fTMVACentralityBins[iReader].size() - 1;
      const int nPtBins = fTMVAPtBins[iReader].size() - 1;

      for(int iCentralityBin=0; iCentralityBin<nCentralityBins; iCentralityBin++)
      {
        for(int iPtBin=0; iPtBin<nPtBins; iPtBin++)
        {
          fTMVAReader[iReader][iCentralityBin][iPtBin] = new TMVA::Reader("Silent");

          TString cutName;
          int firstSymbolOfCutName = 0;      
          unsigned int iCut = 0;
          while(fNtupleCutNames[iReader].Tokenize(cutName,firstSymbolOfCutName,":"))
          {
            fTMVAReader[iReader][iCentralityBin][iPtBin] -> AddVariable( cutName.Data(), &fTMVAParticleParameters[iReader][iCut] );
            iCut++;
            if(iCut == (fTMVAParticleParameters[iReader].size()-1)) break;
          }

          fTMVAReader[iReader][iCentralityBin][iPtBin] -> BookMVA("BDT", fTMVACutFile[iReader][iCentralityBin][iPtBin].Data());
        }
      }
    }
  }

  //Create file with NTuples for cut optimization
  if(fStoreTmvaNTuples)
  {  
    TFile* curFile = gFile;
    TDirectory* curDirectory = gDirectory;
    for(int iNtuple=0; iNtuple<fNNTuples; iNtuple++)
    {
      TString SignalPrefix = "_Signal";
      if(!fProcessSignal) SignalPrefix = "_BG";
      TString currentNTupleFileName = fNtupleNames[iNtuple]+SignalPrefix+TString(".root");
      fNTupleFile[iNtuple] = new TFile(currentNTupleFileName.Data(),"RECREATE");
      fCutsNTuple[iNtuple] = new TNtuple(fNtupleNames[iNtuple].Data(), fNtupleNames[iNtuple].Data(), fNtupleCutNames[iNtuple].Data());
    }
    gFile = curFile;
    gDirectory = curDirectory;
  }

  // fRefmultCorrUtil = CentralityMaker::instance()->getRefMultCorr() ;

  //Initialise the chain with files containing centrality and reaction plane
  if(fFlowAnalysis)
  {
    std::cout << "StKFParticleAnalysisMaker: run flow analysis. Flow file list:"<<std::endl;

    //fFlowChain = new TChain("mTree");
    fFlowChain = new TChain("psi_tree");
    for(unsigned int iFlowFile=0; iFlowFile<fFlowFiles.size(); iFlowFile++)
    {
      std::cout << "      " << fFlowFiles[iFlowFile] << std::endl;
      fFlowChain->Add(fFlowFiles[iFlowFile].Data());
    }

    fFlowChain->SetBranchStatus("*",0);
    //fFlowChain->SetBranchAddress("runid",   &fFlowRunId);   fFlowChain->SetBranchStatus("runid", 1);
    fFlowChain->SetBranchAddress("runnumber",   &fFlowRunId);   fFlowChain->SetBranchStatus("runnumber", 1);
    //fFlowChain->SetBranchAddress("eventid", &fFlowEventId); fFlowChain->SetBranchStatus("eventid", 1);
    fFlowChain->SetBranchAddress("eventid", &fFlowEventId); fFlowChain->SetBranchStatus("eventid", 1);
    //fFlowChain->SetBranchAddress("cent", &fCentrality);  fFlowChain->SetBranchStatus("cent", 1);
    fFlowChain->SetBranchAddress("centnumber", &fCentrality);  fFlowChain->SetBranchStatus("centnumber", 1);
    fFlowChain->SetBranchAddress("FXTMult", &FXTMult); fFlowChain->SetBranchStatus("FXTMult", 1);

    fFlowChain->SetBranchAddress("psi_1_EPD_0", &psi_1_EPD_0);
    fFlowChain->SetBranchStatus("psi_1_EPD_0", 1);
    fFlowChain->SetBranchAddress("psi_1_EPD_1", &psi_1_EPD_1);
    fFlowChain->SetBranchStatus("psi_1_EPD_1", 1);
    fFlowChain->SetBranchAddress("psi_1_EPD_2", &psi_1_EPD_2);
    fFlowChain->SetBranchStatus("psi_1_EPD_2", 1);
    fFlowChain->SetBranchAddress("psi_1_EPD_3", &psi_1_EPD_3);
    fFlowChain->SetBranchStatus("psi_1_EPD_3", 1);
    fFlowChain->SetBranchAddress("psi_1_EPD_4", &psi_1_EPD_4);
    fFlowChain->SetBranchStatus("psi_1_EPD_4", 1);
    fFlowChain->SetBranchAddress("psi_1_EPD_5", &psi_1_EPD_5);
    fFlowChain->SetBranchStatus("psi_1_EPD_5", 1);
    fFlowChain->SetBranchAddress("psi_1_EPD_6", &psi_1_EPD_6);
    fFlowChain->SetBranchStatus("psi_1_EPD_6", 1);
    fFlowChain->SetBranchAddress("psi_1_EPD", &psi_1_EPD);
    fFlowChain->SetBranchStatus("psi_1_EPD", 1);
    fFlowChain->SetBranchAddress("gweight", &gweight);
    fFlowChain->SetBranchStatus("gweight", 1);

    std::cout << "StKFParticleAnalysisMaker: number of entries in the flow chain" << fFlowChain->GetEntries() << std::endl;
    for(int iEntry=0; iEntry<fFlowChain->GetEntries(); iEntry++)
    {
      fFlowChain->GetEvent(iEntry);
      fFlowMap[GetUniqueEventId(fFlowRunId, fFlowEventId)] = iEntry;
    }
  }

  //for 3 GeV centrality
  // mPileupTool = new StPileupUtil();
  // mPileupTool->init();

  file_out = new TFile("ana_tree.root","RECREATE");
  lambda_tree = new TTree("lambda_tree","ana_tree");

  //  lambda_tree->Branch("bVz",&bVz,"bVz/F");
  //  lambda_tree->Branch("bVr",&bVr,"bVr/F");
  //  lambda_tree->Branch("bVzerr",&bVzerr,"bVzerr/F");
  //  lambda_tree->Branch("bVrerr",&bVrerr,"bVrerr/F");

  lambda_tree->Branch("brunId",&brunid,"brunid/I");
  lambda_tree->Branch("beventid",&beventid,"beventid/I");
  lambda_tree->Branch("bVz",&bVz,"bVz/F");
  lambda_tree->Branch("brefmult",&brefmult,"brefmult/I");
  lambda_tree->Branch("btofmult",&btofmult,"btofmult/I");
  lambda_tree->Branch("countrefmult",&countrefmult,"countrefmult/I");

  lambda_tree->Branch("ld_chi2topo",&ld_chi2topo,"ld_chi2topo/F");
  lambda_tree->Branch("ld_chi2primary",&ld_chi2primary,"ld_chi2primary/F");
  lambda_tree->Branch("ld_chi2ndf",&ld_chi2ndf,"ld_chi2ndf/F");
  lambda_tree->Branch("ld_ldl",&ld_ldl,"ld_ldl/F");
  lambda_tree->Branch("ld_l",&ld_l,"ld_l/F");
  lambda_tree->Branch("ld_dl",&ld_dl,"ld_dl/F");

  lambda_tree->Branch("dca_proton",&dca_proton,"dca_proton/F");
  lambda_tree->Branch("dca_pi",&dca_pi,"dca_pi/F");    

  lambda_tree->Branch("chi2primary_proton",&chi2primary_proton,"chi2primary_proton/F");
  lambda_tree->Branch("chi2primary_pi",&chi2primary_pi,"chi2primary_pi/F");

  lambda_tree->Branch("nhits_ld_proton",&nhits_ld_proton,"nhits_ld_proton/I");
  lambda_tree->Branch("nhits_ld_pi",&nhits_ld_pi,"nhits_ld_pi/I");

  lambda_tree->Branch("bparticleid",&bparticleid,"bparticleid/I");
  lambda_tree->Branch("bparticlemass",&bparticlemass,"bparticlemass/F");

  lambda_tree->Branch("bx",&bx,"bx/F");
  lambda_tree->Branch("by",&by,"by/F");
  lambda_tree->Branch("bz",&bz,"bz/F");
  lambda_tree->Branch("bpx",&bpx,"bpx/F");
  lambda_tree->Branch("bpy",&bpy,"bpy/F");
  lambda_tree->Branch("bpz",&bpz,"bpz/F");
  lambda_tree->Branch("bpl",&bpl,"bpl/F");

  lambda_tree->Branch("notbadrun",&notbadrun,"notbadrun/I");
  lambda_tree->Branch("FXTMult2",&FXTMult2,"FXTMult2/I"); 
  lambda_tree->Branch("bmcpx",&bmcpx,"bmcpx/F");
  lambda_tree->Branch("bmcpy",&bmcpy,"bmcpy/F");
  lambda_tree->Branch("bmcpz",&bmcpz,"bmcpz/F");
  lambda_tree->Branch("bmcx",&bmcx,"bmcx/F");
  lambda_tree->Branch("bmcy",&bmcy,"bmcy/F");
  lambda_tree->Branch("bmcz",&bmcz,"bmcz/F");
  lambda_tree->Branch("bmcl",&bmcl,"bmcl/F");
  lambda_tree->Branch("bmcpl",&bmcpl,"bmcpl/F");

  lambda_tree->Branch("ld_bdfvtx",&ld_bdfvtx,"ld_bdfvtx/F");
  lambda_tree->Branch("ld_bdfvtx2",&ld_bdfvtx2,"ld_bdfvtx2/F");

  lambda_tree->Branch("ld_bdfvtx_xy",&ld_bdfvtx_xy,"ld_bdfvtx_xy/F");
  lambda_tree->Branch("ld_bdfvtxdev_xy",&ld_bdfvtxdev_xy,"ld_bdfvtxdev_xy/F");
  lambda_tree->Branch("ld_lifetime",&ld_lifetime,"ld_lifetime/F");

  lambda_tree->Branch("bismc",&bismc,"bismc/I");

  lambda_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
  lambda_tree->Branch("reweight",&reweight,"reweight/D");
  lambda_tree->Branch("cent9",&cent9,"cent9/I");


  if(fFlowAnalysis){
    lambda_tree->Branch("fCentrality",&fCentrality,"fCentrality/I");
    lambda_tree->Branch("psi_1_EPD_0",&psi_1_EPD_0,"psi_1_EPD_0/D");
    lambda_tree->Branch("psi_1_EPD_1",&psi_1_EPD_1,"psi_1_EPD_1/D");
    lambda_tree->Branch("psi_1_EPD_2",&psi_1_EPD_2,"psi_1_EPD_2/D");
    lambda_tree->Branch("psi_1_EPD_3",&psi_1_EPD_3,"psi_1_EPD_3/D");
    lambda_tree->Branch("psi_1_EPD_4",&psi_1_EPD_4,"psi_1_EPD_4/D");
    lambda_tree->Branch("psi_1_EPD_5",&psi_1_EPD_5,"psi_1_EPD_5/D");
    lambda_tree->Branch("psi_1_EPD_6",&psi_1_EPD_6,"psi_1_EPD_6/D");
    lambda_tree->Branch("psi_1_EPD",&psi_1_EPD,"psi_1_EPD/D");
    lambda_tree->Branch("gweight",&gweight,"gweight/F");
    lambda_tree->Branch("FXTMult",&FXTMult,"FXTMult/I");
  }

  lambda_mc_tree = new TTree("lambda_mc_tree","ana_tree");
  lambda_mc_tree->Branch("brunid",&brunid,"brunid/I");
  lambda_mc_tree->Branch("brefmult",&brefmult,"brefmult/I");
  lambda_mc_tree->Branch("bmcparticleid",&bmcparticleid,"bmcparticleid/I");
  lambda_mc_tree->Branch("bmcrawpx",&bmcrawpx,"bmcrawpx/F");
  lambda_mc_tree->Branch("bmcrawpy",&bmcrawpy,"bmcrawpy/F");
  lambda_mc_tree->Branch("bmcrawpz",&bmcrawpz,"bmcrawpz/F");
  lambda_mc_tree->Branch("bmcrawl",&bmcrawl,"bmcrawl/F");
  lambda_mc_tree->Branch("bmcrawpl",&bmcrawpl,"bmcrawpl/F");

  lambda_mc_tree->Branch("countrefmult",&countrefmult,"countrefmult/I");
  lambda_mc_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
  lambda_mc_tree->Branch("reweight",&reweight,"reweight/D");
  lambda_mc_tree->Branch("cent9",&cent9,"cent9/I");

  lambda_mc_tree->Branch("FXTMult2",&FXTMult2,"FXTMult2/I");

  ks_mc_tree = new TTree("ks_mc_tree","ana_tree");
  ks_mc_tree->Branch("brunid",&brunid,"brunid/I");
  ks_mc_tree->Branch("brefmult",&brefmult,"brefmult/I");
  ks_mc_tree->Branch("bmcparticleid",&bmcparticleid,"bmcparticleid/I");
  ks_mc_tree->Branch("bmcrawpx",&bmcrawpx,"bmcrawpx/F");
  ks_mc_tree->Branch("bmcrawpy",&bmcrawpy,"bmcrawpy/F");
  ks_mc_tree->Branch("bmcrawpz",&bmcrawpz,"bmcrawpz/F");
  ks_mc_tree->Branch("bmcrawl",&bmcrawl,"bmcrawl/F");
  ks_mc_tree->Branch("bmcrawpl",&bmcrawpl,"bmcrawpl/F");

  ks_mc_tree->Branch("countrefmult",&countrefmult,"countrefmult/I");
  ks_mc_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
  ks_mc_tree->Branch("reweight",&reweight,"reweight/D");
  ks_mc_tree->Branch("cent9",&cent9,"cent9/I");

  htriton_mc_tree = new TTree("htriton_mc_tree","ana_tree");
  htriton_mc_tree->Branch("brunid",&brunid,"brunid/I");
  htriton_mc_tree->Branch("brefmult",&brefmult,"brefmult/I");
  htriton_mc_tree->Branch("bmcparticleid",&bmcparticleid,"bmcparticleid/I");
  htriton_mc_tree->Branch("bmcrawpx",&bmcrawpx,"bmcrawpx/F");
  htriton_mc_tree->Branch("bmcrawpy",&bmcrawpy,"bmcrawpy/F");
  htriton_mc_tree->Branch("bmcrawpz",&bmcrawpz,"bmcrawpz/F");
  htriton_mc_tree->Branch("bmcrawl",&bmcrawl,"bmcrawl/F");
  htriton_mc_tree->Branch("bmcrawpl",&bmcrawpl,"bmcrawpl/F");

  htriton_mc_tree->Branch("b0mcrawpx",&b0mcrawpx,"b0mcrawpx/F");
  htriton_mc_tree->Branch("b0mcrawpy",&b0mcrawpy,"b0mcrawpy/F");
  htriton_mc_tree->Branch("b0mcrawpz",&b0mcrawpz,"b0mcrawpz/F");
  htriton_mc_tree->Branch("b1mcrawpx",&b1mcrawpx,"b1mcrawpx/F");
  htriton_mc_tree->Branch("b1mcrawpy",&b1mcrawpy,"b1mcrawpy/F");
  htriton_mc_tree->Branch("b1mcrawpz",&b1mcrawpz,"b1mcrawpz/F");
  htriton_mc_tree->Branch("b2mcrawpx",&b2mcrawpx,"b2mcrawpx/F");
  htriton_mc_tree->Branch("b2mcrawpy",&b2mcrawpy,"b2mcrawpy/F");
  htriton_mc_tree->Branch("b2mcrawpz",&b2mcrawpz,"b2mcrawpz/F");

  htriton_mc_tree->Branch("countrefmult",&countrefmult,"countrefmult/I");
  htriton_mc_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
  htriton_mc_tree->Branch("reweight",&reweight,"reweight/D");
  htriton_mc_tree->Branch("cent9",&cent9,"cent9/I");

  htriton_mc_tree->Branch("FXTMult2",&FXTMult2,"FXTMult2/I");

  ks_tree = new TTree("ks_tree","ana_tree");
  ks_tree->Branch("brunId",&brunid,"brunid/I");
  ks_tree->Branch("beventid",&beventid,"beventid/I");
  ks_tree->Branch("bVz",&bVz,"bVz/F");
  ks_tree->Branch("brefmult",&brefmult,"brefmult/I");
  ks_tree->Branch("btofmult",&btofmult,"btofmult/I");
  ks_tree->Branch("ld_chi2topo",&ld_chi2topo,"ld_chi2topo/F");
  ks_tree->Branch("ld_chi2ndf",&ld_chi2ndf,"ld_chi2ndf/F");
  ks_tree->Branch("ld_ldl",&ld_ldl,"ld_ldl/F");
  ks_tree->Branch("ld_l",&ld_l,"ld_l/F");
  ks_tree->Branch("ld_dl",&ld_dl,"ld_dl/F");
  ks_tree->Branch("chi2primary_proton",&chi2primary_proton,"chi2primary_proton/F");
  ks_tree->Branch("chi2primary_pi",&chi2primary_pi,"chi2primary_pi/F");

  ks_tree->Branch("nhits_ld_proton",&nhits_ld_proton,"nhits_ld_proton/I");
  ks_tree->Branch("nhits_ld_pi",&nhits_ld_pi,"nhits_ld_pi/I");

  ks_tree->Branch("bparticleid",&bparticleid,"bparticleid/I");
  ks_tree->Branch("bparticlemass",&bparticlemass,"bparticlemass/F");
  ks_tree->Branch("bpx",&bpx,"bpx/F");
  ks_tree->Branch("bpy",&bpy,"bpy/F");
  ks_tree->Branch("bpz",&bpz,"bpz/F");
  ks_tree->Branch("bmcpx",&bmcpx,"bmcpx/F");
  ks_tree->Branch("bmcpy",&bmcpy,"bmcpy/F");
  ks_tree->Branch("bmcpz",&bmcpz,"bmcpz/F");
  ks_tree->Branch("bismc",&bismc,"bismc/I");
  ks_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
  ks_tree->Branch("reweight",&reweight,"reweight/D");
  ks_tree->Branch("cent9",&cent9,"cent9/I");

  ks_tree->Branch("bx",&bx,"bx/F");
  ks_tree->Branch("by",&by,"by/F");
  ks_tree->Branch("bz",&bz,"bz/F");

  ks_tree->Branch("bmcpx",&bmcpx,"bmcpx/F");
  ks_tree->Branch("bmcpy",&bmcpy,"bmcpy/F");
  ks_tree->Branch("bmcpz",&bmcpz,"bmcpz/F");
  ks_tree->Branch("bmcx",&bmcx,"bmcx/F");
  ks_tree->Branch("bmcy",&bmcy,"bmcy/F");
  ks_tree->Branch("bmcz",&bmcz,"bmcz/F");

  if(fFlowAnalysis){
    ks_tree->Branch("fCentrality",&fCentrality,"fCentrality/I");
    ks_tree->Branch("psi_1_EPD_0",&psi_1_EPD_0,"psi_1_EPD_0/D");
    ks_tree->Branch("psi_1_EPD_1",&psi_1_EPD_1,"psi_1_EPD_1/D");
    ks_tree->Branch("psi_1_EPD_2",&psi_1_EPD_2,"psi_1_EPD_2/D");
    ks_tree->Branch("psi_1_EPD_3",&psi_1_EPD_3,"psi_1_EPD_3/D");
    ks_tree->Branch("psi_1_EPD_4",&psi_1_EPD_4,"psi_1_EPD_4/D");
    ks_tree->Branch("psi_1_EPD_5",&psi_1_EPD_5,"psi_1_EPD_5/D");
    ks_tree->Branch("psi_1_EPD_6",&psi_1_EPD_6,"psi_1_EPD_6/D");
    ks_tree->Branch("psi_1_EPD",&psi_1_EPD,"psi_1_EPD/D");
    ks_tree->Branch("gweight",&gweight,"gweight/F");
    ks_tree->Branch("FXTMult",&FXTMult,"FXTMult/I");
  }


  cascade_tree = new TTree("cascade_tree","ana_tree");
  cascade_tree->Branch("brunid",&brunid,"brunid/I");
  cascade_tree->Branch("beventid",&beventid,"beventid/I");

  cascade_tree->Branch("bVz",&bVz,"bVz/F");
  //  cascade_tree->Branch("bVr",&bVr,"bVr/F");
  //  cascade_tree->Branch("bVzerr",&bVzerr,"bVzerr/F");
  //  cascade_tree->Branch("bVrerr",&bVrerr,"bVrerr/F");

  cascade_tree->Branch("bparticleid",&bparticleid,"bparticleid/I");
  cascade_tree->Branch("bparticlemass",&bparticlemass,"bparticlemass/F");
  //  cascade_tree->Branch("bx",&bx,"bx/F");
  //  cascade_tree->Branch("by",&by,"by/F");
  //  cascade_tree->Branch("bz",&bz,"bz/F");
  cascade_tree->Branch("bpx",&bpx,"bpx/F");
  cascade_tree->Branch("bpy",&bpy,"bpy/F");
  cascade_tree->Branch("bpz",&bpz,"bpz/F");
  //  cascade_tree->Branch("bdl",&bdl,"bdl/F");

  cascade_tree->Branch("brefmult",&brefmult,"brefmult/I");
  cascade_tree->Branch("btofmult",&btofmult,"btofmult/I");

  cascade_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
  cascade_tree->Branch("reweight",&reweight,"reweight/D");
  cascade_tree->Branch("cent9",&cent9,"cent9/I");

  //  cascade_tree->Branch("bbachid",&bbachid,"bbachid/I");
  //  cascade_tree->Branch("bbachpx",&bbachpx,"bbachpx/F");
  //  cascade_tree->Branch("bbachpy",&bbachpy,"bbachpy/F");
  //  cascade_tree->Branch("bbachpz",&bbachpz,"bbachpz/F");
  //  cascade_tree->Branch("bbachmass",&bbachmass,"bbachmass/F");

  //  cascade_tree->Branch("bpionid",&bpionid,"bpionid/I");
  //  cascade_tree->Branch("bpionpx",&bpionpx,"bpionpx/F");
  //  cascade_tree->Branch("bpionpy",&bpionpy,"bpionpy/F");
  //  cascade_tree->Branch("bpionpz",&bpionpz,"bpionpz/F");
  //  cascade_tree->Branch("bpionmass",&bpionmass,"bpionmass/F");

  //  cascade_tree->Branch("bprotonid",&bprotonid,"bprotonid/I");
  //  cascade_tree->Branch("bprotonpx",&bprotonpx,"bprotonpx/F");
  //  cascade_tree->Branch("bprotonpy",&bprotonpy,"bprotonpy/F");
  //  cascade_tree->Branch("bprotonpz",&bprotonpz,"bprotonpz/F");
  //  cascade_tree->Branch("bprotonmass",&bprotonmass,"bprotonmass/F");

  cascade_tree->Branch("xi_chi2topo",&xi_chi2topo,"xi_chi2topo/F");
  cascade_tree->Branch("xi_chi2ndf",&xi_chi2ndf,"xi_chi2ndf/F");
  cascade_tree->Branch("xi_ldl",&xi_ldl,"xi_ldl/F");

  cascade_tree->Branch("xi_ld_chi2topo",&xi_ld_chi2topo,"xi_ld_chi2topo/F");
  cascade_tree->Branch("xi_ld_chi2ndf",&xi_ld_chi2ndf,"xi_ld_chi2ndf/F");
  cascade_tree->Branch("xi_ld_ldl",&xi_ld_ldl,"xi_ld_ldl/F");
  cascade_tree->Branch("xi_ld_l",&xi_ld_l,"xi_ld_l/F");

  cascade_tree->Branch("xi_l",&xi_l,"xi_l/F");
  //cascade_tree->Branch("xi_dl",&xi_dl,"xi_dl/F");

  cascade_tree->Branch("chi2primary_xi_proton",&chi2primary_xi_proton,"chi2primary_xi_proton/F");
  cascade_tree->Branch("chi2primary_xi_pi",&chi2primary_xi_pi,"chi2primary_xi_pi/F");
  cascade_tree->Branch("chi2primary_xi_bach",&chi2primary_xi_bach,"chi2primary_xi_bach/F");
  cascade_tree->Branch("chi2primary_xi_ld",&chi2primary_xi_ld,"chi2primary_xi_ld/F");

  cascade_tree->Branch("notbadrun",&notbadrun,"notbadrun/I");

  cascade_tree->Branch("bmcpx",&bmcpx,"bmcpx/F");
  cascade_tree->Branch("bmcpy",&bmcpy,"bmcpy/F");
  cascade_tree->Branch("bmcpz",&bmcpz,"bmcpz/F");
  cascade_tree->Branch("bismc",&bismc,"bismc/I");


  cascade_mc_tree = new TTree("cascade_mc_tree","ana_tree");
  cascade_mc_tree->Branch("brunid",&brunid,"brunid/I");
  cascade_mc_tree->Branch("brefmult",&brefmult,"brefmult/I");
  cascade_mc_tree->Branch("bmcparticleid",&bmcparticleid,"bmcparticleid/I");
  cascade_mc_tree->Branch("bmcrawpx",&bmcrawpx,"bmcrawpx/F");
  cascade_mc_tree->Branch("bmcrawpy",&bmcrawpy,"bmcrawpy/F");
  cascade_mc_tree->Branch("bmcrawpz",&bmcrawpz,"bmcrawpz/F");

  cascade_mc_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
  cascade_mc_tree->Branch("reweight",&reweight,"reweight/D");
  cascade_mc_tree->Branch("cent9",&cent9,"cent9/I");

  omega_tree = new TTree("omega_tree","ana_tree");
  omega_tree->Branch("brunid",&brunid,"brunid/I");
  // omega_tree->Branch("beventid",&beventid,"beventid/I");

  omega_tree->Branch("bVz",&bVz,"bVz/F");
  // omega_tree->Branch("bVr",&bVr,"bVr/F");
  // omega_tree->Branch("bVzerr",&bVzerr,"bVzerr/F");
  // omega_tree->Branch("bVrerr",&bVrerr,"bVrerr/F");
  omega_tree->Branch("brefmult",&brefmult,"brefmult/I");
  omega_tree->Branch("btofmult",&btofmult,"btofmult/I");

  omega_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
  omega_tree->Branch("reweight",&reweight,"reweight/D");
  omega_tree->Branch("cent9",&cent9,"cent9/I");

  omega_tree->Branch("bparticleid",&bparticleid,"bparticleid/I");
  omega_tree->Branch("bparticlemass",&bparticlemass,"bparticlemass/F");
  // omega_tree->Branch("bx",&bx,"bx/F");
  // omega_tree->Branch("by",&by,"by/F");
  // omega_tree->Branch("bz",&bz,"bz/F");
  omega_tree->Branch("bpx",&bpx,"bpx/F");
  omega_tree->Branch("bpy",&bpy,"bpy/F");
  omega_tree->Branch("bpz",&bpz,"bpz/F");
  // omega_tree->Branch("bdl",&bdl,"bdl/F");

  omega_tree->Branch("om_chi2topo",&om_chi2topo,"om_chi2topo/F");
  omega_tree->Branch("om_chi2ndf",&om_chi2ndf,"om_chi2ndf/F");
  omega_tree->Branch("om_ldl",&om_ldl,"om_ldl/F");

  omega_tree->Branch("om_ld_chi2topo",&om_ld_chi2topo,"om_ld_chi2topo/F");
  omega_tree->Branch("om_ld_chi2ndf",&om_ld_chi2ndf,"om_ld_chi2ndf/F");
  omega_tree->Branch("om_ld_ldl",&om_ld_ldl,"om_ld_ldl/F");
  omega_tree->Branch("om_ld_l",&om_ld_l,"om_ld_l/F");

  omega_tree->Branch("om_l",&om_l,"om_l/F");
  omega_tree->Branch("om_dl",&om_dl,"om_dl/F");

  omega_tree->Branch("chi2primary_om_proton",&chi2primary_om_proton,"chi2primary_om_proton/F");
  omega_tree->Branch("chi2primary_om_pi",&chi2primary_om_pi,"chi2primary_om_pi/F");
  omega_tree->Branch("chi2primary_om_bach",&chi2primary_om_bach,"chi2primary_om_bach/F");
  omega_tree->Branch("chi2primary_om_ld",&chi2primary_om_ld,"chi2primary_om_ld/F");

  omega_tree->Branch("nhits_om_proton",&nhits_om_proton,"nhits_om_proton/I");
  omega_tree->Branch("nhits_om_pi",&nhits_om_pi,"nhits_om_pi/I");
  omega_tree->Branch("nhits_om_bach",&nhits_om_bach,"nhits_om_bach/I");
  omega_tree->Branch("dedx_om_proton",&dedx_om_proton,"dedx_om_proton/F");
  omega_tree->Branch("dedx_om_pi",&dedx_om_pi,"dedx_om_pi/F");
  omega_tree->Branch("dedx_om_bach",&dedx_om_bach,"dedx_om_bach/F");

  omega_tree->Branch("bbachpx",&bbachpx,"bbachpx/F");
  omega_tree->Branch("bbachpy",&bbachpy,"bbachpy/F");
  omega_tree->Branch("bbachpz",&bbachpz,"bbachpz/F");

  omega_tree->Branch("bprotonpx",&bprotonpx,"bprotonpx/F");
  omega_tree->Branch("bprotonpy",&bprotonpy,"bprotonpy/F");
  omega_tree->Branch("bprotonpz",&bprotonpz,"bprotonpz/F");

  omega_tree->Branch("bpionpx",&bpionpx,"bpionpx/F");
  omega_tree->Branch("bpionpy",&bpionpy,"bpionpy/F");
  omega_tree->Branch("bpionpz",&bpionpz,"bpionpz/F");

  omega_tree->Branch("notbadrun",&notbadrun,"notbadrun/I");

  omega_tree->Branch("bmcpx",&bmcpx,"bmcpx/F");
  omega_tree->Branch("bmcpy",&bmcpy,"bmcpy/F");
  omega_tree->Branch("bmcpz",&bmcpz,"bmcpz/F");
  omega_tree->Branch("bismc",&bismc,"bismc/I");


  omega_mc_tree = new TTree("omega_mc_tree","ana_tree");
  omega_mc_tree->Branch("brunid",&brunid,"brunid/I");
  omega_mc_tree->Branch("brefmult",&brefmult,"brefmult/I");
  omega_mc_tree->Branch("bmcparticleid",&bmcparticleid,"bmcparticleid/I");
  omega_mc_tree->Branch("bmcrawpx",&bmcrawpx,"bmcrawpx/F");
  omega_mc_tree->Branch("bmcrawpy",&bmcrawpy,"bmcrawpy/F");
  omega_mc_tree->Branch("bmcrawpz",&bmcrawpz,"bmcrawpz/F");

  omega_mc_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
  omega_mc_tree->Branch("reweight",&reweight,"reweight/D");
  omega_mc_tree->Branch("cent9",&cent9,"cent9/I");

  /*
     omega_tree->Branch("bbachid",&bbachid,"bbachid/I");
     omega_tree->Branch("bbachpx",&bbachpx,"bbachpx/F");
     omega_tree->Branch("bbachpy",&bbachpy,"bbachpy/F");
     omega_tree->Branch("bbachpz",&bbachpz,"bbachpz/F");
     omega_tree->Branch("bbachmass",&bbachmass,"bbachmass/F");

     omega_tree->Branch("bpionid",&bpionid,"bpionid/I");
     omega_tree->Branch("bpionpx",&bpionpx,"bpionpx/F");
     omega_tree->Branch("bpionpy",&bpionpy,"bpionpy/F");
     omega_tree->Branch("bpionpz",&bpionpz,"bpionpz/F");
     omega_tree->Branch("bpionmass",&bpionmass,"bpionmass/F");

     omega_tree->Branch("bprotonid",&bprotonid,"bprotonid/I");
     omega_tree->Branch("bprotonpx",&bprotonpx,"bprotonpx/F");
     omega_tree->Branch("bprotonpy",&bprotonpy,"bprotonpy/F");
     omega_tree->Branch("bprotonpz",&bprotonpz,"bprotonpz/F");
     omega_tree->Branch("bprotonmass",&bprotonmass,"bprotonmass/F");
     */

  htriton_tree = new TTree("htriton_tree","ana_tree");
  htriton_tree->Branch("brunid",&brunid,"brunid/I");
  htriton_tree->Branch("beventid",&beventid,"beventid/I");

  htriton_tree->Branch("bVz",&bVz,"bVz/F");
  // omega_tree->Branch("bVr",&bVr,"bVr/F");
  htriton_tree->Branch("bVzerr",&bVzerr,"bVzerr/F");
  // omega_tree->Branch("bVrerr",&bVrerr,"bVrerr/F"); 
  htriton_tree->Branch("brefmult",&brefmult,"brefmult/I");
  htriton_tree->Branch("btofmult",&btofmult,"btofmult/I");
  htriton_tree->Branch("countrefmult",&countrefmult,"countrefmult/I");

  htriton_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
  htriton_tree->Branch("reweight",&reweight,"reweight/D");
  htriton_tree->Branch("cent9",&cent9,"cent9/I");

  htriton_tree->Branch("bparticleid",&bparticleid,"bparticleid/I");
  htriton_tree->Branch("bparticlemass",&bparticlemass,"bparticlemass/F");
  // omega_tree->Branch("bx",&bx,"bx/F");
  // omega_tree->Branch("by",&by,"by/F");
  // omega_tree->Branch("bz",&bz,"bz/F");
  htriton_tree->Branch("bpx",&bpx,"bpx/F");
  htriton_tree->Branch("bpy",&bpy,"bpy/F");
  htriton_tree->Branch("bpz",&bpz,"bpz/F");
  htriton_tree->Branch("bpl",&bpl,"bpl/F");

  htriton_tree->Branch("chi2primary_he",&chi2primary_he,"chi2primary_he/F");
  htriton_tree->Branch("chi2primary_pi",&chi2primary_pi,"chi2primary_pi/F");
  htriton_tree->Branch("dca_he",&dca_he,"dca_he/F");
  htriton_tree->Branch("dca_pi",&dca_pi,"dca_pi/F");

  htriton_tree->Branch("ht_ldl",&ht_ldl,"ht_ldl/F");
  htriton_tree->Branch("ht_dl",&ht_dl,"ht_dl/F");
  htriton_tree->Branch("ht_l",&ht_l,"ht_l/F");
  htriton_tree->Branch("ht_chi2topo",&ht_chi2topo,"ht_chi2topo/F");
  htriton_tree->Branch("ht_chi2ndf",&ht_chi2ndf,"ht_chi2ndf/F");


  htriton_tree->Branch("nhits_he",&nhits_he,"nhits_he/I");
  htriton_tree->Branch("nhits_pi",&nhits_pi,"nhits_pi/I");

  htriton_tree->Branch("nhitsdedx_he",&nhitsdedx_he,"nhitsdedx_he/I");
  htriton_tree->Branch("nhitsdedx_pi",&nhitsdedx_pi,"nhitsdedx_pi/I");

  htriton_tree->Branch("dedx_he",&dedx_he,"dedx_he/F");
  htriton_tree->Branch("dedx_pi",&dedx_pi,"dedx_pi/F");

  htriton_tree->Branch("ht_chi2",&ht_chi2,"ht_chi2/F");
  htriton_tree->Branch("ht_NDF",&ht_NDF,"ht_NDF/F");

  htriton_tree->Branch("ht_bdfvtx",&ht_bdfvtx,"ht_bdfvtx/F");
  htriton_tree->Branch("ht_bdfvtx2",&ht_bdfvtx2,"ht_bdfvtx2/F");
  htriton_tree->Branch("ht_lifetime",&ht_lifetime,"ht_lifetime/F");

  htriton_tree->Branch("px_pi",&px_pi,"px_pi/F");
  htriton_tree->Branch("py_pi",&py_pi,"py_pi/F");
  htriton_tree->Branch("pz_pi",&pz_pi,"pz_pi/F");
  htriton_tree->Branch("px_he",&px_he,"px_he/F");
  htriton_tree->Branch("py_he",&py_he,"py_he/F");
  htriton_tree->Branch("pz_he",&pz_he,"pz_he/F");

  htriton_tree->Branch("bismc",&bismc,"bismc/I");
  htriton_tree->Branch("bmcpx",&bmcpx,"bmcpx/F");
  htriton_tree->Branch("bmcpy",&bmcpy,"bmcpy/F");
  htriton_tree->Branch("bmcpz",&bmcpz,"bmcpz/F");
  htriton_tree->Branch("bmcl",&bmcl,"bmcl/F");
  htriton_tree->Branch("bmcpl",&bmcpl,"bmcpl/F");

  htriton_tree->Branch("FXTMult2",&FXTMult2,"FXTMult2/I");

  if(fFlowAnalysis){
    htriton_tree->Branch("fCentrality",&fCentrality,"fCentrality/I");
    htriton_tree->Branch("psi_1_EPD_0",&psi_1_EPD_0,"psi_1_EPD_0/D");
    htriton_tree->Branch("psi_1_EPD_1",&psi_1_EPD_1,"psi_1_EPD_1/D");
    htriton_tree->Branch("psi_1_EPD_2",&psi_1_EPD_2,"psi_1_EPD_2/D");
    htriton_tree->Branch("psi_1_EPD_3",&psi_1_EPD_3,"psi_1_EPD_3/D");
    htriton_tree->Branch("psi_1_EPD_4",&psi_1_EPD_4,"psi_1_EPD_4/D");
    htriton_tree->Branch("psi_1_EPD_5",&psi_1_EPD_5,"psi_1_EPD_5/D");
    htriton_tree->Branch("psi_1_EPD_6",&psi_1_EPD_6,"psi_1_EPD_6/D");
    htriton_tree->Branch("psi_1_EPD",&psi_1_EPD,"psi_1_EPD/D");
    htriton_tree->Branch("gweight",&gweight,"gweight/F");
  }


  htriton3_tree = new TTree("htriton3_tree","ana_tree");
  htriton3_tree->Branch("brunid",&brunid,"brunid/I");
  htriton3_tree->Branch("beventid",&beventid,"beventid/I");
  htriton3_tree->Branch("bVz",&bVz,"bVz/F");
  htriton3_tree->Branch("bVx",&bVx,"bVx/F");
  htriton3_tree->Branch("bVy",&bVy,"bVy/F");
  htriton3_tree->Branch("brefmult",&brefmult,"brefmult/I");
  htriton3_tree->Branch("btofmult",&btofmult,"btofmult/I");
  htriton3_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
  htriton3_tree->Branch("reweight",&reweight,"reweight/D");
  htriton3_tree->Branch("cent9",&cent9,"cent9/I");
  htriton3_tree->Branch("bparticleid",&bparticleid,"bparticleid/I");
  htriton3_tree->Branch("bparticlemass",&bparticlemass,"bparticlemass/F");
  htriton3_tree->Branch("bpx",&bpx,"bpx/F");
  htriton3_tree->Branch("bpy",&bpy,"bpy/F");
  htriton3_tree->Branch("bpz",&bpz,"bpz/F");
  htriton3_tree->Branch("ht_chi2",&ht_chi2,"ht_chi2/F");
  htriton3_tree->Branch("ht_ndaughters",&ht_ndaughters,"ht_ndaughters/I");
  htriton3_tree->Branch("ht_NDF",&ht_NDF,"ht_NDF/F");

  htriton3_tree->Branch("chi2primary_proton",&chi2primary_proton,"chi2primary_proton/F");
  htriton3_tree->Branch("chi2primary_pi",&chi2primary_pi,"chi2primary_pi/F");
  htriton3_tree->Branch("chi2primary_d",&chi2primary_d,"chi2primary_d/F");
  htriton3_tree->Branch("ht_ldl",&ht_ldl,"ht_ldl/F");
  htriton3_tree->Branch("ht_dl",&ht_dl,"ht_dl/F");
  htriton3_tree->Branch("ht_l",&ht_l,"ht_l/F");
  htriton3_tree->Branch("ht_chi2topo",&ht_chi2topo,"ht_chi2topo/F");
  htriton3_tree->Branch("ht_chi2ndf",&ht_chi2ndf,"ht_chi2ndf/F");
  htriton3_tree->Branch("bpionpx",&bpionpx,"bpionpx/F");
  htriton3_tree->Branch("bpionpy",&bpionpy,"bpionpy/F");
  htriton3_tree->Branch("bpionpz",&bpionpz,"bpionpz/F");
  htriton3_tree->Branch("bprotonpx",&bprotonpx,"bprotonpx/F");
  htriton3_tree->Branch("bprotonpy",&bprotonpy,"bprotonpy/F");
  htriton3_tree->Branch("bprotonpz",&bprotonpz,"bprotonpz/F");
  htriton3_tree->Branch("bdpx",&bdpx,"bdpx/F");
  htriton3_tree->Branch("bdpy",&bdpy,"bdpy/F");
  htriton3_tree->Branch("bdpz",&bdpz,"bdpz/F");
  htriton3_tree->Branch("bzdeuteron",&bzdeuteron,"bzdeuteron/F");
  htriton3_tree->Branch("bpionnsigma",&bpionnsigma ,"bpionnsigma/F");
  htriton3_tree->Branch("bprotonsigma",&bprotonsigma ,"bprotonsigma/F");
  htriton3_tree->Branch("dca_pion",&dca_pion ,"dca_pion/F");
  htriton3_tree->Branch("dca_proton",&dca_proton ,"dca_proton/F");
  htriton3_tree->Branch("dca_deuteron",&dca_deuteron ,"dca_deuteron/F");

  htriton3_tree->Branch("dca_pion2",&dca_pion2 ,"dca_pion2/F");
  htriton3_tree->Branch("dca_proton2",&dca_proton2 ,"dca_proton2/F");
  htriton3_tree->Branch("dca_deuteron2",&dca_deuteron2 ,"dca_deuteron2/F");

  htriton3_tree->Branch("nhits_pion",&nhits_pion ,"nhits_pion/I");
  htriton3_tree->Branch("nhits_proton",&nhits_proton ,"nhits_proton/I");
  htriton3_tree->Branch("nhits_deuteron",&nhits_deuteron ,"nhits_deuteron/I");
  htriton3_tree->Branch("countrefmult",&countrefmult,"countrefmult/I");
  htriton3_tree->Branch("bdedx",&bdedx,"bdedx/F");
  //htriton3_tree->Branch("btofm2",&btofm2,"btofm2/F");
  //
  htriton3_tree->Branch("v_01_pvdca",&v_01_pvdca,"v_01_pvdca/F");
  htriton3_tree->Branch("v_02_pvdca",&v_02_pvdca,"v_02_pvdca/F");
  htriton3_tree->Branch("v_12_pvdca",&v_12_pvdca,"v_12_pvdca/F");
  htriton3_tree->Branch("ht_bdfvtx",&ht_bdfvtx,"ht_bdfvtx/F");
  htriton3_tree->Branch("ht_bdfvtx2",&ht_bdfvtx2,"ht_bdfvtx2/F");


  htriton3_tree->Branch("bpionm2",&bpionm2,"bpionm2/F");
  htriton3_tree->Branch("bprotonm2",&bprotonm2,"bprotonm2/F");
  htriton3_tree->Branch("bdm2",&bdm2,"bdm2/F");

  htriton3_tree->Branch("bisMix", &bisMix, "bisMix/I");

  htriton3_tree->Branch("bmcpx",&bmcpx,"bmcpx/F");
  htriton3_tree->Branch("bmcpy",&bmcpy,"bmcpy/F");
  htriton3_tree->Branch("bmcpz",&bmcpz,"bmcpz/F");

  htriton3_tree->Branch("b0mcpx",&b0mcpx,"b0mcpx/F");
  htriton3_tree->Branch("b0mcpy",&b0mcpy,"b0mcpy/F");
  htriton3_tree->Branch("b0mcpz",&b0mcpz,"b0mcpz/F");
  htriton3_tree->Branch("b1mcpx",&b1mcpx,"b1mcpx/F");
  htriton3_tree->Branch("b1mcpy",&b1mcpy,"b1mcpy/F");
  htriton3_tree->Branch("b1mcpz",&b1mcpz,"b1mcpz/F");
  htriton3_tree->Branch("b2mcpx",&b2mcpx,"b2mcpx/F");
  htriton3_tree->Branch("b2mcpy",&b2mcpy,"b2mcpy/F");
  htriton3_tree->Branch("b2mcpz",&b2mcpz,"b2mcpz/F");


  htriton3_tree->Branch("bmcl",&bmcl,"bmcl/F");
  htriton3_tree->Branch("bmcpl",&bmcpl,"bmcpl/F");
  htriton3_tree->Branch("bismc",&bismc,"bismc/I");

  htriton3_tree->Branch("v_lambda_mass_0",&v_lambda_mass_0,"v_lambda_mass_0/F");
  htriton3_tree->Branch("v_lambda_ldl_0",&v_lambda_ldl_0,"v_lambda_ldl_0/F");
  htriton3_tree->Branch("v_lambda_chi2primary_0",&v_lambda_chi2primary_0,"v_lambda_chi2primary_0/F");
  htriton3_tree->Branch("v_lambda_l_0",&v_lambda_l_0,"v_lambda_l_0/F");
  htriton3_tree->Branch("mass_01",&mass_01,"mass_01/F");
  htriton3_tree->Branch("mass_02",&mass_02,"mass_02/F");
  htriton3_tree->Branch("mass_12",&mass_12,"mass_12/F");
  htriton3_tree->Branch("v_01_chi2primary",&v_01_chi2primary,"v_01_chi2primary/F");
  htriton3_tree->Branch("v_02_chi2primary",&v_02_chi2primary,"v_02_chi2primary/F");
  htriton3_tree->Branch("v_12_chi2primary",&v_12_chi2primary,"v_12_chi2primary/F");

  htriton3_tree->Branch("v_01_dca",&v_01_dca,"v_01_dca/F");
  htriton3_tree->Branch("v_02_dca",&v_02_dca,"v_02_dca/F");
  htriton3_tree->Branch("v_12_dca",&v_12_dca,"v_12_dca/F");
  htriton3_tree->Branch("v_012_dca",&v_012_dca,"v_012_dca/F");
  htriton3_tree->Branch("v_01_chi2ndf",&v_01_chi2ndf,"v_01_chi2ndf/F");
  htriton3_tree->Branch("v_02_chi2ndf",&v_02_chi2ndf,"v_02_chi2ndf/F");
  htriton3_tree->Branch("v_12_chi2ndf",&v_12_chi2ndf,"v_12_chi2ndf/F");

  he3_tree = new TTree("he3_tree","ana_tree");
  he3_tree->Branch("brunid",&brunid,"brunid/I");
  he3_tree->Branch("bVz",&bVz,"bVz/F");

  he3_tree->Branch("brefmult",&brefmult,"brefmult/I");
  he3_tree->Branch("btofmult",&btofmult,"btofmult/I");

  he3_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
  he3_tree->Branch("reweight",&reweight,"reweight/D");
  he3_tree->Branch("cent9",&cent9,"cent9/I");

  he3_tree->Branch("bparticleid",&bparticleid,"bparticleid/I");
  he3_tree->Branch("bparticlemass",&bparticlemass,"bparticlemass/F");

  he3_tree->Branch("bpx",&bpx,"bpx/F");
  he3_tree->Branch("bpy",&bpy,"bpy/F");
  he3_tree->Branch("bpz",&bpz,"bpz/F");

  he3_tree->Branch("chi2primary_he",&chi2primary_he,"chi2primary_he/F");
  he3_tree->Branch("bnhits",&bnhits,"bnhits/I");
  he3_tree->Branch("bdedx",&bdedx,"bdedx/F");
  he3_tree->Branch("bdca",&bdca,"bdca/F");

  he3_tree->Branch("bmcpx",&bmcpx,"bmcpx/F");
  he3_tree->Branch("bmcpy",&bmcpy,"bmcpy/F");
  he3_tree->Branch("bmcpz",&bmcpz,"bmcpz/F");
  he3_tree->Branch("bismc",&bismc,"bismc/I");

  h4lambda_tree = new TTree("h4lambda_tree","ana_tree");
  h4lambda_tree->Branch("brunid",&brunid,"brunid/I");
  h4lambda_tree->Branch("beventid",&beventid,"beventid/I");
  h4lambda_tree->Branch("bVz",&bVz,"bVz/F");
  h4lambda_tree->Branch("brefmult",&brefmult,"brefmult/I");
  h4lambda_tree->Branch("btofmult",&btofmult,"btofmult/I");
  h4lambda_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
  h4lambda_tree->Branch("reweight",&reweight,"reweight/D");
  h4lambda_tree->Branch("cent9",&cent9,"cent9/I");
  h4lambda_tree->Branch("bparticleid",&bparticleid,"bparticleid/I");
  h4lambda_tree->Branch("bparticlemass",&bparticlemass,"bparticlemass/F");
  h4lambda_tree->Branch("bpx",&bpx,"bpx/F");
  h4lambda_tree->Branch("bpy",&bpy,"bpy/F");
  h4lambda_tree->Branch("bpz",&bpz,"bpz/F");
  h4lambda_tree->Branch("chi2primary_h4",&chi2primary_h4,"chi2primary_h4/F");
  h4lambda_tree->Branch("chi2primary_pi",&chi2primary_pi,"chi2primary_pi/F");
  h4lambda_tree->Branch("hl_ldl",&hl_ldl,"hl_ldl/F");
  h4lambda_tree->Branch("hl_dl",&hl_dl,"hl_dl/F");
  h4lambda_tree->Branch("hl_l",&hl_l,"hl_l/F");
  h4lambda_tree->Branch("hl_chi2topo",&hl_chi2topo,"hl_chi2topo/F");
  h4lambda_tree->Branch("hl_chi2ndf",&hl_chi2ndf,"hl_chi2ndf/F");
  h4lambda_tree->Branch("hl_chi2",&hl_chi2,"hl_chi2/F");
  h4lambda_tree->Branch("hl_NDF",&hl_NDF,"hl_NDF/F");


  he3_mc_tree = new TTree("he3_mc_tree","ana_tree");
  he3_mc_tree->Branch("brunid",&brunid,"brunid/I");
  he3_mc_tree->Branch("brefmult",&brefmult,"brefmult/I");
  he3_mc_tree->Branch("bmcparticleid",&bmcparticleid,"bmcparticleid/I");
  he3_mc_tree->Branch("bmcrawpx",&bmcrawpx,"bmcrawpx/F");
  he3_mc_tree->Branch("bmcrawpy",&bmcrawpy,"bmcrawpy/F");
  he3_mc_tree->Branch("bmcrawpz",&bmcrawpz,"bmcrawpz/F");
  he3_mc_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
  he3_mc_tree->Branch("reweight",&reweight,"reweight/D");
  he3_mc_tree->Branch("cent9",&cent9,"cent9/I");

  h_tree = new TTree("h_tree","ana_tree");
  h_tree->Branch("brunid",&brunid,"brunid/I");
  h_tree->Branch("beventid",&beventid,"beventid/I");
  h_tree->Branch("bVz",&bVz,"bVz/F");
  h_tree->Branch("brefmult",&brefmult,"brefmult/I");
  h_tree->Branch("btofmult",&btofmult,"btofmult/I");
  h_tree->Branch("refmultcor",&refmultcor,"refmultcor/D");
  h_tree->Branch("reweight",&reweight,"reweight/D");
  h_tree->Branch("cent9",&cent9,"cent9/I");
  h_tree->Branch("bparticleid",&bparticleid,"bparticleid/I");
  h_tree->Branch("bparticlemass",&bparticlemass,"bparticlemass/F");
  h_tree->Branch("bpx",&bpx,"bpx/F");
  h_tree->Branch("bpy",&bpy,"bpy/F");
  h_tree->Branch("bpz",&bpz,"bpz/F");
  h_tree->Branch("hl_ldl",&hl_ldl,"hl_ldl/F");
  h_tree->Branch("hl_dl",&hl_dl,"hl_dl/F");
  h_tree->Branch("hl_l",&hl_l,"hl_l/F");
  h_tree->Branch("hl_chi2topo",&hl_chi2topo,"hl_chi2topo/F");
  h_tree->Branch("hl_chi2ndf",&hl_chi2ndf,"hl_chi2ndf/F");

  h_tree->Branch("bpionpx",&bpionpx,"bpionpx/F");
  h_tree->Branch("bpionpy",&bpionpy,"bpionpy/F");
  h_tree->Branch("bpionpz",&bpionpz,"bpionpz/F");
  h_tree->Branch("bprotonpx",&bprotonpx,"bprotonpx/F");
  h_tree->Branch("bprotonpy",&bprotonpy,"bprotonpy/F");
  h_tree->Branch("bprotonpz",&bprotonpz,"bprotonpz/F");
  h_tree->Branch("bdpx",&bdpx,"bdpx/F");
  h_tree->Branch("bdpy",&bdpy,"bdpy/F");
  h_tree->Branch("bdpz",&bdpz,"bdpz/F");
  h_tree->Branch("bzdeuteron",&bzdeuteron,"bzdeuteron/F");
  h_tree->Branch("bpionnsigma",&bpionnsigma ,"bpionnsigma/F");
  h_tree->Branch("bprotonsigma",&bprotonsigma ,"bprotonsigma/F");
  h_tree->Branch("dca_pion",&dca_pion ,"dca_pion/F");
  h_tree->Branch("dca_proton",&dca_proton ,"dca_proton/F");
  h_tree->Branch("dca_deuteron",&dca_deuteron ,"dca_deuteron/F");
  h_tree->Branch("dca_pion2",&dca_pion2 ,"dca_pion2/F");
  h_tree->Branch("dca_proton2",&dca_proton2 ,"dca_proton2/F");
  h_tree->Branch("dca_deuteron2",&dca_deuteron2 ,"dca_deuteron2/F");

  h_tree->Branch("nhits_pion",&nhits_pion ,"nhits_pion/I");
  h_tree->Branch("nhits_proton",&nhits_proton ,"nhits_proton/I");
  h_tree->Branch("nhits_deuteron",&nhits_deuteron ,"nhits_deuteron/I");
  h_tree->Branch("countrefmult",&countrefmult,"countrefmult/I");
  h_tree->Branch("v_01_pvdca",&v_01_pvdca,"v_01_pvdca/F");
  h_tree->Branch("v_02_pvdca",&v_02_pvdca,"v_02_pvdca/F");
  h_tree->Branch("v_12_pvdca",&v_12_pvdca,"v_12_pvdca/F");
  h_tree->Branch("ht_bdfvtx",&ht_bdfvtx,"ht_bdfvtx/F");
  h_tree->Branch("ht_bdfvtx2",&ht_bdfvtx2,"ht_bdfvtx2/F");
  h_tree->Branch("bpionm2",&bpionm2,"bpionm2/F");
  h_tree->Branch("bprotonm2",&bprotonm2,"bprotonm2/F");
  h_tree->Branch("bdm2",&bdm2,"bdm2/F");

  hvtx      = new TH1F("hvtx",    "Vz;Vz(cm);Counts",100,195,205);
  hvtxgood  = new TH1F("hvtxgood","Vz;Vz(cm);Counts",100,195,205);
  hrefmult  = new TH1F("hrefmult", "refmult; hrefmult; N_{evt}", 600,0,600);
  wrefmult  = new TH1F("wrefmult", "refmult; wrefmult; N_{evt}", 600,0,600);
  hCent= new TH1F("hCent", "Centrality; Centrality; N_{evt}", 9,-0.5,8.5);
  hCentWt= new TH1F("hCentWt", "Centrality; Centrality; N_{evt} (Weight)", 9,-0.5,8.5);
  hvtx_xy   = new TH2F("hvtx_xy",  "Vx;Vx(cm);Vy;Vy(cm)",250,-5,5,250,-5,5);

  cout <<""<< endl;

  return kStOK;
}
//________________________________________________________________________________
Int_t StKFParticleAnalysisMaker::InitRun(Int_t runumber) 
{
  //   assert(StPicoDstMaker::instance());
  //   if (StPicoDstMaker::instance()->IOMode() == StPicoDstMaker::ioRead) {
  //TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO Ask Yuri
  //     StPicoDstMaker::instance()->SetStatus("*",0);
  //     const Char_t *ActiveBranches[] = {
  //       "MuEvent"
  //       ,"PrimaryVertices"
  //       ,"PrimaryTracks"
  //       ,"GlobalTracks"
  //       ,"StStMuMcVertex"
  //       ,"StStMuMcTrack"
  //       ,"CovPrimTrack"
  //       ,"CovGlobTrack"
  //       ,"StStMuMcVertex"
  //       ,"StStMuMcTrack"
  //       ,"KFTracks"
  //       ,"KFVertices"
  //       ,"StBTofHit"
  //       ,"StBTofHeader"
  //     }; 
  //     Int_t Nb = sizeof(ActiveBranches)/sizeof(Char_t *);
  //     for (Int_t i = 0; i < Nb; i++) StPicoDstMaker::instance()->SetStatus(ActiveBranches[i],1); // Set Active braches
  //   }
  return StMaker::InitRun(runumber);
}
//_____________________________________________________________________________
void StKFParticleAnalysisMaker::PrintMem(const Char_t *opt)
{
  MemInfo_t info;
  gSystem->GetMemInfo(&info);
  cout << opt 
    << "\tMemory : Total = " << info.fMemTotal 
    << "\tUsed = " << info.fMemUsed
    << "\tFree = " << info.fMemFree
    << "\tSwap Total = " << info.fSwapTotal
    << "\tUsed = " << info.fSwapUsed
    << "\tFree = " << info.fSwapFree << endl;
}
//_____________________________________________________________________________
void StKFParticleAnalysisMaker::BookVertexPlots()
{
  TDirectory *dirs[2] = {0};
  dirs[0] = TDirectory::CurrentDirectory(); assert(dirs[0]);
  dirs[0]->cd();
  if (! dirs[0]->GetDirectory("Particles")) {
    dirs[0]->mkdir("Particles");
  }
  dirs[1] = dirs[0]->GetDirectory("Particles"); assert(dirs[1]);
  dirs[1]->cd();
  PrintMem(dirs[1]->GetPath());

  fStKFParticleInterface = new StKFParticleInterface;
  bool storeMCHistograms = false;
  if(!fIsPicoAnalysis && fProcessSignal) storeMCHistograms = true;
  fStKFParticlePerformanceInterface = new StKFParticlePerformanceInterface(fStKFParticleInterface->GetTopoReconstructor(), storeMCHistograms);
  dirs[0]->cd();
  PrintMem(dirs[1]->GetPath());
}
//_____________________________________________________________________________
Int_t StKFParticleAnalysisMaker::Make()
{  
  if(fIsPicoAnalysis)
  {
    fPicoDst = StPicoDst::instance();
    if(!fPicoDst) return kStOK;
    if(!fPicoDst->event()) return 0;
    brunid = fPicoDst->event()->runId();
  }
  else
  {  
    fMuDst = StMuDst::instance();
    if(!fMuDst) return kStOK;
    else { if(StMuDst::instance()->numberOfPrimaryVertices() == 0 ) return kStOK; }
    if (!fMuDst->event())  return 0;
    brunid = fMuDst->event()->runId();
  }

  int trigger=0;
  if(fIsPicoAnalysis){
    trigger = PassTrigger(fPicoDst, fsnn);
  }
  else{
    trigger+=1; //auto pass trigger because we mudst already has trigger seltcion
  }

  int notbadrun = CheckIfBadRun( brunid, fsnn);

  bool isGoodEvent = true;
  if(trigger==0) 
  {
    // cout <<"no trigger" << endl;  
    isGoodEvent = false; return kStOK;
  }
  else if(notbadrun!=0) { cout <<"bad run" << endl; isGoodEvent = false; return kStOK;}
  //if(isGoodEvent){
  //cout<<"brunid:"<<brunid<<endl;}else{
  //cout<<"brunid2:"<<brunid<<endl;
  //cout<<"trigger:"<<trigger<<" "<<notbadrun<<endl;
  //}

  //_fill_lambda_tree = false;
  _fill_lambda_tree = true;
  //_fill_he3_tree = true;
  _fill_he3_tree = false;
  bool bWriteTree = true;

  if(fIsPicoAnalysis){
    beventid = fPicoDst->event()->eventId();

    const TVector3 picoPV = fPicoDst->event()->primaryVertex();
    const TVector3 picoPVError = fPicoDst->event()->primaryVertexError();

    bVx = picoPV.x();
    bVy = picoPV.y();
    bVr = sqrt(picoPV.x()*picoPV.x()+picoPV.y()*picoPV.y());
    bVz = picoPV.z(); 
    bVrerr = sqrt(picoPVError.x()*picoPVError.x()+picoPVError.y()*picoPVError.y());
    bVxerr = picoPVError.x();
    bVyerr = picoPVError.y();
    bVzerr = picoPVError.z();
  }
  else{
    beventid = fMuDst->event()->eventId();

    //cout<<"brunid:"<<brunid<<endl;
    float bestRank=-1000000;
    int bestPV=-1;
    for(unsigned int iPV=0; iPV<fMuDst->numberOfPrimaryVertices(); iPV++)
    {
      StMuPrimaryVertex *Vtx = fMuDst->primaryVertex(iPV);
      if(!Vtx) continue;
      if (bestRank < Vtx->ranking()) {
        bestRank = Vtx->ranking();
        bestPV = iPV;
      }
      else continue;

      bVr = sqrt(Vtx->position().x()*Vtx->position().x()+Vtx->position().y()*Vtx->position().y());
      bVx = Vtx->position().x();
      bVy = Vtx->position().y();
      bVz = Vtx->position().z();
      bVrerr = sqrt(Vtx->posError().x()*Vtx->posError().x()+Vtx->posError().y()*Vtx->posError().y());
      bVxerr = Vtx->posError().x();
      bVyerr = Vtx->posError().y();
      bVzerr = Vtx->posError().z();
    }
  }

  hvtx->Fill(bVz);
  hvtx_xy->Fill(bVx,bVy);

  Vtxbin = getVtxBin(bVx, bVy, bVz);
  if (Vtxbin<0) isGoodEvent = false;

  //cut on refmult and tofmult, trigger, to cut on good event. this cut is only for real analysis!! 
  //when fStoreTmvaNTuples is true, there is no need to cut on these variables
  if(fIsPicoAnalysis)
  {
    brefmult = fPicoDst->event()->refMult();
    btofmult = fPicoDst->event()->btofTrayMultiplicity();
  }
  else
  {
    brefmult = fMuDst->event()->refMult();
    btofmult = fMuDst->event()->btofTrayMultiplicity();
  }

  if(fIsPicoAnalysis)
  {
    countrefmult = 0;
    FXTMult2 = 0;
    for(unsigned int iTrack = 0; iTrack < fPicoDst->numberOfTracks(); iTrack++)
    {
      const StPicoTrack *ptrk = (StPicoTrack*)fPicoDst->track(iTrack);
      if(! ptrk) continue;
      if(!ptrk->isPrimary())  continue;  // now selecting primary tracks
      FXTMult2++;
      const float dca = ptrk->gDCA( fPicoDst->event()->primaryVertex() ).Mag();
      const int nHitsFit = ptrk->nHitsFit();
      const int nHitsPoss = ptrk->nHitsMax();
      const float quality = (float)nHitsFit/(float)nHitsPoss;
      if( fabs(dca)>3.0 ) continue;
      if( nHitsFit < 15 )  continue;;
      if( quality < 0.52 )  continue;
      countrefmult++;
    }
  }
  else 
  {
    countrefmult = 0;
    FXTMult2 = 0;
    float bestRank=-1000000;
    int bestPV=-1;
    for(unsigned int iPV=0; iPV<fMuDst->numberOfPrimaryVertices(); iPV++)
    {
      StMuPrimaryVertex *Vtx = fMuDst->primaryVertex(iPV);
      if(!Vtx) continue;
      if (bestRank < Vtx->ranking()) {
        bestRank = Vtx->ranking();
        bestPV = iPV;
      }
      else continue;
    }
    if(bestPV!=-1){
      for(unsigned int iTrack = 0; iTrack < fMuDst->numberOfPrimaryTracks(); iTrack++)
      {
        StMuTrack *gTrack = fMuDst->primaryTracks(iTrack);
        if (! gTrack) continue;
        FXTMult2++;
        const int bnHitsFit = gTrack->nHitsFit();
        const int bnHitsPoss = gTrack->nHitsPoss();
        const float bquality = (float)bnHitsFit/(float)bnHitsPoss;
        const float bdca = gTrack->dcaGlobal(bestPV).mag();
        if ( bnHitsFit < 15 ) continue;
        if ( bquality < 0.52 ) continue;
        if ( fabs(bdca) > 3.0 ) continue;
        countrefmult++;
      }
    }
  }

  reweight = 1;
  refmultcor = 0;
  int centralityBin = -1;

  if(fRunCentralityAnalysis)
  {
    if(fIsPicoAnalysis){
      //somehow this tool not work with TFG18n 
      //read centrality from flow tree instead
      // mPileupTool->initEvent(fPicoDst);
      // countrefmult = mPileupTool->get_refMultPrim();
      // centralityBin  = mPileupTool->get_centrality9();  
      // cent9 = centralityBin;  
      // reweight = mPileupTool->get_centralityWeight();
      // if (centralityBin<0) isGoodEvent= false;
    }
    else{
      //no for 3GeV embedding
      //will add manually later

    }
  }

  if(fFlowAnalysis)
  {
    int eventId = -1;
    int runId = -1;
    EPbin=-1;

    if(fIsPicoAnalysis) 
    {
      runId   = fPicoDst->event()->runId();
      eventId = fPicoDst->event()->eventId();
    }
    else
    {
      runId   = fMuDst->event()->runId();
      eventId = fMuDst->event()->eventId();
    }

    long entryId = GetUniqueEventId(runId, eventId);
    //cout<<"entryId:"<<entryId<<" "<<runId<<" "<<eventId<<endl;
    std::map<long,int>::iterator flowMapIterator = fFlowMap.find(entryId);
    if (flowMapIterator != fFlowMap.end())
    {
      fFlowChain->GetEvent(fFlowMap[GetUniqueEventId(runId, eventId)]);       
    }
    centralityBin = fCentrality;
    cent9 = centralityBin;
    reweight = gweight;
    EPbin = getEventPlaneBin(psi_1_EPD_4);
    //cout<<"fFlowRunId:"<<fFlowRunId<<" "<< fFlowEventId<< " "<< fCentrality<<" "<< endl;
    //cout<<"flow:"<<psi_1_EPD_0 <<" "<< psi_1_EPD_1<<" "<<psi_1_EPD_2<<" "<<psi_1_EPD_3<<endl;
    if (EPbin<0) isGoodEvent = false; 
  }

  if (cent9<0) isGoodEvent=false;


  if(isGoodEvent){
    hvtxgood->Fill(bVz);
    // hrefmult->Fill(FXTMult2);
    hrefmult->Fill(countrefmult);
    // cout << centralityBin<<" "<<reweight<< endl;
    hCentWt->Fill(centralityBin, reweight);
    hCent->Fill(centralityBin);
    // wrefmult->Fill(refmultcor);
    wrefmult->Fill(countrefmult, reweight);
  }

  //find max global track index
  int maxGBTrackIndex = -1;
  if(fIsPicoAnalysis)
  {
    for(unsigned int iTrack = 0; iTrack < fPicoDst->numberOfTracks(); iTrack++) 
    {
      StPicoTrack *gTrack = fPicoDst->track(iTrack);
      if (! gTrack) continue;
      int index = gTrack->id();
      if(index > maxGBTrackIndex)
        maxGBTrackIndex = index;
    }
  }
  else
  {
    for(unsigned int iTrack = 0; iTrack < fMuDst->numberOfGlobalTracks(); iTrack++) 
    {
      StMuTrack *gTrack = fMuDst->globalTracks(iTrack);
      if (! gTrack) continue;
      int index = gTrack->id();
      if(index > maxGBTrackIndex)
        maxGBTrackIndex = index;
    }
  }
  vector<KFMCTrack> mcTracks(0);
  vector<int> mcIndices(maxGBTrackIndex+1);
  for(unsigned int iIndex=0; iIndex<mcIndices.size(); iIndex++)
    mcIndices[iIndex] = -1;

  //Process the event
  vector<int> triggeredTracks;
  if(maxGBTrackIndex > 0)
    fStKFParticleInterface->ResizeTrackPidVectors(maxGBTrackIndex+1);
  if(fIsPicoAnalysis)
    isGoodEvent = isGoodEvent && fStKFParticleInterface->ProcessEvent(fPicoDst, triggeredTracks);
  else//b
    //no embedding
    //    isGoodEvent = fStKFParticleInterface->ProcessEvent(fMuDst, triggeredTracks, fProcessSignal);
    //embedding
    isGoodEvent = isGoodEvent && fStKFParticleInterface->ProcessEvent(fMuDst, mcTracks, mcIndices, fProcessSignal);

  // cout <<"test: after process event: " <<isGoodEvent<< endl;
  //collect histograms

  //events removed at this level include no track events and 
  //whatever is in Process event skip level
  //  if (isGoodEvent)  cout <<"is good event" << endl;
  //  else cout <<"not good event" << endl;

  if(isGoodEvent)
  {

    if(fTMVAselection){
      cout<<"TMVA SELECTION!!!!!!!!!!!!! SOMETHING IS WRONG-----------------"<<endl;
    }

    if(fTMVAselection)//check this bool
    {
      for(int iParticle=0; iParticle<fStKFParticlePerformanceInterface->GetNReconstructedParticles(); iParticle++)
      {
        KFParticle particle = fStKFParticleInterface->GetParticles()[iParticle];

        for(int iReader=0; iReader<fNNTuples; iReader++)
        {
          if( abs(particle.GetPDG()) == fNTuplePDG[iReader] )
          {
            GetParticleParameters(iReader, particle);

            const int iTMVACentralityBin = GetTMVACentralityBin(iReader, centralityBin);
            const int iTMVAPtBin = GetTMVAPtBin(iReader, particle.GetPt());

            if(iTMVACentralityBin<0 || iTMVAPtBin<0) 
            {
              fStKFParticleInterface->RemoveParticle(iParticle);
              continue;
            }

            if(fTMVAReader[iReader][iTMVACentralityBin][iTMVAPtBin]->EvaluateMVA("BDT") < fTMVACut[iReader][iTMVACentralityBin][iTMVAPtBin])
              fStKFParticleInterface->RemoveParticle(iParticle);
            /* 
               if(fAnalyseDsPhiPi && abs(fStKFParticleInterface->GetParticles()[iParticle].GetPDG()) == 431)
               {              
               KFParticle phi;
               if(particle.GetPDG() == 431)
               phi += fStKFParticleInterface->GetParticles()[particle.DaughterIds()[0]];
               else
               phi += fStKFParticleInterface->GetParticles()[particle.DaughterIds()[1]];
               phi += fStKFParticleInterface->GetParticles()[particle.DaughterIds()[2]];
               float mass = 0.f, dmass = 0.f;
               phi.GetMass(mass, dmass);
               if( fabs(mass - 1.01946) > 0.015)
               fStKFParticleInterface->RemoveParticle(iParticle);
               }
               */
          }
        }
      }      
    }

    //clean H3L, H4L, Ln, Lnn
    for(int iParticle=0; iParticle<fStKFParticlePerformanceInterface->GetNReconstructedParticles(); iParticle++)
    {
      KFParticle particle = fStKFParticleInterface->GetParticles()[iParticle];
      if( abs(particle.GetPDG())==3003 || abs(particle.GetPDG())==3103 || abs(particle.GetPDG())==3004 || abs(particle.GetPDG())==3005)
      {
        //         if(particle.GetP() < 1.)
        //         {
        //           fStKFParticleInterface->RemoveParticle(iParticle);
        //           continue;
        //         }

        //         if(particle.GetPhi() > -0.8 && particle.GetPhi() < -0.4)
        //         {
        //           fStKFParticleInterface->RemoveParticle(iParticle);
        //           continue;
        //         }

        for(int iD=0; iD<particle.NDaughters(); iD++)
        {
          const int daughterId = particle.DaughterIds()[iD];
          const KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterId];
          //  if (abs(daughter.GetPDG())==211 && daughter.GetP() > 0.5)
          //      fStKFParticleInterface->RemoveParticle(iParticle);
        }
      }
    }


    fStKFParticlePerformanceInterface->SetMCTracks(mcTracks);
    fStKFParticlePerformanceInterface->SetMCIndexes(mcIndices);    
    fStKFParticlePerformanceInterface->SetCentralityBin(centralityBin);
    fStKFParticlePerformanceInterface->SetCentralityWeight(reweight);
    Int_t nevent = 100000;
    fStKFParticlePerformanceInterface->SetPrintEffFrequency(nevent);

    //   cout<<"eventIdcheck:"<<fPicoDst->event()->eventId()<<endl;
    fStKFParticlePerformanceInterface->PerformanceAnalysis();
    //    cout<<"fStKFParticlePerformanceInterface->GetNReconstructedParticles()"<<
    //    fStKFParticlePerformanceInterface->GetNReconstructedParticles()<<endl;


    if(bWriteTree){
      //add mixed event for H3L 103004(in KF particle, first recosntruct ppi, then ppi+d), mix ppi and d 
      // add one more branch, bisMix as mixed event tag 
      // if (fMixEvent && fIsPicoAnalysis && isGoodEvent && EPbin>=0 && centralityBin>=0 && centralityBin<9)
      if (fMixEvent && fIsPicoAnalysis && isGoodEvent && EPbin>=0 && centralityBin>=0 && centralityBin<9 && Vtxbin>=0 && Vtxbin<125)
      {
        StMixEvent newEvent;
        //will make it to be function later
        // makeNewMixEvent(newEvent, dauPDG1, dauPDG2, massrange);
        // cout <<"convert cuttent event" <<endl;
        KFPTrackVector* rTracks=fStKFParticleInterface->GetTopoReconstructor()->GetTracks();
        for (int iTrack= rTracks[0].FirstDeuteron(); iTrack<rTracks[0].LastDeuteron() ;iTrack++)
          // for (int iTrack= 0; iTrack<rTracks[0].Size() ;iTrack++)
        {
          int pdgId = rTracks[0].PDG()[iTrack];

          if (!(abs(pdgId)==1000010020)) continue;
          // cout << rTracks[0].PVIndex()[iTrack]<<endl;
          // if (rTracks[0].PVIndex()[iTrack]>=0) continue; //not secondary track
          // const int_v& trackPVIndex = reinterpret_cast<const  int_v&>(rTracks[0].PVIndex()[iTrack]);
          // const int_m& isTrackSecondary = (trackPVIndex < 0);
          // if (isTrackSecondary.isEmpty()) contiue; 

          // cout << rTracks[0].Id()[iTrack]<<endl;
          KFParticle particle; 
          fStKFParticlePerformanceInterface->GetParticle(particle, rTracks[0].Id()[iTrack]);
          // cout <<particle.DaughterIds()[0] << endl;
          float zdeuteron, chi2prim, dca,m2;
          int nhits;
          FillDaughterInfo(particle, particle.DaughterIds()[0], 1000010020, chi2prim, nhits, dca, zdeuteron, m2 );
          if (nhits<15) continue; 
          newEvent.d_v_chi2prim.push_back(chi2prim);
          newEvent.d_v_dca.push_back(dca);
          newEvent.d_v_nhits.push_back(nhits);
          newEvent.d_v_m2.push_back(m2);
          newEvent.d_v_z.push_back(zdeuteron);
          newEvent.d_v_dEdx.push_back(fStKFParticleInterface->GetdEdX(particle.DaughterIds()[0]));           
          newEvent.d_v.push_back(particle);   
        }
        for (int iParticle=0; iParticle<fStKFParticlePerformanceInterface->GetNReconstructedParticles(); iParticle++)
        { 
          KFParticle particle;
          fStKFParticlePerformanceInterface->GetParticle(particle, iParticle);
          int pdgId = particle.GetPDG();
          double mass = particle.GetMass();
          if (abs(pdgId)==1000010020)
          {
            // newEvent.d_v.push_back(particle);   
            // float zdeuteron, chi2prim, dca,m2;
            // int nhits;
            // FillDaughterInfo(particle, particle.DaughterIds()[0], 1000010020, chi2prim, nhits, dca, zdeuteron, m2 );
            // newEvent.d_v_chi2prim.push_back(chi2prim);
            // newEvent.d_v_dca.push_back(dca);
            // newEvent.d_v_nhits.push_back(nhits);
            // newEvent.d_v_m2.push_back(m2);
            // newEvent.d_v_z.push_back(zdeuteron);
            // newEvent.d_v.push_back(particle);
            // newEvent.d_v_dEdx.push_back(fStKFParticleInterface->GetdEdX(particle.DaughterIds()[0]));

          }
          else if (abs(pdgId)==3012 && mass <1.15  && mass >1.05 )
          {
            if (particle.NDaughters()!=2) {
              cout << " this lambda daughters: "<<particle.NDaughters()<< endl;
              continue; 
            }
            particle.TransportToDecayVertex();
            KFParticle pion = fStKFParticleInterface->GetParticles()[particle.DaughterIds()[0]];
            KFParticle proton = fStKFParticleInterface->GetParticles()[particle.DaughterIds()[1]];
            float sigma, chi2prim, dca, m2;
            int nhits;
            FillDaughterInfo(pion, pion.DaughterIds()[0], 211, chi2prim, nhits, dca, sigma, m2 );
            newEvent.pi_v_chi2prim.push_back(chi2prim);
            newEvent.pi_v_dca.push_back(dca);
            newEvent.pi_v_nhits.push_back(nhits);
            newEvent.pi_v_m2.push_back(m2);
            newEvent.pi_v_sigma.push_back(sigma);
            newEvent.pi_v.push_back(pion);

            FillDaughterInfo(proton, proton.DaughterIds()[0], 2212, chi2prim, nhits, dca, sigma, m2 );
            newEvent.p_v_chi2prim.push_back(chi2prim);
            newEvent.p_v_dca.push_back(dca);
            newEvent.p_v_nhits.push_back(nhits);
            newEvent.p_v_m2.push_back(m2);
            newEvent.p_v_sigma.push_back(sigma);
            newEvent.p_v.push_back(proton);

            newEvent.ppi_v.push_back(particle);
          }
        }
        newEvent.Vx = bVx;
        newEvent.Vy = bVy;
        newEvent.Vz = bVz;

        bool fillBuffer = newEvent.p_v.size()>0 && newEvent.d_v.size()>0;
        if (fillBuffer) {
          // cout << "push current event into buffer"<<endl;
          //fot test
          bisMix=2;
          cout <<newEvent.p_v.size()<<" "<<newEvent.p_v_dca.size() << endl;
            for (int id = 0;id<newEvent.d_v.size();id++){
              for (int iL = 0;iL<newEvent.p_v.size();iL++){
                // cout << " start fill mix tree: d+mix Lambda "<< ie << " lambda idx "<< iL<< endl;
                bool filltree = FillThreeDaughtersMix(newEvent.pi_v[iL],newEvent.p_v[iL],newEvent.d_v[id], newEvent.ppi_v[iL], 0,0,0,0);
                //add daughter further info
                chi2primary_d = newEvent.d_v_chi2prim[id];
                nhits_deuteron = newEvent.d_v_nhits[id];
                dca_deuteron = newEvent.d_v_dca[id];
                bdm2 = newEvent.d_v_m2[id];
                dca_deuteron = newEvent.d_v_dca[id];
                bdedx = newEvent.d_v_dEdx[id];
                bzdeuteron =  newEvent.d_v_z[id];
                 
                chi2primary_pi = newEvent.pi_v_chi2prim[iL];
                bpionnsigma = newEvent.pi_v_sigma[iL];
                nhits_pion = newEvent.pi_v_nhits[iL];
                dca_pion = newEvent.pi_v_dca[iL];
                bpionm2 = newEvent.pi_v_m2[iL];

                chi2primary_proton = newEvent.p_v_chi2prim[iL];
                bprotonsigma = newEvent.p_v_sigma[iL];
                nhits_proton = newEvent.p_v_nhits[iL];
                dca_proton = newEvent.p_v_dca[iL];
                bprotonm2 = newEvent.p_v_m2[iL];    
                if (filltree) 
                {
                  htriton3_tree->Fill();
                } 
              }                
            }
        }
        
        if ( bMixEventBuffer[centralityBin][EPbin][Vtxbin].size()==fMixEventBufferSize && fillBuffer){
          // cout <<"buffer is full, start mix, EPbin: " << EPbin<<" EP angle: "<<psi_1_EPD_4<<" cent: " <<centralityBin<< endl;
          // cout <<"current event d size: " <<  newEvent.d_v.size() << " (ppi): "<< newEvent.p_v.size()<< endl;
          // attention kfparticle pre-cuts
          bisMix = 1;
          bismc = 0;
          for (int ie =0;ie<fMixEventBufferSize;ie++){ 
            // cout <<"mix event d size: " <<  bMixEventBuffer[centralityBin][EPbin][ie].d_v.size() << " Lambda: "<< bMixEventBuffer[centralityBin][EPbin][ie].p_v.size()<< endl;
            float dVx = newEvent.Vx - bMixEventBuffer[centralityBin][EPbin][Vtxbin][ie].Vx;
            float dVy = newEvent.Vy - bMixEventBuffer[centralityBin][EPbin][Vtxbin][ie].Vy;
            float dVz = newEvent.Vz - bMixEventBuffer[centralityBin][EPbin][Vtxbin][ie].Vz;
            // cout <<dVx<<" "<<dVy<<" "<<dVz << endl;

            for (int id = 0;id<bMixEventBuffer[centralityBin][EPbin][Vtxbin][ie].d_v.size();id++){
              for (int iL = 0;iL<newEvent.p_v.size();iL++){
                //will update here!
                bool filltree = FillThreeDaughtersMix(newEvent.pi_v[iL],newEvent.p_v[iL],bMixEventBuffer[centralityBin][EPbin][Vtxbin][ie].d_v[id], newEvent.ppi_v[iL], dVx, dVy, dVz, 1);
                if (filltree) 
                {
                  //add daughter further info
                  nhits_deuteron = bMixEventBuffer[centralityBin][EPbin][Vtxbin][ie].d_v_nhits[id];
                  dca_deuteron = bMixEventBuffer[centralityBin][EPbin][Vtxbin][ie].d_v_dca[id];
                  bdm2 = bMixEventBuffer[centralityBin][EPbin][Vtxbin][ie].d_v_m2[id];
                  chi2primary_d = bMixEventBuffer[centralityBin][EPbin][Vtxbin][ie].d_v_chi2prim[id];
                  bdedx = bMixEventBuffer[centralityBin][EPbin][Vtxbin][ie].d_v_dEdx[id];
                  bzdeuteron =  bMixEventBuffer[centralityBin][EPbin][Vtxbin][ie].d_v_z[id];
                  //     
                  chi2primary_pi = newEvent.pi_v_chi2prim[iL];
                  bpionnsigma = newEvent.pi_v_sigma[iL];
                  nhits_pion = newEvent.pi_v_nhits[iL];
                  dca_pion = newEvent.pi_v_dca[iL];
                  bpionm2 = newEvent.pi_v_m2[iL];

                  chi2primary_proton = newEvent.p_v_chi2prim[iL];
                  bprotonsigma = newEvent.p_v_sigma[iL];
                  nhits_proton = newEvent.p_v_nhits[iL];
                  dca_proton = newEvent.p_v_dca[iL];
                  bprotonm2 = newEvent.p_v_m2[iL];

                  htriton3_tree->Fill(); 
                }
              }                
            }
            for (int id = 0;id<newEvent.d_v.size();id++){
              for (int iL = 0;iL<bMixEventBuffer[centralityBin][EPbin][Vtxbin][ie].p_v.size();iL++){
                bool filltree = FillThreeDaughtersMix(bMixEventBuffer[centralityBin][EPbin][Vtxbin][ie].pi_v[iL],bMixEventBuffer[centralityBin][EPbin][Vtxbin][ie].p_v[iL],newEvent.d_v[id], bMixEventBuffer[centralityBin][EPbin][Vtxbin][ie].ppi_v[iL], dVx, dVy, dVz, 2);
                  //add daughter further info
                  chi2primary_d = newEvent.d_v_chi2prim[id];
                  nhits_deuteron = newEvent.d_v_nhits[id];
                  dca_deuteron = newEvent.d_v_dca[id];
                  bdm2 = newEvent.d_v_m2[id];
                  dca_deuteron = newEvent.d_v_dca[id];
                  bdedx = newEvent.d_v_dEdx[id];
                  bzdeuteron =  newEvent.d_v_z[id];
                  //     
                  chi2primary_pi = bMixEventBuffer[centralityBin][EPbin][Vtxbin][ie].pi_v_chi2prim[iL];
                  bpionnsigma = bMixEventBuffer[centralityBin][EPbin][Vtxbin][ie].pi_v_sigma[iL];
                  nhits_pion = bMixEventBuffer[centralityBin][EPbin][Vtxbin][ie].pi_v_nhits[iL];
                  dca_pion = bMixEventBuffer[centralityBin][EPbin][Vtxbin][ie].pi_v_dca[iL];
                  bpionm2 = bMixEventBuffer[centralityBin][EPbin][Vtxbin][ie].pi_v_m2[iL];

                  chi2primary_proton = bMixEventBuffer[centralityBin][EPbin][Vtxbin][ie].p_v_chi2prim[iL];
                  bprotonsigma = bMixEventBuffer[centralityBin][EPbin][Vtxbin][ie].p_v_sigma[iL];
                  nhits_proton = bMixEventBuffer[centralityBin][EPbin][Vtxbin][ie].p_v_nhits[iL];
                  dca_proton = bMixEventBuffer[centralityBin][EPbin][Vtxbin][ie].p_v_dca[iL];
                  bprotonm2 = bMixEventBuffer[centralityBin][EPbin][Vtxbin][ie].p_v_m2[iL];
                if (filltree) 
                {
                  // cout << " finish fill mix tree: d+mix Lambda "<< ie << " lambda idx "<< iL<< endl;
                  htriton3_tree->Fill();
                } 
              }                
            }
          }
          // cout << "erase the oldest event"<<endl;
          bMixEventBuffer[centralityBin][EPbin][Vtxbin].pop_front();
        }
        if (fillBuffer) {
          // cout << "push current event into buffer"<<endl;
          bMixEventBuffer[centralityBin][EPbin][Vtxbin].push_back(newEvent);
        }
      }

      //check track ID:
      // cout << fStKFParticlePerformanceInterface->GetNReconstructedParticles()<< endl;
      for(int iParticle=0; iParticle<fStKFParticlePerformanceInterface->GetNReconstructedParticles(); iParticle++)
      {
        bisMix = false;
        KFParticle particle;
        //       const KFParticle particle = fStKFParticleInterface->GetParticles()[iParticle];
        bool isMCParticle = fStKFParticlePerformanceInterface->GetParticle(particle, iParticle) && !fIsPicoAnalysis;

        //int beventid;
        //int brunid;
        //int bparticleid;
        //float bparticlemass;
        bismc =0;
        if(isMCParticle) bismc=1;

        bparticleid   = particle.GetPDG();
        bparticlemass = particle.GetMass();

        bx = particle.GetX();
        by = particle.GetY();
        bz = particle.GetZ();
        bpx = particle.GetPx();
        bpy = particle.GetPy();
        bpz = particle.GetPz();
        bdl = particle.GetDecayLength();
        //if(particle.IdTruth()!=0){ 
        //cout<<"IdTruth():"<<particle.IdTruth()<<" "<< particle.IdParentMcVx()  <<endl;
        //         }

        //Lambda
        if(_fill_lambda_tree && isGoodEvent && bparticlemass<1.2){
          if(fabs(particle.GetPDG())==3122 )
            // if(fabs(particle.GetPDG())==3122  || fabs(particle.GetPDG())==3012 )
          {
            KFParticleSIMD tempSIMDParticle(particle);
            float_v l,dl;
            KFParticleSIMD pv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
            tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl);
            ld_ldl = l[0]/dl[0];
            ld_l = l[0];
            ld_dl = dl[0];

            bpl = ld_l/(sqrt(bpx*bpx+bpy*bpy+bpz*bpz)/bparticlemass);

            // cout<<"GetDistanceFromVertex:"<<tempSIMDParticle.GetDistanceFromVertex(pv)<<endl;
            // cout<<"GetLifeTime:"<<tempSIMDParticle.GetLifeTime()<<endl;
            ld_bdfvtx = tempSIMDParticle.GetDistanceFromVertex(pv)[0];
            ld_bdfvtx_xy = tempSIMDParticle.GetDistanceFromVertexXY(pv)[0];
            ld_bdfvtxdev_xy = tempSIMDParticle.GetDeviationFromVertexXY(pv)[0];

            // cout<<"ld_chi2topo:"<<tempSIMDParticle.Chi2()[0]<<"/"<<tempSIMDParticle.NDF()[0]<<endl;
            tempSIMDParticle.SetProductionVertex(pv);
            //cout<<"pv:"<<pv.GetX()<<" "<<pv.GetErrX()<<" "<<pv.GetY()<<" "<<pv.GetErrY()<<pv.GetZ()<<" "<<pv.GetErrZ()<<   endl;
            //cout<<"ld_chi2topo:"<<tempSIMDParticle.GetChi2()[0]<<"/"<<tempSIMDParticle.NDF()[0]<<endl;
            // cout<<"ld_dca:"<<tempSIMDParticle.GetDistanceFromVertex(pv)<<" "<<tempSIMDParticle.GetDeviationFromVertex(pv)<<endl;


            //cout<<"GetLifeTime:"<<tempSIMDParticle.GetLifeTime()<<endl;
            //cout<<"GetDistanceFromVertex2:"<<tempSIMDParticle.GetDistanceFromVertex(pv)<<endl;
            ld_bdfvtx2 = tempSIMDParticle.GetDistanceFromVertex(pv)[0];
            ld_lifetime = tempSIMDParticle.GetLifeTime()[0];


            ld_chi2topo = double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);
            ld_chi2ndf = particle.Chi2()/particle.NDF();
            ld_chi2primary =particle.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex()); 

            for(int iDaughter=0; iDaughter<particle.NDaughters(); iDaughter++)
            {
              int order[4] = {0, 1, 2, 3};
              const int daughterParticleIndex = particle.DaughterIds()[order[iDaughter]];
              KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterParticleIndex];
              if(iDaughter==0){
                //	cout<<"pdg:"<<daughter.GetPDG()<<" "<<particle.GetPDG()<<" "<<iDaughter<<endl;
                chi2primary_pi = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
                nhits_ld_pi = fStKFParticleInterface->GetNHits(daughter.DaughterIds()[0]);
                dca_pi = fStKFParticleInterface->Getdca(daughter.DaughterIds()[0]);
              }
              if(iDaughter==1){
                //						cout<<"pdg:"<<daughter.GetPDG()<<" "<<particle.GetPDG()<<" "<<iDaughter<<endl;
                chi2primary_proton = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
                nhits_ld_proton = fStKFParticleInterface->GetNHits(daughter.DaughterIds()[0]);
                dca_proton = fStKFParticleInterface->Getdca(daughter.DaughterIds()[0]);
              }

            }

            bmcpx=-999;
            bmcpy=-999;
            bmcpz=-999;
            bmcpl=-999;

            //if(isMCParticle && fProcessSignal)
            if(isMCParticle){
              int iMCPart = fStKFParticlePerformanceInterface->GetParticleMCId(iParticle);
              StMuMcTrack *mcTrack = fMuDst->MCtrack(iMCPart);
              bmcpx = mcTrack->Pxyz().x();
              bmcpy = mcTrack->Pxyz().y();
              bmcpz = mcTrack->Pxyz().z();
              bmcidvx = mcTrack->IdVx();
              bmcidvxend = mcTrack->IdVxEnd();
              StMuMcVertex *mcvx = fMuDst->MCvertex(bmcidvx-1);
              StMuMcVertex *mcvxend = fMuDst->MCvertex(bmcidvxend-1);

              //cout<<"particle:"<<particle.X()<<" "<<particle.Y()<<" "<<particle.Z()<<" "<<endl;
              //cout<<"mcpvend:"<<mcvxend->XyzV().x()<<" "<<mcvxend->XyzV().y()<<" "<<mcvxend->XyzV().z()<<endl;
              //cout<<endl;
              //cout<<"particle p:"<<particle.Px()<<" "<<particle.Py()<<" "<<particle.Pz()<<" "<<endl;
              //cout<<"mcpvend:"<<bmcpx<<" "<<bmcpy<<" "<<bmcpz<<endl;
              //cout<<endl;
              //cout<<"pv:"<<pv.X()<<" "<<pv.Y()<<" "<<pv.Z()<<endl;
              //cout<<"mcpv:"<<mcvx->XyzV().x()<<" "<<mcvx->XyzV().y()<<" "<<mcvx->XyzV().z()<<endl;

              bmcx = mcvxend->XyzV().x();
              bmcy = mcvxend->XyzV().y();
              bmcz = mcvxend->XyzV().z(); 
              bx = particle.X();
              by = particle.Y();
              bz = particle.Z();

              bmcl = sqrt( (mcvx->XyzV().x()-mcvxend->XyzV().x())*(mcvx->XyzV().x()-mcvxend->XyzV().x()) + (mcvx->XyzV().y()-mcvxend->XyzV().y())*(mcvx->XyzV().y()-mcvxend->XyzV().y()) + (mcvx->XyzV().z()-mcvxend->XyzV().z())*(mcvx->XyzV().z()-mcvxend->XyzV().z()) ) ;
              bmcpl = bmcl * 1.115683 / sqrt(bmcpx*bmcpx+bmcpy*bmcpy+bmcpz*bmcpz) ;
            }
            lambda_tree->Fill();
          }
        }

        //Kshort
        if(isGoodEvent && particle.GetPDG()==310){
          KFParticleSIMD tempSIMDParticle(particle);
          float_v l,dl;
          KFParticleSIMD pv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
          tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl);
          ld_ldl = l[0]/dl[0];
          ld_l = l[0];
          ld_dl = dl[0];
          tempSIMDParticle.SetProductionVertex(pv);
          ld_chi2topo = double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);
          ld_chi2ndf = particle.Chi2()/particle.NDF();
          for(int iDaughter=0; iDaughter<particle.NDaughters(); iDaughter++)
          {
            int order[4] = {0, 1, 2, 3};
            const int daughterParticleIndex = particle.DaughterIds()[order[iDaughter]];
            KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterParticleIndex];  
            if(iDaughter==0){
              chi2primary_proton = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
              nhits_ld_proton = fStKFParticleInterface->GetNHits(daughter.DaughterIds()[0]);
              dca_proton = fStKFParticleInterface->Getdca(daughter.DaughterIds()[0]);
            }
            if(iDaughter==1){
              chi2primary_pi = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
              nhits_ld_pi = fStKFParticleInterface->GetNHits(daughter.DaughterIds()[0]);
              dca_pi = fStKFParticleInterface->Getdca(daughter.DaughterIds()[0]);
            }
          }
          bmcpx=-999;
          bmcpy=-999;
          bmcpz=-999;
          if(isMCParticle && fProcessSignal){
            int iMCPart = fStKFParticlePerformanceInterface->GetParticleMCId(iParticle);
            StMuMcTrack *mcTrack = fMuDst->MCtrack(iMCPart);
            bmcpx = mcTrack->Pxyz().x();
            bmcpy = mcTrack->Pxyz().y();
            bmcpz = mcTrack->Pxyz().z();

            bmcidvx = mcTrack->IdVx();
            bmcidvxend = mcTrack->IdVxEnd();
            StMuMcVertex *mcvx = fMuDst->MCvertex(bmcidvx-1);
            StMuMcVertex *mcvxend = fMuDst->MCvertex(bmcidvxend-1);

            bmcx = mcvxend->XyzV().x();
            bmcy = mcvxend->XyzV().y();
            bmcz = mcvxend->XyzV().z();
            bx = particle.X();
            by = particle.Y();
            bz = particle.Z();
          }
          ks_tree->Fill();
        }

        //Xi
        if(fabs(particle.GetPDG())==3312 && isGoodEvent){

          KFParticleSIMD tempSIMDParticle(particle);
          float_v l,dl;
          KFParticleSIMD pv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
          tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl);
          xi_ldl = l[0]/dl[0];
          xi_l = l[0];
          xi_dl = dl[0];

          tempSIMDParticle.SetProductionVertex(pv);
          xi_chi2topo = double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);
          xi_chi2ndf = particle.Chi2()/particle.NDF();

          for(int iDaughter=0; iDaughter<particle.NDaughters(); iDaughter++)
          {
            int order[4] = {0, 1, 2, 3};
            const int daughterParticleIndex = particle.DaughterIds()[order[iDaughter]];
            KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterParticleIndex];

            if(iDaughter==0){
              chi2primary_xi_bach = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
            }
            if(iDaughter==1){
              chi2primary_xi_ld = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
              KFParticleSIMD ttempSIMDParticle(daughter);
              float_v tl,tdl;

              KFParticleSIMD tpv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
              ttempSIMDParticle.GetDistanceToVertexLine(tpv, tl, tdl);
              xi_ld_ldl = tl[0]/tdl[0];
              xi_ld_l = tl[0];
              ttempSIMDParticle.SetProductionVertex(tpv);
              xi_ld_chi2topo = double(ttempSIMDParticle.Chi2()[0])/double(ttempSIMDParticle.NDF()[0]);
              xi_ld_chi2ndf = daughter.Chi2()/daughter.NDF();

              for(int jDaughter=0; jDaughter<daughter.NDaughters(); jDaughter++){
                int jorder[4] = {0, 1, 2, 3};
                const int jdaughterParticleIndex = daughter.DaughterIds()[jorder[jDaughter]];
                KFParticle granddaughter = fStKFParticleInterface->GetParticles()[jdaughterParticleIndex];
                if(jDaughter==0){
                  chi2primary_xi_proton = granddaughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
                }
                if(jDaughter==1){
                  chi2primary_xi_pi = granddaughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
                }
              }
            }
          }
          bmcpx=-999;
          bmcpy=-999;
          bmcpz=-999;
          if(isMCParticle && fProcessSignal){
            int iMCPart = fStKFParticlePerformanceInterface->GetParticleMCId(iParticle);
            StMuMcTrack *mcTrack = fMuDst->MCtrack(iMCPart);
            bmcpx = mcTrack->Pxyz().x();
            bmcpy = mcTrack->Pxyz().y();
            bmcpz = mcTrack->Pxyz().z();
          }
          cascade_tree->Fill();
        }


        if(fabs(particle.GetPDG())==3334 && isGoodEvent){

          KFParticleSIMD tempSIMDParticle(particle);
          float_v l,dl;
          KFParticleSIMD pv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
          tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl);
          om_ldl = l[0]/dl[0];
          om_l = l[0];
          om_dl = dl[0];


          tempSIMDParticle.SetProductionVertex(pv);
          om_chi2topo = double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);
          om_chi2ndf = particle.Chi2()/particle.NDF();

          for(int iDaughter=0; iDaughter<particle.NDaughters(); iDaughter++)
          {
            int order[4] = {0, 1, 2, 3};
            const int daughterParticleIndex = particle.DaughterIds()[order[iDaughter]];
            KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterParticleIndex];

            if(iDaughter==0){
              chi2primary_om_bach = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());

              dedx_om_bach = fStKFParticleInterface->GetdEdX(daughter.DaughterIds()[0]);
              nhits_om_bach = fStKFParticleInterface->GetNHits(daughter.DaughterIds()[0]);

              bbachpx = daughter.GetPx();
              bbachpy = daughter.GetPy();
              bbachpz = daughter.GetPz();
            }
            if(iDaughter==1){
              chi2primary_om_ld = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
              KFParticleSIMD ttempSIMDParticle(daughter);
              float_v tl,tdl;

              KFParticleSIMD tpv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
              ttempSIMDParticle.GetDistanceToVertexLine(tpv, tl, tdl);
              om_ld_ldl = tl[0]/tdl[0];
              om_ld_l = tl[0];
              ttempSIMDParticle.SetProductionVertex(tpv);
              om_ld_chi2topo = double(ttempSIMDParticle.Chi2()[0])/double(ttempSIMDParticle.NDF()[0]);
              om_ld_chi2ndf = daughter.Chi2()/daughter.NDF();

              for(int jDaughter=0; jDaughter<daughter.NDaughters(); jDaughter++){
                int jorder[4] = {0, 1, 2, 3};
                const int jdaughterParticleIndex = daughter.DaughterIds()[jorder[jDaughter]];
                KFParticle granddaughter = fStKFParticleInterface->GetParticles()[jdaughterParticleIndex];
                if(jDaughter==0){
                  chi2primary_om_proton = granddaughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());

                  dedx_om_proton = fStKFParticleInterface->GetdEdX(granddaughter.DaughterIds()[0]);
                  nhits_om_proton = fStKFParticleInterface->GetNHits(granddaughter.DaughterIds()[0]);

                  bprotonpx = granddaughter.GetPx();
                  bprotonpy = granddaughter.GetPy();
                  bprotonpz = granddaughter.GetPz();

                }
                if(jDaughter==1){
                  chi2primary_om_pi = granddaughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
                  dedx_om_pi = fStKFParticleInterface->GetdEdX(granddaughter.DaughterIds()[0]);
                  nhits_om_pi = fStKFParticleInterface->GetNHits(granddaughter.DaughterIds()[0]);

                  bpionpx = granddaughter.GetPx();
                  bpionpy = granddaughter.GetPy();
                  bpionpz = granddaughter.GetPz();

                }
              }	
            }
          }
          bmcpx=-999;
          bmcpy=-999;
          bmcpz=-999;
          if(isMCParticle && fProcessSignal){
            int iMCPart = fStKFParticlePerformanceInterface->GetParticleMCId(iParticle);
            StMuMcTrack *mcTrack = fMuDst->MCtrack(iMCPart);
            bmcpx = mcTrack->Pxyz().x();
            bmcpy = mcTrack->Pxyz().y();
            bmcpz = mcTrack->Pxyz().z();
          }
          omega_tree->Fill();
        }

        //he3 tree
        if(fabs(particle.GetPDG())==1000020030 && isGoodEvent && _fill_he3_tree){
          int trackId = particle.DaughterIds()[0];//get the track id
          bdedx = fStKFParticleInterface->GetdEdX(trackId);
          bnhits = fStKFParticleInterface->GetNHits(trackId);
          bdca = fStKFParticleInterface->Getdca(trackId);
          chi2primary_he = particle.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());

          if(isMCParticle){
            int iMCPart = fStKFParticlePerformanceInterface->GetParticleMCId(iParticle);
            StMuMcTrack *mcTrack = fMuDst->MCtrack(iMCPart);
            bmcpx = mcTrack->Pxyz().x();
            bmcpy = mcTrack->Pxyz().y();
            bmcpz = mcTrack->Pxyz().z();
          }else{
            bmcpx=-999;
            bmcpy=-999;
            bmcpz=-999;
          }

          he3_tree->Fill();
        }

        if(  ( fabs(particle.GetPDG())==3004 || fabs(particle.GetPDG())==3005 || fabs(particle.GetPDG())==3103 )  && isGoodEvent){

          KFParticleSIMD tempSIMDParticle(particle);
          float_v l,dl;
          KFParticleSIMD pv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
          tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl);
          ht_ldl = l[0]/dl[0];
          ht_l = l[0];
          ht_dl = dl[0];
          bpl = ht_l/(sqrt(bpx*bpx+bpy*bpy+bpz*bpz)/bparticlemass);
          ht_bdfvtx = tempSIMDParticle.GetDistanceFromVertex(pv)[0];

          tempSIMDParticle.SetProductionVertex(pv);
          ht_bdfvtx2 = tempSIMDParticle.GetDistanceFromVertex(pv)[0];
          ht_lifetime = tempSIMDParticle.GetLifeTime()[0];
          ht_chi2topo = double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);
          ht_chi2 = double(tempSIMDParticle.Chi2()[0]);
          ht_NDF = double(tempSIMDParticle.NDF()[0]);
          ht_chi2ndf = particle.Chi2()/particle.NDF();

          for(int iDaughter=0; iDaughter<particle.NDaughters(); iDaughter++)
          {
            int order[4] = {0, 1, 2, 3};
            const int daughterParticleIndex = particle.DaughterIds()[order[iDaughter]];
            KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterParticleIndex];
            if(iDaughter==0)
            {
              chi2primary_pi = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
              nhits_pi = fStKFParticleInterface->GetNHits(daughter.DaughterIds()[0]);
              nhitsdedx_pi = fStKFParticleInterface->GetNHitsdedx(daughter.DaughterIds()[0]);
              dca_pi = fStKFParticleInterface->Getdca(daughter.DaughterIds()[0]);
              dedx_pi = fStKFParticleInterface->GetdEdX(daughter.DaughterIds()[0]);
              px_pi = daughter.GetPx();
              py_pi = daughter.GetPy();
              pz_pi = daughter.GetPz();
            }
            if(iDaughter==1)
            {
              //			cout<<"daughter:"<<sqrt(daughter.GetPx()*daughter.GetPx()+daughter.GetPy()*daughter.GetPy())<<endl;
              chi2primary_he = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
              nhits_he = fStKFParticleInterface->GetNHits(daughter.DaughterIds()[0]);
              nhitsdedx_he = fStKFParticleInterface->GetNHitsdedx(daughter.DaughterIds()[0]);
              dca_he = fStKFParticleInterface->Getdca(daughter.DaughterIds()[0]);
              dedx_he = fStKFParticleInterface->GetdEdX(daughter.DaughterIds()[0]);
              px_he = daughter.GetPx();
              py_he = daughter.GetPy();
              pz_he = daughter.GetPz();
            }
          }
          bmcpx=-999;
          bmcpy=-999;
          bmcpz=-999;
          bmcl=-999;
          bmcpl=-999;


          if(isMCParticle){
            int iMCPart = fStKFParticlePerformanceInterface->GetParticleMCId(iParticle);
            StMuMcTrack *mcTrack = fMuDst->MCtrack(iMCPart);
            bmcpx = mcTrack->Pxyz().x();
            bmcpy = mcTrack->Pxyz().y();
            bmcpz = mcTrack->Pxyz().z();
            bmcidvx = mcTrack->IdVx();
            bmcidvxend = mcTrack->IdVxEnd();
            StMuMcVertex *mcvx = fMuDst->MCvertex(bmcidvx-1);
            StMuMcVertex *mcvxend = fMuDst->MCvertex(bmcidvxend-1);

            bmcl = sqrt( (mcvx->XyzV().x()-mcvxend->XyzV().x())*(mcvx->XyzV().x()-mcvxend->XyzV().x()) + (mcvx->XyzV().y()-mcvxend->XyzV().y())*(mcvx->XyzV().y()-mcvxend->XyzV().y()) + (mcvx->XyzV().z()-mcvxend->XyzV().z())*(mcvx->XyzV().z()-mcvxend->XyzV().z()) ) ;
            if(fabs(particle.GetPDG())==3004){
              bmcpl = bmcl * 2.99131 / sqrt(bmcpx*bmcpx+bmcpy*bmcpy+bmcpz*bmcpz) ;
            }
            if(fabs(particle.GetPDG())==3005){
              bmcpl = bmcl * 3.9239 / sqrt(bmcpx*bmcpx+bmcpy*bmcpy+bmcpz*bmcpz) ;
            }
          }
          htriton_tree->Fill();
        }

        if((particle.GetPDG()==103004 ||particle.GetPDG()==3007 || particle.GetPDG()==103005 || particle.GetPDG()==3006) && isGoodEvent){

          if((particle.GetMass()>3.05 || particle.GetMass()<2.95)&& particle.GetPDG()==103004) continue;
          if(particle.GetMass()>3.2 && particle.GetPDG()==103005) continue;
          if(particle.GetMass()>5.1 && particle.GetPDG()==3007) continue;

          // cout<<"particle.NDaughters():"<<particle.NDaughters()<<endl;
          KFParticleSIMD tempSIMDParticle(particle);
          float_v l,dl;
          KFParticleSIMD pv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
          tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl);
          ht_ldl = l[0]/dl[0];
          ht_l = l[0];
          ht_dl = dl[0];

          ht_bdfvtx = tempSIMDParticle.GetDistanceFromVertex(pv)[0];
          tempSIMDParticle.SetProductionVertex(pv);
          ht_bdfvtx2 = tempSIMDParticle.GetDistanceFromVertex(pv)[0];

          ht_ndaughters = particle.NDaughters();
          ht_chi2topo = double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);
          ht_chi2 = double(tempSIMDParticle.Chi2()[0]);
          ht_NDF = double(tempSIMDParticle.NDF()[0]);
          ht_chi2ndf = particle.Chi2()/particle.NDF();
          // cout <<particle.NDF() << endl;

          for(int iDaughter=0; iDaughter<particle.NDaughters(); iDaughter++)
          {
            int order[3] = {0, 1, 2};
            if(particle.GetPDG()==3007) {order[1] =  2; order[2] = 1;};

            const int daughterParticleIndex = particle.DaughterIds()[order[iDaughter]];
            KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterParticleIndex];
            // cout << iDaughter<<" "<< daughter.GetPDG()<< endl;
            if(iDaughter==2){
              chi2primary_d = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
              bdpx = daughter.GetPx();
              bdpy = daughter.GetPy();
              bdpz = daughter.GetPz();
              int atrackId = daughter.DaughterIds()[0];//get the track id
              bdedx = fStKFParticleInterface->GetdEdX(atrackId);
              //btofm2 = fStKFParticleInterface->GetTOFm2(atrackId);
              bzdeuteron =  fStKFParticleInterface->Getzdeuteron(atrackId);
              nhits_deuteron = fStKFParticleInterface->GetNHits(daughter.DaughterIds()[0]);
              dca_deuteron = fStKFParticleInterface->Getdca(daughter.DaughterIds()[0]);
              bdm2 = fStKFParticleInterface->Getm2(daughter.DaughterIds()[0]);
            }
            if(iDaughter==0){
              chi2primary_pi = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
              bpionpx = daughter.GetPx();
              bpionpy = daughter.GetPy();
              bpionpz = daughter.GetPz();
              bpionnsigma = fStKFParticleInterface->GetdEdXNSigmaPion(daughter.DaughterIds()[0]);
              nhits_pion = fStKFParticleInterface->GetNHits(daughter.DaughterIds()[0]);
              dca_pion = fStKFParticleInterface->Getdca(daughter.DaughterIds()[0]);
              bpionm2 = fStKFParticleInterface->Getm2(daughter.DaughterIds()[0]);
            }
            if(iDaughter==1){
              chi2primary_proton = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
              bprotonpx = daughter.GetPx();
              bprotonpy = daughter.GetPy();
              bprotonpz = daughter.GetPz();
              bprotonsigma = fStKFParticleInterface->GetdEdXNSigmaProton(daughter.DaughterIds()[0]);
              nhits_proton = fStKFParticleInterface->GetNHits(daughter.DaughterIds()[0]);
              dca_proton = fStKFParticleInterface->Getdca(daughter.DaughterIds()[0]);
              bprotonm2 = fStKFParticleInterface->Getm2(daughter.DaughterIds()[0]);
            }
          }

          KFParticle pion     = fStKFParticleInterface->GetParticles()[particle.DaughterIds()[0]];
          KFParticle proton   = fStKFParticleInterface->GetParticles()[particle.DaughterIds()[1]];
          KFParticle bach     = fStKFParticleInterface->GetParticles()[particle.DaughterIds()[2]];
          // cout <<pion.GetPDG()<<" "<<proton.GetPDG()<<" "<<bach.GetPDG() << endl;


          float m, dm;
          float_v vl, vdl;

          KFParticle v_lambda;
          v_lambda += pion;
          v_lambda += proton;
          v_lambda.GetMass(m,dm);

          v_lambda_mass_0 = m;
          //v_lambda.GetDistanceToVertexLine(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex(), vl, vdl);
          // v_lambda_ldl_0 = vl/vdl;
          v_lambda_chi2primary_0 = v_lambda.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
          KFParticleSIMD tempLambdaSIMD(v_lambda);
          tempLambdaSIMD.GetDistanceToVertexLine(pv, vl, vdl);
          v_lambda_ldl_0 = vl[0]/vdl[0];
          v_lambda_l_0 = vl[0];

          const KFParticle* v_d01[2] = {&pion, &proton};
          const KFParticle* v_d02[2] = {&pion, &bach};
          const KFParticle* v_d12[2] = {&proton, &bach};

          KFParticle v_01;
          KFParticle v_02;
          KFParticle v_12;

          v_01.Construct(v_d01, 2, 0);
          v_01.GetMass(mass_01, mass_01_err);
          v_01_chi2primary = v_01.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
          v_01_pvdca = v_01.GetDistanceFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());

          v_02.Construct(v_d02, 2, 0);
          v_02.GetMass(mass_02, mass_02_err);
          v_02_chi2primary = v_02.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
          v_02_pvdca = v_02.GetDistanceFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());

          v_12.Construct(v_d12, 2, 0);
          v_12.GetMass(mass_12, mass_12_err);
          v_12_chi2primary = v_12.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
          v_12_pvdca = v_12.GetDistanceFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());

          //			float v_01_dca;
          //			float v_02_dca;
          //      float v_12_dca;

          //			float v_01_chi2ndf;
          //			float v_02_chi2ndf;
          //			float v_12_chi2ndf;

          v_01_chi2ndf = v_01.Chi2()/v_01.NDF();
          v_02_chi2ndf = v_02.Chi2()/v_02.NDF();
          v_12_chi2ndf = v_12.Chi2()/v_12.NDF();

          v_01_dca = pion.GetDistanceFromParticleXY(proton);
          v_02_dca = pion.GetDistanceFromParticleXY(bach);
          v_12_dca = proton.GetDistanceFromParticleXY(bach);

          v_012_dca = v_01.GetDistanceFromParticleXY(bach);

          //			cout<<"dca:"<<v_01_dca<<" "<<v_02_dca<<" "<<v_12_dca<<endl;
          //			cout<<"chi2ndf:"<<v_01_chi2ndf<<" "<<v_02_chi2ndf<<" "<<v_12_chi2ndf<<endl;

          dca_pion2 = pion.GetDistanceFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
          dca_proton2 = proton.GetDistanceFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
          dca_deuteron2 = bach.GetDistanceFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());

          // cout<<"dca_pion:"<<dca_pion<<" "<<dca_pion2<<endl;
          // cout<<"dca_proton:"<<dca_proton<<" "<<dca_proton2<<endl;
          // cout<<"dca_deuteron:"<<dca_deuteron-dca_deuteron2<<endl;

          bmcpx=-999;bmcpy=-999;bmcpz=-999;

          b0mcpx=-999;b0mcpy=-999;b0mcpz=-999;
          b1mcpx=-999;b1mcpy=-999;b1mcpz=-999;
          b2mcpx=-999;b2mcpy=-999;b2mcpz=-999;

          bool isMcDeu = fStKFParticlePerformanceInterface->GetParticle(bach, particle.DaughterIds()[2]);
          bool isMcPion = fStKFParticlePerformanceInterface->GetParticle(pion, particle.DaughterIds()[0]);
          bool isMcProton = fStKFParticlePerformanceInterface->GetParticle(proton, particle.DaughterIds()[1]);
          // cout << "check MC"<<isMcPion<<" "<<isMcProton <<" "<< isMcDeu<<endl;

          bismc=0;
          //check if the daughter coming from the same vertex
          if (isMcProton && isMcPion )
          {
            int iMcPion = fStKFParticlePerformanceInterface->GetParticleMCId(particle.DaughterIds()[0]);
            int iMcProton = fStKFParticlePerformanceInterface->GetParticleMCId(particle.DaughterIds()[1]);
            StMuMcTrack *mcProton= fMuDst->MCtrack(iMcProton);
            StMuMcTrack *mcPion = fMuDst->MCtrack(iMcPion);
            // cout <<mcPion->GePid()<<" "<<mcProton->GePid() << endl;
            bool correctDau =mcProton->GePid()==14 && mcPion->GePid()==9; 
            bool isLambdaDau = (mcProton->IdVx() == mcPion->IdVx()) && correctDau;
            // cout << isLambdaDau << " "<<isMcDeu <<endl;

            // 9 -> pion-
            b0mcpx = mcPion->Pxyz().x();
            b0mcpy = mcPion->Pxyz().y();
            b0mcpz = mcPion->Pxyz().z();

            // 14 -> proton
            b1mcpx = mcProton->Pxyz().x();
            b1mcpy = mcProton->Pxyz().y();
            b1mcpz = mcProton->Pxyz().z();

            if (isMcDeu && !isLambdaDau) bismc=-111; //no use rightnow
            if (isMcDeu && isLambdaDau  )
            {
              int iMcDeu = fStKFParticlePerformanceInterface->GetParticleMCId(particle.DaughterIds()[2]);
              StMuMcTrack *mcDeu = fMuDst->MCtrack(iMcDeu);

              // 45 -> DEUTERON
              b2mcpx = mcDeu->Pxyz().x();
              b2mcpy = mcDeu->Pxyz().y();
              b2mcpz = mcDeu->Pxyz().z();


              int idH3L=-1, idLb=-1; 
              //phase space case
              if (mcDeu->IdVx()==mcProton->IdVx() && mcDeu->IdVx()==mcPion->IdVx() &&  mcDeu->GePid()==45 ) 
              {
                //decay from the same particle
                for (Int_t j = 0; j < fMuDst->numberOfMcTracks(); j++) {
                  StMuMcTrack* mcTrack = fMuDst->MCtrack(j);
                  if (mcTrack->GePid()==62053 && mcTrack->IdVxEnd() == mcDeu->IdVx())  
                    // if ( mcTrack->IdVxEnd() == mcDeu->IdVx()) 
                  { 
                    // cout << "check if find the real mother: "<<mcTrack->GePid()<<" pion: "<<mcPion->GePid()<<" mcProton: "<<mcProton->GePid()<<endl; 
                    idH3L= j; 
                    break;
                  }
                }
                if (idH3L>=0) bismc=1;
              } //end of phase space  
              else if ( mcDeu->GePid()==45 )  { //quasi-Lambda decay, search Lambda
                // cout <<"check quasi ";
                for (Int_t j = 0; j < fMuDst->numberOfMcTracks(); j++) 
                {
                  StMuMcTrack* tmpTrack = fMuDst->MCtrack(j);
                  if ( tmpTrack->IdVxEnd() == mcDeu->IdVx() && tmpTrack->GePid()==63053)
                  { 
                    // cout << "check if find the real mother: "<<
                    // tmpTrack->GePid()<<" pion: "<<mcPion->GePid()<<
                    // " mcProton: "<<mcProton->GePid()<<" "<<mcDeu->GePid()<<endl; 
                    idH3L = j;
                  }
                  else if ( tmpTrack->IdVxEnd() == mcProton->IdVx() && tmpTrack->IdVx()==mcDeu->IdVx() && tmpTrack->GePid()==95)
                  { 
                    // cout << "check if find the real Lambda: "<<tmpTrack->GePid()<<endl; 
                    idLb = j;
                  }
                  if (idLb>=0 && idH3L>=0) break;
                }

                if (idH3L>=0 && idLb>=0) bismc=1;
                else bismc=-21;

              } //end of quasi

              if (bismc) {

                StMuMcTrack* mcTrack = fMuDst->MCtrack(idH3L);
                bmcpx = mcTrack->Pxyz().x();
                bmcpy = mcTrack->Pxyz().y();
                bmcpz = mcTrack->Pxyz().z();
                bmcidvx = mcTrack->IdVx();
                bmcidvxend = mcTrack->IdVxEnd();
                StMuMcVertex *mcvx = fMuDst->MCvertex(bmcidvx-1);
                StMuMcVertex *mcvxend = fMuDst->MCvertex(bmcidvxend-1);

                bmcl = sqrt( (mcvx->XyzV().x()-mcvxend->XyzV().x())*(mcvx->XyzV().x()-mcvxend->XyzV().x()) 
                    + (mcvx->XyzV().y()-mcvxend->XyzV().y())*(mcvx->XyzV().y()-mcvxend->XyzV().y()) 
                    + (mcvx->XyzV().z()-mcvxend->XyzV().z())*(mcvx->XyzV().z()-mcvxend->XyzV().z()) ) ;
                if (particle.GetPDG()==103004) 
                {
                  bmcpl = bmcl * 2.99131 / sqrt(bmcpx*bmcpx+bmcpy*bmcpy+bmcpz*bmcpz) ;
                }
                else {
                  bmcpl = bmcl * 4.86824 / sqrt(bmcpx*bmcpx+bmcpy*bmcpy+bmcpz*bmcpz) ;
                }
              } // fill mc parameters
            }
            else if (!isMcDeu) {
              if (isLambdaDau) bismc = -20;
              else bismc=0;
            } // random d from data and Mc H3L->lambda  

          } // end of mc pion and proton 

          htriton3_tree->Fill();

        }

        /*
           if(fabs(particle.GetPDG())==3005 && isGoodEvent){

           KFParticleSIMD tempSIMDParticle(particle);
           float_v l,dl;
           KFParticleSIMD pv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
           tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl);
           hl_ldl = l[0]/dl[0];
           hl_l = l[0];
           hl_dl = dl[0];

           tempSIMDParticle.SetProductionVertex(pv);
           hl_chi2topo = double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);
           hl_chi2 = double(tempSIMDParticle.Chi2()[0]);
           hl_NDF = double(tempSIMDParticle.NDF()[0]);
           hl_chi2ndf = particle.Chi2()/particle.NDF();

           for(int iDaughter=0; iDaughter<particle.NDaughters(); iDaughter++)
           {
           int order[4] = {0, 1, 2, 3};
           const int daughterParticleIndex = particle.DaughterIds()[order[iDaughter]];
           KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterParticleIndex];
           if(iDaughter==0){chi2primary_h4 = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());}
           if(iDaughter==1){chi2primary_pi = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());}

           }

           h4lambda_tree->Fill();
           }
           */
        //if((fabs(particle.GetPDG())==3006 || fabs(particle.GetPDG())==3009  || fabs(particle.GetPDG())==3007) && isGoodEvent)
        if( fabs(particle.GetPDG())==3009 && isGoodEvent)
        {

          KFParticleSIMD tempSIMDParticle(particle);
          float_v l,dl;
          KFParticleSIMD pv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
          tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl);
          hl_l = l[0];
          hl_dl = dl[0];
          hl_ldl = l[0]/dl[0];

          tempSIMDParticle.SetProductionVertex(pv);
          hl_chi2topo = double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);
          hl_chi2ndf = particle.Chi2()/particle.NDF();

          for(int iDaughter=0; iDaughter<particle.NDaughters(); iDaughter++)
          {
            int order[3] = {0, 2, 1};
            const int daughterParticleIndex = particle.DaughterIds()[order[iDaughter]];
            KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterParticleIndex];
            if(iDaughter==2){
              chi2primary_d = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
              bdpx = daughter.GetPx();
              bdpy = daughter.GetPy();
              bdpz = daughter.GetPz();
              int atrackId = daughter.DaughterIds()[0];//get the track id
              bdedx = fStKFParticleInterface->GetdEdX(atrackId);
              bzdeuteron =  fStKFParticleInterface->Getzdeuteron(atrackId);
              nhits_deuteron = fStKFParticleInterface->GetNHits(daughter.DaughterIds()[0]);
              dca_deuteron = fStKFParticleInterface->Getdca(daughter.DaughterIds()[0]);
              bdm2 = fStKFParticleInterface->Getm2(daughter.DaughterIds()[0]);
            }
            if(iDaughter==0){
              chi2primary_pi = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
              bpionpx = daughter.GetPx();
              bpionpy = daughter.GetPy();
              bpionpz = daughter.GetPz();
              bpionnsigma = fStKFParticleInterface->GetdEdXNSigmaPion(daughter.DaughterIds()[0]);
              nhits_pion = fStKFParticleInterface->GetNHits(daughter.DaughterIds()[0]);
              dca_pion = fStKFParticleInterface->Getdca(daughter.DaughterIds()[0]);
              bpionm2 = fStKFParticleInterface->Getm2(daughter.DaughterIds()[0]);
            }
            if(iDaughter==1){
              chi2primary_proton = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
              bprotonpx = daughter.GetPx();
              bprotonpy = daughter.GetPy();
              bprotonpz = daughter.GetPz();
              bprotonsigma = fStKFParticleInterface->GetdEdXNSigmaProton(daughter.DaughterIds()[0]);
              nhits_proton = fStKFParticleInterface->GetNHits(daughter.DaughterIds()[0]);
              dca_proton = fStKFParticleInterface->Getdca(daughter.DaughterIds()[0]);
              bprotonm2 = fStKFParticleInterface->Getm2(daughter.DaughterIds()[0]);
            }
          }

          KFParticle pion     = fStKFParticleInterface->GetParticles()[particle.DaughterIds()[0]];
          KFParticle proton   = fStKFParticleInterface->GetParticles()[particle.DaughterIds()[1]];
          KFParticle bach     = fStKFParticleInterface->GetParticles()[particle.DaughterIds()[2]];
          float m, dm, vl, vdl;

          KFParticle v_lambda;
          v_lambda += pion;
          v_lambda += proton;
          v_lambda.GetMass(m,dm);
          v_lambda_mass_0 = m;
          v_lambda_ldl_0 = vl/vdl;
          v_lambda_chi2primary_0 = v_lambda.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());

          const KFParticle* v_d01[2] = {&pion, &proton};
          const KFParticle* v_d02[2] = {&pion, &bach};
          const KFParticle* v_d12[2] = {&proton, &bach};

          KFParticle v_01;
          KFParticle v_02;
          KFParticle v_12;

          v_01.Construct(v_d01, 2, 0);
          v_01.GetMass(mass_01, mass_01_err);
          v_01_chi2primary = v_01.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
          v_01_pvdca = v_01.GetDistanceFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
          v_02.Construct(v_d02, 2, 0);
          v_02.GetMass(mass_02, mass_02_err);
          v_02_chi2primary = v_02.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
          v_02_pvdca = v_02.GetDistanceFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());

          v_12.Construct(v_d12, 2, 0);
          v_12.GetMass(mass_12, mass_12_err);
          v_12_chi2primary = v_12.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
          v_12_pvdca = v_12.GetDistanceFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
          v_01_chi2ndf = v_01.Chi2()/v_01.NDF();
          v_02_chi2ndf = v_02.Chi2()/v_02.NDF();
          v_12_chi2ndf = v_12.Chi2()/v_12.NDF();

          v_01_dca = pion.GetDistanceFromParticleXY(proton);
          v_02_dca = pion.GetDistanceFromParticleXY(bach);
          v_12_dca = proton.GetDistanceFromParticleXY(bach);

          dca_pion2 = pion.GetDistanceFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
          dca_proton2 = proton.GetDistanceFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
          dca_deuteron2 = bach.GetDistanceFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());

          h_tree->Fill();
        }
      }//loop through paritcles

    }//write tree bool

    //cout<<"fStoreTmvaNTuples:"<<fStoreTmvaNTuples<<" "<< fStKFParticlePerformanceInterface->GetNReconstructedParticles()<<endl;


    if(fStoremctree && isGoodEvent){

      Int_t NoMuMcTracks = fMuDst->numberOfMcTracks();
      for (Int_t k = 0; k < NoMuMcTracks; k++) {

        StMuMcTrack *mcTrack = fMuDst->MCtrack(k);

        bmcparticleid = mcTrack->GePid();
        // cout << bmcparticleid<<endl;
        if(bmcparticleid==40002){
          bmcparticleid = 3334;
        }
        if(bmcparticleid==99999){
          bmcparticleid = 3312;
        }
        if(bmcparticleid==10018 || bmcparticleid==18){
          bmcparticleid = 3122;
        }
        if(bmcparticleid==49){
          bmcparticleid = 1000020030;
        }
        if(bmcparticleid==61053){	
          bmcparticleid = 3004;
        }
        if(bmcparticleid==61055){
          bmcparticleid = 3005;
        }
        if(bmcparticleid==707){
          bmcparticleid = 310;
        }
        if(bmcparticleid==62053 || bmcparticleid==63053){
          bmcparticleid = 103004;
        }
        if(bmcparticleid==63063||bmcparticleid==62063||bmcparticleid==61057){
          bmcparticleid = 3006;
        }
        bmcrawpx = mcTrack->Pxyz().x();
        bmcrawpy = mcTrack->Pxyz().y();
        bmcrawpz = mcTrack->Pxyz().z();		
        bmcrefmult = brefmult;

        if(abs(bmcparticleid)==3334){
          omega_mc_tree->Fill();
        }
        if(abs(bmcparticleid)==3312){
          cascade_mc_tree->Fill();
        }
        if(abs(bmcparticleid)==3122){
          bmcidvx = mcTrack->IdVx();
          bmcidvxend = mcTrack->IdVxEnd();
          if(bmcidvx!=1) continue;
          if(bmcidvxend==0) continue;
          StMuMcVertex *mcvx = fMuDst->MCvertex(bmcidvx-1);
          StMuMcVertex *mcvxend = fMuDst->MCvertex(bmcidvxend-1);
          bmcrawl = sqrt( (mcvx->XyzV().x()-mcvxend->XyzV().x())*(mcvx->XyzV().x()-mcvxend->XyzV().x()) + (mcvx->XyzV().y()-mcvxend->XyzV().y())*(mcvx->XyzV().y()-mcvxend->XyzV().y()) + (mcvx->XyzV().z()-mcvxend->XyzV().z())*(mcvx->XyzV().z()-mcvxend->XyzV().z()) );
          bmcrawpl = bmcrawl* 1.115683 / sqrt(bmcrawpx*bmcrawpx+bmcrawpy*bmcrawpy+bmcrawpz*bmcrawpz) ;
          lambda_mc_tree->Fill();
        }
        if(bmcparticleid==310){
          bmcidvx = mcTrack->IdVx();
          bmcidvxend = mcTrack->IdVxEnd();
          if(bmcidvx!=1) continue;
          if(bmcidvxend==0) continue;
          StMuMcVertex *mcvx = fMuDst->MCvertex(bmcidvx-1);
          StMuMcVertex *mcvxend = fMuDst->MCvertex(bmcidvxend-1);
          bmcrawl = sqrt( (mcvx->XyzV().x()-mcvxend->XyzV().x())*(mcvx->XyzV().x()-mcvxend->XyzV().x()) + (mcvx->XyzV().y()-mcvxend->XyzV().y())*(mcvx->XyzV().y()-mcvxend->XyzV().y()) + (mcvx->XyzV().z()-mcvxend->XyzV().z())*(mcvx->XyzV().z()-mcvxend->XyzV().z()) );
          bmcrawpl = bmcrawl* 0.497611 / sqrt(bmcrawpx*bmcrawpx+bmcrawpy*bmcrawpy+bmcrawpz*bmcrawpz) ;
          ks_mc_tree->Fill();
        }
        if(abs(bmcparticleid)==1000020030){
          he3_mc_tree->Fill();
        }

        if (abs(bmcparticleid)==3004 || abs(bmcparticleid)==3005 || bmcparticleid==103004 || bmcparticleid==3006){
          cout <<"find 103004" << endl;
          bmcidvx = mcTrack->IdVx();
          bmcidvxend = mcTrack->IdVxEnd();
          if(bmcidvx!=1) continue;
          if(bmcidvxend==0) continue;
          StMuMcVertex *mcvx = fMuDst->MCvertex(bmcidvx-1);
          StMuMcVertex *mcvxend = fMuDst->MCvertex(bmcidvxend-1);
          bmcrawl = sqrt( (mcvx->XyzV().x()-mcvxend->XyzV().x())*(mcvx->XyzV().x()-mcvxend->XyzV().x()) + (mcvx->XyzV().y()-mcvxend->XyzV().y())*(mcvx->XyzV().y()-mcvxend->XyzV().y()) + (mcvx->XyzV().z()-mcvxend->XyzV().z())*(mcvx->XyzV().z()-mcvxend->XyzV().z()) );
          if(abs(bmcparticleid)==3004 || bmcparticleid==103004){
            bmcrawpl = bmcrawl* 2.99131/ sqrt(bmcrawpx*bmcrawpx+bmcrawpy*bmcrawpy+bmcrawpz*bmcrawpz) ;
          }
          if(abs(bmcparticleid)==3005){
            bmcrawpl = bmcrawl* 3.924/ sqrt(bmcrawpx*bmcrawpx+bmcrawpy*bmcrawpy+bmcrawpz*bmcrawpz) ;
          }

          if(bmcparticleid==103004){
            for (Int_t j = 0; j < fMuDst->numberOfMcTracks(); j++) {
              StMuMcTrack *dmcTrack = fMuDst->MCtrack(j);

              //	cout<<"dmcTrack->GePid():"<<dmcTrack->GePid()<<" "<< bmcidvxend<<" "<<dmcTrack->IdVx()<<endl;

              if(dmcTrack->IdVx() == bmcidvxend){
                //		cout<<"dmcTrack->GePid():"<<dmcTrack->GePid()<<endl;

                if(dmcTrack->GePid()==45){//deuteron
                  b2mcrawpx = dmcTrack->Pxyz().x();
                  b2mcrawpy = dmcTrack->Pxyz().y();
                  b2mcrawpz = dmcTrack->Pxyz().z();
                }
                else if(dmcTrack->GePid()==95){//"lambda"
                  for (Int_t k = 0; k < fMuDst->numberOfMcTracks(); k++) {
                    StMuMcTrack *ddmcTrack = fMuDst->MCtrack(k);//another loop

                    //		cout<<"ddmcTrack->GePid():"<<ddmcTrack->GePid()<<endl;
                    //		cout<<"ddmcTrack->IdVx():"<<ddmcTrack->IdVx()<<" "<<dmcTrack->IdVxEnd()<<" "<<endl;

                    if(ddmcTrack->IdVx() == dmcTrack->IdVxEnd()){

                      if(ddmcTrack->GePid()==9){//pion
                        b0mcrawpx = ddmcTrack->Pxyz().x();
                        b0mcrawpy = ddmcTrack->Pxyz().y();
                        b0mcrawpz = ddmcTrack->Pxyz().z();
                      }
                      if(ddmcTrack->GePid()==14){//proton
                        b1mcrawpx = ddmcTrack->Pxyz().x();
                        b1mcrawpy = ddmcTrack->Pxyz().y();
                        b1mcrawpz = ddmcTrack->Pxyz().z();
                      }

                    }
                  }		
                }
                else {
                  //the default 3body decay setup
                  if(dmcTrack->GePid()==9){//pion
                    b0mcrawpx = dmcTrack->Pxyz().x();
                    b0mcrawpy = dmcTrack->Pxyz().y();
                    b0mcrawpz = dmcTrack->Pxyz().z();
                  }
                  if(dmcTrack->GePid()==14){//proton
                    b1mcrawpx = dmcTrack->Pxyz().x();
                    b1mcrawpy = dmcTrack->Pxyz().y();
                    b1mcrawpz = dmcTrack->Pxyz().z();
                  }
                }  
              }
            }
          }//loop over daughter	

          htriton_mc_tree->Fill();
        }

      }//mctrack loop

    }



    if(fStoreTmvaNTuples)
    {
      for(int iParticle=0; iParticle<fStKFParticlePerformanceInterface->GetNReconstructedParticles(); iParticle++)
      {
        KFParticle particle;
        bool isMCParticle = fStKFParticlePerformanceInterface->GetParticle(particle, iParticle);

        //       cout<<"isMCParticle:"<<isMCParticle<< " "<< particle.GetPDG()<<endl;      

        //TOBEREVERTED //REVERTED
        if( !( (fProcessSignal && isMCParticle) || (!fProcessSignal && !isMCParticle) ) ) continue;
        //if(isMCParticle) continue;          

        for(int iNTuple=0; iNTuple<fNNTuples; iNTuple++)
        {
          if( particle.GetPDG() == fNTuplePDG[iNTuple] )
          {
            GetParticleParameters(iNTuple, particle);
            double side_mass = particle.GetMass();

            //bool sideband = (side_mass>1.06 && side_mass<1.10) || (side_mass>1.13 && side_mass<1.17);

            bool sideband;
            if(iNTuple==0){ //lambda
              sideband = (side_mass>1.085 && side_mass<1.10) || (side_mass>1.125 && side_mass<1.14);//for consistency with long code
              sideband = true;//for testing
            }
            if(iNTuple==1){ //omega

              //	if(isMCParticle){

              //	  int iMCPart = fStKFParticlePerformanceInterface->GetParticleMCId(iParticle);

              //cout<<"mcPartid:"<<iMCPart<<endl;
              //KFMCParticle &mcPart = vMCParticles[iMCPart];
              //cout<<"mcPart:"<<mcPart.GetPDG()<<endl;
              // mcTracks 
              //StMuMcTrack *mcTrack = fMuDst->MCtrack(iMCPart);
              //double mcpt = mcTrack->Pxyz().x();
              //double mcpt = mcTrack->pT;
              //cout<<"mcPt:"<<mcpt<<" "<< particle.GetPx()<<endl;

              //	}

              //sideband = (side_mass>1.61 && side_mass<1.66) || (side_mass>1.685 && side_mass<1.85);//for testing
              sideband = (side_mass>1.61 && side_mass<1.66) || (side_mass>1.685 && side_mass<1.73);//for consistency with long code

              //TOBEREVERTED
              //  sideband = true;//for testing
              sideband = side_mass<1.78;
              //            cout<<"fProcessSignal:"<<fProcessSignal<<" isMCParticle:"<<isMCParticle<<endl;
            }
            if(iNTuple==2){
              sideband = (side_mass>1.25 && side_mass<1.30) || (side_mass>1.35 && side_mass<1.40);//for consistency with long code
              //TOBEREVERTED
              sideband = true;
            }
            if( fProcessSignal || (!fProcessSignal && sideband) )
              fCutsNTuple[iNTuple]->Fill(fTMVAParticleParameters[iNTuple].data());
          }
        }
      }
    }

  }

  return kStOK;

}

void StKFParticleAnalysisMaker::GetDaughterParameters(const int iReader, int& iDaughterTrack, int& iDaughterParticle, KFParticle& particle)
{
  if(particle.NDaughters() == 1)
  {
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts]   = particle.GetPt();
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+1] = particle.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
    int trackId = particle.DaughterIds()[0];
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+2]   = fStKFParticleInterface->GetdEdXNSigmaPion(trackId);
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+3]   = fStKFParticleInterface->GetdEdXNSigmaKaon(trackId);
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+4]   = fStKFParticleInterface->GetdEdXNSigmaProton(trackId);
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+5]   = fStKFParticleInterface->GetTofNSigmaPion(trackId);
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+6]   = fStKFParticleInterface->GetTofNSigmaKaon(trackId);
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+7]   = fStKFParticleInterface->GetTofNSigmaProton(trackId);

    iDaughterTrack++;
  }
  else if(particle.NDaughters() > 1)
  {
    int order[4] = {0, 1, 2, 3};
    if( particle.GetPDG() == -421 || particle.GetPDG() == -411 || particle.GetPDG() == -431 ||   
        particle.GetPDG() == -429 || particle.GetPDG() == -4122) 
    { 
      order[0] = 1; 
      order[1] = 0; 
    }

    for(int iDaughter=0; iDaughter<particle.NDaughters(); iDaughter++)
    {
      const int daughterParticleIndex = particle.DaughterIds()[order[iDaughter]];
      KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterParticleIndex];
      //set pdg for correct order of cuts
      if(particle.GetPDG() == 521 && daughter.GetPDG() == -1) daughter.SetPDG(-421);
      if(particle.GetPDG() ==-521 && daughter.GetPDG() == -1) daughter.SetPDG( 421);
      if(particle.GetPDG() == 511 && daughter.GetPDG() == -1) daughter.SetPDG(-411);
      if(particle.GetPDG() ==-511 && daughter.GetPDG() == -1) daughter.SetPDG( 411);

      GetDaughterParameters(iReader, iDaughterTrack, iDaughterParticle, daughter);
    }

    fTMVAParticleParameters[iReader][fDaughterNames[iReader].size()*fNTrackTMVACuts + iDaughterParticle*3] = particle.Chi2()/particle.NDF();  

    KFParticleSIMD tempSIMDParticle(particle);
    float_v l,dl;
    KFParticleSIMD pv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
    tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl);
    fTMVAParticleParameters[iReader][fDaughterNames[iReader].size()*fNTrackTMVACuts + iDaughterParticle*3 + 1] = l[0]/dl[0];

    tempSIMDParticle.SetProductionVertex(pv);
    fTMVAParticleParameters[iReader][fDaughterNames[iReader].size()*fNTrackTMVACuts + iDaughterParticle*3 + 2] = 
      double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);

    iDaughterParticle++;
  }
}

void StKFParticleAnalysisMaker::GetParticleParameters(const int iReader, KFParticle& particle)
{
  bool isBMeson = abs(particle.GetPDG()) == 511 || abs(particle.GetPDG()) == 521;
  //   if( !isBMeson ) return;

  int iDaughterTrack = 0;
  int iDaughterParticle = 0;
  GetDaughterParameters(iReader, iDaughterTrack, iDaughterParticle, particle);

  int nDaughterParticleCut = 0;
  if(isBMeson) nDaughterParticleCut += 3;
  nDaughterParticleCut += fDaughterNames[iReader].size()*fNTrackTMVACuts;

  fTMVAParticleParameters[iReader][nDaughterParticleCut]   = particle.Chi2()/particle.NDF();  

  KFParticleSIMD tempSIMDParticle(particle);
  float_v l,dl;
  KFParticleSIMD pv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
  tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl);
  fTMVAParticleParameters[iReader][nDaughterParticleCut + 1] = l[0]/dl[0];

  fTMVAParticleParameters[iReader][nDaughterParticleCut + 2] = l[0];

  tempSIMDParticle.SetProductionVertex(pv);
  fTMVAParticleParameters[iReader][nDaughterParticleCut + 3] = double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);

  if(fIsPicoAnalysis)
    fTMVAParticleParameters[iReader][nDaughterParticleCut + 4] = fPicoDst->event()->refMult();
  else
    fTMVAParticleParameters[iReader][nDaughterParticleCut + 4] = fMuDst->event()->refMult();

  fTMVAParticleParameters[iReader][nDaughterParticleCut + 5] = refmultcor;

  fTMVAParticleParameters[iReader][nDaughterParticleCut + 6] = reweight;

  fTMVAParticleParameters[iReader][nDaughterParticleCut + 7] = cent9;

  fTMVAParticleParameters[iReader][nDaughterParticleCut + 8] = particle.GetMass();

  fTMVAParticleParameters[iReader][nDaughterParticleCut + 9] = particle.GetPx();

  fTMVAParticleParameters[iReader][nDaughterParticleCut + 10] = particle.GetPy();

  fTMVAParticleParameters[iReader][nDaughterParticleCut + 11] = particle.GetPz();

  fTMVAParticleParameters[iReader][nDaughterParticleCut + 12] = bVx;

  fTMVAParticleParameters[iReader][nDaughterParticleCut + 13] = bVy;

  fTMVAParticleParameters[iReader][nDaughterParticleCut + 14] = bVz;

  fTMVAParticleParameters[iReader][nDaughterParticleCut + 15] = bVxerr;

  fTMVAParticleParameters[iReader][nDaughterParticleCut + 16] = bVyerr;

  fTMVAParticleParameters[iReader][nDaughterParticleCut + 17] = bVzerr;

  if(iReader>0){//cascade or omega
    int order[4] = {0, 1, 2, 3};
    const int daughterParticleIndex = particle.DaughterIds()[order[1]];
    KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterParticleIndex];

    fTMVAParticleParameters[iReader][nDaughterParticleCut + 18]   = daughter.Chi2()/daughter.NDF();

    KFParticleSIMD ttempSIMDParticle(daughter);
    float_v tl,tdl;
    KFParticleSIMD tpv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
    ttempSIMDParticle.GetDistanceToVertexLine(tpv, tl, tdl);

    fTMVAParticleParameters[iReader][nDaughterParticleCut + 19] = tl[0]/tdl[0];

    fTMVAParticleParameters[iReader][nDaughterParticleCut + 20] = tl[0];  

    ttempSIMDParticle.SetProductionVertex(tpv); 
    fTMVAParticleParameters[iReader][nDaughterParticleCut + 21] = double(ttempSIMDParticle.Chi2()[0])/double(ttempSIMDParticle.NDF()[0]);

    fTMVAParticleParameters[iReader][nDaughterParticleCut + 22] = particle.GetPDG();

  }
}

Int_t StKFParticleAnalysisMaker::Finish() 
{
  if(fStoreTmvaNTuples)
  {
    TFile* curFile = gFile;
    TDirectory* curDirectory = gDirectory;
    for(int iNtuple=0; iNtuple<fNNTuples; iNtuple++)
    {
      fNTupleFile[iNtuple]->cd();
      fCutsNTuple[iNtuple]->Write();
    }
    gFile = curFile;
    gDirectory = curDirectory;
  }

  //  file_out = new TFile("lambda_tree.root","RECREATE");
  file_out->cd();
  //  lambda_tree->Write();
  file_out->Write();
  delete lambda_tree;
  delete file_out;


  return kStOK;
}

long StKFParticleAnalysisMaker::GetUniqueEventId(const int iRun, const int iEvent) const
{
  long id = 1000000000;
  return id*(iRun%1000) + iEvent;
}

int StKFParticleAnalysisMaker::GetTMVACentralityBin(int iReader, int centrality)
{
  for(unsigned int iBin=0; iBin<fTMVACentralityBins[iReader].size()-1; iBin++)
    if(centrality >= fTMVACentralityBins[iReader][iBin] && centrality < fTMVACentralityBins[iReader][iBin+1])
      return iBin;
  return -1;
}

int StKFParticleAnalysisMaker::GetTMVAPtBin(int iReader, double pt)
{
  for(unsigned int iBin=0; iBin<fTMVAPtBins[iReader].size()-1; iBin++)
    if(pt >= fTMVAPtBins[iReader][iBin] && pt < fTMVAPtBins[iReader][iBin+1])
      return iBin;
  return -1;
}

void StKFParticleAnalysisMaker::SetTMVACentralityBins(int iReader, TString bins)
{
  fTMVACentralityBins[iReader].clear();
  TString value; int firstSymbol = 0;      
  while(bins.Tokenize(value,firstSymbol,":"))
    fTMVACentralityBins[iReader].push_back(value.Atoi());
}

void StKFParticleAnalysisMaker::SetTMVAPtBins(int iReader, TString bins)
{
  fTMVAPtBins[iReader].clear();
  TString value; int firstSymbol = 0;      
  while(bins.Tokenize(value,firstSymbol,":"))
    fTMVAPtBins[iReader].push_back(value.Atof());
}

void StKFParticleAnalysisMaker::SetTMVABins(int iReader, TString centralityBins, TString ptBins)
{
  SetTMVACentralityBins(iReader, centralityBins);
  SetTMVAPtBins(iReader, ptBins);

  const int nCentralityBins = fTMVACentralityBins[iReader].size() - 1;
  const int nPtBins = fTMVAPtBins[iReader].size() - 1;

  fTMVACutFile[iReader].resize(nCentralityBins);
  fTMVACut[iReader].resize(nCentralityBins);
  fTMVAReader[iReader].resize(nCentralityBins);

  for(int iCentralityBin=0; iCentralityBin<nCentralityBins; iCentralityBin++)
  {
    fTMVACutFile[iReader][iCentralityBin].resize(nPtBins);
    fTMVACut[iReader][iCentralityBin].resize(nPtBins);
    fTMVAReader[iReader][iCentralityBin].resize(nPtBins);
  }
}

void StKFParticleAnalysisMaker::FillDaughterInfo(KFParticle& daughter, int trackId, int PDG, float &chi2primary, int &nhits, float &dca,float& sigma, float& bm2)
{
  chi2primary= daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
  if (PDG==211) 
    sigma = fStKFParticleInterface->GetdEdXNSigmaPion(trackId);
  else if (PDG==2212)
    sigma = fStKFParticleInterface->GetdEdXNSigmaProton(trackId);
  else if (PDG==1000010020)
    sigma = fStKFParticleInterface->Getzdeuteron(trackId);
  else sigma=-999;
  nhits= fStKFParticleInterface->GetNHits(trackId);
  dca= fStKFParticleInterface->Getdca(trackId);
  bm2 = fStKFParticleInterface->Getm2(trackId);
}
bool StKFParticleAnalysisMaker::FillThreeDaughtersMix(KFParticle& pion, KFParticle& proton, KFParticle& deuteron, KFParticle& PPi, double dVx, double dVy, double dVz, int mode)
{
  float fDistanceCut=5, fLCut=1;
  float cut[3]={3,10,10};  // ldl, chi2topo, chi2geo

  if (mode == 2){
    PPi.X()=PPi.GetX()+dVx;
    PPi.Y() = PPi.GetY()+dVy;
    PPi.Z() = PPi.GetZ()+dVz;
    
    pion.X() = pion.GetX()+dVx;
    pion.Y() = pion.GetY()+dVy;
    pion.Z() = pion.GetZ()+dVz;

    proton.X() = proton.GetX()+dVx;
    proton.Y() = proton.GetY()+dVy;
    proton.Z() = proton.GetZ()+dVz;
  }
  else if (mode ==1){
    deuteron.X() = deuteron.GetX()+dVx;
    deuteron.Y() = deuteron.GetY()+dVy;
    deuteron.Z() = deuteron.GetZ()+dVz;
  }

  //pip pair
  KFParticle v_lambda(PPi);
  float m, dm;
  float_v vl, vdl;
  v_lambda.GetMass(m,dm);
  v_lambda_mass_0 = m;

  KFParticleSIMD tempLambdaSIMD(v_lambda);
  KFParticleSIMD pv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
  tempLambdaSIMD.GetDistanceToVertexLine(pv, vl, vdl);
  v_lambda_ldl_0 = vl[0]/vdl[0];
  v_lambda_l_0 = vl[0];
  //v_lambda.GetDistanceToVertexLine(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex(), vl, vdl);
  // v_lambda_ldl_0 = vl/vdl;
  v_lambda_chi2primary_0 = v_lambda.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());

  //construct pipd
  KFParticleSIMD tempDeutSIMD(deuteron);
  KFParticleSIMD tempProtonSIMD(proton);
  KFParticleSIMD tempPiSIMD(pion);
  KFParticleSIMD tempSIMDParticle = tempPiSIMD;
  tempSIMDParticle += tempProtonSIMD;
  tempSIMDParticle.SetPDG(deuteron.GetPDG()>0? 103004 : -103004);
  tempDeutSIMD.TransportToPoint(tempSIMDParticle.Parameters());
  tempSIMDParticle += tempDeutSIMD;
  // cout << tempSIMDParticle.NDaughters()<< endl;
  // tempSIMDParticle.SetPDG(deuteron.GetPDG()>0? 103004 : -103004);
  ht_ndaughters = 3;

  if ((tempLambdaSIMD.GetDistanceFromParticle(tempDeutSIMD) < float_v(fDistanceCut)).isEmpty()) return false;


  float_v l,dl;
  float_m isParticleFromVertex;
  tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl, &isParticleFromVertex);
  ht_ldl = l[0]/dl[0];
  ht_l = l[0];
  ht_dl = dl[0];

  // deuteron.TransportToParticle(v_lambda);
  // KFParticle particle=v_lambda;
  // particle.SetPDG(deuteron.GetPDG()>0? 103004 : -103004);
  // particle+= deuteron;
  KFParticle particle;
  tempSIMDParticle.GetKFParticle( particle, 0);
  // cout <<particle.GetMass() << endl;
  ht_chi2ndf = particle.Chi2()/particle.NDF();
  if (particle.GetMass()>3.05 || particle.GetMass()<2.95) return false;

  ht_bdfvtx = tempSIMDParticle.GetDistanceFromVertex(pv)[0];
  KFParticleSIMD tempMotherTopo = tempSIMDParticle;
  tempMotherTopo.SetProductionVertex(pv);
  ht_bdfvtx2 = tempMotherTopo.GetDistanceFromVertex(pv)[0];
  ht_chi2topo = double(tempMotherTopo.Chi2()[0])/double(tempMotherTopo.NDF()[0]);
  ht_chi2 = double(tempMotherTopo.Chi2()[0]);
  ht_NDF = double(tempMotherTopo.NDF()[0]);

  if ((tempSIMDParticle.Chi2()/static_cast<float_v>(tempSIMDParticle.NDF()) < float_v(cut[2])).isEmpty()  ) return false; 
  // cout <<"pass chi2NDF cut" << endl;
  if (KFPMath::Finite(tempSIMDParticle.GetChi2()).isEmpty()) return false;
  // cout <<"pass chi2 finit cut" << endl;
  if ( (tempSIMDParticle.GetChi2() > 0.0f).isEmpty() ) return false;
  // cout <<"pass chi2 >0 cut" << endl;
  if ( (tempSIMDParticle.GetChi2() == tempSIMDParticle.GetChi2()).isEmpty() ) return false; // copy from KFParticleFinder function
  if (((tempMotherTopo.GetChi2()/float_v(tempMotherTopo.GetNDF()))<cut[1]).isEmpty() ) return false;
  // cout <<"pass masstopo NDF cut" << endl;
  if (!(ht_ldl>cut[0])) return false; // not primary
  if (!(ht_l<200.)) return false;
  if (!(ht_l>fLCut)) return false;
  if (isParticleFromVertex.isEmpty()) return false;
  // cout <<"pass cut" << endl;

  bparticleid   = particle.GetPDG();
  bparticlemass = particle.GetMass();
  bx = particle.GetX();
  by = particle.GetY();
  bz = particle.GetZ();
  bpx = particle.GetPx();
  bpy = particle.GetPy();
  bpz = particle.GetPz();
  bdl = particle.GetDecayLength();

  bpionpx = pion.GetPx();
  bpionpy = pion.GetPy();
  bpionpz = pion.GetPz();

  bprotonpx = proton.GetPx();
  bprotonpy = proton.GetPy();
  bprotonpz = proton.GetPz();


  const KFParticle* v_d01[2] = {&pion, &proton};
  const KFParticle* v_d02[2] = {&pion, &deuteron};
  const KFParticle* v_d12[2] = {&proton, &deuteron};

  KFParticle v_01;
  KFParticle v_02;
  KFParticle v_12;

  v_01.Construct(v_d01, 2, 0);
  v_01.GetMass(mass_01, mass_01_err);
  v_01_chi2primary = v_01.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
  v_01_pvdca = v_01.GetDistanceFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());

  v_02.Construct(v_d02, 2, 0);
  v_02.GetMass(mass_02, mass_02_err);
  v_02_chi2primary = v_02.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
  v_02_pvdca = v_02.GetDistanceFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());

  v_12.Construct(v_d12, 2, 0);
  v_12.GetMass(mass_12, mass_12_err);
  v_12_chi2primary = v_12.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
  v_12_pvdca = v_12.GetDistanceFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());

  v_01_chi2ndf = v_01.Chi2()/v_01.NDF();
  v_02_chi2ndf = v_02.Chi2()/v_02.NDF();
  v_12_chi2ndf = v_12.Chi2()/v_12.NDF();

  v_01_dca = pion.GetDistanceFromParticleXY(proton);
  v_02_dca = pion.GetDistanceFromParticleXY(deuteron);
  v_12_dca = proton.GetDistanceFromParticleXY(deuteron);

  v_012_dca = v_01.GetDistanceFromParticleXY(deuteron);

  //			cout<<"dca:"<<v_01_dca<<" "<<v_02_dca<<" "<<v_12_dca<<endl;
  //			cout<<"chi2ndf:"<<v_01_chi2ndf<<" "<<v_02_chi2ndf<<" "<<v_12_chi2ndf<<endl;

  dca_pion2 = pion.GetDistanceFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
  dca_proton2 = proton.GetDistanceFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
  dca_deuteron2 = deuteron.GetDistanceFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
  // chi2primary_d = deuteron.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
  // chi2primary_pi = pion.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
  // chi2primary_proton = proton.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());

  return true;
}
int StKFParticleAnalysisMaker::PassTrigger(StPicoDst* fPicoDst, int fsnn)
{
  int trigger = 0;
  if(fsnn==3){
    if(fPicoDst->event()->isTrigger(620052)) trigger += 1;
  }else if(fsnn==7){
    if(fPicoDst->event()->isTrigger(630052)) trigger += 1;
  }else{
    trigger += 1;
  }
  return trigger;
}
int StKFParticleAnalysisMaker::CheckIfBadRun(int fsnn, int brunid)
{
  int notbadrun = 0;
  if(fsnn==3){
    //benjamin https://drupal.star.bnl.gov/STAR/system/files/Kimelman_3GeV_run_by_run_QA_badRuns.pdf

    if(brunid==19151029) notbadrun+=1;
    if(brunid==19151045) notbadrun+=1;
    if(brunid==19152001) notbadrun+=1;
    if(brunid==19152078) notbadrun+=1;
    if(brunid==19153023) notbadrun+=1;
    if(brunid==19153032) notbadrun+=1;
    if(brunid==19153065) notbadrun+=1;
    if(brunid==19154012) notbadrun+=1;
    if(brunid==19154013) notbadrun+=1;
    if(brunid==19154014) notbadrun+=1;
    if(brunid==19154015) notbadrun+=1;
    if(brunid==19154016) notbadrun+=1;
    if(brunid==19154017) notbadrun+=1;
    if(brunid==19154018) notbadrun+=1;
    if(brunid==19154019) notbadrun+=1;
    if(brunid==19154020) notbadrun+=1;
    if(brunid==19154021) notbadrun+=1;
    if(brunid==19154022) notbadrun+=1;
    if(brunid==19154023) notbadrun+=1;
    if(brunid==19154024) notbadrun+=1;
    if(brunid==19154026) notbadrun+=1;
    if(brunid==19154046) notbadrun+=1;
    if(brunid==19154051) notbadrun+=1;
    if(brunid==19154056) notbadrun+=1;

  }
  return notbadrun;
}
int StKFParticleAnalysisMaker::getVtxBin(float vx, float vy, float vz)
{
  return 0;
  // float vx_low = -0.5;
  // float vx_high = 0.5;
  //
  // float vy_low = -2.5;
  // float vy_high = -1.5;
  //
  // float vz_low = 200.;
  // float vz_high = 201.5;
  //
  // const int binX = 3, binY=3, binZ = 3;
  // float vx_ed[binX+1]={ -0.5,-0.2,0.2,0.5}; 
  // float vy_ed[binY+1]={ -2.5,-2.2,-1.8,-1.5}; 
  // float vz_ed[binZ+1]={ 200.,200.5,201,201.5}; 
  //
  // if (vz>=vz_high || vz<vz_low ) return -1;
  // if (vy>=vy_high || vy<vy_low ) return -1;
  // if (vx>=vx_high || vx<vx_low ) return -1;
  //
  // // int ibX = std::floor((vx-vx_low)/((vx_high-vx_low)/binX));
  // // int ibY = std::floor((vy-vy_low)/((vy_high-vy_low)/binY));
  // // int ibZ = std::floor((vz-vz_low)/((vz_high-vz_low)/binZ));
  // int ibX=-1, ibY=-1,ibZ=-1;
  // for (ibX=0;ibX<binX;ibX++) 
  // {
  //   if (vx_ed[ibX] <= vx && vx_ed[ibX+1]> vx) break;
  // } 
  // for (ibY=0;ibY<binY;ibY++) 
  // {
  //   if (vy_ed[ibY] <= vy && vy_ed[ibY+1]> vy) break;
  // } 
  //
  // for (ibZ=0;ibZ<binZ;ibZ++) 
  // {
  //   if (vz_ed[ibZ] <= vz && vz_ed[ibZ+1]> vz) break;
  // }
  //
  // return ibX+ibY*binX+ibZ*binX*binY;
}
