//*-- Author : Yuri Fisyak 02/02/2016
#include "StKFParticleAnalysisMaker.h"
#include "TDirectory.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TChain.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TSystem.h"
#include <TLorentzVector.h>
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
#include "StRefMultCorr/StRefMultCorr.h"
#include "StRefMultCorr/CentralityMaker.h"
ClassImp(StKFParticleAnalysisMaker);

const Float_t Mass[3] = { //Gev
  0.139570, // pion
  0.493677, // kaon
  0.938272  // proton
};

//________________________________________________________________________________
StKFParticleAnalysisMaker::StKFParticleAnalysisMaker(const char *name) : StMaker(name), fNTrackTMVACuts(0), fIsPicoAnalysis(true), fsnn(-999), fdEdXMode(1), 
fStoreTmvaNTuples(false), fProcessSignal(false), fCollectTrackHistograms(false), fCollectPIDHistograms(false),fTMVAselection(false), fStoremctree(false),
fFlowAnalysis(false), fFlowChain(NULL), fFlowRunId(-1), fFlowEventId(-1), fCentrality(-1), fFlowFiles(), fFlowMap(),
fRunCentralityAnalysis(0), fRefmultCorrUtil(0), fCentralityFile(""), fAnalyseDsPhiPi(false)
{
  memset(mBeg,0,mEnd-mBeg+1);
  
  fNTuplePDG[0] = 3122;
  fNTuplePDG[1] = 3334;
  fNTuplePDG[2] = 3312;
  
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
  
  f = GetTFile();
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
  
  fRefmultCorrUtil = CentralityMaker::instance()->getRefMultCorr() ;
  
  //Initialise the chain with files containing centrality and reaction plane
  if(fFlowAnalysis)
  {
    std::cout << "StKFParticleAnalysisMaker: run flow analysis. Flow file list:"<<std::endl;
    
    fFlowChain = new TChain("psi_tree");
    for(unsigned int iFlowFile=0; iFlowFile<fFlowFiles.size(); iFlowFile++)
    {
      std::cout << "      " << fFlowFiles[iFlowFile] << std::endl;
      fFlowChain->Add(fFlowFiles[iFlowFile].Data());
    }
    
    fFlowChain->SetBranchStatus("*",0);
    fFlowChain->SetBranchAddress("isPileUp", &isPileUp);
    fFlowChain->SetBranchStatus("isPileUp", 1);
    fFlowChain->SetBranchAddress("runnumber", &fFlowRunId);
    fFlowChain->SetBranchStatus("runnumber", 1);
    fFlowChain->SetBranchAddress("eventid", &fFlowEventId);
    fFlowChain->SetBranchStatus("eventid", 1);
    fFlowChain->SetBranchAddress("centnumber", &fCentrality);
    fFlowChain->SetBranchStatus("centnumber", 1);
    fFlowChain->SetBranchAddress("refMultPrim", &refMultPrim);
    fFlowChain->SetBranchStatus("refMultPrim", 1);
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
  
  //  string sRun[] = { "11", "14", "16", "16_sst", "16_nosst", "00", "00" };
  //  mMyFile = new TFile( Form("kfptree_run%s_%d_%d.root",sRun[fRunId].c_str(),mRunNumber,mRunSubId), "recreate","tree and histo",6);
  mMyFile = new TFile("ana_tree.root","RECREATE");
  mMyFile->cd();
  ptree = new TTree( "tree", "tree for strangeness with KFParticle" );
  pEve  = new EventClass(ptree);
  pEve->InitWrite();
//  SetHistograms();
  SetBadRun();
  
  return kStOK;
}
//________________________________________________________________________________
Int_t StKFParticleAnalysisMaker::InitRun(Int_t runumber) 
{
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
  bool bDEBUG = false;
  bool isGoodEvent = false;
  
  StPicoEvent *picoEvent = NULL;
  
  if(fIsPicoAnalysis)
  {
    fPicoDst = StPicoDst::instance();
    if(!fPicoDst) return kStOK;
    picoEvent = (StPicoEvent*)fPicoDst->event();
    if( !picoEvent ) return kStWarn;
  }
  else
  {
    fMuDst = StMuDst::instance();
    if(!fMuDst) return kStOK;
    else { if(StMuDst::instance()->numberOfPrimaryVertices() == 0 ) return kStOK; }
  }
  // ==================================================
  
  if(fIsPicoAnalysis){
    beventid = fPicoDst->event()->eventId();
    brunid = fPicoDst->event()->runId();
    
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
    brunid = fMuDst->event()->runId();
    
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
//  hvtx->Fill(bVz);
  // ==================================================
  int trigger=0;
  if(fIsPicoAnalysis){
    if(fsnn==27){
      if(fPicoDst->event()->isTrigger(610001)) trigger += 1;
      if(fPicoDst->event()->isTrigger(610011)) trigger += 10;
      if(fPicoDst->event()->isTrigger(610021)) trigger += 100;
      if(fPicoDst->event()->isTrigger(610031)) trigger += 1000;
      if(fPicoDst->event()->isTrigger(610041)) trigger += 10000;
      if(fPicoDst->event()->isTrigger(610051)) trigger += 100000;
    }else if(fsnn==3){
      if(fPicoDst->event()->isTrigger(620052)) trigger += 1;
    }else{
      trigger += 1;
    }
  }
  else{
    trigger+=1; //auto pass trigger because we mudst already has trigger seltcion
  }
  
  if(IsBadRun(brunid)) return kStOK;
  //cout<<"trigger:"<<trigger<<endl;
  if(trigger==0) return kStOK;
  // ==================================================
  // additional check
  //cut on refmult and tofmult, trigger, to cut on good event. this cut is only for real analysis!!
  //when fStoreTmvaNTuples is true, there is no need to cut on these variables
  // ==================================================
  int crefmult;
  int ctofmult;
  if(fIsPicoAnalysis)
  {
    crefmult = fPicoDst->event()->refMult();
    ctofmult = fPicoDst->event()->btofTrayMultiplicity();
  }
  else
  {
    crefmult = fMuDst->event()->refMult();
    ctofmult = fMuDst->event()->btofTrayMultiplicity();
  }
  //if(crefmult>(0.42*ctofmult+16) || crefmult<(0.02*ctofmult-5)) isGoodEvent = false;
//  if( fIsPicoAnalysis ){
//    hNBTofvsRM  ->Fill( picoEvent->refMult(),  picoEvent->nBTOFMatch() );
//    hNBTofvsGRM ->Fill( picoEvent->grefMult(), picoEvent->nBTOFMatch() );
//    hTofMvsRM   ->Fill( picoEvent->refMult(),  picoEvent->btofTrayMultiplicity() );
//    hTofMvsGRM  ->Fill( picoEvent->grefMult(), picoEvent->btofTrayMultiplicity() );
//    hRankvsGTrk ->Fill( picoEvent->numberOfGlobalTracks(),          picoEvent->ranking() );
//    hRankvsPTrk ->Fill( picoEvent->refMult()+picoEvent->refMult2(), picoEvent->ranking() );
//    hRankvsNBTof->Fill( picoEvent->nBTOFMatch(),                    picoEvent->ranking() );
//    hGTrkvsPTrk ->Fill( picoEvent->refMult()+picoEvent->refMult2(), picoEvent->numberOfGlobalTracks() );
//    hGTrkvsRM   ->Fill( picoEvent->refMult(), picoEvent->numberOfGlobalTracks() );
//    hGRMvsRM    ->Fill( picoEvent->refMult(), picoEvent->grefMult() );
//
//    hNVpdvsRM   ->Fill( picoEvent->refMult(),              picoEvent->nVpdHitsEast()+picoEvent->nVpdHitsWest() );
//    hNVpdvsGRM  ->Fill( picoEvent->grefMult(),             picoEvent->nVpdHitsEast()+picoEvent->nVpdHitsWest() );
//    hNVpdvsGTrk ->Fill( picoEvent->numberOfGlobalTracks(), picoEvent->nVpdHitsEast()+picoEvent->nVpdHitsWest() );
//    hNBEMCvsRM  ->Fill( picoEvent->refMult(),              picoEvent->nBEMCMatch() );
//    hNBEMCvsGRM ->Fill( picoEvent->grefMult(),             picoEvent->nBEMCMatch() );
//    hNBEMCvsGTrk->Fill( picoEvent->numberOfGlobalTracks(), picoEvent->nBEMCMatch() );
//  }
  // ==================================================
  //work out the primary track count (for FXT)
  // ==================================================
  if(fIsPicoAnalysis)
  {
    countrefmult = 0;
    for(unsigned int iTrack = 0; iTrack < fPicoDst->numberOfTracks(); iTrack++)
    {
      const StPicoTrack *ptrk = (StPicoTrack*)fPicoDst->track(iTrack);
      if(! ptrk) continue;
      if(!ptrk->isPrimary())  continue;  // now selecting primary tracks
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
  //end work out the primary track count (for FXT)
  // ==================================================
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
  //events removed at this level include no track events and
  //whatever is in Process event skip level
  // ==================================================
  vector<int> triggeredTracks;
  if(maxGBTrackIndex > 0)
    fStKFParticleInterface->ResizeTrackPidVectors(maxGBTrackIndex+1);
  
  if(fIsPicoAnalysis)
    isGoodEvent = fStKFParticleInterface->ProcessEvent(fPicoDst, triggeredTracks);
  else//b
    //no embedding
    //isGoodEvent = fStKFParticleInterface->ProcessEvent(fMuDst, triggeredTracks, fProcessSignal);
    //embedding
    isGoodEvent = fStKFParticleInterface->ProcessEvent(fMuDst, mcTracks, mcIndices, fProcessSignal);
  //===================================================================================================
  cent9 = 0;
  reweight = 0;
  refmultcor = 0;
  int centralityBin = -1;

  if(isGoodEvent)
  {
//    hvtxgood->Fill(bVz);
//    hrefmult->Fill(crefmult);
    
    float centralityWeight = 0.;
    if(fRunCentralityAnalysis)
    {
      if(fIsPicoAnalysis){
        
        fRefmultCorrUtil->init(fPicoDst->event()->runId());
        if(! (fRefmultCorrUtil->isBadRun(fPicoDst->event()->runId())) )
        {
          fRefmultCorrUtil->initEvent(fPicoDst->event()->refMult(), fPicoDst->event()->primaryVertex().z(), fPicoDst->event()->ZDCx()) ;
          centralityBin = fRefmultCorrUtil->getCentralityBin9();
          centralityWeight = fRefmultCorrUtil->getWeight();
          
          cent9 = centralityBin;
          reweight = centralityWeight;
        }else{
          isGoodEvent = false;
        }
        
        Bool_t isPileUpEvt = !fRefmultCorrUtil->passnTofMatchRefmultCut(1.*fPicoDst->event()->refMult(), 1.*fPicoDst->event()->nBTOFMatch()); //reject pileup events
        
        if(isPileUpEvt) isGoodEvent = false;
        
        refmultcor = fRefmultCorrUtil->getRefMultCorr();
      }else{
        
        fRefmultCorrUtil->init(fMuDst->event()->runId());
        if(! (fRefmultCorrUtil->isBadRun(fMuDst->event()->runId())) )
        {
          fRefmultCorrUtil->initEvent(fMuDst->event()->refMult(), bVz, fMuDst->event()->runInfo().zdcCoincidenceRate()) ;
          //cout<<"zdcrate:"<<fMuDst->event()->runInfo().zdcCoincidenceRate()<<endl;
          
          centralityBin = fRefmultCorrUtil->getCentralityBin9();
          centralityWeight = fRefmultCorrUtil->getWeight();
          cent9 = centralityBin;
          reweight = centralityWeight;
          //skipt rejecting pileup for mudst
          refmultcor = fRefmultCorrUtil->getRefMultCorr();
        }else{
          isGoodEvent = false;
        }
      }
      
//      if(crefmult>60){
//        wrefmult->Fill(refmultcor);
//      }else{
//        wrefmult->Fill(refmultcor, reweight);
//      }
    }
    if(!isGoodEvent) return kStOK;
    //================================================================================
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
          }
        }
      }
    }
    //================================================================================
    if(fFlowAnalysis){
      long entryId = GetUniqueEventId(brunid, beventid);
      //cout<<"entryId:"<<entryId<<" "<<brunid<<" "<<beventid<<endl;
      std::map<long,int>::iterator flowMapIterator = fFlowMap.find(entryId);
      if (flowMapIterator != fFlowMap.end())
      {
        fFlowChain->GetEvent(fFlowMap[GetUniqueEventId(brunid, beventid)]);
        centralityBin = fCentrality;
//        hrefmult_wPileup->Fill(countrefmult);
//        hrefmult_woPileup->Fill(refMultPrim);
//        hCentrality->Fill(fCentrality);
//        hCentrality_weight->Fill(fCentrality, gweight);
        //cout<<"fFlowRunId:"<<fFlowRunId<<" "<< fFlowEventId<< " "<< fCentrality<<" "<< endl;
        //cout<<"flow:"<<psi_1_EPD_0 <<" "<< psi_1_EPD_1<<" "<<psi_1_EPD_2<<" "<<psi_1_EPD_3<<endl;
      }
    }
    //================================================================================
    centralityWeight = 1;
    fStKFParticlePerformanceInterface->SetMCTracks(mcTracks);
    fStKFParticlePerformanceInterface->SetMCIndexes(mcIndices);
    fStKFParticlePerformanceInterface->SetCentralityBin(centralityBin);
    fStKFParticlePerformanceInterface->SetCentralityWeight(centralityWeight);
    Int_t nevent = 100000;
    fStKFParticlePerformanceInterface->SetPrintEffFrequency(nevent);
    fStKFParticlePerformanceInterface->PerformanceAnalysis();
    //================================================================================
    if(fStoreTmvaNTuples)
    {
      for(int iParticle=0; iParticle<fStKFParticlePerformanceInterface->GetNReconstructedParticles(); iParticle++)
      {
        KFParticle particle;
        bool isMCParticle = fStKFParticlePerformanceInterface->GetParticle(particle, iParticle);
        
        //       cout<<"isMCParticle:"<<isMCParticle<< " "<< particle.GetPDG()<<endl;
        
        //TOBEREVERTED //REVERTED
        if( !( (fProcessSignal && isMCParticle) || (!fProcessSignal && !isMCParticle) ) ) continue;
        
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
              //sideband = (side_mass>1.61 && side_mass<1.66) || (side_mass>1.685 && side_mass<1.85);//for testing
              sideband = (side_mass>1.61 && side_mass<1.66) || (side_mass>1.685 && side_mass<1.73);//for consistency with long code
              
              //TOBEREVERTED
              //sideband = true;//for testing
              sideband = side_mass<1.78;
              //cout<<"fProcessSignal:"<<fProcessSignal<<" isMCParticle:"<<isMCParticle<<endl;
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
  //===================================================================================================
  // track loop
  if(isGoodEvent){
    ftrkIdV0A.clear();
    ftrkIdV0B.clear();
    // Track loop for KFParticle
    unsigned int cur_ntrk = 0;
    for( int iprt=0; iprt<fStKFParticlePerformanceInterface->GetNReconstructedParticles(); iprt++ ){
      
      //      const KFParticle particle = fStKFParticleInterface->GetParticles()[iprt];
      KFParticle particle;
      bool isMCParticle = fStKFParticlePerformanceInterface->GetParticle(particle, iprt);
      //      bismc =0;
      //      if(isMCParticle) bismc=1;
      
      // KS0       =  310
      // Lambda    =  3122
      // LambdaBar = -3122
      // Xi        =  3312
      // Xibar     = -3312
      // Omega     =  3334
      // OmegaBar  = -3334
      // See KFParticle/KFParticleFinder.cxx
      int pdg = particle.GetPDG();
      
      int pid = -1;
      if     ( abs(pdg)==3122 ) pid = 0;
      else if( abs(pdg)==3312 ) pid = 2;
      else if( abs(pdg)==3334 ) pid = 4;
      else if( pdg==310 )       pid = 6;
      if( pid<0 ) continue;
      if( pdg<0 ) pid += 1;
      
      if(abs(pdg)!=3312) continue; // only for Xi analysis
      
      float pt   = particle.GetPt();
      float phi  = particle.GetPhi();
      float eta  = particle.GetEta();
      float rap  = particle.GetRapidity();
      float mass = particle.GetMass();
      float dL   = particle.GetDecayLength();
      float dca  = fStKFParticleInterface->Getdca(particle.DaughterIds()[0]);
      if( mass<1.0e-05 || fabs(phi)<1.0e-05 || fabs(eta)<1.0e-05 ) continue;
      //if( fabs(phi)<1.0e-04 ) cout <<"parent-particle pdg="<< pdg <<" Ndaughter="<< particle.NDaughters() <<" mass="<< mass <<" phi="<< phi << endl;
      
      // mass cut for tree
      bool isMassOK = true;
      if( abs(pdg)==3122 && mass>1.5 ) isMassOK = false;
      if( abs(pdg)==3312 && mass>1.7 ) isMassOK = false;
      if( abs(pdg)==3334 && mass>2.1 ) isMassOK = false;
      if( pdg==310 && mass>0.9 )       isMassOK = false;
      if( !isMassOK ) continue;
      
      // short life time partilce
      float_v l,dl;
      float SIMD_mass, SIMD_pt, SIMD_phi, SIMD_eta, SIMD_rapidity, SIMD_dca, SIMD_dca2D, SIMD_decaylength, SIMD_lifetime, SIMD_pathlength, SIMD_ldl, SIMD_l, SIMD_dl, SIMD_chi2topo, SIMD_chi2ndf, SIMD_pdg;
      if(pid <6) {
        KFParticleSIMD tempSIMDParticle(particle);
        KFParticleSIMD pv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
        tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl);
        SIMD_ldl = l[0]/dl[0];
        SIMD_l = l[0];
        SIMD_dl = dl[0];
        tempSIMDParticle.SetProductionVertex(pv);
        SIMD_chi2topo = double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);
        SIMD_chi2ndf = particle.Chi2()/particle.NDF();
        SIMD_dca = tempSIMDParticle.GetDistanceFromVertex(pv)[0]; //SIMD_dca
        SIMD_decaylength = tempSIMDParticle.GetDecayLength()[0];
        
        if(bDEBUG) {
          if(SIMD_dca > 1.0)    continue;  //xi DCA to primary vertex
          if(SIMD_decaylength < 4.0)   continue; //xi decay length
        }
        
        SIMD_lifetime = tempSIMDParticle.GetLifeTime()[0];
        SIMD_mass = tempSIMDParticle.GetMass()[0];
        SIMD_pt = tempSIMDParticle.GetPt()[0];
        SIMD_pdg = tempSIMDParticle.GetPDG()[0];
        SIMD_phi = tempSIMDParticle.GetPhi()[0];
        SIMD_eta = tempSIMDParticle.GetEta()[0];
        SIMD_rapidity = tempSIMDParticle.GetRapidity()[0];
        SIMD_dca2D = tempSIMDParticle.GetDistanceFromVertexXY(pv)[0];
        SIMD_pathlength = tempSIMDParticle.GetR()[0];  //* distance to the origin
        
        if(bDEBUG) {
          if( SIMD_pt > 5. )    continue;
        }
        
        if(bDEBUG) {
          cout << "mother SIMD_pdg: " << SIMD_pdg << endl;
          cout << "mother SIMD_mass: " << SIMD_mass << ";    SIMD_dca = " << SIMD_dca << "; SIMD_decaylength = " << SIMD_decaylength << ";   SIMD_pt = " << SIMD_pt << ";    SIMD_rapidity = " << SIMD_rapidity << "; SIMD_phi = " << SIMD_phi << endl;
          cout << "Found " << particle.NDaughters() << " daughter particles: " << endl;
        }
      }
      
      // Get daughter information
      float SIMD_A_mass, SIMD_A_pt, SIMD_A_phi, SIMD_A_eta, SIMD_A_rapidity, SIMD_A_dca, SIMD_A_dca2D, SIMD_A_decaylength,  SIMD_A_lifetime, SIMD_A_pathlength, SIMD_A_ldl, SIMD_A_l, SIMD_A_dl, SIMD_A_chi2topo, SIMD_A_chi2ndf, SIMD_A_pdg;
      float chi2primary_A, dL_A, dca_A, nhits_A, rapidity_A, dedx_A, m2_A, pdg_A;
      float chi2primary_B, dca_B, nhits_B, rapidity_B, dedx_B, m2_B, pdg_B;
      float chi2primary_AB, dca_AB, nhits_AB, rapidity_AB, dedx_AB, m2_AB, pdg_AB;
      float chi2primary_AA, dca_AA, nhits_AA, rapidity_AA, dedx_AA, m2_AA, pdg_AA;
      float dgdca, dgdca_db;
      
      float nsigmapion_A, nsigmaproton_A, nsigmaKaon_A;
      float nsigmapion_B, nsigmaproton_B, nsigmaKaon_B;
      float nsigmapion_AA, nsigmaproton_AA, nsigmaKaon_AA;
      float nsigmapion_AB, nsigmaproton_AB, nsigmaKaon_AB;
      
      float massA  = -1.0;
      float massB  = -1.0;
      float massAA = -1.0;
      TVector3 pA( 0., 0., 0. );
      TVector3 pB( 0., 0., 0. );
      TVector3 pAA( 0., 0., 0. );
      TVector3 pAB( 0., 0., 0. );
      int trkIdA = -1;
      int trkIdB = -1;
      
      int targetDg;
      if     ( pid<2 ) targetDg = 2212; // p or pbar
      else if( pid<4 ) targetDg = 3122; // Lambda or LambdaBar
      else if( pid<6 ) targetDg = 3122; // Lambda or LambdaBar
      else             targetDg = 211;  // pion
      
      int Ndaughter = particle.NDaughters();
      if( Ndaughter>2 ){
        cout <<"skip this particle, because of NDaughters = "<< Ndaughter << endl;
        continue;
      }
      
      for( int idg=0;  idg<Ndaughter; idg++ ){
        const int daughterId = particle.DaughterIds()[idg];
        if( daughterId<0 ) continue;
        
        KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterId];
        int dg_pdg = daughter.GetPDG(); // in case of daughter Lambda from Xi or ..., -1 will be returned
        float dg_mass, dg_massE;
        daughter.GetMass(dg_mass,dg_massE);
        
        // select daughter baryons
        if( pid<6 ){
          if( abs(dg_pdg)==targetDg || dg_mass>0.6 ){
            // pi(<- lambda), or lambda
            massA  = dg_mass;
            pA     = TVector3( daughter.GetPx(), daughter.GetPy(), daughter.GetPz() );
            trkIdA = daughter.Id();
            
            // pi(<- lambda)
            chi2primary_A = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
            dca_A = fStKFParticleInterface->Getdca(daughter.DaughterIds()[0]);
            nhits_A = fStKFParticleInterface->GetNHits(daughter.DaughterIds()[0]);
            dL_A = daughter.GetDecayLength();
            //cout << "Ld pion or Ld: " << endl;
            //cout << "daughter.DaughterIds()[0] = " << daughter.DaughterIds()[0] << ", " << daughter.DaughterIds()[1] <<   ";   daughter.GetId() = " << daughter.GetId() << endl;
            //cout << "Getdca(daughter.DaughterIds()[0]) = " << dca_A << ";   Getdca(daughter.GetId()) = " <<           fStKFParticleInterface->Getdca(daughter.GetId()) << endl;
            
            pdg_A = daughter.GetPDG();
            rapidity_A = daughter.GetRapidity();
            nsigmapion_A   = fStKFParticleInterface->GetdEdXNSigmaPion(daughter.DaughterIds()[0]);
            nsigmaproton_A = fStKFParticleInterface->GetdEdXNSigmaProton(daughter.DaughterIds()[0]);
            nsigmaKaon_A   = fStKFParticleInterface->GetdEdXNSigmaKaon(daughter.DaughterIds()[0]);
            dedx_A = fStKFParticleInterface->GetdEdX(daughter.DaughterIds()[0]);
            m2_A = fStKFParticleInterface->GetM2tof(daughter.DaughterIds()[0]);
            if( abs(pdg)!=3122 ){ // except Lambda
              // dA info
              KFParticleSIMD ttempSIMDParticle(daughter);
              KFParticleSIMD tpv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
              ttempSIMDParticle.GetDistanceToVertexLine(tpv, l, dl);
              SIMD_A_ldl = l[0]/dl[0];
              SIMD_A_l = l[0];
              SIMD_A_dl = dl[0];
              
              ttempSIMDParticle.SetProductionVertex(tpv);
              SIMD_A_chi2topo = double(ttempSIMDParticle.Chi2()[0])/double(ttempSIMDParticle.NDF()[0]);
              SIMD_A_chi2ndf = daughter.Chi2()/daughter.NDF();
              
              SIMD_A_dca = ttempSIMDParticle.GetDistanceFromVertex(tpv)[0]; //SIMD_A_dca
              SIMD_A_decaylength = ttempSIMDParticle.GetDecayLength()[0];
              
              if(bDEBUG) {
                if(SIMD_A_dca > 5.0 || SIMD_A_dca < 0.2)    continue;  //Lambda DCA to primary vertex
                if(SIMD_A_decaylength < 3.0)   continue; //Lambda decay length
              }
              
              SIMD_A_lifetime = ttempSIMDParticle.GetLifeTime()[0];
              SIMD_A_pdg = ttempSIMDParticle.GetPDG()[0];
              SIMD_A_mass = ttempSIMDParticle.GetMass()[0];
              SIMD_A_pt = ttempSIMDParticle.GetPt()[0];
              SIMD_A_phi = ttempSIMDParticle.GetPhi()[0];
              SIMD_A_eta = ttempSIMDParticle.GetEta()[0];
              SIMD_A_rapidity = ttempSIMDParticle.GetRapidity()[0];
              SIMD_A_dca2D = ttempSIMDParticle.GetDistanceFromVertexXY(tpv)[0];
              SIMD_A_pathlength = ttempSIMDParticle.GetR()[0];  //* distance to the origin
              
              // Ld daughters
              int Ngrdaughter = daughter.NDaughters();
              if( Ngrdaughter>2 ){
                cout <<"skip this particle, because of NgrDaughters = "<< Ngrdaughter << endl;
                continue;
              }
              
              for( int kdg=0; kdg<Ngrdaughter; kdg++ ){
                const int grdaughterId = daughter.DaughterIds()[kdg];
                if( grdaughterId<0 ) continue;
                
                KFParticle grdaughter = fStKFParticleInterface->GetParticles()[grdaughterId];
                float grdg_mass, grdg_massE;
                grdaughter.GetMass(grdg_mass,grdg_massE);
                
                if( abs(grdaughter.GetPDG())!=2212 && abs(grdaughter.GetPDG())!=211 ) cout <<"idg="<< idg <<" kdg="<< kdg <<" p-pdg="<< pdg <<" dg-pdg="<< dg_pdg <<" grdg_pdg="<< grdaughter.GetPDG() <<" p_mass="<< mass <<" dg_mass="<< dg_mass <<" grdg_mass="<< grdg_mass <<" dg_phi="<< daughter.GetPhi() <<" grdg_pt="<< grdaughter.GetPt() <<" grdg_phi="<< grdaughter.GetPhi() << endl;
                
                if( abs(grdaughter.GetPDG())==2212 ){ // select only protons
                  massAA = grdg_mass;
                  pAA = TVector3( grdaughter.GetPx(), grdaughter.GetPy(), grdaughter.GetPz() );
                  
                  chi2primary_AA = grdaughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
                  dca_AA = fStKFParticleInterface->Getdca(grdaughter.DaughterIds()[0]);
                  nhits_AA = fStKFParticleInterface->GetNHits(grdaughter.DaughterIds()[0]);
                  pdg_AA = grdaughter.GetPDG();
                  rapidity_AA = grdaughter.GetRapidity();
                  nsigmapion_AA   = fStKFParticleInterface->GetdEdXNSigmaPion(grdaughter.DaughterIds()[0]);
                  nsigmaproton_AA = fStKFParticleInterface->GetdEdXNSigmaProton(grdaughter.DaughterIds()[0]);
                  nsigmaKaon_AA   = fStKFParticleInterface->GetdEdXNSigmaKaon(grdaughter.DaughterIds()[0]);
                  dedx_AA = fStKFParticleInterface->GetdEdX(grdaughter.DaughterIds()[0]);
                  m2_AA = fStKFParticleInterface->GetM2tof(grdaughter.DaughterIds()[0]);
                }else{
                  pAB = TVector3( grdaughter.GetPx(), grdaughter.GetPy(), grdaughter.GetPz() );
                  
                  chi2primary_AB = grdaughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
                  dca_AB = fStKFParticleInterface->Getdca(grdaughter.DaughterIds()[0]);
                  nhits_AB = fStKFParticleInterface->GetNHits(grdaughter.DaughterIds()[0]);
                  pdg_AB = grdaughter.GetPDG();
                  rapidity_AB = grdaughter.GetRapidity();
                  nsigmapion_AB   = fStKFParticleInterface->GetdEdXNSigmaPion(grdaughter.DaughterIds()[0]);
                  nsigmaproton_AB = fStKFParticleInterface->GetdEdXNSigmaProton(grdaughter.DaughterIds()[0]);
                  nsigmaKaon_AB   = fStKFParticleInterface->GetdEdXNSigmaKaon(grdaughter.DaughterIds()[0]);
                  dedx_AB = fStKFParticleInterface->GetdEdX(grdaughter.DaughterIds()[0]);
                  m2_AB = fStKFParticleInterface->GetM2tof(grdaughter.DaughterIds()[0]);
                }
              }
              
              //dca between grandaughters
              KFParticle dgA = fStKFParticleInterface->GetParticles()[daughter.DaughterIds()[0]];
              KFParticle dgB = fStKFParticleInterface->GetParticles()[daughter.DaughterIds()[1]];
              dgdca_db = fabs(dgA.GetDistanceFromParticle(dgB));
              
            }// for grand-daughter
          }else{
            // proton or bach(pi, k)
            chi2primary_B = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
            dca_B = fStKFParticleInterface->Getdca(daughter.DaughterIds()[0]);
            nhits_B = fStKFParticleInterface->GetNHits(daughter.DaughterIds()[0]);
            pdg_B = daughter.GetPDG();
            rapidity_B = daughter.GetRapidity();
            nsigmapion_B   = fStKFParticleInterface->GetdEdXNSigmaPion(daughter.DaughterIds()[0]);
            nsigmaproton_B = fStKFParticleInterface->GetdEdXNSigmaProton(daughter.DaughterIds()[0]);
            nsigmaKaon_B   = fStKFParticleInterface->GetdEdXNSigmaKaon(daughter.DaughterIds()[0]);
            dedx_B = fStKFParticleInterface->GetdEdX(daughter.DaughterIds()[0]);
            m2_B = fStKFParticleInterface->GetM2tof(daughter.DaughterIds()[0]);
            
            massB  = dg_mass;
            pB     = TVector3( daughter.GetPx(), daughter.GetPy(), daughter.GetPz() );
            trkIdB = daughter.Id();
          }
        }else{ // if pid>=6
          if( dg_pdg==211 ){
            chi2primary_A = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
            dca_A = fStKFParticleInterface->Getdca(daughter.DaughterIds()[0]);
            nhits_A = fStKFParticleInterface->GetNHits(daughter.DaughterIds()[0]);
            pdg_A = daughter.GetPDG();
            rapidity_A = daughter.GetRapidity();
            nsigmapion_A   = fStKFParticleInterface->GetdEdXNSigmaPion(daughter.DaughterIds()[0]);
            nsigmaproton_A = fStKFParticleInterface->GetdEdXNSigmaProton(daughter.DaughterIds()[0]);
            nsigmaKaon_A   = fStKFParticleInterface->GetdEdXNSigmaKaon(daughter.DaughterIds()[0]);
            dedx_A = fStKFParticleInterface->GetdEdX(daughter.DaughterIds()[0]);
            m2_A = fStKFParticleInterface->GetM2tof(daughter.DaughterIds()[0]);
            
            massA  = dg_mass;
            pA     = TVector3( daughter.GetPx(), daughter.GetPy(), daughter.GetPz() );
            trkIdA = daughter.Id();
          }else{
            chi2primary_B = daughter.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
            dca_B = fStKFParticleInterface->Getdca(daughter.DaughterIds()[0]);
            nhits_B = fStKFParticleInterface->GetNHits(daughter.DaughterIds()[0]);
            pdg_B = daughter.GetPDG();
            rapidity_B = daughter.GetRapidity();
            nsigmapion_B   = fStKFParticleInterface->GetdEdXNSigmaPion(daughter.DaughterIds()[0]);
            nsigmaproton_B = fStKFParticleInterface->GetdEdXNSigmaProton(daughter.DaughterIds()[0]);
            nsigmaKaon_B   = fStKFParticleInterface->GetdEdXNSigmaKaon(daughter.DaughterIds()[0]);
            dedx_B = fStKFParticleInterface->GetdEdX(daughter.DaughterIds()[0]);
            m2_B = fStKFParticleInterface->GetM2tof(daughter.DaughterIds()[0]);
            
            massB  = dg_mass;
            pB     = TVector3( daughter.GetPx(), daughter.GetPy(), daughter.GetPz() );
            trkIdB = daughter.Id();
          }
        }
      }// idg
      
      KFParticle dA = fStKFParticleInterface->GetParticles()[particle.DaughterIds()[0]];
      KFParticle dB = fStKFParticleInterface->GetParticles()[particle.DaughterIds()[1]];
      dgdca = fabs(dA.GetDistanceFromParticle(dB));
      
      // only for Xi
      float mass_fake = -999.0;
      if( abs(pdg)==3312 ){
        TVector3 pV0_fake = pB + pAA;
        float eA_fake = sqrt( pAA.Mag2() + Mass[2]*Mass[2] );
        float eB_fake = sqrt( pB.Mag2() + Mass[0]*Mass[0] );
        float eV0_fake = eA_fake + eB_fake;
        mass_fake = sqrt( eV0_fake*eV0_fake - pV0_fake.Mag2() );
      }
      
      //float pdgMass = KFParticleDatabase::Instance()->GetMass(pdg); // does not work for Lambda (3112), pion mass is set...
      //float pdgE = sqrt( pdgMass*pdgMass + particle.GetP()*particle.GetP() );
      float obsE = sqrt( mass*mass + particle.GetP()*particle.GetP() );
      //cout <<" pdgE="<< pdgE <<" obsE="<< obsE <<" kfpE="<< particle.E() <<" mass="<< mass <<" pdgM="<< pdgMass <<" pdg="<< pdg << endl;
      
      // daughter baryon in the parent rest frame
      float phiA_pRF = -999.0;
      float theA_pRF = -999.0;
      float phiAA_dgRF = -999.0;
      float theAA_dgRF = -999.0;
      if( pid<6 ){
        TLorentzVector mom4d( TVector3( particle.GetPx(), particle.GetPy(), particle.GetPz()), obsE );
        TVector3 pbeta = mom4d.BoostVector();
        
        float eA = sqrt( pA.Mag2() + massA*massA );
        TLorentzVector mom4d_A( pA, eA );
        mom4d_A.Boost( -pbeta );
        phiA_pRF = mom4d_A.Phi();
        theA_pRF = mom4d_A.Theta();
        
        // grdaughter proton in the secondary-Lambda rest frame
        if( pid>=2 ){ // parent= Xi or Omega
          
          float obsE_dg = sqrt( massA*massA + pA.Mag()*pA.Mag() );
          TLorentzVector mom4d_dg( pA, obsE_dg );
          TVector3 dgbeta = mom4d_dg.BoostVector();
          float eAA = sqrt( pAA.Mag2() + massAA*massAA );
          TLorentzVector mom4d_AA( pAA, eAA );
          mom4d_AA.Boost( -dgbeta );
          phiAA_dgRF = mom4d_AA.Phi();
          theAA_dgRF = mom4d_AA.Theta();
        }else{ // parent= Lambda
          // skip here
          ftrkIdV0A.push_back(trkIdA);
          ftrkIdV0B.push_back(trkIdB);
        }
      }// pid<6
      
      //      histVar[pid][ 0]->Fill(mass);
      //      histVar[pid][ 1]->Fill(pt);
      //      histVar[pid][ 2]->Fill(phi);
      //      histVar[pid][ 3]->Fill(eta);
      //      histVar[pid][ 4]->Fill(rap);
      //
      //      histVar[pid][ 5]->Fill(massA);
      //      histVar[pid][ 6]->Fill(pA.Perp());
      //      histVar[pid][ 7]->Fill(pA.Phi());
      //      histVar[pid][ 8]->Fill(pA.Eta());
      //
      //      histVar[pid][ 9]->Fill(massB);
      //      histVar[pid][10]->Fill(pB.Perp());
      //      histVar[pid][11]->Fill(pB.Phi());
      //      histVar[pid][12]->Fill(pB.Eta());
      //
      //      histVar[pid][13]->Fill(massAA);
      //      histVar[pid][14]->Fill(pAA.Perp());
      //      histVar[pid][15]->Fill(pAA.Phi());
      //      histVar[pid][16]->Fill(pAA.Eta());
      
      if( pid<2 ){
        float lamM = 1.115683;
        //float myrap  = log( (sqrt(lamM*lamM+pt*pt*TMath::CosH(eta)*TMath::CosH(eta)) + pt*TMath::SinH(eta)) / sqrt(lamM*lamM+pt*pt) );
        float En = sqrt( lamM*lamM + particle.GetP()*particle.GetP() );
        float myrap = 0.5*log( (En + particle.GetPz()) / (En-particle.Pz()) );
//        histRap[0]->Fill(rap,myrap);
//        if( fabs(mass-lamM)<0.01 ) histRap[1]->Fill(rap,myrap);
      }
      
      // set track info. to tree
      // ------------------------------------------------------
      int trk_stat = pEve->set_crNtrk( cur_ntrk );
      if( trk_stat==-1 ){
        cout << endl <<"Exceeded a maximum of array. There exists missed tracks !!"<< endl << endl;
        cur_ntrk = 0;
        break;
      }
      
      // parent
      pEve->set_pdg (pdg);
      pEve->set_pt  (pt);
      pEve->set_phi (phi);
      pEve->set_eta (eta);
      pEve->set_rap (rap);
      pEve->set_mass(mass);
      pEve->set_dL  (dL);
      pEve->set_dca (dca);
      pEve->set_dgdca(dgdca);
      // add SIMDParticle
      pEve->set_SIMD_ldl      (SIMD_ldl);
      pEve->set_SIMD_l        (SIMD_l);
      pEve->set_SIMD_dl       (SIMD_dl);
      pEve->set_SIMD_chi2topo (SIMD_chi2topo);
      pEve->set_SIMD_chi2ndf  (SIMD_chi2ndf);
      pEve->set_SIMD_dca      (SIMD_dca);
      pEve->set_SIMD_dL       (SIMD_decaylength);
      pEve->set_SIMD_lT       (SIMD_lifetime);
      pEve->set_SIMD_mass     (SIMD_mass);
      pEve->set_SIMD_pt       (SIMD_pt);
      pEve->set_SIMD_phi      (SIMD_phi);
      pEve->set_SIMD_eta      (SIMD_eta);
      pEve->set_SIMD_rap      (SIMD_rapidity);
      pEve->set_SIMD_dca2D    (SIMD_dca2D);
      pEve->set_SIMD_pL       (SIMD_pathlength);
      pEve->set_SIMD_pdg      (SIMD_pdg);
      
      // daughter baryon
      pEve->set_pt_db   (pA.Perp());
      pEve->set_phi_db  (pA.Phi());
      pEve->set_eta_db  (pA.Eta());
      pEve->set_mass_db (massA);
      pEve->set_dca_db  (dca_A);
      pEve->set_dL_db   (dL_A);
      
      pEve->set_acmass  (mass_fake);
      pEve->set_dgdca_db(dgdca_db);
      
      //add more info (for pi<- lambda)
      pEve->set_chi2primaryA  (chi2primary_A);
      pEve->set_nhitsA        (nhits_A);
      pEve->set_rapidityA     (rapidity_A);
      pEve->set_nsigmaA       (0, nsigmapion_A);
      pEve->set_nsigmaA       (1, nsigmaproton_A);
      pEve->set_nsigmaA       (2, nsigmaKaon_A);
      pEve->set_dedxA         (dedx_A);
      pEve->set_m2A           (m2_A);
      pEve->set_pdgA          (pdg_A);
      
      //add SIMDParticle
      pEve->set_SIMD_A_ldl      (SIMD_A_ldl);
      pEve->set_SIMD_A_l        (SIMD_A_l);
      pEve->set_SIMD_A_dl       (SIMD_A_dl);
      pEve->set_SIMD_A_chi2topo (SIMD_A_chi2topo);
      pEve->set_SIMD_A_chi2ndf  (SIMD_A_chi2ndf);
      pEve->set_SIMD_A_dca      (SIMD_A_dca);
      pEve->set_SIMD_A_dL       (SIMD_A_decaylength);
      pEve->set_SIMD_A_lT       (SIMD_A_lifetime);
      pEve->set_SIMD_A_mass     (SIMD_A_mass);
      pEve->set_SIMD_A_pt       (SIMD_A_pt);
      pEve->set_SIMD_A_phi      (SIMD_A_phi);
      pEve->set_SIMD_A_eta      (SIMD_A_eta);
      pEve->set_SIMD_A_rap      (SIMD_A_rapidity);
      pEve->set_SIMD_A_dca2D    (SIMD_A_dca2D);
      pEve->set_SIMD_A_pL       (SIMD_A_pathlength);
      pEve->set_SIMD_A_pdg      (SIMD_A_pdg);
      // add more info
      pEve->set_chi2primaryB  (chi2primary_B);
      pEve->set_dcaB          (dca_B);
      pEve->set_nhitsB        (nhits_B);
      pEve->set_rapidityB     (rapidity_B);
      pEve->set_nsigmaB       (0, nsigmapion_B);
      pEve->set_nsigmaB       (1, nsigmaproton_B);
      pEve->set_nsigmaB       (2, nsigmaKaon_B);
      pEve->set_dedxB         (dedx_B);
      pEve->set_m2B           (m2_B);
      pEve->set_pdgB          (pdg_B);
      
      pEve->set_chi2primaryAA  (chi2primary_AA);
      pEve->set_dcaAA          (dca_AA);
      pEve->set_nhitsAA        (nhits_AA);
      pEve->set_rapidityAA     (rapidity_AA);
      pEve->set_nsigmaAA       (0, nsigmapion_AA);
      pEve->set_nsigmaAA       (1, nsigmaproton_AA);
      pEve->set_nsigmaAA       (2, nsigmaKaon_AA);
      pEve->set_dedxAA         (dedx_AA);
      pEve->set_m2AA           (m2_AA);
      pEve->set_pdgAA          (pdg_AA);
      
      pEve->set_chi2primaryAB  (chi2primary_AB);
      pEve->set_dcaAB          (dca_AB);
      pEve->set_nhitsAB        (nhits_AB);
      pEve->set_rapidityAB     (rapidity_AB);
      pEve->set_nsigmaAB       (0, nsigmapion_AB);
      pEve->set_nsigmaAB       (1, nsigmaproton_AB);
      pEve->set_nsigmaAB       (2, nsigmaKaon_AB);
      pEve->set_dedxAB         (dedx_AB);
      pEve->set_m2AB           (m2_AB);
      pEve->set_pdgAB          (pdg_AB);
      
      // daughter in parent rest frame
      pEve->set_phDBRF ( phiA_pRF );
      pEve->set_thDBRF ( theA_pRF );
      
      // grand daughter in daughter baryon rest frame
      pEve->set_phGDBRF ( phiAA_dgRF );
      pEve->set_thGDBRF ( theAA_dgRF );
      
      // daughter meson / grand daughters info. (truncated)
      pEve->set_ptB  ( static_cast<short>(1000.*pB.Perp())  );
      pEve->set_etaB ( static_cast<short>(1000.*pB.Eta())   );
      pEve->set_phiB ( static_cast<short>(1000.*pB.Phi())   );
      pEve->set_ptAA ( static_cast<short>(1000.*pAA.Perp()) );
      pEve->set_etaAA( static_cast<short>(1000.*pAA.Eta())  );
      pEve->set_phiAA( static_cast<short>(1000.*pAA.Phi())  );
      pEve->set_ptAB ( static_cast<short>(1000.*pAB.Perp()) );
      pEve->set_etaAB( static_cast<short>(1000.*pAB.Eta())  );
      pEve->set_phiAB( static_cast<short>(1000.*pAB.Phi())  );
      
      cur_ntrk++;
      
    }// track loop
    
    // put event info. to tree
    // =====================================
    pEve->set_runid ( brunid );
    pEve->set_vz( bVz );
    pEve->set_vx( bVx );
    pEve->set_vy( bVy );
    pEve->set_zdcRate( fPicoDst->event()->ZDCx() );
    
    pEve->set_refMult( crefmult );
    pEve->set_tofMult( ctofmult );
    pEve->set_countrefmult(countrefmult);
    
    if(fRunCentralityAnalysis){
      pEve->set_trgEff( reweight );
      pEve->set_centid( cent9 );
      pEve->set_refMultcor( refmultcor );
    }
    
    pEve->set_ntrk( cur_ntrk );
    
    if(fFlowAnalysis){
      pEve->set_isPileUp(isPileUp);
      pEve->set_fCentrality(centralityBin);
      pEve->set_psi_1_EPD  (psi_1_EPD);
      pEve->set_psi_1_EPD_0(psi_1_EPD_0);
      pEve->set_psi_1_EPD_1(psi_1_EPD_1);
      pEve->set_psi_1_EPD_2(psi_1_EPD_2);
      pEve->set_psi_1_EPD_3(psi_1_EPD_3);
      pEve->set_psi_1_EPD_4(psi_1_EPD_4);
      pEve->set_psi_1_EPD_5(psi_1_EPD_5);
      pEve->set_psi_1_EPD_6(psi_1_EPD_6);
      pEve->set_gweight(gweight);
    }
    ptree->Fill();
    pEve->InitValue();
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
  //if(_fill_lambda_tree){ld_chi2ndf = particle.Chi2()/particle.NDF();}
  
  KFParticleSIMD tempSIMDParticle(particle);
  float_v l,dl;
  KFParticleSIMD pv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
  tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl);
  fTMVAParticleParameters[iReader][nDaughterParticleCut + 1] = l[0]/dl[0];
  
  fTMVAParticleParameters[iReader][nDaughterParticleCut + 2] = l[0];
  //b
  //if(_fill_lambda_tree){ld_ldl = l[0]/dl[0];}
  
  tempSIMDParticle.SetProductionVertex(pv);
  fTMVAParticleParameters[iReader][nDaughterParticleCut + 3] = double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);
  //if(_fill_lambda_tree){ld_chi2topo = double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);}
  
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

  mMyFile->cd();
  ptree->Write();

//  hvtx->Write();
//  hvtxgood->Write();
//  hrefmult->Write();

//  if(fRunCentralityAnalysis) wrefmult->Write(); //refmultcor
//  if(fFlowAnalysis){
//    hrefmult_wPileup->Write();
//    hrefmult_woPileup->Write();
//    hCentrality->Write();
//    hCentrality_weight->Write();
//  }
//  //  ptree->Write();
//  //  for( int id=0; id<6; id++ ){
//  //    for( int k=0; k<nHVar; k++ ){
//  //      histVar[id][k]->Write();
//  //    }
//  //  }
//  //  for( int k=0; k<2; k++ ) histRap[k]->Write();
//  //
//  //  pileup check
//  hNBTofvsRM  ->Write();
//  hNBTofvsGRM ->Write();
//  hTofMvsRM   ->Write();
//  hTofMvsGRM  ->Write();
//  hRankvsGTrk ->Write();
//  hRankvsPTrk ->Write();
//  hRankvsNBTof->Write();
//  hGTrkvsPTrk ->Write();
//  hGTrkvsRM   ->Write();
//  hGRMvsRM    ->Write();
//  hNVpdvsRM   ->Write();
//  hNVpdvsGRM  ->Write();
//  hNVpdvsGTrk ->Write();
//  hNBEMCvsRM  ->Write();
//  hNBEMCvsGRM ->Write();
//  hNBEMCvsGTrk->Write();
//  //
  mMyFile->Close();
  delete ptree;
  delete mMyFile;
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

Bool_t StKFParticleAnalysisMaker::IsBadRun(Int_t RunId){
  
  vector<Int_t>::iterator iter = std::find( mBadRuns.begin(), mBadRuns.end(), RunId );
  return ( iter != mBadRuns.end() );
}

void StKFParticleAnalysisMaker::SetBadRun(){
  
  ifstream fin;
//  if(fsnn==3)  fin.open("/star/u/yjzhou19/polarization/3GeV_kfptree/jobs/badruns_3.list");
  if(fsnn==3)  fin.open("/star/u/yjzhou19/polarization/3GeV_kfptree/jobs/new_badruns_3.list");
  if(fsnn==27) fin.open("/star/u/yjzhou19/polarization/3GeV_kfptree/jobs/badruns_27.list");
  if( !fin ){
    cout <<"Could not open bad run list... no bad runs"<< endl;
  }
  
  Int_t runid = 0;
  while( fin >> runid ){
    mBadRuns.push_back(runid);
  }
  cout <<"StKFParticleAnalysisMaker::SetBadRun(): #-of-badruns = "<< mBadRuns.size() << endl;
}

void StKFParticleAnalysisMaker::SetHistograms(){
  hvtx      = new TH1F("hvtx",    "Vz;Vz(cm);Counts",400,0.,400);
  hvtxgood  = new TH1F("hvtxgood","Vz;Vz(cm);Counts",400,0.,400);
  hrefmult  = new TH1F("hrefmult", "refmult; hrefmult; N_{evt}", 600,0,600);
  wrefmult  = new TH1F("wrefmult", "refmult; wrefmult; N_{evt}", 600,0,600);
  hrefmult_wPileup  = new TH1F("hrefmult_wPileup", "refmult_wPileup; refmult; N_{evt}", 600,0,600);
  hrefmult_woPileup  = new TH1F("hrefmult_woPileup", "refmult_woPileup; refmult; N_{evt}", 600,0,600);
  hCentrality = new TH1F("hCentrality",";cent bin;N_{evt}", 9, -0.5, 8.5);
  hCentrality_weight = new TH1F("hCentrality_weight",";cent bin;N_{evt}", 9, -0.5, 8.5);
  
  
  for( int id=0; id<7; id++ ){
    int kid = 0;
    histVar[id][kid] = new TH1F( Form("histVar_%d_%d",id,kid), Form("histVar_%d_%d;mass [GeV/c^{2}] ;",id,kid), 500, 1, 2 );      kid = 1;
    histVar[id][kid] = new TH1F( Form("histVar_%d_%d",id,kid), Form("histVar_%d_%d;p_{T} [GeV/c] ;",   id,kid), 100, 0, 10 );     kid = 2;
    histVar[id][kid] = new TH1F( Form("histVar_%d_%d",id,kid), Form("histVar_%d_%d;#phi [rad] ;",      id,kid), 100, -3.2, 3.2 ); kid = 3;
    histVar[id][kid] = new TH1F( Form("histVar_%d_%d",id,kid), Form("histVar_%d_%d;#eta ;",            id,kid), 100, -1.1, 1.1 ); kid = 4;
    histVar[id][kid] = new TH1F( Form("histVar_%d_%d",id,kid), Form("histVar_%d_%d;rapidity ;",        id,kid), 100, -1.1, 1.1 ); kid = 5;
    
    // daughter A
    histVar[id][kid] = new TH1F( Form("histVar_%d_%d",id,kid), Form("histVar_%d_%d;mass [GeV/c^{2}] ;",id,kid), 200, 0.9, 1.4 );  kid = 6;
    histVar[id][kid] = new TH1F( Form("histVar_%d_%d",id,kid), Form("histVar_%d_%d;p_{T} [GeV/c] ;",   id,kid), 100, 0, 10 );     kid = 7;
    histVar[id][kid] = new TH1F( Form("histVar_%d_%d",id,kid), Form("histVar_%d_%d;#phi [rad] ;",      id,kid), 100, -3.2, 3.2 ); kid = 8;
    histVar[id][kid] = new TH1F( Form("histVar_%d_%d",id,kid), Form("histVar_%d_%d;#eta ;",            id,kid), 100, -1.1, 1.1 ); kid = 9;
    
    // daughter B (pi or K)
    histVar[id][kid] = new TH1F( Form("histVar_%d_%d",id,kid), Form("histVar_%d_%d;mass [GeV/c^{2}] ;",id,kid), 100, 0., 1. );    kid = 10;
    histVar[id][kid] = new TH1F( Form("histVar_%d_%d",id,kid), Form("histVar_%d_%d;p_{T} [GeV/c] ;",   id,kid), 100, 0, 10 );     kid = 11;
    histVar[id][kid] = new TH1F( Form("histVar_%d_%d",id,kid), Form("histVar_%d_%d;#phi [rad] ;",      id,kid), 100, -3.2, 3.2 ); kid = 12;
    histVar[id][kid] = new TH1F( Form("histVar_%d_%d",id,kid), Form("histVar_%d_%d;#eta ;",            id,kid), 100, -1.1, 1.1 ); kid = 13;
    
    // grdaughters
    histVar[id][kid] = new TH1F( Form("histVar_%d_%d",id,kid), Form("histVar_%d_%d;mass [GeV/c^{2}] ;",id,kid), 100, 0.8, 1.3 );  kid = 14;
    histVar[id][kid] = new TH1F( Form("histVar_%d_%d",id,kid), Form("histVar_%d_%d;p_{T} [GeV/c] ;",   id,kid), 100, 0, 10 );     kid = 15;
    histVar[id][kid] = new TH1F( Form("histVar_%d_%d",id,kid), Form("histVar_%d_%d;#phi [rad] ;",      id,kid), 100, -3.2, 3.2 ); kid = 16;
    histVar[id][kid] = new TH1F( Form("histVar_%d_%d",id,kid), Form("histVar_%d_%d;#eta ;",            id,kid), 100, -1.1, 1.1 );
  }
  
  // pileup check
  hNBTofvsRM   = new TH2F( Form("hNBTOFvsRM"),  Form("hNBTOFvsRM"),  100, 0, 800, 100, 0, 1200 );
  hNBTofvsGRM  = new TH2F( Form("hNBTOFvsGRM"), Form("hNBTOFvsGRM"), 100, 0, 800, 100, 0, 1200 );
  hTofMvsRM    = new TH2F( Form("hTofvsRM"),  Form("hTofvsRM"),  100, 0, 800, 100, 0, 4000 );
  hTofMvsGRM   = new TH2F( Form("hTofvsGRM"), Form("hTofvsGRM"), 100, 0, 800, 100, 0, 4000 );
  hRankvsGTrk  = new TH2F( Form("hRankvsGTrk"), Form("hRankvsGTrk"),   100, 0, 5000, 50, -10, 5 );
  hRankvsPTrk  = new TH2F( Form("hRankvsPTrk"), Form("hRankvsPTrk"),   100, 0, 1600, 50, -10, 5 );
  hRankvsNBTof = new TH2F( Form("hRankvsNBTOF"), Form("hRankvsNBTOF"), 100, 0, 1200, 50, -10, 5 );
  hGTrkvsPTrk  = new TH2F( Form("hGTrkvsPTrk"), Form("hGTrkvsPTrk"), 100, 0, 1600, 100, 0, 5000 );
  hGTrkvsRM    = new TH2F( Form("hGTrkvsRM"), Form("hGTrkvsRM"), 100, 0, 800, 100, 0, 5000 );
  hGRMvsRM     = new TH2F( Form("hGRMvsRM"), Form("hGRMvsRM"), 100, 0, 800, 100, 0, 800 );
  hNVpdvsRM    = new TH2F( Form("hNVpdvsRM"),   Form("hNVpdvsRM"),   100, 0,  800, 40, 0, 40 );
  hNVpdvsGRM   = new TH2F( Form("hNVpdvsGRM"),   Form("hNVpdvsGRM"),   100, 0,  800, 40, 0, 40 );
  hNVpdvsGTrk  = new TH2F( Form("hNVpdvsGTrk"), Form("hNVpdvsGTrk"), 100, 0, 5000, 40, 0, 40 );
  hNBEMCvsRM   = new TH2F( Form("hBEMCvsRM"), Form("hBEMCvsRM"),     100, 0,  800, 100, 0, 1200 );
  hNBEMCvsGRM  = new TH2F( Form("hBEMCvsGRM"), Form("hBEMCvsGRM"),     100, 0,  800, 100, 0, 1200 );
  hNBEMCvsGTrk = new TH2F( Form("hBEMCvsGTrk"), Form("hBEMCvsGTrk"), 100, 0, 5000, 100, 0, 1200 );
  
  for( int k=0; k<2; k++ ){
    histRap[k] = new TH2F( Form("histRap_%d",k), Form("histRap_%d",k), 200, -1, 1, 200, -1, 1 );
  }
}
