// $Id: StFemtoDstMaker.h,v 1.1 2021/04/19 22:04:11 slan Exp $

#ifndef STAR_StFemtoDstMaker
#define STAR_StFemtoDstMaker
#ifndef StMaker_H
#include "StMaker.h"
#endif
#include "TTree.h"
#include "TFile.h"
class StPicoDst;
class StKFParticleInterface;
class StFemtoDstMaker : public StMaker {
 private:
  StPicoDst *fFemtoDst;
  TTree     *fFemtoTree;
  StKFParticleInterface *fStKFParticleInterface;
  TFile     *fOutFile;
 public: 
  StFemtoDstMaker(const char *name="FemtoDst"): StMaker(name), fFemtoDst(0), fFemtoTree(0), fStKFParticleInterface(0), fOutFile(0) {}
    virtual       ~StFemtoDstMaker();
  virtual Int_t Init();
  virtual Int_t Make();
  virtual Int_t Finish();
  virtual const char *GetCVS() const {
    static const char cvs[]="Tag $Name:  $ $Id: StFemtoDstMaker.h,v 1.1 2021/04/19 22:04:11 slan Exp $ built " __DATE__ " " __TIME__ ; 
    return cvs;
  }
  
  ClassDef(StFemtoDstMaker,0)
};
#endif
// $Log: StFemtoDstMaker.h,v $
// Revision 1.1  2021/04/19 22:04:11  slan
// update v1
//

