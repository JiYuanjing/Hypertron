#ifndef StCuts_h
#define StCuts_h

#include "KFParticle.h"
#include <vector> 
namespace Cuts{
  int const nCentralityBins = 16;
  int const nEPbins = 16;
  int const nVtxbins= 1;
}
struct StMixEvent{
  std::vector <KFParticle> d_v;
  std::vector <float> d_v_chi2prim;
  std::vector <float> d_v_dca;
  std::vector <float> d_v_m2;
  std::vector <float> d_v_z;
  std::vector <float> d_v_dEdx;
  std::vector <int> d_v_nhits;

  // std::vector <KFParticle> Lambda_v;
  std::vector <KFParticle> pi_v;
  std::vector <float> pi_v_chi2prim;
  std::vector <float> pi_v_dca;
  std::vector <float> pi_v_m2;
  std::vector <float> pi_v_sigma;
  std::vector <int> pi_v_nhits;
  std::vector <float> pi_v_mcpx;
  std::vector <float> pi_v_mcpy;
  std::vector <float> pi_v_mcpz;


  std::vector <KFParticle> p_v;
  std::vector <float> p_v_chi2prim;
  std::vector <float> p_v_dca;
  std::vector <float> p_v_m2;
  std::vector <float> p_v_sigma;
  std::vector <int> p_v_nhits;
  std::vector <float> p_v_mcpx;
  std::vector <float> p_v_mcpy;
  std::vector <float> p_v_mcpz;
 
  std::vector <KFParticle> ppi_v;

  float Vx;
  float Vy;
  float Vz;
};
#endif
