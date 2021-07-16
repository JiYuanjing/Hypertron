#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TFile.h"


///weighting function for MC particles
TF1 *levyfit4;
TF1 *levyfit5;
TF1 *levyfit6;
TF1 *levyfit7;
TF1 *levyfit8;
TF1 *t_quadr;
TF1 *t_quadr0;
TF1 *t_quadr1;
TF1 *bolt0;
TF1 *bolt1;
TF1 *bolt2;

//make mc flat distribution
TFile*  fgpt_0;
TH2F*  g_pt_fine_in;


