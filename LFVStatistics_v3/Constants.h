#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <map>
#include <utility>
#include "TRandom.h"
#include "TString.h"
#include "TH1.h"

using namespace std;

//path for root files directory
extern TString PATH;
extern TString PATHnoJets;
extern TString PATHwithJets;
extern TString HISTO_NAME_ME;
extern TString HISTO_NAME_EM;

//mass histogrma variables
extern double COLL_MASS_MIN;
extern double COLL_MASS_MAX;
extern int COLL_MASS_NBINS;
extern int COLL_MASS_MaxSigBin;
extern int COLL_MASS_MinSigBin;

//number of events to produce pdfs
extern int TOY_MC_NUMOFEVENTS;
extern int LIKELIHOOD_PDF_NBINS;
extern double LIKELIHOOD_PDF_MIN;
extern double LIKELIHOOD_PDF_MAX;
extern int LIKELIHOODRATIO_PDF_NBINS;
extern double LIKELIHOODRATIO_PDF_MIN;
extern double LIKELIHOODRATIO_PDF_MAX;
extern int LIKELIHOODRATIOMU_PDF_NBINS;
extern double LIKELIHOODRATIOMU_PDF_MIN;
extern double LIKELIHOODRATIOMU_PDF_MAX;
extern int LIKELIHOOD_PDF_SEED;

extern double GAUSS_WIDTH;
extern double GAUSS_MEAN;

extern TH1D* TEST_GAUSSIAN;
extern TH1D* DATA_GAUSSIAN;

extern unsigned int SEED;
extern TRandom RAND;
extern TH1D* H_DEFAULT;



#define DEBUG_LEVEL 5


void InitExterns();

#endif /* CONSTANTS_H */
