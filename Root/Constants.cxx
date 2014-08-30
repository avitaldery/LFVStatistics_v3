#include "LFVStatistics_v3/Constants.h"

#include "TRandom.h"
#include "TH1.h"
#include "TString.h"


TString PATH;
TString PATHnoJets;
TString PATHwithJets;
TString HISTO_NAME_ME;
TString HISTO_NAME_EM;

double COLL_MASS_MIN;
double COLL_MASS_MAX;
int COLL_MASS_NBINS;
int COLL_MASS_MinSigBin;
int COLL_MASS_MaxSigBin;
int TOY_MC_NUMOFEVENTS;
int LIKELIHOOD_PDF_NBINS;
double LIKELIHOOD_PDF_MIN;
double LIKELIHOOD_PDF_MAX;
int LIKELIHOODRATIO_PDF_NBINS;
double LIKELIHOODRATIO_PDF_MIN;
double LIKELIHOODRATIO_PDF_MAX;
int LIKELIHOODRATIOMU_PDF_NBINS;
double LIKELIHOODRATIOMU_PDF_MIN;
double LIKELIHOODRATIOMU_PDF_MAX;
int LIKELIHOOD_PDF_SEED;
double GAUSS_WIDTH;
double GAUSS_MEAN;
TH1D* TEST_GAUSSIAN;
TH1D* DATA_GAUSSIAN;
TRandom RAND;
unsigned int SEED;
//TH1D* H_DEFAULT;

void InitExterns()
{
	//path for root files and histogram names - CHANGE HERE
	PATH = "../../CutFlowOptimization/14_08_18_20GeVJets/SR_noJets/";
	PATHnoJets = "../../CutFlowOptimization/14_08_18_20GeVJets/SR_noJets/";
	PATHwithJets = "../../CutFlowOptimization/14_08_18_20GeVJets/SR_withJets/";

	HISTO_NAME_ME = "nom/ME_McollHiggs_Unblind";
	HISTO_NAME_EM = "nom/EM_McollHiggs_Unblind";
	//see functions FindEM etc. in utilities.cxx for the complete format

	COLL_MASS_MIN = 0;
	COLL_MASS_MAX = 500;
	COLL_MASS_NBINS = 125;
	COLL_MASS_MinSigBin = 6;
	COLL_MASS_MaxSigBin = 40;

	TOY_MC_NUMOFEVENTS=1000;
	LIKELIHOOD_PDF_NBINS=500;
	LIKELIHOOD_PDF_MIN=-200;
	LIKELIHOOD_PDF_MAX=200;
	LIKELIHOODRATIO_PDF_NBINS=1000;
	LIKELIHOODRATIO_PDF_MIN=-50;
	LIKELIHOODRATIO_PDF_MAX=50;
	LIKELIHOODRATIOMU_PDF_NBINS=550;
	LIKELIHOODRATIOMU_PDF_MIN=-5;
	LIKELIHOODRATIOMU_PDF_MAX=50;
	LIKELIHOOD_PDF_SEED=12;

	GAUSS_WIDTH = 20;
	GAUSS_MEAN = 300;


}
