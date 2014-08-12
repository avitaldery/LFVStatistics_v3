#include "TCanvas.h"
#include "TH1.h"
#include "TString.h"
#include "TFile.h"
#include "TLine.h"
#include "TGraph.h"
#include <iostream>
#include <sstream>

#include "LFVStatistics_v3/Data.h"
#include "LFVStatistics_v3/utilities.h"
#include "LFVStatistics_v3/Toys.h"
#include "LFVStatistics_v3/Constants.h"
#include "LFVStatistics_v3/Likelihood.h"
#include "LFVStatistics_v3/minimization.h"

using namespace utilities;
using namespace Toys;
using namespace Likelihood;
using namespace minimization;

void Sens_limit(double l0, double l1, double Metl, double ll, int Jets)
{
	InitExterns();
	//get base data
	Data d(l0,l1,Metl,ll,Jets);

	double muMax = 10;
	int numMC = 500; int numbins = 100;
	TH1D* h_sens = new TH1D("sens","sens",numbins,0,muMax);

	for (int i=0; i<numMC; i++)
	{
//		cout<<"event # "<< i <<endl;
		if (i % (numMC/10) == 0){ cout << i*(100./numMC) << "%-" << flush;}
		double mu = minimization::GetMuSensitivity_limit(d);
		h_sens->Fill(mu,1./numMC);
	}

	utilities::drawSensitivity(h_sens,numbins,muMax);
}

//void Sens_limit(TH1D* hEM, TH1D* hME, TH1D* hSig)
//{
//	InitExterns();
//	//get base data
//	Data d(hEM,hME,hSig);
//	TH1D* h_b = Likelihood::GetBGEstimation(d);
//
//	double muMax = 10;
//	int numMC = 1000; int numbins = 100;
//	TH1D* h_sens = new TH1D("sens","sens",numbins,0,muMax);
//
//	for (int i=0; i<numMC; i++)
//	{
//		//		cout<<"event # "<< i <<endl;
//		if (i % (numMC/10) == 0){ cout << i*(100./numMC) << "%-" << flush;}
//		Data dRand;
//		dRand = Toys::ToyData(d,h_b);
//		double mu = minimization::GetMuSensitivity_discovery(dRand);
//		h_sens->Fill(mu,1./numMC);
//		dRand.free();
//	}
//
//	utilities::drawSensitivity(h_sens,numbins,muMax);
//}
//
//void Sens_discovery(double l0, double l1, double Metl, double ll)
//{
//	InitExterns();
//	//get base data
//	Data d(l0,l1,Metl,ll);
//
//	TH1D* h_b;
//	h_b= Likelihood::GetBGEstimation(d);
//
//	double muMax = 10;
//	int numMC = 1000; int numbins = 100;
//	TH1D* h_sens = new TH1D("sens","sens",numbins,0,muMax);
//
//	for (int i=0; i<numMC; i++)
//	{
//		//		cout<<"event # "<< i <<endl;
//		if (i % (numMC/10) == 0){ cout << i*(100./numMC) << "%-" << flush;}
//		Data dRand;
//		dRand = Toys::ToyData(d,h_b);
//		double mu = minimization::GetMuSensitivity_discovery(dRand);
//		h_sens->Fill(mu,1./numMC);
//		dRand.free();
//	}
//
//	utilities::drawSensitivity(h_sens,numbins,muMax);
//}
//
//

