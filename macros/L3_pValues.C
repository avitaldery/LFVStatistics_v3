#include "TCanvas.h"
#include "TH1.h"
#include "TString.h"
#include "TFile.h"
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

void L3_pValues(double l0, double l1, double Metl, double ll, int Jets)
{
	InitExterns();
	RAND.SetSeed(5);

	Data d(l0,l1,Metl,ll,Jets);
	TH1D* h_b = Likelihood::GetBGEstimation(d);
	TH1D* h_pValues = new TH1D("pvalues","pvalues",22,0,1.1);
	int numMC = 200;

	for (int i=0; i<numMC; i++)
	{
		if (i % (numMC/10) == 0){ cout << i*(100./numMC) << "%-" << flush;}
		Data dRand;
		dRand = Toys::ToyData(d,h_b);
		TH1D* h_b_rand = Likelihood::GetBGEstimation(dRand);
		TH1D* pdf = new TH1D("pdf","pdf",100000,0,100);
		int numMC2 = 500;
		for (int j=0; j<numMC2; j++)
		{
			Data dToy;
			dToy = Toys::ToyData(dRand,h_b_rand);
			double qZero = Likelihood::qZero(dToy);
			dToy.free();
			pdf->Fill(qZero,1./numMC2);
		}

		double obs = Likelihood::qZero(dRand);
		dRand.free();
//		pdf->Draw();
		double pV = utilities::pValue(pdf,obs);
		delete pdf;
		delete h_b_rand;
		h_pValues->Fill(pV,1./numMC);
	}

	h_pValues->Draw();
	TFile* file = new TFile("f5.root","RECREATE");
	h_pValues->Write("L3_pValues");
	file->Close();
}

