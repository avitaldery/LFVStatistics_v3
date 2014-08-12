#include "TH1.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TLine.h"
#include "TLegend.h"
#include "TMath.h"
#include "TFile.h"
#include "TStyle.h"
#include "TDirectory.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <vector>
#include <stdio.h>
#include "TMinuit.h"

#include "LFVStatistics_v3/minimization.h"
#include "LFVStatistics_v3/Likelihood.h"
#include "LFVStatistics_v3/utilities.h"
#include "LFVStatistics_v3/Data.h"


using namespace Likelihood;
using namespace utilities;

namespace minimization{

TMinuit* minuit=0;


void Lminim(int &npar, double *gin, double &f, double *par, int iflag)
{
	//par[0] is muHat
	npar = npar; gin = gin; iflag = iflag;
	Data *data = (Data*)gMinuit->GetObjectFit();
	double l = Likelihood::LogL(*data,par[0]);
//	delete data;
	f = l;
}

void ThreeSig(int &npar, double *gin, double &f, double *par, int iflag)
{
	//par[0] is mu
	npar = npar; gin = gin; iflag = iflag;
	Data *data = (Data*)gMinuit->GetObjectFit();

	double l = Likelihood::qMu(*data,par[0]);

	f = abs(l-9);
//	cout<<"mu = "<<par[0]<<", f = "<<f<<endl;
}

void TwoSig(int &npar, double *gin, double &f, double *par, int iflag)
{
	//par[0] is mu
	npar = npar; gin = gin; iflag = iflag;
	Data *data = (Data*)minuit->GetObjectFit();
	Data dRandSig;
	dRandSig = Toys::ToyData_Signal(*data,par[0]);
	double l0 = Likelihood::qZero(dRandSig);

	dRandSig.free();
	f = abs(l0-4);
//	cout<<"mu = "<<par[0]<<", l0 = "<<l0<<endl;
}

double GetMuHat(Data d)
{

	TMinuit *gMinuit = new TMinuit(1);
	gMinuit->SetObjectFit(&d);
	gMinuit->SetFCN(Lminim);

	Int_t ierflg = 0; gMinuit->SetPrintLevel(-1);
	// initialize the parameters:
	gMinuit->mnparm(0,"muHat",0,0.0001,-10,10,ierflg);
	gMinuit->SetMaxIterations(500);
	gMinuit->Migrad();

	double muHat,error;
	gMinuit->GetParameter(0,muHat,error);
	delete gMinuit;

	return muHat;
}



double GetMuSensitivity_discovery(Data d)
{
	double mu = 0;
	TMinuit *gMinuit = new TMinuit(1);
	gMinuit->SetObjectFit(&d);
	gMinuit->SetFCN(ThreeSig);

	Int_t ierflg = 0; gMinuit->SetPrintLevel(-1);
	// initialize the parameters:
	double startValue;
	if (d.m_muHat<0){startValue = 0.01;}
	else {startValue = d.m_muHat+0.01;}
	gMinuit->mnparm(0,"mu",startValue,0.01,startValue,10,ierflg);
	gMinuit->SetMaxIterations(500);
	gMinuit->Migrad();

	double error;
	gMinuit->GetParameter(0,mu,error);

//	cout<<"*********************"<<endl;
//	cout<<"muHat = "<<d.m_muHat<<endl;
//	cout<<"mu = "<<mu<<endl;
//	cout<<"*********************"<<endl;

	delete gMinuit;
	return mu;

}

double GetMuSensitivity_limit(Data d)
{
	double mu = 0;
	minuit = new TMinuit(1);
	minuit->SetObjectFit(&d);
	minuit->SetFCN(TwoSig);

	Int_t ierflg = 0; minuit->SetPrintLevel(-1);
	// initialize the parameters:
	double startValue = 0.1;

	minuit->mnparm(0,"mu",startValue,0.001,startValue,10,ierflg);
	minuit->SetMaxIterations(500);
	minuit->Migrad();

	double error;
	minuit->GetParameter(0,mu,error);

//	cout<<"*********************"<<endl;
//	cout<<"mu = "<<mu<<endl;
//	cout<<"*********************"<<endl;

	delete minuit;
	return mu;

}














}
