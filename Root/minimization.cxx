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

void Lminim_B(int &npar, double *gin, double &f, double *par, int iflag)
{
	// par[0] is muHat
	// par[1-nbins] is Bhat
	npar = npar; gin = gin; iflag = iflag;
	Data *data = (Data*)gMinuit->GetObjectFit();
	double l = Likelihood::LogL_B(*data,par);
//	delete data;
	f = l;
}

void ThreeSig(int &npar, double *gin, double &f, double *par, int iflag)
{
	//par[0] is mu
	npar = npar; gin = gin; iflag = iflag;
	Data *data = (Data*)gMinuit->GetObjectFit();

	double l = Likelihood::qMu(*data,par[0]);

	f = abs(l-10.273);//5.1875 //3sigma and 2sigma corresponding values
//	cout<<"mu = "<<par[0]<<", f = "<<f<<endl;
}

void TwoSig(int &npar, double *gin, double &f, double *par, int iflag)
{
	//par[0] is mu
	npar = npar; gin = gin; iflag = iflag;
	Data *data = (Data*)minuit->GetObjectFit();
	double l0sum = 0;
	int iterations = 20;
	for (int i=0; i<iterations; i++){
		Data dRandSig;
		dRandSig = Toys::ToyData_Signal(*data,par[0]);
		l0sum += Likelihood::qZero(dRandSig);
		dRandSig.free();
	}
	double l0 = l0sum/iterations;
	f = abs(l0-5.1875);
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

void GetMuHat_B(Data d, double* muHat_B, double* Errors)
{

	TMinuit *gMinuit = new TMinuit(d.m_nbins+1);
	gMinuit->SetObjectFit(&d);
	gMinuit->SetFCN(Lminim_B);

	Int_t ierflg = 0; gMinuit->SetPrintLevel(-1);
	// initialize the parameters:
	gMinuit->mnparm(0,"muHat",0,0.0001,-10,10,ierflg);
	TString varname;

	for (int j=1;j<=d.m_nbins;j++){
		varname.Form("B%d",j);
		double minV = (d.m_n[j-1] < d.m_m[j-1] ? d.m_n[j-1] : d.m_m[j-1]);
		double maxV = (d.m_n[j-1] < d.m_m[j-1] ? d.m_m[j-1] : d.m_n[j-1]);
		if (minV == maxV){maxV = maxV+1;}
		double startV = (minV+maxV)/2.;
//		cout << "n = " << d.m_n[j-1] << ", m = " << d.m_m[j-1] << ", startV = " << startV << endl;
		gMinuit->mnparm(j,varname,startV,0.0001,minV,maxV,ierflg);
	}

	gMinuit->SetMaxIterations(500);
	gMinuit->Migrad();

	for (int i=0; i<=d.m_nbins; i++){
		gMinuit->GetParameter(i,muHat_B[i],Errors[i]);
//		cout << "GetMuHat_B:: muHat_B[i] = " << muHat_B[i] << endl;
	}

	delete gMinuit;

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
