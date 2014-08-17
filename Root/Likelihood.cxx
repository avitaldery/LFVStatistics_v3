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
#include "TGraph.h"

#include "LFVStatistics_v3/Likelihood.h"
#include "LFVStatistics_v3/Data.h"
#include "LFVStatistics_v3/utilities.h"

namespace Likelihood{

TH1D* GetBGEstimation(Data d)
{	//analytic function of n,m,S,muHat
	delete gDirectory->FindObject("hb");
	TH1D* hb = new TH1D("hb","hb",d.m_nbins,d.m_Bins);
	for (int i=1; i<=d.m_nbins; i++){
		double b = getB(d.m_n[i-1],d.m_m[i-1],d.m_S[i-1],d.m_muHat);
		hb->SetBinContent(i,b);
	}
	utilities::setBErrors(d,hb);
	return hb;
}

double getB(double n,double m,double S,double mu)
{	//analytic function of n,m,S,muHat per bin
	double b;
	double D = n+m;
	double a = D+2*abs(mu)*S;
	if (mu >= 0){
		b = (sqrt(TMath::Power(a,2)-8*mu*S*m)+D-2*mu*S)/4;
	}
	else{
		b = (sqrt(TMath::Power(a,2)+8*mu*S*n)+D+2*mu*S)/4;
	}
	return b;
}


double LogL(Data d,double mu)
//for L3 Discovery
{
	double l=0;
	for (int i=d.m_minBin-1; i<d.m_maxBin; i++){
		//bhathat
		double b = getB(d.m_n[i],d.m_m[i],d.m_S[i],mu);
		if (d.m_n[i]+d.m_m[i] < 0.01){continue;}
		if (mu>=0)
		{
			l -= 2*(TMath::Log(TMath::Poisson(d.m_n[i],b))+TMath::Log(TMath::Poisson(d.m_m[i],b+mu*d.m_S[i])));
			if (l>10000000000) {cout<< "LogL:: i = " << i << ", n = " << d.m_n[i] << ", m = " << d.m_m[i] << ", b = " << b << ", S = " << d.m_S[i]<<endl;}
		}
		else
		{
			l -= 2*(TMath::Log(TMath::Poisson(d.m_n[i],b-mu*d.m_S[i]))+TMath::Log(TMath::Poisson(d.m_m[i],b)));
		}
	}
	return l;
}

double LogL_B(Data d,double* muHat_B)
// for L3 Discovery
{
	double l=0;
	for (int i=d.m_minBin-1; i<d.m_maxBin; i++){
		if (d.m_n[i]+d.m_m[i] < 0.01){continue;}
		if (muHat_B[0]>=0)
		{
			l -= 2*(TMath::Log(TMath::Poisson(d.m_n[i],muHat_B[i+1]))+TMath::Log(TMath::Poisson(d.m_m[i],muHat_B[i+1]+muHat_B[0]*d.m_S[i])));
		}
		else
		{
			l -= 2*(TMath::Log(TMath::Poisson(d.m_n[i],muHat_B[i+1]-muHat_B[0]*d.m_S[i]))+TMath::Log(TMath::Poisson(d.m_m[i],muHat_B[i+1])));
		}
		if (l>10000000000)
		{
			cout<< "LogL_B:: l = " << l << ", i = " << i << ", n = " << d.m_n[i] << ", m = " << d.m_m[i] << ", b = " << muHat_B[i+1] << ", mu = " << muHat_B[0] << ", S = " << d.m_S[i]<<endl;
			break;
		}
	}
	return l;
}



double LogL(Data d,double mu,double muHat)
//for L2 Discovery
{
	double l=0;
	for (int i=d.m_minBin-1; i<d.m_maxBin; i++){
		double bhat = getB(d.m_n[i],d.m_m[i],d.m_S[i],muHat);
		if (d.m_n[i]+d.m_m[i] < 0.01){continue;}
		if(muHat>=0)//if (mu>=0)
		{
			l -= 2*(TMath::Log(TMath::Poisson(d.m_n[i],bhat))+TMath::Log(TMath::Poisson(d.m_m[i],bhat+mu*d.m_S[i])));
		}
		else
		{
			l -= 2*(TMath::Log(TMath::Poisson(d.m_n[i],bhat-mu*d.m_S[i]))+TMath::Log(TMath::Poisson(d.m_m[i],bhat)));
		}
	}
	return l;
}

double qZero(Data d)
{
	double q0 = 0;
	if (1)//d.m_muHat>=0)
	{
		double denom = LogL(d,d.m_muHat);
		double num = LogL(d,0);
		q0 = num - denom;
		if (q0<0) {cout<<"negative q0! num = "<< num << ", denom = " << denom << ", muHat = " << d.m_muHat << endl;}
	}
	else {q0 = 0;}
	return q0;
}

double qZero(Data d, double muHat)
{//for L2
	double q0 = 0;
	if (1)//d.m_muHat>=0)
	{
		double denom = LogL(d,muHat);
		double num = LogL(d,0,muHat);
		q0 = num - denom;
	}
	else {q0 = 0;}
	return q0;
}


double qMu(Data d,double mu)
{
	double q_mu = 0;
	if (1)//mu>=d.m_muHat)
	{
		double denom = LogL(d,d.m_muHat);
		double num = LogL(d,mu);
//		cout << "denom = " << denom << ", num = " << num << endl;
		q_mu = num - denom;
	}
	else {q_mu = 0;}
	return q_mu;
}

double qMu(Data d,double mu, double muHat)
{ // for L2
	double q_mu = 0;
	if (1)//mu>=d.m_muHat)
	{
		double denom = LogL(d,muHat);
		double num = LogL(d,mu,muHat);
//		cout << "denom = " << denom << ", num = " << num << endl;
		q_mu = num - denom;
	}
	else {q_mu = 0;}
	return q_mu;
}

TGraph* GetLambda3Graph(Data d)
{
	int nMus = 100;
	double x[nMus];
	double y[nMus];

	d.setMuHat();

	for(int i=0; i<nMus; i++)
	{
		double mu = i*20./nMus - 10;

		// set denominator value for normalization
		x[i] = mu;
		y[i] = qMu(d,mu);
	}

	TGraph* g = new TGraph(nMus,x,y);
	g->SetName("LogL3");
	g->GetXaxis()->SetTitle("#mu");
	g->GetYaxis()->SetTitle("-2Log[#lambda(#mu)]");
	g->GetYaxis()->SetRangeUser(-1,20);
	return g;
}

TGraph* GetLambda2Graph(Data d)
{
	int nMus = 100;
	double x[nMus];
	double y[nMus];

	d.setMuHat();

	for(int i=0; i<nMus; i++)
	{
		double mu = i*20./nMus - 10;

		// set denominator value for normalization
		x[i] = mu;
		y[i] = qMu(d,mu,d.m_muHat);
	}

	TGraph* g = new TGraph(nMus,x,y);
	g->SetName("LogL2");
	g->GetXaxis()->SetTitle("#mu");
	g->GetYaxis()->SetTitle("-2Log[L2_{#mu}]");
	g->GetYaxis()->SetRangeUser(-1,20);
	return g;
}



















}

