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

void Start(double l0, double l1, double Metl, double ll)
{
	InitExterns();

	Data d(l0,l1,Metl,ll,0);
	cout<<"d.m_muHat = "<<d.m_muHat<<endl;
	TCanvas* c1 = new TCanvas("c1","c1",600,600); c1 = c1;
//	d.DrawEMME();
//	d.m_hsig->Draw("sames");

	double muHatB[d.m_nbins+1];
	double errors[d.m_nbins+1];
	minimization::GetMuHat_B(d,muHatB,errors);

	TH1D* h_b = new TH1D("b","b",d.m_nbins,d.m_Bins);
	for (int i=1; i<=d.m_nbins; i++)
	{
		h_b->SetBinContent(i,muHatB[i]);
	}

	TH1D* h_b_analytic = Likelihood::GetBGEstimation(d);
	h_b->Draw();
	h_b_analytic->SetLineColor(kRed); h_b_analytic->Draw("sames");


	Data d1(l0,l1,Metl,ll,1);
//	TCanvas* c2 = new TCanvas("c2","c2",600,600); c2 = c2;
//	d1.DrawEMME();
//	d1.m_hsig->Draw("sames");


	Data d2(l0,l1,Metl,ll);
//	TCanvas* c3 = new TCanvas("c3","c3",600,600); c3 = c3;
//	d2.DrawEMME();
//	d2.m_hsig->Draw("sames");
//	cout<<"d2: nbins = " << d2.m_nbins << endl;
//	cout<<"Bins[0] = " << d2.m_Bins[0] << ", Bins[1] = " << d2.m_Bins[1] <<endl;
//	cout<<"Bins[44] = " << d2.m_Bins[44] << ", Bins[45] = " << d2.m_Bins[45] <<endl;




	Data dRand;
//	dRand = Toys::ToyData(d.m_hEM,d.m_hsig);
//	dRand = Toys::ToyData_Signal(d,2);
	dRand = Toys::ToyData(d);

//	TCanvas* c2 = new TCanvas("c2","c2",600,600); c2 = c2;
//	dRand.DrawEMME();

//	TCanvas* c3 = new TCanvas("c3","c3",600,600); c3 = c3;
//	cout<<"dRand.m_muHat = "<<dRand.m_muHat<<endl;

//	dRand.setMuHat();
//	dRand.DrawL3Graph();
//	cout<<"muHat = "<<dRand.m_muHat<<endl;

//	TH1D* h_b = minimization::GetBGEstimation(d);

//	TCanvas* c4 = new TCanvas("c4","c4",600,600); c4 = c4;
//	dRand.DrawL2Graph();
//	d.DrawEMME();
//	h_b->SetLineColor(kBlack); h_b->SetLineWidth(2); h_b->Draw("e1");

//	Data dRandSignal;
//	double mu = 3;
//	dRandSignal = Toys::ToyData_Signal(d.m_hEM,d.m_hsig,mu);
//	dRandSignal.DrawL3Graph();
//
//	TCanvas* c5 = new TCanvas("c5","c5",600,600); c5 = c5;
//	Data dRandSignal2;
//	double mu2 = 6;
//	dRandSignal2 = Toys::ToyData_Signal(d.m_hEM,d.m_hsig,mu2);
//	dRandSignal2.DrawL3Graph();
//	dRandSignal.setMuHat();

//	cout<<"muHatsignal = "<<dRandSignal.m_muHat<<endl;
//
//
//	TCanvas* c6 = new TCanvas("c6","c6",600,600); c6 = c6;
//	dRandSignal.DrawL2Graph();

}

