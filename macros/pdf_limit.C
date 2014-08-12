#include "TCanvas.h"
#include "TH1.h"
#include "TString.h"
#include "TFile.h"
#include "TLine.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TF1.h"
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

void pdf_limit(double l0, double l1, double Metl, double ll, int Jets, double mu)
{
	InitExterns();

	Data d(l0,l1,Metl,ll,Jets);

	//mu for testing
//	double mu = 1.4;

	TH1D* h_pdf_zero = new TH1D("q_{#mu}|0","f(q_{#mu}|0)",1000,0,100);
	TH1D* h_pdf_mu = new TH1D("q_{#mu}|#mu","f(q_{#mu}|#mu)",1000,0,100);

	int numMC = 1000;
	//generate f(qMu|mu)
	for (int i=0; i<numMC; i++)
	{
		if (i % (numMC/10) == 0){ cout << i*(100./numMC) << "%-" << flush;}
		Data dRandsig;
		dRandsig= Toys::ToyData_Signal(d,mu);
		double qmu_mu = Likelihood::qMu(dRandsig,mu);
		dRandsig.free();
		h_pdf_mu->Fill(qmu_mu,1./numMC);
	}

	//generate f(qMu|0)
	for (int i=0; i<numMC; i++)
	{
		if (i % (numMC/10) == 0){ cout << i*(100./numMC) << "%-" << flush;}
		Data dRand;
		dRand= Toys::ToyData(d);
		double qmu_zero = Likelihood::qMu(dRand,mu);
		dRand.free();
		h_pdf_zero->Fill(qmu_zero,1./numMC);
	}

	double med = utilities::Median(h_pdf_zero);
//	double obs = Likelihood::qZero(d);
	double pvalue = utilities::pValue(h_pdf_mu,med);
	TLine* line  = new TLine(med,0,med,0.5); line->SetLineColor(kRed);

	cout<<"p-value = "<<pvalue<<endl;

	TF1 *chisquare = new TF1("chisuare","0.05*x^(-0.5)*exp(-x/2)/(2^(0.5))/(TMath::Gamma(0.5))",0.001,10);

	//Plotting and Legend **************************************************
	TCanvas* c2 = new TCanvas("pdfs","pdfs",600,600); c2 = c2;
	h_pdf_zero->SetLineWidth(2); h_pdf_mu->SetLineColor(kGreen); h_pdf_mu->SetLineWidth(2);
	h_pdf_mu->Draw(); h_pdf_zero->Draw("sames");  line->Draw();
	chisquare->Draw("sames");

	TLegend* leg = new TLegend(0.4,0.6,0.9,0.9);
	leg->SetFillColor(kWhite);leg->SetBorderSize(1);leg->SetTextSize(0.05);
	std::ostringstream strs;strs << mu;	std::string s_mu = strs.str();
	const char* c_mu = s_mu.c_str();char formula1[80];char formula2[80];
	strcpy(formula1,"f(q_{#mu}|0)");strcpy(formula2,"f(q_{#mu}|#mu=");
	strcat(formula2,c_mu); strcat(formula2,")");

	leg->AddEntry(h_pdf_zero,formula1,"l");leg->AddEntry(h_pdf_mu,formula2,"l");
	leg->Draw();

}


void pdf_limit(double l0, double l1, double Metl, double ll, double mu)
{
	InitExterns();

	Data d(l0,l1,Metl,ll);

	//mu for testing
//	double mu = 1.4;

	TH1D* h_pdf_zero = new TH1D("q_{#mu}|0","f(q_{#mu}|0)",2000,0,1000);
	TH1D* h_pdf_mu = new TH1D("q_{#mu}|#mu","f(q_{#mu}|#mu)",2000,0,1000);

	int numMC = 1000;
	//generate f(qMu|mu)
	for (int i=0; i<numMC; i++)
	{
		if (i % (numMC/10) == 0){ cout << i*(100./numMC) << "%-" << flush;}
		Data dRandsig;
		dRandsig= Toys::ToyData_Signal(d,mu);
		double qmu_mu = Likelihood::qMu(dRandsig,mu);
		dRandsig.free();
		h_pdf_mu->Fill(qmu_mu,1./numMC);
	}

	//generate f(qMu|0)
	for (int i=0; i<numMC; i++)
	{
		if (i % (numMC/10) == 0){ cout << i*(100./numMC) << "%-" << flush;}
		Data dRand;
		dRand= Toys::ToyData(d);
		double qmu_zero = Likelihood::qMu(dRand,mu);
		dRand.free();
		h_pdf_zero->Fill(qmu_zero,1./numMC);
	}

	double med = utilities::Median(h_pdf_zero);
//	double obs = Likelihood::qZero(d);
	double pvalue = utilities::pValue(h_pdf_mu,med);
	TLine* line  = new TLine(med,0,med,0.5); line->SetLineColor(kRed);

	cout<<"p-value = "<<pvalue<<endl;

	TF1 *chisquare = new TF1("chisuare","0.05*x^(-0.5)*exp(-x/2)/(2^(0.5))/(TMath::Gamma(0.5))",0.001,10);

	//Plotting and Legend **************************************************
	TCanvas* c2 = new TCanvas("pdfs","pdfs",600,600); c2 = c2;
	h_pdf_zero->SetLineWidth(2); h_pdf_mu->SetLineColor(kGreen); h_pdf_mu->SetLineWidth(2);
	h_pdf_mu->Draw(); h_pdf_zero->Draw("sames");  line->Draw();
	chisquare->Draw("sames");

	TLegend* leg = new TLegend(0.4,0.6,0.9,0.9);
	leg->SetFillColor(kWhite);leg->SetBorderSize(1);leg->SetTextSize(0.05);
	std::ostringstream strs;strs << mu;	std::string s_mu = strs.str();
	const char* c_mu = s_mu.c_str();char formula1[80];char formula2[80];
	strcpy(formula1,"f(q_{#mu}|0)");strcpy(formula2,"f(q_{#mu}|#mu=");
	strcat(formula2,c_mu); strcat(formula2,")");

	leg->AddEntry(h_pdf_zero,formula1,"l");leg->AddEntry(h_pdf_mu,formula2,"l");
	leg->Draw();

}

