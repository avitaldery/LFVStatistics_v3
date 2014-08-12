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

void pdf_discovery(double l0, double l1, double Metl, double ll, int Jets, double mu)
{
	InitExterns();

	Data d(l0,l1,Metl,ll,Jets);

	//mu for testing
//	double mu = 2;

	TH1D* h_pdf = new TH1D("q_{0}|0","f(q_{0}|0)",2000,0,200);
	TH1D* h_pdf_sig = new TH1D("q_{0}|#mu","f(q_{0}|#mu)",2000,0,200);

	int numMC = 2000;
	//generate f(q0|mu)
	for (int i=0; i<numMC; i++)
	{
		if (i % (numMC/10) == 0){ cout << i*(100./numMC) << "%-" << flush;}
		Data dRandsig;
		dRandsig= Toys::ToyData_Signal(d,mu);
		//L3
//		double qZerotoysig = Likelihood::qZero(dRandsig);
		//L2
		double qZerotoysig = Likelihood::qZero(dRandsig,dRandsig.m_muHat);
		dRandsig.free();
		h_pdf_sig->Fill(qZerotoysig,1./numMC);
	}

	//generate f(q0|0)
	for (int i=0; i<numMC; i++)
	{
		if (i % (numMC/10) == 0){ cout << i*(100./numMC) << "%-" << flush;}
		Data dRand;
		dRand= Toys::ToyData(d);
		//L3
//		double qZerotoy = Likelihood::qZero(dRand);
		//L2
		double qZerotoy = Likelihood::qZero(dRand,dRand.m_muHat);
		dRand.free();
		h_pdf->Fill(qZerotoy,1./numMC);
	}

	double med = utilities::Median(h_pdf_sig);
	double pvalue = utilities::pValue(h_pdf,med);
	TLine* line  = new TLine(med,0,med,0.5); line->SetLineColor(kRed);

	cout<<"p-value = "<<pvalue<<endl;

	TF1 *chisquare = new TF1("chisuare","0.01*x^(-0.5)*exp(-x/2)/(2^(0.5))/(TMath::Gamma(0.5))",0.001,10);

	//Plotting and Legend **************************************************
	TCanvas* c2 = new TCanvas("pdfs","pdfs",600,600); c2 = c2;
	h_pdf->SetLineWidth(2); h_pdf_sig->SetLineColor(kGreen); h_pdf_sig->SetLineWidth(2);
	h_pdf->Draw(); h_pdf_sig->Draw("sames"); line->Draw();
	chisquare->Draw("sames");

	TLegend* leg = new TLegend(0.4,0.6,0.9,0.9);
	leg->SetFillColor(kWhite);leg->SetBorderSize(1);leg->SetTextSize(0.05);
	std::ostringstream strs;strs << mu;	std::string s_mu = strs.str();
	const char* c_mu = s_mu.c_str();char formula1[80];char formula2[80];
	strcpy(formula1,"f(q_{0}|0)");strcpy(formula2,"f(q_{0}|#mu=");
	strcat(formula2,c_mu); strcat(formula2,")");

	leg->AddEntry(h_pdf,formula1,"l");leg->AddEntry(h_pdf_sig,formula2,"l");
	leg->Draw();

}

void pdf_discovery(double l0, double l1, double Metl, double ll, double mu)
{
	InitExterns();

	Data d(l0,l1,Metl,ll);

	//mu for testing
//	double mu = 2;

	TH1D* h_pdf = new TH1D("q_{0}|0","f(q_{0}|0)",2000,0,200);
	TH1D* h_pdf_sig = new TH1D("q_{0}|#mu","f(q_{0}|#mu)",2000,0,200);

	int numMC = 10000;
	//generate f(q0|mu)
	for (int i=0; i<numMC; i++)
	{
		if (i % (numMC/10) == 0){ cout << i*(100./numMC) << "%-" << flush;}
		Data dRandsig;
		dRandsig= Toys::ToyData_Signal(d,mu);
		//L3
//		double qZerotoysig = Likelihood::qZero(dRandsig);
		//L2
		double qZerotoysig = Likelihood::qZero(dRandsig,dRandsig.m_muHat);
		dRandsig.free();
		h_pdf_sig->Fill(qZerotoysig,1./numMC);
	}

	//generate f(q0|0)
	for (int i=0; i<numMC; i++)
	{
		if (i % (numMC/10) == 0){ cout << i*(100./numMC) << "%-" << flush;}
		Data dRand;
		dRand= Toys::ToyData(d);
		//L3
//		double qZerotoy = Likelihood::qZero(dRand);
		//L2
		double qZerotoy = Likelihood::qZero(dRand,dRand.m_muHat);
		dRand.free();
		h_pdf->Fill(qZerotoy,1./numMC);
	}

	double med = utilities::Median(h_pdf_sig);
	double pvalue = utilities::pValue(h_pdf,med);
	TLine* line  = new TLine(med,0,med,0.5); line->SetLineColor(kRed);

	cout<<"p-value = "<<pvalue<<endl;

	TF1 *chisquare = new TF1("chisuare","0.01*x^(-0.5)*exp(-x/2)/(2^(0.5))/(TMath::Gamma(0.5))",0.001,10);

	//Plotting and Legend **************************************************
	TCanvas* c2 = new TCanvas("pdfs","pdfs",600,600); c2 = c2;
	h_pdf->SetLineWidth(2); h_pdf_sig->SetLineColor(kGreen); h_pdf_sig->SetLineWidth(2);
	h_pdf->Draw(); h_pdf_sig->Draw("sames"); line->Draw();
	chisquare->Draw("sames");

	TLegend* leg = new TLegend(0.4,0.6,0.9,0.9);
	leg->SetFillColor(kWhite);leg->SetBorderSize(1);leg->SetTextSize(0.05);
	std::ostringstream strs;strs << mu;	std::string s_mu = strs.str();
	const char* c_mu = s_mu.c_str();char formula1[80];char formula2[80];
	strcpy(formula1,"f(q_{0}|0)");strcpy(formula2,"f(q_{0}|#mu=");
	strcat(formula2,c_mu); strcat(formula2,")");

	leg->AddEntry(h_pdf,formula1,"l");leg->AddEntry(h_pdf_sig,formula2,"l");
	leg->Draw();

}

