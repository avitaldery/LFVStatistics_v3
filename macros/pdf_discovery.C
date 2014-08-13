#include "TCanvas.h"
#include "TH1.h"
#include "TString.h"
#include "TFile.h"
#include "TLine.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TMath.h"
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

	int numbins = 1000; double pdfmax = 200;
	TH1D* h_pdf = new TH1D("q_{0}|0","f(q_{0}|0)",numbins,0,pdfmax);
	TH1D* h_pdf_sig = new TH1D("q_{0}|#mu","f(q_{0}|#mu)",numbins,0,pdfmax);

	int numMC = 10000;
	//generate f(q0|mu)
	for (int i=0; i<numMC; i++)
	{
		if (i % (numMC/5) == 0){ cout << i*(50./numMC) << "%-" << flush;}
		Data dRandsig;
		dRandsig= Toys::ToyData_Signal(d,mu);
		//L3
		double qZerotoysig = Likelihood::qZero(dRandsig);
		//L2
//		double qZerotoysig = Likelihood::qZero(dRandsig,dRandsig.m_muHat);
		dRandsig.free();
		h_pdf_sig->Fill(qZerotoysig,1./numMC);
	}

	//generate f(q0|0)
	for (int i=0; i<numMC; i++)
	{
		if (i % (numMC/5) == 0){ cout << 50+i*(50./numMC) << "%-" << flush;}
		Data dRand;
		dRand= Toys::ToyData(d);
		//L3
		double qZerotoy = Likelihood::qZero(dRand);
		//L2
//		double qZerotoy = Likelihood::qZero(dRand,dRand.m_muHat);
		dRand.free();
		h_pdf->Fill(qZerotoy,1./numMC);
	}

	double med = utilities::Median(h_pdf_sig);
	double pvalue = utilities::pValue(h_pdf,med);
	TLine* line  = new TLine(med,0,med,0.5); line->SetLineColor(kRed);

	cout<<"p-value = "<<pvalue<<endl;

	TF1 *chisquare = new TF1("chisquare","0.5*x^(-0.5)*exp(-x/2)/((2*pi)^(0.5))*(TMath::Gamma(0.5))",0.001,40);
	TH1D* half_chisquare  = new TH1D("chisquare","chisquare",numbins,0,pdfmax);
	double x1 = half_chisquare->GetXaxis()->GetBinCenter(1);
	double chis1 = 0.5 + 0.25*TMath::Power(x1,-0.5)*exp(-x1/2)/(TMath::Power(TMath::TwoPi(),0.5)*TMath::Gamma(0.5));
	half_chisquare->SetBinContent(1,chis1);
	for (int j=2; j<numbins; j++)
	{
		double x = half_chisquare->GetXaxis()->GetBinCenter(j);
		double chis = 0.25*TMath::Power(x,-0.5)*exp(-x/2)/(TMath::Power(TMath::TwoPi(),0.5)*TMath::Gamma(0.5));
		half_chisquare->SetBinContent(j,chis);
	}
	//Plotting and Legend **************************************************
	TCanvas* c2 = new TCanvas("pdfs","pdfs",600,600); c2 = c2;
	h_pdf->SetLineWidth(2); h_pdf_sig->SetLineColor(kGreen); h_pdf_sig->SetLineWidth(2);
	h_pdf->Draw(); h_pdf_sig->Draw("sames"); line->Draw();
	half_chisquare->Draw("sames");
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
	int numbins = 1000; double pdfmax = 200;
	TH1D* h_pdf = new TH1D("q_{0}|0","f(q_{0}|0)",numbins,0,pdfmax);
	TH1D* h_pdf_sig = new TH1D("q_{0}|#mu","f(q_{0}|#mu)",numbins,0,pdfmax);

	int numMC = 10000;
	//generate f(q0|mu)
	for (int i=0; i<numMC; i++)
	{
		if (i % (numMC/5) == 0){ cout << i*(50./numMC) << "%-" << flush;}
		Data dRandsig;
		dRandsig= Toys::ToyData_Signal(d,mu);
		//L3
		double qZerotoysig = Likelihood::qZero(dRandsig);
		//L2
//		double qZerotoysig = Likelihood::qZero(dRandsig,dRandsig.m_muHat);
		dRandsig.free();
		h_pdf_sig->Fill(qZerotoysig,1./numMC);
	}

	//generate f(q0|0)
	for (int i=0; i<numMC; i++)
	{
		if (i % (numMC/5) == 0){ cout << 50+i*(50./numMC) << "%-" << flush;}
		Data dRand;
		dRand= Toys::ToyData(d);
		//L3
		double qZerotoy = Likelihood::qZero(dRand);
		//L2
//		double qZerotoy = Likelihood::qZero(dRand,dRand.m_muHat);
		dRand.free();
		h_pdf->Fill(qZerotoy,1./numMC);
	}

	double med = utilities::Median(h_pdf_sig);
	double pvalue = utilities::pValue(h_pdf,med);
	TLine* line  = new TLine(med,0,med,0.5); line->SetLineColor(kRed);

	cout<<"p-value = "<<pvalue<<endl;

//	TF1 *chisquare = new TF1("chisuare","0.1*x^(-0.5)*exp(-x/2)/(2^(0.5))/(TMath::Gamma(0.5))",0.001,10);
	TH1D* h_chisquare  = new TH1D("chisquare","chisquare",numbins,0,pdfmax);
	for (int j=1; j<=numbins; j++)
	{
		double x = h_chisquare->GetXaxis()->GetBinCenter(j);
		double chis = 0.5*TMath::Power(x,-0.5)*exp(-x/2)/(TMath::Power(2,0.5)*TMath::Gamma(0.5));
		h_chisquare->SetBinContent(j,chis);
	}
	//Plotting and Legend **************************************************
	TCanvas* c2 = new TCanvas("pdfs","pdfs",600,600); c2 = c2;
	h_pdf->SetLineWidth(2); h_pdf_sig->SetLineColor(kGreen); h_pdf_sig->SetLineWidth(2);
	h_pdf->Draw(); h_pdf_sig->Draw("sames"); line->Draw();
	h_chisquare->Draw("sames");

	TLegend* leg = new TLegend(0.4,0.6,0.9,0.9);
	leg->SetFillColor(kWhite);leg->SetBorderSize(1);leg->SetTextSize(0.05);
	std::ostringstream strs;strs << mu;	std::string s_mu = strs.str();
	const char* c_mu = s_mu.c_str();char formula1[80];char formula2[80];
	strcpy(formula1,"f(q_{0}|0)");strcpy(formula2,"f(q_{0}|#mu=");
	strcat(formula2,c_mu); strcat(formula2,")");

	leg->AddEntry(h_pdf,formula1,"l");leg->AddEntry(h_pdf_sig,formula2,"l");
	leg->Draw();

}

