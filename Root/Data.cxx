
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

#include "LFVStatistics_v3/Data.h"
#include "LFVStatistics_v3/utilities.h"
#include "LFVStatistics_v3/Likelihood.h"
#include "LFVStatistics_v3/minimization.h"

using namespace std;
using namespace utilities;
using namespace Likelihood;
using namespace minimization;

//CONSTRUCTORS
Data::Data(TH1D* hEM, TH1D* hME, TH1D* hsig)
{
	m_numSR = 1;
	setEM(hEM);
	setME(hME);
	setSig(hsig);
//	setSigHTM(hsig);
//	setSigHTE(hsig);
	setEM2(NULL);
	setME2(NULL);
	m_nbins = m_numSR*(hEM->GetXaxis()->GetNbins());
	m_minBin = 1;
	m_maxBin = m_nbins;
	setN();
	setM();
	setS();
	setBins();
//	m_muHat = 0;
	setBlindHistos();
	setMuHat();

}

Data::Data(TH1D* hEM, TH1D* hME, TH1D* hsig, int minBin, int maxBin)
{
	m_numSR = 1;
	setEM(hEM);
	setME(hME);
	setSig(hsig);
//	setSigHTM(hsig);
//	setSigHTE(hsig);
	setEM2(NULL);
	setME2(NULL);
	m_nbins = m_numSR*(hEM->GetXaxis()->GetNbins());
	m_minBin = minBin;
	m_maxBin = maxBin;
	setN();
	setM();
	setS();
	setBins();
	m_muHat = 0;
	setBlindHistos();
	setMuHat();
}

Data::Data(double l0, double l1, double Metl, double ll, double delPt, double Metl0)
{
	m_numSR = 1;
	TH1D* h_EM  = utilities::FindEM(l0,l1,Metl,ll,delPt,Metl0);
	TH1D* h_ME = utilities::FindME(l0,l1,Metl,ll,delPt,Metl0);
	TH1D* h_sig = utilities::FindSig(l0,l1,Metl,ll,delPt,Metl0);
//	TH1D* h_sigHTM = utilities::FindSig_HTM(l0,l1,Metl,ll,delPt,Metl0);
//	TH1D* h_sigHTE = utilities::FindSig_HTE(l0,l1,Metl,ll,delPt,Metl0);


	setEM(h_EM);
	setME(h_ME);
	setSig(h_sig);
//	setSigHTM(h_sigHTM);
//	setSigHTE(h_sigHTE);
	setEM2(NULL);
	setME2(NULL);
	m_nbins = m_numSR*(h_EM->GetXaxis()->GetNbins());
	m_minBin = 1;
	m_maxBin = m_nbins;
	setN();
	setM();
	setS();
	setBins();
	m_muHat = 0;
	setBlindHistos();
	setMuHat();
}


Data::Data(double l0, double l1, double Metl, double ll, int Jets)
{
	m_numSR = 1;
	TH1D* h_EM  = utilities::FindEM(l0,l1,Metl,ll,Jets);
	TH1D* h_ME = utilities::FindME(l0,l1,Metl,ll,Jets);
	TH1D* h_sig = utilities::FindSig(l0,l1,Metl,ll,Jets);

	setEM(h_EM);
	setME(h_ME);
	setSig(h_sig);
	setEM2(NULL);
	setME2(NULL);
	m_nbins = m_numSR*(h_EM->GetXaxis()->GetNbins());
	m_minBin = 1;
	m_maxBin = m_nbins;
	setN();
	setM();
	setS();
	setBins();
	m_muHat = 0;
	setBlindHistos();
	setMuHat();
}

Data::Data(double l0, double l1, double Metl, double ll)
{
	m_numSR = 2;
	TH1D* h_EM  = utilities::FindEM(l0,l1,Metl,ll,0);
	TH1D* h_ME = utilities::FindME(l0,l1,Metl,ll,0);
	TH1D* h_sig = utilities::FindSig(l0,l1,Metl,ll,0);
	TH1D* h_EM1  = utilities::FindEM(l0,l1,Metl,ll,1);
	TH1D* h_ME1 = utilities::FindME(l0,l1,Metl,ll,1);
	TH1D* h_sig1 = utilities::FindSig(l0,l1,Metl,ll,1);

	setEM(h_EM);
	setME(h_ME);
	setEM2(h_EM1);
	setME2(h_ME1);

	m_nbins = m_numSR*(h_EM->GetXaxis()->GetNbins());
	setBins();

	setEM(h_EM,h_EM1);
	setME(h_ME,h_ME1);
	setSig(h_sig,h_sig1);

	m_minBin = 1;
	m_maxBin = m_nbins;
	setN();
	setM();
	setS();

	m_muHat = 0;
	setBlindHistos();
	setMuHat();
}

Data::Data(double l0, double l1, double Metl, double ll,
		double l0_2, double l1_2, double Metl_2, double ll_2)
{
	m_numSR = 2;
	TH1D* h_EM  = utilities::FindEM(l0,l1,Metl,ll,0);
	TH1D* h_ME = utilities::FindME(l0,l1,Metl,ll,0);
	TH1D* h_sig = utilities::FindSig(l0,l1,Metl,ll,0);
	TH1D* h_EM1  = utilities::FindEM(l0_2,l1_2,Metl_2,ll_2,1);
	TH1D* h_ME1 = utilities::FindME(l0_2,l1_2,Metl_2,ll_2,1);
	TH1D* h_sig1 = utilities::FindSig(l0_2,l1_2,Metl_2,ll_2,1);

	setEM(h_EM);
	setME(h_ME);
	setEM2(h_EM1);
	setME2(h_ME1);

	m_nbins = m_numSR*(h_EM->GetXaxis()->GetNbins());
	setBins();

	setEM(h_EM,h_EM1);
	setME(h_ME,h_ME1);
	setSig(h_sig,h_sig1);

	m_minBin = 1;
	m_maxBin = m_nbins;
	setN();
	setM();
	setS();

	m_muHat = 0;
	setBlindHistos();
	setMuHat();
}

Data::Data(double l0, double l1, double Metl, double ll, double delPt, double Metl0,
				double l0_2, double l1_2, double Metl_2, double ll_2, double delPt_2, double Metl0_2)
{
	m_numSR = 2;
	TH1D* h_EM  = utilities::FindEM(l0,l1,Metl,ll,delPt,Metl0,0);
	TH1D* h_ME = utilities::FindME(l0,l1,Metl,ll,delPt,Metl0,0);
	TH1D* h_sig = utilities::FindSig(l0,l1,Metl,ll,delPt,Metl0,0);
	TH1D* h_EM1  = utilities::FindEM(l0_2,l1_2,Metl_2,ll_2,delPt_2,Metl0_2,1);
	TH1D* h_ME1 = utilities::FindME(l0_2,l1_2,Metl_2,ll_2,delPt_2,Metl0_2,1);
	TH1D* h_sig1 = utilities::FindSig(l0_2,l1_2,Metl_2,ll_2,delPt_2,Metl0_2,1);

	setEM(h_EM);
	setME(h_ME);
	setEM2(h_EM1);
	setME2(h_ME1);

	m_nbins = m_numSR*(h_EM->GetXaxis()->GetNbins());
	setBins();

	setEM(h_EM,h_EM1);
	setME(h_ME,h_ME1);
	setSig(h_sig,h_sig1);

	m_minBin = 1;
	m_maxBin = m_nbins;
	setN();
	setM();
	setS();

	m_muHat = 0;
	setBlindHistos();
	setMuHat();
}


Data::Data(double l0, double l1, double Metl, double ll, int Jets,
		double l0_2, double l1_2, double Metl_2, double ll_2, int Jets_2)
{
	m_numSR = 2;
	TH1D* h_EM  = utilities::FindEM(l0,l1,Metl,ll,Jets);
	TH1D* h_ME = utilities::FindME(l0,l1,Metl,ll,Jets);
	TH1D* h_sig = utilities::FindSig(l0,l1,Metl,ll,Jets);
	TH1D* h_EM1  = utilities::FindEM(l0_2,l1_2,Metl_2,ll_2,Jets_2);
	TH1D* h_ME1 = utilities::FindME(l0_2,l1_2,Metl_2,ll_2,Jets_2);
	TH1D* h_sig1 = utilities::FindSig(l0_2,l1_2,Metl_2,ll_2,Jets_2);

	setEM(h_EM);
	setME(h_ME);
	setEM2(h_EM1);
	setME2(h_ME1);

	m_nbins = m_numSR*(h_EM->GetXaxis()->GetNbins());
	setBins();

	setEM(h_EM,h_EM1);
	setME(h_ME,h_ME1);
	setSig(h_sig,h_sig1);

	m_minBin = 1;
	m_maxBin = m_nbins;
	setN();
	setM();
	setS();

	m_muHat = 0;
	setBlindHistos();
	setMuHat();
}

Data::Data(Data &d)
{
	setEM(d.m_hEM);
	setME(d.m_hME);
	setSig(d.m_hsig);
//	setSigHTM(d.m_hsigHTM);
//	setSigHTE(d.m_hsigHTE);
	setEMBlind(d.m_hEMblind);
	setMEBlind(d.m_hMEblind);
	setEM2(NULL);
	setME2(NULL);
	m_nbins = d.m_nbins;
	m_minBin = d.m_minBin;
	m_maxBin = d.m_maxBin;
	setN(d.m_n);
	setM(d.m_m);
	setS(d.m_S);
	setBins(d.m_Bins);
	m_muHat = 0;
//	setBlindHistos();
	setMuHat(d.m_muHat);
}

//Data::~Data()
//{
//
//}

void Data::free()
{
	if (m_hEM) {delete m_hEM;}
	if (m_hME) {delete m_hME;}
	if (m_hEMblind) {delete m_hEMblind;}
	if (m_hMEblind) {delete m_hMEblind;}
	if (m_SR2_hEM) {delete m_SR2_hEM;}
	if (m_SR2_hME) {delete m_SR2_hME;}
//	if (m_hsigHTM) {delete m_hsigHTM;}
//	if (m_hsigHTE) {delete m_hsigHTE;}
	delete[] m_Bins;
	delete[] m_n;
	delete[] m_m;
	delete[] m_S;
}

//Likelihood Graphs
void Data::DrawL2Graph()
{
	TGraph* p = Likelihood::GetLambda2Graph(*this); p->Draw("AC*");
	TLine* line  = new TLine(-5,9,5,9); line->SetLineColor(kRed); line->Draw();
	TLine* vline  = new TLine(0,-1,0,20); vline->SetLineColor(kBlack);vline->SetLineStyle(2); vline->Draw();
}

void Data::DrawL3Graph()
{
	TGraph* p = Likelihood::GetLambda3Graph(*this);	p->Draw("AC*");
	TLine* line  = new TLine(-5,9,5,9); line->SetLineColor(kRed); line->Draw();
	TLine* vline  = new TLine(0,-1,0,20); vline->SetLineColor(kBlack);vline->SetLineStyle(2); vline->Draw();
}

void Data::setBlindHistos()
{
	TString tempEM = m_hEM->GetName();
	TString nameEM = tempEM + "_blind";
	TString tempME = m_hME->GetName();
	TString nameME = tempME + "_blind";
//	delete gDirectory->FindObject(nameEM);
//	delete gDirectory->FindObject(nameME);

	m_hEMblind = new TH1D(nameEM,nameEM,m_nbins,m_Bins);
	m_hMEblind = new TH1D(nameME,nameME,m_nbins,m_Bins);
	int bin1 = m_hEM->GetXaxis()->FindBin(100);
	int bin2 = m_hEM->GetXaxis()->FindBin(150);

	for (int j=1; j<bin1; j++){
		m_hEMblind->SetBinContent(j,m_hEM->GetBinContent(j));
		m_hMEblind->SetBinContent(j,m_hME->GetBinContent(j));
	}
	for (int j=bin2; j<=m_nbins; j++){
		m_hEMblind->SetBinContent(j,m_hEM->GetBinContent(j));
		m_hMEblind->SetBinContent(j,m_hME->GetBinContent(j));
	}
	m_hEMblind->SetLineColor(kGreen+2);
	m_hMEblind->SetLineColor(kBlue);
}


void Data::setMuHat()
{
	m_muHat = minimization::GetMuHat(*this);
}

void Data::setMuHat(double muHat)
{
	m_muHat = muHat;
}

void Data::appendSignalME(double strength)
{
	m_hME->Add(m_hsig,strength);
	setMuHat();
	setM();
}

void Data::setEM(TH1D* hEM1, TH1D* hEM2)
{
	m_hEM = new TH1D("hEM2","hEM2",m_nbins,m_Bins);
	for (int i=1; i<=m_nbins/m_numSR; i++)
	{
		m_hEM->SetBinContent(i,hEM1->GetBinContent(i));
	}
	for (int i=m_nbins/m_numSR+1; i<=m_nbins; i++)
	{
		m_hEM->SetBinContent(i,hEM2->GetBinContent(i-m_nbins/m_numSR));

	}
}

void Data::setME(TH1D* hME1, TH1D* hME2)
{
	m_hME = new TH1D("hME2","hME2",m_nbins,m_Bins);
	for (int i=1; i<=m_nbins/m_numSR; i++)
	{
		m_hME->SetBinContent(i,hME1->GetBinContent(i));
	}
	for (int i=m_nbins/m_numSR+1; i<=m_nbins; i++)
	{
		m_hME->SetBinContent(i,hME2->GetBinContent(i-m_nbins/m_numSR));

	}
}

void Data::setSig(TH1D* hs1, TH1D* hs2)
{
	m_hsig = new TH1D("hsig2","hsig2",m_nbins,m_Bins);
	for (int i=1; i<=m_nbins/m_numSR; i++)
	{
		m_hsig->SetBinContent(i,hs1->GetBinContent(i));
	}
	for (int i=m_nbins/m_numSR+1; i<=m_nbins; i++)
	{
		m_hsig->SetBinContent(i,hs2->GetBinContent(i-m_nbins/m_numSR));

	}
}


