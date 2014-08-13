
#include "TCanvas.h"
#include "TH1.h"
#include "TString.h"
#include "TFile.h"
#include "TLine.h"
#include "TGraph.h"
#include <iostream>
#include <sstream>
#include "TStyle.h"
#include <fstream>

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


void Sens_discovery(TString datafile, TString signalfile)
{
	InitExterns();
	//get histos
	TFile* f = new TFile(datafile);
	TFile* f_sig = new TFile(signalfile);
	TH1D* h_EM =(TH1D*)f->Get(HISTO_NAME_EM);
	TH1D* h_ME =(TH1D*)f->Get(HISTO_NAME_ME);
	TH1D* h_sig =(TH1D*)f_sig->Get(HISTO_NAME_ME);

	Data d(h_EM,h_ME,h_sig);
	TH1D* h_b = Likelihood::GetBGEstimation(d);

	utilities::PrintSideBandProbabilities(d.m_hEM,d.m_hME,h_b);

	double muMax = 10;
	int numMC = 10000; int numbins = 100;
	TH1D* h_sens = new TH1D("sens","sens",numbins,0,muMax);

	for (int i=0; i<numMC; i++)
	{
		//		cout<<"event # "<< i <<endl;
		if (i % (numMC/10) == 0){ cout << i*(100./numMC) << "%-" << flush;}
		Data dRand;
		dRand = Toys::ToyData(d,h_b);
		double mu = minimization::GetMuSensitivity_discovery(dRand);
		h_sens->Fill(mu,1./numMC);
		dRand.free();
	}

	// write to file
	std::ofstream SensFile;
	SensFile.open ("SensFile.txt",ios_base::app);
	SensFile.close();
	// draw brazil plot
	utilities::drawSensitivity(h_sens,"SensFile.txt",numbins,muMax);
}



void Sens_discovery(double l0, double l1, double Metl, double ll, int Jets)
{
	InitExterns();
	//get base data
	Data d(l0,l1,Metl,ll,Jets);

	TH1D* h_b = Likelihood::GetBGEstimation(d);

	utilities::PrintSideBandProbabilities(d.m_hEM,d.m_hME,h_b);

	double muMax = 10;
	int numMC = 10000; int numbins = 100;
	TH1D* h_sens = new TH1D("sens","sens",numbins,0,muMax);

	for (int i=0; i<numMC; i++)
	{
//		cout<<"event # "<< i <<endl;
		if (i % (numMC/10) == 0){ cout << i*(100./numMC) << "%-" << flush;}
		Data dRand;
		dRand = Toys::ToyData(d,h_b);
		double mu = minimization::GetMuSensitivity_discovery(dRand);
		h_sens->Fill(mu,1./numMC);
		dRand.free();
	}

	// write to file
	char temp1[50];
	sprintf(temp1,",{ %2.0f , %2.0f , %1.1f , %1.1f , %d , ",l0,l1,Metl,ll,Jets);
	std::ofstream SensFile;
	SensFile.open ("SensFile.txt",ios_base::app);
	SensFile << temp1 ; SensFile.close();
	// draw brazil plot
	utilities::drawSensitivity(h_sens,"SensFile.txt",numbins,muMax);
}


void Sens_discovery(double l0, double l1, double Metl, double ll)
{
	InitExterns();
	//get base data
	Data d(l0,l1,Metl,ll);

	TH1D* h_b;
	h_b= Likelihood::GetBGEstimation(d);

	utilities::PrintSideBandProbabilities(d.m_hEM,d.m_hME,h_b);

	double muMax = 10;
	int numMC = 10000; int numbins = 100;
	TH1D* h_sens = new TH1D("sens","sens",numbins,0,muMax);

	for (int i=0; i<numMC; i++)
	{
		//		cout<<"event # "<< i <<endl;
		if (i % (numMC/10) == 0){ cout << i*(100./numMC) << "%-" << flush;}
		Data dRand;
		dRand = Toys::ToyData(d,h_b);
		double mu = minimization::GetMuSensitivity_discovery(dRand);
		h_sens->Fill(mu,1./numMC);
		dRand.free();
	}

	// write to file
	char temp1[50];
	sprintf(temp1,",{ %2.0f , %2.0f , %1.1f , %1.1f , ",l0,l1,Metl,ll);
	std::ofstream SensFile;
	SensFile.open ("SensFile.txt",ios_base::app);
	SensFile << temp1 ; SensFile.close();
	// draw brazil plot
	utilities::drawSensitivity(h_sens,"SensFile.txt",numbins,muMax);

}


void Sens_discovery(double l0, double l1, double Metl, double ll,
		double l0_2,double l1_2, double Metl_2, double ll_2)
{
	InitExterns();
	//get base data
	Data d(l0,l1,Metl,ll,l0_2,l1_2,Metl_2,ll_2);

	TH1D* h_b;
	h_b= Likelihood::GetBGEstimation(d);

	utilities::PrintSideBandProbabilities(d.m_hEM,d.m_hME,h_b);

	double muMax = 10;
	int numMC = 10000; int numbins = 100;
	TH1D* h_sens = new TH1D("sens","sens",numbins,0,muMax);

	for (int i=0; i<numMC; i++)
	{
		//		cout<<"event # "<< i <<endl;
		if (i % (numMC/10) == 0){ cout << i*(100./numMC) << "%-" << flush;}
		Data dRand;
		dRand = Toys::ToyData(d,h_b);
		double mu = minimization::GetMuSensitivity_discovery(dRand);
		h_sens->Fill(mu,1./numMC);
		dRand.free();
	}

	// write to file
	char temp1[50];
	sprintf(temp1,",{ %2.0f , %2.0f , %1.1f , %1.1f , ",l0,l1,Metl,ll);
	std::ofstream SensFile;
	SensFile.open ("SensFile.txt",ios_base::app);
	SensFile << temp1 ; SensFile.close();
	// draw brazil plot
	utilities::drawSensitivity(h_sens,"SensFile.txt",numbins,muMax);

}

void Sens_discovery(double l0, double l1, double Metl, double ll,int Jets,
		double l0_2,double l1_2, double Metl_2, double ll_2,int Jets_2)
{
	InitExterns();
	//get base data
	Data d(l0,l1,Metl,ll,Jets,l0_2,l1_2,Metl_2,ll_2,Jets_2);

	TH1D* h_b;
	h_b= Likelihood::GetBGEstimation(d);

	utilities::PrintSideBandProbabilities(d.m_hEM,d.m_hME,h_b);

	double muMax = 10;
	int numMC = 10000; int numbins = 100;
	TH1D* h_sens = new TH1D("sens","sens",numbins,0,muMax);

	for (int i=0; i<numMC; i++)
	{
		//		cout<<"event # "<< i <<endl;
		if (i % (numMC/10) == 0){ cout << i*(100./numMC) << "%-" << flush;}
		Data dRand;
		dRand = Toys::ToyData(d,h_b);
		double mu = minimization::GetMuSensitivity_discovery(dRand);
		h_sens->Fill(mu,1./numMC);
		dRand.free();
	}

	// write to file
	char temp1[50];
	sprintf(temp1,",{ %2.0f , %2.0f , %1.1f , %1.1f , ",l0,l1,Metl,ll);
	std::ofstream SensFile;
	SensFile.open ("SensFile.txt",ios_base::app);
	SensFile << temp1 ; SensFile.close();
	// draw brazil plot
	utilities::drawSensitivity(h_sens,"SensFile.txt",numbins,muMax);

}



void Sens_discovery(TH1D* hEM, TH1D* hME, TH1D* hSig)
{
	InitExterns();
	//get base data
	Data d(hEM,hME,hSig);
	TH1D* h_b = Likelihood::GetBGEstimation(d);

	utilities::PrintSideBandProbabilities(d.m_hEM,d.m_hME,h_b);

	double muMax = 10;
	int numMC = 10000; int numbins = 100;
	TH1D* h_sens = new TH1D("sens","sens",numbins,0,muMax);

	for (int i=0; i<numMC; i++)
	{
		//		cout<<"event # "<< i <<endl;
		if (i % (numMC/10) == 0){ cout << i*(100./numMC) << "%-" << flush;}
		Data dRand;
		dRand = Toys::ToyData(d,h_b);
		double mu = minimization::GetMuSensitivity_discovery(dRand);
		h_sens->Fill(mu,1./numMC);
		dRand.free();
	}

	// write to file
	std::ofstream SensFile;
	SensFile.open ("SensFile.txt",ios_base::app);
	SensFile.close();
	utilities::drawSensitivity(h_sens,"SensFile.txt",numbins,muMax);
}



