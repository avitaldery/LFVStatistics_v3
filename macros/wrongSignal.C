
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


void wrongSignal(TString signalfile)
{
	InitExterns();
	//get histos
	TFile* f_sig = new TFile(signalfile);
	TH1D* h_sig =(TH1D*)f_sig->Get(HISTO_NAME_ME);
	TH1D* h_wrongsig =(TH1D*)f_sig->Get(HISTO_NAME_EM);

	double wrong = h_wrongsig->Integral()/h_sig->Integral();

	cout << "wrong signal percentage = " << wrong << endl;

}


void wrongSignal(double l0, double l1, double Metl, double ll, double delPt, double Metl0)
{
	InitExterns();
	//get histos
	TH1D* h_sig = utilities::FindSig(l0,l1,Metl,ll,delPt,Metl0);
	TH1D* h_wrongsig = utilities::FindWrongSig(l0,l1,Metl,ll,delPt,Metl0);

	double wrong = h_wrongsig->Integral()/h_sig->Integral();

	// write to file
	char temp1[50];
	sprintf(temp1,",{ %2.0f , %2.0f , %1.1f , %1.1f , %1.1f , %1.1f , %4.4f ,",l0,l1,Metl,ll,delPt,Metl0,wrong);
	std::ofstream WrongSigFile;
	WrongSigFile.open ("WrongSignal.txt",ios_base::app);
	WrongSigFile << temp1 << endl; WrongSigFile.close();

}

