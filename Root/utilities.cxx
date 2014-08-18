
#include "TH1.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TLine.h"
#include "TLegend.h"
#include "TMath.h"
#include "TFile.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <vector>
#include <stdio.h>
#include "TMinuit.h"

#include "LFVStatistics_v3/utilities.h"
#include "LFVStatistics_v3/Constants.h"

namespace utilities{


TH1D* FindME(double l0, double l1, double Metl, double ll, double delPt, double Metl0)
{
    char temp[50];
    sprintf(temp,"%2.0f_%2.0f_%1.1f_%1.1f_%1.1f_%1.1f",l0,l1,Metl,ll,delPt,Metl0);

    TFile* f=NULL;
    if( 1){ f = new TFile(PATH+"LFV."+temp+".root");}
//    if( Jets==1 ){ f = new TFile(PATH+"LFV_1jet."+temp+".root");}

    if (f == NULL){return NULL;}
    TH1D* h_ME =(TH1D*)f->Get(HISTO_NAME_ME);
    if (h_ME==NULL){cout<<"NULL histo"<<endl;}
    return h_ME;
}

TH1D* FindME(double l0, double l1, double Metl, double ll, int Jets)
{
    char temp[50];
    sprintf(temp,"%2.0f_%2.0f_%1.1f_%1.1f",l0,l1,Metl,ll);

    TFile* f=NULL;
    if( Jets==0 ){ f = new TFile(PATH+"LFV."+temp+".root");}
    if( Jets==1 ){ f = new TFile(PATH+"LFV_1jet."+temp+".root");}

    if (f == NULL){return NULL;}
    TH1D* h_ME =(TH1D*)f->Get(HISTO_NAME_ME);
    if (h_ME==NULL){cout<<"NULL histo"<<endl;}
    return h_ME;
}

TH1D* FindEM(double l0, double l1, double Metl, double ll, double delPt, double Metl0)
{
    char temp[50];
    sprintf(temp,"%2.0f_%2.0f_%1.1f_%1.1f_%1.1f_%1.1f",l0,l1,Metl,ll,delPt,Metl0);

    TFile* f=NULL;
    if( 1 ){ f = new TFile(PATH+"LFV."+temp+".root");}
//    if( Jets==1 ){ f = new TFile(PATH+"LFV_1jet."+temp+".root");}
    if (f == NULL){return NULL;}
    TH1D* h_EM =(TH1D*)f->Get(HISTO_NAME_EM);
    if (h_EM==NULL){cout<<"NULL histo"<<endl;}
    return h_EM;
}



TH1D* FindEM(double l0, double l1, double Metl, double ll, int Jets)
{
    char temp[50];
    sprintf(temp,"%2.0f_%2.0f_%1.1f_%1.1f",l0,l1,Metl,ll);

    TFile* f=NULL;
    if( Jets==0 ){ f = new TFile(PATH+"LFV."+temp+".root");}
    if( Jets==1 ){ f = new TFile(PATH+"LFV_1jet."+temp+".root");}
    if (f == NULL){return NULL;}
    TH1D* h_EM =(TH1D*)f->Get(HISTO_NAME_EM);
    if (h_EM==NULL){cout<<"NULL histo"<<endl;}
    return h_EM;
}

TH1D* FindSig(double l0, double l1, double Metl, double ll, double delPt, double Metl0)
{
    char temp[50];
    sprintf(temp,"%2.0f_%2.0f_%1.1f_%1.1f_%1.1f_%1.1f",l0,l1,Metl,ll,delPt,Metl0);

    TFile* fsig=NULL;
    if( 1 ){ fsig = new TFile(PATH+"HTM_"+temp+".root");}
//    if( Jets==1 ){ fsig = new TFile(PATH+"HTM_"+temp+".LFV_1jet.root");}

    if (fsig == NULL){return NULL;}
    TH1D* h_sig =(TH1D*)fsig->Get(HISTO_NAME_ME);
    if (h_sig==NULL){cout<<"NULL histo"<<endl;}

    return h_sig;
}


TH1D* FindSig(double l0, double l1, double Metl, double ll, int Jets)
{
    char temp[50];
    sprintf(temp,"%2.0f_%2.0f_%1.1f_%1.1f",l0,l1,Metl,ll);

    TFile* fsig=NULL;
    if( Jets==0 ){ fsig = new TFile(PATH+"HTM_"+temp+".root");}
    if( Jets==1 ){ fsig = new TFile(PATH+"HTM_"+temp+".LFV_1jet.root");}

    if (fsig == NULL){return NULL;}
    TH1D* h_sig =(TH1D*)fsig->Get(HISTO_NAME_ME);
    if (h_sig==NULL){cout<<"NULL histo"<<endl;}

    return h_sig;
}


void setBErrors(Data d,TH1D* h_b)
{
	for (int i=1; i<=d.m_nbins; i++)
	{
		double D = d.m_n[i-1] + d.m_m[i-1];
		double diff;
		if (d.m_muHat >= 0){diff = d.m_n[i-1] - d.m_m[i-1];}
		else {diff = d.m_m[i-1] - d.m_n[i-1];}
		double a = D*D-4*d.m_muHat*d.m_S[i-1]+4*d.m_muHat*d.m_muHat*d.m_S[i-1]*d.m_S[i-1];
		double sigmai = sqrt((D*sqrt(a)+D*D-2*d.m_muHat*d.m_S[i-1]*diff)/(8*sqrt(a)));
		if (isnan(sigmai)){sigmai=0.5*sqrt(d.m_n[i-1]+d.m_m[i-1]);}
		h_b->SetBinError(i,sigmai);
	}
}

double pValue(TH1D* h_pdf, double value)
{
	int bin_i = h_pdf->GetXaxis()->FindBin(value);
	int bin_end = h_pdf->GetXaxis()->GetNbins();

	if ((h_pdf->Integral(0.,bin_end)-1) > 0.001){cout<<"PDF Integral is not 1!"<<endl;}

	return h_pdf->Integral(bin_i,bin_end);
}

double Median(const TH1D * h1)
{
   int n = h1->GetXaxis()->GetNbins();
   std::vector<double>  x(n);
   h1->GetXaxis()->GetCenter( &x[0] );
   const double * y = h1->GetArray();
   // exclude underflow/overflows from bin content array y
   return TMath::Median(n, &x[0], &y[1]);
}

void drawSensitivity(TH1D* h_vanilla,TString filename,int nbins,double maxMu)
{
	std::ofstream SensFile;
	SensFile.open (filename,ios_base::app);

	TCanvas* csens = new TCanvas("sensitivity","sensitivity",600,600); csens = csens;

	h_vanilla->Draw();

	Double_t quantile;
	Double_t quantile_1sigma;
	Double_t quantile_1sigma2;
	Double_t quantile_2sigma;
	Double_t quantile_2sigma2;
	const Double_t prob = 0.5;
	const Double_t prob2 = 0.683;
	const Double_t prob3 = 0.317;
	const Double_t prob4 = 0.954;
	const Double_t prob5 = 0.046;
	h_vanilla->GetQuantiles(1,&quantile,&prob);
	h_vanilla->GetQuantiles(1,&quantile_1sigma,&prob2);
	h_vanilla->GetQuantiles(1,&quantile_1sigma2,&prob3);
	h_vanilla->GetQuantiles(1,&quantile_2sigma,&prob4);
	h_vanilla->GetQuantiles(1,&quantile_2sigma2,&prob5);

	TH1D* h_green = new TH1D("greenquantile","greenquantile",nbins,0,maxMu);
	int bin_med = h_green->GetXaxis()->FindBin(quantile);
	int bin_1sig = h_green->GetXaxis()->FindBin(quantile_1sigma);

	for (int k=bin_med;k<=bin_1sig;k++){
		h_green->SetBinContent(k,h_vanilla->GetBinContent(k));
	}
	h_green->SetFillColor(kGreen);
	h_green->Draw("sames");

	TH1D* h_green2 = new TH1D("greenquantile2","greenquantile2",nbins,0,maxMu);
	int bin_1sig2 = h_green2->GetXaxis()->FindBin(quantile_1sigma2);

	for (int k=bin_1sig2;k<=bin_med;k++){
		h_green2->SetBinContent(k,h_vanilla->GetBinContent(k));
	}
	h_green2->SetFillColor(kGreen);
	h_green2->Draw("sames");
	TH1D* h_yellow = new TH1D("yellowquantile","yellowquantile",nbins,0,maxMu);
	int bin_2sig = h_yellow->GetXaxis()->FindBin(quantile_2sigma);

	for (int k=bin_1sig;k<=bin_2sig;k++){
		h_yellow->SetBinContent(k,h_vanilla->GetBinContent(k));
	}
	h_yellow->SetFillColor(kYellow);
	h_yellow->Draw("sames");
	TH1D* h_yellow2 = new TH1D("yellowquantile2","yellowquantile2",nbins,0,maxMu);
	int bin_2sig2 = h_yellow2->GetXaxis()->FindBin(quantile_2sigma2);

	for (int k=bin_2sig2;k<=bin_1sig2;k++){
		h_yellow2->SetBinContent(k,h_vanilla->GetBinContent(k));
	}
	h_yellow2->SetFillColor(kYellow);
	h_yellow2->Draw("sames");

	//print values
	cout<<"median value = "<<quantile<<endl;
	cout<<"1sigma values = "<<quantile_1sigma2<<", "<<quantile_1sigma<<endl;
	cout<<"2sigma values = "<<quantile_2sigma2<<", "<<quantile_2sigma<<endl;

	 // write results to file
	SensFile << quantile << "," << quantile_1sigma2 << "," << quantile_1sigma << "," << quantile_2sigma2 << "," << quantile_2sigma << "}" << endl;

	SensFile.close();
}

void drawSensitivity(TH1D* h_vanilla,int nbins,double maxMu)
{
	TCanvas* csens = new TCanvas("sensitivity","sensitivity",600,600); csens = csens;

	h_vanilla->Draw();

	Double_t quantile, quantile_1sigma, quantile_1sigma2, quantile_2sigma, quantile_2sigma2;
	const Double_t prob = 0.5;
	const Double_t prob2 = 0.683;
	const Double_t prob3 = 0.317;
	const Double_t prob4 = 0.954;
	const Double_t prob5 = 0.046;
	h_vanilla->GetQuantiles(1,&quantile,&prob);
	h_vanilla->GetQuantiles(1,&quantile_1sigma,&prob2);
	h_vanilla->GetQuantiles(1,&quantile_1sigma2,&prob3);
	h_vanilla->GetQuantiles(1,&quantile_2sigma,&prob4);
	h_vanilla->GetQuantiles(1,&quantile_2sigma2,&prob5);

	TH1D* h_green = new TH1D("greenquantile","greenquantile",nbins,0,maxMu);
	int bin_med = h_green->GetXaxis()->FindBin(quantile);
	int bin_1sig = h_green->GetXaxis()->FindBin(quantile_1sigma);

	for (int k=bin_med;k<=bin_1sig;k++){
		h_green->SetBinContent(k,h_vanilla->GetBinContent(k));
	}
	h_green->SetFillColor(kGreen);
	h_green->Draw("sames");

	TH1D* h_green2 = new TH1D("greenquantile2","greenquantile2",nbins,0,maxMu);
	int bin_1sig2 = h_green2->GetXaxis()->FindBin(quantile_1sigma2);

	for (int k=bin_1sig2;k<=bin_med;k++){
		h_green2->SetBinContent(k,h_vanilla->GetBinContent(k));
	}
	h_green2->SetFillColor(kGreen);
	h_green2->Draw("sames");
	TH1D* h_yellow = new TH1D("yellowquantile","yellowquantile",nbins,0,maxMu);
	int bin_2sig = h_yellow->GetXaxis()->FindBin(quantile_2sigma);

	for (int k=bin_1sig;k<=bin_2sig;k++){
		h_yellow->SetBinContent(k,h_vanilla->GetBinContent(k));
	}
	h_yellow->SetFillColor(kYellow);
	h_yellow->Draw("sames");
	TH1D* h_yellow2 = new TH1D("yellowquantile2","yellowquantile2",nbins,0,maxMu);
	int bin_2sig2 = h_yellow2->GetXaxis()->FindBin(quantile_2sigma2);

	for (int k=bin_2sig2;k<=bin_1sig2;k++){
		h_yellow2->SetBinContent(k,h_vanilla->GetBinContent(k));
	}
	h_yellow2->SetFillColor(kYellow);
	h_yellow2->Draw("sames");

	//print values
	cout<<"median value = "<<quantile<<endl;
	cout<<"1sigma values = "<<quantile_1sigma2<<", "<<quantile_1sigma<<endl;
	cout<<"2sigma values = "<<quantile_2sigma2<<", "<<quantile_2sigma<<endl;

}

double PrintSideBandProbabilities(TH1D* h_EM,TH1D* h_ME,TH1D* h_B)
{
	int nbins = h_ME->GetXaxis()->GetNbins();
	int sigbin1 = h_ME->GetXaxis()->FindBin(100);
	int sigbin2 = h_ME->GetXaxis()->FindBin(150);
	double sigEM[nbins],sigME[nbins];
	double sigmaSqrEM = 0;
	double sigmaSqrME = 0;
	double maxSigME = 0;
	double maxSigEM = 0;
	int N = 0;
	int maxBinEM=0;
	int maxBinME=0;
	for (int bin = 1; bin<sigbin1; bin++){
		double sig = h_B->GetBinError(bin);
		if (sig>0){
			double b = h_B->GetBinContent(bin);
			double n = h_EM->GetBinContent(bin);
			double m = h_ME->GetBinContent(bin);
//			cout<<"sigma bin = "<<(n-b)/sig<<endl;
			sigEM[N] = (n-b)/sig;
			sigME[N] = (m-b)/sig;
			sigmaSqrEM += sigEM[N]*sigEM[N];
			sigmaSqrME += sigME[N]*sigME[N];
			if (maxSigEM < abs(sigEM[N])) {maxSigEM = abs(sigEM[N]); maxBinEM = bin;}
			if (maxSigME < abs(sigME[N])) {maxSigME = abs(sigME[N]); maxBinME = bin;}
			N++;
		}
	}
	for (int bin = sigbin2; bin<nbins; bin++){
		double sig = h_B->GetBinError(bin);
		if (sig>0){
			double b = h_B->GetBinContent(bin);
			double n = h_EM->GetBinContent(bin);
			double m = h_ME->GetBinContent(bin);
//			cout<<"sigma bin = "<<(n-b)/sig<<endl;
			sigEM[N] = (n-b)/sig;
			sigME[N] = (m-b)/sig;
			sigmaSqrEM += sigEM[N]*sigEM[N];
			sigmaSqrME += sigME[N]*sigME[N];
			if (maxSigEM < abs(sigEM[N])) {maxSigEM = abs(sigEM[N]); maxBinEM = bin;}
			if (maxSigME < abs(sigME[N])) {maxSigME = abs(sigME[N]); maxBinME = bin;}
			N++;
		}
	}

	sigmaSqrEM = sigmaSqrEM/(N);
	sigmaSqrME = sigmaSqrME/(N);

	double sigma = (sigmaSqrEM > sigmaSqrME ? sigmaSqrEM : sigmaSqrME);

	double maxSig = maxSigEM > maxSigME ? maxSigEM : maxSigME;
	int maxBin = maxSigEM > maxSigME ? maxBinEM : maxBinME;

	cout<<"sigma EM = "<<sqrt(sigmaSqrEM)<<endl;
	cout<<"sigma ME = "<<sqrt(sigmaSqrME)<<endl;
	cout<<"maximum single bin deviation = "<<maxSig<<" in bin "<<maxBin<<endl;

	return sigma;
}


}
