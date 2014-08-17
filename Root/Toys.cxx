#include "TFile.h"
#include "TH1.h"
#include "TRandom.h"

#include "LFVStatistics_v3/Constants.h"

#include "LFVStatistics_v3/Toys.h"
#include "LFVStatistics_v3/Data.h"
#include "LFVStatistics_v3/Likelihood.h"

namespace Toys{

Data ToyData(TH1D* bkg,TH1D* hsig)
{
	TH1D* h_brand = drawFromHistoGaus(bkg,"brand");
	TH1D* hEM = drawFromHisto(h_brand,"toyEM");
	TH1D* hME = drawFromHisto(h_brand,"toyME");
	Data d(hEM,hME,hsig);
	return d;
}

Data ToyData(Data d)
{
	TH1D* h_b = Likelihood::GetBGEstimation(d);
	TH1D* hEM = drawFromHisto2(h_b,"toyEM");
	TH1D* hME = drawFromHisto2(h_b,"toyME");
	delete h_b;
	Data newd(hEM,hME,d.m_hsig);
	return newd;
}

Data ToyData(Data d, TH1D* h_b)
{
	TH1D* hEM = drawFromHisto2(h_b,"toyEM");
	TH1D* hME = drawFromHisto2(h_b,"toyME");
	Data newd(hEM,hME,d.m_hsig);
	cout<< "bla4" << endl;
	return newd;
}


Data ToyData_Signal(TH1D* bkg,TH1D* hsig, double mu)
{
	TH1D* hEM = drawFromHisto(bkg,"toyEM");
	TH1D* hME = drawFromHistoPlusSignal(bkg,hsig,mu,"toyME");
	Data d(hEM,hME,hsig);
	return d;
}

Data ToyData_Signal(Data d, double mu)
{
	TH1D* h_b = Likelihood::GetBGEstimation(d);
	TH1D* hEM = NULL;
	TH1D* hME = NULL;
	if (mu>=0)
	{
		hEM = drawFromHisto2(h_b,"toyEM");
	}
	else { hME = drawFromHisto2(h_b,"toyME");}
	h_b->Add(d.m_hsig,mu);
	if (mu>=0)
	{
		hME = drawFromHisto2(h_b,"toyME");
	}
	else { hEM = drawFromHisto2(h_b,"toyEM");}
	delete h_b;
	Data newd(hEM,hME,d.m_hsig);
	return newd;
}


TH1D* drawFromHisto(TH1D* h_source,TString name)
{
	int nbins = h_source->GetXaxis()->GetNbins();
	const TArrayD* binsarray = h_source->GetXaxis()->GetXbins();
	const Double_t * bins = binsarray->GetArray();

	delete gDirectory->FindObject(name);
	TH1D* h_new = new TH1D(name,name,nbins,bins);
	for (int i=1; i<=nbins; i++){
		if (h_source->GetBinContent(i)<0.0001){h_new->SetBinContent(i,h_source->GetBinContent(i));}
		else{
			double v = RAND.Poisson(h_source->GetBinContent(i));
			if (v<0){v = 0;}
			h_new->SetBinContent(i,v);
		}
	}

	return h_new;
}

TH1D* drawFromHistoPlusSignal(TH1D* h_source, TH1D* hsig, double mu, TString name)
{
	int nbins = h_source->GetXaxis()->GetNbins();
	const TArrayD* binsarray = h_source->GetXaxis()->GetXbins();
	const Double_t * bins = binsarray->GetArray();

	delete gDirectory->FindObject(name);
	TH1D* h_new = new TH1D(name,name,nbins,bins);
	for (int i=1; i<=nbins; i++){
		double sum = h_source->GetBinContent(i)+mu*hsig->GetBinContent(i);
		if (sum<0.0001){h_new->SetBinContent(i,sum);}
		else{
			double v = RAND.Poisson(sum);
			if (v<0){v = 0;}
			h_new->SetBinContent(i,v);
		}
	}

	return h_new;
}

TH1D* drawFromHistoGaus(TH1D* h_source, TString name, double sigma)
{
	int nbins = h_source->GetXaxis()->GetNbins();
	const TArrayD* binsarray = h_source->GetXaxis()->GetXbins();
	const Double_t * bins = binsarray->GetArray();
	delete gDirectory->FindObject(name);
	TH1D* h_new = new TH1D(name,name,nbins,bins);

	for (int i=1; i<=nbins; i++){
		double n = h_source->GetBinContent(i);
		double err = sigma * (h_source->GetBinError(i));
		double v = RAND.Gaus(n,err);//sqrt(0.5n));
		if (v<0){v = 0;}
		h_new->SetBinContent(i,v);
	}
	return h_new;
}

TH1D* drawFromHisto2(TH1D* h_source, TString name, double sigma)
{	//do the gaus and poisson draws at once
	int nbins = h_source->GetXaxis()->GetNbins();
	const TArrayD* binsarray = h_source->GetXaxis()->GetXbins();const Double_t * bins = binsarray->GetArray();
//	delete gDirectory->FindObject(name);
	std::ostringstream strs;strs << RAND.Gaus(600,40)+0.1*RAND.Gaus(80,4)+0.01*RAND.Gaus(584,12);	std::string srand = strs.str();

	TH1D* h_new = new TH1D(name+srand,name+srand,nbins,bins);

	for (int i=1; i<=nbins; i++){
		double n = h_source->GetBinContent(i);
		double err = sigma * (h_source->GetBinError(i));
		//get gaussian distributed around b
		double vGaus = RAND.Gaus(n,err);//sqrt(0.5n));
		if (vGaus<0){vGaus = 0;}
		double vPoiss = RAND.Poisson(vGaus);
		if (vPoiss<0){vPoiss = 0;}
		h_new->SetBinContent(i,vPoiss);

	}

	return h_new;
}





}
