/*
 * Toys.h
 *
 *  Created on: Jul 31, 2014
 *      Author: avitald
 */

#ifndef TOYS_H_
#define TOYS_H_

#include "TFile.h"
#include "TH1.h"
#include "TRandom.h"

namespace Toys
{

Data ToyData(TH1D* bkg,TH1D* hsig); //generate random toy data from bkg histo
Data ToyData(Data d); //generate random toy data from data object (finds BG estimation)
Data ToyData(Data d, TH1D* h_b);
Data ToyData_Signal(TH1D* bkg,TH1D* hsig, double muE);
Data ToyData_Signal(Data d, double mu);
TH1D* drawFromHisto(TH1D* h_source,TString name);
TH1D* drawFromHistoPlusSignal(TH1D* h_source, TH1D* hsig, double mu, TString name);
TH1D* drawFromHistoGaus(TH1D* h_source, TString name, double sigma = 1);
TH1D* drawFromHisto2(TH1D* h_source, TString name, double sigma=1);

}


#endif /* TOYS_H_ */
