/*
 * utilities.h
 *
 *  Created on: Jul 31, 2014
 *      Author: avitald
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include "TH1.h"

namespace utilities{

TH1D* FindME(double l0, double l1, double Metl, double ll, int Jets);
TH1D* FindEM(double l0, double l1, double Metl, double ll, int Jets);
TH1D* FindSig(double l0, double l1, double Metl, double ll, int Jets);

double PrintSideBandProbabilities(TH1D* h_EM,TH1D* h_ME,TH1D* h_B);

TH1D* getBlindHisto(TH1D* h_source,TString name = "default");

double pValue(TH1D* h_pdf, double value);
double Median(const TH1D * h1);

void drawSensitivity(TH1D* h_vanilla,TString filename,int nbins,double maxMu);

}


#endif /* UTILITIES_H_ */
