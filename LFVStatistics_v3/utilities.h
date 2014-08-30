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

TH1D* FindME(double l0, double l1, double Metl, double ll, double delPt, double Metl0);
TH1D* FindEM(double l0, double l1, double Metl, double ll, double delPt, double Metl0);
TH1D* FindSig_HTE(double l0, double l1, double Metl, double ll, double delPt, double Metl0);
TH1D* FindSig_HTM(double l0, double l1, double Metl, double ll, double delPt, double Metl0);
TH1D* FindSig_HTE_ME(double l0, double l1, double Metl, double ll, double delPt, double Metl0);
TH1D* FindSig_HTM_EM(double l0, double l1, double Metl, double ll, double delPt, double Metl0);
TH1D* FindSig(double l0, double l1, double Metl, double ll, double delPt, double Metl0);
TH1D* FindWrongSig(double l0, double l1, double Metl, double ll, double delPt, double Metl0);

TH1D* FindME(double l0, double l1, double Metl, double ll, int Jets);
TH1D* FindEM(double l0, double l1, double Metl, double ll, int Jets);
TH1D* FindSig(double l0, double l1, double Metl, double ll, int Jets);

TH1D* FindEM(double l0, double l1, double Metl, double ll, double delPt, double Metl0, int Jets);
TH1D* FindME(double l0, double l1, double Metl, double ll, double delPt, double Metl0, int Jets);
TH1D* FindSig(double l0, double l1, double Metl, double ll, double delPt, double Metl0, int Jets);

double PrintSideBandProbabilities(TH1D* h_EM,TH1D* h_ME,TH1D* h_B);

TH1D* getBlindHisto(TH1D* h_source,TString name = "default");

double pValue(TH1D* h_pdf, double value);
double Median(const TH1D * h1);

void drawSensitivity(TH1D* h_vanilla,TString filename,int nbins,double maxMu);
void drawSensitivity(TH1D* h_vanilla,int nbins,double maxMu);

void drawPoly(double* param, double* bins, int nbins);
TH1D* getPoly(double* param, double* bins, int nbins);


}


#endif /* UTILITIES_H_ */
