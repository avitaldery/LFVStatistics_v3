

#ifndef LIKELIHOOD_H_
#define LIKELIHOOD_H_

#include "TGraph.h"


namespace Likelihood
{

double LogL(Data d,double mu); //L3
double LogL_poly(Data d,double mu,double* muHat_poly);
double LogL_B(Data d,double* muHat_B);
double LogL(Data d,double mu,double muHat); //L2

TGraph* GetLambda3Graph(Data d);
TGraph* GetLambda2Graph(Data d);

double qZero(Data d);
double qZero(Data d, double muHat); //L2

double qMu(Data d,double mu);
double qMu(Data d,double mu, double muHat); //L@

double getB(double n,double m,double S,double mu); //analytic expression for b in one bin
double getB(double n,double m,double S,double mu,double f);

TH1D* GetBGEstimation(Data d);
TH1D* GetBGEstimation(Data d,double* coeff);// with polynomial function





}
#endif /* LIKELIHOOD_H_ */
