

#ifndef LIKELIHOOD_H_
#define LIKELIHOOD_H_

#include "TGraph.h"


namespace Likelihood
{

double LogL(Data d,double mu); //LogL for L3 discovery
double LogL(Data d,double mu,double muHat); //LogL for L2 discovery
double LogL_mu(Data d,double mu); //LogL for L3 limit (with mu>=muhat in the definitions)

TGraph* GetLambda3Graph(Data d);
TGraph* GetLambda2Graph(Data d);

double LogL3(Data d, double mu, double muHat);
double LogL2(Data d, double mu, double muHat);

double qZero(Data d);
double qZero(Data d, double muHat); //L2

double qMu(Data d,double mu);

double getB(double n,double m,double S,double mu); //analytic expression for b in one bin
TH1D* GetBGEstimation(Data d);




}
#endif /* LIKELIHOOD_H_ */
