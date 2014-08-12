/*
 * minimization.h
 *
 *  Created on: Jul 31, 2014
 *      Author: avitald
 */

#ifndef MINIMIZATION_H_
#define MINIMIZATION_H_

namespace minimization{

double GetMuHat(Data d);
void Lminim(int &npar, double *gin, double &f, double *par, int iflag); //function for GetMuHat

double GetMuSensitivity_discovery(Data d);
double GetMuSensitivity_limit(Data d);
void ThreeSig(int &npar, double *gin, double &f, double *par, int iflag); //function for GetMuSensitivity_discovery
void TwoSig(int &npar, double *gin, double &f, double *par, int iflag); //function for GetMuSensitivity_limit


}
#endif /* MINIMIZATION_H_ */
