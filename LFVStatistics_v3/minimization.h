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
void GetMuHat_B(Data d, double* muHat_B, double* Errors);
void GetMuHat_Poly(Data d, double* muHat_poly, double* Errors);
void Lminim(int &npar, double *gin, double &f, double *par, int iflag); //function for GetMuHat
void Lminim_poly(int &npar, double *gin, double &f, double *par, int iflag);
void Lminim_B(int &npar, double *gin, double &f, double *par, int iflag);

double GetMuSensitivity_discovery(Data d);
double GetMuSensitivity_limit(Data d);
void ThreeSig(int &npar, double *gin, double &f, double *par, int iflag); //function for GetMuSensitivity_discovery
void TwoSig(int &npar, double *gin, double &f, double *par, int iflag); //function for GetMuSensitivity_limit


}
#endif /* MINIMIZATION_H_ */
