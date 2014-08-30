
#ifndef DATA_H_
#define DATA_H_

#include "TH1.h"
#include "TGraph.h"


class Data	: public TObject
{
public:
	//CONSTRUCTORS
	//default constructor (sets muHat to ZERO)
	Data(): m_hEM(NULL), m_hME(NULL), m_hsig(NULL),  m_hEMblind(NULL), m_hMEblind(NULL), m_SR2_hEM(NULL),
		m_SR2_hME(NULL), m_Bins(NULL), m_muHat(0) {}
//	~Data();
	Data(TH1D* h1,TH1D* h2,TH1D* hsig);
	Data(TH1D* hEM, TH1D* hME, TH1D* hsig, int minBin, int maxBin);
	//copy constructor
	Data(Data&);
	//constructor from cut values
	Data(double l0, double l1, double Metl, double ll, double delPt, double Metl0);
	Data(double l0, double l1, double Metl, double ll, int Jets);
	Data(double l0, double l1, double Metl, double ll);
	Data(double l0, double l1, double Metl, double ll,
			double l0_2, double l1_2, double Metl_2, double ll_2);
	Data(double l0, double l1, double Metl, double ll, double delPt, double Metl0,
				double l0_2, double l1_2, double Metl_2, double ll_2, double delPt_2, double Metl0_2);
	Data(double l0, double l1, double Metl, double ll,int Jets,
				double l0_2, double l1_2, double Metl_2, double ll_2,int Jets_2);
	void setEM(TH1D* hEM){ m_hEM = hEM; m_hEM->SetLineColor(kGreen+2);}
	void setME(TH1D* hME){ m_hME = hME; m_hME->SetLineColor(kBlue);}
	void setEMBlind(TH1D* hEMblind){ m_hEMblind = hEMblind;}
	void setMEBlind(TH1D* hMEblind){ m_hMEblind = hMEblind;}
//	void setSigHTM(TH1D* hsigHTM){m_hsigHTM = hsigHTM;}
//	void setSigHTE(TH1D* hsigHTE){m_hsigHTE = hsigHTE;}
//	void setSigHTM_EM(TH1D* hsigHTMEM){m_hsigHTM_EM = hsigHTMEM;}
//	void setSigHTE_ME(TH1D* hsigHTEME){m_hsigHTE_ME = hsigHTEME;}
	void setSig(TH1D* hsig){ m_hsig = hsig;}
	void setEM(TH1D* hEM1, TH1D* hEM2);
	void setME(TH1D* hME1, TH1D* hME2);
	void setSig(TH1D* hs1, TH1D* hs2);
	void setEM2(TH1D* hEM2){ m_SR2_hEM = hEM2;}
	void setME2(TH1D* hME2){ m_SR2_hME = hME2;}
	void setN(){m_n = new double[m_nbins]; for (int i=0; i<m_nbins; i++){m_n[i]=m_hEM->GetBinContent(i+1);}}
	void setM(){m_m = new double[m_nbins]; for (int i=0; i<m_nbins; i++){m_m[i]=m_hME->GetBinContent(i+1);}}
	void setS(){m_S = new double[m_nbins]; for (int i=0; i<m_nbins; i++){m_S[i]=m_hsig->GetBinContent(i+1);}}
	void setN(double* N){m_n = new double[m_nbins]; memcpy(m_n,N,(m_nbins)*sizeof(double));}
	void setM(double* M){m_m = new double[m_nbins]; memcpy(m_m,M,(m_nbins)*sizeof(double));}
	void setS(double* S){m_S = new double[m_nbins]; memcpy(m_S,S,(m_nbins)*sizeof(double));}
	void setBins(){
		m_Bins = new double[m_nbins+1]; const TArrayD* array = m_hEM->GetXaxis()->GetXbins();
		const Double_t *bins = array->GetArray();
		memcpy(m_Bins,bins,(m_nbins/m_numSR+1)*sizeof(double));
		if(m_numSR==2){double bins2[m_nbins/m_numSR+1];
			for (int i=0; i<m_nbins/m_numSR; i++){
				bins2[i] = bins[m_nbins/m_numSR]+bins[i+1];
			}
			memcpy(m_Bins+m_nbins/m_numSR+1,bins2,(m_nbins/m_numSR)*sizeof(double));
		}
	}
	void setBins(double* bins){m_Bins = new double[m_nbins+1]; memcpy(m_Bins,bins,(m_nbins+1)*sizeof(double));}
	void setBlindHistos();
	//do the minimization to calculate muHat.
	//  Called from all constructors besides the default and the copy constructor
	void setMuHat();
	//copy muHat. the copy constructor uses this
	void setMuHat(double muHat);

	void DrawEMME(){ m_hMEblind->Draw("E1"); m_hEMblind->Draw("E1 sames");}
	void DrawL3Graph();
	void DrawL2Graph();

	void appendSignalME(double strength);

	void free();

private:
	int m_numSR;
	TH1D* m_hEM;
	TH1D* m_hME;
	//TODO: make the object hold two signal histos and use them accordingly
	TH1D* m_hsig;
//	TH1D* m_hsigHTM;
//	TH1D* m_hsigHTE;
//	TH1D* m_hsigHTM_EM;
//	TH1D* m_hsigHTE_ME;
	TH1D* m_hEMblind;
	TH1D* m_hMEblind;
	TH1D* m_SR2_hEM;
	TH1D* m_SR2_hME;
	int m_nbins;
	int m_minBin;
	int m_maxBin;
	double* m_Bins;
	double* m_n;
	double* m_m;
	double* m_S;
	double* m_SR2_n;
	double* m_SR2_m;
	double* m_SR2_S;
	double m_muHat;


};






#endif /* DATA_H_ */
