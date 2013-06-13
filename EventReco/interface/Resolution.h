#ifndef Resolution_h
#define Resolution_h

// IPHC headers
#include "NTFormat/interface/NTEvent.h"
#include "NTFormat/interface/NTGenParticle.h"

// ROOT headers
#include <TTree.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TMath.h>
#include <TF1.h>
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnSimplex.h"
#include "Minuit2/MnScan.h"
#include "Minuit2/FunctionMinimum.h"
#include "TFitterMinuit.h"

// STL headers
#include <iostream>
#include <memory>
#include <vector>
#include <assert.h>
#include <math.h>

class Resolution
{
	
  // -------------------------------------------------------------
  //                       method members
  // -------------------------------------------------------------
 
 public :

	//Constructor
	Resolution();

	//Destructor
	~Resolution();

	double GetSigmaPtJet(const TLorentzVector& jet);

	std::vector<float> GetVecSigmaPtJets(const std::vector<IPHCTree::NTJet>& jets);
 
	double GetChi2(const std::vector<IPHCTree::NTJet>& jets);

	//double GetSelChi2(const TTbarMetSelection* sel);

	static double fc2 (double c1, double m12, double m22, double m02);

	static double fchi2 (double c1, double pt1, double sigma1, double pt2,
								double sigma2, double m12, double m22, double m02);

	static void minuitFunction(int&, double* , double &result, double* par, int);

	//!\Accessors require execution of GetChi2 Method !

	//Accessor to Sigma reconstructed top
	const double GetSigmaRecoTop()
  {return smtop_;}

	//Accessor Mass reconstructed top
	const double GetMRecoTop()
	{return mtop_;}	
	
	//Accessor TLV reconstructed top
	const TLorentzVector GetTLVRecoTop()
	{return hadT_;}
	

 private :
	
	std::vector<float> vecSigmaPtJets;
	
	//methode chi2
	vector<float> btag;
	vector<float> sigma_jets;	
	int n_jets;
	vector<int> v_i, v_j;
	vector<double> v_k1, v_k2;
	TLorentzVector hadW;
	double p1;
	double c1;
	double c2;
	int n_btag;
	double chi2min;
	double pt_b;
	int k, l;
	int nwb;
	double pt_w1, pt_w2;
	double massW;
	TLorentzVector hadT;
	double massT;
	double pt_w;
	double sigma_w2;
	double smw2;
	double pt_t;
	double sigma_t2;
	double smtop2;
	double c_chi2;
	TLorentzVector TLVhadT;
	double smtop_;
	double mtop_;
	TLorentzVector hadT_;

//    TTbarMetSelection* sel;
	IPHCTree::NTMonteCarlo mcInfo;
	std::vector<IPHCTree::NTGenParticle> MCParticles;
	std::vector<int> TopId;
	std::vector<TLorentzVector> TopHad;
	double Chi2topHad;
	double Chi2;

	TF1 sigma;
	TF1 sigma_void;
	TFitterMinuit minimizer;
	TFitterMinuit minimizer_void;

};
#endif
