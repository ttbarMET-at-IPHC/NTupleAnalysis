#ifndef StopAnaReco_h
#define StopAnaReco_h

// IPHC headers
#include "NTFormat/interface/NTEvent.h"

// ROOT headers
#include <TTree.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TH1.h>
#include <TGraph.h>
#include <TMath.h>

// STL headers
#include <iostream>
#include <memory>
#include <vector>
#include <stdlib.h>

class StopAnaReco
{

  // -------------------------------------------------------------
  //                       method members
  // -------------------------------------------------------------

 public :

	//Constructor
	StopAnaReco();

	StopAnaReco(TH1F* HS,TH1F* HB);

	StopAnaReco(const IPHCTree::NTMET& NTMET,
						  const std::vector<IPHCTree::NTJet>& NTJet,
						  const std::vector<IPHCTree::NTElectron>& NTElectron,
						  const std::vector<IPHCTree::NTMuon>& NTMuon);

	//Destructor
	~StopAnaReco();


	bool IsSSM(const IPHCTree::NTMET& MET,const TLorentzVector& obj) const;

	pair<float,float> HT() const;

	float GetHTSSM() const;

	float GetHTOSM() const;

	float GetHT() const;

	float GetMt();

	//Normalisation de l'histogramme
	void Normalize(TH1F* h1);

	//Overlap de deux histogrammes
	double GetOverlap(TH1F* HS,TH1F* HB);

	//void ComputeEff(TH1F* HistoS,TH1F* HistoB);

	//std::vector<TGraph*> GetComputeEff(TH1F* HistoS,TH1F* HistoB);
	std::vector<TGraph*> ComputeEff(TH1F* HistoS,TH1F* HistoB);

	//Calcul du max de Fsigma
	double GetMaxCutEff(TH1F* HistoS,TH1F* HistoB);

	double GetK(const double& Sigma,
							const double& Luminosity,
							const int& NofEvents) const;

	TH1F* GetHistoAdd(std::vector<TH1F*>& HistoVect);

	void NormalizeKfactor(TH1F* h1,const double Kfactor);

	TH1F* GetHistoNorm(const std::vector<TH1F*>& HistoVect,
										 const std::vector<Double_t>& KVect);

	double GetDeltaPhi();

  // -------------------------------------------------------------
  //                        data members
  // -------------------------------------------------------------

 private :

	TLorentzVector obj;

	std::vector<IPHCTree::NTJet> NTJet;
	std::vector<IPHCTree::NTElectron> NTElectron;
	std::vector<IPHCTree::NTMuon> NTMuon;
	IPHCTree::NTMET NTMET;

	std::vector<IPHCTree::NTJet> Jet;
	std::vector<IPHCTree::NTElectron> Electron;
	std::vector<IPHCTree::NTMuon> Muon;
	IPHCTree::NTMET MET;

	float HTSSM;
	float HTOSM;

	std::vector<TH1F*> HistoVect;
	std::vector<Double_t> KVect;
	TH1F* KHisto;
	TH1F* mHistoVect;

	float Mt;

	TH1F* HS;
	TH1F* HB;
	TH1F* Histo;

	TH1F* HistoS;
	TH1F* HistoB;
	float overLap;

	TGraph* Cut;
	TGraph* Rap;
	std::vector<TGraph*> GraphVect;
	int max;
	double* x;
	double* Fsig;
	double* EpsS;
	double* EpsB;

	TGraph Graph;
	float m;

	double Sigma;
	double Luminosity;
	double NofEvents;

	double Kfactor;

	TH1F* h1;
	TH1F* h2;
	TH1F* h3;

	int i;
	int j;
	int k;
};

#endif
