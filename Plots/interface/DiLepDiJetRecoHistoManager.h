#ifndef DiLepDiJetRecoHistoManager_h
#define DiLepDiJetRecoHistoManager_h

#include "NTFormat/interface/NTEvent.h"

#include "Plots/interface/HistoManager.h"


using namespace IPHCTree;

class DiLepDiJetRecoHistoManager: public HistoManager{

  public:
	DiLepDiJetRecoHistoManager();
	~DiLepDiJetRecoHistoManager();

	//Initialisation methods

	void CreateHistos(); /** Create a bunch of standard histos */

	//Fill methods
	void Fill(NTEvent* event, const vector<NTMuon>& candMuon, const vector<NTElectron>& candElec, const vector<NTJet>& candJet, const int& maxSelStep, const int& iChannel, const int& iDataset, const float& weight);
	void FillSelStep(NTEvent* event, const vector<NTMuon>& candMuon, const vector<NTElectron>& candElec, const vector<NTJet>& candJet, const int& iSelStep, const int& iChannel, const int& iDataset, const float& weight);


};

#endif
