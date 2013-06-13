#ifndef DiLepRecoHistoManager_h
#define DiLepRecoHistoManager_h

#include "NTFormat/interface/NTEvent.h"

#include "Plots/interface/HistoManager.h"


using namespace IPHCTree;

class DiLepRecoHistoManager: public HistoManager{

  public:
	DiLepRecoHistoManager();
	~DiLepRecoHistoManager();

	//Initialisation methods

	void CreateHistos(); /** Create a bunch of standard histos */

	//Fill methods
	void Fill(NTEvent* event, const vector<NTMuon>& candMuon, const vector<NTElectron>& candElec, const int& maxSelStep, const int& iChannel, const int& iDataset, const float& weight);
	void FillSelStep(NTEvent* event, const vector<NTMuon>& candMuon, const vector<NTElectron>& candElec, const int& iSelStep, const int& iChannel, const int& iDataset, const float& weight);


};

#endif
