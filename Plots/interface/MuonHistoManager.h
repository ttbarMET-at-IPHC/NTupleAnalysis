#ifndef MuonHistoManager_h
#define MuonHistoManager_h

#include "NTFormat/interface/NTMuon.h"
#include "Plots/interface/HistoManager.h"


using namespace IPHCTree;

class MuonHistoManager: public HistoManager{

  public:
	MuonHistoManager();
	~MuonHistoManager();

	//Initialisation methods

	void CreateHistos(); /** Create a bunch of standard histos */

	//Fill methods
	void Fill(const vector<NTMuon>& muons, const int& maxSelStep, const int& iChannel, const int& iDataset, const float& weight);
	void FillSelStep(const vector<NTMuon>& muons, const int& iSelStep, const int& iChannel, const int& iDataset, const float& weight);


};

#endif
