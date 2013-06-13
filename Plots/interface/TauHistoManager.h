#ifndef TauHistoManager_h
#define TauHistoManager_h

#include "NTFormat/interface/NTTau.h"
#include "NTFormat/interface/NTVertex.h"
#include "NTFormat/interface/NTEvent.h"
#include "Plots/interface/HistoManager.h"

using namespace IPHCTree;

class TauHistoManager: public HistoManager{

  public:
	TauHistoManager();
	~TauHistoManager();

	//Initialisation methods

	void CreateHistos(); /** Create a bunch of standard histos */

	//Fill methods
	void Fill(NTMET met, const vector<NTTau>& taus, const vector<NTVertex>& vertices, const int& maxSelStep, const int& iChannel, const int& iDataset, const float& weight);
	void FillSelStep(NTMET met, const vector<NTTau>& taus, const vector<NTVertex>& vertices,const int& iSelStep, const int& iChannel, const int& iDataset, const float& weight);


};

#endif
