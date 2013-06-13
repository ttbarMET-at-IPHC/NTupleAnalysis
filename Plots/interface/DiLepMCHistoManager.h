#ifndef DiLepMCHistoManager_h
#define DiLepMCHistoManager_h

#include "NTFormat/interface/NTEvent.h"
#include "Plots/interface/HistoManager.h"


using namespace IPHCTree;

class DiLepMCHistoManager: public HistoManager
{

  public:
	DiLepMCHistoManager();
	~DiLepMCHistoManager();

	//Initialisation methods

	void CreateHistos(); /** Create a bunch of standard histos */

	//Fill methods
	void Fill(NTEvent* event, const int& maxSelStep, const int& iChannel, const int& iDataset, const float& weight);
	void FillSelStep(NTEvent* event, const int& iSelStep, const int& iChannel, const int& iDataset, const float& weight);


};

#endif
