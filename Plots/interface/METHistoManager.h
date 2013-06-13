#ifndef METHistoManager_h
#define METHistoManager_h

#include "NTFormat/interface/NTMET.h"

#include "Plots/interface/HistoManager.h"


using namespace IPHCTree;

class METHistoManager: public HistoManager{

  public:
	METHistoManager();
	~METHistoManager();

	//Initialisation methods

	void CreateHistos(); /** Create a bunch of standard histos */

	//Fill methods
	void Fill(NTMET met, const int& maxSelStep, const int& iChannel, const int& iDataset, const float& weight);
	void FillSelStep(NTMET met, const int& iSelStep, const int& iChannel, const int& iDataset, const float& weight);


};

#endif
