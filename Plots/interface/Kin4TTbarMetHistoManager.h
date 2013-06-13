#ifndef Kin4TTbarMetHistoManager_h
#define Kin4TTbarMetHistoManager_h

#include "NTFormat/interface/NTEvent.h"
#include "Selection/interface/TTbarMetSelection.h"
#include "Plots/interface/HistoManager.h"
#include "TMath.h"



using namespace IPHCTree;

class Kin4TTbarMetHistoManager: public HistoManager{

  public:
	Kin4TTbarMetHistoManager();
	~Kin4TTbarMetHistoManager();

	//Initialisation methods

	void CreateHistos(); /** Create a bunch of standard histos */

	//Fill methods

	void Fill(const TTbarMetSelection& sel, const int& maxSelStep, const int& iChannel, const int& iDataset, const float& weight);
	void FillSelStep(const TTbarMetSelection& sel, const int& iSelStep, const int& iChannel, const int& iDataset, const float& weight);


};

#endif
