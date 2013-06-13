#ifndef TTbarMetHistoManager_h
#define TTbarMetHistoManager_h

#include "NTFormat/interface/NTEvent.h"

#include "Selection/interface/TTbarMetSelection.h"
#include "Tools/interface/Dataset.h"
#include "ElectronHistoManager.h"
#include "MuonHistoManager.h"
#include "JetHistoManager.h"
#include "BJetHistoManager.h"
#include "METHistoManager.h"
#include "NofEventsHistoManager.h"
#include "Kin4TTbarMetHistoManager.h"

#include <TH1F.h>

using namespace IPHCTree;

class TTbarMetHistoManager{

  public:
	TTbarMetHistoManager();
	~TTbarMetHistoManager();

 	void LoadDatasets(vector<Dataset> datasets);    
 	void LoadSelectionSteps(vector<string> selectionSteps);
 	void LoadChannels(vector<string> channels);

        void Clear();	
  	void CreateHistos();	

	void Fill(const TTbarMetSelection& sel, NTEvent* event, const vector<NTMuon>& candMuon, const vector<NTElectron>& candElec, const int& maxSelStep, const int& iChannel, const int& iDataset, const float& weight);
        void Fill(const TTbarMetSelection& sel, const int& maxSelStep, const int& iChannel, const int& iDataset, const float& weight);
	//void FillwBweight(const DiLeptonSelection& sel, NTEvent* event, const vector<NTMuon>& candMuon, const vector<NTElectron>& candElec, const int& maxSelStep, const int& iChannel, const int& iDataset, const float& weight, const vector<float>& weightb);
	
	void Compute();
	void ComputeForProof();
	void Write(TFile* file);


  private:

	ElectronHistoManager elecHistos;
	MuonHistoManager muHistos;
	JetHistoManager jetHistos;
	BJetHistoManager bjetHistos;
	METHistoManager metHistos;
	NofEventsHistoManager nofEvtHistos;
	Kin4TTbarMetHistoManager  kinHistos;
};

#endif
