#ifndef NofEventsHistoManager_h
#define NofEventsHistoManager_h

#include "HistoManager.h"

class NofEventsHistoManager: public HistoManager{
  public:
	NofEventsHistoManager();
	~NofEventsHistoManager();

	//Initialisation method
	void CreateHistos(); /** Create a bunch of standard histos */

        //Fill method
	void Fill(const int& maxSelStep, const int& iChannel, const int& iDataset, const float& weight);
	
};

#endif
