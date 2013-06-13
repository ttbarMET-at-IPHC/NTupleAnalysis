#include "../interface/NofEventsHistoManager.h"
#include <sstream>


NofEventsHistoManager::NofEventsHistoManager(){
}

NofEventsHistoManager::~NofEventsHistoManager(){
}


void NofEventsHistoManager::CreateHistos(){
	AddHisto(string("NofEvents"),string("#(Events)"),string("Nof Events"),SelectionSteps.size(),0,SelectionSteps.size());	
}

void NofEventsHistoManager::Fill(const int& maxSelStep, const int& iChannel, const int& iDataset, const float& weight){
	if(!Check(iChannel, iDataset)) return;
	for(unsigned int i=0;i<SelectionSteps.size();i++){
		if(maxSelStep>=(int)i) Histos[0][iChannel][0][iDataset].Fill(i, weight); // Fill 1 histo for all selection steps
	}	
}
