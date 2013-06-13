#include "Plots/interface/EventShapesHistoManager.h"




EventShapesHistoManager::EventShapesHistoManager(){
}

EventShapesHistoManager::~EventShapesHistoManager(){
}


void EventShapesHistoManager::CreateHistos(){
	AddHisto(string("Aplanarity"),string("Aplanarity"),string("Aplanarity"),100,0,1.);
	AddHisto(string("Sphericity"),string("Sphericity"),string("Sphericity"),100,0,1);
	AddHisto(string("Circularity"),string("Circularity"),string("Circularity"),100,0,1);
	AddHisto(string("Isotropy"),string("Isotropy"),string("Isotropy"),100,0,1);
	AddHisto(string("C"),string("C"),string("C"),100,0,1);
	AddHisto(string("D"),string("D"),string("D"),100,0,1);
	AddHisto(string("Ht"),string("HT"),string("H_{T}"),200,0,1000);
	AddHisto(string("H"),string("H"),string("H"),200,0,1000);
	AddHisto(string("sqrts"),string("sqrt(s)"),string("sqrt(s)"),200,0,1000);
	AddHisto(string("NNOutput"),string("NNOutput"),string("NNOutput"), 100, -1.25, 1.5 );
	AddHisto(string("M3"),string("M3"),string("M3"), 500, 0, 500 );
}

void EventShapesHistoManager::Fill(const vector<NTJet>& candJet, const int& maxSelStep, const int& iChannel, const int& iDataset, const float& weight, const float& NNOutput){
	if(!Check(iChannel, iDataset)) return;
	for(unsigned int i=0;i<SelectionSteps.size();i++){
	  if(maxSelStep>=(int)i) FillSelStep(candJet, i, iChannel, iDataset, weight, NNOutput);	
	}	
}

void EventShapesHistoManager::FillSelStep(const vector<NTJet>& candJet, const int& iSelStep, const int& iChannel, const int& iDataset, const float& weight, const float& NNOutput){
        if(!Check(iChannel, iSelStep, iDataset, 0) ) return;
	EventShapes evShape(candJet);
	Histos[0][iChannel][iSelStep][iDataset].Fill(evShape.aplanarity(),weight);
	Histos[1][iChannel][iSelStep][iDataset].Fill(evShape.sphericity(),weight);
	Histos[2][iChannel][iSelStep][iDataset].Fill(evShape.circularity(),weight);
	Histos[3][iChannel][iSelStep][iDataset].Fill(evShape.isotropy(),weight);
	Histos[4][iChannel][iSelStep][iDataset].Fill(evShape.C(),weight);
	Histos[5][iChannel][iSelStep][iDataset].Fill(evShape.D(),weight);
	Histos[6][iChannel][iSelStep][iDataset].Fill(evShape.HT(),weight);
	Histos[7][iChannel][iSelStep][iDataset].Fill(evShape.H(),weight);
	Histos[8][iChannel][iSelStep][iDataset].Fill(evShape.sqrt_s(),weight);
	Histos[9][iChannel][iSelStep][iDataset].Fill(NNOutput,weight);
	Histos[10][iChannel][iSelStep][iDataset].Fill(evShape.M3(candJet),weight);
}
