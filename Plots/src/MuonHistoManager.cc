#include "Plots/interface/MuonHistoManager.h"




MuonHistoManager::MuonHistoManager(){
}

MuonHistoManager::~MuonHistoManager(){
}


void MuonHistoManager::CreateHistos(){
	AddHisto(string("Multiplicity"),string("#(muons)"),string("Nof muons"),10,0,10);	
	AddHisto(string("Pt"),string("Pt(muon)"),string("p_{T}(muon)"),50,0,100);	
	AddHisto(string("Eta"),string("Eta(muon)"),string("#eta(muon)"),60,-3,3);	
	AddHisto(string("Phi"),string("Phi(muon)"),string("#phi(muon)"),64,-3.2,3.2);	
	AddHisto(string("d0"),string("d0(muon)"),string("d0(muon)"),50,0,2);	
	AddHisto(string("charge"),string("charge(muon)"),string("charge(muon)"),5,-2,2);	
	AddHisto(string("trackIso"),string("trackIso(muon)"),string("trackIso(muon)"),50,0,10);	
	AddHisto(string("caloIso"),string("caloIso(muon)"),string("caloIso(muon)"),50,0,10);	
	AddHisto(string("hcalIso"),string("hcalIso(muon)"),string("hcalIso(muon)"),50,0,10);	
	AddHisto(string("relIso"),string("relIso(muon)"),string("relIso(muon)"),50,0,10);	
	AddHisto(string("Pt1"),string("Pt(muon1)"),string("p_{T}(muon1)"),50,0,100);	
	AddHisto(string("Eta1"),string("Eta(muon1)"),string("#eta(muon1)"),60,-3,3);	
	AddHisto(string("Pt2"),string("Pt(muon2)"),string("p_{T}(muon2)"),50,0,100);	
	AddHisto(string("Eta2"),string("Eta(muon2)"),string("#eta(muon2)"),60,-3,3);	
}

void MuonHistoManager::Fill(const vector<NTMuon>& muons, const int& maxSelStep, const int& iChannel, const int& iDataset, const float& weight){
	if(!Check(iChannel, iDataset)) return;
	for(unsigned int i=0;i<SelectionSteps.size();i++){
		if(maxSelStep>=(int)i) FillSelStep(muons, i , iChannel, iDataset, weight);	
	}	
}

void MuonHistoManager::FillSelStep(const vector<NTMuon>& muons, const int& iSelStep, const int& iChannel, const int& iDataset, const float& weight){
        if(!Check(iChannel, iSelStep, iDataset, 0) ) return;
	Histos[0][iChannel][iSelStep][iDataset].Fill(muons.size(),weight);
	for(int i=0;i<(int) muons.size();i++){
		//very important:
		//respect the order of the function CreateHistos to fill the histograms
		Histos[1][iChannel][iSelStep][iDataset].Fill(muons[i].p4.Pt(),weight);
		Histos[2][iChannel][iSelStep][iDataset].Fill(muons[i].p4.Eta(),weight);
		Histos[3][iChannel][iSelStep][iDataset].Fill(muons[i].p4.Phi(),weight);
		Histos[4][iChannel][iSelStep][iDataset].Fill(muons[i].D0,weight);
		Histos[5][iChannel][iSelStep][iDataset].Fill(muons[i].charge,weight);
		Histos[6][iChannel][iSelStep][iDataset].Fill(muons[i].isolation["Trk03"],weight);
		Histos[7][iChannel][iSelStep][iDataset].Fill(muons[i].isolation["ECalo03"],weight);
		Histos[8][iChannel][iSelStep][iDataset].Fill(muons[i].isolation["HCalo03"],weight);
		//Histos[[9]iChannel][iSelStep][iDataset].Fill(muons[i].RelIso03(),weight);
		NTMuon mu = muons[i];
		float relIso = mu.RelIso03();
		Histos[9][iChannel][iSelStep][iDataset].Fill(relIso,weight);
	}
	if(muons.size()>0){
		Histos[10][iChannel][iSelStep][iDataset].Fill(muons[0].p4.Pt(),weight);
		Histos[11][iChannel][iSelStep][iDataset].Fill(muons[0].p4.Eta(),weight);
	}
	if(muons.size()>1){
		Histos[12][iChannel][iSelStep][iDataset].Fill(muons[1].p4.Pt(),weight);
		Histos[13][iChannel][iSelStep][iDataset].Fill(muons[1].p4.Eta(),weight);
	}
}


