#include "Plots/interface/DiLepDiJetRecoHistoManager.h"
#include "Selection/interface/SSDiLeptonSelection.h"



DiLepDiJetRecoHistoManager::DiLepDiJetRecoHistoManager(){
}

DiLepDiJetRecoHistoManager::~DiLepDiJetRecoHistoManager(){
}


void DiLepDiJetRecoHistoManager::CreateHistos(){
	AddHisto(string("DiLepMass"),string("DiLepMass"),string("Di-lepton mass [GeV]"),100,0,800);	
	AddHisto(string("DiLepJetMass"),string("DiLepJetMass"),string("Di-lepton + jet mass [GeV]"),100,0,800);	
	AddHisto(string("DiLepDiJetMass"),string("DiLepDiJetMass"),string("Di-lepton + di-jet mass [GeV]"),100,0,800);	
}

void DiLepDiJetRecoHistoManager::Fill(NTEvent* event, const vector<NTMuon>& candMuon, const vector<NTElectron>& candElec, const vector<NTJet>& candJet, const int& maxSelStep, const int& iChannel, const int& iDataset, const float& weight){
	if(!Check(iChannel, iDataset)) return;
	for(unsigned int i=0;i<SelectionSteps.size();i++){
		if(maxSelStep>=(int)i) FillSelStep(event, candMuon, candElec, candJet, i , iChannel, iDataset, weight);	
	}	
}

void DiLepDiJetRecoHistoManager::FillSelStep(NTEvent* event, const vector<NTMuon>& candMuon, const vector<NTElectron>& candElec, const vector<NTJet>& candJet, const int& iSelStep, const int& iChannel, const int& iDataset, const float& weight){
        if(!Check(iChannel, iSelStep, iDataset, 0) ) return;
	SSDiLeptonSelection sel;
	float mass = sel.DiLeptonMass(candMuon,candElec);
	/*float mt = sel.DiLeptonMass(candMuon,candElec);
   	float DPhi = -999;
	bool isSameSign = false;
	if(candMuon.size()==2) {
		DPhi = candMuon[0].p4.Phi()-candMuon[1].p4.Phi();	
		if(candMuon[0].charge == candMuon[1].charge) isSameSign = true;
	}
	if(candElec.size()==2){
		DPhi = candElec[0].p4.Phi()-candElec[1].p4.Phi();	
		if(candElec[0].charge == candElec[1].charge) isSameSign = true;
	}
	if(candMuon.size()==1 && candElec.size()==1){
		DPhi = candMuon[0].p4.Phi()-candElec[0].p4.Phi();
		if(candMuon[0].charge == candElec[0].charge) isSameSign = true;
	}*/
	TLorentzVector DiLepton = sel.DiLeptonCand(candMuon,candElec);
	if(DiLepton.E()!=0) Histos[0][iChannel][iSelStep][iDataset].Fill(mass,weight);
	if(candJet.size()>0 && DiLepton.E()!=0){
		TLorentzVector llj;
		llj = DiLepton + candJet[0].p4;
		Histos[1][iChannel][iSelStep][iDataset].Fill(llj.M(),weight);
	}
	if(candJet.size()>1 && DiLepton.E()!=0){
		TLorentzVector lljj;
		lljj = DiLepton + candJet[0].p4 + candJet[1].p4;
		Histos[2][iChannel][iSelStep][iDataset].Fill(lljj.M(),weight);
	}
	
}


