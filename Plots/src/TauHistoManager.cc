#include "Plots/interface/TauHistoManager.h"




TauHistoManager::TauHistoManager(){
}

TauHistoManager::~TauHistoManager(){
}


void TauHistoManager::CreateHistos(){
	AddHisto(string("taulmutiplicity"),string("#(taus)"),string("Nof taus"),10,0,10);	
	AddHisto(string("Pt"),string("Pt(tau)"),string("p_{T}(tau)"),50,0,100);	
	AddHisto(string("Eta"),string("Eta(tau)"),string("#eta(tau)"),60,-3,3);	
	AddHisto(string("Phi"),string("Phi(tau)"),string("#phi(tau)"),64,-3.2,3.2);	
	AddHisto(string("d0"),string("d0(tau)"),string("d0(tau)"),50,0,2);	
	AddHisto(string("dz"),string("dz(tau,vtx)"),string("dz(tau,vtx)"),100,-20,20);	
	AddHisto(string("charge"),string("charge(tau)"),string("charge(tau)"),5,-2,2);
	AddHisto(string("leadTrackPt"),string("leadTrackPt(tau)"),string("leadTrackPt(tau)"),50,0,100);
	AddHisto(string("byLooseIsolation"),string("byLooseIsolation(tau)"),string("byLooseIsolation(tau)"),50,-1,1);
	AddHisto(string("byMediumIsolation"),string("byMediumIsolation(tau)"),string("byMediumIsolation(tau)"),50,-1,1);
	AddHisto(string("byTightIsolation"),string("byTightIsolation(tau)"),string("byTightIsolation(tau)"),50,-1,1);	
	AddHisto(string("numSigConeTracks"),string("numSigConeTracks(tau)"),string("numSigConeTracks(tau)"),20,0,20);
	AddHisto(string("numIsoConeTracks"),string("numIsoConeTracks(tau)"),string("numIsoConeTracks(tau)"),20,0,20);
        AddHisto(string("isolationPFChargedHadrCandsPtSum"),string("isolationPFChargedHadrCandsPtSum(tau)"),string("isolationPFChargedHadrCandsPtSum(tau)"),10,0,1.);
	AddHisto(string("isolationPFGammaCandsEtSum"),string("isolationPFGammaCandsEtSum(tau)"),string("isolationPFGammaCandsEtSum(tau)"),10,0,1.);
	AddHisto(string("emFraction"),string("emFraction(tau)"),string("emFraction(tau)"),20,0,1.);
	AddHisto(string("againstElectronLoose"),string("againstElectronLoose(tau)"),string("againstElectronLoose(tau)"),20,0,1.);
        AddHisto(string("againstMuonLoose"),string("againstMuonLoose(tau)"),string("againstMuonLoose(tau)"),20,0,1.);
        AddHisto(string("q*Eta"),string("q*Eta(tau)"),string("#eta(tau)*q"),20,-3,3);
	AddHisto(string("MT"),string("MT"),string("MT(W)"),40,0,160);
}

void TauHistoManager::Fill(NTMET met, const vector<NTTau>& taus, const vector<NTVertex>& vertices, const int& maxSelStep, const int& iChannel, const int& iDataset, const float& weight){
	if(!Check(iChannel, iDataset)) return;
	for(unsigned int i=0;i<SelectionSteps.size();i++){
		if(maxSelStep>=(int)i) FillSelStep(met, taus, vertices, i, iChannel, iDataset, weight);	
	}	
}

void TauHistoManager::FillSelStep(NTMET met, const vector<NTTau>& taus, const vector<NTVertex>& vertices, const int& iSelStep, const int& iChannel, const int& iDataset, const float& weight){
        if(!Check(iChannel, iSelStep, iDataset, 0) ) return;
	Histos[0][iChannel][iSelStep][iDataset].Fill(taus.size(),weight);
	float MTW = -9999;
	for(int i=0;i<(int) taus.size();i++){
	        NTTau taui = taus[i];
          MTW =sqrt (pow(taui.p4.Et() + met.met(),2) - pow(taui.p4.Px() + met.p2.Px(),2) - pow(taui.p4.Py() + met.p2.Py(),2));
	        //very important:
		//respect the order of the function CreateHistos to fill the histograms
		Histos[1][iChannel][iSelStep][iDataset].Fill(taus[i].p4.Pt(),weight);
		Histos[2][iChannel][iSelStep][iDataset].Fill(taus[i].p4.Eta(),weight);
		Histos[3][iChannel][iSelStep][iDataset].Fill(taus[i].p4.Phi(),weight);
		Histos[4][iChannel][iSelStep][iDataset].Fill(taus[i].D0,weight);
		Histos[5][iChannel][iSelStep][iDataset].Fill(taus[i].vertex.Z()-vertices[0].p3.Z(),weight);
		Histos[6][iChannel][iSelStep][iDataset].Fill(taus[i].charge,weight);
		Histos[7][iChannel][iSelStep][iDataset].Fill(taus[i].leadTrackPt,weight);
		Histos[8][iChannel][iSelStep][iDataset].Fill(taui.GetDiscriminator("byLooseIsolation"),weight);
		Histos[9][iChannel][iSelStep][iDataset].Fill(taui.GetDiscriminator("byMediumIsolation"),weight);
		Histos[10][iChannel][iSelStep][iDataset].Fill(taui.GetDiscriminator("byTightIsolation"),weight);		
		Histos[11][iChannel][iSelStep][iDataset].Fill(taus[i].numSigConeTracks,weight);
		Histos[12][iChannel][iSelStep][iDataset].Fill(taus[i].numIsoConeTracks,weight);
		Histos[13][iChannel][iSelStep][iDataset].Fill(taus[i].isolationPFChargedHadrCandsPtSum,weight);
		Histos[14][iChannel][iSelStep][iDataset].Fill(taus[i].isolationPFGammaCandsEtSum,weight);
		Histos[15][iChannel][iSelStep][iDataset].Fill(taus[i].emFraction,weight);
		Histos[16][iChannel][iSelStep][iDataset].Fill(taui.GetDiscriminator("againstElectronLoose"),weight);
		Histos[17][iChannel][iSelStep][iDataset].Fill(taui.GetDiscriminator("againstMuonLoose"),weight);
	        Histos[18][iChannel][iSelStep][iDataset].Fill((taus[i].p4.Eta())*(taus[i].charge),weight);
		Histos[19][iChannel][iSelStep][iDataset].Fill(MTW,weight);}
}


