#include "Plots/interface/SSDiLepAnaHistoManager.h"



SSDiLepAnaHistoManager::SSDiLepAnaHistoManager(){
}

SSDiLepAnaHistoManager::~SSDiLepAnaHistoManager(){
}

void SSDiLepAnaHistoManager::LoadDatasets(vector<Dataset> datasets){
	elecHistos.LoadDatasets(datasets);
	muHistos.LoadDatasets(datasets);
	jetHistos.LoadDatasets(datasets);
	metHistos.LoadDatasets(datasets);
	recoHistos.LoadDatasets(datasets);
	mcHistos.LoadDatasets(datasets);
	bjetHistos.LoadDatasets(datasets);
	lljjrecoHistos.LoadDatasets(datasets);
	nofEvtHistos.LoadDatasets(datasets);
}    
 	

void SSDiLepAnaHistoManager::LoadSelectionSteps(vector<string> selectionSteps){
	elecHistos.LoadSelectionSteps(selectionSteps);
	muHistos.LoadSelectionSteps(selectionSteps);
	jetHistos.LoadSelectionSteps(selectionSteps);
	metHistos.LoadSelectionSteps(selectionSteps);
	recoHistos.LoadSelectionSteps(selectionSteps);
	mcHistos.LoadSelectionSteps(selectionSteps);
	bjetHistos.LoadSelectionSteps(selectionSteps);
	lljjrecoHistos.LoadSelectionSteps(selectionSteps);
	nofEvtHistos.LoadSelectionSteps(selectionSteps);
}
	
 
void SSDiLepAnaHistoManager::LoadChannels(vector<string> channels){
	elecHistos.LoadChannels(channels);
	muHistos.LoadChannels(channels);
	jetHistos.LoadChannels(channels);
	metHistos.LoadChannels(channels);
	recoHistos.LoadChannels(channels);
	mcHistos.LoadChannels(channels);
	bjetHistos.LoadChannels(channels);
	lljjrecoHistos.LoadChannels(channels);
	nofEvtHistos.LoadChannels(channels);
}

void SSDiLepAnaHistoManager::Clear(){
	elecHistos.Clear();
	muHistos.Clear();
	jetHistos.Clear();
	metHistos.Clear();
	recoHistos.Clear();
	mcHistos.Clear();
	bjetHistos.Clear();
	lljjrecoHistos.Clear();
	nofEvtHistos.Clear();
}


void SSDiLepAnaHistoManager::CreateHistos(){
	elecHistos.CreateHistos();
	muHistos.CreateHistos();
	jetHistos.CreateHistos();
	metHistos.CreateHistos();
	recoHistos.CreateHistos();
	mcHistos.CreateHistos();
	bjetHistos.CreateHistos();
	lljjrecoHistos.CreateHistos();
	nofEvtHistos.CreateHistos();
}

void SSDiLepAnaHistoManager::Fill(const SSDiLeptonSelection& sel, NTEvent* event, const vector<NTMuon>& candMuon, const vector<NTElectron>& candElec, const int& maxSelStep, const int& iChannel, const int& iDataset, const float& weight){
        elecHistos.Fill(candElec, maxSelStep, iChannel, iDataset, weight);
        muHistos.Fill(candMuon, maxSelStep, iChannel, iDataset, weight);
        jetHistos.Fill(sel.GetJetsForAna(), maxSelStep, iChannel, iDataset, weight);
        metHistos.Fill(sel.GetMET(), maxSelStep, iChannel, iDataset, weight);
        recoHistos.Fill(event, candMuon, candElec, maxSelStep, iChannel, iDataset, weight);
        mcHistos.Fill(event, maxSelStep, iChannel, iDataset, weight);
        bjetHistos.Fill(sel.GetBJetsForAna().size(), maxSelStep, iChannel, iDataset, weight, weight);
        lljjrecoHistos.Fill(event, candMuon, candElec, sel.GetJetsForAna(), maxSelStep, iChannel, iDataset, weight);
	nofEvtHistos.Fill(maxSelStep, iChannel, iDataset, weight);
}

void SSDiLepAnaHistoManager::Fill(const SSDiLeptonSelection& sel, const int& maxSelStep, const int& iChannel, const int& iDataset, const float& weight){
        elecHistos.Fill(sel.GetSelectedElectrons(), maxSelStep, iChannel, iDataset, weight);
        muHistos.Fill(sel.GetSelectedMuons(), maxSelStep, iChannel, iDataset, weight);
        jetHistos.Fill(sel.GetJetsForAna(), maxSelStep, iChannel, iDataset, weight);
        metHistos.Fill(sel.GetMET(), maxSelStep, iChannel, iDataset, weight);
        bjetHistos.Fill(sel.GetBJetsForAna().size(), maxSelStep, iChannel, iDataset, weight, weight);
	nofEvtHistos.Fill(maxSelStep, iChannel, iDataset, weight);
}



void SSDiLepAnaHistoManager::Compute(){
	elecHistos.MergeChannels();
	muHistos.MergeChannels();
	jetHistos.MergeChannels();
	metHistos.MergeChannels();
	recoHistos.MergeChannels();
	mcHistos.MergeChannels();
	bjetHistos.MergeChannels();
	lljjrecoHistos.MergeChannels();
	nofEvtHistos.MergeChannels();

	//MCStack
	elecHistos.DoMCStack();
	muHistos.DoMCStack();
	jetHistos.DoMCStack();
	metHistos.DoMCStack();
	recoHistos.DoMCStack();
	mcHistos.DoMCStack();
	bjetHistos.DoMCStack();
	lljjrecoHistos.DoMCStack();
	nofEvtHistos.DoMCStack();

	//MCDatasets
	elecHistos.MergeMCDatasets();
	muHistos.MergeMCDatasets();
	jetHistos.MergeMCDatasets();
	metHistos.MergeMCDatasets();
	recoHistos.MergeMCDatasets();
	mcHistos.MergeMCDatasets();
	bjetHistos.MergeMCDatasets();
	lljjrecoHistos.MergeMCDatasets();
	nofEvtHistos.MergeMCDatasets();

	//Cut by cut
	elecHistos.PlotsCutByCut();
	muHistos.PlotsCutByCut();
	jetHistos.PlotsCutByCut();
	metHistos.PlotsCutByCut();
	recoHistos.PlotsCutByCut();
	mcHistos.PlotsCutByCut();
	bjetHistos.PlotsCutByCut();
	lljjrecoHistos.PlotsCutByCut();
	nofEvtHistos.PlotsCutByCut();

}


void SSDiLepAnaHistoManager::ComputeForProof(){
        // Only 1 dataset filled with Proof
        elecHistos.MergeChannels();
        muHistos.MergeChannels();
        jetHistos.MergeChannels();
        metHistos.MergeChannels();
        recoHistos.MergeChannels();
        mcHistos.MergeChannels();
        bjetHistos.MergeChannels();
	lljjrecoHistos.MergeChannels();
	nofEvtHistos.MergeChannels();

        //Cut by cut not done, because Proof can not merge Canvas 

}


void SSDiLepAnaHistoManager::Write(TFile* file){
        TDirectory *fdir = file->GetDirectory("SSDiLepAnaPlots");
        if(!fdir) fdir = file->mkdir("SSDiLepAnaPlots");
	TDirectory * dir = 0;
        dir = fdir->GetDirectory("NofEvents");
	if(!dir) dir = fdir->mkdir("NofEvents");
	nofEvtHistos.Write(dir);
        dir = fdir->GetDirectory("Electrons");
	if(!dir) dir = fdir->mkdir("Electrons");
	elecHistos.Write(dir);
        dir = fdir->GetDirectory("Muons");
        if(!dir) dir = fdir->mkdir("Muons");
        muHistos.Write(dir);
        dir = fdir->GetDirectory("Jets");
        if(!dir) dir = fdir->mkdir("Jets");
        jetHistos.Write(dir);
        dir = fdir->GetDirectory("MET");
        if(!dir) dir = fdir->mkdir("MET");
        metHistos.Write(dir);
        dir = fdir->GetDirectory("RECO");
        if(!dir) dir = fdir->mkdir("RECO");
        recoHistos.Write(dir);
        lljjrecoHistos.Write(dir);
        dir = fdir->GetDirectory("MC");
        if(!dir) dir = fdir->mkdir("MC");
        mcHistos.Write(dir);
        dir = fdir->GetDirectory("BJets");
        if(!dir) dir = fdir->mkdir("BJets");
	bjetHistos.Write(dir);
	dir = 0;
	delete dir;
}

