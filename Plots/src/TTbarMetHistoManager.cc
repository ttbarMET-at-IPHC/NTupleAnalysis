#include "Plots/interface/TTbarMetHistoManager.h"



TTbarMetHistoManager::TTbarMetHistoManager(){
}

TTbarMetHistoManager::~TTbarMetHistoManager(){
}

void TTbarMetHistoManager::LoadDatasets(vector<Dataset> datasets){
	elecHistos.LoadDatasets(datasets);
	muHistos.LoadDatasets(datasets);
	jetHistos.LoadDatasets(datasets);
	metHistos.LoadDatasets(datasets);
	bjetHistos.LoadDatasets(datasets);
	nofEvtHistos.LoadDatasets(datasets);
	kinHistos.LoadDatasets(datasets);
}    
 	

void TTbarMetHistoManager::LoadSelectionSteps(vector<string> selectionSteps){
	elecHistos.LoadSelectionSteps(selectionSteps);
	muHistos.LoadSelectionSteps(selectionSteps);
	jetHistos.LoadSelectionSteps(selectionSteps);
	metHistos.LoadSelectionSteps(selectionSteps);
	bjetHistos.LoadSelectionSteps(selectionSteps);
	nofEvtHistos.LoadSelectionSteps(selectionSteps);
	kinHistos.LoadSelectionSteps(selectionSteps);
}
	
 
void TTbarMetHistoManager::LoadChannels(vector<string> channels){
	elecHistos.LoadChannels(channels);
	muHistos.LoadChannels(channels);
	jetHistos.LoadChannels(channels);
	metHistos.LoadChannels(channels);
	bjetHistos.LoadChannels(channels);
	nofEvtHistos.LoadChannels(channels);
	kinHistos.LoadChannels(channels);
}

void TTbarMetHistoManager::Clear(){
	elecHistos.Clear();
	muHistos.Clear();
	jetHistos.Clear();
	metHistos.Clear();
	bjetHistos.Clear();
	nofEvtHistos.Clear();
	kinHistos.Clear();
}


void TTbarMetHistoManager::CreateHistos(){
	elecHistos.CreateHistos();
	muHistos.CreateHistos();
	jetHistos.CreateHistos();
	metHistos.CreateHistos();
	bjetHistos.CreateHistos();
	nofEvtHistos.CreateHistos();
	kinHistos.CreateHistos();
}

void TTbarMetHistoManager::Fill(const TTbarMetSelection& sel, NTEvent* event, const vector<NTMuon>& candMuon, const vector<NTElectron>& candElec, const int& maxSelStep, const int& iChannel, const int& iDataset, const float& weight){
        elecHistos.Fill(candElec, maxSelStep, iChannel, iDataset, weight);
        muHistos.Fill(candMuon, maxSelStep, iChannel, iDataset, weight);
        jetHistos.Fill(sel.GetJetsForAna(), maxSelStep, iChannel, iDataset, weight);
        metHistos.Fill(sel.GetMET(), maxSelStep, iChannel, iDataset, weight);
        bjetHistos.Fill(sel.GetBJetsForAna().size(), maxSelStep, iChannel, iDataset, weight, weight);
	nofEvtHistos.Fill(maxSelStep, iChannel, iDataset, weight);
	kinHistos.Fill(sel,maxSelStep, iChannel, iDataset, weight);
}

void TTbarMetHistoManager::Fill(const TTbarMetSelection& sel, const int& maxSelStep, const int& iChannel, const int& iDataset, const float& weight){
        elecHistos.Fill(sel.GetElectronsForAna(), maxSelStep, iChannel, iDataset, weight);
        muHistos.Fill(sel.GetMuonsForAna(), maxSelStep, iChannel, iDataset, weight);
        jetHistos.Fill(sel.GetJetsForAna(), maxSelStep, iChannel, iDataset, weight);
        metHistos.Fill(sel.GetMET(), maxSelStep, iChannel, iDataset, weight);
        bjetHistos.Fill(sel.GetBJetsForAna().size(), maxSelStep, iChannel, iDataset, weight, weight);
	nofEvtHistos.Fill(maxSelStep, iChannel, iDataset, weight);
	kinHistos.Fill(sel,maxSelStep, iChannel, iDataset, weight);
}



void TTbarMetHistoManager::Compute(){
	elecHistos.MergeChannels();
	muHistos.MergeChannels();
	jetHistos.MergeChannels();
	metHistos.MergeChannels();
	bjetHistos.MergeChannels();
	nofEvtHistos.MergeChannels();
	kinHistos.MergeChannels();

	//MCStack
	elecHistos.DoMCStack();
	muHistos.DoMCStack();
	jetHistos.DoMCStack();
	metHistos.DoMCStack();
	bjetHistos.DoMCStack();
	nofEvtHistos.DoMCStack();
	kinHistos.DoMCStack();

	//MCDatasets
	elecHistos.MergeMCDatasets();
	muHistos.MergeMCDatasets();
	jetHistos.MergeMCDatasets();
	metHistos.MergeMCDatasets();
	bjetHistos.MergeMCDatasets();
	nofEvtHistos.MergeMCDatasets();
	kinHistos.MergeMCDatasets();

	//Cut by cut
	elecHistos.PlotsCutByCut();
	muHistos.PlotsCutByCut();
	jetHistos.PlotsCutByCut();
	metHistos.PlotsCutByCut();
	bjetHistos.PlotsCutByCut();
	nofEvtHistos.PlotsCutByCut();
	kinHistos.PlotsCutByCut();

}


void TTbarMetHistoManager::ComputeForProof(){
        // Only 1 dataset filled with Proof
        elecHistos.MergeChannels();
        muHistos.MergeChannels();
        jetHistos.MergeChannels();
        metHistos.MergeChannels();
        bjetHistos.MergeChannels();
	nofEvtHistos.MergeChannels();
	kinHistos.MergeChannels();

        //Cut by cut not done, because Proof can not merge Canvas 

}


void TTbarMetHistoManager::Write(TFile* file){
        TDirectory *fdir = file->GetDirectory("TTbarMetAnaPlots");
        if(!fdir) fdir = file->mkdir("TTbarMetAnaPlots");
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
        dir = fdir->GetDirectory("BJets");
        if(!dir) dir = fdir->mkdir("BJets");
	bjetHistos.Write(dir);
        dir = fdir->GetDirectory("KinPlots");
        if(!dir) dir = fdir->mkdir("KinPlots");
	kinHistos.Write(dir);
	dir = 0;
	delete dir;
}

