#include "Plots/interface/Kin4TTbarMetHistoManager.h"
#include "Selection/interface/TTbarMetSelection.h"



Kin4TTbarMetHistoManager::Kin4TTbarMetHistoManager(){
}

Kin4TTbarMetHistoManager::~Kin4TTbarMetHistoManager(){
}


void Kin4TTbarMetHistoManager::CreateHistos(){
        AddHisto(string("HT"),           string("HT"),           string("H_{T} [GeV]"),              100,0.0,1000.0);
        AddHisto(string("M3"),           string("M3"),           string("M^{top had} [GeV]"),        100,0.0,1000.0);
        AddHisto(string("Mtl"),          string("Mtl"),          string("M^{top lep} [GeV]"),        100,0.0,1000.0);
        AddHisto(string("MTw"),          string("MTw"),          string("M_{T}^{W} [GeV]"),          100,0.0,1000.0);
        AddHisto(string("PT3"),          string("PT3"),          string("P_{T}^{top had} [GeV]"),    100,0.,1000.);
        AddHisto(string("PTtl"),         string("PTtl"),         string("P_{T}^{top lep} [GeV]"),    100,0.,1000.);
        AddHisto(string("PTw"),          string("PTw"),          string("P_{T}^{W} [GeV]"),          100,0.0,1000.0);
        AddHisto(string("Dphi_l_met"),   string("Dphi_l_met"),   string("#Delta #Phi (l,MET)"),       30,-1.*TMath::Pi(),TMath::Pi());
        AddHisto(string("Dphi_l_4thjet"),string("Dphi_l_4thjet"),string("#Delta #Phi (l,4^{th}jet)"), 30,-1.*TMath::Pi(),TMath::Pi());
        AddHisto(string("Dphi_2tops"),   string("Dphi_2tops"),   string("#Delta #Phi (tops)"),        30,-1.*TMath::Pi(),TMath::Pi());
        AddHisto(string("Deta_l_thad"),  string("Deta_l_thad"),  string("#Delta #eta (l,top had)"),   50,0.,5.);
        AddHisto(string("Deta_l_4thjet"),string("Deta_l_4thjet"),string("#Delta #eta (l,4^{th}jet)"), 50,0.,5.);
        AddHisto(string("MET"),          string("MET"),          string("MET [GeV]"),                 20,0.0,500.0);
}

void Kin4TTbarMetHistoManager::Fill(const TTbarMetSelection& sel, const int& maxSelStep, const int& iChannel, const int& iDataset, const float& weight){
	if(!Check(iChannel, iDataset)) return;
        // DO NOT FILL BEFORE Njet SELECTION
	for(unsigned int i=4;i<SelectionSteps.size();i++){
		if(maxSelStep>=(int)i) FillSelStep(sel, i , iChannel, iDataset, weight);	
	}	
}

void Kin4TTbarMetHistoManager::FillSelStep(const TTbarMetSelection& sel, const int& iSelStep, const int& iChannel, const int& iDataset, const float& weight){
        if(!Check(iChannel, iSelStep, iDataset, 0) ) return;

/*
        cout << " ht " << sel.HT() 
             << " mtw " << sel.MT_wleptonic() 
             << " ptw " << sel.PT_wleptonic() 
             << " pttl " << sel.PT_topleptonic() << endl;
*/

	Histos[0][iChannel][iSelStep][iDataset].Fill(sel.HT()            ,weight);
	Histos[1][iChannel][iSelStep][iDataset].Fill(sel.M3()            ,weight);
	Histos[2][iChannel][iSelStep][iDataset].Fill(sel.M_topleptonic() ,weight);
	Histos[3][iChannel][iSelStep][iDataset].Fill(sel.MT_wleptonic()  ,weight);
	Histos[4][iChannel][iSelStep][iDataset].Fill(sel.PT_tophad()     ,weight);
	Histos[5][iChannel][iSelStep][iDataset].Fill(sel.PT_topleptonic(),weight);
	Histos[6][iChannel][iSelStep][iDataset].Fill(sel.PT_wleptonic()  ,weight);
	Histos[7][iChannel][iSelStep][iDataset].Fill(sel.Dphi_lmet()  ,weight);
	Histos[8][iChannel][iSelStep][iDataset].Fill(sel.Dphi_ljet4()  ,weight);
	Histos[9][iChannel][iSelStep][iDataset].Fill(sel.Dphi_tops()  ,weight);
	Histos[10][iChannel][iSelStep][iDataset].Fill(sel.Deta_lth()  ,weight);
	Histos[11][iChannel][iSelStep][iDataset].Fill(sel.Deta_ljet4()  ,weight);
	Histos[12][iChannel][iSelStep][iDataset].Fill(sel.Met() ,weight);
}


