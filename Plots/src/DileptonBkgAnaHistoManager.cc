#include "Plots/interface/DileptonBkgAnaHistoManager.h"




DileptonBkgAnaHistoManager::DileptonBkgAnaHistoManager(){
}

DileptonBkgAnaHistoManager::~DileptonBkgAnaHistoManager(){
}


void DileptonBkgAnaHistoManager::AddHisto_(string name, string title, string xaxis, const int& nbins, const float& min, const float& max)
{
	histoNames.push_back(name);	
	AddHisto(name,title,xaxis,nbins,min,max);
}

void DileptonBkgAnaHistoManager::AddHisto2D_(string name, string title, string xaxis, const int& nxbins, const float& xmin, const float& xmax, string yaxis, const int& nybins, const float& ymin, const float& ymax) 
{
	histo2DNames.push_back(name);
	AddHisto2D(name,title,xaxis,nxbins,xmin,xmax,yaxis,nybins,ymin,ymax);
}

void DileptonBkgAnaHistoManager::CreateHistos(){
  AddHisto_(string("DecayCase"),        string("Decay case  (0=loste/mu, 1=hadronictau, 2=leptonictau, 3=2taus)"),	string("Case"),4,-0.5,3.5);
 
  //todo
  AddHisto_(string("FoundMCLepton_DeltaRSel"),string("Delta R between NOTLOST MC lepton and SELECTED lepton"),string("DeltaR"),50,0,2);
  
  AddHisto_(string("LostMCLepton_WhenIsItLost"),     string("Lost lepton was lost before selection (0) or after (1)"),			string("Case"),2,-0.5,1.5);
  AddHisto_(string("LostMCLepton_Eta"),    string("Lost MC lepton : Eta"),	                                 string("Eta"),   50,-5,5);
  AddHisto_(string("LostMCLepton_EtaInt"), string("Eta of lost MC lepton ; <1.4 ; <1.6 ; <2.5 ; <inf"),        string("AbsEta"),4,-0.5,3.5);
  AddHisto_(string("LostMCLepton_Pt"),     string("Lost MC lepton : Pt"),	                                     string("Pt"),    50,0,200);
  AddHisto_(string("LostMCLepton_Type"),   string("Lost MC lepton : type (1=el, 2=mu, negative if from tau)"), string("Type"),  6,-3.5,2.5);
  
  AddHisto_(string("PFMC_MatchSeemsOK"),string("Match MC<->PFCand : doesMatchSeemsOK (0 : no ; 1 : yes)"),	string("Answer"),2,-0.5,1.5);
  AddHisto_(string("PFMC_DeltaPtRel"),  string("Match MC<->PFCand : DeltaPt/Pt"),		string("DeltaPt/Pt"),50,-1.2,1.2);
  AddHisto_(string("PFMC_DeltaEta"),    string("Match MC<->PFCand : DeltaEta"),			string("DeltaEta"),50,-3,3);
  AddHisto_(string("PFMC_DeltaR"),      string("Match MC<->PFCand : DeltaR"),			string("DeltaR"),50,0,0.5);
  AddHisto_(string("PFMC_Pt"),          string("Match MC<->PFCand : Pt(pfcand)"),		string("Pt"),50,0,150);
  AddHisto_(string("PFMC_Eta"),         string("Eta of matched pf candidate wrt lost MC lepton"),string("Eta"),50,0,200);
  AddHisto_(string("PFMC_EtaInt"),      string("EtaInt of matched pf candidate wrt lost MC lepton ; <1.4 ; <1.6 ; <2.5 ; <inf"),string("EtaInt"),4,-0.5,3.5);
  AddHisto_(string("PFMC_dz"),          string("Match MC<->PFCand : fabs(dz(pfcand))"),	string("dz"),50,0,0.1);
  AddHisto_(string("PFMC_iso"),         string("Match MC<->PFCand : reliso(pfcand)"),	string("reliso"),50,0,0.2);
  
  AddHisto_(string("PATMC_MatchSeemsOK"),string("Match MC<->PAT : doesMatchSeemsOK (0 : no ; 1 : yes)"),	string("Answer"),2,-0.5,1.5);
  AddHisto_(string("PATMC_DeltaPtRel"),	string("Match MC<->PAT : DeltaPt/Pt"),			string("DeltaPt/Pt"),50,-1.2,1.2);
  AddHisto_(string("PATMC_DeltaEta"),	string("Match MC<->PAT : DeltaEta"),			string("DeltaEta"),50,-3,3);
  AddHisto_(string("PATMC_DeltaR"),		string("Match MC<->PAT : DeltaR"),				string("DeltaR"),50,0,0.5);
  AddHisto_(string("PATMC_Pt"),			string("Match MC<->PAT : Pt"),					string("Pt"),50,0,150);
  AddHisto_(string("PATMC_FailingStep"),	string("Match MC<->PAT : Failing step during selection (negative if muon)"),string("Pt"),30,-14.5,15.5);
 
  AddHisto_(string("HadrTauMC_Pt"),string("Pt of MC hadronic tau"),string("Pt"),50,0,150);
  AddHisto_(string("HadrTauMC_Eta"),string("Eta of MC hadronic tau"),string("Eta"),50,-5,5);
  AddHisto_(string("HadrTauMC_EtaInt"),string("Eta of MC hadronic tau ; <1.4 ; <1.6 ; <2.5 ; <inf"),string("AbsEta"),4,-0.5,3.5);
  
  AddHisto_(string("HadrTauMC_decayModeCompact"),string("Hadronic tau - Compact decay mode (10*charged+neutral)"),string("mode"),56,-0.5,55.5);
  AddHisto_(string("HadrTauMC_NneutralDecay"),string("Hadronic tau - Number of neutral decay (ex. photons)"),string("n"),6,-0.5,5.5);
  AddHisto_(string("HadrTauMC_NchargedDecay"),string("Hadronic tau - Number of charged decay"),string("n"),6,-0.5,5.5);
  AddHisto_(string("HadrTauMC_Nphotons"),string("Hadronic tau - Number of photons"),string("n"),6,-0.5,5.5);
  AddHisto_(string("HadrTauMC_leptonFlavor"),string("Hadronic tau - Lepton flavor (expect 0)"),string("n"),3,-0.5,2.5);

  AddHisto2D_(string("LostMCLepton_PtvsEta"),string("Pt of lost MC lepton (vs Eta)"),string("Pt"),50,0,200,string("Eta"),50,-5,5);
  AddHisto2D_(string("LostMCLepton_PtvsEtaInt"),string("Pt of lost MC lepton (vs Eta Int)"),string("Pt"),50,0,200,string("Eta"),4,-0.5,3.5);

  AddHisto2D_(string("HadrTauMC_PtvsEta"),string("Pt of MC hadronic tau (vs Eta)"),string("Pt"),50,0,200,string("Eta"),50,-5,5);
  AddHisto2D_(string("HadrTauMC_PtvsEtaInt"),string("Pt of MC hadronic tau (vs Eta Int)"),string("Pt"),50,0,200,string("Eta"),4,-0.5,3.5);

  AddHisto_(string("HadrTauPFMatch_DeltaR"),string("dR between TauReco and matched pfcandidate"),string("DeltaR"),50,0,0.5);
  AddHisto2D_(string("HadrTauPFMatch_ChargedDecayvsPt"),string("Number of charged decay vs Pt"),string("ChargedDecay"),6,-0.5,5.5,string("Pt"),50,0,200);
  AddHisto2D_(string("HadrTauPFMatch_ChargedDecayvsdz"),string("Number of charged decay vs dz"),string("ChargedDecay"),6,-0.5,5.5,string("dz"),50,0,0.1);
  
  AddHisto2D_(string("HadrTauPFMatch_ChargedDecayvsIso03"),string("Number of charged decay vs Iso03"),string("ChargedDecay"),6,-0.5,5.5,string("Iso"),200,0,2);
  AddHisto2D_(string("HadrTauPFMatch_ChargedDecayvsIso05"),string("Number of charged decay vs Iso05"),string("ChargedDecay"),6,-0.5,5.5,string("Iso"),200,0,2);
  AddHisto2D_(string("HadrTauPFMatch_ChargedDecayvsIso05no015"),string("Number of charged decay vs Iso05no015"),string("ChargedDecay"),6,-0.5,5.5,string("Iso"),200,0,2);
  AddHisto2D_(string("HadrTauPFMatch_ChargedDecayvsNTrack015"),string("Number of charged decay vs NTrack015"),string("ChargedDecay"),6,-0.5,5.5,string("NTracks"),20,-0.5,19.5);
  AddHisto2D_(string("HadrTauPFMatch_ChargedDecayvsNTrack03"),string("Number of charged decay vs NTrack03"),string("ChargedDecay"),6,-0.5,5.5,string("NTracks"),20,-0.5,19.5);
  AddHisto2D_(string("HadrTauPFMatch_ChargedDecayvsNTrack05"),string("Number of charged decay vs NTrack05"),string("ChargedDecay"),6,-0.5,5.5,string("NTracks"),20,-0.5,19.5);
  AddHisto2D_(string("HadrTauPFMatch_ChargedDecayvsNTrack015Pt1GeV"),string("Number of charged decay vs NTrack015Pt1GeV"),string("ChargedDecay"),6,-0.5,5.5,string("NTracks"),20,-0.5,19.5);
  AddHisto2D_(string("HadrTauPFMatch_ChargedDecayvsNTrack03Pt1GeV"),string("Number of charged decay vs NTrack03Pt1GeV"),string("ChargedDecay"),6,-0.5,5.5,string("NTracks"),20,-0.5,19.5);
  AddHisto2D_(string("HadrTauPFMatch_ChargedDecayvsNTrack05Pt1GeV"),string("Number of charged decay vs NTrack05Pt1GeV"),string("ChargedDecay"),6,-0.5,5.5,string("NTracks"),20,-0.5,19.5);

  AddHisto_(string("HadrTauPFMath_1ChargedDecay_WhyNotVetoed"),string("(1 charged) Reason why not vetoed : 0 = pT, 1 = dz, 2 = iso"),string("Reason"),3,-0.5,2.5);
  AddHisto_(string("HadrTauPFMath_3OrMoreChargedDecay_WhyNotVetoed"),string("(3+ charged decay) Reason why not vetoed, 0 = pT, 1 = dz, 2 = iso"),string("Reason"),3,-0.5,2.5);


  AddHisto2D_(string("HadrTauPFMatch_1ChargedIso03vsIso05no015"),string("(1 charged) Iso03 vs Iso05no015"),string("Iso03"),200,0,2,string("Iso05No015"),200,0,2);
  AddHisto2D_(string("HadrTauPFMatch_3OrMoreChargedIso03vsIso05no015"),string("(3 charged) Iso03 vs Iso05no015"),string("Iso03"),200,0,2,string("Iso05No015"),200,0,2);


}

void DileptonBkgAnaHistoManager::Fill(std::vector<std::string> dataName, std::vector<float> dataValue, const int& maxSelStep, const int& iChannel, const int& iDataset, const float& weight)
{
	if(!Check(iChannel, iDataset)) return;
	for(unsigned int i=0;i<SelectionSteps.size();i++)
		if(maxSelStep>=(int)i) FillSelStep(dataName,dataValue, i , iChannel, iDataset, weight);	
}


void DileptonBkgAnaHistoManager::Fill2D( std::vector<std::string> dataName, std::vector<float> dataValueX, std::vector<float> dataValueY, const int& maxSelStep, const int& iChannel, const int& iDataset, const float& weight)
{
	if(!Check(iChannel, iDataset)) return;
	for(unsigned int i=0;i<SelectionSteps.size();i++)
	{
		if(maxSelStep>=(int)i) Fill2DSelStep(dataName,dataValueX,dataValueY, i , iChannel, iDataset, weight);
	}
}

void DileptonBkgAnaHistoManager::FillSelStep(std::vector<std::string> dataName, std::vector<float> dataValue, const int& iSelStep, const int& iChannel, const int& iDataset, const float& weight){

		if(!Check(iChannel, iSelStep, iDataset, 0) ) return;

		for (unsigned int i = 0 ; i < dataName.size(); i++)
		{
			int histoIndex = -1;
			for (unsigned int j = 0 ; j < histoNames.size() ; j++)
			{
				if (histoNames[j] == dataName[i].c_str()) histoIndex = j;
			}

			if (histoIndex == -1) cout << "Warning : histo " << dataName[i] << " not found." << endl;
			else
			Histos[histoIndex][iChannel][iSelStep][iDataset].Fill(dataValue[i],weight);
			
		}
}

void DileptonBkgAnaHistoManager::Fill2DSelStep(std::vector<std::string> dataName,std::vector<float> dataValueX, std::vector<float> dataValueY, const int& iSelStep, const int& iChannel, const int& iDataset, const float& weight){

		if(!Check2D(iChannel, iSelStep, iDataset, 0) ) return;

		for (unsigned int i = 0 ; i < dataName.size() ; i++)
		{
			int histoIndex = -1;
			for (unsigned int j = 0 ; j < histo2DNames.size() ; j++)
			{
				if (histo2DNames[j] == dataName[i].c_str()) histoIndex = j;
			}

			if (histoIndex == -1) cout << "Warning : histo " << dataName[i] << " not found." << endl;
			else Histos2D[histoIndex][iChannel][iSelStep][iDataset].Fill(dataValueX[i],dataValueY[i],weight);
		}
}
