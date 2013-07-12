#include <iomanip>
#include <iostream>
#include <time.h>
#include "NTFormat/interface/NTEvent.h"

#include "Selection/interface/TTbarMetSelection.h"
#include "Selection/interface/SelectionTable.h"
#include "Tools/interface/Dataset.h"
#include "Tools/interface/AnalysisEnvironmentLoader.h"
#include "Plots/interface/TTbarMetHistoManager.h"
#include "Tools/interface/PUWeighting.h"
#include "Tools/interface/LumiReweightingStandAlone.h"
#include "Tools/interface/JetCorrector.h"

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>

using namespace IPHCTree;
using namespace std;

#include "Tools/interface/TableScrewdriver.h"
#include "eventListSynchro.h"
#include "miscTools.h"
#include "DileptonBackgroundStudy.h"


// #########################################################################""
// #########################################################################""
// #########################################################################""

typedef struct 
{
	Float_t iso_02;
	Float_t iso_03;
	Float_t iso_03_pT1GeV;
	Float_t iso_04;
	Float_t iso_05;
	Float_t iso_05_no01;
	Float_t iso_05_no015;
	Float_t iso_05_no015_pT1GeV;
	Float_t iso_05_no02;
	Float_t iso_045_no015;
	Float_t iso_055_no015;
	Float_t iso_05_no015_pTEqualIso02;
	
	Float_t hadronicMC_dRmeanDecay;
	Float_t hadronicMC_dRmaxDecay;
	Float_t hadronicMC_dRmeanChargedDecay;
	Float_t hadronicMC_dRmaxChargedDecay;
	Float_t hadronicMC_pTminChargedDecay;
	Float_t hadronicMC_pTmeanChargedDecay;
	Float_t hadronicMC_pTmaxChargedDecay;

	Float_t bestIso03_iso_03;
	Float_t bestIso03_iso_05_no015;
	Float_t bestIso03_iso_05_no01;
	Float_t bestIso03_iso_05_ntrk_01;
	Float_t bestIso03_iso_05_ntrk_015;
	Float_t bestIso03_iso_05_ntrk_015_pT1GeV;

	Float_t bestIso05no015_iso_03;
	Float_t bestIso05no015_iso_05_no015;
	Float_t bestIso05no015_iso_05_no01;
	Float_t bestIso05no015_iso_05_ntrk_01;
	Float_t bestIso05no015_iso_05_ntrk_015;
	Float_t bestIso05no015_iso_05_ntrk_015_pT1GeV;

	Float_t bestIso05no01_iso_03;
	Float_t bestIso05no01_iso_05_no015;
	Float_t bestIso05no01_iso_05_no01;
	Float_t bestIso05no01_iso_05_ntrk_01;
	Float_t bestIso05no01_iso_05_ntrk_015;
	Float_t bestIso05no01_iso_05_ntrk_015_pT1GeV;

	Float_t eventInfo_isSemilep;
	Float_t eventInfo_isDilep;
	Float_t eventInfo_haveTau;
	Float_t eventInfo_haveHadronicTau;
	Float_t eventInfo_haveHadronicTau3;

} 
dileptonData;


void resetDileptonData(dileptonData* dileptonDataForTree_)
{
	(*dileptonDataForTree_).eventInfo_isSemilep = 0.0;
	(*dileptonDataForTree_).eventInfo_isDilep = 0.0;
	(*dileptonDataForTree_).eventInfo_haveTau = 0.0;
	(*dileptonDataForTree_).eventInfo_haveHadronicTau = 0.0;
	(*dileptonDataForTree_).eventInfo_haveHadronicTau3 = 0.0;

	(*dileptonDataForTree_).iso_02 = 10.0;
	(*dileptonDataForTree_).iso_03 = 10.0;
	(*dileptonDataForTree_).iso_03_pT1GeV = 10.0;
	(*dileptonDataForTree_).iso_04 = 10.0;
	(*dileptonDataForTree_).iso_05 = 10.0;
	(*dileptonDataForTree_).iso_05_no01 = 10.0;
	(*dileptonDataForTree_).iso_05_no015 = 10.0;
	(*dileptonDataForTree_).iso_05_no015_pT1GeV = 10.0;
	(*dileptonDataForTree_).iso_05_no02 = 10.0;
	(*dileptonDataForTree_).iso_045_no015 = 10.0;
	(*dileptonDataForTree_).iso_055_no015 = 10.0;

	(*dileptonDataForTree_).bestIso03_iso_03 = 10.0;
	(*dileptonDataForTree_).bestIso03_iso_05_no015 = 10.0;
	(*dileptonDataForTree_).bestIso03_iso_05_no01 = 10.0;
	(*dileptonDataForTree_).bestIso03_iso_05_ntrk_01 = -1.0;
	(*dileptonDataForTree_).bestIso03_iso_05_ntrk_015 = -1.0;
	(*dileptonDataForTree_).bestIso03_iso_05_ntrk_015_pT1GeV = -1.0;

	(*dileptonDataForTree_).bestIso05no015_iso_03 = 10.0;
	(*dileptonDataForTree_).bestIso05no015_iso_05_no015 = 10.0;
	(*dileptonDataForTree_).bestIso05no015_iso_05_no01 = 10.0;
	(*dileptonDataForTree_).bestIso05no015_iso_05_ntrk_01 = -1.0;
	(*dileptonDataForTree_).bestIso05no015_iso_05_ntrk_015 = -1.0;
	(*dileptonDataForTree_).bestIso05no015_iso_05_ntrk_015_pT1GeV = -1.0;

	(*dileptonDataForTree_).bestIso05no01_iso_03 = 10.0;
	(*dileptonDataForTree_).bestIso05no01_iso_05_no015 = 10.0;
	(*dileptonDataForTree_).bestIso05no01_iso_05_no01 = 10.0;
	(*dileptonDataForTree_).bestIso05no01_iso_05_ntrk_01 = -1.0;
	(*dileptonDataForTree_).bestIso05no01_iso_05_ntrk_015 = -1.0;
	(*dileptonDataForTree_).bestIso05no01_iso_05_ntrk_015_pT1GeV = -1.0;

	(*dileptonDataForTree_).iso_05_no015_pTEqualIso02 = 10.0;

	(*dileptonDataForTree_).hadronicMC_dRmeanDecay = -1.0;
	(*dileptonDataForTree_).hadronicMC_dRmaxDecay = -1.0;
	(*dileptonDataForTree_).hadronicMC_dRmeanChargedDecay = -1.0;
	(*dileptonDataForTree_).hadronicMC_dRmaxChargedDecay = -1.0;
	(*dileptonDataForTree_).hadronicMC_pTminChargedDecay = -1.0;
	(*dileptonDataForTree_).hadronicMC_pTmeanChargedDecay = -1.0;
	(*dileptonDataForTree_).hadronicMC_pTmaxChargedDecay = -1.0;

}
// #########################################################################""
// #########################################################################""
// #########################################################################""

//testStruct myStruct;

dileptonData dileptonDataForTree;

int main (int argc, char *argv[])
{
  cout << endl;
  cout << "   ,-----------------------," << endl;
  cout << "   |   Starting analysis   |" << endl;
  cout << "   `-----------------------`" << endl;
  cout << endl;

  // ##################################
  // #	Some parameters for the user  #
  // ##################################
 
  bool doHistoManager = false;		// False will not produce histos
  bool doLatexOutput = true;		// False will not produce latex table
  bool doPileUpReweigthing = false;	// False will not use pileup reweighting

  // (can be overwritten by argv[1])
  string defaultConfiguration("../../config/TTbarMETAnalysis.xml");
  // (can be overwritten by argv[2])
  string defaultLatexTableName("CutflowTable");
  // (can be overwritten by argv[3])
  string defaultRootOutputName("TTbarMETanalysis.root");

  // ############################
  // #	Initializing variables  #
  // ############################
  
  INFO1_MSG << "Initializing variables..." << endl;
 
  vector < Dataset > datasets;
  TTbarMetSelection sel;
  
  float Luminosity = 0;
  float LumiError = 0.;
  int DataType = 0;				// DataType : 0: MC - 1: Data - 2 Data & MC
  int verbosity = -1;


  reweight::LumiReWeighting *LumiWeights;
  IPHCTree::NTEvent * event = 0;


  std::vector<std::string> dileptonCol;
  dileptonCol.push_back("l+jets");
  dileptonCol.push_back("ll");
  dileptonCol.push_back("llHadrTau");
  dileptonCol.push_back("llHadrTau3+");

  std::vector<std::string> dileptonRow;
  dileptonRow.push_back("baseline");
  dileptonRow.push_back("MET");
  dileptonRow.push_back("stdVeto");
  dileptonRow.push_back("stdVetoEff");
  dileptonRow.push_back("vetoTest0");
  dileptonRow.push_back("vetoTest0Eff");
  dileptonRow.push_back("vetoTest1");
  dileptonRow.push_back("vetoTest1Eff");
  dileptonRow.push_back("vetoTest2");
  dileptonRow.push_back("vetoTest2Eff");
  dileptonRow.push_back("vetoTest3");
  dileptonRow.push_back("vetoTest3Eff");

  // Reading parameters from argv
  // 	-> configuration file
  string xmlFileName;
  if (argc > 1) xmlFileName = string(argv[1]);
  else          xmlFileName = defaultConfiguration;
  // 	-> root output name
  string rootOutputName;
  if (argc > 2) rootOutputName = string(argv[2]);
  else 			rootOutputName = defaultRootOutputName;
  // 	-> latex table name
  string latexTableName;
  if (argc > 3) latexTableName = string(argv[3]);
  else latexTableName = defaultLatexTableName;
  
  // #############################
  // # 	 Loading configuration   #
  // #############################

  cout << endl;
  INFO1_MSG << "Loading configuration..." << endl;
  cout << "        (config : " << xmlFileName << ")" << endl;
  
 
  AnalysisEnvironmentLoader anaEL (xmlFileName);
  anaEL.LoadSamples (datasets);	// now the list of datasets written in the xml file is known
  anaEL.LoadSelection (sel);	// now the parameters for the selection are given to the selection // no specific TTbarMET parameters
  anaEL.LoadGeneralInfo (DataType, Luminosity, LumiError, verbosity);

  // #########################################
  // # 	 Creating tables and histomanagers   #
  // #########################################

  cout << endl;
  INFO1_MSG << "Creating tables and histomanagers..." << endl;

  // Stuff
  
  TFile dileptonTreeFile("dileptonTree.root","RECREATE","Tree for Dilepton background study");
  TTree *dileptonTree = new TTree("dileptonTree","Tree for Dilepton background study");
  
  dileptonTree->Branch("dileptonData",&dileptonDataForTree,"iso_02:iso_03:iso_03_pT1GeV:iso_04:iso_05:iso_05_no01:iso_05_no015:iso_05_no015_pT1GeV:iso_05_no02:iso_045_no015:iso_055_no015:iso_05_no015_pTEqualIso02:hadronicMC_dRmeanDecay:hadronicMC_dRmaxDecay:hadronicMC_dRmeanChargedDecay:hadronicMC_dRmaxChargedDecay:hadronicMC_pTminChargedDecay:hadronicMC_pTmeanChargedDecay:hadronicMC_pTmaxChargedDecay:bestIso03_iso_03:bestIso03_iso_05_no015:bestIso03_iso_05_no01:bestIso03_iso_05_ntrk_01:bestIso03_iso_05_ntrk_015:bestIso03_iso_05_ntrk_015_pT1GeV:bestIso05no015_iso_03:bestIso05no015_iso_05_no015:bestIso05no015_iso_05_no01:bestIso05no015_iso_05_ntrk_01:bestIso05no015_iso_05_ntrk_015:bestIso05no015_iso_05_ntrk_015_pT1GeV:bestIso05no01_iso_03:bestIso05no01_iso_05_no015:bestIso05no01_iso_05_no01:bestIso05no01_iso_05_ntrk_01:bestIso05no01_iso_05_ntrk_015:bestIso05no01_iso_05_ntrk_015_pT1GeV:eventInfo_isSemilep:eventInfo_isDilep:eventInfo_haveTau:eventInfo_haveHadronicTau:eventInfo_haveHadronicTau3");




  TableScrewdriver dileptonTable(dileptonCol,dileptonRow);
  dileptonTable.PrintTable();
  
  
  // Tables
  
  SelectionTable selTable_e (sel.GetCutList (), datasets, string ("e"));
  SelectionTable selTable_mu (sel.GetCutList (), datasets, string ("µ"));
  SelectionTable selTable_l (sel.GetCutList (), datasets, string ("lepton"));
  
  // Histomanagers
  
  TTbarMetHistoManager histoManager;
  if(doHistoManager)
  {
  	histoManager.LoadDatasets(datasets);
  	histoManager.LoadSelectionSteps(sel.GetCutList ());
  	histoManager.LoadChannels(sel.GetChannelList ());
  	histoManager.CreateHistos();
  }

  // HistoManager specific to Dilepton background analysis
  
  DileptonBkgAnaHistoManager histoManagerDileptonBkg;
  if(doHistoManager)
  {
  	histoManagerDileptonBkg.LoadDatasets (datasets);
  	histoManagerDileptonBkg.LoadSelectionSteps (sel.GetCutList ());
  	histoManagerDileptonBkg.LoadChannels (sel.GetChannelList ());
  	histoManagerDileptonBkg.CreateHistos ();
  }


  // ##############################
  // # 	 Printing general infos   #
  // ##############################
  
  cout << endl;
  cout << "   -----------------------------" << endl;
  cout << "     Verbosity mode : " << verbosity << endl;
  cout << "     Luminosity : " << Luminosity << endl;
  cout << "     DataType (whole config) : ";
  switch (DataType) 
  {
  	case 0: { cout << "MC" << endl;   break; }
  	case 1: { cout << "Data" << endl; break; }
  	case 2: { cout << "Data & MC" << endl; break; }
  	default: { cout << "Unknwon" << endl; break; }
  }
  cout << "   -----------------------------" << endl;
  
  // #########################################
  // # 	 Start loop over the dataset infos   #
  // #########################################
  if (verbosity > 0)
  {
  	cout << endl;
  	cout << "   ,---------------------------------," << endl;
  	cout << "   |   Starting loop over datasets   |" << endl;
  	cout << "   `---------------------------------`" << endl;
  	cout << endl;
  }

  int nEvents_tot;
  for (unsigned int datasetId = 0; datasetId < datasets.size (); datasetId++) 
  {

	  // ########################
	  // # 	 Load the dataset   #
	  // ########################
   
	INFO1_MSG << "Loading next dataset..." << endl;

    datasets[datasetId].eventTree()->SetBranchAddress ("NTEvent", &event);
    unsigned int nEvents = static_cast<unsigned int>(datasets[datasetId].eventTree()->GetEntries ());
    nEvents_tot = nEvents;

    if (verbosity > 2)
	{
		cout << endl;
  		cout << "         [ Dataset n°" << datasetId+1 << " ]" << endl;
		cout << "         " << datasets[datasetId].Name() << endl;
		cout << endl;
  		INFO1_MSG << "   NEvents total    : " << nEvents << endl;
  		INFO1_MSG << "NEvents to run over : " << datasets[datasetId].NofEvtsToRunOver() << endl;
		cout << endl;
	}

	// ########################
	// #   Load the pile-up   #
	// ########################
	
	INFO1_MSG << "Loading pile-up informations..." << endl;
   
    // PU from JL's code
	if (datasets[datasetId].isData() == false) {
//   if(datasets[datasetId].Name() == "DY1" || datasets[datasetId].Name() == "signal" ){
//      string datafile = "/opt/sbg/data/data1/cms/ccollard/CMSSW/CMSSW_4_2_8_patch7/src/NTuple/NTupleAnalysis/macros/data/PUdata.root";
//      string mcfile   = "/opt/sbg/data/data1/cms/ccollard/CMSSW/CMSSW_4_2_8_patch7/src/NTuple/NTupleAnalysis/macros/data/PU3DMC_Fall11.root";
//
//      LumiWeights    = new reweight::LumiReWeighting(mcfile, datafile, "histoMCPU", "pileup" );
//    }
//    else{
      string datafile = "/opt/sbg/data/data1/cms/ccollard/CMSSW/CMSSW_4_2_8_patch7/src/NTuple/NTupleAnalysis/macros/data/default73.5mb.root";
      string mcfile   = "/opt/sbg/data/data1/cms/ccollard/CMSSW/CMSSW_4_2_8_patch7/src/NTuple/NTupleAnalysis/macros/data/PU3DMC.root";

      LumiWeights    = new reweight::LumiReWeighting(mcfile, datafile, "histoMCPU", "pileup" );
      LumiWeights->weight3D_init( 1. );

//    }
    }
 
	// ############################
	// #   Loop over the events   #
	// ############################
			  
	if (verbosity > 0)
	{
 		cout << endl;
 	 	cout << "   ,-------------------------------," << endl;
 	 	cout << "   |   Starting loop over events   |" << endl;
  		cout << "   `-------------------------------`" << endl;
  		cout << endl;
	}

    for (unsigned int ievt = 0; ievt < datasets[datasetId].NofEvtsToRunOver(); ievt++)
    {

	  // Display some information

      if (verbosity > 3)
      {
		        INFO1_MSG <<  "run: " << event->general.runNb;
		        std::cout << " lumi: " << event->general.lumiblock;
      		    std::cout << " event: " << event->general.eventNb;
		  		std::cout << std::endl;
      }

      if (ievt % 100000 == 0) printProgressBar(ievt,datasets[datasetId].NofEvtsToRunOver());

	  // Load the event
	   
	  datasets[datasetId].eventTree()->GetEntry(ievt);
	  IPHCTree::NTTransient::InitializeAfterReading(event);
	  int eventId = event->general.eventNb;
	  sel.LoadEvent(event);

	  // Get the weight for pile-up reweighting

	  float weight = 1.;
	  if (doPileUpReweigthing)
	  {
	  		if(datasets[datasetId].isData() == false) weight = datasets[datasetId].NormFactor()*Luminosity; //if Data , weight = 1

      		float weightpu=1.;
      		if(datasets[datasetId].isData() == false) 
			{ 
				// MC
  	//      if( datasets[datasetId].Name() == "DY1" || datasets[datasetId].Name() == "signal")  
	//      {
	//          double ave_npu = (event->pileup.before_npu+event->pileup.intime_npu+event->pileup.after_npu)/3.;
	//          weightpu = LumiWeights->ITweight3BX(ave_npu);
	//      }
	//      else 
	//      {
         	 	weightpu = LumiWeights->weight3D(event->pileup.before_npu, event->pileup.intime_npu, event->pileup.after_npu);
	//      }
	//      weight *= weightpu; //if Data , weight = 1
      		}
	//      else 
	//      { 
	//      	// DATA
	//         JEC_L2L3Residuals.ApplyCorrections(event); // n'appliquer la correction que pour les donnees
	//      }
	  }


	  // Apply full selection
	  int triggerME = 0;
      int selLastStep = sel.doFullSelection (&(datasets[datasetId]), string("all"),&triggerME);


	  // Get trigger info from selection
      bool trigger_e = false;
      bool trigger_mu = false;
	  if (triggerME % 1 == 1) trigger_e  = true;
	  if (triggerME    >= 10) trigger_mu = true; 

	  // Also get the lepton type
      int lep = sel.GetLeptonType();


      // Fill the table

      		//  1. No Selection 
      selTable_e.Fill( datasetId, 0, weight);
      selTable_mu.Fill(datasetId, 0, weight);
      selTable_l.Fill( datasetId, 0, weight);

      		//  2. Trigger
      if (selLastStep > 0) 
	  {
		selTable_l.Fill( datasetId, 1, weight);
		if (trigger_e)  selTable_e.Fill( datasetId, 1, weight);
		if (trigger_mu) selTable_mu.Fill(datasetId, 1, weight);
      }

      		//  3. Rest of the table
	  for (unsigned int i = 2; i < (sel.GetCutList()).size() ; i++)
	  {
        if (selLastStep >= (int) i && lep==0)  selTable_e.Fill( datasetId, i, weight);
        if (selLastStep >= (int) i && lep==1)  selTable_mu.Fill(datasetId, i, weight);
        if (selLastStep >= (int) i) 		   selTable_l.Fill( datasetId, i, weight);
      }

      // Fill the histo
	  
      int selLastStep_e = 0;
      int selLastStep_mu = 0;
      if (trigger_e) 
	  {
		  selLastStep_e = 1;
      	  if (lep == 0) selLastStep_e=selLastStep;
	  }
      if (trigger_mu) 
	  {
		  selLastStep_mu = 1;
      	  if (lep == 1) selLastStep_mu=selLastStep;
	  }

      histoManager.Fill(sel, event, sel.GetMuonsForAna(), sel.GetElectronsForAna(), selLastStep_e, 0, datasetId, weight);
      histoManager.Fill(sel, event, sel.GetMuonsForAna(), sel.GetElectronsForAna(), selLastStep_mu, 1, datasetId, weight);
      //histoManager.Fill(sel, event, sel.GetMuonsForAna(), sel.GetElectronsForAna(), selLastStep, 2, datasetId, weight);

	  //dileptonBackgroundStudy(&sel,trigger_e,trigger_mu,lep,selLastStep,&histoManagerDileptonBkg,datasetId,weight);

	  // Calculate efficiency for semi-leptonic and di-leptonic
		{
 		
			const IPHCTree::NTMonteCarlo mcInfo = *(sel.GetPointer2MC());
	
  			int TMEME =  mcInfo.TMEME; 
  			int  MEME =  TMEME % 10000; 
  			int   EME =   MEME % 1000; 
  			int    ME =    EME % 100; 
 			int     E =     ME % 10;

  			int nTau      = TMEME / 10000;
  			int nMuFromTau = MEME / 1000;
  			int nElFromTau  = EME / 100;
  			int nMuon        = ME / 10;
			int nElec         = E / 1;

			bool foundVeto0 = false;
			bool foundVeto1 = false;
			bool foundVeto2 = false;
			bool foundVeto3 = false;



	  		resetDileptonData(&dileptonDataForTree);
			if (selLastStep >= 5)
			{

						TLorentzVector lepton_p;
						// Get selected muon for the check
						if (sel.GetMuonsForAna().size()==1) lepton_p = (sel.GetMuonsForAna()[0]).p4;
						else                                lepton_p = (sel.GetElectronsForAna()[0]).p4;



			if (nTau > 0)
			{
			  DEBUG_MSG << " NTauMC = " << nTau;
			  if (nMuFromTau + nElFromTau == 0) cout << " hadr ";
			  std::vector<IPHCTree::NTTau> localTaus = *(sel.GetPointer2Taus());
			  for(unsigned int i=0;i<localTaus.size();i++)
			  {
				if ( lepton_p.DeltaR(localTaus[i].p4) < 0.1) continue;
				
				//DEBUG_MSG << "tipi 0" << endl;

				if (localTaus[i].leadTrackPt <= 5) continue;
	
				//DEBUG_MSG << "tipi 1" << endl;
				
				if ((localTaus[i].ID["byLooseIsolation"]  != 1)
				 && (localTaus[i].ID["byMediumIsolation"] != 1)
				 && (localTaus[i].ID["byTightIsolation"]  != 1)
				 && (localTaus[i].ID["byLooseCombinedIsolationDeltaBetaCorr"] != 1)
				 && (localTaus[i].ID["byMediumCombinedIsolationDeltaBetaCorr"] != 1)
				 && (localTaus[i].ID["byTightCombinedIsolationDeltaBetaCorr"]  != 1)) continue;
				
				if(fabs(localTaus[i].p4.Eta()) >= 2.5)              continue;

				if(fabs(localTaus[i].p4.Eta())<1.566 && fabs(localTaus[i].p4.Eta())>1.4442) continue;
				if(localTaus[i].p4.Pt() < 5)                      continue;

				//if ( fabs( localTaus[i].vertex.Z() - sel.GetSelectedVertex()[0].p3.Z() )  > cfg.TauVertexMatchThr_ ) continue;
				if ( fabs( localTaus[i].D0)  >= 0.04 )     continue;


			  	cout << " ; found";
			  }
			cout << endl;
			}




						// Output vector
						std::vector<IPHCTree::NTPFCandidate> vetotracks = sel.GetPFCandidates();
						// Loop over pfcandidates tracks
						for(unsigned int i=0 ; i < vetotracks.size() ; i++)
						{
							if ((vetotracks[i].p4.Pt() > 10)
							&& (fabs(vetotracks[i].dz_firstGoodVertex) < 0.05))
							{
								TLorentzVector vetoTrack_p = vetotracks[i].p4;
								// Check pfcandidate doesnt match the selected lepton
								if (lepton_p.DeltaR(vetoTrack_p) > 0.1)
								{
									float thePt = vetotracks[i].p4.Pt();
								
									if (dileptonDataForTree.iso_03 > vetotracks[i].others["iso_03"]/thePt)
									{
										dileptonDataForTree.iso_03 = vetotracks[i].others["iso_03"]/thePt;
										dileptonDataForTree.bestIso03_iso_03                 = vetotracks[i].others["iso_03"]/thePt;
										dileptonDataForTree.bestIso03_iso_05_no015           = vetotracks[i].others["iso_05_no015"]/thePt;
										dileptonDataForTree.bestIso03_iso_05_no01            = vetotracks[i].others["iso_05_no01"]/thePt;
										dileptonDataForTree.bestIso03_iso_05_ntrk_01         = (float) vetotracks[i].others["ntrk_01"]; 
										dileptonDataForTree.bestIso03_iso_05_ntrk_015        = (float) vetotracks[i].others["ntrk_015"];
										dileptonDataForTree.bestIso03_iso_05_ntrk_015_pT1GeV = (float) vetotracks[i].others["ntrk_015_pT1GeV"];
									}
									if (dileptonDataForTree.iso_05_no01 > vetotracks[i].others["iso_05_no01"]/thePt)
									{
										dileptonDataForTree.iso_05_no01 = vetotracks[i].others["iso_05_no01"]/thePt;
										dileptonDataForTree.bestIso05no01_iso_03                 = vetotracks[i].others["iso_03"]/thePt;
										dileptonDataForTree.bestIso05no01_iso_05_no015           = vetotracks[i].others["iso_05_no015"]/thePt;
										dileptonDataForTree.bestIso05no01_iso_05_no01            = vetotracks[i].others["iso_05_no01"]/thePt;
										dileptonDataForTree.bestIso05no01_iso_05_ntrk_01         = (float) vetotracks[i].others["ntrk_01"]; 
										dileptonDataForTree.bestIso05no01_iso_05_ntrk_015        = (float) vetotracks[i].others["ntrk_015"];
										dileptonDataForTree.bestIso05no01_iso_05_ntrk_015_pT1GeV = (float) vetotracks[i].others["ntrk_015_pT1GeV"];
									}
									if (dileptonDataForTree.iso_05_no015 > vetotracks[i].others["iso_05_no015"]/thePt)
									{
										dileptonDataForTree.iso_05_no015 = vetotracks[i].others["iso_05_no015"]/thePt;
										dileptonDataForTree.bestIso05no015_iso_03                 = vetotracks[i].others["iso_03"]/thePt;
										dileptonDataForTree.bestIso05no015_iso_05_no015           = vetotracks[i].others["iso_05_no015"]/thePt;
										dileptonDataForTree.bestIso05no015_iso_05_no01            = vetotracks[i].others["iso_05_no01"]/thePt;
										dileptonDataForTree.bestIso05no015_iso_05_ntrk_01         = (float) vetotracks[i].others["ntrk_01"]; 
										dileptonDataForTree.bestIso05no015_iso_05_ntrk_015        = (float) vetotracks[i].others["ntrk_015"];
										dileptonDataForTree.bestIso05no015_iso_05_ntrk_015_pT1GeV = (float) vetotracks[i].others["ntrk_015_pT1GeV"];
									}
									if (dileptonDataForTree.iso_02 > vetotracks[i].others["iso_02"]/thePt) dileptonDataForTree.iso_02 = vetotracks[i].others["iso_02"]/thePt;
									if (dileptonDataForTree.iso_04 > vetotracks[i].others["iso_04"]/thePt) dileptonDataForTree.iso_04 = vetotracks[i].others["iso_04"]/thePt;
									if (dileptonDataForTree.iso_05 > vetotracks[i].others["iso_05"]/thePt) dileptonDataForTree.iso_05 = vetotracks[i].others["iso_05"]/thePt;
									if (dileptonDataForTree.iso_03_pT1GeV > vetotracks[i].others["iso_03_pT1GeV"]/thePt) 
										dileptonDataForTree.iso_03_pT1GeV = vetotracks[i].others["iso_03_pT1GeV"]/thePt;
									if (dileptonDataForTree.iso_05_no015_pT1GeV > vetotracks[i].others["iso_05_no015_pT1GeV"]/thePt)
										dileptonDataForTree.iso_05_no015_pT1GeV = vetotracks[i].others["iso_05_no015_pT1GeV"]/thePt;
									if (dileptonDataForTree.iso_05_no02 > vetotracks[i].others["iso_05_no02"]/thePt) 
										dileptonDataForTree.iso_05_no02 = vetotracks[i].others["iso_05_no02"]/thePt;
									if (dileptonDataForTree.iso_045_no015 > vetotracks[i].others["iso_045_no015"]/thePt)
										dileptonDataForTree.iso_045_no015 = vetotracks[i].others["iso_045_no015"]/thePt;
									if (dileptonDataForTree.iso_055_no015 > vetotracks[i].others["iso_055_no015"]/thePt)
										dileptonDataForTree.iso_055_no015 = vetotracks[i].others["iso_055_no015"]/thePt;

									
								}
							}
						}

						// Loop over pfcandidates tracks
						for(unsigned int i=0 ; i < vetotracks.size() ; i++)
						{
							if ((vetotracks[i].p4.Pt() + vetotracks[i].others["iso_02"] > 10)
							&& (fabs(vetotracks[i].dz_firstGoodVertex) < 0.05))
							{
								TLorentzVector vetoTrack_p = vetotracks[i].p4;
								float thePt = vetotracks[i].p4.Pt() + vetotracks[i].others["iso_02"];
								// Check pfcandidate doesnt match the selected lepton
								if (lepton_p.DeltaR(vetoTrack_p) > 0.1)
								{
									if (dileptonDataForTree.iso_05_no015_pTEqualIso02 > vetotracks[i].others["iso_05_no015"]/thePt)
										dileptonDataForTree.iso_05_no015_pTEqualIso02 = vetotracks[i].others["iso_05_no015"]/thePt;
								}
							}
						}
			
			string channelType(""); 

			// Semi-leptonic case
			if (nTau + nMuon + nElec == 1)
			{
				dileptonDataForTree.eventInfo_isSemilep = 1.0;
				channelType = string("l+jets");
				if (nTau >= 1)
					dileptonDataForTree.eventInfo_haveTau = 1.0;
			}
			// Di-leptonic case
			else if (nTau + nMuon + nElec >= 2)	
			{
				dileptonDataForTree.eventInfo_isDilep = 1.0;
				channelType = string("ll");
				if (nTau >= 1)
				{
					// Di-leptonic with a tau
					dileptonDataForTree.eventInfo_haveTau = 1.0;
					if (nMuFromTau + nElFromTau <= 0)
					{
						// Di-leptonic with a hadronic tau
						dileptonDataForTree.eventInfo_haveHadronicTau = 1.0;
						
						std::vector<IPHCTree::NTGenParticle> genParticles = mcInfo.genParticles; 
						IPHCTree::NTGenParticle theGenTau;
						int genTauIndex = -1;
						bool foundTau = false;

						for (int i = 0 ; i < genParticles.size() ; i++)
						{
							IPHCTree::NTGenParticle particle = genParticles[i];
							if (fabs(particle.id) != 15) continue; 
							if (!foundTau) { theGenTau = particle; genTauIndex = i; }
							else { DEBUG_MSG << "Warning : second tau found in the event" << endl; }
						}
					
						int PKKPPL_ = theGenTau.decayMode;
						int  KKPPL_ = PKKPPL_ % 100000;
						int   KPPL_ =  KKPPL_ % 10000;
						int    PPL_ =   KPPL_ % 1000;
						int     PL_ =    PPL_ % 100;
						int      L_ =     PL_ % 10;

						int numberOfPhotons   = PKKPPL_ / 100000;
						int numberOfKCharged  =  KKPPL_ / 10000;
						int numberOfKZero     =   KPPL_ / 1000;
						int numberOfPiCharged =    PPL_ / 100;
						int numberOfPiZero    =     PL_ / 10;
						int leptonFlavor      =      L_ / 1;

						int numberOfCharged   = numberOfKCharged + numberOfPiCharged;
						if (numberOfCharged >= 3) dileptonDataForTree.eventInfo_haveHadronicTau3 = 1.0;

								// Compute dR mean between decays products
								
								// first we find the highest pt charged decay
								int maxPt = -1.0;
								IPHCTree::NTGenParticle maxPtChargedDecay;
								for (int i = 0 ; i < genParticles.size() ; i++)
								{
									IPHCTree::NTGenParticle particle = genParticles[i];
									if (particle.motherIndex_ == genTauIndex)
									if (((abs(particle.id) == 211) || (abs(particle.id) == 321)) && (particle.p4.Pt() > maxPt))
									{
										maxPt = particle.p4.Pt();
										maxPtChargedDecay = particle;
									}
								}
							
								// then compute the sumDr
								float sumDr = 0.0;
								float DrMax = -1.0;
								int nDecay = 1;
								
								float sumDrCharged = 0.0;
								float DrMaxCharged = -1.0;
								float sumPtCharged = maxPtChargedDecay.p4.Pt();
								float pTminCharged = maxPtChargedDecay.p4.Pt();
								float pTmaxCharged = maxPtChargedDecay.p4.Pt();
								int nDecayCharged = 1;

								for (int i = 0 ; i < genParticles.size() ; i++)
								{
									IPHCTree::NTGenParticle particle = genParticles[i];
									if (particle.p4 == maxPtChargedDecay.p4) continue;
									if ((particle.motherIndex_ == genTauIndex) && (abs(particle.id) != 16))
									{
										nDecay++;
										
										// For all decays
										sumDr += maxPtChargedDecay.p4.DeltaR(particle.p4);
										
										if (maxPtChargedDecay.p4.DeltaR(particle.p4) > DrMax)
											DrMax = maxPtChargedDecay.p4.DeltaR(particle.p4);

										// ChargedDecay
										if (((abs(particle.id) == 211) || (abs(particle.id) == 321)))
										{
											nDecayCharged++;
											sumPtCharged += particle.p4.Pt();
											sumDrCharged += maxPtChargedDecay.p4.DeltaR(particle.p4);
											if (particle.p4.Pt() > pTmaxCharged)
												pTmaxCharged = particle.p4.Pt();
											if (particle.p4.Pt() < pTminCharged)
												pTminCharged = particle.p4.Pt();
											if (maxPtChargedDecay.p4.DeltaR(particle.p4) > DrMaxCharged)
												DrMaxCharged = maxPtChargedDecay.p4.DeltaR(particle.p4);
										}
									}

								}

								// Fill tree info
								dileptonDataForTree.hadronicMC_dRmeanDecay = sumDr / nDecay;
								dileptonDataForTree.hadronicMC_dRmaxDecay = DrMax;
								dileptonDataForTree.hadronicMC_dRmeanChargedDecay = sumDrCharged / nDecayCharged;
								dileptonDataForTree.hadronicMC_dRmaxChargedDecay = DrMaxCharged;
								dileptonDataForTree.hadronicMC_pTminChargedDecay = pTminCharged;
								dileptonDataForTree.hadronicMC_pTmeanChargedDecay = sumPtCharged / nDecayCharged;
								dileptonDataForTree.hadronicMC_pTmaxChargedDecay = pTmaxCharged;

					}

				}
				
			}
		
				if (channelType != string(""))
				{
						dileptonTable.Fill(channelType,"baseline",1.0);
					if (selLastStep >= 5)
					{
						dileptonTable.Fill(channelType,"MET"      ,1.0);
						dileptonTable.Fill(channelType,"vetoTest0",(int) foundVeto0);
						dileptonTable.Fill(channelType,"vetoTest1",(int) foundVeto1);
						dileptonTable.Fill(channelType,"vetoTest2",(int) foundVeto2);
						dileptonTable.Fill(channelType,"vetoTest3",(int) foundVeto3);
					}
					if (selLastStep >= 6) 
						dileptonTable.Fill(channelType,"stdVeto",1.0);
				}

				dileptonTree->Fill();
			}
	  	}

	  } // end of loop over evts
  }		// end of loop over datasets 

  dileptonTree->Print();
  dileptonTreeFile.Write();
  dileptonTreeFile.Close();
  dileptonTable.PrintTable();

  // ################################
  // #   Computations after loops   #
  // ################################
	
  cout << endl;
  cout << "   ,-------------------------------," << endl;
  cout << "   |   Compute histos and tables   |" << endl;
  cout << "   `-------------------------------`" << endl;
  cout << endl;

  // Compute and save histograms
  
  if(!doHistoManager) INFO1_MSG << "HistoManagers [skipped]" << endl;
  else
  {	
	  INFO1_MSG << "HistoManagers ..." << endl;
	  
      // Standard histograms

	  histoManager.Compute ();
	  
	  string orootfilename = rootOutputName;
	  TFile *fout = new TFile (orootfilename.c_str(), "RECREATE");
	  histoManager.Write (fout);
	  fout->Close ();  
	  
	  // DileptonBkg histograms
	  
	  histoManagerDileptonBkg.MergeChannels();
	  histoManagerDileptonBkg.DoMCStack();
	  histoManagerDileptonBkg.MergeMCDatasets();
	  histoManagerDileptonBkg.PlotsCutByCut();
	  
	  string orootfilenameDileptonBkg;
	  orootfilenameDileptonBkg = string("TTbarMETanalysis_DileptonBkgAna.root");
	  TFile *fout2 = new TFile (orootfilenameDileptonBkg.c_str(), "RECREATE");
	  histoManagerDileptonBkg.Write (fout2);
	  fout2->Close();

	  string orootfilenameDileptonBkg2D;
	  orootfilenameDileptonBkg2D = string("TTbarMETanalysis_DileptonBkgAna2D.root");
	  TFile *fout3 = new TFile (orootfilenameDileptonBkg2D.c_str(), "RECREATE");
	  histoManagerDileptonBkg.Write2D (fout3);
	  fout3->Close();

	  // Clear
  
	  histoManager.Clear ();
	  histoManagerDileptonBkg.Clear ();

	  delete fout;
	  delete fout2;
	  delete fout3;

  }
  
  // Write and compile tables

  if (!doLatexOutput) INFO1_MSG << "Tables [skipped]" << endl;
  {
	  INFO1_MSG << "Tables ..." << endl;
	
	  selTable_e.TableCalculator();
	  selTable_mu.TableCalculator();
	  selTable_l.TableCalculator();

	  makeSelectionTable(latexTableName,selTable_e,selTable_mu,selTable_l);

  }

  cout << endl;
  cout << "   ,-----------------------," << endl;
  cout << "   |   Program completed   |" << endl;
  cout << "   `-----------------------`" << endl;
  cout << endl;

  return (0);
}


