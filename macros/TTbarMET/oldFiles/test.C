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

#include "interface/Table.h" 
#include "interface/SonicScrewdriver.h" 
using namespace theDoctor;

bool noRegionSelector()  { return true; }
bool noChannelSelector() { return true; }

float* jetMassForSelector;
float* jetPtForSelector;
bool jetMassSelector() { if ((*jetMassForSelector > 60) && (*jetMassForSelector < 130)) return true; else return false; }
bool highPtSelector() { if (*jetPtForSelector > 200) return true; else return false; }
//bool overlapRemovalSelector() { if (*jetdRlepForSelector > 200) return true; else return false; }

// #########################################################################""
// #########################################################################""
// #########################################################################""

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
  bool doLatexOutput = false;		// False will not produce latex table
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
  
  //cout << "Initializing variables..." << endl;
 
  vector < Dataset > datasets;
  TTbarMetSelection sel;
  
  float Luminosity = 0;
  float LumiError = 0.;
  int DataType = 0;				// DataType : 0: MC - 1: Data - 2 Data & MC
  int verbosity = -1;


  reweight::LumiReWeighting *LumiWeights;
  IPHCTree::NTEvent * event = 0;

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
  cout << "Loading configuration..." << endl;
  cout << "        (config : " << xmlFileName << ")" << endl;
  
  AnalysisEnvironmentLoader anaEL (xmlFileName);
  anaEL.LoadSamples (datasets);	// now the list of datasets written in the xml file is known
  anaEL.LoadSelection (sel);	// now the parameters for the selection are given to the selection // no specific TTbarMET parameters
  anaEL.LoadGeneralInfo (DataType, Luminosity, LumiError, verbosity);

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
   
	cout << "Loading next dataset..." << endl;

    datasets[datasetId].eventTree()->SetBranchAddress ("NTEvent", &event);
    unsigned int nEvents = static_cast<unsigned int>(datasets[datasetId].eventTree()->GetEntries ());
    nEvents_tot = nEvents;

    if (verbosity > 2)
	{
		cout << endl;
  		cout << "         [ Dataset n°" << datasetId+1 << " ]" << endl;
		cout << "         " << datasets[datasetId].Name() << endl;
		cout << endl;
        cout << "  Before skimming   : " << datasets[datasetId].getNSkimmedEvent() << endl;
  		cout << "  Afeter skimming   : " << nEvents << endl;
  		cout << "NEvents to run over : " << datasets[datasetId].NofEvtsToRunOver() << endl;
		cout << endl;
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

          int totalNumberOfJets = 0;
          int totalNumberOfJetsTagged = 0;

          int totalNumberOfEvents = 0;
          int totalNumberOfEventsRejected = 0;

        for (unsigned int ievt = 0; ievt < datasets[datasetId].NofEvtsToRunOver(); ievt++)
        {
          //if (ievt % 10000 == 0) printProgressBar(ievt,datasets[datasetId].NofEvtsToRunOver());
          if (ievt % 10000 == 0) cout << "i = "<< ievt << endl;

          // Load the event
          datasets[datasetId].eventTree()->GetEntry(ievt);
          IPHCTree::NTTransient::InitializeAfterReading(event);
          int eventId = event->general.eventNb;
          sel.LoadEvent(event);

    int triggerME = 0;
    int selection_lastStep = sel.doFullSelection(&(datasets[datasetId]),string("all"),&triggerME);
    if (selection_lastStep < 5) continue;

          // Loop on W-candidate collection
          totalNumberOfEvents++;
          bool foundGoodJet = false;
          std::vector<IPHCTree::NTJet> WCand = sel.GetHeavyTagJets();
          for (unsigned int i = 0 ; i < WCand.size() ; i++)
          {

             float Wcand_pT      = WCand[i].p4.Pt();
             float Wcand_mass    = WCand[i].p4.M();
             
             if (Wcand_pT > 200)
             {
                 totalNumberOfJets++;
                 if ((Wcand_mass > 60) && (Wcand_mass < 130))
                 {
                    totalNumberOfJetsTagged++;
                    foundGoodJet = true;
                 }

             }
                 
          }
          if (foundGoodJet == false)
            totalNumberOfEventsRejected++;

          if (ievt % 1000 == 0)
          {
            cout << "mistag. eff = " << ((float)totalNumberOfJetsTagged)/((float)totalNumberOfJets) << endl;
            cout << "event reject. eff = " << ((float)totalNumberOfEventsRejected)/((float)totalNumberOfEvents) << endl;
          }


	  } // end of loop over evts
  
  }		// end of loop over datasets 


  cout << "   > Saving plots..." << endl;

  cout << endl;
  cout << "   ,-----------------------," << endl;
  cout << "   |   Program completed   |" << endl;
  cout << "   `-----------------------`" << endl;
  cout << endl;

  return (0);
}

