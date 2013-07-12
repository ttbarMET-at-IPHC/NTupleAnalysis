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

#include "eventListSynchro.h"
#include "miscTools.h"
          
#include "sonicScrewdriver2/interface/SonicScrewdriver.h"
#include "sonicScrewdriver2/interface/TableScrew.h"

bool noRegionSelector()  { return true; }
bool noChannelSelector() { return true; }

float* jetMassForSelector;
bool jetMassSelector() { if ((*jetMassForSelector > 60) && (*jetMassForSelector < 130)) return true; else return false; }

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
  
  INFO1_MSG << "Initializing variables..." << endl;
 
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


  
    // ############################
	// #    Config screwdriver    #
	// ############################

  // Create the screwdriver
  SonicScrewdriver mySonic;

  // Add variables
  float Wcand_pT;               mySonic.AddVariable("pT",   "p_{T}" ,               "GeV", 20,  0, 600, &Wcand_pT   );
  float Wcand_eta;              mySonic.AddVariable("eta",  "#eta",                 "",    12, -3, 3,   &Wcand_eta  );
  float Wcand_mass;             mySonic.AddVariable("mass", "Jet mass",             "GeV", 50,  0, 500, &Wcand_mass );
  float Wcand_dRgen;            mySonic.AddVariable("dR",   "#DeltaR(genW)",        "",    20,  0, 5,   &Wcand_dRgen);
  float Wcand_nSubjet;          mySonic.AddVariable("nSubjet",       "# of subjets",            "",    20,  -0.5, 19.5, &Wcand_nSubjet);
  float Wcand_nSubjet_pT5;      mySonic.AddVariable("nSubjet_pT5",   "# of subjets",            "",    20,  -0.5, 19.5, &Wcand_nSubjet_pT5);
  float Wcand_Msj12;            mySonic.AddVariable("Msj12",         "M(sj1,sj2)",              "GeV", 20,  0,    100,  &Wcand_Msj12);
  float Wcand_leadSubPt;        mySonic.AddVariable("leadSubPt",     "p_{T}(sj1)",              "GeV", 20,  0,    100,  &Wcand_leadSubPt);
  float Wcand_fPt_sj1;          mySonic.AddVariable("fPt_sj1",       "p_{T}(sj1) / p_{T}(fat)", "GeV", 11,  0,    1.1,  &Wcand_fPt_sj1);
  float Wcand_fPt_sj2;          mySonic.AddVariable("fPt_sj2",       "p_{T}(sj2) / p_{T}(fat)", "GeV",  11,  0,    1.1,  &Wcand_fPt_sj2);
  float Wcand_fPt_sj12;         mySonic.AddVariable("fPt_sj12",      "(p_{T}(sj1) + p_{T}(sj2)) / p_{T}(fat)", "GeV",    11,  0,     1.1,  &Wcand_fPt_sj12);
  jetMassForSelector = &Wcand_mass;

  // Add processClass
  mySonic.AddProcessClass("fakes",      "fakes",                "background",COLORPLOT_BLUE);
  mySonic.AddProcessClass("matched05",  "matched #DeltaR<0.5",  "background",COLORPLOT_GREEN);
  mySonic.AddProcessClass("matched01",  "matched #DeltaR<0.1",  "background",COLORPLOT_ORANGE);

  // Add region and channels
  mySonic.AddRegion("noSelection","",&noRegionSelector);
        mySonic.ScheduleVariablesForRegion("noSelection","all");
  
  mySonic.AddRegion("60 < mass(jet) < 130","",&jetMassSelector);
        mySonic.ScheduleVariablesForRegion("60 < mass(jet) < 130","all");
  
  mySonic.AddChannel("noChannel","",&noChannelSelector);


  // Create histos and schedule plots
  mySonic.Create1DHistos();
  mySonic.SchedulePlots("1DSuperpRenorm");


  // -----------------------------------------------
  // Second screwdriver for efficiencies computation
  // -----------------------------------------------

  SonicScrewdriver mySonicEff;
  SonicScrewdriver mySonicMiseff;

  mySonicEff.AddVariable("genPt",  "True p_{T} of W", "GeV", 20,  0, 600);
  mySonicEff.AddVariable("genEta", "True #eta of W",  "",    12, -3, 3  );
  mySonicMiseff.AddVariable("candPt",  "p_{T} of non-matched W-cand.", "GeV", 20,  0, 600 );
  mySonicMiseff.AddVariable("candEta", "#eta of non-matched W-cand.",  "",    12, -3, 3   );
  
  mySonicEff.AddProcessClass("genW",              "genW",                "background",COLORPLOT_GREEN);
  mySonicEff.AddProcessClass("matchedRecoW",      "matchedRecoW",        "background",COLORPLOT_BLUE );
  mySonicMiseff.AddProcessClass("unmatchedCandW", "unmatchedCandW",      "background",COLORPLOT_ORANGE);
  mySonicMiseff.AddProcessClass("unmatchedRecoW", "unmatchedRecoW",      "background",COLORPLOT_RED   );

  mySonicEff.AddRegion("noSelection","",&noRegionSelector);                     mySonicMiseff.AddRegion("noSelection","",&noRegionSelector);
        mySonicEff.ScheduleVariablesForRegion("noSelection","all");                 mySonicMiseff.ScheduleVariablesForRegion("noSelection","all");
  mySonicEff.AddRegion("60 < mass(jet) < 130","",&jetMassSelector);             mySonicMiseff.AddRegion("60 < mass(jet) < 130","",&jetMassSelector);
        mySonicEff.ScheduleVariablesForRegion("60 < mass(jet) < 130","all");        mySonicMiseff.ScheduleVariablesForRegion("60 < mass(jet) < 130","all");

  mySonicEff.AddChannel("noChannel","",&noChannelSelector);                     mySonicMiseff.AddChannel("noChannel","",&noChannelSelector);

  mySonicEff.Create1DHistos();
  mySonicMiseff.Create1DHistos();
  mySonicEff.SchedulePlots(   "TaggingEfficiency","genObjects=genW,matchedRecoObjects=matchedRecoW"                        );
  mySonicMiseff.SchedulePlots("MistaggingEfficiency","unmatchedCandObjects=unmatchedCandW,unmatchedRecoObjects=unmatchedRecoW");
  
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
        INFO1_MSG << "  Before skimming   : " << datasets[datasetId].getNSkimmedEvent() << endl;
  		INFO1_MSG << "  Afeter skimming   : " << nEvents << endl;
  		INFO1_MSG << "NEvents to run over : " << datasets[datasetId].NofEvtsToRunOver() << endl;
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

        for (unsigned int ievt = 0; ievt < datasets[datasetId].NofEvtsToRunOver(); ievt++)
        {
          if (ievt % 10000 == 0) printProgressBar(ievt,datasets[datasetId].NofEvtsToRunOver());

          // Load the event
          datasets[datasetId].eventTree()->GetEntry(ievt);
          IPHCTree::NTTransient::InitializeAfterReading(event);
          int eventId = event->general.eventNb;
          sel.LoadEvent(event);

          // Find hadronic W (assume there's only one)
          const IPHCTree::NTMonteCarlo mcInfo = *(sel.GetPointer2MC());
          vector<IPHCTree::NTGenParticle> MCParticles = mcInfo.genParticles;
          int genW_index = -1;
          for (unsigned int i = 0 ; i < MCParticles.size() ; i++)	
          {
              if (MCParticles[i].motherIndex_ == -1) continue;
              if ( (abs(MCParticles[MCParticles[i].motherIndex_].id) == 24)
               && (abs(MCParticles[i].id) <= 5) )
              {
                  genW_index = MCParticles[i].motherIndex_;
                  break;
              }
          }
          
          // Loop on W-candidate collection
          bool foundMatchedRecoW = false;
          float foundMatchedRecoW_jetMass = -1.0;

          std::vector<IPHCTree::NTJet> WCand = sel.GetHeavyTagJets();
          for (unsigned int i = 0 ; i < WCand.size() ; i++)
          {

             Wcand_pT      = WCand[i].p4.Pt();
             Wcand_eta     = WCand[i].p4.Eta();
             Wcand_mass    = WCand[i].p4.M();
             Wcand_dRgen   = WCand[i].p4.DeltaR(MCParticles[genW_index].p4);
             Wcand_nSubjet = WCand[i].subjets.size();

             Wcand_nSubjet_pT5 = 0;
             TLorentzVector subjet1;
             TLorentzVector subjet2;
             for (int j = 0 ; j < Wcand_nSubjet ; j++)
             {
                if (WCand[i].subjets[j].p4.Pt() > 5)   Wcand_nSubjet_pT5 += 1;
                     if (WCand[i].subjets[j].p4.Pt() > subjet1.Pt()) subjet1 = WCand[i].subjets[j].p4;
                else if (WCand[i].subjets[j].p4.Pt() > subjet2.Pt()) subjet2 = WCand[i].subjets[j].p4;
             }
               
             Wcand_leadSubPt = subjet1.Pt();
             Wcand_Msj12     = (subjet1 + subjet2).M();
             Wcand_fPt_sj1   = subjet1.Pt() / WCand[i].p4.Pt();
             Wcand_fPt_sj2   = subjet2.Pt() / WCand[i].p4.Pt();
             Wcand_fPt_sj12  = (subjet1.Pt() + subjet2.Pt()) / WCand[i].p4.Pt();

                  if (Wcand_dRgen < 0.1) mySonic.AutoFillProcessClass("matched01");
             else if (Wcand_dRgen < 0.5) mySonic.AutoFillProcessClass("matched05");
             else                        mySonic.AutoFillProcessClass("fakes");
                  
             // Matching for (mis)tag efficiencies
             if (Wcand_dRgen < 0.1) 
             {
                foundMatchedRecoW = true; 
                foundMatchedRecoW_jetMass = Wcand_mass;
             }
             else
             {
                float Wcand_mass_bkp = Wcand_mass;
                Wcand_mass = 80;
                mySonicMiseff.Fill("candPt", "unmatchedCandW",MCParticles[genW_index].p4.Pt());
                mySonicMiseff.Fill("candEta","unmatchedCandW",MCParticles[genW_index].p4.Eta());
                
                Wcand_mass = Wcand_mass_bkp;
                mySonicMiseff.Fill("candPt", "unmatchedRecoW",MCParticles[genW_index].p4.Pt());
                mySonicMiseff.Fill("candEta","unmatchedRecoW",MCParticles[genW_index].p4.Eta());
             }

          }

          // Fill tagging-efficiency info
          {
                Wcand_mass = 80;
                mySonicEff.Fill("genPt", "genW",MCParticles[genW_index].p4.Pt());
                mySonicEff.Fill("genEta","genW",MCParticles[genW_index].p4.Eta());

                if (foundMatchedRecoW == true)
                {
                  Wcand_mass = foundMatchedRecoW_jetMass;
                  mySonicEff.Fill("genPt", "matchedRecoW",MCParticles[genW_index].p4.Pt());
                  mySonicEff.Fill("genEta","matchedRecoW",MCParticles[genW_index].p4.Eta());
                }
          }

	  } // end of loop over evts
  
  }		// end of loop over datasets 


  cout << "   > Making plots..." << endl;
  mySonic.MakePlots();
  mySonicEff.MakePlots();
  mySonicMiseff.MakePlots();
  cout << "   > Saving plots..." << endl;
  mySonic.WritePlots("plots/firstLookAtWTagging/","");
  mySonicEff.WritePlots("plots/firstLookAtWTagging/","");
  mySonicMiseff.WritePlots("plots/firstLookAtWTagging/","");
                                                                //"exportPngAndEps=true");

  cout << endl;
  cout << "   ,-----------------------," << endl;
  cout << "   |   Program completed   |" << endl;
  cout << "   `-----------------------`" << endl;
  cout << endl;

  return (0);
}

