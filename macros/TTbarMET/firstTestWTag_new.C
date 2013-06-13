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
  // # 	 Creating tables and histomanagers   #
  // #########################################

  cout << endl;
  cout << "Creating tables and histomanagers..." << endl;


  
    // ############################
	// #    Config screwdriver    #
	// ############################

  // Create the screwdriver
  SonicScrewdriver mySonic;

  // Add variables
  float Wcand_pT;               mySonic.AddVariable("pT",   "p_{T}" ,               "GeV", 20,  0, 600, &Wcand_pT   );
  float Wcand_eta;              mySonic.AddVariable("eta",  "#eta",                 "",    12, -3, 3,   &Wcand_eta  );
  float Wcand_mass;             mySonic.AddVariable("mass", "Jet mass",             "GeV", 50,  0, 500, &Wcand_mass );
  float Wcand_dRgen;            mySonic.AddVariable("dR",   "#DeltaR(genW)",        "",    40,  0, 4,   &Wcand_dRgen);
  float Wcand_dRlep;            mySonic.AddVariable("dRlep","#DeltaR(lep)",         "",    40,  0, 4,   &Wcand_dRlep);
  float Wcand_dRbgen;           mySonic.AddVariable("dRbgen","#DeltaR(gen. b)",     "",    40,  0, 4,   &Wcand_dRbgen);
  float Wcand_dRlgen;           mySonic.AddVariable("dRlgen","#DeltaR(gen. l)",     "",    40,  0, 4,   &Wcand_dRlgen);
  float Wcand_nSubjet;          mySonic.AddVariable("nSubjet",       "# of subjets",            "",    20,  -0.5, 19.5, &Wcand_nSubjet);
  float Wcand_nSubjet_pT5;      mySonic.AddVariable("nSubjet_pT5",   "# of subjets",            "",    20,  -0.5, 19.5, &Wcand_nSubjet_pT5);
  float Wcand_Msj12;            mySonic.AddVariable("Msj12",         "M(sj1,sj2)",              "GeV", 20,  0,    100,  &Wcand_Msj12);
  float Wcand_leadSubPt;        mySonic.AddVariable("leadSubPt",     "p_{T}(sj1)",              "GeV", 20,  0,    100,  &Wcand_leadSubPt);
  float Wcand_fPt_sj1;          mySonic.AddVariable("fPt_sj1",       "p_{T}(sj1) / p_{T}(fat)", "GeV", 11,  0,    1.1,  &Wcand_fPt_sj1);
  float Wcand_fPt_sj2;          mySonic.AddVariable("fPt_sj2",       "p_{T}(sj2) / p_{T}(fat)", "GeV",  11,  0,    1.1,  &Wcand_fPt_sj2);
  float Wcand_fPt_sj12;         mySonic.AddVariable("fPt_sj12",      "(p_{T}(sj1) + p_{T}(sj2)) / p_{T}(fat)", "GeV",    11,  0,     1.1,  &Wcand_fPt_sj12);
  jetMassForSelector = &Wcand_mass;
  jetPtForSelector = &Wcand_pT;

  // Add processClass
  mySonic.AddProcessClass("fakes",      "fakes",                "background",COLORPLOT_BLUE);
  mySonic.AddProcessClass("matched05",  "matched #DeltaR<0.5",  "background",COLORPLOT_GREEN);
  mySonic.AddProcessClass("matched01",  "matched #DeltaR<0.1",  "background",COLORPLOT_ORANGE);

  // Add region and channels
  mySonic.AddRegion("noSelection","",&noRegionSelector);
        mySonic.ScheduleVariablesForRegion("noSelection","all");
  
  mySonic.AddRegion("60 < mass(jet) < 130","",&jetMassSelector);
        mySonic.ScheduleVariablesForRegion("60 < mass(jet) < 130","all");
  
  mySonic.AddRegion("pT(jet) > 200","",&highPtSelector);
        mySonic.ScheduleVariablesForRegion("pT(jet) > 200","all");
  
  mySonic.AddChannel("noChannel","",&noChannelSelector);


  // Create histos and schedule plots
  mySonic.Create1DHistos();
  mySonic.Add2DHisto("pT","mass");
  mySonic.Add2DHisto("dRlep","mass");
  mySonic.Add2DHisto("dRlgen","mass");
  mySonic.Add2DHisto("dRbgen","mass");
  mySonic.Add2DHisto("dRlgen","dRbgen");

  mySonic.SchedulePlots("1DSuperpRenorm");
  mySonic.SchedulePlots("2D");


  // -----------------------------------------------
  // Second screwdriver for efficiencies computation
  // -----------------------------------------------

  SonicScrewdriver mySonicEff;
  SonicScrewdriver mySonicMiseff;

  mySonicEff.AddVariable("genPt",  "True p_{T} of W", "GeV", 20,  0, 600);
  mySonicEff.AddVariable("genEta", "True #eta of W",  "",    12, -3, 3  );
  mySonicEff.AddVariable("flag",   "Have a matched reco W", "", 2,  -0.5,1.5);

  mySonicMiseff.AddVariable("recoPt",      "p_{T}(reco jet)", "GeV", 20,  0, 600);
  mySonicMiseff.AddVariable("selectedPt",  "p_{T}(selected jet)", "GeV", 20,  0, 600);
  mySonicMiseff.AddVariable("flag_"  ,  "Is tagged",     "",    2,-0.5, 1.5);

  mySonicEff.AddProcessClass("ttbarSemilept",  "ttbarSemilept",                "background",COLORPLOT_GREEN);
  mySonicMiseff.AddProcessClass("ttbarSemilept_",  "ttbarSemilept_",                "background",COLORPLOT_GREEN);

  mySonicEff.AddRegion("noSelection","",&noRegionSelector);                     mySonicMiseff.AddRegion("noSelection","",&noRegionSelector);
        mySonicEff.ScheduleVariablesForRegion("noSelection","all");                 mySonicMiseff.ScheduleVariablesForRegion("noSelection","all");
  mySonicEff.AddRegion("60 < mass(jet) < 130","",&jetMassSelector);             mySonicMiseff.AddRegion("60 < mass(jet) < 130","",&jetMassSelector);
        mySonicEff.ScheduleVariablesForRegion("60 < mass(jet) < 130","all");        mySonicMiseff.ScheduleVariablesForRegion("60 < mass(jet) < 130","all");

  mySonicEff.AddChannel("noChannel","",&noChannelSelector);                     mySonicMiseff.AddChannel("noChannel","",&noChannelSelector);

  mySonicEff.Create1DHistos();
  mySonicEff.Add2DHisto("genPt","flag");
  mySonicEff.SchedulePlots("1DSuperpRenorm");
  mySonicEff.SchedulePlots("2D");
  mySonicEff.SchedulePlots("2DProjectedTo1D","varX=genPt,varY=flag,projectionType=mean,labelY=TaggingEfficiency");
  
  mySonicMiseff.Create1DHistos();
  mySonicMiseff.Add2DHisto("recoPt","flag_");
  mySonicMiseff.SchedulePlots("1DSuperpRenorm");
  mySonicMiseff.SchedulePlots("1DStack");
  mySonicMiseff.SchedulePlots("2D");
  mySonicMiseff.SchedulePlots("2DProjectedTo1D","varX=recoPt,varY=flag_,projectionType=mean,labelY=MistaggingEfficiency");
  
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

        for (unsigned int ievt = 0; ievt < datasets[datasetId].NofEvtsToRunOver(); ievt++)
        {
          //if (ievt % 10000 == 0) printProgressBar(ievt,datasets[datasetId].NofEvtsToRunOver());
          if (ievt % 10000 == 0) cout << "i = "<< ievt << endl;

          // Load the event
          datasets[datasetId].eventTree()->GetEntry(ievt);
          IPHCTree::NTTransient::InitializeAfterReading(event);
          int eventId = event->general.eventNb;
          sel.LoadEvent(event);


          // Apply selection
    int triggerME = 0;
    int selection_lastStep = sel.doFullSelection(&(datasets[datasetId]),string("all"),&triggerME);
//    if (selection_lastStep < 5) continue;

    TLorentzVector lepton_p;
    int lepton_type;
    if (sel.GetMuonsForAna().size()==1)
    {
        lepton_p = (sel.GetMuonsForAna()[0]).p4;
        lepton_type = 2;
    }
    else
    {
        lepton_p = (sel.GetElectronsForAna()[0]).p4;
        lepton_type = 1;
    }

          // Find hadronic W (assume there's only one)
          const IPHCTree::NTMonteCarlo mcInfo = *(sel.GetPointer2MC());
          vector<IPHCTree::NTGenParticle> MCParticles = mcInfo.genParticles;
          TLorentzVector lepton_gen1;
          TLorentzVector lepton_gen2;
          TLorentzVector b_gen1;
          TLorentzVector b_gen2;
          int genW_index = -1;
          for (unsigned int i = 0 ; i < MCParticles.size() ; i++)	
          {
              if (MCParticles[i].motherIndex_ == -1) continue;
              if ( (abs(MCParticles[MCParticles[i].motherIndex_].id) == 24)
               && (abs(MCParticles[i].id) <= 5) )
              {
                  genW_index = MCParticles[i].motherIndex_;
                  //break;
              }
               
              if (((abs(MCParticles[i].id) == 11) || (abs(MCParticles[i].id) == 13) || (abs(MCParticles[i].id) == 15))
                  && (abs(MCParticles[MCParticles[i].motherIndex_].id) == 24))
              {
                       if (lepton_gen1 == TLorentzVector(0,0,0,0)) lepton_gen1 = MCParticles[i].p4;
                  else if (lepton_gen2 == TLorentzVector(0,0,0,0)) lepton_gen2 = MCParticles[i].p4;
              }
              if (abs(MCParticles[i].id) == 5)
              {
                       if (b_gen1 == TLorentzVector(0,0,0,0)) b_gen1 = MCParticles[i].p4;
                  else if (b_gen2 == TLorentzVector(0,0,0,0)) b_gen2 = MCParticles[i].p4;
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
             if (genW_index != -1) Wcand_dRgen   = WCand[i].p4.DeltaR(MCParticles[genW_index].p4);
             else Wcand_dRgen = 999.0;
             
             Wcand_dRlep   = WCand[i].p4.DeltaR(lepton_p);
             Wcand_dRlgen  = min(WCand[i].p4.DeltaR(lepton_gen1),WCand[i].p4.DeltaR(lepton_gen2));
             Wcand_dRbgen  = min(WCand[i].p4.DeltaR(b_gen1),WCand[i].p4.DeltaR(b_gen2));

                  if (Wcand_dRgen < 0.1) mySonic.AutoFillProcessClass("matched01");
             else if (Wcand_dRgen < 0.5) mySonic.AutoFillProcessClass("matched05");
             else                        mySonic.AutoFillProcessClass("fakes");
                  
             // Matching for (mis)tag efficiencies
             if ((Wcand_dRgen < 0.1) && (Wcand_mass > 60) && (Wcand_mass < 130) && (Wcand_pT > 200))
             {
                foundMatchedRecoW = true; 
                foundMatchedRecoW_jetMass = Wcand_mass;
             }
             
             if (Wcand_dRgen > 0.1)
             {
/*
                 cout << "   i = " << i;
                 cout << " ; pT = " << Wcand_pT;
                 cout << " ; mass = " << Wcand_mass;
                 cout << endl;
*/

                 float flag = 0;
                 if ((Wcand_mass > 60) && (Wcand_mass < 130) && (Wcand_pT > 200) && (Wcand_dRlep > 0.5))
                    flag = 1;
                 else flag = 0;

                 mySonicMiseff.Fill("recoPt","ttbarSemilept_",Wcand_pT);
                
                 if (flag == 1) mySonicMiseff.Fill("selectedPt","ttbarSemilept_",Wcand_pT);
                 else mySonicMiseff.Fill("selectedPt","ttbarSemilept_",-1.0);
                 mySonicMiseff.Fill("flag_","ttbarSemilept_",flag);
                 mySonicMiseff.Fill("recoPt","flag_","ttbarSemilept_",Wcand_pT,flag);
             }

          }

          // Fill tagging-efficiency info
          if (genW_index != -1) 
          {
                Wcand_mass = 80;

                mySonicEff.Fill("genPt", "ttbarSemilept",MCParticles[genW_index].p4.Pt());
                mySonicEff.Fill("genEta","ttbarSemilept",MCParticles[genW_index].p4.Eta());

                float flag = 0;
                if (foundMatchedRecoW == true) flag = 1;
                else flag = 0;
                

                mySonicEff.Fill("flag","ttbarSemilept",flag);
                mySonicEff.Fill("genPt","flag","ttbarSemilept",MCParticles[genW_index].p4.Pt(),flag);
          }

	  } // end of loop over evts
  
  }		// end of loop over datasets 


  cout << "   > Making plots..." << endl;
  mySonic.MakePlots();
  mySonicEff.MakePlots();
  mySonicMiseff.MakePlots();
  cout << "   > Saving plots..." << endl;
  mySonic.WritePlots("plots/firstLookAtWTagging/","");
  mySonicEff.WritePlots("plots/firstLookAtWTaggingTag/","");
  mySonicMiseff.WritePlots("plots/firstLookAtWTaggingMistag/","");
                                                                //"exportPngAndEps=true");

  cout << endl;
  cout << "   ,-----------------------," << endl;
  cout << "   |   Program completed   |" << endl;
  cout << "   `-----------------------`" << endl;
  cout << endl;

  return (0);
}

