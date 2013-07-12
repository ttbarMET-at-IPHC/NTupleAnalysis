
#include <math.h>
#include <iomanip>
#include <iostream>
#include <time.h>
using namespace std;

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

//#include <TH1.h>
//#include <TH2.h>
//#include <TCanvas.h>
//#include <TMarker.h>


#include "sonicScrewdriver/interface/SonicScrewdriver.h"
#include "sonicScrewdriver/interface/TableScrew.h"
#include "Tools/interface/LumiReweightingStandAlone.h"
using namespace reweight;

typedef struct 
{
    // Dataset
	Float_t dataset_nameHash;
    Float_t dataset_trueNumberOfEvents;
    
    // Event
    
    // Pileup
	Float_t trueNumberOfInteractions;
	Float_t numberOfVertices;

    // Multiplicity of jets, bTags
    Float_t nJets;
    Float_t nBtags;

    // Lepton info
    Float_t lep_pT;
    Float_t lep_eta;
    Float_t lep_flavor;

    // Variables of intereset 
    Float_t MET;
    Float_t MT;
    
    Float_t MET_lep_dphi;
    Float_t HT;
    Float_t HT_ratio;
    Float_t M3;
    Float_t hadrChi2;

    // Vetos
    Float_t isolatedTrackVeto;
    Float_t tauVeto;
    Float_t tauVeto_loose;
    
    // More info on tau vetoed
    
        // Kinematic
    Float_t goodTauPt;
    Float_t goodTauEta;
          
    Float_t goodTau_dRbJets;

        // Matching with MC leptons
    Float_t matchingTau_genElec;
    Float_t matchingTau_genMuon;
    Float_t matchingTau_genTau;
   
    // Monte-Carlo about leptons
    Float_t mcInfo_nLep;
    Float_t mcInfo_lep1_pdgid;
    Float_t mcInfo_lep1_eta;
    Float_t mcInfo_lep1_phi;
    Float_t mcInfo_lep1_pt;
    Float_t mcInfo_lep2_pdgid;
    Float_t mcInfo_lep2_eta;
    Float_t mcInfo_lep2_phi;
    Float_t mcInfo_lep2_pt;

} 
eventInfo;



#define FOLDER_MICROTUPLES "/opt/sbg/data/data4/cms/aaubin/tauVetoValidation_beta/store/microTuples318/"

#define ELECTRON_CHANNEL 1.0
#define MUON_CHANNEL 2.0

#define LEPTON_CHANNEL ELECTRON_CHANNEL
//#define LEPTON_CHANNEL MUON_CHANNEL


LumiReWeighting loadPileUp(string dataPileUp);
float triggerWeight(int flavor, float pt, float eta);



int main()
{
  cout << endl;
  cout << "   ,--------------------------," << endl;
  cout << "   |   Starting plot making   |" << endl;
  cout << "   `--------------------------`" << endl;
  cout << endl;

  gStyle->SetFrameLineWidth(2);

  // Create a container for the data
  eventInfo event;

  // Create the screwdriver
  SonicScrewdriver mySonic;

  // Add the variables
  
  // ###################
  // ##   Variables   ##
  // ###################
  
  mySonic.AddVariable("pTlep","p_{T}(lepton)","GeV",15,0,300,&(event.lep_pT));

  mySonic.AddVariable("Nvtx",                   "Number of good vertices",      "",40,-0.5,39.5,    &(event.numberOfVertices));
  mySonic.AddVariable("Njets",                  "# of jets",                    "",9,-0.5,8.5,      &(event.nJets)          );

  mySonic.AddVariable("MET_preTrack",           "MET",                          "GeV",30,50,500,    0,  "logY=true");
  mySonic.AddVariable("MET_CR5Track",           "MET",                          "GeV",30,50,500,    0,  "logY=true");
  mySonic.AddVariable("MET_preTau",             "MET",                          "GeV",30,50,500,    0,  "logY=true");
  mySonic.AddVariable("MET_postTau",            "MET",                          "GeV",30,50,500,    0,  "logY=true");
  mySonic.AddVariable("MET_CR5Tau",             "MET",                          "GeV",30,50,500,    0,  "logY=true");
  
  mySonic.AddVariable("MT_preTrack",            "MT",                           "GeV",17,-10,500,   0,  "logY=true");
  mySonic.AddVariable("MT_CR5Track",            "MT",                           "GeV",17,-10,500,   0,  "logY=true");
  mySonic.AddVariable("MT_preTau",              "MT",                           "GeV",17,-10,500,   0,  "logY=true");
  mySonic.AddVariable("MT_postTau",             "MT",                           "GeV",17,-10,500,   0,  "logY=true");
  mySonic.AddVariable("MT_CR5Tau",              "MT",                           "GeV",17,-10,500,   0,  "logY=true");
  
  mySonic.AddVariable("MTpeak_preTrack",        "MT_{peak} (preTrack)",         "GeV",10,10,110,    0              );
  mySonic.AddVariable("MTpeak_CR5Track",        "MT_{peak} (CR5Track)",         "GeV",10,10,110,    0              );
  mySonic.AddVariable("MTpeak_preTau",          "MT_{peak} (preTau)",           "GeV",10,10,110,    0              );
  mySonic.AddVariable("MTpeak_postTau",         "MT_{peak} (postTau)",          "GeV",10,10,110,    0              );
  mySonic.AddVariable("MTpeak_CR5Tau",          "MT_{peak} (CR5Tau)",           "GeV",10,10,110,    0              );
  
  mySonic.AddVariable("pT_goodTau",             "p_{T}(#tau-vetoed)",           "GeV",12,0,240,     0,  "logY=true");
  
  mySonic.AddVariable("MT_preTau_region1",      "MT",                           "GeV",20,0,500,     0,  "logY=true");
  mySonic.AddVariable("MT_postTau_region1",     "MT",                           "GeV",20,0,500,     0,  "logY=true");
  
  mySonic.AddVariable("MT_preTau_region2",      "MT",                           "GeV",20,0,500,     0,  "logY=true");
  mySonic.AddVariable("MT_postTau_region2",     "MT",                           "GeV",20,0,500,     0,  "logY=true");
  
  mySonic.AddVariable("DPhi_preTau",            "#Delta#phi(l,MET) (before tau veto)",  "",     10,-3.14,3.14);
  mySonic.AddVariable("HT_preTau",                     "HT (before tau veto)",                 "GeV",  20,0,500     );
  mySonic.AddVariable("HT_ratio_preTau",               "HT_{ratio} (before tau veto)",         "",     10,0,1       );
  mySonic.AddVariable("M3_preTau",                     "M_{3} (before tau veto)",              "GeV",  20,0,800     );
  //mySonic.AddVariable("hadrChi2","hadr.#Chi^{2} (before tau veto)","",10,-3.14,3.14,0);
  
  mySonic.AddVariable("DPhi_preTau_MTpeak",     "#Delta#phi(l,MET) (before tau veto, MTpeak)",  "",     10,-3.14,3.14);
  mySonic.AddVariable("HT_MTpeak",              "HT (before tau veto, MTpeak)",                 "GeV",  20,0,500     );
  mySonic.AddVariable("HT_ratio_MTpeak",        "HT_{ratio} (before tau veto, MTpeak)",         "",     10,0,1       );
  mySonic.AddVariable("M3_MTpeak",              "M_{3} (before tau veto, MTpeak)",              "GeV",  20,0,800     );


  // ######################
  // ##   ProcessClass   ##
  // ######################

  // ttbar -> l + jets
  mySonic.AddProcessClass("ttbarSemiLept",    "t#bar{t} #rightarrow l^{#pm} + jets",     "background",COLORPLOT_RED);
      mySonic.AddDataset("ttbarSemiLept",     "ttbarSemiLept",100,234*(0.67*0.33*2) );
  
  // ttbar -> ll
  mySonic.AddProcessClass("ttbarFullLept",    "t#bar{t} #rightarrow l^{+}l^{-}",         "background",COLORPLOT_CYAN);
      mySonic.AddDataset("ttbarFullLept",     "ttbarFullLept",100,234*(0.33*0.33)   );
 
  // Rare
  mySonic.AddProcessClass("rare",              "rare",                                   "background",COLORPLOT_GREEN2);
      
      mySonic.AddDataset("DY4Jets",           "rare",100,27.59    );
      mySonic.AddDataset("WGstarToLNu2E",     "rare",100,5.873    );
      mySonic.AddDataset("WWJetsTo2L2Nu",     "rare",100,5.8123   );
      mySonic.AddDataset("TBZToLL",           "rare",100,0.0114   );         
      mySonic.AddDataset("TTGJets",           "rare",100,2.166    );         
      mySonic.AddDataset("TTWJets",           "rare",100,0.23     );
      mySonic.AddDataset("TTWWJets",          "rare",100,0.002037 );
      mySonic.AddDataset("TTZJets",           "rare",100,0.2057   );
      mySonic.AddDataset("WZJetsTo3LNu",      "rare",100,1.0575   );
      mySonic.AddDataset("WZJetsTo2L2Q",      "rare",100,2.206    );
      mySonic.AddDataset("WGstarToLNu2Mu",    "rare",100,1.914    );
      mySonic.AddDataset("WGstarToLNu2Tau",   "rare",100,0.336    ); 
      mySonic.AddDataset("ZZJetsTo2L2Nu",     "rare",100,0.365    );      
      mySonic.AddDataset("WWWJets",           "rare",100,0.08058  );
      mySonic.AddDataset("WWZNoGstarJets",    "rare",100,0.05798  );
      mySonic.AddDataset("WWGJets",           "rare",100,0.528    ); 
      mySonic.AddDataset("WZZNoGstarJets",    "rare",100,0.01698  );
      mySonic.AddDataset("ZZJetsTo2L2Q",      "rare",100,2.4487   );
      mySonic.AddDataset("ZZJetsTo4L",        "rare",100,0.1769   );
      mySonic.AddDataset("ZZZNoGstarJets",    "rare",100,0.0055269);

  // Singletop
  mySonic.AddProcessClass("Singletop",         "single top",                            "background",COLORPLOT_BLUE);

      mySonic.AddDataset("Singletop_T_t",     "Singletop",100,56.4);
      mySonic.AddDataset("Singletop_Tbar_tW", "Singletop",100,11.1);
      mySonic.AddDataset("Singletop_Tbar_s",  "Singletop",100,1.76);
      mySonic.AddDataset("Singletop_T_s",     "Singletop",100,3.79);
      mySonic.AddDataset("Singletop_T_tW",    "Singletop",100,11.1);
      mySonic.AddDataset("Singletop_Tbar_t",  "Singletop",100,30.7);

  // Wjets
  mySonic.AddProcessClass("W+jets",             "W+jets",                               "background",COLORPLOT_MAGENTA);
  
      mySonic.AddDataset("W2Jets",             "W+jets",100,2159);
      mySonic.AddDataset("W3Jets",             "W+jets",100,640);
      mySonic.AddDataset("W4Jets",             "W+jets",100,264);


  // Data
  string lep_channel;
  LumiReWeighting pileUpWeights;

  mySonic.AddProcessClass("data","data","data",COLORPLOT_PINK);
  if (LEPTON_CHANNEL == ELECTRON_CHANNEL)
  {  
        lep_channel = "e-channel";
  
            /*
            mySonic.AddDataset("Elec_A", "data",1,386);
            mySonic.AddDataset("Elec_B", "data",1,841);
            mySonic.AddDataset("Elec_C1","data",1,426);
            */ 

            /* New prod */
            mySonic.AddDataset("Elec_A",   "data",1, 764);
            mySonic.AddDataset("Elec_B-1", "data",1, 2574/3.0);
            mySonic.AddDataset("Elec_B-2", "data",1, 2574/3.0);
            mySonic.AddDataset("Elec_B-3", "data",1, 2574/3.0);
            mySonic.AddDataset("Elec_C1",  "data",1, 478);
            mySonic.AddDataset("Elec_C2-1","data",1, 2430/2.0);
            mySonic.AddDataset("Elec_C2-2","data",1, 2430/2.0);
            mySonic.AddDataset("Elec_D-1", "data",1, 2092/2.0);
            mySonic.AddDataset("Elec_D-2", "data",1, 2092/2.0);
            /**/

        // Pile-up
        //pileUpWeights = loadPileUp("pileUpElecABC1.root");
       pileUpWeights = loadPileUp("pileUpElecABC1C2D.root");
  }
  else if (LEPTON_CHANNEL == MUON_CHANNEL)
  {  
        lep_channel = "#mu-channel";

            /*
            mySonic.AddDataset("Muon_A",   "data",1,669   );
            mySonic.AddDataset("Muon_B-1", "data",1,3040/2);
            mySonic.AddDataset("Muon_B-2", "data",1,3040/2);
            mySonic.AddDataset("Muon_C1",  "data",1,482   );
            */

            /* New prod */
            mySonic.AddDataset("Muon_A",   "data",1,499     );
            mySonic.AddDataset("Muon_B-1", "data",1,1783/2.0);
            mySonic.AddDataset("Muon_B-2", "data",1,1783/2.0);
            mySonic.AddDataset("Muon_C1",  "data",1,482     );
            mySonic.AddDataset("Muon_C2-1","data",1,2765/3.0);
            mySonic.AddDataset("Muon_C2-2","data",1,2765/3.0);
            mySonic.AddDataset("Muon_C2-3","data",1,2765/3.0);
            mySonic.AddDataset("Muon_D-1", "data",1,2802/3.0);
            mySonic.AddDataset("Muon_D-2", "data",1,2802/3.0);
            mySonic.AddDataset("Muon_D-3", "data",1,2802/3.0);
            /**/
        
        // Pile-up
        //pileUpWeights = loadPileUp("pileUpMuonABC1.root");
        pileUpWeights = loadPileUp("pileUpMuonABC1C2D.root");
  }
  
 
  mySonic.Add1DHistoAll();


  mySonic.SchedulePlots("MCDataComparison");
  mySonic.SchedulePlots("1DSuperpRenorm");

  // ###################################
  // ## Loop over datasets and events ##
  // ###################################

  vector<string> datasetsList;
  mySonic.GetDatasetList(&datasetsList);

  for (unsigned int i = 0 ; i < datasetsList.size() ; i++)
  {
      string currentDataset = datasetsList[i];
      string currentProcessClass = mySonic.GetProcessClass(currentDataset);

      TFile f((string(FOLDER_MICROTUPLES)+currentDataset+".root").c_str());
      TTree* theTree = (TTree*) f.Get("microTuple"); 
      theTree->SetBranchAddress("events",&event);
      
      cout << "Running on dataset : " << currentDataset << " (" << theTree->GetEntries();
     
      // ##############################
      //      Loop over events
      // ##############################
      
      float weight_lumi = 0.0;
      
      for (int j = 0 ; j < theTree->GetEntries() ; j++)
      {
        theTree->GetEntry(j);
        
        if (j == 0) 
        {    
           mySonic.SetTrueNumberOfEvents(currentDataset,event.dataset_trueNumberOfEvents);
           cout << " events from " << event.dataset_trueNumberOfEvents << " before skimming)" << endl;
           weight_lumi = mySonic.GetDatasetLumiWeight(currentDataset);
        }
        
        // ##############################
        // ## Apply selection on event ##
        // ##############################
        
        if (event.lep_flavor != LEPTON_CHANNEL) continue;
        if (event.nJets < 4.0) continue;
        if (event.nBtags < 1.0) continue;

        // ######################################
        // ##   Compute weight for the event   ##
        // ######################################
       
        // Pile-up
        float weight_pu;
        if (mySonic.GetProcessClassType(currentProcessClass) == "background")
            weight_pu = pileUpWeights.ITweight((int) event.trueNumberOfInteractions);
        else //if (mySonic.GetProcessClassType(currentProcessClass) == "data")
            weight_pu = 1.0;

        // Trigger
        float weight_trigger;
        if (mySonic.GetProcessClassType(currentProcessClass) == "background")
            weight_trigger = triggerWeight(event.lep_flavor,event.lep_pT,event.lep_eta);
        else
            weight_trigger = 1.0;

        // Total weight
        float weight = weight_lumi * weight_pu * weight_trigger;
        
        // #########################
        // ##   Fill histograms   ##
        // #########################

        // Get region flag
        
        // MT peak
        bool inMTpeak = false;
        if ((event.MT >= 50) && (event.MT <= 80))
            inMTpeak = true;

        // Region 1 and 2
        bool inRegion1 = false;
        bool inRegion2 = false;
//        if       ((event.nJets <= 3.0)  && (abs(event.MET_lep_dphi) > 2) && (event.HT_ratio < 0.2)) inRegion1 = true;
//        else if  ((event.nJets >= 4.0) &&         inMTpeak              ) inRegion2 = true;

        if (inMTpeak) mySonic.Fill("MTpeak_preTrack",currentProcessClass,event.MT,weight);
        
        if (event.isolatedTrackVeto != 1.0) continue;

        // Hard-coded values for SFpre_track and SFpost_track

        #ifdef APPLY_SCALE_FACTORS_TRACK
            if (LEPTON_CHANNEL == ELECTRON_CHANNEL)
            {
                if (currentProcessClass == "ttbarFullLept")      weight *= 1.01;
                else if ((currentProcessClass == "ttbarSemiLept") 
                      || (currentProcessClass == "W+jets")
                      || (currentProcessClass == "Singletop"))   weight *= 1.0;
            }
            else if (LEPTON_CHANNEL == MUON_CHANNEL)
            {
                if (currentProcessClass == "ttbarFullLept")      weight *= 1.09;
                else if ((currentProcessClass == "ttbarSemiLept") 
                      || (curren:tProcessClass == "W+jets")
                      || (currentProcessClass == "Singletop"))   weight *= 1.08;
            }
        #endif

        mySonic.AutoFillProcessClass(currentProcessClass,weight);
        
        mySonic.Fill("MET_preTau",currentProcessClass,event.MET,weight);
        mySonic.Fill("MT_preTau" ,currentProcessClass,event.MT,weight);
       
        mySonic.Fill("DPhi_preTau",currentProcessClass,event.MET_lep_dphi,weight);
        mySonic.Fill("HT_preTau",currentProcessClass,event.HT,weight);
        mySonic.Fill("HT_ratio_preTau",currentProcessClass,event.HT_ratio,weight);
        mySonic.Fill("M3_preTau",currentProcessClass,event.M3,weight);

        if (inMTpeak)
        {
            mySonic.Fill("MTpeak_preTau",currentProcessClass,event.MT,weight);
            
            mySonic.Fill("DPhi_preTau_MTpeak",currentProcessClass,event.MET_lep_dphi,weight);
            mySonic.Fill("HT_MTpeak",currentProcessClass,event.HT,weight);
            mySonic.Fill("HT_ratio_MTpeak",currentProcessClass,event.HT_ratio,weight);
            mySonic.Fill("M3_MTpeak",currentProcessClass,event.M3,weight);
        }
        
        if (inRegion1) mySonic.Fill("MT_preTau_region1",currentProcessClass,event.MT,weight);
        if (inRegion2) mySonic.Fill("MT_preTau_region2",currentProcessClass,event.MT,weight);

        if (event.tauVeto == 1.0)
        {
            mySonic.Fill("MET_CR5Tau",currentProcessClass,event.MET,weight);
            mySonic.Fill("MT_CR5Tau" ,currentProcessClass,event.MT,weight);
            if (inMTpeak) mySonic.Fill("MTpeak_CR5Tau" ,currentProcessClass,event.MT,weight);

        
            mySonic.Fill("pT_goodTau",currentProcessClass,event.goodTauPt,weight);

                // Match
                float minDrMatch = 999.0;
                string channelMatch = "fake";

                if ((event.matchingTau_genElec < 0.15) && (event.matchingTau_genElec < minDrMatch)) 
                { channelMatch = "matchElec"; minDrMatch = event.matchingTau_genElec; }
                if ((event.matchingTau_genMuon < 0.15) && (event.matchingTau_genMuon < minDrMatch)) 
                { channelMatch = "matchMuon"; minDrMatch = event.matchingTau_genMuon; }
                if ((event.matchingTau_genTau < 0.15) && (event.matchingTau_genTau < minDrMatch)) 
                { channelMatch = "matchTau"; minDrMatch = event.matchingTau_genTau; }
              
        }
        else
        {
            // Hard-coded value for SFpre_tau and SFpost_tau
/*
             if (currentProcessClass == "ttbarFullLept") weight *= 1;
        else if ((currentProcessClass == "ttbarSemiLept") 
              || (currentProcessClass == "W+jets")
              || (currentProcessClass == "Singletop"))   weight *= 0.99716;
*/
 
            mySonic.Fill("MET_postTau",currentProcessClass,event.MET,weight);
            mySonic.Fill("MT_postTau" ,currentProcessClass,event.MT,weight);

            if (inMTpeak)  mySonic.Fill("MTpeak_postTau",currentProcessClass,event.MT,weight);
            if (inRegion1) mySonic.Fill("MT_postTau_region1",currentProcessClass,event.MT,weight);
            if (inRegion2) mySonic.Fill("MT_postTau_region2",currentProcessClass,event.MT,weight);
         }
      }
     
      f.Close();
  }

  // #######################
  // ## Tables for yields ##
  // #######################

  vector<string> processClassesList;
  mySonic.GetProcessClassList(&processClassesList);
  processClassesList.push_back("total MC");
  processClassesList.push_back("data");

  std::vector<std::string> region;
  region.push_back("pre-track-veto");
  region.push_back("pre-tau-veto");
  region.push_back("post-veto");
  region.push_back("CR5tau");

  TableScrew yieldsMTPeak(region,processClassesList);

  FigureScrew tmpYield;

  FigureScrew yieldTotalPreTrack = 0.0;
  FigureScrew yieldTotalPreTau = 0.0;
  FigureScrew yieldTotalPost = 0.0;
  FigureScrew yieldTotalCR5tau = 0.0;

  for (int i = 0 ; i < processClassesList.size() ; i++)
  {
     if (processClassesList[i] == "total MC") continue;

     tmpYield = mySonic.GetYieldAndError("MTpeak_preTrack",processClassesList[i]);
     yieldsMTPeak.Set("pre-track-veto",processClassesList[i],tmpYield);
     if (processClassesList[i] != "data") yieldTotalPreTrack += tmpYield;

     tmpYield = mySonic.GetYieldAndError("MTpeak_preTau",processClassesList[i]);
     yieldsMTPeak.Set("pre-tau-veto",processClassesList[i],tmpYield);
     if (processClassesList[i] != "data") yieldTotalPreTau += tmpYield;

     tmpYield = mySonic.GetYieldAndError("MTpeak_postTau",processClassesList[i]);
     yieldsMTPeak.Set("post-veto",processClassesList[i],tmpYield);
     if (processClassesList[i] != "data") yieldTotalPost += tmpYield;
     
     tmpYield = mySonic.GetYieldAndError("MTpeak_CR5Tau",processClassesList[i]);
     yieldsMTPeak.Set("CR5tau",processClassesList[i],tmpYield);
     if (processClassesList[i] != "data") yieldTotalCR5tau += tmpYield;
  }
 
  yieldsMTPeak.Set("pre-track-veto","total MC",yieldTotalPreTrack);
  yieldsMTPeak.Set("pre-tau-veto","total MC",yieldTotalPreTau);
  yieldsMTPeak.Set("post-veto","total MC",yieldTotalPost);
  yieldsMTPeak.Set("CR5tau","total MC",yieldTotalCR5tau);
 
  cout.precision(5);
  yieldsMTPeak.PrintTable("withError");
 
  FigureScrew SFpre_track = (yieldsMTPeak.Get("pre-track-veto","data")     - yieldsMTPeak.Get("pre-track-veto","rare"))
                    / (yieldsMTPeak.Get("pre-track-veto","total MC") - yieldsMTPeak.Get("pre-track-veto","rare"));

  FigureScrew SFpost_track=    (yieldsMTPeak.Get("pre-tau-veto","data")     - yieldsMTPeak.Get("pre-tau-veto","rare") - SFpre_track * yieldsMTPeak.Get("pre-tau-veto","ttbarFullLept"))     
                      /  (yieldsMTPeak.Get("pre-tau-veto","total MC") - yieldsMTPeak.Get("pre-tau-veto","rare") - yieldsMTPeak.Get("pre-tau-veto","ttbarFullLept"));

  FigureScrew SFpre_tau = (yieldsMTPeak.Get("pre-tau-veto","data")     - yieldsMTPeak.Get("pre-tau-veto","rare")) 
                 / ((yieldsMTPeak.Get("pre-tau-veto","total MC") - yieldsMTPeak.Get("pre-tau-veto","rare") - yieldsMTPeak.Get("pre-tau-veto","ttbarFullLept"))*SFpost_track  
                    + yieldsMTPeak.Get("pre-tau-veto","ttbarFullLept") * SFpre_track
                   );
  
  FigureScrew SFpost_tau = (yieldsMTPeak.Get("post-veto","data")     - yieldsMTPeak.Get("post-veto","rare") - SFpre_track * SFpre_tau * yieldsMTPeak.Get("post-veto","ttbarFullLept"))    
                 /  ((yieldsMTPeak.Get("post-veto","total MC") - yieldsMTPeak.Get("post-veto","rare") - yieldsMTPeak.Get("post-veto","ttbarFullLept"))*SFpost_track );
  
   cout << "SFpre_track  = " << SFpre_track.Print()  << endl;
   cout << "SFpost_track = " << SFpost_track.Print() << endl;
   cout << endl;
   cout << "SFpre_tau = " << SFpre_tau.Print()   << endl;
   cout << "SFpost_tau = " << SFpost_tau.Print() << endl;

  FigureScrew remaining =  (yieldsMTPeak.Get("post-veto","data")         - yieldsMTPeak.Get("post-veto","rare")) 
                 /  (( yieldsMTPeak.Get("post-veto","total MC")    - yieldsMTPeak.Get("post-veto","rare") - yieldsMTPeak.Get("post-veto","ttbarFullLept"))*SFpost_track*SFpost_tau 
                     + yieldsMTPeak.Get("post-veto","ttbarFullLept")*SFpre_track*SFpre_tau
                     );

   cout << endl;
   cout << "remaining = " << remaining.Print() << endl;
 




   // Apply scale factors


      // Apply tayId-scaleFactor on dilepton
      FigureScrew tauID_scaleFactor(1.0,0.07);
      mySonic.ApplyScaleFactor("MTpeak_postTau","ttbarFullLept",tauID_scaleFactor);
      mySonic.ApplyScaleFactor("MTpeak_CR5Tau","ttbarFullLept",  tauID_scaleFactor);

      // Apply pre-track SF on dilepton
      mySonic.ApplyScaleFactor("MTpeak_postTau","ttbarFullLept",SFpre_track);
      mySonic.ApplyScaleFactor("MTpeak_CR5Tau","ttbarFullLept",  SFpre_track);
  
      // Apply post-track SF on 1-lep backgrounds
      mySonic.ApplyScaleFactor("MTpeak_postTau","ttbarSemiLept",SFpost_track);
      mySonic.ApplyScaleFactor("MTpeak_postTau","Singletop",    SFpost_track);
      mySonic.ApplyScaleFactor("MTpeak_postTau","W+jets",       SFpost_track);
      mySonic.ApplyScaleFactor("MTpeak_CR5Tau",  "ttbarSemiLept",SFpost_track);
      mySonic.ApplyScaleFactor("MTpeak_CR5Tau",  "Singletop",    SFpost_track);
      mySonic.ApplyScaleFactor("MTpeak_CR5Tau",  "W+jets",       SFpost_track);
 
      // Here : scale factor for tau's


   //
   //   Method 2
   //
  
  std::vector<std::string> reducedProcessClassList;
  reducedProcessClassList.push_back("ll");
  reducedProcessClassList.push_back("others");
  reducedProcessClassList.push_back("total MC");
  reducedProcessClassList.push_back("data");

  std::vector<std::string> regionMethod2;
  regionMethod2.push_back("pre-veto1");
  regionMethod2.push_back("post-veto1");
  regionMethod2.push_back("eff1");
  regionMethod2.push_back("pre-veto2");
  regionMethod2.push_back("post-veto2");
  regionMethod2.push_back("eff2");

  TableScrew tableMethod2(regionMethod2,reducedProcessClassList);

  string region1pre("MT_preTau_region1");
  string region1post("MT_postTau_region1");
  
  string region2pre("MT_preTau_region2");
  string region2post("MT_postTau_region2");

  FigureScrew yieldPre1 = 0;
  FigureScrew yieldPost1 = 0;
  FigureScrew yieldPre2 = 0;
  FigureScrew yieldPost2 = 0;
  
  FigureScrew yieldPre1TotMC = 0;
  FigureScrew yieldPost1TotMC = 0;
  FigureScrew yieldPre2TotMC = 0;
  FigureScrew yieldPost2TotMC = 0;

  tmpYield = mySonic.GetYieldAndError(region1pre, "ttbarFullLept"); tableMethod2.Set("pre-veto1" ,"ll",tmpYield); yieldPre1TotMC +=tmpYield;
  tmpYield = mySonic.GetYieldAndError(region1post,"ttbarFullLept"); tableMethod2.Set("post-veto1","ll",tmpYield); yieldPost1TotMC+=tmpYield;
  tmpYield = mySonic.GetYieldAndError(region2pre, "ttbarFullLept"); tableMethod2.Set("pre-veto2" ,"ll",tmpYield); yieldPre2TotMC +=tmpYield;
  tmpYield = mySonic.GetYieldAndError(region2post,"ttbarFullLept"); tableMethod2.Set("post-veto2","ll",tmpYield); yieldPost2TotMC+=tmpYield;

  tmpYield = mySonic.GetYieldAndError(region1pre, "ttbarSemiLept") 
           + mySonic.GetYieldAndError(region1pre, "Singletop") 
           + mySonic.GetYieldAndError(region1pre, "W+jets") 
           + mySonic.GetYieldAndError(region1pre, "rare");     tableMethod2.Set("pre-veto1" ,"others",tmpYield); yieldPre1TotMC +=tmpYield;
  tmpYield = mySonic.GetYieldAndError(region1post,"ttbarSemiLept") 
           + mySonic.GetYieldAndError(region1post,"Singletop") 
           + mySonic.GetYieldAndError(region1post,"W+jets") 
           + mySonic.GetYieldAndError(region1post,"rare");    tableMethod2.Set("post-veto1","others",tmpYield); yieldPost1TotMC +=tmpYield;
  tmpYield = mySonic.GetYieldAndError(region2pre, "ttbarSemiLept") 
           + mySonic.GetYieldAndError(region2pre, "Singletop") 
           + mySonic.GetYieldAndError(region2pre, "W+jets") 
           + mySonic.GetYieldAndError(region2pre, "rare");     tableMethod2.Set("pre-veto2" ,"others",tmpYield); yieldPre2TotMC +=tmpYield;
  tmpYield = mySonic.GetYieldAndError(region2post,"ttbarSemiLept") 
           + mySonic.GetYieldAndError(region2post,"Singletop") 
           + mySonic.GetYieldAndError(region2post,"W+jets") 
           + mySonic.GetYieldAndError(region2post,"rare");    tableMethod2.Set("post-veto2","others",tmpYield); yieldPost2TotMC +=tmpYield;

  tableMethod2.Set("pre-veto1" ,"total MC",yieldPre1TotMC);
  tableMethod2.Set("post-veto1","total MC",yieldPost1TotMC);
  tableMethod2.Set("pre-veto2" ,"total MC",yieldPre2TotMC);
  tableMethod2.Set("post-veto2","total MC",yieldPost2TotMC);

  tmpYield = mySonic.GetYieldAndError(region1pre, "data"); tableMethod2.Set("pre-veto1" ,"data",tmpYield);
  tmpYield = mySonic.GetYieldAndError(region1post,"data"); tableMethod2.Set("post-veto1","data",tmpYield);
  tmpYield = mySonic.GetYieldAndError(region2pre, "data"); tableMethod2.Set("pre-veto2" ,"data",tmpYield);
  tmpYield = mySonic.GetYieldAndError(region2post,"data"); tableMethod2.Set("post-veto2","data",tmpYield);

  for (unsigned int i = 0 ; i  < reducedProcessClassList.size() ; i++)
  {
    tableMethod2.Set("eff1",reducedProcessClassList[i],  tableMethod2.Get("post-veto1",reducedProcessClassList[i]) 
                                                       / tableMethod2.Get("pre-veto1" ,reducedProcessClassList[i]) );
    tableMethod2.Set("eff2",reducedProcessClassList[i],  tableMethod2.Get("post-veto2",reducedProcessClassList[i])
                                                       / tableMethod2.Get("pre-veto2" ,reducedProcessClassList[i]) );
  }

  tableMethod2.PrintTable();

  FigureScrew SF_1 = tableMethod2.Get("pre-veto1","data") / tableMethod2.Get("pre-veto1","total MC");
  FigureScrew SF_2 = tableMethod2.Get("pre-veto2","data") / tableMethod2.Get("pre-veto2","total MC");

  //FigureScrew N_1_data = tableMethod2.Get("post-veto1","data"); 
  //FigureScrew N_2_data = tableMethod2.Get("post-veto2","data"); 
 
  FigureScrew N_1_data = tableMethod2.Get("post-veto1","total MC") * SF_1; 
  FigureScrew N_2_data = tableMethod2.Get("post-veto2","total MC") * SF_2; 
  
  /*
  FigureScrew SF_1(1.0,0.0);
  FigureScrew SF_2(1.0,0.0);

  FigureScrew N_1_data = (tableMethod2.Get("pre-veto1","others") * 0.97 + tableMethod2.Get("pre-veto1","ll") * 0.85 ) * SF_1;
  FigureScrew N_2_data = (tableMethod2.Get("pre-veto2","others") * 0.97 + tableMethod2.Get("pre-veto2","ll") * 0.85 ) * SF_2; 
  */

  FigureScrew N_1_ll  = tableMethod2.Get("pre-veto1","ll")     * SF_1; 
  FigureScrew N_1_oth = tableMethod2.Get("pre-veto1","others") * SF_1; 
  
  FigureScrew N_2_ll  = tableMethod2.Get("pre-veto2","ll")     * SF_2; 
  FigureScrew N_2_oth = tableMethod2.Get("pre-veto2","others") * SF_2; 

  cout << "(1) N_ll  = " << N_1_ll.Print() << endl;
  cout << "(1) N_oth = " << N_1_oth.Print() << endl;

  cout << endl;

  cout << "(1) fraction ll = " << (N_1_ll/(N_1_ll + N_1_oth)).Print() << endl;
  cout << "(1) fraction others = " << (N_1_oth/(N_1_ll + N_1_oth)).Print() << endl;
  
  cout << endl;
  
  cout << "(2) N_ll  = " << N_2_ll.Print() << endl;
  cout << "(2) N_oth = " << N_2_oth.Print() << endl;
  
  cout << endl;
  
  cout << "(2) fraction ll = " << (N_2_ll/(N_2_ll + N_2_oth)).Print() << endl;
  cout << "(2) fraction others = " << (N_2_oth/(N_2_ll + N_2_oth)).Print() << endl;

  FigureScrew eps_others_fromData = ( N_2_data - (N_2_ll / N_1_ll)   * N_1_data ) / ( N_2_oth - (N_2_ll / N_1_ll)   * N_1_oth );
  FigureScrew eps_ll_fromData     = ( N_2_data - (N_2_oth / N_1_oth) * N_1_data ) / ( N_2_ll -  (N_2_oth / N_1_oth) * N_1_ll  );

  cout << "-----------------" << endl;
  
  cout << "(1) N_data = " << N_1_data.Print() << endl;
  cout << "(2) N_data = " << N_2_data.Print() << endl;
  
  cout << endl;
  cout << "eps_others_data = " << eps_others_fromData.Print() << endl;
  cout << "eps_ll_data = " << eps_ll_fromData.Print() << endl;
  cout << endl;

  cout << "SF_others(2) = " << (eps_others_fromData/tableMethod2.Get("eff2","others")).Print() << endl;
  cout << "SF_ll(2) = " << (eps_ll_fromData/tableMethod2.Get("eff2","ll")).Print() << endl;

  // #######################
  // ##   Compute plots   ##
  // #######################

  

  std::ostringstream lumiInString;
  lumiInString << mySonic.GetLumi();
  string infoText("CMS Internal      #sqrt{s} = 8 TeV, ");
  infoText += string("L = ")+lumiInString.str()+string(" pb^{-1}, ");
  infoText += lep_channel;

  cout << endl;
  cout << "   > Making plots..." << endl;
  mySonic.MakePlots();
  cout << "   > Saving plots..." << endl;
  //mySonic.WritePlots("plots/dataMCComparison.root",infoText,"exportPngAndEps=true");
  mySonic.WritePlots("plots/dataMCComparison.root",infoText);
  cout << endl;
  
  
  cout << endl;
  cout << "   ,-----------------------," << endl;
  cout << "   |   Program completed   |" << endl;
  cout << "   `-----------------------`" << endl;
  cout << endl;


 



  return (0);

}




LumiReWeighting loadPileUp(string dataPileUp)
{

 Double_t Summer2012_S10[60] = {2.560E-06, 5.239E-06, 1.420E-05, 5.005E-05, 1.001E-04, 2.705E-04, 1.999E-03,
 6.097E-03, 1.046E-02, 1.383E-02, 1.685E-02, 2.055E-02, 2.572E-02, 3.262E-02, 4.121E-02, 4.977E-02, 5.539E-02,
 5.725E-02, 5.607E-02, 5.312E-02, 5.008E-02, 4.763E-02, 4.558E-02, 4.363E-02, 4.159E-02, 3.933E-02, 3.681E-02,
 3.406E-02, 3.116E-02, 2.818E-02, 2.519E-02, 2.226E-02, 1.946E-02, 1.682E-02, 1.437E-02, 1.215E-02, 1.016E-02,
 8.400E-03, 6.873E-03, 5.564E-03, 4.457E-03, 3.533E-03, 2.772E-03, 2.154E-03, 1.656E-03, 1.261E-03, 9.513E-04,
 7.107E-04, 5.259E-04, 3.856E-04, 2.801E-04, 2.017E-04, 1.439E-04, 1.017E-04, 7.126E-05, 4.948E-05, 3.405E-05,
 2.322E-05, 1.570E-05, 5.005E-06};

 vector<float> mc_vect;
 vector<float> data_vect;


 TFile *filepuest = new TFile(dataPileUp.c_str(),"READ");
 TH1F* npu_dat = (TH1F*) filepuest->Get("pileup");

 float norm_MC = 0.0;
 float norm_data = 0.0;

 for(int i=0; i<60; ++i) 
 {
    mc_vect.push_back(Summer2012_S10[i]);
    norm_MC += Summer2012_S10[i];
 }
 for (int i=0; i<60; i++)
 {
    data_vect.push_back(npu_dat->GetBinContent(i+1));
    norm_data += npu_dat->GetBinContent(i+1);
 }
  
  filepuest->Close();

  return reweight::LumiReWeighting(mc_vect,data_vect);
}




float triggerWeight(int flavor, float pt, float eta)
{

  //electron efficiencies
  if (flavor == 1) 
  {
    if ( fabs(eta)<1.5) 
    {
      if ( pt>=20 && pt<22 )   return 0.00;
      if ( pt>=22 && pt<24 )   return 0.00;
      if ( pt>=24 && pt<26 )   return 0.00;
      if ( pt>=26 && pt<28 )   return 0.08;
      if ( pt>=28 && pt<30 )   return 0.61;
      if ( pt>=30 && pt<32 )   return 0.86;
      if ( pt>=32 && pt<34 )   return 0.88;
      if ( pt>=34 && pt<36 )   return 0.90;
      if ( pt>=36 && pt<38 )   return 0.91;
      if ( pt>=38 && pt<40 )   return 0.92;
      if ( pt>=40 && pt<50 )   return 0.94;
      if ( pt>=50 && pt<60 )   return 0.95;
      if ( pt>=60 && pt<80 )   return 0.96;
      if ( pt>=80 && pt<100 )  return 0.96;
      if ( pt>=100 && pt<150 ) return 0.96;
      if ( pt>=150 && pt<200 ) return 0.97;
      if ( pt>=200 )           return 0.97;
    } 
    else if ( fabs(eta)>=1.5 && fabs(eta)<2.1) 
    {
      if ( pt>=20 && pt<22 )   return 0.00;
      if ( pt>=22 && pt<24 )   return 0.00;
      if ( pt>=24 && pt<26 )   return 0.02;
      if ( pt>=26 && pt<28 )   return 0.18;
      if ( pt>=28 && pt<30 )   return 0.50;
      if ( pt>=30 && pt<32 )   return 0.63;
      if ( pt>=32 && pt<34 )   return 0.68;
      if ( pt>=34 && pt<36 )   return 0.70;
      if ( pt>=36 && pt<38 )   return 0.72;
      if ( pt>=38 && pt<40 )   return 0.74;
      if ( pt>=40 && pt<50 )   return 0.76;
      if ( pt>=50 && pt<60 )   return 0.77;
      if ( pt>=60 && pt<80 )   return 0.78;
      if ( pt>=80 && pt<100 )  return 0.80;
      if ( pt>=100 && pt<150 ) return 0.79;
      if ( pt>=150 && pt<200 ) return 0.76;
      if ( pt>=200 )           return 0.81;
    }
  } 
  //muon efficiencies
  else if (flavor == 2) 
  {
    if ( fabs(eta)<0.8 ) 
    {
      if (pt>=20 && pt<22)   return 0.00;     
      if (pt>=22 && pt<24)   return 0.03;      
      if (pt>=24 && pt<26)   return 0.87;
      if (pt>=26 && pt<28)   return 0.90;
      if (pt>=28 && pt<30)   return 0.91;
      if (pt>=30 && pt<32)   return 0.91;
      if (pt>=32 && pt<34)   return 0.92;
      if (pt>=34 && pt<36)   return 0.93;
      if (pt>=36 && pt<38)   return 0.93;
      if (pt>=38 && pt<40)   return 0.93;
      if (pt>=40 && pt<50)   return 0.94;
      if (pt>=50 && pt<60)   return 0.95;
      if (pt>=60 && pt<80)   return 0.95;
      if (pt>=80 && pt<100)  return 0.94;
      if (pt>=100 && pt<150) return 0.94;
      if (pt>=150 && pt<200) return 0.93;
      if (pt>=200)           return 0.92;
    } 
    else if ( fabs(eta)>=0.8 && fabs(eta)<1.5 ) 
    {
      if (pt>=20 && pt<22)   return 0.00;
      if (pt>=22 && pt<24)   return 0.05;
      if (pt>=24 && pt<26)   return 0.78;
      if (pt>=26 && pt<28)   return 0.81;
      if (pt>=28 && pt<30)   return 0.81;
      if (pt>=30 && pt<32)   return 0.81;
      if (pt>=32 && pt<34)   return 0.82;
      if (pt>=34 && pt<36)   return 0.82;
      if (pt>=36 && pt<38)   return 0.83;
      if (pt>=38 && pt<40)   return 0.83;
      if (pt>=40 && pt<50)   return 0.84;
      if (pt>=50 && pt<60)   return 0.84;
      if (pt>=60 && pt<80)   return 0.84;
      if (pt>=80 && pt<100)  return 0.84;
      if (pt>=100 && pt<150) return 0.84;
      if (pt>=150 && pt<200) return 0.84;
      if (pt>=200)           return 0.82;
    } 
    else if ( fabs(eta)>=1.5 && fabs(eta)<2.1 ) 
    {
      if (pt>=20 && pt<22)   return 0.00;
      if (pt>=22 && pt<24)   return 0.11;
      if (pt>=24 && pt<26)   return 0.76;
      if (pt>=26 && pt<28)   return 0.78;
      if (pt>=28 && pt<30)   return 0.79;
      if (pt>=30 && pt<32)   return 0.80;
      if (pt>=32 && pt<34)   return 0.80;
      if (pt>=34 && pt<36)   return 0.81;
      if (pt>=36 && pt<38)   return 0.81;
      if (pt>=38 && pt<40)   return 0.82;
      if (pt>=40 && pt<50)   return 0.82;
      if (pt>=50 && pt<60)   return 0.83;
      if (pt>=60 && pt<80)   return 0.83;
      if (pt>=80 && pt<100)  return 0.83;
      if (pt>=100 && pt<150) return 0.83;
      if (pt>=150 && pt<200) return 0.82;
      if (pt>=200)           return 0.82;
    }
  }

  return 1.;

}



TH1F* ratioForMETBeforeTauVeto;
TH1F* ratioForMETAfterTauVeto;

TH1F* ratioForMTBeforeTauVeto;
TH1F* ratioForMTAfterTauVeto;

float yieldDataBeforeTauVeto;
float yieldMCBeforeTauVeto;

float yieldDataAfterTauVeto;
float yieldMCAfterTauVeto;





  // --------- Ratio of ratio

  /*
    TF1* unity = new TF1("unity","1",-1000,1000);
    unity->SetLineColor(kBlack);
    unity->SetLineStyle(1);
    unity->SetLineWidth(1);

 

  TCanvas* canvas_ratioOfRatioMET = new TCanvas("ratioOfRatioMET","",800,200);
  ratioForMETAfterTauVeto->Divide(ratioForMETBeforeTauVeto);
  ratioForMETAfterTauVeto->Draw();
  ratioForMETAfterTauVeto->GetYaxis()->SetTitle("with/without #tau veto");
  ratioForMETAfterTauVeto->GetXaxis()->SetLabelSize(1.0);
  ratioForMETAfterTauVeto->GetXaxis()->SetLabelSize(0.015);
  ratioForMETAfterTauVeto->GetXaxis()->SetTitle("MET");

  unity->Draw("SAME");
  TFile filerrMET("plots/rrMET.root","RECREATE");
  canvas_ratioOfRatioMET->Write();


  TCanvas* canvas_ratioOfRatioMT = new TCanvas("ratioOfRatioMT","",800,200);
  ratioForMTAfterTauVeto->Divide(ratioForMTBeforeTauVeto);
  ratioForMTAfterTauVeto->Draw();
  ratioForMTAfterTauVeto->GetYaxis()->SetTitle("with/without #tau veto");
  ratioForMTAfterTauVeto->GetXaxis()->SetLabelSize(1.0);
  ratioForMTAfterTauVeto->GetXaxis()->SetLabelSize(0.015);
  ratioForMTAfterTauVeto->GetXaxis()->SetTitle("MET");
  unity->Draw("SAME");
  TFile filerrMT("plots/rrMT.root","RECREATE");
  canvas_ratioOfRatioMT->Write();

  cout << "before ; data =" << yieldDataBeforeTauVeto << " ; MC = " << yieldMCBeforeTauVeto << endl;
  cout << " ratio  = " << yieldDataBeforeTauVeto /  yieldMCBeforeTauVeto << endl;
  cout << " after ; data =" << yieldDataAfterTauVeto <<  " ; MC = " << yieldMCAfterTauVeto  << endl;
  cout << " ratio  = " << yieldDataAfterTauVeto /  yieldMCAfterTauVeto << endl;

  cout << " scale factor (after/before) =" << (yieldDataAfterTauVeto /  yieldMCAfterTauVeto) / (yieldDataBeforeTauVeto /  yieldMCBeforeTauVeto) << endl;
    */


