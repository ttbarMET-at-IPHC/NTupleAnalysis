
#include <iomanip>
#include <iostream>
#include <time.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <time.h>

#include <TFile.h>
#include <TMarker.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLorentzVector.h>

using namespace std;

#include "MicroTuple_Format_WTag0430.h" 
microEvent* myEventPointer;

#include "interface/Table.h" 
#include "interface/SonicScrewdriver.h" 
using namespace theDoctor;

#define FOLDER_MICROTUPLES "microTuples_WTag0430-2/" 


// #########################################################################
//                          Region selectors
// #########################################################################

bool WTagged(float mass, float pT)
{
    if ((pT > 0) && (mass > 60) && (mass < 130)) return true;
    else return false;
}

bool inclusiveChannelSelector() { return true; }

bool Selection_presel() { return true; }

bool Selection_MTcutOnly() 
{ 
    if (myEventPointer->MT                >= 120)
         return true;
    else return false;
}

bool Selection_baseline()
{
    if ( 
          (myEventPointer->nJets             >= 4)
    &&    (myEventPointer->isolatedTrackVeto == 1)
    &&    (myEventPointer->tauVeto           == 1)
    &&    (myEventPointer->MET               >= 100)
    &&    (myEventPointer->MT                >= 120)
       )
         return true;
    else return false;
}

bool Selection_lowdM()
{
    if ( 
          (myEventPointer->DPhi_MetJets > 0.8)
       && (myEventPointer->HadronicChi2 < 5.0)
       )
         return Selection_baseline();
    else return false;

}

bool Selection_highdM()
{
    if (myEventPointer->MT2W > 200) return Selection_lowdM();
    else                            return false;
}

bool Selection_WtagSelection()
{
    bool foundWtaggedJet = false;
    for (int i = 0 ; i < 10 ; i++)
    {
        if (WTagged(myEventPointer->recoW_Mass[i],myEventPointer->recoW_pT[i]))
        { foundWtaggedJet = true; break; }
    }

    if ( 
          (myEventPointer->nJets             >= 3)
    &&    (foundWtaggedJet)
    &&    (myEventPointer->isolatedTrackVeto == 1)
    &&    (myEventPointer->tauVeto           == 1)
    &&    (myEventPointer->MET               >= 200)
    &&    (myEventPointer->MT                >= 120)
    &&    (myEventPointer->DPhi_MetJets      > 0.8)
       )
         return true;
    else return false;
}

bool Selection_WtagSelection_noW()
{
    if ( 
          (myEventPointer->nJets             >= 3)
    &&    (myEventPointer->isolatedTrackVeto == 1)
    &&    (myEventPointer->tauVeto           == 1)
    &&    (myEventPointer->MET               >= 200)
    &&    (myEventPointer->MT                >= 120)
    &&    (myEventPointer->DPhi_MetJets      > 0.8)
       )
         return true;
    else return false;
}


// #########################################################################
//                          Others tools/stuff
// #########################################################################

float stopCrossSection(float);

// #########################################################################
//                              Main function
// #########################################################################

int main (int argc, char *argv[])
{

  cout << endl;
  cout << "   ┌──────────────────────────────┐  " << endl;
  cout << "   │   Starting plot generation   │  " << endl;
  cout << "   └──────────────────────────────┘  " << endl; 
  cout << endl;

  // ####################
  // ##   Init tools   ##
  // ####################
  
	 // Create a sonic Screwdriver
 	 SonicScrewdriver mySonic;

  	 // Create a container for the event
	 microEvent myEvent;
	 myEventPointer = &myEvent;

  // ##########################
  // ##   Create Variables   ##
  // ##########################

 	 mySonic.AddVariable("mStop",          "m_{#tilde{t}}",           "GeV",    29,87.5,812.5,    &(myEvent.mStop),        "");
 	 mySonic.AddVariable("mNeutralino",    "m_{#chi^0}",              "GeV",    17,-12.5,412.5,   &(myEvent.mNeutralino),  "");
 	 
     // AK8 collection 
     float nJetsAK8;            mySonic.AddVariable("nJetsAK8",       "# of AK8 jets (p_{T}>20)",      "", 10,-0.5,9.5,   &(nJetsAK8),         "");
     
     // W-tagging
     mySonic.AddVariable("genW_Pt",        "p_{T}(gen. W)",           "GeV",    30,0,900,         &(myEvent.genW_pT),      "");
     float genW_recoFlag;       mySonic.AddVariable("genW_recoFlag",  "gen W has been reconstructed",  "",  2,-0.5,1.5,   &(genW_recoFlag),    "");
                                mySonic.AddVariable("genW_recoFlag2", "gen W has been reconstructed",  "",  2,-0.5,1.5,   &(genW_recoFlag),    "");
     float nJetsAK8Tagged;      mySonic.AddVariable("nJetsAK8Tagged", "# of AK8 W-tagged jets",        "",  4,-0.5,3.5,   &(nJetsAK8Tagged),   "");
     
     // W-mistagging
   //  mySonic.AddVariable("fakeW_Pt",       "p_{T}(fake reco. W)",           "GeV",    30,0,600);
   //  mySonic.AddVariable("fakeW_tagFlag",  "Fake W cand. has been tagged",  "GeV",    30,0,600);
 	 
  // #########################################################
  // ##   Create ProcessClasses (and associated datasets)   ##
  // #########################################################

      // ttbar -> l + jets
      //mySonic.AddProcessClass("ttbarSemiLept",    "t#bar{t} #rightarrow l^{#pm} + jets",     "background",COLORPLOT_RED);
      mySonic.AddProcessClass("1lepTop",          "1l top",     "background",kRed-7);

          mySonic.AddDataset("ttbarSemiLept",     "1lepTop",100,234*(0.67*0.33*2) );

     // Signal
     //mySonic.AddProcessClass("signal", "signal",  "signal",COLORPLOT_AZURE);
      //    mySonic.AddDataset("signal", "signal", 100,1);

  // ##########################
  // ##    Create Regions    ##
  // ##########################

     mySonic.AddRegion("presel",         "Preselection",             &Selection_presel);     mySonic.ScheduleVariablesForRegion("presel","all");
     mySonic.AddRegion("MTCutOnly",      "MT cut only",              &Selection_MTcutOnly);  mySonic.ScheduleVariablesForRegion("MTCutOnly","all");
     mySonic.AddRegion("baseline",       "Baseline",                 &Selection_baseline);   mySonic.ScheduleVariablesForRegion("baseline","all");
 
//     mySonic.AddRegion("W-tagSelection_noW", "W-tag selection (noW)",&Selection_WtagSelection_noW); mySonic.ScheduleVariablesForRegion("W-tagSelection_noW","all");
//     mySonic.AddRegion("W-tagSelection",     "W-tag selection",      &Selection_WtagSelection);     mySonic.ScheduleVariablesForRegion("W-tagSelection","all");


  // ##########################
  // ##   Create Channels    ##
  // ##########################
   
	 mySonic.AddChannel("inclusiveChannel","",&inclusiveChannelSelector);

  // ########################################
  // ##       Create histograms and        ##
  // ##  schedule type of plots to produce ##
  // ########################################
  
     mySonic.SetLumi(20000);

	 // Create histograms
  	 mySonic.Create1DHistos();

  	 mySonic.Add2DHisto("mStop","mNeutralino");
  	 mySonic.Add2DHisto("genW_Pt","genW_recoFlag");
  	 mySonic.Add2DHisto("genW_Pt","genW_recoFlag2");
  	// mySonic.Add2DHisto("fakeW_Pt","fakeW_tagFlag");

  	 mySonic.Add3DHisto("mStop","mNeutralino","genW_Pt");

  	 // Schedule plots
  	 mySonic.SchedulePlots("1DSuperpRenorm");
  	 mySonic.SchedulePlots("1DStack");
  	 mySonic.SchedulePlots("2D");
     
     mySonic.SchedulePlots("2DProjectedTo1D","varX=genW_Pt,varY=genW_recoFlag,projectionType=mean,labelY=TaggingEfficiency");
 //  mySonic.SchedulePlots("2DProjectedTo1D","varX=fakeW_Pt,varY=fakeW_tagFlag,projectionType=mean,labelY=MistaggingEfficiency");

     mySonic.SchedulePlots("3DProjectedTo2D","varX=mStop,varY=mNeutralino,varZ=genW_Pt,projectionType=mean,labelZ=mean(genW_Pt)");
  
  // ########################################
  // ##       Run over the datasets        ##
  // ########################################

  TH2F* signalInitialNumberOfEvents = new TH2F("signalInitialNumberOfEvents_","",29,87.5,812.5,17,-12.5,412.5);
  vector<string> datasetsList;
  mySonic.GetDatasetList(&datasetsList);

  cout << "   > Running on dataset : " << endl;

  for (unsigned int d = 0 ; d < datasetsList.size() ; d++)
  {
     string currentDataset = datasetsList[d];
     string currentProcessClass = mySonic.GetProcessClass(currentDataset);

     // Open the tree
     TFile f((string(FOLDER_MICROTUPLES)+currentDataset+".root").c_str());
     TTree* theTree = (TTree*) f.Get("microTuple"); 
     theTree->SetBranchAddress("microEvents",&myEvent);
    
     // Read the number of event for signal (for each datapoint)
     if (currentDataset == "signal")
     {
         TTree* theTree2 = (TTree*) f.Get("signalInitialNumberOfEvents"); 
         signalPoint mySignalPoint;
         theTree2->SetBranchAddress("signalPoint",&mySignalPoint);
         for (int i = 0 ; i < theTree2->GetEntries() ; i++)
         {
            theTree2->GetEntry(i);
            signalInitialNumberOfEvents->Fill(mySignalPoint.mStop,mySignalPoint.mNeutralino);
         }
     }

     cout << "                    " << currentDataset << " (" << theTree->GetEntries();

  // ########################################
  // ##        Run over the events         ##
  // ########################################

  float weight_lumi = 0.0;
  for (int i = 0 ; i < theTree->GetEntries() ; i++)
  {

	  // Get the i-th entry
      theTree->GetEntry(i);

      if (i == 0) 
      {    
           mySonic.SetTrueNumberOfEvents(currentDataset,myEvent.trueNumberOfEvents);
           cout << " events from " << myEvent.trueNumberOfEvents << " before skimming)" << endl;
           weight_lumi = mySonic.GetDatasetLumiWeight(currentDataset);
      }

      float weight = weight_lumi;
      // Weight computation for signal
      if (currentDataset == "signal")
      {
         weight = mySonic.GetLumi() * stopCrossSection(myEvent.mStop) 
                 / signalInitialNumberOfEvents->GetBinContent(  
                          signalInitialNumberOfEvents->FindBin(myEvent.mStop,myEvent.mNeutralino)  
                                                             );
         weight *= 0.9; // Trigger efficiency (trigger were not included for selection on signal when microtupling)
      }






       weight = 1.0;






      // Dirty management of had. chi2 overflow
      if (myEvent.HadronicChi2 >= 19.99) myEvent.HadronicChi2 = 19.99;


      // Tagging/Mistagging efficiencies
      nJetsAK8       = 0;       
      nJetsAK8Tagged = 0;       
      genW_recoFlag  = 0;
      for (int j = 0 ; j < 10 ; j++)
      {
        if (myEventPointer->recoW_pT[j] == 0) break;
        else nJetsAK8++;

        if (WTagged(myEvent.recoW_Mass[j],myEvent.recoW_pT[j]))
        {
            nJetsAK8Tagged++;
            float dEta = myEvent.genW_Eta - myEvent.recoW_Eta[j]; 
            float dPhi = myEvent.genW_Phi - myEvent.recoW_Phi[j]; 
            float deltaR = sqrt(dEta*dEta + dPhi*dPhi);
            if (deltaR < 0.1) genW_recoFlag = 1;
        }
      }


	  // Fill all the variables with autoFill-mode activated
	  mySonic.AutoFillProcessClass(currentProcessClass,weight);

  } }

  // ###################################
  // ##   Make plots and write them   ##
  // ###################################
 
  std::ostringstream lumiInString;
  lumiInString << mySonic.GetLumi() / 1000;
  string infoText("CMS Internal      #sqrt{s} = 8 TeV, ");
  infoText += string("L = ")+lumiInString.str()+string(" fb^{-1}, ");
  
  cout << endl;
  cout << "   > Making plots..." << endl;
  mySonic.MakePlots();
  cout << "   > Saving plots..." << endl;
  mySonic.WritePlots("plots/WTagPerformances/",infoText); //"exportPngAndEps=true"

  cout << endl;
  cout << "   ┌──────────────────────────────┐  " << endl;
  cout << "   │   Plot generation completed  │  " << endl;
  cout << "   └──────────────────────────────┘  " << endl; 
  cout << endl;

  return (0);
}



float stopCrossSection(float inputMass)
{
         if (abs(inputMass - 170) <= 5) return 42.6441;
    else if (abs(inputMass - 180) <= 5) return 31.8695;
    else if (abs(inputMass - 190) <= 5) return 24.1585;
    else if (abs(inputMass - 200) <= 5) return 18.5245;
    else if (abs(inputMass - 210) <= 5) return 14.3201;
    else if (abs(inputMass - 220) <= 5) return 11.1808;
    else if (abs(inputMass - 230) <= 5) return 8.78125;
    else if (abs(inputMass - 240) <= 5) return 6.96892;
    else if (abs(inputMass - 250) <= 5) return 5.57596;
    else if (abs(inputMass - 260) <= 5) return 4.48773;
    else if (abs(inputMass - 270) <= 5) return 3.63085;
    else if (abs(inputMass - 280) <= 5) return 2.95613;
    else if (abs(inputMass - 290) <= 5) return 2.42299;
    else if (abs(inputMass - 300) <= 5) return 1.99608;
    else if (abs(inputMass - 310) <= 5) return 1.64956;
    else if (abs(inputMass - 320) <= 5) return 1.3733;
    else if (abs(inputMass - 330) <= 5) return 1.14277;
    else if (abs(inputMass - 340) <= 5) return 0.959617;
    else if (abs(inputMass - 350) <= 5) return 0.807323;
    else if (abs(inputMass - 360) <= 5) return 0.681346;
    else if (abs(inputMass - 370) <= 5) return 0.576882;
    else if (abs(inputMass - 380) <= 5) return 0.489973;
    else if (abs(inputMass - 390) <= 5) return 0.4176;
    else if (abs(inputMass - 400) <= 5) return 0.35683;
    else if (abs(inputMass - 410) <= 5) return 0.305512;
    else if (abs(inputMass - 420) <= 5) return 0.262683;
    else if (abs(inputMass - 430) <= 5) return 0.226367;
    else if (abs(inputMass - 440) <= 5) return 0.195812;
    else if (abs(inputMass - 450) <= 5) return 0.169668;
    else if (abs(inputMass - 460) <= 5) return 0.147492;
    else if (abs(inputMass - 470) <= 5) return 0.128326;
    else if (abs(inputMass - 480) <= 5) return 0.112241;
    else if (abs(inputMass - 490) <= 5) return 0.0977878;
    else if (abs(inputMass - 500) <= 5) return 0.0855847;
    else if (abs(inputMass - 510) <= 5) return 0.0751004;
    else if (abs(inputMass - 520) <= 5) return 0.0660189;
    else if (abs(inputMass - 530) <= 5) return 0.0580348;
    else if (abs(inputMass - 540) <= 5) return 0.0511747;
    else if (abs(inputMass - 550) <= 5) return 0.0452067;
    else if (abs(inputMass - 560) <= 5) return 0.0399591;
    else if (abs(inputMass - 570) <= 5) return 0.0354242;
    else if (abs(inputMass - 580) <= 5) return 0.0313654;
    else if (abs(inputMass - 590) <= 5) return 0.0279395;
    else if (abs(inputMass - 600) <= 5) return 0.0248009;
    else if (abs(inputMass - 610) <= 5) return 0.0220672;
    else if (abs(inputMass - 620) <= 5) return 0.0196331;
    else if (abs(inputMass - 630) <= 5) return 0.0175075;
    else if (abs(inputMass - 640) <= 5) return 0.0155809;
    else if (abs(inputMass - 650) <= 5) return 0.0139566;
    else if (abs(inputMass - 660) <= 5) return 0.0125393;
    else if (abs(inputMass - 670) <= 5) return 0.0112223;
    else if (abs(inputMass - 680) <= 5) return 0.0100516;
    else if (abs(inputMass - 690) <= 5) return 0.0090306;
    else if (abs(inputMass - 700) <= 5) return 0.0081141;
    else if (abs(inputMass - 710) <= 5) return 0.00730084;
    else if (abs(inputMass - 720) <= 5) return 0.00656729;
    else if (abs(inputMass - 730) <= 5) return 0.00591771;
    else if (abs(inputMass - 740) <= 5) return 0.00532605;
    else if (abs(inputMass - 750) <= 5) return 0.00480639;
    else if (abs(inputMass - 760) <= 5) return 0.00433688;
    else if (abs(inputMass - 770) <= 5) return 0.00391839;
    else if (abs(inputMass - 780) <= 5) return 0.00354211;
    else if (abs(inputMass - 790) <= 5) return 0.00320476;
    else if (abs(inputMass - 800) <= 5) return 0.00289588;
    else if (abs(inputMass - 810) <= 5) return 0.0026184;
    else if (abs(inputMass - 820) <= 5) return 0.00237168;
    else if (abs(inputMass - 830) <= 5) return 0.00214607;
    else if (abs(inputMass - 840) <= 5) return 0.00195172;
    else if (abs(inputMass - 850) <= 5) return 0.00176742;
    else if (abs(inputMass - 860) <= 5) return 0.00160403;
    else if (abs(inputMass - 870) <= 5) return 0.00145772;
    else if (abs(inputMass - 880) <= 5) return 0.00132077;
    else if (abs(inputMass - 890) <= 5) return 0.00120568;
    else if (abs(inputMass - 900) <= 5) return 0.00109501;
    else if (abs(inputMass - 910) <= 5) return 0.000996193;
    else if (abs(inputMass - 920) <= 5) return 0.000907494;
    else if (abs(inputMass - 930) <= 5) return 0.000826533;
    else if (abs(inputMass - 940) <= 5) return 0.000753768;
    else if (abs(inputMass - 950) <= 5) return 0.000687022;
    else if (abs(inputMass - 960) <= 5) return 0.000626876;
    else if (abs(inputMass - 970) <= 5) return 0.000571551;
    else if (abs(inputMass - 980) <= 5) return 0.000522495;
    else if (abs(inputMass - 990) <= 5) return 0.000476255;
    else if (abs(inputMass - 1000) <= 5) return 0.000435488;
    else return 0.0;
}
