
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

#include "interface/SonicScrewdriver.h" 
using namespace theDoctor;

typedef struct 
{

    // Dataset
    Float_t dataset_trueNumberOfEvents;
    
    // Jets, bTags info
    Float_t nJets;
    Float_t nBtags;

    Float_t jet_pT[10];
    Float_t jet_eta[10];
    Float_t jet_phi[10];
    Float_t jet_bTag[10];

    // Lepton info
    Float_t lep_pT;
    Float_t lep_eta;
    Float_t lep_phi;
    Float_t lep_flavor;

    // Variables of intereset 
    Float_t MET;
    Float_t MT;
    
    // Vetos
    Float_t isolatedTrackVeto;
    Float_t tauVeto;

    // Generated hadronic W
    Float_t genW_pT;
    Float_t genW_eta;
    Float_t genW_phi;

    // W-tagged jet candidates
    Float_t recoW_pT[10];
    Float_t recoW_eta[10];
    Float_t recoW_phi[10];
    Float_t recoW_mass[10];

} 
EventStruct;




// #########################################################################""
// #########################################################################""

EventStruct* myEventPointer;

bool baselineSelector();
bool noChannelSelector();

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
	  EventStruct myEvent;
	  myEventPointer = &myEvent;

  // ##########################
  // ##   Create Variables   ##
  // ##########################

 	 mySonic.AddVariable("MET", "MET",   "GeV", 25,0,500,  &(myEvent.MET)   );
 	 mySonic.AddVariable("MT",  "M_{T}", "GeV", 25,0,500,   &(myEvent.MT)    );
 	 
     float recoW_bestMass;
     mySonic.AddVariable("mass(Wcand)",  "M(W_{cand})", "GeV", 20,0,200,   &(recoW_bestMass)    );
     float deltaR_bjet;
     mySonic.AddVariable("deltaR(b-jet)",  "#DeltaR(W_{cand},b-jet)", "GeV", 20,0,5,   &(deltaR_bjet)    );
 	 
  // #########################################################
  // ##   Create ProcessClasses (and associated datasets)   ##
  // #########################################################

 	 mySonic.AddProcessClass("dilepton", "t#bar{t}#rightarrow ll",  "background",COLORPLOT_AZURE);
  
  // ##########################
  // ##    Create Regions    ##
  // ##########################

     mySonic.AddRegion("baseline","Baseline",&baselineSelector);
   	 	mySonic.ScheduleVariablesForRegion("baseline","all");
     
     mySonic.AddRegion("selection","Selection",&baselineSelector);
   	 	mySonic.ScheduleVariablesForRegion("selection","all");
       
  // ##########################
  // ##   Create Channels    ##
  // ##########################
   
	 mySonic.AddChannel("noChannel","",&noChannelSelector);

  // ########################################
  // ##       Create histograms and        ##
  // ##  schedule type of plots to produce ##
  // ########################################
  
	 // Create histograms
  	 mySonic.Create1DHistos();

  	 // Schedule plots
  	 mySonic.SchedulePlots("1DSuperpRenorm");

  // ########################################
  // ##        Run over the events         ##
  // ########################################

  	 // Open the tree
  	 TFile* f = new TFile("microTuples_WTag_423/test.root");
  	 TTree* theTree;
  	 f->GetObject("microTuple",theTree);
  	 theTree->SetBranchAddress("events",&myEvent);
 
  	 // Loop over the events

  int eventSelected = 0;
  int eventTotal = 0;
  for (int i = 0 ; i < theTree->GetEntries() ; i++)
  {
	  // Get the i-th entry
      theTree->GetEntry(i);
    
      //cout << " nJets, nBtags = " << myEvent.nJets << ", " << myEvent.nBtags << endl;
 
      if (myEvent.nBtags <= 0) continue;

      eventTotal++;
      deltaR_bjet = -1;

      recoW_bestMass = 0;
      for (int j = 0 ; j < 10 ; j++)
      {
            if ((myEvent.recoW_mass[j] < 60) || (myEvent.recoW_mass[j] > 130)) continue;
           
            recoW_bestMass = myEvent.recoW_mass[j];
      
            for (int k = 0 ; k < 10 ; k++)
            {
                if (myEvent.jet_bTag[k] != 1.0) continue;

                float dPhi = myEvent.jet_phi[k] - myEvent.recoW_phi[j];
                float dEta = myEvent.jet_eta[k] - myEvent.recoW_eta[j];
                float deltaR_test = sqrt(dPhi*dPhi + dEta*dEta);

                if (deltaR_test < deltaR_bjet) deltaR_bjet = deltaR_test;
            }
            eventSelected++;
	        mySonic.AutoFillProcessClass("dilepton");
            break;
      }


	  // Fill all the variables with autoFill-mode activated
	  mySonic.AutoFillProcessClass("dilepton");

  }

    cout << " eventTotal    : " << eventTotal    << endl; 
    cout << " eventSelected : " << eventSelected << endl; 
    
    cout << " eff : " << (100 * eventSelected) / eventTotal << endl; 


  // ###################################
  // ##   Make plots and write them   ##
  // ###################################
 
  
  mySonic.MakePlots();
  mySonic.WritePlots("plots/WTagging/");
  

  cout << endl;
  cout << "   ┌──────────────────────────────┐  " << endl;
  cout << "   │   Plot generation completed  │  " << endl;
  cout << "   └──────────────────────────────┘  " << endl; 
  cout << endl;
  
  return (0);

}

bool baselineSelector()
{
	return true;
}


bool noChannelSelector()
{
	return true;
}

