
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

float* fakeW_dRlep_ptr;

bool WTagged(float mass, float pT)
{
    if ((pT > 200) && (mass > 60) && (mass < 130) && (*fakeW_dRlep_ptr > 0.7))
         return true;
    //if ((mass > 60) && (mass < 130)) return true;
    else return false;
}

float deltaR(float genEta, float genPhi, float recoEta, float recoPhi)
{
    float dEta = genEta - recoEta;
    float dPhi = TVector2::Phi_mpi_pi(genPhi - recoPhi);
    float dR = sqrt(dEta*dEta + dPhi*dPhi);
    return dR;
}

// #########################################################################
//                          Region selectors
// #########################################################################

bool inclusiveChannelSelector() { return true; }

bool Selection_presel() { return true; }

bool Selection_dRlep_07() 
{
    if (*fakeW_dRlep_ptr > 0.7)
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
 	 SonicScrewdriver mySonicTagEff;
 	 SonicScrewdriver mySonicMistagEff;

  	 // Create a container for the event
	 microEvent myEvent;
	 myEventPointer = &myEvent;

  // ##########################
  // ##   Create Variables   ##
  // ##########################

     // W-tagging
                                mySonicTagEff.AddVariable("genW_Pt",        "p_{T}(gen. W)",              "GeV",  30,0,900,     &(myEvent.genW_pT),        "");
     float recoW_tagged;        mySonicTagEff.AddVariable("genW_recoFlag",  "gen W has been reconstructed",  "",  2,-0.5,1.5,   &(genW_recoFlag),  "");
     float recoW_pT;            mySonicTagEff.AddVariable("recoW_Pt",       "p_{T}(reco W)",              "GeV",  15,0,900,     &(recoW_pT),      "");
     float recoW_dR;            mySonicTagEff.AddVariable("recoW_dR",       "#DeltaR(gen W,reco W)",      "",     20,0,2,       &(recoW_dR),      "");
     float deltaPtRel_genReco;  mySonicTagEff.AddVariable("deltaPtRel",     "#Deltap_{T}(gen W,reco W)/p_{T}(gen W)",  "",     20,-0.1,0.1,       &(deltaPtRel_genReco),      "");

     // W-mistagging
     float fakeW_pT;      mySonicMistagEff.AddVariable("fakeW_Pt",       "p_{T}(fake. W)",                 "GeV",  15,0,900,   &(fakeW_pT),      "");
     float fakeW_Mass;    mySonicMistagEff.AddVariable("fakeW_Mass",     "mass(fake. W)",                  "GeV",  20,0,200,   &(fakeW_Mass),    "");
     float fakeW_tagFlag; mySonicMistagEff.AddVariable("fakeW_tagFlag",  "Fake W cand. has been tagged",   "",     2,-0.5,1.5, &(fakeW_tagFlag), "");
     float fakeW_dR;      mySonicMistagEff.AddVariable("fakeW_dR",       "#DeltaR(gen W,fake W)",          "",     20,0,2,     &(fakeW_dR),      "");
     float fakeW_dRlep;   mySonicMistagEff.AddVariable("fakeW_dRlep",    "#DeltaR(lepton,fake W)",         "",     20,0,2,     &(fakeW_dRlep),      "");
     fakeW_dRlep_ptr = &fakeW_dRlep;

  // #########################################################
  // ##   Create ProcessClasses (and associated datasets)   ##
  // #########################################################

      mySonicTagEff.AddProcessClass("process","process","background",COLORPLOT_AZURE);
      mySonicMistagEff.AddProcessClass("process","process","background",COLORPLOT_AZURE);

  // ##########################
  // ##    Create Regions    ##
  // ##########################

     mySonicTagEff.AddRegion("presel",         "Preselection",             &Selection_presel);  mySonicTagEff.ScheduleVariablesForRegion("presel","all");
     mySonicMistagEff.AddRegion("presel",      "Preselection",             &Selection_presel);  mySonicMistagEff.ScheduleVariablesForRegion("presel","all");
     mySonicMistagEff.AddRegion("dRlep>0.7",   "#DeltaR(lep,jet)>0.7",     &Selection_dRlep_07);  mySonicMistagEff.ScheduleVariablesForRegion("dRlep>0.7","all");

  // ##########################
  // ##   Create Channels    ##
  // ##########################
   
	 mySonicTagEff.AddChannel("inclusiveChannel","",&inclusiveChannelSelector);
	 mySonicMistagEff.AddChannel("inclusiveChannel","",&inclusiveChannelSelector);

  // ########################################
  // ##       Create histograms and        ##
  // ##  schedule type of plots to produce ##
  // ########################################
  
	 // Create histograms
  	 mySonicTagEff.Create1DHistos();
  	 mySonicMistagEff.Create1DHistos();

  	 mySonicTagEff.Add2DHisto("genW_Pt","recoW_dR");
  	 mySonicTagEff.Add2DHisto("genW_Pt","genW_recoFlag");
  	 mySonicTagEff.Add2DHisto("recoW_Pt","genW_recoFlag");
  	 mySonicMistagEff.Add2DHisto("fakeW_Pt","fakeW_tagFlag");
  	 mySonicMistagEff.Add2DHisto("fakeW_dR","fakeW_Mass");
  	 mySonicMistagEff.Add2DHisto("fakeW_dRlep","fakeW_Mass");
  	 mySonicMistagEff.Add2DHisto("fakeW_Pt","fakeW_Mass");

  	 // Schedule plots
  	 mySonicTagEff.SchedulePlots("1DSuperpRenorm");
  	 mySonicMistagEff.SchedulePlots("1DSuperpRenorm");
  	 mySonicTagEff.SchedulePlots("2D");
  	 mySonicMistagEff.SchedulePlots("2D");
     
     mySonicTagEff.SchedulePlots("2DProjectedTo1D","varX=genW_Pt,varY=genW_recoFlag,projectionType=mean,labelY=TaggingEfficiency");
     mySonicTagEff.SchedulePlots("2DProjectedTo1D","varX=recoW_Pt,varY=genW_recoFlag,projectionType=mean,labelY=TaggingEfficiency");
     mySonicMistagEff.SchedulePlots("2DProjectedTo1D","varX=fakeW_Pt,varY=fakeW_tagFlag,projectionType=mean,labelY=MistaggingEfficiency");

  // ########################################
  // ##       Run over the datasets        ##
  // ########################################

     // Open the tree
     TFile f((string(FOLDER_MICROTUPLES)+"ttbarSemiLept.root").c_str());
     TTree* theTree = (TTree*) f.Get("microTuple"); 
     theTree->SetBranchAddress("microEvents",&myEvent);
    
  // ########################################
  // ##        Run over the events         ##
  // ########################################

  for (int i = 0 ; i < theTree->GetEntries() ; i++)
  {
	  // Get the i-th entry
      theTree->GetEntry(i);

      genW_recoFlag = 0;
      recoW_pT      = -1;
      recoW_dR      = -1;
      deltaPtRel_genReco = -1;

      for (int j = 0 ; j < 10 ; j++)
      {
        if (myEvent.recoW_pT[j] == 0.0) break;

        bool tagged  = WTagged(myEvent.recoW_Mass[j],myEvent.recoW_pT[j]);
        float dR = deltaR(myEvent.genW_Eta,
                          myEvent.genW_Phi,
                          myEvent.recoW_Eta[j],
                          myEvent.recoW_Phi[j]);
 
        float dRlep = deltaR(myEvent.lep_eta,
                          myEvent.lep_phi,
                          myEvent.recoW_Eta[j],
                          myEvent.recoW_Phi[j]);
            
        if (dR > 1.99) dR = 1.99;

        if (dR < 0.1)
        {
            if (tagged) 
            {   
                recoW_tagFlag = 1;
            }
            else 
            {
                recoW_tagFlag = 0;
            }
            recoW_pT = myEvent.recoW_pT[j];
            recoW_dR = dR; 
            deltaPtRel_genReco = (myEvent.genW_pT - myEvent.recoW_pT[j])/ myEvent.genW_pT;
	        mySonicTagEff.AutoFillProcessClass("process");
        }
        else
        {
            fakeW_pT      = myEvent.recoW_pT[j];
            fakeW_Mass    = myEvent.recoW_Mass[j];
            fakeW_dR      = dR; 
            fakeW_dRlep   = dRlep; 
            
            if (tagged) fakeW_tagFlag = 1.0;
            else        fakeW_tagFlag = 0.0;
	        mySonicMistagEff.AutoFillProcessClass("process");

        }
      }
  } 

  // ###################################
  // ##   Make plots and write them   ##
  // ###################################
 
  cout << endl;
  cout << "   > Making plots..." << endl;
  mySonicTagEff.MakePlots();
  mySonicMistagEff.MakePlots();
  cout << "   > Saving plots..." << endl;
  mySonicTagEff.WritePlots("plots/WTagPerfTag/",""); //"exportPngAndEps=true"
  mySonicMistagEff.WritePlots("plots/WTagPerfMisTag/",""); //"exportPngAndEps=true"

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
