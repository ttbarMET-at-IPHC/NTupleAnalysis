
#include <cmath>
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

bool WTagged(float mass, float pT, float dRlep)
{
    if ((mass > 60) && (mass < 130) && (dRlep > 0.6))
         return true;
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

     float genW_pT;                   mySonic.AddVariable("PtGen",                   "p_{T}(gen W)",                          "GeV",    24,0,600,  &(genW_pT),                  "");
     float genW_hasMatchedAK8;        mySonic.AddVariable("genW_hasMatchedAK8",        "gen W has a matched AK8 jet",           "",     2,-0.5,1.5,  &(genW_hasMatchedAK8),       "");
     float genW_hasMatchedTaggedAK8;  mySonic.AddVariable("genW_hasMatchedTaggedAK8",  "gen W has a matched tagged AK8 jet",    "",     2,-0.5,1.5,  &(genW_hasMatchedTaggedAK8), "");
     float AK8jet_mass;   mySonic.AddVariable("mass",       "mass(AK8 jet)",                         "GeV",    20,0,200,  &(AK8jet_mass),    "");
     float AK8jet_pT;     mySonic.AddVariable("Pt",         "p_{T}(AK8 jet)",                        "GeV",    24,0,600,  &(AK8jet_pT),      "");
     float AK8jet_dR;     mySonic.AddVariable("dR",         "#DeltaR(AK8 jet,gen W)",                   "",      40,0,2,  &(AK8jet_dR),      "");
     float AK8jet_dRlep;  mySonic.AddVariable("dRlep",      "#DeltaR(AK8 jet,sel. l)",                  "",      40,0,2,  &(AK8jet_dRlep),   "");
     float AK8jet_tagged; mySonic.AddVariable("tagged",     "W-tag flag",                               "",  2,-0.5,1.5,  &(AK8jet_tagged),  "");
     float deltaPtRel;    mySonic.AddVariable("deltaPtRel", "p_{T}(AK8 jet) - p_{T}(gen W) / p_{T}(gen W)",  "", 20,-0.3,0.3,  &(deltaPtRel),     "");
  
  // #########################################################
  // ##   Create ProcessClasses (and associated datasets)   ##
  // #########################################################

      mySonic.AddProcessClass("matched","matched #DeltaR<0.1","background",COLORPLOT_ORANGE);
      mySonic.AddProcessClass("fakes",  "fakes",  "background",COLORPLOT_AZURE);
      mySonic.AddProcessClass("gen",    "gen",    "data",COLORPLOT_GREEN);

  // ##########################
  // ##    Create Regions    ##
  // ##########################

     mySonic.AddRegion("presel","Preselection",&Selection_presel);  

  // ##########################
  // ##   Create Channels    ##
  // ##########################
   
	 mySonic.AddChannel("inclusiveChannel","",&inclusiveChannelSelector);

  // ########################################
  // ##       Create histograms and        ##
  // ##  schedule type of plots to produce ##
  // ########################################
  
	 // Create histograms
  	 mySonic.Create1DHistos();

  	 mySonic.Add2DHisto("Pt","tagged");
  	 mySonic.Add2DHisto("Pt","mass");
  	 mySonic.Add2DHisto("dR","mass");
  	 mySonic.Add2DHisto("dRlep","mass");
     mySonic.Add2DHisto("PtGen","genW_hasMatchedAK8");
     mySonic.Add2DHisto("PtGen","genW_hasMatchedTaggedAK8");

  	 // Schedule plots
  	 mySonic.SchedulePlots("1DSuperpRenorm");
  	 mySonic.SchedulePlots("2D");
    
     // (Mis)tagging efficiency plots
     mySonic.SchedulePlots("2DProjectedTo1D","varX=Pt,varY=tagged,projectionType=mean,labelY=Efficiency");

     mySonic.SchedulePlots("2DProjectedTo1D","varX=PtGen,varY=genW_hasMatchedAK8,projectionType=mean,labelY=hasMatched");
     mySonic.SchedulePlots("2DProjectedTo1D","varX=PtGen,varY=genW_hasMatchedTaggedAK8,projectionType=mean,labelY=Efficiency");

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

     int numberOfAK8JetHighPtMatched = 0;
     int numberOfAK8JetHighPtMatchedTagged = 0;

  for (int i = 0 ; i < theTree->GetEntries() ; i++)
  {
	  // Get the i-th entry
      theTree->GetEntry(i);

      genW_pT = -1;
      genW_hasMatchedAK8 = -1;
      genW_hasMatchedTaggedAK8 = -1;

      bool foundMatchedAK8 = false;
      bool foundMatchedTaggedAK8 = false;
      for (int j = 0 ; j < 10 ; j++)
      {
        if (myEvent.recoW_pT[j] == 0.0) break;

         AK8jet_mass = myEvent.recoW_Mass[j];
         AK8jet_pT   = myEvent.recoW_pT[j] ;
         AK8jet_dR   = deltaR(myEvent.genW_Eta,
                              myEvent.genW_Phi,
                              myEvent.recoW_Eta[j],
                              myEvent.recoW_Phi[j]);
         if (AK8jet_dR > 1.99) AK8jet_dR = 1.99;

         AK8jet_dRlep = deltaR(myEvent.lep_eta,
                               myEvent.lep_phi,
                               myEvent.recoW_Eta[j],
                               myEvent.recoW_Phi[j]);
         if (AK8jet_dRlep > 1.99) AK8jet_dRlep = 1.99;

         AK8jet_tagged = WTagged(AK8jet_mass,AK8jet_pT,AK8jet_dRlep);
         deltaPtRel = (AK8jet_pT - myEvent.genW_pT) / myEvent.genW_pT;

         if (AK8jet_dR < 0.1) 
         {
             mySonic.AutoFillProcessClass("matched");
             foundMatchedAK8 = true;
             if (AK8jet_tagged >= 1.0) foundMatchedTaggedAK8 = true;
         
             if (AK8jet_pT >= 200)
             {
                numberOfAK8JetHighPtMatched++;
                if (AK8jet_tagged >= 1.0) numberOfAK8JetHighPtMatchedTagged++;
             }
         }
         else mySonic.AutoFillProcessClass("fakes");

      }

      genW_pT = myEvent.genW_pT;
      if (foundMatchedAK8) genW_hasMatchedAK8 = 1;
      else                 genW_hasMatchedAK8 = 0;

      if (foundMatchedTaggedAK8) genW_hasMatchedTaggedAK8 = 1;
      else                       genW_hasMatchedTaggedAK8 = 0;
      
      mySonic.AutoFillProcessClass("gen");

  } 

  cout << "numberOfAK8JetHighPtMatched = " << numberOfAK8JetHighPtMatched << endl;
  cout << "numberOfAK8JetHighPtMatchedTagged = " <<numberOfAK8JetHighPtMatchedTagged << endl;
  cout << "ratio = " << (float) numberOfAK8JetHighPtMatched / (float) numberOfAK8JetHighPtMatchedTagged << endl;


  // ###################################
  // ##   Make plots and write them   ##
  // ###################################
 
  cout << "   > Making plots..." << endl;
  mySonic.MakePlots();
  cout << "   > Saving plots..." << endl;
  mySonic.WritePlots("plots/WTagPerformances/","","exportPngAndEps=true");

  cout << endl;
  cout << "   ┌──────────────────────────────┐  " << endl;
  cout << "   │   Plot generation completed  │  " << endl;
  cout << "   └──────────────────────────────┘  " << endl; 
  cout << endl;

  return (0);

}


