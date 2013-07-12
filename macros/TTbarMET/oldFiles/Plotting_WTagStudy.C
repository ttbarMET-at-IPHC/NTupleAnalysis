
#include <iomanip>
#include <iostream>
#include <time.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <time.h>

#include <TFile.h>
#include <TRandom.h>
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


bool WTagged(float mass, float pT, float dR_lepton)
{
    if ((pT > 200) && (mass > 60) && (mass < 130) && (dR_lepton > 0.6)) return true;
    else return false;
}

TRandom* r0 = new TRandom();

bool WTagged_Randomize_tag(float mass, float pT, float dR_lepton)
{    if ((pT > 200) && (r0->Rndm() < 0.1)) return true;    else return false;   }
bool WTagged_Randomize_tag(float mass, float pT, float dR_lepton)
{    if ((pT > 200) && (r0->Rndm() < 0.1*1.1)) return true;    else return false;   }
bool WTagged_Randomize_tag(float mass, float pT, float dR_lepton)
{    if ((pT > 200) && (r0->Rndm() < 0.1*0.9)) return true;    else return false;   }

bool WTagged_Randomize_mistag(float mass, float pT, float dR_lepton)
{    if ((pT > 200) && (r0->Rndm() < 0.1)) return true;    else return false;   }
bool WTagged_Randomize_mistag_plus1Sigma(float mass, float pT, float dR_lepton)
{    if ((pT > 200) && (r0->Rndm() < 0.1*1.1)) return true;    else return false;   }
bool WTagged_Randomize_mistag_minus1Sigma(float mass, float pT, float dR_lepton)
{    if ((pT > 200) && (r0->Rndm() < 0.1*0.9)) return true;    else return false;   }


bool WTagged_Randomize_mistag(float mass, float pT, float dR_lepton)



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

bool Selector_presel() { 
    return true;
}


/*
bool Selector_baseline()
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

bool Selector_lowdM(float cutMET)
{
    if ( 
          (myEventPointer->DPhi_MetJets > 0.8   )
       && (myEventPointer->HadronicChi2 < 5.0   )
       && (myEventPointer->MET          > cutMET)
       )
         return Selector_baseline();
    else return false;

}

bool Selector_highdM(float cutMET)
{
    if (myEventPointer->MT2W > 200) return Selector_lowdM(cutMET);
    else                            return false;
}

bool Selector_lowdM_200()  { return Selector_lowdM(200); }    
bool Selector_highdM_200() { return Selector_highdM(200); }
*/


bool Selector_highdMSignal()
{
  
    if ((myEventPointer->mStop == -1)
     &&  (myEventPointer->mNeutralino == -1))
        return true;
    else if (
        (myEventPointer->mStop >= 700)
     && (myEventPointer->mStop <= 725)
     && (myEventPointer->mNeutralino >= 0)
     && (myEventPointer->mNeutralino <= 25))
      return true;
    else
      return false; 

}

bool Selector_3jets()
{
    if ( 
          (myEventPointer->nJets             >= 3)
    &&    (myEventPointer->isolatedTrackVeto == 1)
    &&    (myEventPointer->tauVeto           == 1)
    &&    (myEventPointer->MET               >= 200)
    &&    (myEventPointer->MT                >= 120)
    &&    (myEventPointer->DPhi_MetJets      > 0.8)
    &&    (myEventPointer->MT2W              > 200)
       )
         return Selector_highdMSignal();
    else return false;
}

bool Selector_3jetsExactly()         { if (myEventPointer->nJets == 3) { return Selector_3jets(); } else return false; }
bool Selector_4jets()                { if (myEventPointer->nJets >= 4) { return Selector_3jets(); } else return false; }

bool Selector_AtLeastOneWTag()
{
    for (int i = 0 ; i < 10 ; i++)
        if (WTagged(myEventPointer->recoW_Mass[i],       myEventPointer->recoW_pT[i],
                    deltaR(myEventPointer->lep_eta,      myEventPointer->lep_phi,
                           myEventPointer->recoW_Eta[i], myEventPointer->recoW_Phi[i])
                    ))
        return true;

    return false;
}

bool Selector_3jetsPlusWtag()        { if (Selector_AtLeastOneWTag())    { return Selector_3jets();        } else return false; }
bool Selector_3jetsExactlyPlusWtag() { if (Selector_AtLeastOneWTag())    { return Selector_3jetsExactly(); } else return false; }
bool Selector_4jetsPlusWtag()        { if (Selector_AtLeastOneWTag())    { return Selector_4jets();        } else return false; }

bool Selector_3jetsPlusVetoWtag()        { if (!Selector_AtLeastOneWTag())   { return Selector_3jets();        } else return false; }
bool Selector_3jetsExactlyPlusVetoWtag() { if (!Selector_AtLeastOneWTag())   { return Selector_3jetsExactly(); } else return false; }
bool Selector_4jetsPlusVetoWtag()        { if (!Selector_AtLeastOneWTag())   { return Selector_4jets();        } else return false; }

bool Selector_HighPtGenW()
{
    for (int i = 0 ; i < 10 ; i++)
        if (myEventPointer->genW_pT > 200) return true;

    return false;
}
bool Selector_3jetsPlusHighPtGenW()    { if (!Selector_HighPtGenW())   { return Selector_3jets(); } else return false; }

// #########################################################################
//                          Others tools/stuff
// #########################################################################

float stopCrossSection(float);
void fillMCSignalTable(SonicScrewdriver* mySonic, vector<string> region, vector<string> process, Table* table);

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

 	 mySonic.AddVariable("MET",            "MET",                     "GeV",    15,50,500,   &(myEvent.MET),           "logY=true");
 	 mySonic.AddVariable("MT",             "MT",                      "GeV",    17,0,510,    &(myEvent.MT),            "logY=true");
 	 mySonic.AddVariable("DPhi_MetJets",   "#Delta#Phi(MET,j_{1,2})", "rad",    16,0,3.2,    &(myEvent.DPhi_MetJets),  "logY=true");
 	 mySonic.AddVariable("HadronicChi2",   "had. #Chi^{2}",           "",       20,0,20,     &(myEvent.HadronicChi2),      "");
 	 mySonic.AddVariable("MT2W",           "M_{T2}^{W}",              "GeV",    20,0,500,    &(myEvent.MT2W),              "");
 	 //mySonic.AddVariable("mStop",        "m_{#tilde{t}}",           "GeV",    29,87.5,812.5,    &(myEvent.mStop),        "");
 	 //mySonic.AddVariable("mNeutralino",  "m_{#chi^0}",              "GeV",    17,-12.5,412.5,   &(myEvent.mNeutralino),  "");
 	 mySonic.AddVariable("mStop",          "m_{#tilde{t}}",           "GeV",    14,112.5,812.5,   &(myEvent.mStop),        "");
 	 mySonic.AddVariable("mNeutralino",    "m_{#chi^{0}}",              "GeV",    8,-12.5,387.5,    &(myEvent.mNeutralino),  "");
 	 mySonic.AddVariable("nJets",          "# of jets",               "",       5,2.5,7.5,        &(myEvent.nJets),        "");
     mySonic.AddVariable("genW_Pt",        "p_{T}(gen. W)",           "GeV",    30,0,600,         &(myEvent.genW_pT),      "");
 	 
     // AK8 collection 
     float nJetsAK8;            mySonic.AddVariable("nJetsAK8",        "# of AK8 jets (p_{T}>20)",       "", 10,-0.5,9.5,   &(nJetsAK8),          "");
     float nJetsAK8_tagged;     mySonic.AddVariable("nJetsAK8_tagged", "# of AK8 W-tagged jets",         "", 4,-0.5,3.5,    &(nJetsAK8_tagged),   "");
     float nJetsAK8_highPt;     mySonic.AddVariable("nJetsAK8_highPt", "# of AK8 jets (p_{T}>200)",      "", 6,-0.5,5.5,    &(nJetsAK8_highPt),   "");
     float nJetsAK8_matched;    mySonic.AddVariable("nJetsAK8_matched","# of AK8 jets matched to gen W", "", 3,-0.5,2.5,    &(nJetsAK8_matched),  "");
     
     float highPtGenWFlag;      mySonic.AddVariable("highPtGenWFlag",  "Event has an hadronic W with pT > 200",  "", 2,-0.5,1.5,    &(highPtGenWFlag),          "");
     float WtagFlag;            mySonic.AddVariable("WtagFlag",        "Event has an AK8 jet W-tagged",  "", 2,-0.5,1.5,    &(WtagFlag),          "");
     
  // #########################################################
  // ##   Create ProcessClasses (and associated datasets)   ##
  // #########################################################
      
      // ttbar -> l + jets
      //mySonic.AddProcessClass("ttbarSemiLept",    "t#bar{t} #rightarrow l^{#pm} + jets",     "background",COLORPLOT_RED);
      mySonic.AddProcessClass("1lepTop",          "1l top",     "background",kRed-7);

          mySonic.AddDataset("ttbarSemiLept",     "1lepTop",100,234*(0.67*0.33*2) );
          mySonic.AddDataset("Singletop_T_t",     "1lepTop",100,56.4);
          mySonic.AddDataset("Singletop_T_tW",    "1lepTop",100,11.1);
          mySonic.AddDataset("Singletop_T_s",     "1lepTop",100,3.79);
          mySonic.AddDataset("Singletop_Tbar_t",  "1lepTop",100,30.7);
       // mySonic.AddDataset("Singletop_Tbar_tW", "Singletop",100,11.1);
          mySonic.AddDataset("Singletop_Tbar_s",  "1lepTop",100,1.76);    
 
      // ttbar -> ll
      mySonic.AddProcessClass("ttbarFullLept",    "t#bar{t} #rightarrow l^{+}l^{-}",         "background",kCyan-3);

          mySonic.AddDataset("ttbarFullLept",     "ttbarFullLept",100,234*(0.33*0.33)   );

      // Wjets
      mySonic.AddProcessClass("W+jets",             "W+jets",                               "background",kOrange-2);
      
          mySonic.AddDataset("W2Jets",             "W+jets",100,2159);
          mySonic.AddDataset("W3Jets",             "W+jets",100,640);
          mySonic.AddDataset("W4Jets",             "W+jets",100,264);

      // Rare
      mySonic.AddProcessClass("rare",              "rare",                                   "background",kMagenta-5);
          
          mySonic.AddDataset("DY4Jets",           "rare",100,27.59    );
          mySonic.AddDataset("WGstarToLNu2E",     "rare",100,5.873    );
          mySonic.AddDataset("WWJetsTo2L2Nu",     "rare",100,5.8123   );
          mySonic.AddDataset("TBZ",               "rare",100,0.0114   );         
       // mySonic.AddDataset("TTG",               "rare",100,2.166    );         
          mySonic.AddDataset("TTW",               "rare",100,0.23     );
          mySonic.AddDataset("TTWW",              "rare",100,0.002037 );
          mySonic.AddDataset("TTZ",               "rare",100,0.2057   );
       // mySonic.AddDataset("WZJetsTo3LNu",      "rare",100,1.0575   );
          mySonic.AddDataset("WZJetsTo2L2Q",      "rare",100,2.206    );
          mySonic.AddDataset("WGstarToLNu2Mu",    "rare",100,1.914    );
          mySonic.AddDataset("WGstarToLNu2Tau",   "rare",100,0.336    ); 
       // mySonic.AddDataset("ZZJetsTo2L2Nu",     "rare",100,0.365    );      
          mySonic.AddDataset("WWWJets",           "rare",100,0.08058  );
          mySonic.AddDataset("WWZNoGstarJets",    "rare",100,0.05798  );
          mySonic.AddDataset("WWGJets",           "rare",100,0.528    ); 
       // mySonic.AddDataset("WZZNoGstarJets",    "rare",100,0.01698  );
       // mySonic.AddDataset("ZZJetsTo2L2Q",      "rare",100,2.4487   );
       // mySonic.AddDataset("ZZJetsTo4L",        "rare",100,0.1769   );
          mySonic.AddDataset("ZZZNoGstarJets",    "rare",100,0.0055269);
     
     // Signal
     mySonic.AddProcessClass("signal", "signal",  "signal",COLORPLOT_AZURE);
          mySonic.AddDataset("signal", "signal", 100,1);

  // ##########################
  // ##    Create Regions    ##
  // ##########################

     mySonic.AddRegion("presel",                "Preselection",             &Selector_presel);
     mySonic.AddRegion("baseline",              "Baseline",                 &Selector_baseline);
 
     // Reference
     mySonic.AddRegion("lowdM_MET>200",         "low DeltaM | MET > 200",   &Selector_lowdM_200);
     mySonic.AddRegion("highdM_MET>200",        "high DeltaM | MET > 200",  &Selector_highdM_200);

     mySonic.AddRegion("highdMSignal",          "highdMSignal",             &Selector_highdMSignal);
     
     // >= 3 jets
     mySonic.AddRegion(">=3jets",               "",                         &Selector_3jets);
     mySonic.AddRegion(">=3jets+Wtag",          "",                         &Selector_3jetsPlusWtag);
     mySonic.AddRegion(">=3jets+vetoWtag",      "",                         &Selector_3jetsPlusVetoWtag);
     
     // = 3 jets
     mySonic.AddRegion("=3jets",                "",                         &Selector_3jetsExactly);
     mySonic.AddRegion("=3jets+Wtag",           "",                         &Selector_3jetsExactlyPlusWtag);
     mySonic.AddRegion("=3jets+vetoWtag",       "",                         &Selector_3jetsExactlyPlusVetoWtag);
     
     // >= 4 jets
     mySonic.AddRegion(">=4jets",               "",                         &Selector_4jets);
     mySonic.AddRegion(">=4jets+Wtag",          "",                         &Selector_4jetsPlusWtag);
     mySonic.AddRegion(">=4jets+vetoWtag",      "",                         &Selector_4jetsPlusVetoWtag);
     
     // for info..
     mySonic.AddRegion(">=3jets+highPtGenW",    "",                         &Selector_3jetsPlusHighPtGenW);
     
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
     mySonic.Add2DHisto("nJetsAK8_highPt","nJetsAK8_matched");
  	 
     mySonic.Add3DHisto("mStop","mNeutralino","WtagFlag");
     mySonic.Add3DHisto("mStop","mNeutralino","highPtGenWFlag");
     mySonic.Add3DHisto("mStop","mNeutralino","genW_Pt");
     mySonic.Add3DHisto("mStop","mNeutralino","nJets");

  	 // Schedule plots
  	 mySonic.SchedulePlots("1DSuperpRenorm");
  	 mySonic.SchedulePlots("1DStack");
  	 mySonic.SchedulePlots("2D");
     
     //mySonic.SchedulePlots("2DProjectedTo1D","varX=genW_Pt,varY=genW_recoFlag,projectionType=mean,labelY=TaggingEfficiency");
     //mySonic.SchedulePlots("2DProjectedTo1D","varX=fakeW_Pt,varY=fakeW_tagFlag,projectionType=mean,labelY=MistaggingEfficiency");

     mySonic.SchedulePlots("3DProjectedTo2D","varX=mStop,varY=mNeutralino,varZ=genW_Pt,projectionType=mean,labelZ=mean(genW_Pt)");
     mySonic.SchedulePlots("3DProjectedTo2D","varX=mStop,varY=mNeutralino,varZ=highPtGenWFlag,projectionType=mean,labelZ=frac(genW_highPt)");
     mySonic.SchedulePlots("3DProjectedTo2D","varX=mStop,varY=mNeutralino,varZ=nJets,projectionType=mean,labelZ=mean(nJets)");
     mySonic.SchedulePlots("3DProjectedTo2D","varX=mStop,varY=mNeutralino,varZ=WtagFlag,projectionType=mean,labelZ=W-Tagging_Efficiency");
//     mySonic.SchedulePlots("3DProjectedTo2D","varX=mStop,varY=mNeutralino,varZ=genW_Pt,projectionType=fracSupTo200,labelZ=frac(genW_Pt>200)");
  
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
//         if (!(Selector_highdMSignal())) continue;
         weight = mySonic.GetLumi() * stopCrossSection(myEvent.mStop) 
                 / signalInitialNumberOfEvents->GetBinContent(  
                          signalInitialNumberOfEvents->FindBin(myEvent.mStop,myEvent.mNeutralino)  
                                                             );
         weight *= 0.9; // Trigger efficiency (trigger were not included for selection on signal when microtupling)
         weight *= 0.25; // Divide yields by 4, accounting for bin merging (cf declaration of mStop and mNeutralino variables)
      }

      // Dirty management of had. chi2 overflow
      if (myEvent.HadronicChi2 >= 19.99) myEvent.HadronicChi2 = 19.99;

      if (myEvent.genW_pT > 200) highPtGenWFlag = 1.0;
      else highPtGenWFlag = 0.0;

      // Tagging/Mistagging efficiencies
      WtagFlag = 0.0;
      nJetsAK8        = 0;       
      nJetsAK8_tagged = 0;       
      nJetsAK8_highPt = 0;
      nJetsAK8_matched = 0;
      for (int j = 0 ; j < 10 ; j++)
      {
        if (myEventPointer->recoW_pT[j] == 0.0) break;
        else nJetsAK8++;

        if (myEventPointer->recoW_pT[j] > 200) nJetsAK8_highPt++;
        float dRlep = deltaR(myEvent.lep_eta,
                             myEvent.lep_phi,
                             myEvent.recoW_Eta[j],
                             myEvent.recoW_Phi[j]);
        float dR = deltaR(myEvent.genW_Eta,
                          myEvent.genW_Phi,
                          myEvent.recoW_Eta[j],
                          myEvent.recoW_Phi[j]);
        bool tagged;
        if (dR < 0.1) tagged = WTagged(myEvent.recoW_Mass[j],myEvent.recoW_pT[j],dRlep);
        else          tagged = WTagged_Randomize(myEvent.recoW_Mass[j],myEvent.recoW_pT[j],dRlep);
        if (tagged) 
        {
           nJetsAK8_tagged++;
           WtagFlag = 1.0;
        }
        if ((myEventPointer->recoW_pT[j] > 200) && (dR < 0.1)) nJetsAK8_matched++;
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
  mySonic.WritePlots("plots/WTagStudies/",infoText);  //"exportPngAndEps=true");

  cout << endl;
  cout << "   ┌──────────────────────────────┐  " << endl;
  cout << "   │   Plot generation completed  │  " << endl;
  cout << "   └──────────────────────────────┘  " << endl; 
  cout << endl;
 
  cout << endl;
  cout << "   ┌──────────────────────────────┐  " << endl;
  cout << "   │  Starting tables generation  │  " << endl;
  cout << "   └──────────────────────────────┘  " << endl; 
  cout << endl;

      vector<string> processClassList;
      mySonic.GetProcessClassList(&processClassList);
      processClassList.pop_back();
      processClassList.push_back("total");
      processClassList.push_back("signal");

      vector<string> referenceSelections;
      referenceSelections.push_back("lowdM_MET>200");
      referenceSelections.push_back("highdM_MET>200");
      
      vector<string> threeJetsSelections;
      threeJetsSelections.push_back(">=3jets");
      threeJetsSelections.push_back(">=3jets+Wtag");
      threeJetsSelections.push_back(">=3jets+vetoWtag");
      
      vector<string> threeJetsExactlySelections;
      threeJetsExactlySelections.push_back("=3jets");
      threeJetsExactlySelections.push_back("=3jets+Wtag");
      threeJetsExactlySelections.push_back("=3jets+vetoWtag");
      
      vector<string> fourJetsSelections;
      fourJetsSelections.push_back(">=4jets");
      fourJetsSelections.push_back(">=4jets+Wtag");
      fourJetsSelections.push_back(">=4jets+vetoWtag");

      Table reference(referenceSelections,processClassList);
      Table threeJets(threeJetsSelections,processClassList);
      Table threeJetsExactly(threeJetsExactlySelections,processClassList);
      Table fourJets(fourJetsSelections,processClassList);
      
      fillMCSignalTable(&mySonic,referenceSelections,processClassList,&reference);
      fillMCSignalTable(&mySonic,threeJetsSelections,processClassList,&threeJets);
      fillMCSignalTable(&mySonic,threeJetsExactlySelections,processClassList,&threeJetsExactly);
      fillMCSignalTable(&mySonic,fourJetsSelections,processClassList,&fourJets);

      reference.PrintTable();
      threeJets.PrintTable();
      threeJetsExactly.PrintTable();
      fourJets.PrintTable();

    //
    // ROOT files for Eric
    //

    TH2F* signal2D_exactly3jets;
    TH2F* signal2D_exactly3jetsPlusWtag;
    TH2F* signal2D_exactly3jetsPlusVetoWtag;
    TH2F* signal2D_atLeast4jets;
    TH2F* signal2D_atLeast4jetsPlusWtag;
    TH2F* signal2D_atLeast4jetsPlusVetoWtag;
    TH2F* signal2D_atLeast3jets;

    vector<Histo2D>* vectorOfHisto2D = mySonic.Get2DHistoList();
    for (unsigned int i = 0 ; i < vectorOfHisto2D->size() ; i++)
    {
        Histo2D currentHisto = (*vectorOfHisto2D)[i];
          
        if ((currentHisto.getProcessClassTag() == "signal")
         && (currentHisto.getVariableXTag() == "mStop"	  )
         && (currentHisto.getVariableYTag() == "mNeutralino"))
        {
                if (currentHisto.getRegionTag() == "=3jets")           signal2D_exactly3jets             = currentHisto.getClone();
                if (currentHisto.getRegionTag() == "=3jets+Wtag")      signal2D_exactly3jetsPlusWtag     = currentHisto.getClone();
                if (currentHisto.getRegionTag() == "=3jets+vetoWtag")  signal2D_exactly3jetsPlusVetoWtag = currentHisto.getClone();
                if (currentHisto.getRegionTag() == ">=4jets")          signal2D_atLeast4jets             = currentHisto.getClone();
                if (currentHisto.getRegionTag() == ">=4jets+Wtag")     signal2D_atLeast4jetsPlusWtag     = currentHisto.getClone();
                if (currentHisto.getRegionTag() == ">=4jets+vetoWtag") signal2D_atLeast4jetsPlusVetoWtag = currentHisto.getClone();
                if (currentHisto.getRegionTag() == ">=3jets")          signal2D_atLeast3jets             = currentHisto.getClone();
        }
    }

    unsigned int nBinsXSignal2D = signal2D_exactly3jets->GetNbinsX();
    unsigned int nBinsYSignal2D = signal2D_exactly3jets->GetNbinsY();
    
    for (unsigned int i = 1 ; i <= nBinsXSignal2D ; i++)
    for (unsigned int j = 1 ; j <= nBinsXSignal2D ; j++)
    {
       // Get current signal point
       
       float mStopMean       = signal2D_exactly3jets->GetXaxis()->GetBinCenter(i);
       float mNeutralinoMean = signal2D_exactly3jets->GetYaxis()->GetBinCenter(j);

       std::ostringstream fileName;
       fileName << "scenarios_mStop" << mStopMean << "_mNeutralino" << mNeutralinoMean << ".root";

       // Open file
       
       TFile scenariosSave(fileName.str().c_str(),"RECREATE");

       // Create histo

       TH1F* oneBox_AtLeast4Jets__bkg = new TH1F("1box_>=4Jets__bkg","",4,0.5,4.5);
       TH1F* oneBox_AtLeast4Jets__sig = new TH1F("1box_>=4Jets__sig","",4,0.5,4.5);
     
       TH1F* oneBox_AtLeast3Jets__bkg = new TH1F("1box_>=3Jets__bkg","",4,0.5,4.5);
       TH1F* oneBox_AtLeast3Jets__sig = new TH1F("1box_>=3Jets__sig","",4,0.5,4.5);
      
       TH1F* twoBox_Exactly3Jets_vs_AtLeast4Jets__bkg = new TH1F("2box_=3Jets_vs_>=4Jets__bkg","",4,0.5,4.5);
       TH1F* twoBox_Exactly3Jets_vs_AtLeast4Jets__sig = new TH1F("2box_=3Jets_vs_>=4Jets__sig","",4,0.5,4.5);
       
       TH1F* fourBox_Exactly3Jets_AtLeast4Jets_vs_0Wtag_1Wtag__bkg = new TH1F("4box_=3Jets_>=4Jets_vs_0Wtag_1Wtag__bkg","",4,0.5,4.5);
       TH1F* fourBox_Exactly3Jets_AtLeast4Jets_vs_0Wtag_1Wtag__sig = new TH1F("4box_=3Jets_>=4Jets_vs_0Wtag_1Wtag__sig","",4,0.5,4.5);
      
       // Fill

       oneBox_AtLeast4Jets__bkg->Fill(1.0,fourJets.Get(">=4jets","total").value());
       oneBox_AtLeast4Jets__sig->Fill(1.0,signal2D_atLeast4jets->GetBinContent(i,j));

       oneBox_AtLeast3Jets__bkg->Fill(1.0,threeJets.Get(">=3jets","total").value());
       oneBox_AtLeast3Jets__sig->Fill(1.0,signal2D_atLeast3jets->GetBinContent(i,j));

       twoBox_Exactly3Jets_vs_AtLeast4Jets__bkg->Fill(1.0,threeJetsExactly.Get("=3jets","total").value());
       twoBox_Exactly3Jets_vs_AtLeast4Jets__sig->Fill(1.0,signal2D_exactly3jets->GetBinContent(i,j));
       twoBox_Exactly3Jets_vs_AtLeast4Jets__bkg->Fill(2.0,fourJets.Get(">=4jets","total").value());
       twoBox_Exactly3Jets_vs_AtLeast4Jets__sig->Fill(2.0,signal2D_atLeast4jets->GetBinContent(i,j));

       fourBox_Exactly3Jets_AtLeast4Jets_vs_0Wtag_1Wtag__bkg->Fill(1.0,threeJetsExactly.Get("=3jets+vetoWtag","total").value());
       fourBox_Exactly3Jets_AtLeast4Jets_vs_0Wtag_1Wtag__sig->Fill(1.0,signal2D_exactly3jetsPlusVetoWtag->GetBinContent(i,j));
       fourBox_Exactly3Jets_AtLeast4Jets_vs_0Wtag_1Wtag__bkg->Fill(2.0,threeJetsExactly.Get("=3jets+Wtag","total").value());
       fourBox_Exactly3Jets_AtLeast4Jets_vs_0Wtag_1Wtag__sig->Fill(2.0,signal2D_exactly3jetsPlusWtag->GetBinContent(i,j));
       fourBox_Exactly3Jets_AtLeast4Jets_vs_0Wtag_1Wtag__bkg->Fill(3.0,fourJets.Get(">=4jets+vetoWtag","total").value());
       fourBox_Exactly3Jets_AtLeast4Jets_vs_0Wtag_1Wtag__sig->Fill(3.0,signal2D_atLeast4jetsPlusVetoWtag->GetBinContent(i,j));
       fourBox_Exactly3Jets_AtLeast4Jets_vs_0Wtag_1Wtag__bkg->Fill(4.0,fourJets.Get(">=4jets+Wtag","total").value());
       fourBox_Exactly3Jets_AtLeast4Jets_vs_0Wtag_1Wtag__sig->Fill(4.0,signal2D_atLeast4jetsPlusWtag->GetBinContent(i,j));

       // Write

       oneBox_AtLeast4Jets__bkg->Write();
       oneBox_AtLeast4Jets__sig->Write();
     
       oneBox_AtLeast3Jets__bkg->Write();
       oneBox_AtLeast3Jets__sig->Write();
      
       twoBox_Exactly3Jets_vs_AtLeast4Jets__bkg->Write();
       twoBox_Exactly3Jets_vs_AtLeast4Jets__sig->Write();
       
       fourBox_Exactly3Jets_AtLeast4Jets_vs_0Wtag_1Wtag__bkg->Write();
       fourBox_Exactly3Jets_AtLeast4Jets_vs_0Wtag_1Wtag__sig->Write();
      
       // Close

       scenariosSave.Close();
    }

  cout << endl;
  cout << "   ┌───────────────────────────────┐  " << endl;
  cout << "   │  Tables generation completed  │  " << endl;
  cout << "   └───────────────────────────────┘  " << endl; 
  cout << endl;

  cout << endl;
  cout << "   ┌────────────────────────┐  " << endl;
  cout << "   │  Computing ratio plot  │  " << endl;
  cout << "   └────────────────────────┘  " << endl; 
  cout << endl;

  // #########################
  // ##   Make ratio plot   ##
  // #########################
 
  TH2F* boostedSelectionSignalYield;
  TH2F* standardSelectionSignalYield;
  TH2F* ratioSelectionSignalYield;

  //vector<Histo2D>* vectorOfHisto2D = mySonic.Get2DHistoList();
  for (unsigned int i = 0 ; i < vectorOfHisto2D->size() ; i++)
  {
      Histo2D currentHisto = (*vectorOfHisto2D)[i];
      
      if ((currentHisto.getProcessClassTag() == "signal")
       && (currentHisto.getVariableXTag() == "mStop"	  )
       && (currentHisto.getVariableYTag() == "mNeutralino"))
      {
              if (currentHisto.getRegionTag() == "highdM_MET>200")
                  standardSelectionSignalYield = currentHisto.getClone();
         else if (currentHisto.getRegionTag() == ">=3jets+Wtag")
         {              
                  boostedSelectionSignalYield = currentHisto.getClone();
                  ratioSelectionSignalYield   = currentHisto.getClone();
         }
      }
  }

  ratioSelectionSignalYield->Divide(standardSelectionSignalYield);

  Figure boostedSelectionBackgroundYield = threeJets.Get(">=3jets+Wtag","total");
  Figure standardSelectionBackgroundYield = reference.Get("highdM_MET>200","total");

  Figure ratioSelectionBackgroundYield = standardSelectionBackgroundYield / boostedSelectionBackgroundYield;

  ratioSelectionSignalYield->Scale(sqrt(ratioSelectionBackgroundYield.value()));

  TFile f("WTagStudyRatio.root","RECREATE");
  boostedSelectionSignalYield->Write();
  standardSelectionSignalYield->Write();
  ratioSelectionSignalYield->Write();
  f.Close();

  TFile f2("signalInitialNumberOfEvents.root","RECREATE");
  signalInitialNumberOfEvents->Write();
  f2.Close();

  // ####################################
  // ##   Save raw plot for ext. use   ##
  // ####################################

  TFile f3("multiplicityAK8_highPt_vs_matched.root","RECREATE");
  TH2F* clone = 0;
  vectorOfHisto2D = mySonic.Get2DHistoList();
  for (unsigned int i = 0 ; i < vectorOfHisto2D->size() ; i++)
  {
      Histo2D currentHisto = (*vectorOfHisto2D)[i];
      
      if ((currentHisto.getRegionTag() == ">=3jets")
       && (currentHisto.getVariableXTag() == "nJetsAK8_highPt" )
       && (currentHisto.getVariableYTag() == "nJetsAK8_matched"))
      {
         string process = currentHisto.getProcessClassTag();
         clone = currentHisto.getClone();
         clone->Write(process.c_str());
      }
  }
  f3.Close();


  return (0);

}

void fillMCSignalTable(SonicScrewdriver* mySonic, vector<string> region, vector<string> process, Table* table)
{
    string varUsedToGetYields = "HadronicChi2";
    string channelUsedToGetYields = "inclusiveChannel";

    for (unsigned int r = 0 ; r < region.size()          ; r++)
    {
        Figure tmpTotal(0.0,0.0);
        for (unsigned int p = 0 ; p < process.size() ; p++)
        {
            if (process[p] == "total") continue;
            table->Set(region[r],
                      process[p],
                      mySonic->GetYieldAndError(varUsedToGetYields,
                                               process[p],
                                               region[r],
                                               channelUsedToGetYields));

            if (process[p] != "signal")
                tmpTotal += mySonic->GetYieldAndError(varUsedToGetYields,
                                                     process[p],
                                                     region[r],
                                                     channelUsedToGetYields);
        }
        table->Set(region[r],"total",tmpTotal);
    }

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
