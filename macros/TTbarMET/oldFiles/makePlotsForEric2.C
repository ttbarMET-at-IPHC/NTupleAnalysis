

#include <iomanip>
#include <iostream>

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TMath.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <iomanip>


using namespace std;



float simulateWtagRequirement(float e_tag, float e_mis, int nMatched, int nNotMatched);

int main()
{
    TFile f("tmp.root","UPDATE");
  	
    TH3F* signal2D_exactly3jets;
    TH3F* signal2D_atLeast4jets;
    TH3F* signal2D_atLeast3jets;
    TH1F* totalSM_exactly3jets;
    TH1F* totalSM_atLeast4jets;
    TH1F* totalSM_atLeast3jets;

    f.GetObject("variableX:mStop|variableY:mNeutralino|variableZ:nJetsAK8_matched_notMatched|processclass:signal|region:=3jets|channel:inclusiveChannel",signal2D_exactly3jets);
    f.GetObject("variableX:mStop|variableY:mNeutralino|variableZ:nJetsAK8_matched_notMatched|processclass:signal|region:>=4jets|channel:inclusiveChannel",signal2D_atLeast4jets);
    f.GetObject("variableX:mStop|variableY:mNeutralino|variableZ:nJetsAK8_matched_notMatched|processclass:signal|region:>=3jets|channel:inclusiveChannel",signal2D_atLeast3jets);
    f.GetObject("variable:nJetsAK8_matched_notMatched|processclass:1lepTop|region:=3jets|channel:inclusiveChannel",totalSM_exactly3jets);
    f.GetObject("variable:nJetsAK8_matched_notMatched|processclass:1lepTop|region:>=4jets|channel:inclusiveChannel",totalSM_atLeast4jets);
    f.GetObject("variable:nJetsAK8_matched_notMatched|processclass:1lepTop|region:>=3jets|channel:inclusiveChannel",totalSM_atLeast3jets);




    unsigned int nBinsXSignal2D = signal2D_exactly3jets->GetNbinsX();
    unsigned int nBinsYSignal2D = signal2D_exactly3jets->GetNbinsY();
    
    for (unsigned int i = 1 ; i <= nBinsXSignal2D ; i++)
    for (unsigned int j = 1 ; j <= nBinsYSignal2D ; j++)
    {
       // Get current signal point
       float mStopMean       = signal2D_exactly3jets->GetXaxis()->GetBinCenter(i);
       float mNeutralinoMean = signal2D_exactly3jets->GetYaxis()->GetBinCenter(j);

       cout << endl;
       cout << "mStop, mNeutralino = " << mStopMean << "," << mNeutralinoMean << endl;
       cout << endl;

       float exactly3jets_signal = 0;
       float exactly3jets_backgr = 0;
       float atLeast3jets_signal = 0;
       float atLeast3jets_backgr = 0;
       float atLeast4jets_signal = 0;
       float atLeast4jets_backgr = 0;

       float exactly3jets_Wtag_signal = 0;                           float exactly3jets_vetoWtag_signal = 0;
       float exactly3jets_Wtag_backgr = 0;                           float exactly3jets_vetoWtag_backgr = 0;
       float atLeast3jets_Wtag_signal = 0;                           float atLeast3jets_vetoWtag_signal = 0;
       float atLeast3jets_Wtag_backgr = 0;                           float atLeast3jets_vetoWtag_backgr = 0;
       float atLeast4jets_Wtag_signal = 0;                           float atLeast4jets_vetoWtag_signal = 0;
       float atLeast4jets_Wtag_backgr = 0;                           float atLeast4jets_vetoWtag_backgr = 0;

       unsigned int nBinsZSignal2D = signal2D_exactly3jets->GetNbinsZ();
       for (unsigned int k = 1 ; k < nBinsZSignal2D ; k++)
       {    
           int nMatched_nNotMatched = signal2D_exactly3jets->GetZaxis()->GetBinCenter(k);
           int nMatched = nMatched_nNotMatched / 10;
           int nNotMatched = nMatched_nNotMatched % 10;

           float p = simulateWtagRequirement(0.8,0.23,nMatched,nNotMatched);

           exactly3jets_signal += signal2D_exactly3jets->GetBinContent(i,j,k) - signal2D_exactly3jets->GetBinError(i,j,k);
           exactly3jets_backgr += totalSM_exactly3jets->GetBinContent(k)      + totalSM_exactly3jets->GetBinError(k);
           atLeast3jets_signal += signal2D_atLeast3jets->GetBinContent(i,j,k) - signal2D_atLeast3jets->GetBinError(i,j,k);
           atLeast3jets_backgr += totalSM_atLeast3jets->GetBinContent(k)      + totalSM_atLeast3jets->GetBinError(k);
           atLeast4jets_signal += signal2D_atLeast4jets->GetBinContent(i,j,k) - signal2D_atLeast4jets->GetBinError(i,j,k);
           atLeast4jets_backgr += totalSM_atLeast4jets->GetBinContent(k)      + totalSM_atLeast4jets->GetBinError(k);

           exactly3jets_Wtag_signal += p * exactly3jets_signal;  exactly3jets_vetoWtag_signal += (1 - p) * exactly3jets_signal;
           exactly3jets_Wtag_backgr += p * exactly3jets_backgr;  exactly3jets_vetoWtag_backgr += (1 - p) * exactly3jets_backgr;
           atLeast3jets_Wtag_signal += p * atLeast3jets_signal;  atLeast3jets_vetoWtag_signal += (1 - p) * atLeast3jets_signal;
           atLeast3jets_Wtag_backgr += p * atLeast3jets_backgr;  atLeast3jets_vetoWtag_backgr += (1 - p) * atLeast3jets_backgr;
           atLeast4jets_Wtag_signal += p * atLeast4jets_signal;  atLeast4jets_vetoWtag_signal += (1 - p) * atLeast4jets_signal;
           atLeast4jets_Wtag_backgr += p * atLeast4jets_backgr;  atLeast4jets_vetoWtag_backgr += (1 - p) * atLeast4jets_backgr;
       }

       std::ostringstream fileName;
       fileName << "scenarios_mStop" << mStopMean << "_mNeutralino" << mNeutralinoMean << ".root";

       // Open file
       
       TFile scenariosSave(fileName.str().c_str(),"RECREATE");

       // Create histo

       TH1F* oneBox_AtLeast4Jets__bkg = new TH1F("1box_>=4Jets__bkg","",4,0.5,4.5);
       TH1F* oneBox_AtLeast4Jets__sig = new TH1F("1box_>=4Jets__sig","",4,0.5,4.5);

       TH1F* oneBox_AtLeast4JetsOr3JetsPlus1WTag__bkg = new TH1F("1box_>=4JetsOr3JetsPlus1WTag__bkg","",4,0.5,4.5);
       TH1F* oneBox_AtLeast4JetsOr3JetsPlus1WTag__sig = new TH1F("1box_>=4JetsOr3JetsPlus1WTag__sig","",4,0.5,4.5);

       TH1F* twoBox_AtLeast4Jets_Exactly3Jets__bkg = new TH1F("2box_>=4Jets_=3Jets__bkg","",4,0.5,4.5);
       TH1F* twoBox_AtLeast4Jets_Exactly3Jets__sig = new TH1F("2box_>=4Jets_=3Jets__sig","",4,0.5,4.5);
       
       TH1F* twoBox_AtLeast4Jets_Exactly3JetsPlus1WTag__bkg = new TH1F("2box_>=4Jets_=3JetsPlus1WTag__bkg","",4,0.5,4.5);
       TH1F* twoBox_AtLeast4Jets_Exactly3JetsPlus1WTag__sig = new TH1F("2box_>=4Jets_=3JetsPlus1WTag__sig","",4,0.5,4.5);
       
       TH1F* threeBox_AtLeast4JetsPlus0WTag_AtLeast4JetsPlus1Wtag_Exactly3JetsPlus1WTag__bkg = new TH1F("3box_>=4JetsPlus0WTag_>=4JetsPlus1WTag_=3JetsPlus1WTag__bkg","",4,0.5,4.5);
       TH1F* threeBox_AtLeast4JetsPlus0WTag_AtLeast4JetsPlus1Wtag_Exactly3JetsPlus1WTag__sig = new TH1F("3box_>=4JetsPlus0WTag_>=4JetsPlus1WTag_=3JetsPlus1WTag__sig","",4,0.5,4.5);
       
       TH1F* fourBox__bkg = new TH1F("4box__bkg","",4,0.5,4.5);
       TH1F* fourBox__sig = new TH1F("4box__sig","",4,0.5,4.5);
      
       // Fill

       oneBox_AtLeast4Jets__bkg->Fill(1.0,atLeast4jets_backgr);
       oneBox_AtLeast4Jets__sig->Fill(1.0,atLeast4jets_signal);

       oneBox_AtLeast4JetsOr3JetsPlus1WTag__bkg->Fill(1.0,atLeast4jets_backgr+exactly3jets_Wtag_backgr);
       oneBox_AtLeast4JetsOr3JetsPlus1WTag__sig->Fill(1.0,atLeast4jets_signal+exactly3jets_Wtag_signal);

       twoBox_AtLeast4Jets_Exactly3Jets__bkg->Fill(1.0,atLeast4jets_backgr);
       twoBox_AtLeast4Jets_Exactly3Jets__sig->Fill(1.0,atLeast4jets_signal);
       twoBox_AtLeast4Jets_Exactly3Jets__bkg->Fill(2.0,exactly3jets_backgr);
       twoBox_AtLeast4Jets_Exactly3Jets__sig->Fill(2.0,exactly3jets_signal);
       
       twoBox_AtLeast4Jets_Exactly3JetsPlus1WTag__bkg->Fill(1.0,atLeast4jets_backgr);
       twoBox_AtLeast4Jets_Exactly3JetsPlus1WTag__sig->Fill(1.0,atLeast4jets_signal);
       twoBox_AtLeast4Jets_Exactly3JetsPlus1WTag__bkg->Fill(2.0,exactly3jets_Wtag_backgr);
       twoBox_AtLeast4Jets_Exactly3JetsPlus1WTag__sig->Fill(2.0,exactly3jets_Wtag_signal);
 
       threeBox_AtLeast4JetsPlus0WTag_AtLeast4JetsPlus1Wtag_Exactly3JetsPlus1WTag__bkg->Fill(1.0,atLeast4jets_vetoWtag_backgr);
       threeBox_AtLeast4JetsPlus0WTag_AtLeast4JetsPlus1Wtag_Exactly3JetsPlus1WTag__sig->Fill(1.0,atLeast4jets_vetoWtag_signal);
       threeBox_AtLeast4JetsPlus0WTag_AtLeast4JetsPlus1Wtag_Exactly3JetsPlus1WTag__bkg->Fill(2.0,atLeast4jets_Wtag_backgr);
       threeBox_AtLeast4JetsPlus0WTag_AtLeast4JetsPlus1Wtag_Exactly3JetsPlus1WTag__sig->Fill(2.0,atLeast4jets_Wtag_signal);
       threeBox_AtLeast4JetsPlus0WTag_AtLeast4JetsPlus1Wtag_Exactly3JetsPlus1WTag__bkg->Fill(3.0,exactly3jets_Wtag_backgr);
       threeBox_AtLeast4JetsPlus0WTag_AtLeast4JetsPlus1Wtag_Exactly3JetsPlus1WTag__sig->Fill(3.0,exactly3jets_Wtag_signal);
      
       fourBox__bkg->Fill(1.0,atLeast4jets_vetoWtag_backgr);
       fourBox__sig->Fill(1.0,atLeast4jets_vetoWtag_signal);
       fourBox__bkg->Fill(2.0,atLeast4jets_Wtag_backgr);
       fourBox__sig->Fill(2.0,atLeast4jets_Wtag_signal);
       fourBox__bkg->Fill(3.0,exactly3jets_vetoWtag_backgr);
       fourBox__sig->Fill(3.0,exactly3jets_vetoWtag_signal);
       fourBox__bkg->Fill(4.0,exactly3jets_Wtag_backgr);
       fourBox__sig->Fill(4.0,exactly3jets_Wtag_signal);

      // Write

       oneBox_AtLeast4Jets__bkg->Write();
       oneBox_AtLeast4Jets__sig->Write();

       oneBox_AtLeast4JetsOr3JetsPlus1WTag__bkg->Write();
       oneBox_AtLeast4JetsOr3JetsPlus1WTag__sig->Write();

       twoBox_AtLeast4Jets_Exactly3Jets__bkg->Write();
       twoBox_AtLeast4Jets_Exactly3Jets__sig->Write();
       
       twoBox_AtLeast4Jets_Exactly3JetsPlus1WTag__bkg->Write();
       twoBox_AtLeast4Jets_Exactly3JetsPlus1WTag__sig->Write();
       
       threeBox_AtLeast4JetsPlus0WTag_AtLeast4JetsPlus1Wtag_Exactly3JetsPlus1WTag__bkg->Write();
       threeBox_AtLeast4JetsPlus0WTag_AtLeast4JetsPlus1Wtag_Exactly3JetsPlus1WTag__sig->Write();
       
       fourBox__bkg->Write();
       fourBox__sig->Write();
      
      // Close

       scenariosSave.Close();

    }
}


float simulateWtagRequirement(float e_tag, float e_mis, int nMatched, int nNotMatched)
{
    if (nMatched >= 2) 
    {
        cout << "Error : nMatched >= 2" << endl; 
        return 0;
    }

    if (nMatched == 0)
    {
        return 
            TMath::Binomial(nNotMatched,1) *    e_mis    * pow(1 - e_mis , nNotMatched-1)  // P( 1W-tag | 0matched + n NotMatched)
          + TMath::Binomial(nNotMatched,2) * e_mis*e_mis * pow(1 - e_mis , nNotMatched-2); // P( 2W-tag | 0matched + n NotMatched)
    }
    else
    {
        return
                                                e_tag                  * pow(1 - e_mis , nNotMatched)
          + TMath::Binomial(nNotMatched,1) * (1 - e_tag) *   e_mis     * pow(1 - e_mis , nNotMatched-1)  // P(1W-tag | 1matched + nNotMatched)
          + TMath::Binomial(nNotMatched,1) *    e_tag    *   e_mis     * pow(1 - e_mis , nNotMatched-1)
          + TMath::Binomial(nNotMatched,2) * (1 - e_tag) * e_mis*e_mis * pow(1 - e_mis , nNotMatched-2); // P(2W-tag | 1matched + nNotMatched)
    }
}
