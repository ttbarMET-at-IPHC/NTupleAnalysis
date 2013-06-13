
#include <math.h>
#include <iomanip>
#include <iostream>
#include <time.h>

#include <TFile.h>
#include <TMarker.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>

using namespace std;

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

    // Variables for analysis 
    Float_t MET;
    Float_t MET_lep_dphi;
    Float_t MT;

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





TH1F* ratioForMETBeforeTauVeto;
TH1F* ratioForMETAfterTauVeto;

TH1F* ratioForMTBeforeTauVeto;
TH1F* ratioForMTAfterTauVeto;

float yieldDataBeforeTauVeto;
float yieldMCBeforeTauVeto;

float yieldDataAfterTauVeto;
float yieldMCAfterTauVeto;



#define FOLDER_MICROTUPLES "/opt/sbg/data/data4/cms/aaubin/tauVetoValidation_beta/store/microTuples314/"

#define ELECTRON_CHANNEL 1.0
#define MUON_CHANNEL 2.0

//#define LEPTON_CHANNEL ELECTRON_CHANNEL
#define LEPTON_CHANNEL MUON_CHANNEL


LumiReWeighting loadPileUp(string dataPileUp);
float triggerWeight(int flavor, float pt, float eta);

int main()
{

   //
   //   Method 2
   //
  
  std::vector<std::string> channelListMethod2;
  channelListMethod2.push_back("ll");
  channelListMethod2.push_back("others");
  channelListMethod2.push_back("total MC");
  channelListMethod2.push_back("data");

  std::vector<std::string> regionMethod2;
  regionMethod2.push_back("pre-veto1");
  regionMethod2.push_back("post-veto1");
  regionMethod2.push_back("eff1");
  regionMethod2.push_back("pre-veto2");
  regionMethod2.push_back("post-veto2");
  regionMethod2.push_back("eff2");

  TableScrew tableMethod2(regionMethod2,channelListMethod2);

  string region1pre("MTpeak-preTauVeto-region1");
  string region1post("MTpeak-postTauVeto-region1");
  
  string region2pre("MTpeak-preTauVeto-region2");
  string region2post("MTpeak-postTauVeto-region2");

  FigureScrew yieldPre1 = 0;
  FigureScrew yieldPost1 = 0;
  FigureScrew yieldPre2 = 0;
  FigureScrew yieldPost2 = 0;
  
  FigureScrew yieldPre1TotMC = 0;
  FigureScrew yieldPost1TotMC = 0;
  FigureScrew yieldPre2TotMC = 0;
  FigureScrew yieldPost2TotMC = 0;

  FigureScrew tmpYield;

  float f  = 1;
  float fe = 1;
     
  tmpYield = FigureScrew(  246*f,1*fe);  tableMethod2.Set("pre-veto1" ,"ll",tmpYield); yieldPre1TotMC +=tmpYield;
  tmpYield = FigureScrew(  215*f,1*fe);  tableMethod2.Set("post-veto1","ll",tmpYield); yieldPost1TotMC+=tmpYield;
  tmpYield = FigureScrew(  938*f,2*fe);  tableMethod2.Set("pre-veto2" ,"ll",tmpYield); yieldPre2TotMC +=tmpYield;
  tmpYield = FigureScrew(  835*f,2*fe);  tableMethod2.Set("post-veto2","ll",tmpYield); yieldPost2TotMC+=tmpYield;

  tmpYield = FigureScrew( 4395*f,21*fe); tableMethod2.Set("pre-veto1" ,"others",tmpYield); yieldPre1TotMC  +=tmpYield;
  tmpYield = FigureScrew( 4332*f,21*fe); tableMethod2.Set("post-veto1","others",tmpYield); yieldPost1TotMC +=tmpYield;
  tmpYield = FigureScrew( 5445*f,18*fe); tableMethod2.Set("pre-veto2" ,"others",tmpYield); yieldPre2TotMC  +=tmpYield;
  tmpYield = FigureScrew( 5337*f,18*fe); tableMethod2.Set("post-veto2","others",tmpYield); yieldPost2TotMC +=tmpYield;

  tableMethod2.Set("pre-veto1" ,"total MC",yieldPre1TotMC);
  tableMethod2.Set("post-veto1","total MC",yieldPost1TotMC);
  tableMethod2.Set("pre-veto2" ,"total MC",yieldPre2TotMC);
  tableMethod2.Set("post-veto2","total MC",yieldPost2TotMC);

  tmpYield = FigureScrew(5374*f,73*fe); tableMethod2.Set("pre-veto1" ,"data",tmpYield);
  tmpYield = FigureScrew(5262*f,72*fe); tableMethod2.Set("post-veto1","data",tmpYield);
  tmpYield = FigureScrew(6806*f,75*fe); tableMethod2.Set("pre-veto2" ,"data",tmpYield);
  tmpYield = FigureScrew(6493*f,74*fe); tableMethod2.Set("post-veto2","data",tmpYield);

  for (unsigned int i = 0 ; i  < channelListMethod2.size() ; i++)
  {
    tableMethod2.Set("eff1",channelListMethod2[i],  tableMethod2.Get("post-veto1",channelListMethod2[i]) / tableMethod2.Get("pre-veto1",channelListMethod2[i]) );
    tableMethod2.Set("eff2",channelListMethod2[i],  tableMethod2.Get("post-veto2",channelListMethod2[i]) / tableMethod2.Get("pre-veto2",channelListMethod2[i]) );
  }

  tableMethod2.PrintTable();

  FigureScrew SF_1 = tableMethod2.Get("pre-veto1","data") / tableMethod2.Get("pre-veto1","total MC");
  FigureScrew SF_2 = tableMethod2.Get("pre-veto2","data") / tableMethod2.Get("pre-veto2","total MC");

  
  FigureScrew N_1_data = tableMethod2.Get("post-veto1","data"); 
  FigureScrew N_2_data = tableMethod2.Get("post-veto2","data"); 


  FigureScrew N_1_ll = tableMethod2.Get("pre-veto1","ll") * SF_1; 
  FigureScrew N_2_ll = tableMethod2.Get("pre-veto2","ll") * SF_2; 
  
  FigureScrew N_1_oth = tableMethod2.Get("pre-veto1","others") * SF_1; 
  FigureScrew N_2_oth = tableMethod2.Get("pre-veto2","others") * SF_2; 


  FigureScrew eps_others_fromData = ( N_2_data - (N_2_ll / N_1_ll)   * N_1_data ) / ( N_2_oth - (N_2_ll / N_1_ll)   * N_1_oth );
  FigureScrew eps_ll_fromData     = ( N_2_data - (N_2_oth / N_1_oth) * N_1_data ) / ( N_2_ll -  (N_2_oth / N_1_oth) * N_1_ll  );
  cout << "eps_others_data = " << eps_others_fromData.Print() << endl;
  cout << "eps_ll_data = " << eps_ll_fromData.Print() << endl;



  FigureScrew a(1.0,1.0);
  FigureScrew b(2.0,1.0);

  cout << "a-b : " << (a-b).Print() << endl;
  cout << "a+b : " << (a+b).Print() << endl;
  cout << "a*b : " << (a*b).Print() << endl;
  cout << "a/b : " << (a/b).Print() << endl;

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




