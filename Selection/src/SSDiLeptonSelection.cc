#include "Selection/interface/SSDiLeptonSelection.h"
using namespace std;


// ----------------------------------------------------------------------------
// Default constructor
// ----------------------------------------------------------------------------
SSDiLeptonSelection::SSDiLeptonSelection ()
{
  selCode_        = -1;
  MinMassCut_     = 0.;
  ZMassWindow_    = std::make_pair(9999., 0.);
  METCuts_        = std::make_pair(9999., 0.);
  btagAlgo_       = -1;
  btagDiscriCut_  = -999.;
  NofBtagJets_    = 0;

  //Fill the table of cuts
  cuts_.push_back ("All");
  cuts_.push_back ("All dilept");
  cuts_.push_back ("Trigger");
  cuts_.push_back ("DiLepton SS pair");
  cuts_.push_back ("Z mass veto");
  cuts_.push_back ("NJets cut");
  cuts_.push_back ("MET cuts");
  cuts_.push_back ("NbtagJets cut");

  //Fill Channels
  channels_.push_back (string ("ee_ss"));
  channels_.push_back (string ("emu_ss"));
  channels_.push_back (string ("mumu_ss"));
  channels_.push_back (string ("ee_os"));
  channels_.push_back (string ("emu_os"));
  channels_.push_back (string ("mumu_os"));
}


// ----------------------------------------------------------------------------
// Copy constructor
// ----------------------------------------------------------------------------
SSDiLeptonSelection::SSDiLeptonSelection(const SSDiLeptonSelection & s):
                                                                 Selection (s)
{
  selCode_       = s.selCode_;
  MinMassCut_    = s.MinMassCut_;
  METCuts_       = s.METCuts_;
  ZMassWindow_   = s.ZMassWindow_;
  cuts_          = s.cuts_;
  channels_      = s.channels_;
  btagAlgo_      = s.btagAlgo_;
  btagDiscriCut_ = s.btagDiscriCut_;
  NofBtagJets_   = s.NofBtagJets_;
}


// ----------------------------------------------------------------------------
// SetParameters
// ----------------------------------------------------------------------------
void SSDiLeptonSelection::SetParameters(float MinValue,
                                        std::pair<float,float> METCuts,
                                        std::pair<float,float> ZMassWindow,
                                        int btagAlgo,
                                        float btagDiscriCut,
                                        int NofBtagJets)
{
  MinMassCut_    = MinValue;
  METCuts_       = METCuts;
  ZMassWindow_   = ZMassWindow;
  btagAlgo_      = btagAlgo;
  btagDiscriCut_ = btagDiscriCut;
  NofBtagJets_   = NofBtagJets;
}



bool SSDiLeptonSelection::GetLeptonPair (const std::vector < IPHCTree::NTMuon >& muon_in, 
                                         const std::vector < IPHCTree::NTElectron >& elec_in, 
                                         std::vector < NTMuon > &muon_out,
                                         std::vector < NTElectron > &elec_out,
                                         std::string & CandPairType, 
					 bool isForMM, float iso1_in, float iso2_in)
{
  // isForMM -> apply different isolation criteria on both lepton

  //important: reset the out collections
  muon_out.clear ();
  elec_out.clear ();

  // Electron combination
  float sum_pT_ee = 0.;
  bool pass_elec = false;
  int ie1 = -1;
  int ie2 = -1;
  if (elec_in.size () >= 2)
  {
    for (unsigned int i = 0; i < elec_in.size (); i++) 
    {
      for (unsigned int j = i + 1; j < elec_in.size (); j++) 
      {
        if (pass_elec) continue;
        if ((elec_in[i].charge == elec_in[j].charge))
	if( !isForMM || (isForMM && 
	 ((elec_in[i].RelIso03PF()<iso1_in && elec_in[j].RelIso03PF()<iso2_in) ||
	 (elec_in[i].RelIso03PF()<iso2_in && elec_in[j].RelIso03PF()<iso1_in))))
        {
          pass_elec = true;
          sum_pT_ee = elec_in[i].p4.Pt () + elec_in[j].p4.Pt ();
          ie1 = i;
          ie2 = j;
        }
      }
    }
    if(!pass_elec) // No SS pair found
    if( !isForMM || (isForMM && 
     ((elec_in[0].RelIso03PF()<iso1_in && elec_in[1].RelIso03PF()<iso2_in) ||
     (elec_in[0].RelIso03PF()<iso2_in && elec_in[1].RelIso03PF()<iso1_in))))
    { 
      sum_pT_ee = elec_in[0].p4.Pt () + elec_in[1].p4.Pt ();
      ie1=0;
      ie2=1;
    }
  }

  // Muon combination
  float sum_pT_mumu = 0.;
  bool pass_muon = false;
  int imu1 = -1;
  int imu2 = -1;
  if (muon_in.size () >= 2)
  {
    for (unsigned int i = 0; i < muon_in.size (); i++)
    {
      for (unsigned int j = i + 1; j < muon_in.size (); j++) 
      {
        if (pass_muon)  continue;
        if ((muon_in[i].charge == muon_in[j].charge))
	if( !isForMM || (isForMM && 
	((muon_in[i].RelIso03PF()<iso1_in && muon_in[j].RelIso03PF()<iso2_in) ||
	(muon_in[i].RelIso03PF()<iso2_in && muon_in[j].RelIso03PF()<iso1_in))))
        {
          pass_muon = true;
          sum_pT_mumu = muon_in[i].p4.Pt () + muon_in[j].p4.Pt ();
          imu1 = i;
          imu2 = j;
        }
      }
    }
    if(!pass_muon) // No SS pair found
    if( !isForMM || (isForMM && 
    ((muon_in[0].RelIso03PF()<iso1_in && muon_in[1].RelIso03PF()<iso2_in) ||
    (muon_in[0].RelIso03PF()<iso2_in && muon_in[1].RelIso03PF()<iso1_in))))
    { 
      sum_pT_mumu = muon_in[0].p4.Pt () + muon_in[1].p4.Pt ();
      imu1=0;
      imu2=1;
    }
  }


  // ElecMu combination
  float sum_pT_emu_start = 0.;
  float sum_pT_emu = 0.;
  bool pass_emu = false;
  int je1 = -1;
  int jmu2 = -1;
  if (muon_in.size () >= 1 && elec_in.size () >= 1) {
    for (unsigned int i = 0; i < muon_in.size (); i++) {
      for (unsigned int j = 0; j < elec_in.size (); j++) {
        if ((muon_in[i].charge == elec_in[j].charge))
	if( !isForMM || (isForMM && 
	( (muon_in[i].RelIso03PF()<iso1_in && elec_in[j].RelIso03PF()<iso2_in))))
        {
	  pass_emu = true;
          sum_pT_emu = muon_in[i].p4.Pt () + elec_in[j].p4.Pt ();
          if (sum_pT_emu > sum_pT_emu_start) {
            sum_pT_emu_start = sum_pT_emu;
            je1 = j;
            jmu2 = i;
          }
        }
      }
    }
    if(!pass_emu) // No SS pair found
    if( !isForMM || (isForMM && 
    ( (muon_in[0].RelIso03PF()<iso1_in && elec_in[0].RelIso03PF()<iso2_in))))
    {
      sum_pT_emu = muon_in[0].p4.Pt () + elec_in[0].p4.Pt ();
      je1=0;
      jmu2=0;
    }
  }

  if(pass_elec && !pass_muon) 
  { 
    elec_out.push_back (elec_in[ie1]);
    elec_out.push_back (elec_in[ie2]);
  }
  if(!pass_elec && pass_muon) 
  { 
    muon_out.push_back (muon_in[imu1]);
    muon_out.push_back (muon_in[imu2]);
  }
  if(pass_elec && pass_muon)
  {
    float sum[2] = { sum_pT_ee, sum_pT_mumu};
    int sortedIndices[2];
    TMath::Sort (2, sum, sortedIndices);
    if (sortedIndices[0] == 0 && sum_pT_ee != 0.) {
      elec_out.push_back (elec_in[ie1]);
      elec_out.push_back (elec_in[ie2]);
    }
    else if (sortedIndices[0] == 1 && sum_pT_mumu != 0.) {
      muon_out.push_back (muon_in[imu1]);
      muon_out.push_back (muon_in[imu2]);
    }    
  }
  if(!pass_elec && !pass_muon && pass_emu)
  {
    elec_out.push_back (elec_in[je1]);
    muon_out.push_back (muon_in[jmu2]);
    
  }
  
  // No SS pair
  if(!pass_elec && !pass_muon && !pass_emu)
  {
    float sum[3] = { sum_pT_ee, sum_pT_mumu, sum_pT_emu };
    int sortedIndices[3];
    TMath::Sort (3, sum, sortedIndices);
    if (sortedIndices[0] == 0 && sum_pT_ee != 0.) {
      elec_out.push_back (elec_in[ie1]);
      elec_out.push_back (elec_in[ie2]);
    }
    else if (sortedIndices[0] == 1 && sum_pT_mumu != 0.) {
      muon_out.push_back (muon_in[imu1]);
      muon_out.push_back (muon_in[imu2]);
    }
    else if (sortedIndices[0] == 2 && sum_pT_emu != 0.) {
      elec_out.push_back (elec_in[je1]);
      muon_out.push_back (muon_in[jmu2]);
    }
  }



  if (elec_out.size () + muon_out.size () == 2) {
    if (muon_out.size () == 2) {
      if(pass_muon) CandPairType = "mumu_ss";
      else CandPairType = "mumu_os";
    }
    if (elec_out.size () == 2) {
      if(pass_elec) CandPairType = "ee_ss";
      else CandPairType = "ee_os";
    }
    if (muon_out.size () == 1 && elec_out.size () == 1) {
      if(pass_emu) CandPairType = "emu_ss";
      else CandPairType = "emu_os";
    }
  }
  else CandPairType="false";


  if (CandPairType == "ee_ss" || CandPairType == "emu_ss" || CandPairType == "mumu_ss"
       || CandPairType == "ee_os" || CandPairType == "emu_os" || CandPairType == "mumu_os")
    return true;
  else
    return false;

}



bool SSDiLeptonSelection::GetLeptonPair (std::vector < NTMuon > &muon_out, std::vector < NTElectron > &elec_out, string & CandPairType)
{
  return GetLeptonPair (GetSelectedMuons (), GetSelectedElectrons (), muon_out, elec_out, CandPairType, false);
}

bool SSDiLeptonSelection::GetLeptonPairForMM (std::vector < NTMuon > &muon_out, std::vector < NTElectron > &elec_out, string & CandPairType, float iso1_in, float iso2_in)
{
  return GetLeptonPair (GetSelectedMuonsNoIso (), GetSelectedElectronsNoIso (), muon_out, elec_out, CandPairType, true, iso1_in, iso2_in);
}


bool SSDiLeptonSelection::TestIsolationOfPair (float iso1_in, float iso2_in, std::vector < NTMuon > muon_in, std::vector < IPHCTree::NTElectron > elec_in)
{

  bool pass_elec = false;
  if (elec_in.size () == 2) {
    for (unsigned int i = 0; i < elec_in.size (); i++) {
      for (unsigned int j = i + 1; j < elec_in.size (); j++) {
	if (pass_elec)
	  continue;
	if ((elec_in[i].RelIso03PF () < iso1_in && elec_in[j].RelIso03PF () < iso2_in) || (elec_in[i].RelIso03PF () < iso2_in && elec_in[j].RelIso03PF () < iso1_in)) {
	  pass_elec = true;
	}
      }
    }
  }



  bool pass_muon = false;
  if (muon_in.size () == 2) {
    for (unsigned int i = 0; i < muon_in.size (); i++) {
      for (unsigned int j = i + 1; j < muon_in.size (); j++) {
	if (pass_muon)
	  continue;
	if ((muon_in[i].RelIso03PF () < iso1_in && muon_in[j].RelIso03PF () < iso2_in) || (muon_in[i].RelIso03PF () < iso2_in && muon_in[j].RelIso03PF () < iso1_in)) {
	  pass_muon = true;
	}
      }
    }
  }


  bool pass_emu = false;
  if ((muon_in.size () + elec_in.size ()) == 2) {
    for (unsigned int i = 0; i < muon_in.size (); i++) {
      for (unsigned int j = 0; j < elec_in.size (); j++) {
	if (pass_emu)
	  continue;
	if (muon_in[i].RelIso03PF () < iso1_in && elec_in[j].RelIso03PF () < iso2_in) {
	  pass_emu = true;
	}
      }
    }
  }

  if (pass_elec || pass_muon || pass_emu)
    return true;
  else
    return false;

}


TLorentzVector SSDiLeptonSelection::DiLeptonCand (const std::vector < NTMuon > &muons_cand, const std::vector < NTElectron > &electrons_cand)
{
  TLorentzVector DiLepton;
  if (muons_cand.size () + electrons_cand.size () != 2)
    return DiLepton;
  for (unsigned int i = 0; i < muons_cand.size (); i++)
    DiLepton += muons_cand[i].p4;
  for (unsigned int i = 0; i < electrons_cand.size (); i++)
    DiLepton += electrons_cand[i].p4;
  return DiLepton;
}

float SSDiLeptonSelection::DiLeptonMass (const std::vector < NTMuon > &muons_cand, const std::vector < NTElectron > &electrons_cand)
{
  return DiLeptonCand (muons_cand, electrons_cand).M ();
}

float SSDiLeptonSelection::DiLeptonMT (const std::vector < NTMuon > &muons_cand, const std::vector < NTElectron > &electrons_cand)
{
  return DiLeptonCand (muons_cand, electrons_cand).Mt ();
}

bool SSDiLeptonSelection::DiLeptonMassCut (float MinValue, pair < float, float >ZMassWindow, const std::vector < NTMuon > &muons_cand, const std::vector < NTElectron > &electrons_cand,
					 string channelName)
{


  string leptonpairname;
  if (channelName != string ("") && channelName != string ("*") && channelName != string ("allChannels")) {
    leptonpairname = channelName;
  }
  else {
    if (electrons_cand.size () == 2)
      leptonpairname = "ee";
    if (muons_cand.size () == 2)
      leptonpairname = "mumu";
    if (muons_cand.size () == 1 && electrons_cand.size () == 1)
      leptonpairname = "emu";
  }

  bool iresult = true;
  float mass = DiLeptonMass (muons_cand, electrons_cand);

  // reject low DY mass;
  if (mass < MinValue)
    iresult = false;
  // Z veto;
  if (leptonpairname != "emu" && ZMassWindow.first < mass && mass < ZMassWindow.second)
    iresult = false;

  return iresult;

}

bool SSDiLeptonSelection::DiLeptonMassCut (const std::vector < NTMuon > &muons_cand, const std::vector < NTElectron > &electrons_cand, string channelName)
{
  return DiLeptonMassCut (MinMassCut_, ZMassWindow_, muons_cand, electrons_cand, channelName);
}


int SSDiLeptonSelection::doFullSelection (Dataset * dataset, string channelName, bool print, 
          bool isForMM, float iso1_in, float iso2_in, 
	  bool applyJES, float JESParam, bool applyEES, float EESParam, bool applyMES, float MESParam, bool applyJER, float JERFactor, bool applyMETS, float METScale)
{
 
  //clear object collections
  jetsAna.clear();
  bjetsAna.clear();
  electronsAna.clear();
  muonsAna.clear();

  bool applyChannel = false;
  if (channelName != string ("") && channelName != string ("*") && channelName != string ("allChannels"))
    applyChannel = true;
  //boolean for the selection step: true = pass the cut
  bool step_trigger = false;
  bool step_pair = false;
  bool step_Zveto = false;
  bool step_jets = false;
  bool step_met = false;
  bool step_bjets = false;

  TString dump;
  ostringstream runNumber_oss;
  ostringstream eventNumber_oss;
  runNumber_oss << getRunNumber();
  eventNumber_oss << getEventNumber();
//  dump = channel_.ChannelName () + string (" | ") + runNumber_oss.str () + string (" | ") + eventNumber_oss.str () + string (" | ");
  //dump+=string("")+runNumberS;  

  //double METEMu = METCuts_.first;
  double METLL = METCuts_.second;
  std::vector < NTMuon > muon_cand;
  std::vector < NTElectron > elec_cand;
  pairType_ = "";
  dimass_=0.;

  int FinalStep = 0;
  //Step 1        MC-match
    FinalStep++;
    //Step 2        Trigger
    if (passTriggerSelection (dataset, channelName)) {
      step_trigger = true;
    //Step 3        Dilepton pair choice
    // Rajouter selections avec iso
    if ((!isForMM && (GetLeptonPair (GetSelectedMuons(applyMES, MESParam), GetSelectedElectrons(applyEES, EESParam), muon_cand, elec_cand, pairType_) == true) 
	&& (!applyChannel || (applyChannel && pairType_ == channelName))) ||
	(isForMM && (GetLeptonPair(GetSelectedMuonsNoIso(), GetSelectedElectronsNoIso(), muon_cand, elec_cand, pairType_, true, 10000, 10000)) == true 
	&& (!applyChannel || (applyChannel && pairType_ == channelName)) && TestIsolationOfPair(iso1_in, iso2_in, muon_cand, elec_cand))
	) {
      step_pair = true;
      muonsAna = muon_cand;
      electronsAna = elec_cand;

      float lep1PtxCharge = 0;
      float lep2PtxCharge = 0;
      float lep1RelIso = 0;
      float lep2RelIso = 0;
      ostringstream lep1PtxCharge_oss;
      ostringstream lep2PtxCharge_oss;
      ostringstream lep1RelIso_oss;
      ostringstream lep2RelIso_oss;
      //if (pairType_ == "mumu") {
      if (pairType_.find("mumu")!=string::npos) {
	lep1PtxCharge = muon_cand[0].p4.Pt () * muon_cand[0].charge;
	lep2PtxCharge = muon_cand[1].p4.Pt () * muon_cand[1].charge;
	lep1RelIso = muon_cand[0].RelIso03PF ();
	lep2RelIso = muon_cand[1].RelIso03PF ();
      }
      //if (pairType_ == "ee") {
      if (pairType_.find("ee")!=string::npos) {
	lep1PtxCharge = elec_cand[0].p4.Pt () * elec_cand[0].charge;
	lep2PtxCharge = elec_cand[1].p4.Pt () * elec_cand[1].charge;
	lep1RelIso = elec_cand[0].RelIso03PF ();
	lep2RelIso = elec_cand[1].RelIso03PF ();
      }
      //if (pairType_ == "emu") {
      if (pairType_.find("emu")!=string::npos) {
	if (elec_cand[0].p4.Pt () > muon_cand[0].p4.Pt ()) {

	  lep1PtxCharge = elec_cand[0].p4.Pt () * elec_cand[0].charge;
	  lep2PtxCharge = muon_cand[0].p4.Pt () * muon_cand[0].charge;
	  lep1RelIso = elec_cand[0].RelIso03PF ();
	  lep2RelIso = muon_cand[0].RelIso03PF ();
	}
	else {
	  lep2PtxCharge = elec_cand[0].p4.Pt () * elec_cand[0].charge;
	  lep1PtxCharge = muon_cand[0].p4.Pt () * muon_cand[0].charge;
	  lep2RelIso = elec_cand[0].RelIso03PF ();
	  lep1RelIso = muon_cand[0].RelIso03PF ();
	}
      }
      dimass_ = DiLeptonMass (muon_cand, elec_cand);
      ostringstream dimass_oss;
      dimass_oss << dimass_;
      lep1PtxCharge_oss << lep1PtxCharge;
      lep2PtxCharge_oss << lep2PtxCharge;
      lep1RelIso_oss << lep1RelIso;
      lep2RelIso_oss << lep2RelIso;
      dump += lep1PtxCharge_oss.str () + "," + lep2PtxCharge_oss.str () + " | " + lep1RelIso_oss.str () + "," + lep2RelIso_oss.str () + " | " + dimass_oss.str () + " | ";
      //Step 4     Z mass veto 
      if (DiLeptonMassCut (muon_cand, elec_cand, "")) {
	step_Zveto = true;
      }
      //Step 5    Minimal jet multiplicity 
      vector < NTJet > SelectedJets = GetSelectedJets (muon_cand, elec_cand, applyJES, JESParam, applyJER, JERFactor);
      jetsAna = SelectedJets;

      if (SelectedJets.size () >= 2) {
	step_jets = true;
      }
      //Step 6  MET cuts
      if (GetSelectedMET(applyJES, JESParam, applyMETS, METScale).met() < METLL) {
	  step_met = true;
      }
      //Step 7 btagging
      vector < NTJet > btagjets;
      vector < float >btagDiscri;
      ostringstream jet1pt_oss;
      ostringstream jet2pt_oss;
      ostringstream jet1bdisc_oss;
      ostringstream jet2bdisc_oss;
      ostringstream met_oss;
      if(step_jets){
      	jet1pt_oss << SelectedJets[0].p4.Pt ();
      	jet2pt_oss << SelectedJets[1].p4.Pt ();
      	jet1bdisc_oss << SelectedJets[0].bTag["trackCountingHighEffBJetTags"];
      	jet2bdisc_oss << SelectedJets[1].bTag["trackCountingHighEffBJetTags"];
     	 met_oss << GetMET().met();
      	dump += jet1pt_oss.str () + " , " + jet2pt_oss.str () + " | " + met_oss.str () + " | " + jet1bdisc_oss.str () + " , " + jet2bdisc_oss.str () + " | " + pairType_;
      }

      for (unsigned int j = 0; j < SelectedJets.size (); j++) {
	switch (btagAlgo_) {
	case 0:
	  if (SelectedJets[j].bTag["trackCountingHighEffBJetTags"] >= btagDiscriCut_) {
	    btagjets.push_back (SelectedJets[j]);
	    btagDiscri.push_back (SelectedJets[j].bTag["trackCountingHighEffBJetTags"]);
	  }
	  break;
	case 1:
	  if (SelectedJets[j].bTag["combinedSecondaryVertexBJetTags"] >= btagDiscriCut_) {
	    btagjets.push_back (SelectedJets[j]);
	    btagDiscri.push_back (SelectedJets[j].bTag["combinedSecondaryVertexBJetTags"]);
	  }
	  break;
	case 2:
	  if (SelectedJets[j].bTag["jetProbabilityBJetTags"] >= btagDiscriCut_) {
	    btagjets.push_back (SelectedJets[j]);
	    btagDiscri.push_back (SelectedJets[j].bTag["jetProbabilityBJetTags"]);
	  }
	  break;
	default:
	  cerr << "btagAlgo doesn't exist !" << endl;
	  break;
	}
      }
      bjetsAna = btagjets;
      if ((int) btagjets.size () == NofBtagJets_) {
	//FinalStep++;
	step_bjets = true;
      }
    }
  }

  //Compute a collection of jets and b-jets and a weight even if the selection of pair fails
  if(!step_pair){
      vector < NTJet > SelectedJets = GetSelectedJets (GetSelectedMuons(), GetSelectedElectrons(), applyJES, JESParam);
      jetsAna = SelectedJets;
      vector < NTJet > btagjets;
      vector < float >btagDiscri;
      for (unsigned int j = 0; j < SelectedJets.size (); j++) {
	switch (btagAlgo_) {
	case 0:
	  if (SelectedJets[j].bTag["trackCountingHighEffBJetTags"] >= btagDiscriCut_) {
	    btagjets.push_back (SelectedJets[j]);
	    btagDiscri.push_back (SelectedJets[j].bTag["trackCountingHighEffBJetTags"]);
	  }
	  break;
	case 1:
	  if (SelectedJets[j].bTag["combinedSecondaryVertexBJetTags"] >= btagDiscriCut_) {
	    btagjets.push_back (SelectedJets[j]);
	    btagDiscri.push_back (SelectedJets[j].bTag["combinedSecondaryVertexBJetTags"]);
	  }
	  break;
	case 2:
	  if (SelectedJets[j].bTag["jetProbabilityBJetTags"] >= btagDiscriCut_) {
	    btagjets.push_back (SelectedJets[j]);
	    btagDiscri.push_back (SelectedJets[j].bTag["jetProbabilityBJetTags"]);
	  }
	  break;
	default:
	  cerr << "btagAlgo doesn't exist !" << endl;
	  break;
	}
      }
      bjetsAna = btagjets;
  }


  //compute FinalStep
  if (step_trigger) {
    FinalStep++;
    if (step_pair) {
      FinalStep++;
      if (step_Zveto) {
	FinalStep++;
	if (step_jets) {
	  FinalStep++;
	  if (step_met) {
	    FinalStep++;
	    if (step_bjets) {
	      FinalStep++;
	      //the event is selected
	      if (print)
		cout << dump << endl;
	    }
	  }
	}
      }
    }
  }
  selCode_ = step_trigger*maskTriggerCut + step_pair*maskPairCut + 
	step_Zveto * maskZvetoCut + step_jets * maskJetCut + 
	step_met * maskMETCut + step_bjets * maskBjetsCut;

  return FinalStep;
}



int SSDiLeptonSelection::FillTable (SelectionTable & selTable,
                                    Dataset * dataset,
                                    int idataset, float weight)
{

  int sel = doFullSelection (dataset, selTable.Channel (), false);	// true-> has to be modified !!
  for (unsigned int i = 0; i < cuts_.size () + 1; i++)
    if (sel >= (int) i)
      selTable.Fill (idataset, i, weight);
  return sel;
}




//bool SSDiLeptonSelection::passTriggerSelection(string datasetName, string channelName){
bool SSDiLeptonSelection::passTriggerSelection (Dataset * dataset, string channelName)
{

  //bool match_HLT_Ele10_LW_L1R_recoEl = eleHLTMatch;
//cout << "start TS\n";
  ////cout << eleHLTMatch << endl;
  bool passEl = false;
  bool passMu = false;
  bool passElMu = false;
  int skim = -1;

  string datasetName = dataset->Name ();
//   if (datasetName == "MuData")
//     skim = 0;
//   if (datasetName == "EGData")
//     skim = 1;
// to be compatible with MyCutFlow and others
  if (datasetName.compare(0,6,"DataMu")==0)
    skim = 0;
  if (datasetName.compare(0,6,"DataEG")==0)
    skim = 1;
//  if (datasetName == "DataMuEG" || datasetName == "DataEGMu")
  if ( datasetName.compare(0,8,"DataMuEG")==0 || datasetName.compare(0,8,"DataEGMu")==0 )
    skim = 2;
  
  
  //cout << " datasetName " << datasetName << endl;
  
  
  //  if(runNumber == 1){
  if (!dataset->isData ()) 	//MC
  {
    //GetPointer2Trigger()->Dump();
    passMu = (GetPointer2Trigger()->IsFired("HLT_DoubleMu6_v1") // Summer11
               || GetPointer2Trigger()->IsFired("HLT_DoubleMu6_v8")); // Fall11
    passEl = (GetPointer2Trigger()->IsFired("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2") // Summer11
               || GetPointer2Trigger()->IsFired("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8")); // Fall11
    passElMu = (GetPointer2Trigger()->IsFired("HLT_Mu8_Ele17_CaloIdL_v2") || GetPointer2Trigger()->IsFired("HLT_Mu10_Ele10_CaloIdL_v3") // Summer11
               || GetPointer2Trigger()->IsFired("HLT_Mu8_Ele17_CaloIdL_v9")); // Fall11
  }
  else
  {
    // DATA --> Taken from TopDileptonRefAnalysis2010Pass6
    // WARNING : COULD SOMEONE CHECK THE < AND <= FOR THE RUN NUMBER????

    // ERIC TO DO
    /*
      if ( 160329 <= getRunNumber() && getRunNumber() < 164237 && 
      GetPointer2Trigger()->IsFired("HLT_DoubleMu7");
      else if ( 165085 <= getRunNumber() && 
      GetPointer2Trigger()->IsFired("HLT_Mu13_Mu8_vleMu7");
      ( triggerList[i].first.compare(0,14,"HLT_Mu13_Mu8_v")==0  && triggerList[i].second == true)  )
      passMu = true;

      if (( triggerList[i].first.compare(0,50,"HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL")==0  && triggerList[i].second == true)  )
      passEl = true;

      if (160329 <= getRunNumber()  &&
      ( ( triggerList[i].first.compare(0,21,"HLT_Mu17_Ele8_CaloIdL")==0  && triggerList[i].second == true)
      || ( triggerList[i].first.compare(0,21,"HLT_Mu8_Ele17_CaloIdL")==0  && triggerList[i].second == true) ) )
      passElMu = true;
    */
	
  }

  //if (channelName == string ("ee")) 
  if (channelName.find("ee")!=string::npos) 
  {
    if (passEl)
      return true;
    else
      return false;
  }
//  if (channelName == string ("mumu")) 
  if (channelName.find("mumu")!=string::npos) 
  {
    if (passMu)
      return true;
    else
      return false;
  }
//  if (channelName == string ("emu")) 
  if (channelName.find("emu")!=string::npos) 
  {
    bool thereturn = false;
    //     if (  (passEl && skim==1 ) || ( passMu && skim==0 && !passEl ) ) thereturn = true;
    //     if ( skim == -1 && (passEl  || passMu ) ) thereturn = true;
    if ( skim == -1 &&  passElMu) return true;
    if ( skim == 2 && (passElMu) )
      thereturn = true;
    return thereturn;
  }
  if (channelName == string ("") || channelName == string ("*") || channelName == string ("allChannels")) {
    if (passEl || passMu)
      return true;
    return false;
  }

  return false;

}




float SSDiLeptonSelection::GetMinValueMassCut ()
{
  return MinMassCut_;
}

pair < float, float >SSDiLeptonSelection::GetMETCut ()
{
  return METCuts_;
}

pair < float, float >SSDiLeptonSelection::GetZmassWindowCut ()
{
  return ZMassWindow_;
}

     int SSDiLeptonSelection::GetbtagAlgo () const
     {
       return btagAlgo_;
     }

     float SSDiLeptonSelection::GetbtagDiscriCut () const
     {
       return btagDiscriCut_;
     }

     int SSDiLeptonSelection::GetNofBtagJetsCut () const
     {
       return NofBtagJets_;
     }


bool SSDiLeptonSelection::passBtagSelection(const NTJet & jet) const
{
  return (getBtagDiscr(jet) >= btagDiscriCut_);
}

double SSDiLeptonSelection::getBtagDiscr(const NTJet & jet) const
{
  switch (btagAlgo_) {
  case 0:
    return jet.bTag["trackCountingHighEffBJetTags"];
    break;
  case 1:
    return jet.bTag["combinedSecondaryVertexBJetTags"];
    break;
  case 2:
    return jet.bTag["jetProbabilityBJetTags"];
    break;
  default:
    cerr << "btagAlgo doesn't exist !" << endl;
    return -10000000;
  }
}



double SSDiLeptonSelection::getLeptonScaleFactor(double pt1, double eta1, double pt2, double eta2, string channel){
  double the_getScaleFactor = 0;
  
  if(pt1 > 100) pt1 = 99;
  if(fabs(eta1) > 2.5) eta1 = 2.4;
  
  if(pt2 > 100) pt2 = 99;
  if(fabs(eta2) > 2.5) eta2 = 2.4;
  //if(channel == "ee"){
  if(channel.find("ee")!=string::npos){
     int binx1 = getScaleFactEl()->GetXaxis()->FindBin( pt1 );
     int biny1 = getScaleFactEl()->GetYaxis()->FindBin( fabs(eta1) );
     
     int binx2 = getScaleFactEl()->GetXaxis()->FindBin( pt2 );
     int biny2 = getScaleFactEl()->GetYaxis()->FindBin( fabs(eta2) );
     
     the_getScaleFactor = getScaleFactEl()->GetBinContent( binx1, biny1 )*getScaleFactEl()->GetBinContent( binx2, biny2 );

  }
  
  //if(channel == "mumu"){
  if(channel.find("mumu")!=string::npos){
     int binx1 = getScaleFactMu()->GetXaxis()->FindBin( pt1 );
     int biny1 = getScaleFactMu()->GetYaxis()->FindBin( fabs(eta1) );
     int binx2 = getScaleFactMu()->GetXaxis()->FindBin( pt2 );
     int biny2 = getScaleFactMu()->GetYaxis()->FindBin( fabs(eta2) );
     
     the_getScaleFactor = getScaleFactMu()->GetBinContent( binx1, biny1 )*getScaleFactMu()->GetBinContent( binx2, biny2 );
  }

  //if(channel == "emu"){
  if(channel.find("emu")!=string::npos){
  
     int binx1 = getScaleFactEl()->GetXaxis()->FindBin( pt1 );
     int biny1 = getScaleFactEl()->GetYaxis()->FindBin( fabs(eta1) );
     
     int binx2 = getScaleFactMu()->GetXaxis()->FindBin( pt2 );
     int biny2 = getScaleFactMu()->GetYaxis()->FindBin( fabs(eta2) );
     
     the_getScaleFactor = getScaleFactEl()->GetBinContent( binx1, biny1 )*getScaleFactMu()->GetBinContent( binx2, biny2 );
  }
  
  
  return the_getScaleFactor;
}





double SSDiLeptonSelection::getLeptonScaleFactorError(double pt1, double eta1, double pt2, double eta2, string channel){
  double the_getScaleFactor = 0;
  
  if(pt1 > 100) pt1 = 99;
  if(fabs(eta1) > 2.5) eta1 = 2.4;
  
  if(pt2 > 100) pt2 = 99;
  if(fabs(eta2) > 2.5) eta2 = 2.4;
  
  //if(channel == "ee"){
  if(channel.find("ee")!=string::npos){
     int binx1 = getScaleFactEl()->GetXaxis()->FindBin( pt1 );
     int biny1 = getScaleFactEl()->GetYaxis()->FindBin( fabs(eta1) );
     
     int binx2 = getScaleFactEl()->GetXaxis()->FindBin( pt2 );
     int biny2 = getScaleFactEl()->GetYaxis()->FindBin( fabs(eta2) );
     
     the_getScaleFactor = getScaleFactEl()->GetBinError( binx1, biny1 )*getScaleFactEl()->GetBinContent( binx2, biny2 );

  }
  
//  if(channel == "mumu"){
  if(channel.find("mumu")!=string::npos){
  
     int binx1 = getScaleFactMu()->GetXaxis()->FindBin( pt1 );
     int biny1 = getScaleFactMu()->GetYaxis()->FindBin( fabs(eta1) );
     
     int binx2 = getScaleFactMu()->GetXaxis()->FindBin( pt2 );
     int biny2 = getScaleFactMu()->GetYaxis()->FindBin( fabs(eta2) );
     
     the_getScaleFactor = getScaleFactMu()->GetBinError( binx1, biny1 )*getScaleFactMu()->GetBinContent( binx2, biny2 );
  }

  //if(channel == "emu"){
  if(channel.find("emu")!=string::npos){
  
     int binx1 = getScaleFactEl()->GetXaxis()->FindBin( pt1 );
     int biny1 = getScaleFactEl()->GetYaxis()->FindBin( fabs(eta1) );
     
     int binx2 = getScaleFactMu()->GetXaxis()->FindBin( pt2 );
     int biny2 = getScaleFactMu()->GetYaxis()->FindBin( fabs(eta2) );
     
     the_getScaleFactor = getScaleFactEl()->GetBinError( binx1, biny1 )*getScaleFactMu()->GetBinContent( binx2, biny2 );
  }
  
  
  return the_getScaleFactor;
}












