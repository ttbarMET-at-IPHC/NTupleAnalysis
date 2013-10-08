#include "Selection/interface/TTbarMetSelection.h"
using namespace std;

#define DEBUG_MSG std::cout << "DEBUG (" << __FILE__ << ", l." << __LINE__ << ") "

// ----------------------------------------------------------------------------
// Default constructor
// ----------------------------------------------------------------------------
TTbarMetSelection::TTbarMetSelection ()
{
  btagAlgo_       = -1;
  btagDiscriCut_  = -999.;
  NofBtagJets_    = 0;

  //Fill the table of cuts
  cuts_.push_back ("All");
  cuts_.push_back ("Trigger");
  cuts_.push_back ("1 Lepton");
  cuts_.push_back ("4 jets");
  cuts_.push_back ("1 btag");
  cuts_.push_back ("MET");
  cuts_.push_back ("isolated track veto");
  cuts_.push_back ("tau veto");

  //Fill Channels
  channels_.push_back (string ("e"));
  channels_.push_back (string ("mu"));
//  channels_.push_back (string ("lepton"));
}


// ----------------------------------------------------------------------------
// Copy constructor
// ----------------------------------------------------------------------------
TTbarMetSelection::TTbarMetSelection(const TTbarMetSelection & s):
                                                                 Selection (s)
{
  cuts_          = s.cuts_;
  channels_      = s.channels_;
  btagAlgo_      = s.btagAlgo_;
  btagDiscriCut_ = s.btagDiscriCut_;
  NofBtagJets_   = s.NofBtagJets_;
}


// ----------------------------------------------------------------------------
// SetParameters
// ----------------------------------------------------------------------------
void TTbarMetSelection::SetParameters(  int btagAlgo,
                                        float btagDiscriCut,
                                        int NofBtagJets)
{
  btagAlgo_      = btagAlgo;
  btagDiscriCut_ = btagDiscriCut;
  NofBtagJets_   = NofBtagJets;
}



int TTbarMetSelection::doFullSelection (Dataset * dataset, string channelName, int* triggerME, 
      bool applyJES, float JESParam, bool applyEES, float EESParam, bool applyMES, float MESParam, bool applyJER, float JERFactor, bool applyMETS, float METScale)
{
 
  //clear object collections
  jetsAna.clear();
  bjetsAna.clear();
  electronsAna.clear();
  muonsAna.clear();

  //boolean for the selection step: true = pass the cut
  bool step_trigger    = false;
  bool step_trigger_e  = false;
  bool step_trigger_mu = false;
  bool step_1lepton    = false;
  bool step_jets       = false;
  bool step_bjets1     = false;
  bool step_met        = false;
  bool step_vetoTrack  = false;
  bool step_vetoTau    = false;

  //Collection of selected objects

  std::vector<IPHCTree::NTMuon>     goodMuons     = GetSUSYstopGoodMuons();
  std::vector<IPHCTree::NTElectron> goodElectrons = GetSUSYstopGoodElectrons(goodMuons);

  muonsAna     = GetSUSYstopSelectedMuons();
  electronsAna = GetSUSYstopSelectedElectrons(goodMuons);


  NTMET themet = GetSUSYstopType1PhiMET(dataset->isData(),goodMuons,goodElectrons);
  jetsAna      = GetSUSYstopSelectedJets(dataset->isData(),goodMuons,goodElectrons);


  int FinalStep = 0;

  // ##########################
  // #  Step 1        Trigger #
  // # + check rho value < 40 #
  // ##########################
  
  step_trigger_e = passTriggerSelection(dataset, string("e") );
  step_trigger_mu = passTriggerSelection(dataset, string("mu") );
 

  bool step_rho = false;
  if ((getRho() >= 0.0) && (getRho() < 40.0)) step_rho = true;

  // Check if step trigger is OK
  if ((step_trigger_e || step_trigger_mu) && (step_rho)) step_trigger = true;
  if (step_trigger) {

      // ##########################
      // #  Step 2   Lepton (iso) #
      // ##########################
  
     if (electronsAna.size()==1 && step_trigger_e && muonsAna.size()==0)
     {
         step_1lepton = true;
         LeptonType="e";
     }
     else if (muonsAna.size()==1 && step_trigger_mu && electronsAna.size()==0) 
     {
         step_1lepton = true;
         LeptonType="mu";
     }

     if (step_1lepton) 
     {

       TLorentzVector lepton_p;
       float lepton_charge;
       if (muonsAna.size()==1) 
       {
           lepton_p      = muonsAna[0].p4;
           lepton_charge = muonsAna[0].charge;
       }
       else
       {
           lepton_p      = electronsAna[0].p4;
           lepton_charge = electronsAna[0].charge;
       }

       // ##########################
       // #  Step 3       4 jets   #
       // ##########################

       //!! Switch to 3 jets to be able to run selection with 3 jets only in the plotter

       if (jetsAna.size () >= 3) 
       {
           step_jets = true;
           FillKinematicP4(muonsAna, electronsAna, themet, jetsAna);
    
       // ##########################
       // #  Step 4        btag    #
       // ##########################

       //cas CSV ou JP : recommandation de BTV  
       // CSVL wp = 0.244 ; CSVM wp = 0.679 ; CSVT wp = 0.898 combinedSecondaryVertexBJetTags 
       // JPL  wp = 0.275 ; JPM  wp = 0.545 ; JPT  wp = 0.790 jetProbabilityBJetTags
        
       for (unsigned int j = 0; j < jetsAna.size (); j++)
       {
           // Get discriminant
           float discr = jetsAna[j].bTag["combinedSecondaryVertexBJetTags"];
       
           // Apply CSV medium working point
           if (discr >= 0.679) bjetsAna.push_back (jetsAna[j]);
       }

       if (bjetsAna.size() >= 1)
       {
           step_bjets1 = true;

           // ##########################
           // #  Step 5        MET     #
           // ##########################

           if (themet.met() > 50) 
           {
               step_met = true;

           // ##################################
           // #  Step 6  Isolated track veto   #
           // ##################################
           
            if (GetSUSYstopIsolatedTrackVeto(lepton_p,lepton_charge))
            {
              step_vetoTrack = true;
                   
              if (GetSUSYstopTauVeto(lepton_p,lepton_charge))
              {
                 step_vetoTau = true;
              } // step 7
            } // step 6
          } // step 5
        } // step 4
      } // step 3
    } // step 2
  } // step 1

  //compute FinalStep
  if (step_trigger)   { FinalStep++;
  if (step_1lepton)   { FinalStep++;
  if (step_jets)      { FinalStep++;
  if (step_bjets1)    { FinalStep++;
  if (step_met)       { FinalStep++;
  if (step_vetoTrack) { FinalStep++;
  if (step_vetoTau)   { FinalStep++;
    } } } } } } }

  if (triggerME != 0)
  {
      (*triggerME) = 0;
      if (step_trigger_e)  (*triggerME) += 1;
      if (step_trigger_mu) (*triggerME) += 10;
  }

  return FinalStep;
}

bool TTbarMetSelection::GetSUSYstopIsolatedTrackVeto(TLorentzVector lepton_p, float lepton_charge) const
{
    // Input vector
    std::vector<IPHCTree::NTPFCandidate> vetotracks = GetPFCandidates();

    // Loop over pfcandidates
    for(unsigned int i=0 ; i < vetotracks.size() ; i++)
    {
        bool passCuts = false;

        float pfCandId = vetotracks[i].others["id"];

        if (abs(pfCandId) == 11)
        {
            if ((vetotracks[i].others["gsfPt"] > 5)
            && (fabs(vetotracks[i].others["gsfdz"]) < 0.05)
            && (vetotracks[i].trackIso / vetotracks[i].p4.Pt() < 0.2))    
            passCuts = true;
        }
        else if (abs(pfCandId == 13))
        {
            if ((vetotracks[i].p4.Pt() > 5)
            && (fabs(vetotracks[i].dz_firstGoodVertex) < 0.05)
            && (vetotracks[i].trackIso / vetotracks[i].p4.Pt() < 0.2))    
            passCuts = true;
        }
        else
        {
            if ((vetotracks[i].p4.Pt() > 10)
            && (fabs(vetotracks[i].dz_firstGoodVertex) < 0.05)
            && (vetotracks[i].trackIso / vetotracks[i].p4.Pt() < 0.1)
            && (lepton_charge != vetotracks[i].others["charge_fromID"]))
            passCuts = true;
        }

        if (passCuts == false) continue;

        // Check pfcandidate doesnt match the selected lepton
        // + apply opposite charge req.
        TLorentzVector vetoTrack_p = vetotracks[i].p4;
        if (lepton_p.DeltaR(vetoTrack_p) < 0.1) continue;

        return false;
    }

    return true;
}

bool TTbarMetSelection::GetSUSYstopTauVeto(TLorentzVector lepton_p, float lepton_charge) const
{
    std::vector<IPHCTree::NTTau> localTaus = (*GetPointer2Taus());
    for (unsigned int i = 0 ; i < localTaus.size() ; i++)
    {
        // Reject tau candidates with pT < 20 GeV
        if (localTaus[i].p4.Pt() < 20) continue;
        // Reject tau candidates with same charge than selected lepton
        if (localTaus[i].charge == lepton_charge) continue;
        // Reject tau candidates not satisfying IDs
        if (localTaus[i].ID["decayModeFinding"] != 1.0) continue;
        // Reject tau candidates too close from selected lepton
        if (localTaus[i].p4.DeltaR(lepton_p) < 0.4) continue;
        // Apply ID
        if (localTaus[i].ID["byMediumIsolationMVA2"] != 1.0) continue;
      
        return false; 
    }

    return true;
}


void TTbarMetSelection::FillKinematicP4(const std::vector < IPHCTree::NTMuon >& muon_kin,
                      const std::vector < IPHCTree::NTElectron >& elec_kin,
                      const IPHCTree::NTMET& met_kin,
                      const std::vector < IPHCTree::NTJet >& jet_kin) {

      // 0. Reset
      TLorentzVector resetvector(0.,0.,0.,0.);
      the_lepton=resetvector;
      all_hadronic=resetvector;
      top_hadronic=resetvector;
      the_met=resetvector;
      the_lepton=resetvector;
      the_4thjet=resetvector;
      w_leptonic=resetvector;
      top_leptonic=resetvector;
      

      // 1. The Lepton
           if (muon_kin.size()==1) the_lepton=muon_kin[0].p4;
      else if (elec_kin.size()==1) the_lepton=elec_kin[0].p4;
      
      // 2. The leptonic W
      TLorentzVector met_initial(met_kin.p2.Px(),met_kin.p2.Py(),0.,met_kin.p2.Mod());
      the_met=met_initial;
      w_leptonic = the_lepton + the_met;

      // 3. Total hadronic activity
      for (UInt_t ind=0;ind<jet_kin.size();ind++) all_hadronic+=jet_kin[ind].p4;

      // 4. The hadronic Top (from M3)
      // 5. The 4th jet (after M3)

      float test_pt=0;
      float test_pt_2=0;

      for (UInt_t ind=0;       ind <jet_kin.size(); ind++ )
      for (UInt_t ind1=ind+1;  ind1<jet_kin.size(); ind1++)
      for (UInt_t ind2=ind1+1; ind2<jet_kin.size(); ind2++) 
      {
        TLorentzVector combi_test;
        combi_test =  jet_kin[ind].p4;
        combi_test += jet_kin[ind1].p4;
        combi_test += jet_kin[ind2].p4;
        if (combi_test.Pt()>test_pt) 
        {
          test_pt=combi_test.Pt();
          top_hadronic=combi_test;
          for (UInt_t ind3=0; ind3<jet_kin.size();ind3++) 
          {
            if (ind3!=ind && ind3!=ind1 && ind3!=ind2 && jet_kin[ind3].p4.Pt()>test_pt_2) 
            {
               test_pt_2=jet_kin[ind3].p4.Pt();
               the_4thjet=jet_kin[ind3].p4;
            }
          }
        }
      }
           
      //5. The leptonic Top
      top_leptonic = w_leptonic + the_4thjet;


}



int TTbarMetSelection::FillTable (SelectionTable & selTable,
                                    Dataset * dataset,
                                    int idataset, float weight)
{

  int sel = doFullSelection (dataset, selTable.Channel ());    // true-> has to be modified !!
  for (unsigned int i = 0; i < cuts_.size () + 1; i++)
    if (sel >= (int) i)
      selTable.Fill (idataset, i, weight);
  return sel;
}

bool TTbarMetSelection::passTriggerSelection (Dataset * dataset, string channelName)
{

  bool passEl = false;
  bool passMu = false;

  string datasetName = dataset->Name ();
 
    // ######################################
    // #   According to CMS AN-2012/379     #
    // ######################################

            // #########################
            // #        Muon           #
            // #########################

            std::vector<IPHCTree::NTTriggerPathType> myPaths0;
            GetPointer2Trigger()->GetSubTable("HLT_IsoMu24_v*",myPaths0);
            for (unsigned int i=0;i<myPaths0.size();i++) {
              if (myPaths0[i].fired==1) {
                passMu=true;
                if (myPaths0[i].prescale>1) cout << " warning TRIGGER " << myPaths0[i].name 
                                                 << " is PRESCALED with a factor " << myPaths0[i].prescale 
                                                 << endl;
              }
            }

            std::vector<IPHCTree::NTTriggerPathType> myPaths1;
            GetPointer2Trigger()->GetSubTable("HLT_IsoMu24_eta2p1_v*",myPaths1);
            for (unsigned int i=0;i<myPaths1.size();i++) {
              if (myPaths1[i].fired==1) {
                passMu=true;
                if (myPaths1[i].prescale>1) cout << " warning TRIGGER " << myPaths1[i].name 
                                                 << " is PRESCALED with a factor " << myPaths1[i].prescale 
                                                 << endl;
              }
            }

            // #########################
            // #      Electron         #
            // #########################

            std::vector<IPHCTree::NTTriggerPathType> myPaths2;
            GetPointer2Trigger()->GetSubTable("HLT_Ele27_WP80_v*",myPaths2);
            for (unsigned int i=0;i<myPaths2.size();i++) {
                if (myPaths2[i].fired==1) {
                passEl=true;
                if (myPaths2[i].prescale>1) cout << " warning TRIGGER " << myPaths2[i].name 
                                                 << " is PRESCALED with a factor " << myPaths2[i].prescale 
                                                 << endl;
              }
            }
            
  if (channelName.find("e") !=string::npos) return passEl; 
  if (channelName.find("mu")!=string::npos) return passMu; 
  
  return false;

}



int TTbarMetSelection::GetbtagAlgo () const
{
  return btagAlgo_;
}

float TTbarMetSelection::GetbtagDiscriCut () const
{
  return btagDiscriCut_;
}

int TTbarMetSelection::GetNofBtagJetsCut () const
{
  return NofBtagJets_;
}

// ##########################################################
// #                   Objects selection                    #
// ##########################################################


// ----------------------------------------------------------------------------
//   Good muons selection
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTMuon> TTbarMetSelection::GetSUSYstopGoodMuons() const
{

  // Container for output 
  std::vector<IPHCTree::NTMuon> selectedMuons;

  // Get Muons
  std::vector<IPHCTree::NTMuon> localMuons;
  localMuons = *GetPointer2Muons();

  // Loop over muons
  for(unsigned int i=0;i<localMuons.size();i++)
  {

    if (!localMuons[i].isPFMuon)               continue;
    if (!localMuons[i].isGlobalMuon)           continue;
    if (localMuons[i].p4.Pt()        < 10)     continue;
    if (fabs(localMuons[i].p4.Eta()) > 2.4)    continue;

    float pfIsoCharged = localMuons[i].isolation["PF03Char"];
    float pfIsoNeutral = localMuons[i].isolation["PF03Neut"];
    float pfIsoPhoton  = localMuons[i].isolation["PF03Phot"];
    float pfIsoPU      = localMuons[i].isolation["PF03PU"];
    float relIso =  ( pfIsoCharged + max((float) 0.0,(float) (pfIsoNeutral + pfIsoPhoton- 0.5*pfIsoPU )) ) / localMuons[i].p4.Pt(); 
    if (relIso >= 0.15)                        continue;

    if (localMuons[i].Chi2 > 10)               continue;
    
    if (localMuons[i].NValidHits <= 0)         continue;
    if (localMuons[i].numMatchedStations <= 1) continue;   
    if (localMuons[i].pixelHits <= 0)          continue;

    if (localMuons[i].numTrackerLayersWithMeasurement <= 5) continue;
    if (localMuons[i].dxy_vertex >= 0.02)      continue;
    if (localMuons[i].dz_vertex >= 0.5)        continue;

    selectedMuons.push_back(localMuons[i]);
  }
  std::sort(selectedMuons.begin(),selectedMuons.end(),HighestPt());
  return selectedMuons;
}

// ----------------------------------------------------------------------------
//   Complete muons selection
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTMuon> TTbarMetSelection::GetSUSYstopSelectedMuons() const
{

  // Container for output
  std::vector<IPHCTree::NTMuon> selectedMuons;

  // Get muons
  std::vector<IPHCTree::NTMuon> muons = GetSUSYstopGoodMuons();

  // Loop over muons
  for(unsigned int i=0;i<muons.size();i++)
  {
    
    // Pt and Eta
    if (muons[i].p4.Pt() < 25)   continue;
    if (muons[i].p4.Eta() > 2.1) continue;

    // Reco - PF matching
    if (fabs(muons[i].bestMatch_pT - muons[i].p4.Pt()) >= 10) continue;

    // Absolute isolation
    float pfIsoCharged = muons[i].isolation["PF03Char"];
    float pfIsoNeutral = muons[i].isolation["PF03Neut"];
    float pfIsoPhoton  = muons[i].isolation["PF03Phot"];
    float pfIsoPU      = muons[i].isolation["PF03PU"];
    float absIso       = pfIsoCharged + max(0., pfIsoNeutral + pfIsoPhoton- 0.5*pfIsoPU); 
    if (absIso > 5) continue;

    selectedMuons.push_back(muons[i]);
  }
  
  std::sort(selectedMuons.begin(),selectedMuons.end(),HighestPt());
  return selectedMuons;
}

// ----------------------------------------------------------------------------
//   Good electrons selection
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTElectron> TTbarMetSelection::GetSUSYstopGoodElectrons(std::vector<IPHCTree::NTMuon> goodMuons) const
{

  // Output vector
  std::vector<IPHCTree::NTElectron> selectedElectrons;
  
  // Get new electrons
  std::vector<IPHCTree::NTElectron> localElectrons;
  localElectrons = *GetPointer2Electrons();
 
  // Loop over electrons
  for(unsigned int i=0;i<localElectrons.size();i++)
  {
      // Check muon overlap
      bool foundMuonOverlap = false;
      for (unsigned int j = 0 ; j < goodMuons.size() ; j++)
      {
          if (localElectrons[i].p4.DeltaR(goodMuons[j].p4) < 0.1) foundMuonOverlap = true;
      }
      if (foundMuonOverlap) continue;

      // p_T
      if (localElectrons[i].p4.Pt()        < 10)       continue;

      // eta cuts
      if (fabs(localElectrons[i].p4.Eta()) >= 2.4)     continue;
      if (  (fabs(localElectrons[i].etaSuperCluster) >= 1.4442) 
         && (fabs(localElectrons[i].etaSuperCluster) <= 1.566)) continue;

      // abs(deltaEta) and abs(deltaPhi)
      if ((localElectrons[i].isEB) && (fabs(localElectrons[i].deltaEtaSuperClusterTrackAtVtx) >= 0.004))    continue;
      if ((localElectrons[i].isEE) && (fabs(localElectrons[i].deltaEtaSuperClusterTrackAtVtx) >= 0.007))    continue;
      
      if ((localElectrons[i].isEB) && (fabs(localElectrons[i].deltaPhiSuperClusterTrackAtVtx) >= 0.06))     continue;
      if ((localElectrons[i].isEE) && (fabs(localElectrons[i].deltaPhiSuperClusterTrackAtVtx) >= 0.03))     continue;
      
      // sigmaIetaIeta
       if ((localElectrons[i].isEB) && (localElectrons[i].sigmaIetaIeta >= 0.01))    continue;
       if ((localElectrons[i].isEE) && (localElectrons[i].sigmaIetaIeta >= 0.03))    continue;

      // hadOverEM
      if ((localElectrons[i].isEB) && (localElectrons[i].hadronicOverEm >= 0.12))    continue;
      if ((localElectrons[i].isEE) && (localElectrons[i].hadronicOverEm >= 0.10))    continue;

      // dxy, dz
      if (localElectrons[i].dxy_vertex >= 0.02)           continue;
      if (localElectrons[i].dz_vertex  >= 0.1)            continue;

      // fabs(1/E-1/pin)
      float overE_m_overPin = fabs( (1.0 - localElectrons[i].eSuperClusterOverP) / localElectrons[i].EmEnergy_  );
      if (overE_m_overPin >= 0.05) continue;

      // Rel Iso
      float chargedIso  = localElectrons[i].isolation["RA4Charg"];
      float photonIso   = localElectrons[i].isolation["RA4Photo"];
      float neutralIso  = localElectrons[i].isolation["RA4Neutr"];
      float rho_relIso  = localElectrons[i].isolation["rho"];
      float Aeff_relIso = localElectrons[i].isolation["Aeff"];
      float relIso = (chargedIso + max((float) 0.0, (float) (photonIso + neutralIso - rho_relIso * Aeff_relIso)))/localElectrons[i].p4.Pt();
      
      if (localElectrons[i].isEB)
        if (relIso > 0.15) continue;

      if ((localElectrons[i].isEE) && (localElectrons[i].p4.Pt() >= 20))
        if (relIso > 0.15) continue;

      if ((localElectrons[i].isEE) && (localElectrons[i].p4.Pt() < 20))
        if (relIso > 0.10) continue;

      // Conversion rejection
      if (localElectrons[i].passConversionVeto == false) continue;
      if (localElectrons[i].missingHits > 1)            continue;

      // Add to selected electrons
      selectedElectrons.push_back(localElectrons[i]);
  }

  // Return output vector after sorting
  std::sort(selectedElectrons.begin(),selectedElectrons.end(),HighestPt());
  return selectedElectrons;
}

// ----------------------------------------------------------------------------
//   Complete electron selection
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTElectron> TTbarMetSelection::GetSUSYstopSelectedElectrons(std::vector<IPHCTree::NTMuon> goodMuons) const
{
  std::vector<IPHCTree::NTElectron> selectedElectrons;
  
  std::vector<IPHCTree::NTElectron> electrons = GetSUSYstopGoodElectrons(goodMuons);

  for(unsigned int i=0;i<electrons.size();i++)
  {
    
    // Pt and Eta
    if (electrons[i].p4.Pt() < 30) continue;
    if (fabs(electrons[i].etaSuperCluster) >= 1.4442)     continue;

    // Absolute isolation
    float chargedIso = electrons[i].isolation["RA4Charg"];
    float photonIso = electrons[i].isolation["RA4Photo"];
    float neutralIso = electrons[i].isolation["RA4Neutr"];
    float rho_relIso = electrons[i].isolation["rho"];
    float Aeff_relIso = electrons[i].isolation["Aeff"];
    float absIso = chargedIso + max((float) 0.0,(float) (photonIso + neutralIso - rho_relIso * Aeff_relIso));
    if (absIso > 5) continue;

    // E/Pin
    if (electrons[i].eSuperClusterOverP > 4) continue;
    
    // PF - Reco matching
    if (fabs(electrons[i].bestMatch_pT - electrons[i].p4.Pt()) > 10) continue;

    selectedElectrons.push_back(electrons[i]);
  }
 
  std::sort(selectedElectrons.begin(),selectedElectrons.end(),HighestPt());
  return selectedElectrons;
}

// ----------------------------------------------------------------------------
//   Jets selection
// ----------------------------------------------------------------------------
/*
void TTbarMetSelection::InitSUSYstopJEC(string tag)
{
    // Create the JetCorrectorParameter objects, the order does not matter.
    JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters((tag+"_L1FastJet_AK5PF.txt").c_str());
    JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters((tag+"_L2Relative_AK5PF.txt").c_str());
    JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters((tag+"_L3Absolute_AK5PF.txt").c_str());
    //JetCorrectorParameters *ResJetPar = new JetCorrectorParameters((tag+"_L2L3Residual_AK5PF.txt").c_str()); 

    //  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
    vector<JetCorrectorParameters>* vPar = new vector<JetCorrectorParameters>;
    vPar->push_back(*L1JetPar);
    vPar->push_back(*L2JetPar);
    vPar->push_back(*L3JetPar);
    //vPar->push_back(*ResJetPar);

    JetCorrector = new FactorizedJetCorrector(*vPar);
}

void TTbarMetSelection::CorrectSUSYstopJets(int DataType, std::vector<IPHCTree::NTJet> scaledJets) const
{
    
    for(unsigned int i=0;i<scaledJets.size();i++)
    {
        JetCorrector->setJetEta(scaledJets[i].p4.Eta());
        JetCorrector->setJetPt(scaledJets[i].p4.Pt());
        JetCorrector->setJetE(scaledJets[i].p4.E());
        JetCorrector->setJetA(scaledJets[i].others["jetArea"]);
        JetCorrector->setRho(getRho());
        JetCorrector->setNPV(GetVertex().size());
    
        vector<float> factors = JetCorrector->getSubCorrections();

        TLorentzVector p4uncorr = scaledJets[i].p4 * scaledJets[i].others["corr_Uncorr"];

        DEBUG_MSG << " i = " << i
                << " ; eta = " << p4uncorr.Eta() 
                << " ; jetA = " << scaledJets[i].others["jetArea"]
                << " ; ID = " << scaledJets[i].ID["LOOSE"]
                << " ; pTuncorr = " << p4uncorr.Pt() 
                << " ; pTold = " << scaledJets[i].p4.Pt()
                << " ; pTnew = " << (p4uncorr * factors[2]).Pt()
                << endl;

        DEBUG_MSG  << "0) L1         = " << scaledJets[i].others["corr_L1FastJet"]
                            << " L2  = " << scaledJets[i].others["corr_L2Relative"]
                            << " L3  = " << scaledJets[i].others["corr_L3Absolute"]
                                         << endl;

        DEBUG_MSG  << "1) factors[0] = " << factors[0]  
                            << " [1] = " << factors[1]
                            << " [2] = " << factors[2]
                            << " [3] = " << factors[3] << endl; */ /*
    }
}
*/

std::vector<IPHCTree::NTJet> TTbarMetSelection::GetSUSYstopSelectedJets(
            int DataType,
            const std::vector<IPHCTree::NTMuon>& muon_cand,
            const std::vector<IPHCTree::NTElectron>& elec_cand) const
{
  // Container for output
  std::vector<IPHCTree::NTJet> selectedJets;

  // Get scaled jets
  std::vector<IPHCTree::NTJet> scaledJets;
  scaledJets = *GetPointer2Jets();

  //CorrectSUSYstopJets(DataType,scaledJets);

  for(unsigned int i=0;i<scaledJets.size();i++)
  {
    
    //cout << "jet i = " << i << " ; Pt = " << scaledJets[i].p4.Pt() * scaledJets[i].others["corr_L3Absolute"] << " ; Eta = " << scaledJets[i].p4.Eta() << endl;

    // Loose ID
    if (scaledJets[i].ID["LOOSE"] != 1.) continue;
    if (scaledJets[i].ID["PU_IDTight5x"] != 1.) continue;

    // Eta and Pt cuts
         if ((DataType == 0) && (fabs(scaledJets[i].p4.Eta()) >= 2.4 || scaledJets[i].p4.Pt() * scaledJets[i].others["corr_L3Absolute"]  < 30)) 
        continue; 
    else if ((DataType == 1) && (fabs(scaledJets[i].p4.Eta()) >= 2.4 || scaledJets[i].p4.Pt() * scaledJets[i].others["corr_L2L3Residual"]  < 30)) 
        continue;
    
    //cout << "  > pass Pt/Eta cuts" << endl;

    // Overlap removal
    double deltaRmu = 10000;
    double deltaRel = 10000;
    
    for(unsigned int imu=0; imu< muon_cand.size(); imu++)
    {
      double deltaR = scaledJets[i].p4.DeltaR(muon_cand[imu].p4);
      if(deltaR < deltaRmu) deltaRmu = deltaR;
    }
    
    for(unsigned int iel=0; iel< elec_cand.size(); iel++)
    {
      double deltaR = scaledJets[i].p4.DeltaR(elec_cand[iel].p4);
      if(deltaR < deltaRel) deltaRel = deltaR;
    }
    
    if( deltaRmu > 0.4  && deltaRel > 0.4)
    {
        //cout << "  > pass deltaR" << endl;
        //cout << "  > OK !" << endl;
                         selectedJets.push_back(scaledJets[i]);
     }
  }
  std::sort(selectedJets.begin(),selectedJets.end(),HighestPt());
  return selectedJets;
}

// ----------------------------------------------------------------------------
//   Met computation
// ----------------------------------------------------------------------------

IPHCTree::NTMET TTbarMetSelection::GetSUSYstopType1MET(
            int DataType,
            const std::vector<IPHCTree::NTMuon>& muon_cand,
            const std::vector<IPHCTree::NTElectron>& elec_cand) const
{


  // Container for output
  std::vector<IPHCTree::NTJet> selectedJets;

  // Get scaled jets
  std::vector<IPHCTree::NTJet> scaledJets;
  scaledJets = *GetPointer2Jets();

  // # #           Get the jets                #

  for(unsigned int i=0;i<scaledJets.size();i++)
  {
    
    if (scaledJets[i].ID["LOOSE"] != 1.) continue; 
    
    if ((DataType == 0) && (fabs(scaledJets[i].p4.Eta()) >= 4.7 
                        || scaledJets[i].p4.Pt() * scaledJets[i].others["corr_L3Absolute"]  < 10)) 
        continue; 
    else if ((DataType == 1) && (fabs(scaledJets[i].p4.Eta()) >= 4.7 
                             || scaledJets[i].p4.Pt() * scaledJets[i].others["corr_L2L3Residual"]  < 10)) 
        continue;

    double deltaRmu = 10000;
    double deltaRel = 10000;
    
    for(unsigned int imu=0; imu< muon_cand.size(); imu++)
    {
      double deltaR = scaledJets[i].p4.DeltaR(muon_cand[imu].p4);
      if(deltaR < deltaRmu) deltaRmu = deltaR;
    }
    
    for(unsigned int iel=0; iel< elec_cand.size(); iel++)
    {
      double deltaR = scaledJets[i].p4.DeltaR(elec_cand[iel].p4);
      if(deltaR < deltaRel) deltaRel = deltaR;
    }
    
    if( deltaRmu > 0.4  && deltaRel > 0.4)
                         selectedJets.push_back(scaledJets[i]);
  }

  // # # Actually compute the type1 met from pfMet #

  NTMET rawPfMET = GetSelectedMET(false,1.0,false,1.0);

  float metx = rawPfMET.p2.Mod() * cos( rawPfMET.p2.Phi() );
  float mety = rawPfMET.p2.Mod() * sin( rawPfMET.p2.Phi() );

  for (unsigned int i = 0 ; i < selectedJets.size() ; i++)
  {
    
    float l1Corr = selectedJets[i].others["corr_L1FastJet"];
    float lastCorr = 999.0;
    if (DataType == 0)
        lastCorr = selectedJets[i].others["corr_L3Absolute"];
    else if (DataType == 1)
        lastCorr = selectedJets[i].others["corr_L2L3Residual"];
        
    metx += selectedJets[i].p4.Px() * (l1Corr - lastCorr);
    mety += selectedJets[i].p4.Py() * (l1Corr - lastCorr);

  }

  rawPfMET.p2.Set(metx,mety);
  
  return rawPfMET;
}

IPHCTree::NTMET TTbarMetSelection::GetSUSYstopType1PhiMET(
            int DataType,
            const std::vector<IPHCTree::NTMuon>& muon_cand,
            const std::vector<IPHCTree::NTElectron>& elec_cand) const
{
      NTMET the_type1met_ = GetSUSYstopType1MET(DataType,muon_cand,elec_cand);
      NTMET the_type1phimet_ = the_type1met_;

      int Nvtx = GetVertex().size();

      //DEBUG_MSG << "Nvtx = " << Nvtx << endl;;
      float metx = the_type1met_.p2.Px();
      float mety = the_type1met_.p2.Py();

      float metx_phiCorr = 0.0;
      float mety_phiCorr = 0.0;

      // MC corrections
      if (DataType == 0)
      {
          metx_phiCorr = metx - (+0.1166 + 0.0200*Nvtx);
          mety_phiCorr = mety - (+0.2764 - 0.1280*Nvtx);
      }
      // Data corrections
      else if (DataType == 1)
      {
          metx_phiCorr = metx - (+0.2661 + 0.3217*Nvtx);
          mety_phiCorr = mety - (-0.2251 - 0.1747*Nvtx);
      }

      the_type1phimet_.p2.Set(metx_phiCorr,mety_phiCorr);
      
      return the_type1phimet_;
}


