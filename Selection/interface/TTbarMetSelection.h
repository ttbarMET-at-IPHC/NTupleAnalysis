#ifndef TTbarMetSelection_h
#define TTbarMetSelection_h

// IPHC headers
#include "Selection/interface/Selection.h"
#include "Selection/interface/SelectionTable.h"
#include "Tools/interface/Dataset.h"
//#include "JEC/interface/JetCorrectorParameters.h"
//#include "JEC/interface/FactorizedJetCorrector.h"

#include "EventReco/interface/StopAnaReco.h"
#include "EventReco/interface/Resolution.h"
#include "EventReco/interface/Mt2Com_bisect.h"

// system include files
#include <memory>
#include <vector>
#include <sstream>

/**
	Steps of the selection: (integer returned by doFullSelection() or FillTable(...))
	- Step 1        Trigger
	- Step 2  	Lepton
	- Step 3 	Minimal jet multiplicity 
	- Step 4 	btagging cuts

*/

class TTbarMetSelection: public Selection
{

  // -------------------------------------------------------------
  //                       method members
  // -------------------------------------------------------------
 public:

  /*
  void CorrectSUSYstopJets(int DataType, std::vector<IPHCTree::NTJet> scaledJets) const;
  void InitSUSYstopJEC(string tag);
  */

  //! Constructor without argument
  TTbarMetSelection();

  //! Copy constructor
  TTbarMetSelection(const TTbarMetSelection &);

  //! Destructor
  ~TTbarMetSelection()
  { }

  //! Mutator for parameters
  void SetParameters(int btagAlgo,
                     float btagDiscriCut,
                     int NofBtagJets);

  //! Accessor to cut list
  std::vector<std::string> GetCutList() const {return cuts_;}

  //! Accessor to channel list
  std::vector<std::string> GetChannelList() const {return channels_;}

  //! Accessor to Channel
  signed int GetChannel(const std::string& LeptonType) const
  {
    for (unsigned int i=0; i<channels_.size(); i++)
    {
      if (channels_[i]==LeptonType) return static_cast<signed int>(i);
    }
    return -999;
  }

  //! LoadEvent 
  bool LoadEvent(const NTEvent * event)
  { return Event::LoadEvent(event); }

  //! Accessor to jetsAna
  const std::vector<IPHCTree::NTJet>& GetJetsForAna() const
  { return jetsAna; }

  //! Accessor to bjetsAna
  const std::vector<IPHCTree::NTJet>& GetBJetsForAna() const
  { return bjetsAna; }

  //! Accessor to electronsAna
  const std::vector<IPHCTree::NTElectron>& GetElectronsForAna() const
  { return electronsAna; }

  //! Accessor to muonsAna
  const std::vector<IPHCTree::NTMuon>& GetMuonsForAna() const
  { return muonsAna; }


  void FillKinematicP4(const std::vector < IPHCTree::NTMuon >& muon_kin,
                      const std::vector < IPHCTree::NTElectron >& elec_kin,
                      const IPHCTree::NTMET& met_kin,
                      const std::vector < IPHCTree::NTJet >& jet_kin);

  int FillTable(SelectionTable& selTable, Dataset* dataset, int idataset, float weight); /** Fill the selectionTable according to the result of doFullSelection  for an 
                                                                                             event of weight "weight" of a given dataset idataset - Returns the integer of doFullSelection() */


  //
  // Objects selection
  //
  
  bool GetSUSYstopIsolatedTrackVeto(TLorentzVector lepton_p, float lepton_charge) const;
  bool GetSUSYstopTauVeto(TLorentzVector lepton_p, float lepton_charge) const;

  std::vector<IPHCTree::NTMuon> GetSUSYstopGoodMuons() const;
  std::vector<IPHCTree::NTMuon> GetSUSYstopSelectedMuons() const;
  
  std::vector<IPHCTree::NTElectron> GetSUSYstopGoodElectrons(
		    std::vector<IPHCTree::NTMuon> goodMuons) const;
  std::vector<IPHCTree::NTElectron> GetSUSYstopSelectedElectrons(
		    std::vector<IPHCTree::NTMuon> goodMuons) const;
  
  std::vector<IPHCTree::NTJet> GetSUSYstopSelectedJets(
            int DataType,
			const std::vector<IPHCTree::NTMuon>& muon_cand,
            const std::vector<IPHCTree::NTElectron>& elec_cand) const;

  IPHCTree::NTMET GetSUSYstopType1MET(
    		    int DataType,
				const std::vector<IPHCTree::NTMuon>& muon_cand,
	            const std::vector<IPHCTree::NTElectron>& elec_cand) const;

  IPHCTree::NTMET GetSUSYstopType1PhiMET(
    	        int DataType,
				const std::vector<IPHCTree::NTMuon>& muon_cand,
    			const std::vector<IPHCTree::NTElectron>& elec_cand) const;

  /**
   * return a integer which correspond to the last step that the event passes in the selection 
   * - possibility to check if the candPair correspond to the correct channel 
   *  & to cut at gen level for a given channel 
   * -  compute also the weight associated to btag
   */

  int doFullSelection(Dataset* dataset, string channelName=string(""), int* triggerME = 0, 
	  bool applyJES = false, float JESParam = 1., bool applyEES = false, float EESParam = 1., bool applyMES = false, float MESParam = 1., bool applyJER = false, float JERFactor = 0., bool applyMETS = false, float METScale = 1.);


  bool passTriggerSelection(Dataset * dataset, string channelName);
  

  int GetbtagAlgo() const;
  float GetbtagDiscriCut() const;
  int GetNofBtagJetsCut() const; 


  /**
   * Returns the 32-bit selection cote. It is only filled during the doFullSelection method.
   * Each bit (from the LSB) represents the
   * outcome of a particular filter step. It can be probed with the masks,
   * e.g. GetSelCode()&SSDiLeptonSelection::maskZvetoCut tells whether the Z-veto cut 
   * was successful.
   */

  /** Invariant mass of the two leptons.
   *  It is only filled during the doFullSelection method.
   */
  float HT() const             { return all_hadronic.Pt();}
  float M3() const             { return top_hadronic.M();}
  float M_topleptonic() const  { return top_leptonic.M();}
  float MT_wleptonic() const   { return sqrt( 2.* the_lepton.Pt() * the_met.Pt() *(1. - cos( the_lepton.Phi() - the_met.Phi()) )) ;}
  float PT_tophad() const      { return top_hadronic.Pt();}
  float PT_topleptonic() const { return top_leptonic.Pt();}
  float PT_wleptonic()  const  { return w_leptonic.Pt();}
  float Dphi_lmet() const      { return the_lepton.DeltaPhi(the_met);}
  float Dphi_ljet4() const     { return the_lepton.DeltaPhi(the_4thjet);}
  float Dphi_tops() const      { return top_hadronic.DeltaPhi(top_leptonic);}
  float Deta_lth() const       { 
                                 float deta_e_tophad= the_lepton.Eta() - top_hadronic.Eta();
                                 if (deta_e_tophad<0.) deta_e_tophad*=-1.;
                                 return deta_e_tophad; }
  float Deta_ljet4() const     { 
      float deta_e_jet4= the_lepton.Eta() - the_4thjet.Eta();
      if (deta_e_jet4<0.) deta_e_jet4*=-1.;
     return deta_e_jet4; }
  float Met() const            { return the_met.Pt(); }

  int GetLeptonType() {
    if (LeptonType=="e") return 0;
    else if (LeptonType=="mu") return 1;
    else return -1;
  }
      
  float HT_ratio() const      
  { 
      float HT_onTheSideOfMET = 0;
      float HT_total = 0;
      for (unsigned int i = 0 ; i < jetsAna.size() ; i++)
      {
          if (abs(the_met.DeltaPhi(jetsAna[i].p4)) < 3.1415/2.0) 
              HT_onTheSideOfMET += jetsAna[i].p4.Pt();
            
          HT_total += jetsAna[i].p4.Pt();
      }

      return HT_onTheSideOfMET / HT_total;
  }
      
  float HadronicChi2(bool runningOnData) const
  {
    Resolution Chi2;
    float value = Chi2.GetChi2(GetJetsForAna(),runningOnData);
    return value;
  }

  float MT2W() const
  {
    Mt2Com_bisect Mt2;
    float value = Mt2.calculateMT2w(GetJetsForAna(),
                                     GetBJetsForAna(),
                                     GetMuonsForAna(),
                                     GetElectronsForAna(),
                                     the_met.Vect().XYvector(),
                                     "MT2w");
    return value;
  }
  
  float DPhi_MET_leadingJets() const
  {
    TLorentzVector firstLeadingJet;
    TLorentzVector secondLeadingJet;

    for (unsigned int j = 0 ; j < jetsAna.size() ; j++)
    {
       if (jetsAna[j].p4.Pt() > firstLeadingJet.Pt())
           firstLeadingJet = jetsAna[j].p4;

       else if (jetsAna[j].p4.Pt() > secondLeadingJet.Pt())
           secondLeadingJet = jetsAna[j].p4;
    }

    if (abs(the_met.DeltaPhi(firstLeadingJet)) < abs(the_met.DeltaPhi(secondLeadingJet)))
        return abs(the_met.DeltaPhi(firstLeadingJet));
    else
        return abs(the_met.DeltaPhi(secondLeadingJet));
  }

    
  // -------------------------------------------------------------
  //                        data members
  // -------------------------------------------------------------
 private:
 
  // parameters
  std::vector<std::string> cuts_;
  std::vector<std::string> channels_;
  int                      btagAlgo_;
  float                    btagDiscriCut_;
  int                      NofBtagJets_;


  TLorentzVector all_hadronic;
  TLorentzVector top_hadronic;
  TLorentzVector the_met;
  TLorentzVector the_lepton;
  TLorentzVector the_4thjet;  
  TLorentzVector w_leptonic;
  TLorentzVector top_leptonic;
  string LeptonType;
  
  

  //Objects for analysis
  std::vector<IPHCTree::NTElectron> electronsAna;
  std::vector<IPHCTree::NTMuon>     muonsAna;
  std::vector<IPHCTree::NTJet>      jetsAna;
  std::vector<IPHCTree::NTJet>      bjetsAna;

  //FactorizedJetCorrector* JetCorrector;
};

#endif
