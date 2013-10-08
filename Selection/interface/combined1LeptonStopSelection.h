#ifndef combined1LeptonStopSelection_h
#define combined1LeptonStopSelection_h

// IPHC headers
#include "Selection/interface/Selection.h"
#include "Selection/interface/SelectionTable.h"
#include "Tools/interface/Dataset.h"

#include "EventReco/interface/StopAnaReco.h"
#include "EventReco/interface/Resolution.h"
#include "EventReco/interface/Mt2Com_bisect.h"

// system include files
#include <memory>
#include <vector>
#include <sstream>

class combined1LeptonStopSelection: public Selection
{

  // -------------------------------------------------------------
  //                       method members
  // -------------------------------------------------------------
 public:

  //! Constructor without argument
  combined1LeptonStopSelection();

  //! Copy constructor
  combined1LeptonStopSelection(const combined1LeptonStopSelection &);

  //! Destructor
  ~combined1LeptonStopSelection()
  { }

  //! LoadEvent 
  bool LoadEvent(const NTEvent * event) { return Event::LoadEvent(event); }

  //! Accessor to jetsAna
  const std::vector<IPHCTree::NTJet>& GetJetsForAna() const { return jetsAna; }

  //! Accessor to bjetsAna
  const std::vector<IPHCTree::NTJet>& GetBJetsForAna() const { return bjetsAna; }

  void FillKinematicP4(const std::vector < IPHCTree::NTMuon >& muon_kin,
                       const std::vector < IPHCTree::NTElectron >& elec_kin,
                       const IPHCTree::NTMET& met_kin,
                       const std::vector < IPHCTree::NTJet >& jet_kin);

  bool GetSUSYstopIsolatedTrackVeto(TLorentzVector lepton_p, float lepton_charge) const;
  bool GetSUSYstopTauVeto(TLorentzVector lepton_p, float lepton_charge) const;

  std::vector<IPHCTree::NTMuon> GetSUSYstopGoodMuons()     const;
  std::vector<IPHCTree::NTMuon> GetSUSYstopSelectedMuons() const;
  
  std::vector<IPHCTree::NTElectron> GetSUSYstopGoodElectrons    (std::vector<IPHCTree::NTMuon> goodMuons) const;
  std::vector<IPHCTree::NTElectron> GetSUSYstopSelectedElectrons(std::vector<IPHCTree::NTMuon> goodMuons) const;
  
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

  bool doFullSelection(Dataset* dataset);
  bool checkPathHasBeenFired(string path);
  bool passTrigger(string channel);

  int GetbtagAlgo() const;
  float GetbtagDiscriCut() const;
  int GetNofBtagJetsCut() const; 

  float M3() const             { return top_hadronic.M();}
  float M_topleptonic() const  { return top_leptonic.M();}
  float MT_wleptonic() const   { return sqrt( 2.* the_leading_lepton.Pt() * the_met.Pt() *(1. - cos( the_leading_lepton.Phi() - the_met.Phi()) )) ;}
  float PT_tophad() const      { return top_hadronic.Pt();}
  float PT_topleptonic() const { return top_leptonic.Pt();}
  float PT_wleptonic()  const  { return w_leptonic.Pt();}
  float Dphi_lmet() const      { return the_leading_lepton.DeltaPhi(the_met);}
  float Dphi_ljet4() const     { return the_leading_lepton.DeltaPhi(the_4thjet);}
  float Dphi_tops() const      { return top_hadronic.DeltaPhi(top_leptonic);}
  float Deta_lth() const      
  { 
      float deta_e_tophad= the_leading_lepton.Eta() - top_hadronic.Eta();
      if (deta_e_tophad<0.) deta_e_tophad*=-1.;
      return deta_e_tophad; 
  }
  float Deta_ljet4() const     
  { 
      float deta_e_jet4= the_leading_lepton.Eta() - the_4thjet.Eta();
      if (deta_e_jet4<0.) deta_e_jet4*=-1.;
      return deta_e_jet4; 
  }
  float Met() const            { return the_met.Pt(); }

  float HT() const      
  { 
      float HT_ = 0;
      for (unsigned int i = 0 ; i < jetsAna.size() ; i++)
          HT_ += jetsAna[i].p4.Pt();

      return HT_;
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
                                    the_leading_lepton,
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

  float M3b()
  {
      TLorentzVector sum;
      if(jetsAna.size() == 3)
      {
          sum = jetsAna[0].p4 + jetsAna[1].p4 + jetsAna[2].p4;
      }
      else
      {
          // check which jet is closest to lepton, then take other 3
          double dphimin = 99.;
          int index_closest_jet = -1;
          for(int i=0; i < 4; i++)
          {
              double dphi = the_leading_lepton.DeltaPhi(jetsAna[i].p4);
              if (dphi < dphimin)
              {
                  dphimin = dphi;
                  index_closest_jet = i;
              }
          }

          for(int i=0; i<4; i++)
          {
              if(i!=index_closest_jet)
                  sum = sum + jetsAna[i].p4;
          }

      }    

      return sum.M();    

  }






  TLorentzVector getTheLeadingLepton()      { return the_leading_lepton; };
  int            getTheLeadingLeptonPDGId() { return the_leading_lepton_pdgid; };
  TLorentzVector getTheSecondLepton()       { return the_second_lepton; };
  int            getTheSecondLeptonPDGId()  { return the_second_lepton_pdgid; };


  int getTheNumberOfSelectedLeptons() { return numberOfSelectedLeptons; }


  // -------------------------------------------------------------
  //                        data members
  // -------------------------------------------------------------
 private:
 

  TLorentzVector all_hadronic;
  TLorentzVector top_hadronic;
  TLorentzVector the_met;
  TLorentzVector the_leading_lepton;
  TLorentzVector the_second_lepton;
  TLorentzVector the_4thjet;  
  TLorentzVector w_leptonic;
  TLorentzVector top_leptonic;
  int the_leading_lepton_pdgid;
  int the_second_lepton_pdgid;
  
  //Objects for analysis
  std::vector<IPHCTree::NTElectron> electronsAna;
  std::vector<IPHCTree::NTMuon>     muonsAna;
  std::vector<IPHCTree::NTJet>      jetsAna;
  std::vector<IPHCTree::NTJet>      bjetsAna;

  //FactorizedJetCorrector* JetCorrector;
  int numberOfSelectedLeptons;
};

#endif
