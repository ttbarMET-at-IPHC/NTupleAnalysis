#ifndef SSDiLeptonSelection_h
#define SSDiLeptonSelection_h

// IPHC headers
#include "Selection/interface/Selection.h"
#include "Selection/interface/SelectionTable.h"
#include "Tools/interface/Dataset.h"

// system include files
#include <memory>
#include <vector>
#include <sstream>

/**
	Steps of the selection: (integer returned by doFullSelection() or FillTable(...))
	- Step 1        MC-match
	- Step 2        Trigger
	- Step 3  	Dilepton pair choice
	- Step 4 	Z mass veto 
	- Step 5 	Minimal jet multiplicity 
	- Step 6 	MET cuts
	- Step 7 	btagging cuts

*/

class SSDiLeptonSelection: public Selection
{

  // -------------------------------------------------------------
  //                       method members
  // -------------------------------------------------------------
 public:

  //! Constructor without argument
  SSDiLeptonSelection();

  //! Copy constructor
  SSDiLeptonSelection(const SSDiLeptonSelection &);

  //! Destructor
  ~SSDiLeptonSelection()
  { }

  //! Mutator for parameters
  void SetParameters(float MinValue,
                     std::pair<float,float> METCuts,
                     std::pair<float,float> ZMassWindow,
                     int btagAlgo,
                     float btagDiscriCut,
                     int NofBtagJets);

  //! Accessor to cut list
  std::vector<std::string> GetCutList() const {return cuts_;}

  //! Accessor to channel list
  std::vector<std::string> GetChannelList() const {return channels_;}

  //! Accessor to Channel
  signed int GetChannel(const std::string& CandPairType) const
  {
    for (unsigned int i=0; i<channels_.size(); i++)
    {
      if (channels_[i]==CandPairType) return static_cast<signed int>(i);
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

  bool GetLeptonPair(const std::vector<IPHCTree::NTMuon>& muon_in,const std::vector<IPHCTree::NTElectron>& elec_in, std::vector<IPHCTree::NTMuon>& muon_out,std::vector<IPHCTree::NTElectron>& elec_out,  string& CandPairType, bool isForMM = false, float iso1_in = -1., float iso2_in = -1.);
  /** muon_out & elec_out will be filled with the 2 di-leptons candidates \n Returns true if a lepton pair is found \n
      CandPairType = "ee" or "emu" or "mumu" or "false"
	*/
  bool GetLeptonPair(std::vector<IPHCTree::NTMuon>& muon_out,std::vector<IPHCTree::NTElectron>& elec_out,  string& CandPairType); /** Idem with default GetSelected Electrons & Muons as input*/
  bool GetLeptonPairForMM(std::vector<IPHCTree::NTMuon>& muon_out,std::vector<IPHCTree::NTElectron>& elec_out,  string& CandPairType, float iso1_in = -1., float iso2_in = -1.); /** Idem with default GetSelectedNoIso Electrons & Muons as input*/


  int FillTable(SelectionTable& selTable, Dataset* dataset, int idataset, float weight); /** Fill the selectionTable according to the result of doFullSelection  for an 
                                                                                             event of weight "weight" of a given dataset idataset - Returns the integer of doFullSelection() */


  bool TestIsolationOfPair(float iso1_in , float iso2_in, std::vector<IPHCTree::NTMuon> muon_in,std::vector<IPHCTree::NTElectron> elec_in); /** It allows to test the isolation of the pair of leptons selected by GetLeptonPairForMM (or by GetLeptonPair).*/


  float DiLeptonMass(const std::vector<IPHCTree::NTMuon>& muons_cand, const std::vector<IPHCTree::NTElectron>& electrons_cand); /** Return the mass of the di-lepton candidate */
  float DiLeptonMT(const std::vector<IPHCTree::NTMuon>& muons_cand, const std::vector<IPHCTree::NTElectron>& electrons_cand); /** Return MT of the di-lepton candidate */
  TLorentzVector DiLeptonCand(const std::vector<IPHCTree::NTMuon>& muons_cand, const std::vector<IPHCTree::NTElectron>& electrons_cand); /** Return the p4 of the di-lepton candidate */
  bool DiLeptonMassCut(float MinValue, pair<float,float> ZMassWindow, const std::vector<IPHCTree::NTMuon>& muons_cand, const std::vector<IPHCTree::NTElectron>& electrons_cand,string channelName); /** Return true if dilepton candidate passes the  mass cuts */
  bool DiLeptonMassCut(const std::vector<IPHCTree::NTMuon>& muons_cand, const std::vector<IPHCTree::NTElectron>& electrons_cand, string channelName); /** Idem with defaults mass cuts */

  /**
   * return a integer which correspond to the last step that the event passes in the selection 
   * - possibility to check if the candPair correspond to the correct channel 
   *  & to cut at gen level for a given channel 
   * -  compute also the weight associated to btag
   */

  int doFullSelection(Dataset* dataset, string channelName=string(""), bool print = false, 
          bool isForMM = false , float iso1_in = -1., float iso2_in = -1., 
	  bool applyJES = false, float JESParam = 1., bool applyEES = false, float EESParam = 1., bool applyMES = false, float MESParam = 1., bool applyJER = false, float JERFactor = 0., bool applyMETS = false, float METScale = 1.);


  bool passTriggerSelection(Dataset* dataset, string channelName = string(""));
     
  float GetMinValueMassCut();
  pair<float,float> GetMETCut();
  pair<float,float> GetZmassWindowCut();
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
  int GetSelCode() const {return selCode_;};

  static const int maskTriggerCut = 0x1;
  static const int maskPairCut = 0x2;
  static const int maskZvetoCut = 0x4;
  static const int maskJetCut = 0x8;
  static const int maskMETCut = 0x10;
  static const int maskBjetsCut = 0x20;

  /**
   * Dilepton type. It is only filled during the doFullSelection method.
   */
  string pairType() const {return pairType_;}

  /** Invariant mass of the two leptons.
   *  It is only filled during the doFullSelection method.
   */
  float dileptonMass() const {return dimass_;}

  /**
   * Tells whether a particular jet passes the b-tagging requirement
   */
  bool passBtagSelection(const IPHCTree::NTJet & jet) const;

  /**
   *  Returns the b-tagging discriminant of a particular jet according to
   *  the algorithm specified in the config file.
   */

  double getBtagDiscr(const IPHCTree::NTJet & jet) const;
      
      
  double getLeptonScaleFactor(double pt1, double eta1, double pt2, double eta2, string channel);
  double getLeptonScaleFactorError(double pt1, double eta1, double pt2, double eta2, string channel);


  // -------------------------------------------------------------
  //                        data members
  // -------------------------------------------------------------
 private:
 
  // parameters
  float                    MinMassCut_;
  std::pair<float,float>   ZMassWindow_;
  std::pair<float,float>   METCuts_;
  std::vector<std::string> cuts_;
  std::vector<std::string> channels_;
  int                      btagAlgo_;
  float                    btagDiscriCut_;
  int                      NofBtagJets_;
  std::string              pairType_;

  //! invariant mass of the two leptons
  float dimass_;

  //! 32-bit outcome of the full selection
  int selCode_;

  //Objects for analysis
  std::vector<IPHCTree::NTElectron> electronsAna;
  std::vector<IPHCTree::NTMuon>     muonsAna;
  std::vector<IPHCTree::NTJet>      jetsAna;
  std::vector<IPHCTree::NTJet>      bjetsAna;

};

#endif
