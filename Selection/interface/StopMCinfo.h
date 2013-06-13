#ifndef StopMCInfo_h
#define StopMCInfo_h

// IPHC headers
#include "NTFormat/interface/NTEvent.h"

#include <vector>
#include <string>

class StopMCinfo 
{

 public:

 enum ttChannel{unknown,fullyHad,eJets,muJets,tauJets,ee,emu,mumu,etau,mutau,tautau};

  StopMCinfo(); 
  ~StopMCinfo(){;}

 void LoadEvent(IPHCTree::NTEvent* event);
 void Reset();

 // Accessor to booleans

 bool isSUSYEvent() { return isSUSYEvent_; }
 bool isStop2topChi2Event() { return isStop2topChiEvent_; }
 bool isStop2bCharginoEvent() { return isStop2bCharginoEvent_; }

 ttChannel GetChannel() { return ttbarChannel_; }
 std::string PrintChannel()const;
 Int_t TMEME() { return TMEME_; }
 Long_t SCNTWMEME() {return SCNTWMEME_;}

 //To be completeed
 bool MassSelector(const float& mStop, const float& mNeutralino);
 bool MassRangeSelector(const float& mStopMin, const float& mStopMax, const float& mNeutralinoMin, const float& mNeutralinMax);

 //Output
 std::string Bool2Text(const bool& value) const;
 void Print() const;

 //Access to masses
 Float_t GetStopMass()const {return mStop_;}
 Float_t GetNeutralinoMass()const {return mNeutralino_;}
 Float_t GetCharginoMass()const {return mChargino_;}

 //Add accessors to the objects
 vector<IPHCTree::NTGenParticle*> GetStops() { return stop_; }
 vector<IPHCTree::NTGenParticle*> GetCharginos() { return charginos_; }
 vector<IPHCTree::NTGenParticle*> GetNeutralinos() { return neutralinos_; }
 vector<IPHCTree::NTGenParticle*> GetWs() { return W_; }
 vector<IPHCTree::NTGenParticle*> GetHadronicWs() { return HadronicW_; }
 vector<IPHCTree::NTGenParticle*> GetTops() { return top_; }
 vector<IPHCTree::NTGenParticle*> GetHadronicTops() { return HadronicTop_; }
 vector<IPHCTree::NTGenParticle*> GetLeptonsFromW() { return lepton_; }
 vector<IPHCTree::NTGenParticle*> GetNeutrinosFromW() { return neutrino_; }
 vector<IPHCTree::NTGenParticle*> GetHadronicTauFromW() { return hadronicTau_; }
 vector<IPHCTree::NTGenParticle*> GetLeptonicTauDecay() { return TauDecay_; }

 private:

  // Characterization of the event
  bool isSUSYEvent_;
  bool isStop2topChiEvent_;
  bool isStop2bCharginoEvent_;

  // Decay identification
  ttChannel ttbarChannel_;
  Int_t TMEME_;
  Long_t SCNTWMEME_;
  //S: number of stop
  //C: number of charginos
  //N: number of neutralinos
  //W: number of Ws

  // Masses
  Float_t mStop_;
  Float_t mNeutralino_;
  Float_t mChargino_;

  // Pointer to event
  IPHCTree::NTEvent* event_;

  // List of particles
  vector<IPHCTree::NTGenParticle*> stop_;
  vector<IPHCTree::NTGenParticle*> charginos_;
  vector<IPHCTree::NTGenParticle*> neutralinos_;
  vector<IPHCTree::NTGenParticle*> top_;
  vector<IPHCTree::NTGenParticle*> HadronicTop_;
  vector<IPHCTree::NTGenParticle*> W_;
  vector<IPHCTree::NTGenParticle*> HadronicW_;
  vector<IPHCTree::NTGenParticle*> lepton_;
  vector<IPHCTree::NTGenParticle*> neutrino_;
  vector<IPHCTree::NTGenParticle*> hadronicTau_;
  vector<IPHCTree::NTGenParticle*> TauDecay_;


};

#endif
