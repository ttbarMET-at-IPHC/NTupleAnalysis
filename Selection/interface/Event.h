#ifndef IPHC_Event_h
#define IPHC_Event_h

// IPHC headers
#include "NTFormat/interface/NTEvent.h"

// STL headers
#include <memory>
#include <vector>
#include <string>

//! \class Event
class Event
{

  // -------------------------------------------------------------
  //                       method members
  // -------------------------------------------------------------
 public:

  //! Constructor without argument
  Event() 
  { Reset(); }

  //! Copy constructor
  Event(const Event&);

  //! Destructor
  ~Event()
  {}
   
  //! Clear collections but parameters are cleared
  void Reset();

  //! Collections of objects are updated
  bool LoadEvent(const IPHCTree::NTEvent*event);


  // ------------- accessor to collection ( copy ) ---------------

  //! Get a copy of the vertex collection
  std::vector<IPHCTree::NTVertex> GetVertex() const
  { return (*vertices_); }

  //! Get a copy of the track collection
  std::vector<IPHCTree::NTTrack> GetTracks() const
  { return (*tracks_); }

  //! Get a copy of the pfcandidate collection
  std::vector<IPHCTree::NTPFCandidate> GetPFCandidates() const
  { return (*pfcandidates_); }

  //! Get a copy of the photon collection
  std::vector<IPHCTree::NTPhoton> GetPhotons() const
  { return (*photons_); }

  //! Get a copy of the jet collection
  std::vector<IPHCTree::NTJet> GetJets() const
  { return (*jets_); }
  
  //! Get a copy of the heavyTag jets collection
  std::vector<IPHCTree::NTJet> GetHeavyTagJets() const
  { return (*heavyTagJets_); }

  //! Get a copy of the MET collection
  IPHCTree::NTMET GetMET() const
  { return (*met_); }

  //! Get a copy of the electron collection
  std::vector<IPHCTree::NTElectron> GetElectrons() const
  { return (*electrons_); }

  //! Get a copy of the muon collection
  std::vector<IPHCTree::NTMuon> GetMuons() const   
  { return (*muons_); }

  //! Get a copy of the tau collection
  std::vector<IPHCTree::NTTau> GetTaus() const   
  { return (*taus_); }

  //! Get a copy of the GenTaus collection
  std::vector<TLorentzVector> GetGenTaus() const
  { return mc_->Generatedtaus; }

  //! Get a copy of the GenATaus collection
  std::vector<TLorentzVector> GetGenATaus() const
  { return mc_->GeneratedAtaus; }

  // ------------- accessor to collection (pointer) --------------

  //! Get a const pointer to the vertex collection
  const std::vector<IPHCTree::NTVertex>* GetPointer2Vertex() const
  { return vertices_; }

  //! Get a const pointer to the track collection
  const std::vector<IPHCTree::NTTrack>* GetPointer2Tracks() const
  { return tracks_; }

  //! Get a const pointer to the pfcandidate collection
  const std::vector<IPHCTree::NTPFCandidate>* GetPointer2PFCandidates() const
  { return pfcandidates_; }

  //! Get a const pointer to the photon collection
  const std::vector<IPHCTree::NTPhoton>* GetPointer2Photons() const
  { return photons_; }

  //! Get a const pointer to the jet collection
  const std::vector<IPHCTree::NTJet>* GetPointer2Jets() const
  { return jets_; }
  
  //! Get a const pointer to the heavyTag jet collection
  const std::vector<IPHCTree::NTJet>* GetPointer2HeavyTagJets() const
  { return heavyTagJets_; }

  //! Get a const pointer to the MET collection
  const IPHCTree::NTMET* GetPointer2MET() const
  { return met_; }

  //! Get a const pointer to the electron collection
  const std::vector<IPHCTree::NTElectron>* GetPointer2Electrons() const
  { return electrons_; }

  //! Get a const pointer to the muon collection
  const std::vector<IPHCTree::NTMuon>* GetPointer2Muons() const   
  { return muons_; }

  //! Get a const pointer to the tau collection
  const std::vector<IPHCTree::NTTau>* GetPointer2Taus() const   
  { return taus_; }

  //! Get a const pointer to Monte Carlo data
  const IPHCTree::NTMonteCarlo* GetPointer2MC() const   
  { return mc_; }

  //! Get a const pointer to Monte Carlo data
  const IPHCTree::NTMonteCarlo* mc() const   
  { return mc_; }

  //! Get a const pointer to Trigger data
  const IPHCTree::NTTrigger* GetPointer2Trigger() const   
  { return trigger_; }

  // ------------- mutator to collection labels ------------------

  //! Initialize the JetMet collection label
  void SetJetMetCollectionLabel(const std::string& label)
  { JetMetType_=label; }
  
  //! Initialize the HeavyTagJet collection label
  void SetHeavyTagJetCollectionLabel(const std::string& label)
  { HeavyTagJetType_=label; }

  //! Initialize the photon collection label
  void SetPhotonCollectionLabel(const std::string& label)
  { PhotonType_=label; }

  //! Initialize the electron collection label
  void SetElectronCollectionLabel(const std::string& label)
  { ElectronType_=label; }

  //! Initialize the muon collection label
  void SetMuonCollectionLabel(const std::string& label)
  { MuonType_=label; }

  //! Initialize the tau collection label
  void SetTauCollectionLabel(const std::string& label)
  { TauType_=label; }

  //! Initialize the track collection label
  void SetTrackCollectionLabel(const std::string& label)
  { TrackType_=label; }

  //! Initialize the pfcandidate collection label
  void SetPFCandidateCollectionLabel(const std::string& label)
  { PFCandidateType_=label; }

  //! Initialize the vertex collection label
  void SetVertexCollectionLabel(const std::string& label)
  { VertexType_=label; }


  // ------------- acessor to collection labels ------------------

  //! Access to the JetMet collection label
  std::string GetJetMetCollectionLabel() const
  { return JetMetType_; }
  
  //! Access to the HeavyTagJet collection label
  std::string GetHeavyTagJetCollectionLabel() const
  { return HeavyTagJetType_; }

  //! Access to the photon collection label
  std::string GetPhotonCollectionLabel() const
  { return PhotonType_; }

  //! Access to the electron collection label
  std::string GetElectronCollectionLabel() const
  { return ElectronType_; }

  //! Access to the muon collection label
  std::string GetMuonCollectionLabel() const
  { return MuonType_; }

  //! Access to the tau collection label
  std::string GetTauCollectionLabel() const
  { return TauType_; }

  //! Access to the track collection label
  std::string GetTrackCollectionLabel() const
  { return TrackType_; }

  //! Access to the pfcandidate collection label
  std::string GetPFCandidateCollectionLabel() const
  { return PFCandidateType_; }

  //! Access to the vertex collection label
  std::string GetVertexCollectionLabel() const
  { return VertexType_; }

  // ------------- mutator to enable/disable collection ------------------

  //! Initialize the JetMet collection label
  void EnableJetMetCollection()
  { JetMetEnabled_=true; }
  
  //! Initialize the HeavyTagJet collection label
  void EnableHeavyTagJetCollection()
  { HeavyTagJetEnabled_=true; }

  //! Initialize the photon collection label
  void EnablePhotonCollection()
  { PhotonEnabled_=true; }

  //! Initialize the electron collection label
  void EnableElectronCollection()
  { ElectronEnabled_=true; }

  //! Initialize the muon collection label
  void EnableMuonCollection()
  { MuonEnabled_=true; }

  //! Initialize the tau collection label
  void EnableTauCollection()
  { TauEnabled_=true; }

  //! Initialize the track collection label
  void EnableTrackCollection()
  { TrackEnabled_=true; }

  //! Initialize the pfcandidate collection label
  void EnablePFCandidateCollection()
  { PFCandidateEnabled_=true; }

  //! Initialize the vertex collection label
  void EnableVertexCollection()
  { VertexEnabled_=true; }

  //! Initialize the JetMet collection label
  void DisableJetMetCollection()
  { JetMetEnabled_=false; }
  
  //! Initialize the HeavyTagJet collection label
  void DisableHeavyTagJetCollection()
  { HeavyTagJetEnabled_=false; }

  //! Initialize the photon collection label
  void DisablePhotonCollection()
  { PhotonEnabled_=false; }

  //! Initialize the electron collection label
  void DisableElectronCollection()
  { ElectronEnabled_=false; }

  //! Initialize the muon collection label
  void DisableMuonCollection()
  { MuonEnabled_=false; }

  //! Initialize the tau collection label
  void DisableTauCollection()
  { TauEnabled_=false; }

  //! Initialize the track collection label
  void DisableTrackCollection()
  { TrackEnabled_=false; }

  //! Initialize the pfcandidate collection label
  void DisablePFCandidateCollection()
  { PFCandidateEnabled_=false; }

  //! Initialize the vertex collection label
  void DisableVertexCollection()
  { VertexEnabled_=false; }

  UInt_t getRunNumber()   const
  {return general_->runNb;}

  UInt_t getEventNumber() const
  {return general_->eventNb;}

  UInt_t getNpu() const
  {return pileup_->intime_npu;}

  UInt_t getTnpv() const
  {return pileup_->Tnpv;}

  float getRho() const
  {return pileup_->rho_PUUE_dens;}

  //  double getEleHLTMatch(){return eleHLTMatch;};

  // -------------------------------------------------------------
  //                        data members
  // -------------------------------------------------------------
 private:

  // Pointer to stored data
  const IPHCTree::NTGeneral*               general_;
  const IPHCTree::NTMonteCarlo*            mc_;
  const IPHCTree::NTTrigger*               trigger_;
  const IPHCTree::NTPileUp*                pileup_;
  const std::vector<IPHCTree::NTJet>*      jets_;
  const std::vector<IPHCTree::NTJet>*      heavyTagJets_;
  const IPHCTree::NTMET*                   met_;
  const std::vector<IPHCTree::NTPhoton>*   photons_;
  const std::vector<IPHCTree::NTElectron>* electrons_;
  const std::vector<IPHCTree::NTMuon>*     muons_;
  const std::vector<IPHCTree::NTTau>*      taus_;
  const std::vector<IPHCTree::NTTrack>*    tracks_; 
  const std::vector<IPHCTree::NTPFCandidate>* pfcandidates_; 
  const std::vector<IPHCTree::NTVertex>*   vertices_; 

  // Collection name initialized by XML configurator 
  std::string PhotonType_;
  std::string JetMetType_;
  std::string HeavyTagJetType_;
  std::string ElectronType_;
  std::string MuonType_;
  std::string TauType_;
  std::string TrackType_;
  std::string VertexType_;
  std::string PFCandidateType_;

  // Collection Enabled
  bool PhotonEnabled_;
  bool JetMetEnabled_;
  bool HeavyTagJetEnabled_;
  bool ElectronEnabled_;
  bool MuonEnabled_;
  bool TauEnabled_;
  bool TrackEnabled_;
  bool VertexEnabled_;
  bool PFCandidateEnabled_;

  // Static empty collection if a collection 
  // is not stored in the NTuple
  static const std::vector<IPHCTree::NTJet>      empty_jets_;
  static const std::vector<IPHCTree::NTJet>      empty_heavyTagJets_;
  static const std::vector<IPHCTree::NTPhoton>   empty_photons_;
  static const std::vector<IPHCTree::NTElectron> empty_electrons_;
  static const std::vector<IPHCTree::NTMuon>     empty_muons_;
  static const std::vector<IPHCTree::NTTau>      empty_taus_;
  static const std::vector<IPHCTree::NTTrack>    empty_tracks_; 
  static const std::vector<IPHCTree::NTPFCandidate>    empty_pfcandidates_; 
  static const std::vector<IPHCTree::NTVertex>   empty_vertices_; 
};

#endif
