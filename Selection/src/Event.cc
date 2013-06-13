#include "Selection/interface/Event.h"


// ----------------------------------------------------------------------------
// Copy constructor
// ----------------------------------------------------------------------------
Event::Event(const Event& evt)
{
  general_   = evt.general_;
  mc_        = evt.mc_;
  pileup_    = evt.pileup_;
  trigger_   = evt.trigger_;
  jets_      = evt.jets_;
  heavyTagJets_   = evt.heavyTagJets_;
  met_       = evt.met_;
  photons_   = evt.photons_;
  electrons_ = evt.electrons_;
  muons_     = evt.muons_;
  taus_      = evt.taus_;
  tracks_    = evt.tracks_; 
  vertices_  = evt.vertices_; 
 
  PhotonType_   = evt.PhotonType_;
  JetMetType_   = evt.JetMetType_;
  ElectronType_ = evt.ElectronType_;
  MuonType_     = evt.MuonType_;
  TauType_      = evt.TauType_;
  TrackType_    = evt.TrackType_;
  VertexType_   = evt.VertexType_;

  PhotonEnabled_   = evt.PhotonEnabled_;
  JetMetEnabled_   = evt.JetMetEnabled_;
  ElectronEnabled_ = evt.ElectronEnabled_;
  MuonEnabled_     = evt.MuonEnabled_;
  TauEnabled_      = evt.TauEnabled_;
  TrackEnabled_    = evt.TrackEnabled_;
  VertexEnabled_   = evt.VertexEnabled_;

}


// ----------------------------------------------------------------------------
// Reset
// ----------------------------------------------------------------------------
void Event::Reset()
{
  general_         = 0;
  mc_              = 0;
  pileup_          = 0;
  trigger_         = 0;
  photons_         = 0;
  jets_            = 0;
  heavyTagJets_    = 0;
  met_             = 0;
  electrons_       = 0;
  muons_           = 0;
  taus_            = 0;
  tracks_          = 0; 
  pfcandidates_    = 0; 
  vertices_        = 0; 
 
  PhotonType_      = "";
  JetMetType_      = "";
  ElectronType_    = "";
  MuonType_        = "";
  TauType_         = "";
  TrackType_       = "";
  PFCandidateType_ = "";
  VertexType_      = "";

  PhotonEnabled_   = false;
  JetMetEnabled_   = false;
  ElectronEnabled_ = false;
  MuonEnabled_     = false;
  TauEnabled_      = false;
  TrackEnabled_    = false;
  PFCandidateEnabled_    = false;
  VertexEnabled_   = false;

}


// ----------------------------------------------------------------------------
// LoadEvent
// ----------------------------------------------------------------------------
bool Event::LoadEvent(const IPHCTree::NTEvent* evt)
{
  
  // safety : input pointer is null ?
  if (evt==0)
  {
    std::cerr << "The pointer to NTuple event is null !!!!" << std::endl;
    return false; 
  }

  // Declaring indicator for the success of event reading
  bool success=true;

  // general info
  general_ = &(evt->general);
  mc_      = &(evt->mc);
  trigger_ = &(evt->trigger);
  pileup_  = &(evt->pileup);

  // get only one photon collection 
  if (PhotonEnabled_)
  {
    photons_ = evt->photons.GetCollection(PhotonType_);
    if(photons_==0)
    {
      success=false;
      std::cerr<<"The photon collection called '" + PhotonType_ + "' is not found !"<<std::endl;
    }
  }

  // get only one Jet collection 
  if (JetMetEnabled_)
  {
    jets_ = evt->jets.GetCollection(JetMetType_);
    if(jets_==0)
    {
      success=false;
      std::cerr<<"The jet collection called '" + JetMetType_ + "' is not found !"<<std::endl;
      std::set<std::string> names;
      evt->jets.GetCollectionList(names);
      std::cerr << "Available collections are : "; 
      for (std::set<std::string>::const_iterator it=names.begin();it!=names.end();it++)
      {
        std::cerr << *it;
        if (it!=names.begin()) std::cerr << " , ";
      }
      std::cerr << std::endl;
    }
  }
  // get only one Met collection 
  if (JetMetEnabled_)
  {
    const std::vector<IPHCTree::NTMET>* mets = 
      evt->met.GetCollection(JetMetType_);
    if(mets==0)
    {
      success=false;
      std::cerr<<"The MET collection called '" + JetMetType_ + "' is not found !"<<std::endl;
      std::set<std::string> names;
      evt->met.GetCollectionList(names);
      std::cerr << "Available collections are : "; 
      for (std::set<std::string>::const_iterator it=names.begin();it!=names.end();it++)
      {
        std::cerr << *it;
        if (it!=names.begin()) std::cerr << " , ";
      }
      std::cerr << std::endl;
    }
    else
    {
    // met_ = &((*mets)[0]);
      met_ = &((*mets)[1]);
    }
  }

  // get one HeavyTagJet collection 
  if (HeavyTagJetEnabled_)
  {
    heavyTagJets_ = evt->jets.GetCollection(HeavyTagJetType_);
    if(heavyTagJets_==0)
    {
      success=false;
      std::cerr<<"The heavyTagJet collection called '" + HeavyTagJetType_ + "' is not found !"<<std::endl;
      std::set<std::string> names;
      evt->jets.GetCollectionList(names);
      std::cerr << "Available collections are : "; 
      for (std::set<std::string>::const_iterator it=names.begin();it!=names.end();it++)
      {
        std::cerr << *it;
        if (it!=names.begin()) std::cerr << " , ";
      }
      std::cerr << std::endl;
    }
  }

  // get only one electron collection 
  if (ElectronEnabled_)
  {
    electrons_ = evt->electrons.GetCollection(ElectronType_);
    if(electrons_==0)
    {
      success=false;
      std::cerr<<"The electron collection called '" + ElectronType_ + "' is not found !"<<std::endl;
      std::set<std::string> names;
      evt->electrons.GetCollectionList(names);
      std::cerr << "Available collections are : "; 
      for (std::set<std::string>::const_iterator it=names.begin();it!=names.end();it++)
      {
        std::cerr << *it;
        if (it!=names.begin()) std::cerr << " , ";
      }
      std::cerr << std::endl;
    }
  }

  // get only one muon collection 
  if (MuonEnabled_)
  {
    muons_ = evt->muons.GetCollection(MuonType_);
    if(muons_==0)
    {
      success=false;
      std::cerr<<"The muon collection called '" + MuonType_ + "' is not found !"<<std::endl;
      std::set<std::string> names;
      evt->muons.GetCollectionList(names);
      std::cerr << "Available collections are : "; 
      for (std::set<std::string>::const_iterator it=names.begin();it!=names.end();it++)
      {
        std::cerr << *it;
        if (it!=names.begin()) std::cerr << " , ";
      }
      std::cerr << std::endl;
    }
  }

  // get only one tau collection 
  if (TauEnabled_)
  {
    taus_ = evt->taus.GetCollection(TauType_);
    if(taus_==0)
    {
      success=false;
      std::cerr<<"The tau collection called '" + TauType_ + "' is not found !"<<std::endl;
      std::set<std::string> names;
      evt->taus.GetCollectionList(names);
      std::cerr << "Available collections are : "; 
      for (std::set<std::string>::const_iterator it=names.begin();it!=names.end();it++)
      {
        std::cerr << *it;
        if (it!=names.begin()) std::cerr << " , ";
      }
      std::cerr << std::endl;
    }
  }

  // get only one track collection 
  if (TrackEnabled_)
  {
    tracks_ = evt->tracks.GetCollection(TrackType_); 
    if(tracks_==0)
    {
      success=false;
      std::cerr<<"The track collection called '" + TrackType_ + " is not found !"<<std::endl;
      std::set<std::string> names;
      evt->tracks.GetCollectionList(names);
      std::cerr << "Available collections are : "; 
      for (std::set<std::string>::const_iterator it=names.begin();it!=names.end();it++)
      {
        std::cerr << *it;
        if (it!=names.begin()) std::cerr << " , ";
      }
      std::cerr << std::endl;
    }
  }

  // get only one pfcandidate collection 
  if (PFCandidateEnabled_)
  {
    pfcandidates_ = evt->pfcandidates.GetCollection(PFCandidateType_); 
    if(pfcandidates_==0)
    {
      success=false;
      std::cerr<<"The pfcandidate collection called '" + PFCandidateType_ + " is not found !"<<std::endl;
      std::set<std::string> names;
      evt->pfcandidates.GetCollectionList(names);
      std::cerr << "Available collections are : "; 
      for (std::set<std::string>::const_iterator it=names.begin();it!=names.end();it++)
      {
        std::cerr << *it;
        if (it!=names.begin()) std::cerr << " , ";
      }
      std::cerr << std::endl;
    }
  }

  // get only one vertex collection 
  if (VertexEnabled_)
  {
	vertices_ = evt->vertices.GetCollection(VertexType_); 
    if(vertices_==0)
    {
      success=false;
      std::cerr<<"The vertex collection called '" + VertexType_ + "' is not found !"<<std::endl;
      std::set<std::string> names;
      evt->vertices.GetCollectionList(names);
      std::cerr << "Available collections are : "; 
      for (std::set<std::string>::const_iterator it=names.begin();it!=names.end();it++)
      {
        std::cerr << *it;
        if (it!=names.begin()) std::cerr << " , ";
      }
      std::cerr << std::endl;
    }
  }

  return success;

}
