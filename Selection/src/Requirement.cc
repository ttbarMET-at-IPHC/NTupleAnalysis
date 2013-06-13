#include "Selection/interface/Requirement.h"


// ----------------------------------------------------------------------------
// Constructor without arguments
// ----------------------------------------------------------------------------
Requirement::Requirement()
{
  // Jet requirements
  JetPtThreshold_     = 0.;
  JetEtaThreshold_    = 999.;

  // Muon requirements
  MuonPtThreshold_    = 0.;
  MuonEtaThreshold_   = 999.;
  MuonRelIso_         = 999.;
  MuonNofValidHits_   = 0;
  MuonNofValidTrHits_ = 10;
  MuonD0Cut_          = 0.;
  MuonNormChi2_       = 0;
  MuonVertexMatchThr_ = 999.;

  // Electron requirements
  ElectronPtThreshold_    = 0.;
  ElectronEtaThreshold_   = 999.;
  ElectronRelIso_         = 999.;
  ElectronD0Cut_          = 0.;
  ElectronVertexMatchThr_ = 999.;
  ElectronETSCThr_        = 0;
  DRemuThr_               = 999.;

  // Tau requirements
  TauPtThreshold_   = 0.;
  TauEtaThreshold_  = 999.;
  TauLeadTrkPtCut_  = 0;
  TauVertexMatchThr_= 999.;

  // Vertex requirements
  VertexNdofThr_    = 0;
  VertexZThr_       = 999.;
  VertexRhoThr_     = 999.;

  // MET
  METThreshold_     = 0;
}
