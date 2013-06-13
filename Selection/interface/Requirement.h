#ifndef IPHC_Requirement_h
#define IPHC_Requirement_h

// STL headers
#include <memory>
#include <vector>
#include <string>

class Selection;
class AnalysisEnvironmentLoader;

//! \class Requirement
class Requirement
{

  // -------------------------------------------------------------
  //                       method members
  // -------------------------------------------------------------
 public:

  friend class Selection;
  friend class AnalysisEnvironmentLoader;

  //! Constructor without argument
  Requirement();

  //! Destructor
  ~Requirement()
  {}

  //! Set electron requirements
  void SetElectronRequirements(float PtThr, float EtaThr,
                               float RelIsoThr, float D0Thr,
                               float VertexMatchThr, float ElectronETSCThr,
                               float DRemuThr)
  {
    ElectronPtThreshold_    = PtThr;
    ElectronEtaThreshold_   = EtaThr;
    ElectronRelIso_         = RelIsoThr;
    ElectronD0Cut_          = D0Thr;
    ElectronVertexMatchThr_ = VertexMatchThr;
    ElectronETSCThr_        = ElectronETSCThr;
    DRemuThr_               = DRemuThr;
  }

  //! Set muon requirements
  void SetMuonRequirements(float PtThr, float EtaThr,
                           float RelIsoThr, float D0Thr,
                           float VertexMatchThr, float NValidHitsThr,
                           float NValidTkHitsThr, float Chi2Thr)
  {
    MuonPtThreshold_     = PtThr;
    MuonEtaThreshold_    = EtaThr;
    MuonRelIso_          = RelIsoThr;
    MuonD0Cut_           = D0Thr;
    MuonVertexMatchThr_  = VertexMatchThr;
    MuonNofValidTrHits_  = NValidTkHitsThr;
    MuonNofValidHits_    = NValidHitsThr;
    MuonNormChi2_        = Chi2Thr;
  }

  //! Set tau requirements
  void SetTauRequirements(float PtThr, float EtaThr,
                          float VertexMatchThr, float LeadTrkPtThr)
  {
    TauPtThreshold_    = PtThr;
    TauEtaThreshold_   = EtaThr;
    TauVertexMatchThr_ = VertexMatchThr;
    TauLeadTrkPtCut_   = LeadTrkPtThr;
  }

  //! Set jet requirements
  void SetJetRequirements(float PtThr, float EtaThr)
  {
    JetPtThreshold_=PtThr;
    JetEtaThreshold_=EtaThr;
  }

  //! Set MET requirements
  void SetMETRequirements(float PtThr)
  {
    METThreshold_=PtThr;
  }

  //! Set vertex requirements
  void SetVertexRequirements(float VertexNdofThr, float VertexZThr,
                             float VertexRhoThr)
  {
    VertexNdofThr_ = VertexNdofThr;
    VertexZThr_    = VertexZThr;
    VertexRhoThr_  = VertexRhoThr;
  }
   

  // -------------------------------------------------------------
  //                        data members
  // -------------------------------------------------------------
 private:

  //! Muon requirements
  float MuonPtThreshold_;
  float MuonEtaThreshold_;
  float MuonRelIso_;
  float MuonNofValidTrHits_;
  float MuonNofValidHits_;
  float MuonD0Cut_;
  float MuonNormChi2_;
  float MuonVertexMatchThr_;

  //! Electron requirements
  float ElectronPtThreshold_;
  float ElectronEtaThreshold_;
  float ElectronRelIso_;
  float ElectronD0Cut_;
  float ElectronVertexMatchThr_;
  float ElectronETSCThr_;
  float DRemuThr_;

  //! Tau requirements
  float TauPtThreshold_;
  float TauEtaThreshold_;
  float TauLeadTrkPtCut_;
  float TauVertexMatchThr_;

  //! Vertex requirements
  float VertexNdofThr_;
  float VertexZThr_;
  float VertexRhoThr_;

  //! Jet requirements
  float JetPtThreshold_;
  float JetEtaThreshold_;

  //! MET
  float METThreshold_;
};

#endif
