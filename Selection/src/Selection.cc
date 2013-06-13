#include "Selection/interface/Selection.h"
#include "Tools/interface/FileExists.h"

// ----------------------------------------------------------------------------
// Global functions
// ----------------------------------------------------------------------------
void InitVectorOfWeight(std::vector<float>& weightb)
{
  weightb.push_back (1.);
  weightb.push_back (0.);
  weightb.push_back (0.);
  weightb.push_back (0.);
  weightb.push_back (0.);
}

void LoadBWeight(Selection& sel, 
                 std::vector<float>& weightb,
                 const std::vector<NTJet>& theselJets)
{
  std::vector<float> weight_temp = 
           sel.GetSFBweight().GetWeigth4BSel(sel.GetMethodb(),
                                             sel.GetSystb(),
                                             theselJets);
  weightb[0] = weight_temp[0];	//weight of the event
  weightb[1] = weight_temp[1];	//proba 0 jet
  weightb[2] = weight_temp[2];	//proba 1 jet
  weightb[3] = weight_temp[3];	//proba 2 jets
  weightb[4] = weight_temp[4];	//proba at least 3 jets
}

// ----------------------------------------------------------------------------
// Constructor without arguments
// ----------------------------------------------------------------------------
Selection::Selection()
{
  // Pile-Up weights
  PUWeights.clear();

  //! Scale factors
  scaleFactEl = 0;
  scaleFactMu = 0;

  // BTag variables
  flag_btagweight_=0;
  methodb_=0;
  systb_=0;

  // JES
  histo_jesunc_ = 0;
}

// ----------------------------------------------------------------------------
// Constructor without arguments
// ----------------------------------------------------------------------------
Selection::Selection(const Selection & s) : Event(s)
{
  // Pile-Up weights
  PUWeights=s.PUWeights;

  // Requirements
  cfg = s.cfg;

  // BTag variables
  flag_btagweight_        = s.flag_btagweight_;
  systb_                  = s.systb_;
  methodb_                = s.methodb_;
  sfb_                    = s.sfb_;

  // JES & Scale factors
  if (s.histo_jesunc_) histo_jesunc_  = (TH2F*) s.histo_jesunc_->Clone();
  if (s.scaleFactEl)   scaleFactEl	  = (TH2F*) s.scaleFactEl->Clone();
  if (s.scaleFactMu)   scaleFactMu	  = (TH2F*) s.scaleFactMu->Clone();
}

    
// ----------------------------------------------------------------------------
// GetSelectedVertex
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTVertex> Selection::GetSelectedVertex() const
{
  // Create container for output
  std::vector<IPHCTree::NTVertex> selectedVertex;

  // Get pointer to data stored in MiniTree
  const std::vector<IPHCTree::NTVertex>* vertices = GetPointer2Vertex();
  if (vertices==0) return selectedVertex;

  // Loop over stored objects
	for(unsigned int i=0;i<vertices->size();i++)
  {
    const IPHCTree::NTVertex* myvertex = &((*vertices)[i]);
		if (myvertex->isFake)                        continue;
		if (myvertex->ndof<=cfg.VertexNdofThr_)      continue;
		if (fabs(myvertex->p3.Z())> cfg.VertexZThr_) continue;
		if (myvertex->Rho()>= cfg.VertexRhoThr_)     continue;
		selectedVertex.push_back(*myvertex);
	}
	return selectedVertex;
} 
    


// ----------------------------------------------------------------------------
// GetCleanJets
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTJet> Selection::GetCleanJets(
                       std::vector<IPHCTree::NTJet> jet_cand, 
		       std::vector<IPHCTree::NTMuon> muon_cand, 
		       std::vector<IPHCTree::NTElectron> elec_cand) const
{

  std::vector<IPHCTree::NTJet> cleanJets;
 
  for(unsigned int i=0;i<jet_cand.size();i++){
      if(fabs(jet_cand[i].p4.Eta())> cfg.JetEtaThreshold_ || 
           jet_cand[i].p4.Pt()<cfg.JetPtThreshold_) continue;
      
      double deltaRmu = 10000;
      double deltaRel = 10000;
      
      for(unsigned int imu=0; imu< muon_cand.size(); imu++){
	double deltaR = jet_cand[i].p4.DeltaR(muon_cand[imu].p4);
	if(deltaR < deltaRmu) deltaRmu = deltaR;
      }
      
      for(unsigned int iel=0; iel< elec_cand.size(); iel++){
	double deltaR = jet_cand[i].p4.DeltaR(elec_cand[iel].p4);
	if(deltaR < deltaRel) deltaRel = deltaR;
      }
      
      if( deltaRmu > 0.4  && deltaRel > 0.4) cleanJets.push_back(jet_cand[i]);
      
  }
  std::sort(cleanJets.begin(),cleanJets.end(),HighestPt());
  return cleanJets;

}


// ----------------------------------------------------------------------------
// GetScaledJets
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTJet> Selection::GetScaledJets(float scale) const
{
  std::vector<IPHCTree::NTJet> newJets = GetJets();
  for (unsigned int i=0; i<newJets.size(); i++)
  {
    float newscale = 0;
    float eta  = newJets[i].p4.Eta();
    float pt   = newJets[i].p4.Pt();
    if (pt>=300.) pt = 300;
    if (scale<1.0)
    {
	//DEBUG_MSG << endl;
      newscale = (1.- histo_jesunc_->GetBinContent(
                           histo_jesunc_->FindBin(eta,pt)));

	  //DEBUG_MSG << endl;
    }
    else
    {
	//DEBUG_MSG << endl;
      newscale = (1.+ histo_jesunc_->GetBinContent(
                           histo_jesunc_->FindBin(eta,pt)));
	//DEBUG_MSG << endl;
    }
	//DEBUG_MSG << endl;
    newJets[i].p4.SetPxPyPzE(newscale*newJets[i].p4.Px(),
                             newscale*newJets[i].p4.Py(),
                             newscale*newJets[i].p4.Pz(),
                             newscale*newJets[i].p4.E());
  }
    
  return newJets;
}


// ----------------------------------------------------------------------------
// GetSmearedJets
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTJet>  Selection::GetSmearedJets(
                           const vector<IPHCTree::NTJet>& injets,
                           float jetResFactor) const
{
  std::vector<IPHCTree::NTJet> newJets = GetJets();
  for (unsigned int i=0; i<newJets.size(); i++)
  {
    double gen_pt  = newJets[i].p4Gen.Pt(); //to be changed
    double jet_pt  = newJets[i].p4.Pt();
    double deltapt = (jet_pt - gen_pt) * jetResFactor;
    double ptscale = ((jet_pt + deltapt) / jet_pt);
    if (ptscale<0) ptscale = 0;
    newJets[i].p4.SetPxPyPzE(ptscale*newJets[i].p4.Px(),
                             ptscale*newJets[i].p4.Py(),
                             ptscale*newJets[i].p4.Pz(),
                             ptscale*newJets[i].p4.E());
  }
  return newJets;
}


// ----------------------------------------------------------------------------
// GetScaleMET
// ----------------------------------------------------------------------------
IPHCTree::NTMET Selection::GetScaledMET(float scale) const
{
  // Create container for output
  IPHCTree::NTMET newMET;

  // Get pointer to data stored in MiniTree
  const IPHCTree::NTMET* met  = GetPointer2MET();
  const std::vector<IPHCTree::NTJet>* jets = GetPointer2Jets();
  if (jets==0 || met==0) return newMET;

  double missetX = 0.;
  double missetY = 0.;
  for (unsigned int i=0; i<jets->size(); i++)
  {
    missetX += (1-scale) * (*jets)[i].p4.Px();
    missetY += (1-scale) * (*jets)[i].p4.Py();
  }

  newMET = *met;
  newMET.p2.Set(met->p2.Px() + missetX,
                met->p2.Py() + missetY);
  return newMET;
}


// ----------------------------------------------------------------------------
// GetSmearedMET
// ----------------------------------------------------------------------------
IPHCTree::NTMET Selection::GetSmearedMET(
                         const std::vector<IPHCTree::NTJet>& injets,
                         float jetResFactor) const
{
  // Create container for output
  IPHCTree::NTMET newMET;

  // Get pointer to data stored in MiniTree
  const IPHCTree::NTMET* met  = GetPointer2MET();
  if (met==0) return newMET;

  double missetX = 0.;
  double missetY = 0.;
  for (unsigned int i=0; i<injets.size(); i++)
  {
    double gen_pt = injets[i].p4Gen.Pt(); //to be changed:
                                          //should be the pt of genjet !
    double jet_pt = injets[i].p4.Pt();
    double deltapt = (jet_pt - gen_pt) * jetResFactor;
    double ptscale = ((jet_pt + deltapt) / jet_pt);
    if(ptscale <0 ) ptscale = 0;
    missetX = (ptscale-1)*injets[i].p4.Px();
    missetY = (ptscale-1)*injets[i].p4.Py();
  }
  newMET = *met;
  newMET.p2.Set(met->p2.Px() + missetX,
                met->p2.Py() + missetY);
  return newMET;
}


// ----------------------------------------------------------------------------
// GetUnclusScaleMET
// ----------------------------------------------------------------------------
IPHCTree::NTMET Selection::GetUnclusScaledMET(bool applyUnclusScale,
                                              float scale) const
{
  // Create container for output
  IPHCTree::NTMET newMET;

  double missetX = 0.;
  double missetY = 0.;

  const IPHCTree::NTMET* met  = GetPointer2MET();
  if (met!=0)
  {
    missetX = met->p2.Px();
    missetY = met->p2.Py();
  }

  const std::vector<IPHCTree::NTJet>* jets  = GetPointer2Jets();
  if (jets!=0) for (unsigned int i=0; i<jets->size(); i++)
  {
    missetX -= (*jets)[i].p4.Px();
    missetY -= (*jets)[i].p4.Py();
  }

  const std::vector<IPHCTree::NTElectron>* electrons  = GetPointer2Electrons();
  for (unsigned int i=0; i<electrons->size(); i++)
  {
    missetX -= (*electrons)[i].p4.Px();
    missetY -= (*electrons)[i].p4.Py();
  }

  const std::vector<IPHCTree::NTMuon>* muons  = GetPointer2Muons();
  for (unsigned int i=0; i<muons->size(); i++)
  {
    missetX -= (*muons)[i].p4.Px();
    missetY -= (*muons)[i].p4.Py();
  }

  const std::vector<IPHCTree::NTTau>* taus = GetPointer2Taus();
  for (unsigned int i=0; i<taus->size(); i++)
  {
    missetX -= (*taus)[i].p4.Px();
    missetY -= (*taus)[i].p4.Py();
  }

  newMET = *met;
  newMET.p2.Set(met->p2.Px()+(scale-1)*missetX,
                met->p2.Py()+(scale-1)*missetY);
  return newMET;
}


// ----------------------------------------------------------------------------
// GetMET
// ----------------------------------------------------------------------------
IPHCTree::NTMET Selection::GetSelectedMET(bool applyJES, float scale, 
                                          bool applyJER, float ResFactor) const
{
	// WARNING: force to not apply the JES 
  // it's due to JES/JEC correction scenario in PAT
	applyJES = false;
	if(applyJES)
  {
    if(applyJER)
    {
			return GetSmearedMET(GetScaledJets(scale),ResFactor);
		}
		return GetScaledMET(scale); 
	}
  if(applyJER)
  {
		return GetSmearedMET(GetJets(), ResFactor);
	}
	else return GetMET();
}


//*******************************************************
//tools for Type1 MET calculated ourselves, on the flight.
//*******************************************************

IPHCTree::NTMET Selection::GetScaledType1MET(IPHCTree::NTMET &themet, float scale) const
{
    double missetX = 0.;
    double missetY = 0.;

    const std::vector<IPHCTree::NTJet>* jets  = GetPointer2Jets();
    for (unsigned int i=0; i<jets->size(); i++){
      missetX += (1-scale) * (*jets)[i].p4.Px();
      missetY += (1-scale) * (*jets)[i].p4.Py();
    }
    IPHCTree::NTMET newMET = themet;
    newMET.p2.Set( themet.p2.Px()+missetX, themet.p2.Py()+missetY);
    return newMET;
}

IPHCTree::NTMET  Selection::GetSmearedType1MET(IPHCTree::NTMET &themet, vector<IPHCTree::NTJet> injets, float jetResFactor) const
{
    double missetX = 0.;
    double missetY = 0.;

    for (unsigned int i=0; i<injets.size(); i++){
      double gen_pt = injets[i].p4Gen.Pt(); //to be changed: should be the pt of genjet !
      double jet_pt = injets[i].p4.Pt();
      double deltapt = (jet_pt - gen_pt) * jetResFactor;
      double ptscale = ((jet_pt + deltapt) / jet_pt);
      if(ptscale <0 ) ptscale = 0;
      missetX   = (ptscale-1)*injets[i].p4.Px();
      missetY   = (ptscale-1)*injets[i].p4.Py();
    }

    const IPHCTree::NTMET* met  = GetPointer2MET();
    IPHCTree::NTMET newMET = *met;
    newMET.p2.Set( themet.p2.Px()+missetX, themet.p2.Py()+missetY);
    return newMET;
}


IPHCTree::NTMET Selection::GetType1MET( vector<IPHCTree::NTJet> injets, bool applyJES, float scale, bool applyJER, float ResFactor) const{
	//WARNING: force to not apply the JES - it's due to JES/JEC correction scenario in PAT
	
	double missetX = 0.;
        double missetY = 0.;
 
        //JLA Info corrJES is not available in NTuples
	/*for (unsigned int i=0; i<injets.size(); i++){
	  if(injets[i].p4.Pt() < 10 ) continue;
          missetX   += - injets[i].corrJES*injets[i].p4.Px() + injets[i].p4.Px();
          missetY   += - injets[i].corrJES*injets[i].p4.Py() + injets[i].p4.Py();
        }*/
	std::cout<<" JLA WARNING : GetType1MET, corrJES not available."<<std::endl;
		
        const IPHCTree::NTMET* met  = GetPointer2MET();
	IPHCTree::NTMET themet = *met;
        themet.p2.Set( themet.p2.Px()+missetX, themet.p2.Py()+missetY);
	  
	
	if(applyJES){
        	if(applyJER){
			return GetSmearedType1MET(themet, GetScaledJets(scale),ResFactor );
		}
		return GetScaledType1MET(themet, scale); 
	}
        if(applyJER){
		return GetSmearedType1MET(themet, GetJets(), ResFactor);
	}
	else return themet;
}



IPHCTree::NTMET Selection::GetScaledType1METWithJER(vector<IPHCTree::NTJet> injets, bool applyJES, float scale, bool applyJER, float ResFactor) const{
    
    double missetX = 0.;
    double missetY = 0.;
 
    IPHCTree::NTMET themet = GetType1MET(injets, applyJES, scale, applyJER, ResFactor);
    
    double MCkfac[4];
    MCkfac[0] = 1.06177;
    MCkfac[1] = 1.08352;
    MCkfac[2] = 1.02911;
    MCkfac[3] = 1.15288;
    
    for (unsigned int i=0; i<injets.size(); i++){
	  if(injets[i].p4.Pt() < 10 ) continue;
	  double  DeltaPt = injets[i].p4.Pt() - injets[i].p4Gen.Pt();
	  double kfac = 0;
	  if(fabs(injets[i].p4.Eta()) < 1.1                                   ) kfac = MCkfac[0];
	  if(fabs(injets[i].p4.Eta()) > 1.1 && fabs(injets[i].p4.Eta()) < 1.7 ) kfac = MCkfac[1];
	  if(fabs(injets[i].p4.Eta()) > 1.7 && fabs(injets[i].p4.Eta()) < 2.3 ) kfac = MCkfac[2];
	  if(fabs(injets[i].p4.Eta()) >2.3                                    ) kfac = MCkfac[3];
	  
	  double pt_smear = injets[i].p4.Pt() + (kfac-1)*DeltaPt;
	   
	  double pt_ratio = pt_smear/injets[i].p4.Pt();
	  
          missetX   += -pt_ratio*injets[i].p4.Px() + injets[i].p4.Px();
          missetY   += -pt_ratio*injets[i].p4.Py() + injets[i].p4.Py();
          //missetX   += - injets[i].corrJES*injets[i].p4.Px() + injets[i].p4.Px();
          //missetY   += - injets[i].corrJES*injets[i].p4.Py() + injets[i].p4.Py();
    }
    
    
    IPHCTree::NTMET newMET = themet;  
    newMET.p2.Set( themet.p2.Px()+missetX, themet.p2.Py()+missetY);
    return newMET;
}


IPHCTree::NTMET Selection::GetGenMET(std::vector<IPHCTree::WDecaysMC> &wAndDecays)
{
    IPHCTree::NTMET theGenMET ;  
    TLorentzVector p4GenMET;
    for(unsigned int i = 0; i < wAndDecays.size(); i++){
       p4GenMET += (wAndDecays)[i].p4_Neu_gen;
    }
    
    theGenMET.p2.Set( -p4GenMET.Px(), -p4GenMET.Py());
    return theGenMET;
}

//*******************************************************


// ----------------------------------------------------------------------------
// GetSelectedJets
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTJet> Selection::GetSelectedJets(
            const std::vector<IPHCTree::NTMuon>& muon_cand,
            const std::vector<IPHCTree::NTElectron>& elec_cand,
            bool applyJES, float scale, 
            bool applyJER, float ResFactor) const
{
  // Containers for output
  std::vector<IPHCTree::NTJet> selectedJets;
  std::vector<IPHCTree::NTJet> scaledJets;

  // Get scaled jets
  if(applyJES) scaledJets = GetScaledJets(scale);
  else if (GetPointer2Jets()!=0) scaledJets = *GetPointer2Jets();

  // Apply JER
  if(applyJER) scaledJets = GetSmearedJets(scaledJets, ResFactor);
  for(unsigned int i=0;i<scaledJets.size();i++)
  {
    if (fabs(scaledJets[i].p4.Eta())> cfg.JetEtaThreshold_ ||
        scaledJets[i].p4.Pt()<cfg.JetPtThreshold_) continue;
	
    if (!scaledJets[i].ID["TIGHT"]) continue;
    
        
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
  std::sort(selectedJets.begin(),selectedJets.end(),HighestPt());
  return selectedJets;
}


// ----------------------------------------------------------------------------
// GetSelectedJets
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTJet> Selection::GetSelectedJets(
                  const std::vector<IPHCTree::NTMuon>& muon_cand,
                  const std::vector<IPHCTree::NTElectron>& elec_cand,
                  const std::vector<IPHCTree::NTTau>& tau_cand,
                  bool applyJES, float scale, 
                  bool applyJER, float ResFactor) const
{

  // Create containers for output
  std::vector<IPHCTree::NTJet> selectedJets;
  std::vector<IPHCTree::NTJet> scaledJets;

  // Get scaled jets
  if(applyJES) scaledJets = GetScaledJets(scale);
  else if (GetPointer2Jets()!=0) scaledJets = *GetPointer2Jets();

  // Apply JER
  if(applyJER) scaledJets = GetSmearedJets(scaledJets, ResFactor);

  for(unsigned int i=0;i<scaledJets.size();i++)
  {
    if (fabs(scaledJets[i].p4.Eta()) > cfg.JetEtaThreshold_ ||
        scaledJets[i].p4.Pt() < cfg.JetPtThreshold_) continue;
    
    double deltaRmu  = 10000;
    double deltaRel  = 10000;
    double deltaRtau = 10000;
    
    for(unsigned int imu=0; imu< muon_cand.size(); imu++)
    {
      double deltaR = scaledJets[i].p4.DeltaR(muon_cand[imu].p4);
      if(deltaR < deltaRmu) deltaRmu = deltaR;
    }
    
    for(unsigned int iel=0; iel< elec_cand.size(); iel++){
      double deltaR = scaledJets[i].p4.DeltaR(elec_cand[iel].p4);
      if(deltaR < deltaRel) deltaRel = deltaR;
    }
    
    for(unsigned int itau=0; itau< tau_cand.size(); itau++)
    {
      double deltaR = scaledJets[i].p4.DeltaR(tau_cand[itau].p4);
      if(deltaR < deltaRtau) deltaRtau = deltaR;
    }
    
    if( deltaRmu > 0.4  && deltaRel > 0.4 && deltaRtau > 0.4) 
      selectedJets.push_back(scaledJets[i]);
  }
  std::sort(selectedJets.begin(),selectedJets.end(),HighestPt());
  return selectedJets;
}

// ----------------------------------------------------------------------------
// GetSelectedJets
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTJet> Selection::GetSelectedJets(
                  const std::vector<IPHCTree::NTTau>& tau_cand, 
		  bool applyJES, float scale, 
		  bool applyJER, float ResFactor) const
{
  // Container for output
  std::vector<IPHCTree::NTJet> selectedJets;
  std::vector<IPHCTree::NTJet> scaledJets;

  // Get scaled jets
  if(applyJES) scaledJets = GetScaledJets(scale);
  else if (GetPointer2Jets()!=0) scaledJets = *GetPointer2Jets();
  if(applyJER) scaledJets = GetSmearedJets(scaledJets, ResFactor);
 
  for(unsigned int i=0;i<scaledJets.size();i++)
  {
    if(fabs(scaledJets[i].p4.Eta()) > cfg.JetEtaThreshold_ || 
        scaledJets[i].p4.Pt() < cfg.JetPtThreshold_) continue;
   
    double deltaRtau = 10000;
    
    for(unsigned int itau=0; itau< tau_cand.size(); itau++)
    {
      double deltaR = scaledJets[i].p4.DeltaR(tau_cand[itau].p4);
      if(deltaR < deltaRtau) deltaRtau = deltaR;  
    }
    
    if(deltaRtau > 0.5) selectedJets.push_back(scaledJets[i]);
	
  }
  std::sort(selectedJets.begin(),selectedJets.end(),HighestPt());
  return selectedJets;
}

// ----------------------------------------------------------------------------
// GetSelectedJets
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTJet> Selection::GetSelectedJets(
                               bool applyJES, float scale, 
                               bool applyJER, float ResFactor) const
{
  // Container for output
  std::vector<IPHCTree::NTJet> selectedJets;
  std::vector<IPHCTree::NTJet> scaledJets;

  // Get scaled jets
  if(applyJES) scaledJets = GetScaledJets(scale);
  else if (GetPointer2Jets()!=0) scaledJets = *GetPointer2Jets();
  if(applyJER) scaledJets = GetSmearedJets(scaledJets, ResFactor);

  for(unsigned int i=0;i<scaledJets.size();i++)
  {
    if (fabs(scaledJets[i].p4.Eta())> cfg.JetEtaThreshold_ ||
        scaledJets[i].p4.Pt()< cfg.JetPtThreshold_) continue;   
    selectedJets.push_back(scaledJets[i]);
	
  }
  std::sort(selectedJets.begin(),selectedJets.end(),HighestPt());
  return selectedJets;
}


// ----------------------------------------------------------------------------
// GetSelectedJets
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTJet> Selection::GetSelectedJets(
                                float PtThr, float EtaThr,
                                bool applyJES, float scale,
                                bool applyJER, float ResFactor) const
{
  // Container for output
  std::vector<IPHCTree::NTJet> selectedJets;
  std::vector<IPHCTree::NTJet> scaledJets;

  // Get scaled jets
  if(applyJES) scaledJets = GetScaledJets(scale);
  else if (GetPointer2Jets()!=0) scaledJets = *GetPointer2Jets();
  if(applyJER) scaledJets = GetSmearedJets(scaledJets, ResFactor);

  for(unsigned int i=0;i<scaledJets.size();i++)
  {
   	if (fabs(scaledJets[i].p4.Eta())<EtaThr &&
        scaledJets[i].p4.Pt()>PtThr)
      selectedJets.push_back(scaledJets[i]);
  }
  std::sort(selectedJets.begin(),selectedJets.end(),HighestPt());
  return selectedJets;
}

// ----------------------------------------------------------------------------
// GetSelectedJetsForLJets
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTJet> Selection::GetSelectedJetsForLJets(
                               bool applyJES, float scale, 
                               bool applyJER, float ResFactor) const
{
  // Container for output
  std::vector<IPHCTree::NTJet> selectedJets;
  std::vector<IPHCTree::NTJet> scaledJets;

  // Get scaled jets
  if(applyJES) scaledJets = GetScaledJets(scale);
  else if (GetPointer2Jets()!=0) scaledJets = *GetPointer2Jets();
  if(applyJER) scaledJets = GetSmearedJets(scaledJets, ResFactor);

  for(unsigned int i=0;i<scaledJets.size();i++)
  {
    if (scaledJets[i].ID["LOOSE"] != 1.) continue;
    if (fabs(scaledJets[i].p4.Eta())> cfg.JetEtaThreshold_ ||
        scaledJets[i].p4.Pt()< cfg.JetPtThreshold_) continue;   
    selectedJets.push_back(scaledJets[i]);
	
  }
  std::sort(selectedJets.begin(),selectedJets.end(),HighestPt());
  return selectedJets;
}

// ----------------------------------------------------------------------------
// GetSelectedBJets
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTJet> Selection::GetSelectedBJets(
                    const std::vector<IPHCTree::NTJet>& SelectedJets,
                    const int& algo, const float & discricut ) const
{
  std::vector<IPHCTree::NTJet> btagjets;
  for(unsigned int j=0;j<SelectedJets.size();j++)
  {
    if (passBtagSelection(SelectedJets[j], algo, discricut)) 
      btagjets.push_back(SelectedJets[j]);
  }
  return btagjets;
}

bool Selection::passBtagSelection(const NTJet & jet, const int& algo, const float & discricut) const
{
  return (sfb_.getBtagDiscr(jet, algo) >= discricut);
}

// ----------------------------------------------------------------------------
// GetScaledElectrons
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTElectron> Selection::GetScaledElectrons
                                                        (float scale) const
{
  // Container for output
  std::vector<IPHCTree::NTElectron> newElectrons;

  // Get Electrons  
  if (GetPointer2Electrons()!=0) newElectrons = *GetPointer2Electrons();
  else return newElectrons;
     
  for (unsigned int i=0; i<newElectrons.size(); i++)
  {
    if(scale > 1)
    {
      if(newElectrons[i].isEB == 1) scale = 1.005;
      else scale = 1.0275;
    }
    newElectrons[i].p4.SetPxPyPzE(scale*newElectrons[i].p4.Px(),
                                  scale*newElectrons[i].p4.Py(),
                                  scale*newElectrons[i].p4.Pz(),
                                  scale*newElectrons[i].p4.E());
  }
  return newElectrons;
}


// ----------------------------------------------------------------------------
// GetSmearedElectrons
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTElectron> Selection::GetSmearedElectrons
                                                        (float resol) const
{
  // Container for output
  std::vector<IPHCTree::NTElectron> newElectrons;

  // Get Electrons  
  if (GetPointer2Electrons()!=0) newElectrons = *GetPointer2Electrons();
  else return newElectrons;

  // Loop over Electrons
  for (unsigned int i=0; i<newElectrons.size(); i++)
  {
    if(resol > 0)
    {
      if (newElectrons[i].isEB == 1) resol = 0.25;
      else resol = 0.34;
    }
      
    double deltaR = -1;
    double gen_pt = -1; //to be changed

    // Look through W decays
    for(unsigned int igenEl=0; igenEl<mc()->wAndDecays.size(); igenEl++)
    {
      if(abs(mc()->wAndDecays[igenEl].mcLepId) != 11  ) continue;
      deltaR = mc()->wAndDecays[igenEl].p4_Lep_gen.DeltaR(newElectrons[i].p4);
      if(deltaR<0.1) gen_pt = mc()->wAndDecays[igenEl].p4_Lep_gen.Pt();
    }

    // Look through Z decays
    for(unsigned int igenEl=0; igenEl<mc()->zAndDecays.size(); igenEl++)
      {
        if (mc()->zAndDecays[igenEl].Lep1_pdgID == 11 || 
            mc()->zAndDecays[igenEl].Lep1_pdgID == -11)
        {
          deltaR = mc()->zAndDecays[igenEl].p4_Lep1_gen.DeltaR
                                                   (newElectrons[i].p4);
          if(deltaR<0.1) gen_pt = mc()->zAndDecays[igenEl].p4_Lep1_gen.Pt();
        }
        if (mc()->zAndDecays[igenEl].Lep2_pdgID == 11 || 
            mc()->zAndDecays[igenEl].Lep2_pdgID == -11)
        {
          deltaR = mc()->zAndDecays[igenEl].p4_Lep2_gen.DeltaR
                                                    (newElectrons[i].p4);
          if(deltaR<0.1) gen_pt = mc()->zAndDecays[igenEl].p4_Lep2_gen.Pt();
        }
    }
      
    if( gen_pt > 0)
    {
      double ele_pt = newElectrons[i].p4.Pt();
      double deltapt = (ele_pt - gen_pt) * resol;
      double ptscale = ((ele_pt + deltapt) / ele_pt);
      if(ptscale <0 ) ptscale = 0;
      newElectrons[i].p4.SetPxPyPzE(ptscale*newElectrons[i].p4.Px(),
                                    ptscale*newElectrons[i].p4.Py(),
                                    ptscale*newElectrons[i].p4.Pz(),
                                    ptscale*newElectrons[i].p4.E() );
    }
  }

  return newElectrons;
}


// ----------------------------------------------------------------------------
// GetSelectedElectronsNoIso
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTElectron> Selection::GetSelectedElectronsNoIso(
                                           float PtThr, float EtaThr,
                                           bool applyLES, float scale,
                         	                 bool applyLER, float resol) const
{
  // Containers for output
  std::vector<IPHCTree::NTElectron> selectedElectrons;
  std::vector<IPHCTree::NTElectron> localElectrons;
  
  // Get new electrons
  if (applyLES ) { localElectrons = GetScaledElectrons(scale); }
  else { if(!applyLER) localElectrons = *GetPointer2Electrons();}
  
  if(applyLER ) { localElectrons = GetSmearedElectrons(resol); }
  else{ if(!applyLES) localElectrons = *GetPointer2Electrons();}
  
  // Loop over electrons
  for(unsigned int i=0;i<localElectrons.size();i++)
  {
    //bool hadId = localElectrons[i].hadId(
    //    static_cast<unsigned int>(localElectrons[i].ID["simpleEleId90relIso"]) & 0x1
    //                            );
				
    bool hadId = localElectrons[i].hadId(
        static_cast<unsigned int>(localElectrons[i].ID["simpleEleId90cIso"]) & 0x1
                                );
    
    //useless 
    if (!localElectrons[i].isGsfElectron) continue; 
    if (!hadId)                           continue;
   
    // TO DO 
     //   if ( localElectrons[i].nlost >=2 || 
     //        ( fabs(localElectrons[i].deltaCotTheta)<0.02 && 
     //          fabs(localElectrons[i].deltaDistance)<0.02)     )   continue;
    if ( localElectrons[i].nLost >=2 || 
             ( fabs(localElectrons[i].deltaCotTheta)<0.02 && 
               fabs(localElectrons[i].deltaDistance)<0.02)     )   continue;

    if(!localElectrons[i].isEcalDriven)                        continue;
    if(fabs(localElectrons[i].D0)       >=cfg.ElectronD0Cut_)  continue; 
    if(fabs(localElectrons[i].p4.Eta()) >=EtaThr)              continue;
    if(localElectrons[i].p4.Pt()        <=PtThr)               continue;
    if(localElectrons[i].ET_SC          <=cfg.ElectronETSCThr_)continue;
     
    // DR computed between localElectrons candidates
    // and all Global or Tracker Muon
    bool SharedCone = false;
    const std::vector<IPHCTree::NTMuon>* muons = GetPointer2Muons();
    for(unsigned int j=0;j<muons->size();j++)
    {
      if( (*muons)[j].isGlobalMuon || 
          (*muons)[j].isTrackerMuon   )
      {
        if(localElectrons[i].p4.DeltaR((*muons)[j].p4)<0.1)
        {
          SharedCone = true;
          break;
        }
      }
    }
    if(SharedCone)                      continue;
    if(GetSelectedVertex().size() == 0) continue;
    if ( fabs( localElectrons[i].vertex.Z() - 
         GetSelectedVertex()[0].p3.Z() ) > cfg.ElectronVertexMatchThr_ )
        continue;
    selectedElectrons.push_back(localElectrons[i]);
  }
  std::sort(selectedElectrons.begin(),selectedElectrons.end(),HighestPt());
  return selectedElectrons;
}


// ----------------------------------------------------------------------------
// GetSelectedElectrons
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTElectron> Selection::GetSelectedElectrons(
                               float PtThr, float EtaThr,
                               float ElectronRelIso, bool applyLES,
                               float scale, bool applyLER , float resol) const
{
  std::vector<IPHCTree::NTElectron> selectedElectrons;

  std::vector<IPHCTree::NTElectron> electrons = 
              GetSelectedElectronsNoIso(PtThr, EtaThr, applyLES,
                                        scale, applyLER, resol);

  for(unsigned int i=0;i<electrons.size();i++)
  {
    if ( static_cast<double>(electrons[i].RelIso03PF()) 
                                        > ElectronRelIso) continue;
    selectedElectrons.push_back(electrons[i]);
  }
  std::sort(selectedElectrons.begin(),selectedElectrons.end(),HighestPt());
  return selectedElectrons;
}


// ----------------------------------------------------------------------------
// GetSelectedElectronIsoNonId
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTElectron> Selection::GetSelectedElectronsIsoNonID(
                                bool applyLES, float scale,
                                bool applyLER , float resol) const
{
  std::vector<IPHCTree::NTElectron> selectedElectrons;
  std::vector<IPHCTree::NTElectron> localElectrons;
  
  if(applyLES) localElectrons = GetScaledElectrons(scale); 
  else localElectrons = *GetPointer2Electrons();
  
  if(applyLER) localElectrons = GetSmearedElectrons(resol); 
  else localElectrons = *GetPointer2Electrons();
   
  for(unsigned int i=0;i<localElectrons.size();i++)
  {
    if(localElectrons[i].p4.Pt()        < cfg.ElectronPtThreshold_ ) continue;
    if(fabs(localElectrons[i].p4.Eta()) > cfg.ElectronEtaThreshold_) continue;
    if(static_cast<double>(localElectrons[i].RelIso03PF()) 
                                            > cfg.ElectronRelIso_)  continue;
    selectedElectrons.push_back(localElectrons[i]);
  }
  std::sort(selectedElectrons.begin(),selectedElectrons.end(),HighestPt());
  return selectedElectrons;
}

// ----------------------------------------------------------------------------
// GetSelectedElectronNoIsoNonId
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTElectron> Selection::GetSelectedElectronsNoIsoNonID(
                      bool applyLES, float scale, 
		      bool applyLER , float resol ) const
{
  
  std::vector<IPHCTree::NTElectron> selectedElectrons;
  std::vector<IPHCTree::NTElectron> localElectrons;
  
  if(applyLES) localElectrons = GetScaledElectrons(scale); 
  else localElectrons = *GetPointer2Electrons();
  
  if(applyLER) localElectrons = GetSmearedElectrons(resol); 
  else localElectrons = *GetPointer2Electrons();
 
  for(unsigned int i=0;i<localElectrons.size();i++)
  {
    if(localElectrons[i].p4.Pt()        < cfg.ElectronPtThreshold_ ) continue;
    if(fabs(localElectrons[i].p4.Eta()) > cfg.ElectronEtaThreshold_) continue;
    if(!localElectrons[i].isGsfElectron) continue; 
    
    selectedElectrons.push_back(localElectrons[i]);
  }
  std::sort(selectedElectrons.begin(),selectedElectrons.end(),HighestPt());
  return selectedElectrons;
}

// ----------------------------------------------------------------------------
// GetSelectedElectronNoIso
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTElectron> Selection::GetSelectedElectronsNoIso(
                                 bool applyLES, float scale,
                                 bool applyLER , float resol ) const
{
	return GetSelectedElectronsNoIso(cfg.ElectronPtThreshold_,
                                   cfg.ElectronEtaThreshold_,
                                   applyLES, scale);
}


// ----------------------------------------------------------------------------
// GetSelectedElectron
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTElectron> Selection::GetSelectedElectrons(
                              bool applyLES, float scale,
                              bool applyLER, float resol ) const
{
	return GetSelectedElectrons(cfg.ElectronPtThreshold_,
                              cfg.ElectronEtaThreshold_,
                              cfg.ElectronRelIso_,
                              applyLES, scale,
                              applyLER, resol);
}


// ----------------------------------------------------------------------------
// GetSelectedElectronForLJets
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTElectron> Selection::GetSelectedElectronsForLJets(
//                              const std::vector<IPHCTree::NTJet>& SelectedJets,
                              bool applyLES, float scale,
                              bool applyLER, float resol ) const
{
        return GetSelectedElectronsForLJets(cfg.ElectronPtThreshold_,
                              cfg.ElectronEtaThreshold_,
                              cfg.ElectronRelIso_,
                              applyLES, scale,
                              applyLER, resol);
} 


std::vector<IPHCTree::NTElectron> Selection::GetSelectedElectronsNoIsoForLJets(
                               const std::vector<IPHCTree::NTJet>& SelectedJets,
                                           float PtThr, float EtaThr,
                                           bool applyLES, float scale,
                         	                 bool applyLER, float resol) const
{
  // Containers for output
  std::vector<IPHCTree::NTElectron> selectedElectrons;
  std::vector<IPHCTree::NTElectron> localElectrons;
  
  // Get new electrons
  if (applyLES ) { localElectrons = GetScaledElectrons(scale); }
  else { if(!applyLER) localElectrons = *GetPointer2Electrons();}
  
  if(applyLER ) { localElectrons = GetSmearedElectrons(resol); }
  else{ if(!applyLES) localElectrons = *GetPointer2Electrons();}
  
  // Loop over electrons
  for(unsigned int i=0;i<localElectrons.size();i++)
  {
    bool hadId = localElectrons[i].hadId(
        static_cast<unsigned int>(localElectrons[i].ID["cicHyperTight1MC"]) & 0x1   
                                );
  
    if (!localElectrons[i].isGsfElectron) continue; 
    if(localElectrons[i].p4.Pt()        <=PtThr)               continue;
    if(fabs(localElectrons[i].p4.Eta()) >=EtaThr)              continue;
    // missing cut on ETA_SC PATElectron.superCluster()->eta())
    // excluding EB-EE transition region 1.4442 < abs(eta_sc) < 1.5660
    if(fabs(localElectrons[i].D0)       >=cfg.ElectronD0Cut_)  continue; 

    // dR(e,jet)>0.3 where jet is any jet passing the jet requirements 
    /*
    double deltaRel  = 10000;
    for(unsigned int ijet=0; ijet< SelectedJets.size(); ijet++){
      double deltaR = SelectedJets[ijet].p4.DeltaR(localElectrons[i].p4);
      if(deltaR < deltaRel) deltaRel = deltaR;
      if (deltaR<0.3) cout << " ze jet " << SelectedJets[ijet].p4.Pt() << " eta " << SelectedJets[ijet].p4.Eta() << " phi " << SelectedJets[ijet].p4.Phi()
                           << " ze el " << localElectrons[i].p4.Pt() << " eta " << localElectrons[i].p4.Eta() << " phi " << localElectrons[i].p4.Phi() 
                           << " dR " << deltaR << endl;
    }
    if( deltaRel <= 0.3) continue;
    */

    if(localElectrons[i].nLost>0) continue;
    if( fabs(localElectrons[i].deltaCotTheta)<0.02 or
               fabs(localElectrons[i].deltaDistance)<0.02     )   continue;
    if (!hadId)                  continue;

    // plus utilise on dirait...
    //if(!localElectrons[i].isEcalDriven)                        continue;
    //if(localElectrons[i].ET_SC          <=cfg.ElectronETSCThr_)continue;
    //if(GetSelectedVertex().size() == 0) continue;
    //if ( fabs( localElectrons[i].vertex.Z() - 
    //     GetSelectedVertex()[0].p3.Z() ) > cfg.ElectronVertexMatchThr_ )
    //    continue;
     
    // DR computed between localElectrons candidates
    // and all Global or Tracker Muon
    /*
    bool SharedCone = false;
    const std::vector<IPHCTree::NTMuon>* muons = GetPointer2Muons();
    for(unsigned int j=0;j<muons->size();j++)
    {
      if( (*muons)[j].isGlobalMuon || 
          (*muons)[j].isTrackerMuon   )
      {
        if(localElectrons[i].p4.DeltaR((*muons)[j].p4)<0.1)
        {
          SharedCone = true;
          break;
        }
      }
    }
    if(SharedCone)                      continue;  // encore utilise ou pas? top projection???
    */
	
    selectedElectrons.push_back(localElectrons[i]);
  }
  std::sort(selectedElectrons.begin(),selectedElectrons.end(),HighestPt());
  return selectedElectrons;
}

// ----------------------------------------------------------------------------
// GetScaledMuons
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTMuon> Selection::GetScaledMuons(float scale) const
{
  std::vector<IPHCTree::NTMuon> newMuons = *GetPointer2Muons();
  for (unsigned int i=0; i<newMuons.size(); i++)
  {
    newMuons[i].p4.SetPxPyPzE(scale*newMuons[i].p4.Px(),
                              scale*newMuons[i].p4.Py(),
                              scale*newMuons[i].p4.Pz(),
                              scale*newMuons[i].p4.E());
  }
  return newMuons;
}
  

// ----------------------------------------------------------------------------
// GetScaledMuonsNoIso
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTMuon> Selection::GetSelectedMuonsNoIso(
                                 float PtThr, float EtaThr,
                                 bool applyLES, float scale) const
{
  // Containers for output
  std::vector<IPHCTree::NTMuon> selectedMuons;
  vector<IPHCTree::NTMuon> localMuons;

  // Get muons
  if (applyLES) localMuons = GetScaledMuons(scale); 
  else localMuons = *GetPointer2Muons();

  for(unsigned int i=0;i<localMuons.size();i++)
  {
    if (!localMuons[i].isGlobalMuon)  continue; // isGlobalMuon
    if (!localMuons[i].isTrackerMuon) continue; // isTrackerMuon
    if (localMuons[i].Chi2           >= cfg.MuonNormChi2_)       continue;
    if (localMuons[i].NTrValidHits   <= cfg.MuonNofValidTrHits_) continue;
    if (localMuons[i].NValidHits     <= cfg.MuonNofValidHits_  ) continue;
    if (fabs(localMuons[i].D0Inner)  >= cfg.MuonD0Cut_)          continue;
    if (fabs(localMuons[i].p4.Eta()) >= EtaThr)                  continue;
    if (localMuons[i].p4.Pt()        <  PtThr)                   continue;
    if (GetSelectedVertex().size() == 0) continue;
    if ( fabs( localMuons[i].vertex.Z() - 
         GetSelectedVertex()[0].p3.Z() ) > cfg.MuonVertexMatchThr_ ) continue;
    selectedMuons.push_back(localMuons[i]);
  }

  std::sort(selectedMuons.begin(),selectedMuons.end(),HighestPt());
  return selectedMuons;
}


// ----------------------------------------------------------------------------
// GetSelectedMuons
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTMuon> Selection::GetSelectedMuons(
                                  float PtThr, float EtaThr,
                                  float MuonRelIso, bool applyLES,
                                  float scale) const
{
  // Containers for output
  std::vector<IPHCTree::NTMuon> selectedMuons;
  std::vector<IPHCTree::NTMuon> muons = 
                 GetSelectedMuonsNoIso(PtThr,EtaThr,applyLES,scale);

  // Loop over muons 
  for(unsigned int i=0;i<muons.size();i++)
  {
    if ( static_cast<double>(muons[i].RelIso03PF()) > MuonRelIso) continue;
    selectedMuons.push_back(muons[i]);
  }
  std::sort(selectedMuons.begin(),selectedMuons.end(),HighestPt());
  return selectedMuons;
}


// ----------------------------------------------------------------------------
// GetSelectedMuonIsoNonId
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTMuon> Selection::GetSelectedMuonIsoNonID(
                                        bool applyLES,
                                        float scale) const
{
  // Containers for output 
  std::vector<IPHCTree::NTMuon> selectedMuons;
  vector<IPHCTree::NTMuon> localMuons;

  // Get muons
  if(applyLES) localMuons = GetScaledMuons(scale); 
  else localMuons = *GetPointer2Muons();

  // Loop over muons
  for(unsigned int i=0;i<localMuons.size();i++)
  {
    if( localMuons[i].p4.Pt()              < cfg.MuonPtThreshold_)  continue;
    if( fabs(localMuons[i].p4.Eta())       > cfg.MuonEtaThreshold_) continue;
    if((double) localMuons[i].RelIso03PF() > cfg.MuonRelIso_)       continue;
    selectedMuons.push_back(localMuons[i]);
  }

  // Sort muons according to PT rank
  std::sort(selectedMuons.begin(),selectedMuons.end(),HighestPt());
  return selectedMuons;
}


// ----------------------------------------------------------------------------
// GetSelectedMuonsNoIso
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTMuon> Selection::GetSelectedMuonsNoIso(
                                   bool applyLES, float scale) const{
	return GetSelectedMuonsNoIso(cfg.MuonPtThreshold_,
                               cfg.MuonEtaThreshold_,
                               applyLES, scale);
}


// ----------------------------------------------------------------------------
// GetSelectedMuons
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTMuon> Selection::GetSelectedMuons(
                                   bool applyLES, float scale) const
{
	return GetSelectedMuons(cfg.MuonPtThreshold_,
                          cfg.MuonEtaThreshold_,
                          cfg.MuonRelIso_, applyLES, scale);
}


// ----------------------------------------------------------------------------
// GetSelectedMuonsNoIsoForLJets
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTMuon> Selection::GetSelectedMuonsNoIsoForLJets(
                               const std::vector<IPHCTree::NTJet>& SelectedJets,
                                   float PtThr, float EtaThr,
                                   bool applyLES, float scale) const
{
  // originally developped by Camille & Jeremy
  // modified by Caroline for TTbarMET in lepton+jets+met

  // Containers for output 
  std::vector<IPHCTree::NTMuon> selectedMuons;
  std::vector<IPHCTree::NTMuon> localMuons;

  // Get Muons
  if(applyLES) localMuons = GetScaledMuons(scale); 
  else localMuons = *GetPointer2Muons();

  // Loop over muons
  for(unsigned int i=0;i<localMuons.size();i++)
  {

    if(!localMuons[i].isGlobalMuon)  continue;
    if(!localMuons[i].isTrackerMuon) continue;   // pas en 2012
    if(localMuons[i].p4.Pt()        <PtThr)                      continue;
    if(fabs(localMuons[i].p4.Eta()) >= EtaThr)                   continue;
    if(localMuons[i].Chi2           >= cfg.MuonNormChi2_)        continue;
    // if(localMuons[i].NTrValidHits   <= cfg.MuonNofValidTrHits_)  continue;  // remplace trackerLayersWithMeasurement() ?
    // missing cut on muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 8  // 2012 >5
    // missing cut on muon.innerTrack()->hitPattern().pixelLayersWithMeasurement() >=1    // 2012 >0 
    // missing cut on numberOfMatchedStations() >1                                        // 2012 >1
    if(localMuons[i].NValidHits     <= cfg.MuonNofValidHits_  )  continue;
    if(fabs(localMuons[i].D0Inner)  >= cfg.MuonD0Cut_)           continue;
    if(GetSelectedVertex().size() == 0)  continue;
    if ( fabs( localMuons[i].vertex.Z() - 
      GetSelectedVertex()[0].p3.Z()) > cfg.MuonVertexMatchThr_ ) continue;
    // dR(μ,jet)>0.3 where jet is any jet passing the jet requirements
    /*
    double deltaRmu  = 10000;
    for(unsigned int ijet=0; ijet< SelectedJets.size(); ijet++){
      double deltaR = SelectedJets[ijet].p4.DeltaR(localMuons[i].p4);
      if(deltaR < deltaRmu) deltaRmu = deltaR;
      if (deltaR<0.3) cout << " ze jet " << SelectedJets[ijet].p4.Pt() << " eta " << SelectedJets[ijet].p4.Eta() << " phi " << SelectedJets[ijet].p4.Phi()
                           << " ze mu " << localMuons[i].p4.Pt() << " eta " << localMuons[i].p4.Eta() << " phi " << localMuons[i].p4.Phi() 
                           << " dR " << deltaR << endl;
    }
    if( deltaRmu <= 0.3) continue;
    */

    selectedMuons.push_back(localMuons[i]);
  }
  std::sort(selectedMuons.begin(),selectedMuons.end(),HighestPt());
  return selectedMuons;
}

// ----------------------------------------------------------------------------
// GetVetoMuonsForLJets
// ----------------------------------------------------------------------------

std::vector<IPHCTree::NTMuon> Selection::GetVetoMuonsForLJets(
                              bool applyLES, float scale) const
{       

  // for 2012 ana!

  // Containers for output 
  std::vector<IPHCTree::NTMuon> selectedMuons;
  std::vector<IPHCTree::NTMuon> localMuons;

  // Get Muons
  if(applyLES) localMuons = GetScaledMuons(scale); 
  else localMuons = *GetPointer2Muons();

  // Loop over muons
  for(unsigned int i=0;i<localMuons.size();i++)
  {
    if(!localMuons[i].isGlobalMuon && !localMuons[i].isTrackerMuon)  continue;
    if(localMuons[i].p4.Pt()        <10.)                      continue;
    if(fabs(localMuons[i].p4.Eta()) >= 2.5)                   continue;
    if((localMuons[i].RelIso03PF()) > 0.20)                   continue;


/*  // autres cuts aussi ou pas????

    if(localMuons[i].Chi2           >= cfg.MuonNormChi2_)        continue;
    // if(localMuons[i].NTrValidHits   <= cfg.MuonNofValidTrHits_)  continue;  // remplace trackerLayersWithMeasurement() ?
    // missing cut on muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 8
    // missing cut on muon.innerTrack()->hitPattern().pixelLayersWithMeasurement() >=1
    // missing cut on numberOfMatchedStations() >1
    if(localMuons[i].NValidHits     <= cfg.MuonNofValidHits_  )  continue;
    if(fabs(localMuons[i].D0Inner)  >= cfg.MuonD0Cut_)           continue;
    if(GetSelectedVertex().size() == 0)  continue;
    if ( fabs( localMuons[i].vertex.Z() - 
      GetSelectedVertex()[0].p3.Z()) > cfg.MuonVertexMatchThr_ ) continue;
    // dR(μ,jet)>0.3 where jet is any jet passing the jet requirements   ???? top projection???
*/

    selectedMuons.push_back(localMuons[i]);
  }
  std::sort(selectedMuons.begin(),selectedMuons.end(),HighestPt());
  return selectedMuons;


}


std::vector<IPHCTree::NTElectron> Selection::GetVetoElectronsForLJets(
                              bool applyLES, float scale,
                              bool applyLER, float resol ) const
{
  // Containers for output
  std::vector<IPHCTree::NTElectron> selectedElectrons;
  std::vector<IPHCTree::NTElectron> localElectrons;
  
  // Get new electrons
  if (applyLES ) { localElectrons = GetScaledElectrons(scale); }
  else { if(!applyLER) localElectrons = *GetPointer2Electrons();}
  
  if(applyLER ) { localElectrons = GetSmearedElectrons(resol); }
  else{ if(!applyLES) localElectrons = *GetPointer2Electrons();}
  
  // Loop over electrons
  for(unsigned int i=0;i<localElectrons.size();i++)
  {
    bool hadId = localElectrons[i].hadId(
        static_cast<unsigned int>(localElectrons[i].ID["cicLooseMC"]) & 0x1   
                                );
  
    if (!localElectrons[i].isGsfElectron) continue; 
    if (!hadId)                  continue;
    if(localElectrons[i].p4.Pt()        <=15)               continue;
    if(fabs(localElectrons[i].p4.Eta()) >=2.5)              continue;
    // missing cut on ETA_SC PATElectron.superCluster()->eta())
    // excluding EB-EE transition region 1.4442 < abs(eta_sc) < 1.5660
    if ( localElectrons[i].RelIso03PF() > 0.2)              continue;


    // if(fabs(localElectrons[i].D0)       >=cfg.ElectronD0Cut_)  continue; 

    // dR(e,jet)>0.3 where jet is any jet passing the jet requirements  ??? top projection? 

    // DR computed between localElectrons candidates
    // and all Global or Tracker Muon
    bool SharedCone = false;
    const std::vector<IPHCTree::NTMuon>* muons = GetPointer2Muons();
    for(unsigned int j=0;j<muons->size();j++)
    {
      if( (*muons)[j].isGlobalMuon || 
          (*muons)[j].isTrackerMuon   )
      {
        if(localElectrons[i].p4.DeltaR((*muons)[j].p4)<0.1)
        {
          SharedCone = true;
          break;
        }
      }
    }
    if(SharedCone)                      continue;  // encore utilise ou pas? top projection???

    selectedElectrons.push_back(localElectrons[i]);
  }
  std::sort(selectedElectrons.begin(),selectedElectrons.end(),HighestPt());
  return selectedElectrons;


} 
// ----------------------------------------------------------------------------
// GetSelectedLooseMuonForMuJets
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTMuon> Selection::GetSelectedLooseMuonsForMuJets(
                                bool applyLES, float scale) const
{
  // Container for output
  std::vector<IPHCTree::NTMuon> selectedLooseMuons;

  // Get muons 
  vector<IPHCTree::NTMuon> localMuons;
  if(applyLES) localMuons = GetScaledMuons(scale); 
  else localMuons = *GetPointer2Muons();

  // Loop over muons
  for(unsigned int i=0; i<localMuons.size(); i++)
  {
    if ( fabs(localMuons[i].p4.Eta()) < 2.5 &&
         localMuons[i].p4.Pt() > 10         &&
         static_cast<double>(localMuons[i].RelIso03PF()) > 0.2)     
      selectedLooseMuons.push_back(localMuons[i]);
  }
  return selectedLooseMuons;
}


// ----------------------------------------------------------------------------
// GetSelectedLooseElectronsForMuJets
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTElectron> 
Selection::GetSelectedLooseElectronsForMuJets
                                (bool applyLES, float scale) const
{
  // Container for output
  std::vector<IPHCTree::NTElectron> selectedLooseElectrons;

  // Get electrons
  std::vector<IPHCTree::NTElectron> localElectrons;
  if(applyLES) localElectrons = GetScaledElectrons(scale); 
  else localElectrons = *GetPointer2Electrons();

  // Loop over electrons
  for(unsigned int i=0;i<localElectrons.size();i++)
  {
    if ( fabs(localElectrons[i].p4.Eta()) < 2.5 &&
         localElectrons[i].p4.Pt() > 15 &&
         static_cast<double>(localElectrons[i].RelIso03PF()) > 0.2)     
      selectedLooseElectrons.push_back(localElectrons[i]);
  }
  return selectedLooseElectrons;
}


// ----------------------------------------------------------------------------
// GetScaledTaus
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTTau> Selection::GetScaledTaus(float scale) const
{
  std::vector<IPHCTree::NTTau> newTaus = *GetPointer2Taus();
  for (unsigned int i=0; i<newTaus.size(); i++)
  {
    newTaus[i].p4.SetPxPyPzE(scale*newTaus[i].p4.Px(),
                             scale*newTaus[i].p4.Py(),
                             scale*newTaus[i].p4.Pz(),
                             scale*newTaus[i].p4.E() );
  }
  return newTaus;
}
  

// ----------------------------------------------------------------------------
// GetScaledTaus
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTTau> Selection::GetSelectedTaus(
                                float PtThr, float EtaThr, bool applyLES,
                                float scale, int isoLevel, bool antiLep) const
{
  // Contaienr for output
  std::vector<IPHCTree::NTTau> selectedTaus;

  // Get taus
  std::vector<IPHCTree::NTTau> localTaus;
  if(applyLES) localTaus = GetScaledTaus(scale);
  else localTaus = *GetPointer2Taus();

  // Loop over taus   
  for(unsigned int i=0;i<localTaus.size();i++)
  {
  
    if ( localTaus[i].leadTrackPt <= cfg.TauLeadTrkPtCut_ ) continue;
    if ( isoLevel ==0 &&  localTaus[i].ID["byLooseIsolation"]  != 1)  continue;
    if ( isoLevel ==1 &&  localTaus[i].ID["byMediumIsolation"] != 1)  continue;
    if ( isoLevel ==2 &&  localTaus[i].ID["byTightIsolation"]  != 1)  continue;
    if(antiLep)
    {
      if ( localTaus[i].ID["againstElectronTight"] != 1) continue;
      if ( localTaus[i].ID["againstMuonTight"] != 1)     continue;
    }
    if(fabs(localTaus[i].p4.Eta())>=EtaThr)              continue;
    if(fabs(localTaus[i].p4.Eta())<1.566 && fabs(localTaus[i].p4.Eta())>1.4442) continue;
     if(localTaus[i].p4.Pt()<PtThr)                      continue;

    if ( fabs( localTaus[i].vertex.Z() - GetSelectedVertex()[0].p3.Z() )
          > cfg.TauVertexMatchThr_ ) continue;
    if( fabs( localTaus[i].D0)  >= 0.04 )     continue;
     selectedTaus.push_back(localTaus[i]);
  }

  // Sort taus according to PT rank
  std::sort(selectedTaus.begin(),selectedTaus.end(),HighestPt());
  return selectedTaus;
}


// ----------------------------------------------------------------------------
// GetTaus
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTTau> Selection::GetTaus(
                      const std::vector<IPHCTree::NTMuon>& muon_cand,
                      const std::vector<IPHCTree::NTElectron>& elec_cand) const
{
  
	// Container for output  
  std::vector<IPHCTree::NTTau> cleanedTaus;

  // Loop over taus
  const std::vector<IPHCTree::NTTau>* taus = GetPointer2Taus();
  for(unsigned int i=0;i<taus->size();i++)
  {
    //IPHCTree::NTTau tau = taus[i];
    double deltaRmu = 10000;
    double deltaRel = 10000;
     
    for(unsigned int imu=0; imu<muon_cand.size(); imu++)
    {
      double deltaR = (*taus)[i].p4.DeltaR(muon_cand[imu].p4);
      if(deltaR < deltaRmu) deltaRmu = deltaR;
    }
    
    for(unsigned int iel=0; iel<elec_cand.size(); iel++)
    {
      double deltaR = (*taus)[i].p4.DeltaR(elec_cand[iel].p4);
      if(deltaR < deltaRel) deltaRel = deltaR;
    }
    
    if( deltaRmu > 0.4  && deltaRel > 0.4)
                       cleanedTaus.push_back((*taus)[i]); 
  }

  // Sort taus according to PT rank
  std::sort(cleanedTaus.begin(),cleanedTaus.end(),HighestPt());
  return cleanedTaus;
}


// ----------------------------------------------------------------------------
// GetTaus
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTTau> Selection::GetSelectedTaus(
                               bool applyLES, float scale,
                               int isoLevel, bool antiLep) const
{
	return GetSelectedTaus(cfg.TauPtThreshold_, cfg.TauEtaThreshold_,
                         applyLES, scale, isoLevel, antiLep);
}


// ----------------------------------------------------------------------------
// GetSelectedTaus
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTTau> Selection::GetSelectedTaus(
                  const std::vector<IPHCTree::NTMuon>& muon_cand, 
                  const std::vector<IPHCTree::NTElectron>& elec_cand,
                  float PtThr, float EtaThr,
                  bool applyLES, float scale,
                  int isoLevel, bool antiLep) const
{
  // Container for output
  std::vector<IPHCTree::NTTau> cleanTaus;

  // Get taus
  std::vector<IPHCTree::NTTau> selectedTaus = 
            GetSelectedTaus(PtThr,EtaThr,applyLES,scale,isoLevel,antiLep);

  // Loop over taus
  for(unsigned int i=0; i< selectedTaus.size(); i++)
  {
    double deltaRmu = 10000;
    double deltaRel = 10000;
     
    for(unsigned int imu=0; imu< muon_cand.size(); imu++)
    {
      double deltaR = selectedTaus[i].p4.DeltaR(muon_cand[imu].p4);
      if(deltaR < deltaRmu) deltaRmu = deltaR;
    }
     
    for(unsigned int iel=0; iel< elec_cand.size(); iel++)
    {
      double deltaR = selectedTaus[i].p4.DeltaR(elec_cand[iel].p4);
      if(deltaR < deltaRel) deltaRel = deltaR;
    }
     
    if( deltaRmu > 0.4  && deltaRel > 0.4)
                     cleanTaus.push_back(selectedTaus[i]);
  }
  return cleanTaus;
}


// ----------------------------------------------------------------------------
// GetSelectedTaus
// ----------------------------------------------------------------------------
std::vector<IPHCTree::NTTau> Selection::GetSelectedTaus(
                    const std::vector<IPHCTree::NTMuon>& muon_cand,
                    const std::vector<IPHCTree::NTElectron>& elec_cand,
                    bool applyLES, float scale,
                    int isoLevel, bool antiLep) const
{
  return GetSelectedTaus(muon_cand, elec_cand, 
                         cfg.TauPtThreshold_, cfg.TauEtaThreshold_,
                         applyLES, scale, isoLevel, antiLep);
}


// ----------------------------------------------------------------------------
// IsAnEventSelected
// ----------------------------------------------------------------------------
bool Selection::isAnEventSelected(unsigned int nElectrons,
                                  unsigned int nMuons,
                                  unsigned int nTaus,
                                  unsigned int nJets)
{
  std::vector<IPHCTree::NTElectron> electrons = GetSelectedElectrons();
  std::vector<IPHCTree::NTMuon>     muons     = GetSelectedMuons();

	if(electrons.size()                        < nElectrons) return false;
	if(muons.size()                            < nMuons)     return false;
	if(GetSelectedTaus().size()                < nTaus)      return false;
	if(GetSelectedJets(muons,electrons).size() < nJets)      return false;
	return true;
}


// ----------------------------------------------------------------------------
// Weight for B-Tag
// ----------------------------------------------------------------------------
void Selection::InitSFBWeight (int flagb, int methodb,
                               int systb, int btagAlgo_,
                               float btagDiscriCut_, int NofBtagJets_)
{
  flag_btagweight_ = flagb;
  methodb_         = methodb;
  systb_           = systb;
  sfb_.SFBinit (btagAlgo_, btagDiscriCut_, NofBtagJets_);
  if (flagb>0)
  {
    sfb_.LoadInfo("efficacite_btag.root");
    sfb_.LoadInfo2("eff_from_ttmadgraph_summer11_multiwp.root");
    sfb_.InitAlgoAndWP(btagAlgo_, btagDiscriCut_);
  }
}

void Selection::ReInitSFBWeight (int flagb, int methodb,
                                 int systb, int btagAlgo_,
				 float btagDiscriCut_, int NofBtagJets_)
{
  flag_btagweight_ = flagb;
  methodb_ = methodb;
  systb_ = systb;

  sfb_.SFBinit (btagAlgo_, btagDiscriCut_, NofBtagJets_);
  if (flagb>0) {
    //sfb_.LoadInfo();
    //sfb_.LoadInfo2();
    sfb_.InitAlgoAndWP(btagAlgo_, btagDiscriCut_);
  }
}

// ----------------------------------------------------------------------------
// ResetParameters4Bweight
// ----------------------------------------------------------------------------
void Selection::ResetParameters4Bweight (int flagb, int methodb, int systb)
{
  flag_btagweight_ = flagb;
  methodb_         = methodb;
  systb_           = systb;

  if (sfb_.GetHistoSFB()==0)
  {
    sfb_.LoadInfo("efficacite_btag.root");
    sfb_.LoadInfo2("eff_from_ttmadgraph_summer11_multiwp.root");
  }
}


// ----------------------------------------------------------------------------
// generate_flat10_weights
// ----------------------------------------------------------------------------
std::vector<double> Selection::generate_flat10_weights(
                                                  TH1D* data_npu_estimated)
{
  // distri verte
  // see SimGeneral/MixingModule/python/mix_E7TeV_FlatDist10_2011EarlyData_inTimeOnly_cfi.py; copy and paste from there: 
  const double npu_probs[25] = {0.0698146584, 0.0698146584, 0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,
                                0.0630151648,0.0526654164,0.0402754482,0.0292988928,0.0194384503,0.0122016783,0.007207042,0.004003637,0.0020278322,
                                0.0010739954,0.0004595759,0.0002229748,0.0001028162,4.58337152809607E-05 };


  // see SimGeneral/MixingModule/python/mix_E7TeV_Summer_2011_50ns_PoissonOOT.py for Summer11
  // const double npu_probs[25] = {0.0400676665,0.040358009,0.0807116334,0.0924154156,0.0924154156,0.0924154156,0.0924154156,0.0924154156,0.0870356742,0.0767913175,0.0636400516,0.0494955563,0.036223831,0.0249767088,0.0162633216,0.0099919945,0.0058339324,0.0032326433,0.0017151846,0.0008505404,0.0004108859,0.0001905137,0.0000842383,0.000034939,0.0000142801};

  /*
  // distri bleue
  // info taken on ttbar mc (spring11 madgraph - no skim) : /opt/sbg/data/data1/cms/ccollard/CMSSW/fichier_root2011/pu_in_mc.C  
  const double npu_probs[25] = {0.0795275, 0.0761536, 0.069365, 0.0730743, 0.0699924, 0.075654, 0.076423, 0.0684544, 0.0732846, 0.0748051, 0.0669443, 0.0582711, 0.0484596, 0.0348441, 0.0204753, 0.0157868, 0.00964813, 0.00537668, 0.00266835, 0.000791556, 0, 0, 0, 0, 0};
  */

  std::vector<double> result(25);
  double s = 0.0;
  for(int inpu=0; inpu<25; ++inpu)
  {
    double npu_estimated = data_npu_estimated->GetBinContent
      (data_npu_estimated->GetXaxis()->FindBin(inpu));

    if (npu_probs[inpu]>0) result[inpu] = npu_estimated / npu_probs[inpu]; 
    else 
    {
      result[inpu] = 0.;
    }
    s += npu_estimated;
  }

  // normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
  for(int inpu=0; inpu<25; ++inpu) result[inpu] /= s;

  return result;
}


// ----------------------------------------------------------------------------
// generatePUWeight
// ----------------------------------------------------------------------------
void Selection::GeneratePUWeight(std::string filename)
{
	TFile file(filename.c_str(),"READ");
	TH1D* histo = (TH1D*) file.Get("pileup");
	// get TH1D* histo from the root file created by estimatePileupD above
	PUWeights = generate_flat10_weights(histo);
}


// ----------------------------------------------------------------------------
// getPUWeight
// ----------------------------------------------------------------------------
float Selection::GetPUWeight()
{
	// get the number of true pileup events from MC in n_pu
	if(getNpu()<PUWeights.size())
  {
    return PUWeights[getNpu()];
	}
	else
  { 
    // should not happen as we have a weight for
    // all simulated n_pu multiplicities!
    cerr<<"n_pu too big"<<endl;
    return 1.;
	}
}


// ----------------------------------------------------------------------------
// LoadElScaleFactors
// ----------------------------------------------------------------------------
void Selection::LoadElScaleFactors()
{
  string fileName(getenv( "CMSSW_BASE" )+string("/src/NTuple/NTupleAnalysis/macros/data/electronSF.root"));
  std::cout<<"Reading the ElScaleFactors file "<<fileName<<endl;
  fexists(fileName, true);
  TFile *f_Data_El = new TFile(fileName.c_str());
  f_Data_El->cd();
  scaleFactEl = (TH2F*)gROOT->FindObject("thehistSF")->Clone("");
  scaleFactEl->SetDirectory(0);
  f_Data_El->Close();
  delete f_Data_El;
} 


// ----------------------------------------------------------------------------
// LoadMuScaleFactors
// ----------------------------------------------------------------------------
void Selection::LoadMuScaleFactors()
{
  string fileName(getenv( "CMSSW_BASE" )+string("/src/NTuple/NTupleAnalysis/macros/data/muonSF.root"));
  std::cout<<"Reading the MuScaleFactors file "<<fileName<<endl;
  fexists(fileName, true);
  TFile *f_Data_Mu = new TFile(fileName.c_str());
  f_Data_Mu->cd();
  scaleFactMu = (TH2F*)gROOT->FindObject("thehistSF")->Clone("");
  scaleFactMu->SetDirectory(0);
  f_Data_Mu->Close();
  delete f_Data_Mu;
} 


// ----------------------------------------------------------------------------
// InitJESUnc
// ----------------------------------------------------------------------------
void Selection::InitJESUnc()
{
  string fileName(getenv( "CMSSW_BASE" )+string("/src/NTuple/NTupleAnalysis/macros/data/JESUncMC.root"));
  std::cout<<"Reading the JES Uncertainty file "<<fileName<<endl;
  fexists(fileName, true);
  TFile* f1 = new TFile(fileName.c_str());
  f1->cd();
  histo_jesunc_    = (TH2F*) f1->Get("JESUncHisto")->Clone("JESUncHisto") ;
  std::cout<<"Entries: "<< histo_jesunc_->GetEntries()<<std::endl;
  histo_jesunc_->SetDirectory(0); // 
  f1->Close();
  delete f1;
}


