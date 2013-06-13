//Corrig√e par Thomas DESCHLER.

#include "EventReco/interface/Resolution.h"

Resolution::Resolution() {}

Resolution::~Resolution() {}

// ----------------------------------------------------------------------------
// Global functions
// ----------------------------------------------------------------------------

//Calcul du sigma sur la masse du jet
double Resolution::GetSigmaPtJet(const TLorentzVector& jet)
{
  sigma=sigma_void;
  TF1 sigma("sigma","sqrt(((TMath::Sign(1,[0])*sq([0]/x))+(sq([1])*(x^([3]-1))))+sq([2]))",0,10);

  if (jet.Eta() > -9.9 && jet.Eta() <= -4)
  {
    sigma.SetParameter(0,1.41584);
    sigma.SetParameter(1,0.209477);
    sigma.SetParameter(2,0.588872);
  }
  if (jet.Eta() > -4 && jet.Eta() <= -3.5)
  {
    sigma.SetParameter(0,1.65966);
    sigma.SetParameter(1,0.223683);
    sigma.SetParameter(2,0.60873);
  }
  if (jet.Eta() > -3.5 && jet.Eta() <= -3)
  {
    sigma.SetParameter(0,2.81978);
    sigma.SetParameter(1,0.272373);
    sigma.SetParameter(2,0.579396);
  }
  if (jet.Eta() > -3 && jet.Eta() <= -2.5)
  {
    sigma.SetParameter(0,2.56933);
    sigma.SetParameter(1,0.305802);
    sigma.SetParameter(2,0.398929);
  }
  if (jet.Eta() > -2.5 && jet.Eta() <= -2)
  {
    sigma.SetParameter(0,1.04792);
    sigma.SetParameter(1,0.466763);
    sigma.SetParameter(2,0.193137);
  }
  if (jet.Eta() > -2 && jet.Eta() <= -1.5)
  {
    sigma.SetParameter(0,-1.12329);
    sigma.SetParameter(1,0.657891);
    sigma.SetParameter(2,0.139595);
  }
  if (jet.Eta() > -1.5 && jet.Eta() <= -1)
  {
    sigma.SetParameter(0,-0.561649);
    sigma.SetParameter(1,0.420293);
    sigma.SetParameter(2,0.392398);
  }
  if (jet.Eta() > -1 && jet.Eta() <= -0.5)
  {
    sigma.SetParameter(0,-0.499735);
    sigma.SetParameter(1,0.336391);
    sigma.SetParameter(2,0.430689);
  }
  if (jet.Eta() > -0.5 && jet.Eta() <= 0.5)
  {
    sigma.SetParameter(0,-0.349206);
    sigma.SetParameter(1,0.297831);
    sigma.SetParameter(2,0.471121);
  }
  if (jet.Eta() > 0.5 && jet.Eta() <= 1)
  {
    sigma.SetParameter(0,-0.499735);
    sigma.SetParameter(1,0.336391);
    sigma.SetParameter(2,0.430689);
  }
  if (jet.Eta() > 1 && jet.Eta() <= 1.5)
  {
    sigma.SetParameter(0,-0.561649);
    sigma.SetParameter(1,0.420293);
    sigma.SetParameter(2,0.392398);
  }
  if (jet.Eta() > 4 && jet.Eta() <= 9.9)
  {
    sigma.SetParameter(0,1.41584);
    sigma.SetParameter(1,0.209477);
    sigma.SetParameter(2,0.588872);
  }
  if (jet.Eta() > 3.5 && jet.Eta() <= 4)
  {
    sigma.SetParameter(0,1.65966);
    sigma.SetParameter(1,0.223683);
    sigma.SetParameter(2,0.60873);
  }
  if (jet.Eta() > 3 && jet.Eta() <= 3.5)
  {
    sigma.SetParameter(0,2.81978);
    sigma.SetParameter(1,0.272373);
    sigma.SetParameter(2,0.579396);
  }
  if (jet.Eta() > 2.5 && jet.Eta() <= 3)
  {
    sigma.SetParameter(0,2.56933);
    sigma.SetParameter(1,0.305802);
    sigma.SetParameter(2,0.398929);
  }
  if (jet.Eta() > 2 && jet.Eta() <= 2.5)
  {
    sigma.SetParameter(0,1.04792);
    sigma.SetParameter(1,0.466763);
    sigma.SetParameter(2,0.193137);
  }
  if (jet.Eta() > 1.5 && jet.Eta() <= 2)
  {
    sigma.SetParameter(0,-1.12329);
    sigma.SetParameter(1,0.657891);
    sigma.SetParameter(2,0.139595);
  }
  
  return sigma.Eval(jet.Pt());
}

// Get vector<btag_discri>
//Les probl√®mes de memory leaks ne viennent pas d'ici.
std::vector<float> Resolution::GetVecSigmaPtJets(const std::vector<IPHCTree::NTJet>& jets) 
{
  vecSigmaPtJets.clear();

  for(unsigned int i=0; i<jets.size(); i++)
    vecSigmaPtJets.push_back(Resolution::GetSigmaPtJet(jets[i].p4));
  
  return vecSigmaPtJets;
}

double Resolution::GetChi2(const std::vector<IPHCTree::NTJet>& jets)
{
  v_i.clear();
  v_j.clear();
  v_k1.clear();
  v_k2.clear();
  btag.clear();

  sigma_jets = Resolution::GetVecSigmaPtJets(jets); // !!! CAUSE LES MEMORY LEAKS !!!

  for (unsigned int i = 0; i < jets.size(); i++)
    btag.push_back(jets[i].bTag["combinedSecondaryVertexBJetTags"]);


  static const float BTAG_MED = 0.679;
  static const float PDG_TOP_MASS = 173.5;
  static const float PDG_W_MASS = 80.385;
  assert(jets.size() == sigma_jets.size());

  //check at most first 6 jets
  n_jets = jets.size();
  if (n_jets>6) n_jets = 6;
  //consider at least 3 jets
  if (n_jets<3) return 99999.;

  for ( int i=0; i<n_jets; ++i )
  for ( int j=i+1; j<n_jets; ++j )
  {
      minimizer=minimizer_void;

      //  W
      hadW = jets[i].p4 + jets[j].p4;

      //  W Mass Constraint.
      p1 = -1;

      minimizer.ExecuteCommand("SET PRINTOUT", &p1, 1);
      minimizer.SetFCN(&(Resolution::minuitFunction));
      minimizer.SetParameter(0 , "c1"     , 1.1             , 1 , 0 , 0);
      minimizer.SetParameter(1 , "pt1"    , 1.0             , 1 , 0 , 0);
      minimizer.SetParameter(2 , "sigma1" , sigma_jets[i]   , 1 , 0 , 0);
      minimizer.SetParameter(3 , "pt2"    , 1.0             , 1 , 0 , 0);
      minimizer.SetParameter(4 , "sigma2" , sigma_jets[j]   , 1 , 0 , 0);
      minimizer.SetParameter(5 , "m12"    , jets[i].p4.M2() , 1 , 0 , 0);
      minimizer.SetParameter(6 , "m22"    , jets[j].p4.M2() , 1 , 0 , 0);
      minimizer.SetParameter(7 , "m02"    , hadW.M2()    , 1 , 0 , 0);

      for (unsigned int k = 1; k < 8; k++)
        minimizer.FixParameter(k);

      minimizer.ExecuteCommand("SIMPLEX", 0, 0);
      minimizer.ExecuteCommand("MIGRAD", 0, 0);

      c1 = minimizer.GetParameter(0);
      if (c1!=c1) {
        cout<<"[PartonCombinatorics::recoHadronicTop] ERROR: c1 parameter is NAN! Skipping this parton combination"<<endl;
        continue;
      }

      c2 = Resolution::fc2(c1, jets[i].p4.M2(), jets[j].p4.M2(), hadW.M2());

      //     * W Mass check :)
      //     *  Never trust a computer you can't throw out of a window.
      //     *  - Steve Wozniak
      //    
      //     *  Never trust somebody you can't throw out of a window.
      //     *  - Me

      v_i.push_back(i);
      v_j.push_back(j);
      v_k1.push_back(c1);
      v_k2.push_back(c2);
  }

  minimizer=minimizer_void;

  //Apply b-consistency requirement
  n_btag = 0;
  for(int i=0;i<n_jets;i++)
    if(btag.at(i)>BTAG_MED) n_btag++;

  chi2min = 99999.;

  //consider b-jet in leading 3 jets
  for(int b=0;b<n_jets;++b)
  {    
    //if not tagged, consider only 3 leading jets
    if( btag.at(b)<BTAG_MED && b>2 ) continue;

    //require b-tagging if have more than 1 b-tag
    if( n_btag>1 && btag.at(b) < BTAG_MED ) continue;

    pt_b = jets[b].p4.Pt();
      
    for (unsigned int w=0;w<v_i.size();++w)
    {
      k = v_i[w];
      l = v_j[w];

      if ( k==b || l==b ) continue;

      //count number of b-tagged Ws
      nwb = 0;
      if(k<btag.size())
        if (btag.at(k) > BTAG_MED) nwb++;

      if(l<btag.size())
        if (btag.at(l) > BTAG_MED) nwb++;

      //no btagged jets in W if have few btags
      if ( n_btag<3  && nwb>0 ) continue;

      //In 3 b-tag case, allow for 1 W jet to be tagged
      // If have more b-tags then btagging information not useful
      if ( n_btag==3 && nwb>1 ) continue;
 
      pt_w1 = jets[k].p4.Pt();
      pt_w2 = jets[l].p4.Pt();
      
      //  W Mass.
      hadW = jets[k].p4 + jets[l].p4;
      massW = hadW.M();
      
      c1 = v_k1[w];
      c2 = v_k2[w];
      
      // Top Mass.
      hadT = (jets[k].p4 * c1) + (jets[l].p4 * c2) + jets[b].p4;
      massT = hadT.M();
      pt_w = hadW.Pt();
      sigma_w2 = pow(pt_w1*sigma_jets[k], 2)+ pow(pt_w2*sigma_jets[l], 2);
      smw2 = (1. + 2.*pow(pt_w,2)/pow(massW,2))*sigma_w2;
      pt_t = hadT.Pt();
      sigma_t2 = pow(c1*pt_w1*sigma_jets[k],2)+ pow(c2*pt_w2*sigma_jets[l],2)+ pow(pt_b*sigma_jets[b],2);
      smtop2 = (1. + 2.*pow(pt_t,2)/pow(massT,2))*sigma_t2;    
      c_chi2 = pow(massT-PDG_TOP_MASS, 2)/smtop2+ pow(massW-PDG_W_MASS, 2)/smw2;

      if (c_chi2<chi2min)
      {
      chi2min = c_chi2;
      hadT_ = hadT;
      smtop_ = sqrt(smtop2);
      mtop_ = massT;
      }
    }
  }
  return chi2min;
}
/*
double Resolution::GetSelChi2(const TTbarMetSelection* sel)
{
  mcInfo=*(sel->GetPointer2MC());
  MCParticles=mcInfo.genParticles;
  TopId.clear();
  TopHad.clear();
  Chi2topHad=-9999;

  for (int i=0;i<MCParticles.size();i++)      
    if ((abs(MCParticles[i].id) == 1)
      ||(abs(MCParticles[i].id) == 2)
      ||(abs(MCParticles[i].id) == 3)
      ||(abs(MCParticles[i].id) == 4)
      ||(abs(MCParticles[i].id) == 5))
      if(MCParticles[i].motherIndex_!=-1)
        if(abs(MCParticles[MCParticles[i].motherIndex_].id)==24)
          if(MCParticles[MCParticles[i].motherIndex_].motherIndex_!=-1)
          {
            TopId.push_back(MCParticles[MCParticles[i].motherIndex_].motherIndex_);
            Chi2=Resolution::GetChi2(sel->GetJetsForAna());
          }

  for(int i=0;i<TopId.size();i++)
    TopHad.push_back(MCParticles[TopId[i]].p4);

  if(TopHad.size()!=0)
    for(int i=0;i<TopHad.size();i++)
      if(TopHad[i].DeltaR(Resolution::GetTLVRecoTop())<0.3) Chi2topHad=Chi2;

  return Chi2topHad;
}
*/
double Resolution::fc2 (double c1, double m12, double m22, double m02)
{
  static const float PDG_W_MASS = 80.385;
  double a = m22;
  double b = (m02 - m12 - m22) * c1;
  double c = m12 * c1 * c1 - PDG_W_MASS * PDG_W_MASS;
  double num = -1. * b + sqrt(b * b - 4 * a * c);
  double den = 2 * a;

  return (num/den);
}

double Resolution::fchi2 (double c1, double pt1, double sigma1, double pt2, double
sigma2,double m12, double m22, double m02)
{
   double rat1 = pt1 * (1 - c1) / sigma1;
   double rat2 = pt2 * (1 - Resolution::fc2(c1, m12, m22, m02)) / sigma2;

   return ( rat1 * rat1 + rat2 * rat2);
}

//--------------------------------------------------------------------
void Resolution::minuitFunction(int& titi, double* tata, double &result, double* par, int toto)
{
   result=Resolution::fchi2(par[0], par[1], par[2], par[3], par[4], par[5], par[6],par[7]);
}

