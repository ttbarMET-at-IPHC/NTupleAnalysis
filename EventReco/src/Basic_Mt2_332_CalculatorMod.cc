#include "EventReco/interface/Basic_Mt2_332_CalculatorMod.h"

/*#include "Mt2/interface/Mt2Units.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnSimplex.h"
#include "Minuit2/MnScan.h"
#include "Minuit2/FunctionMinimum.h"
#include <fstream>*/

using ROOT::Minuit2::MnUserParameters;
using ROOT::Minuit2::MnMigrad;
using ROOT::Minuit2::MnSimplex;
using ROOT::Minuit2::MnScan;
using ROOT::Minuit2::FunctionMinimum;

namespace Mt2
{

	double Basic_Mt2_332_CalculatorMod::Pt(const TVector2& vect) const
	{
		return sqrt(vect.TVector2::Px()*vect.TVector2::Px()+vect.TVector2::Py()*vect.TVector2::Py());
	}

  /*
    
  stealing parts of files in
  
  http://isscvs.cern.ch/cgi-bin/cvsweb.cgi/offline/PhysicsAnalysis/SUSYPhys/SUSYPhysUser/src/?cvsroot=atlas
  
  like
  
  http://isscvs.cern.ch/cgi-bin/cvsweb.cgi/offline/PhysicsAnalysis/SUSYPhys/SUSYPhysUser/src/SusyUserExampleHistTool.cxx?rev=1.4;content-type=text%2Fplain;cvsroot=atlas
  
  and
  
  http://isscvs.cern.ch/cgi-bin/cvsweb.cgi/~checkout~/offline/PhysicsAnalysis/SUSYPhys/SUSYPhysUser/src/mT2Fcn.cxx?rev=1.1;content-type=text%2Fplain;cvsroot=atlas
  
  BUT CRUCIALLY:

  then modifying them to remove the assumption about massless visible particles!
  */



  double Basic_Mt2_332_CalculatorMod::mt2_332(const TLorentzVector& visA, 
																				   const TLorentzVector& visB,
																					 const TVector2& ptmiss, 
																					 const double mEachInvisible) 
	{
    const double mSq = mt2_332_Sq(visA, visB, ptmiss, mEachInvisible);

    if (mSq>=0) {
      return sqrt(mSq);
    } else {
      // something went wrong -- maybe just rounding errors;
      return sqrt(fabs(mSq));
    }
  }


  double Basic_Mt2_332_CalculatorMod::mt2_332_Sq(const TLorentzVector& visA,
																							const TLorentzVector& visB,
																							const TVector2& ptmiss,
																							const double mEachInvisible)
	{
    mT2Fcn theFCN(ptmiss.Px(),
									ptmiss.Py(),
									mEachInvisible,
									visA.Px(),
									visA.Py(),
									visA.M(),
									visB.Px(),
									visB.Py(),
									visB.M());

    const double massScale = (
			      //ptmiss.Pt() +
			      Pt(ptmiss) +
			      mEachInvisible +
			      visA.Pt()  +
			      visA.M() +
			      visB.Pt() +
			      visB.M()
			      )/6.0;
    // DANG! Try to get rid of Minuit output:
    //std::ofstream    DANG_log("/dev/null");
    //std::streambuf * DANG_save = std::cerr.rdbuf();
    //if (DANG_log.is_open()) {
    //  std::cerr.rdbuf(DANG_log.rdbuf());
    // }

    double guessx = 0.5*(ptmiss.Px());
    double guessy = 0.5*(ptmiss.Py());

		//MnUserParameters est definit dans Minuit2   
    MnUserParameters upar;
    upar.Add("etx", guessx, 0.02*massScale); 
    upar.Add("ety", guessy, 0.02*massScale);
    const int highQuality=2;    

    // Usually migrad produces the best minumum.
    // But when the minimum is in a fold, migrad can fail badly.
    // On the fold, simplex does well.  We therefore do both separately
    // and record the answer of the one that did best.
    
    // Further to the above notes, it now seems that by choosing the massScale sensibly, and by making the "tolerance" (the second argument to the call that extracts the FunctionMinimum below) much smaller, the simplex algorithm seems to work so well that we don't need the migrad algorithm.  This is good news, as the migrad algorithm produces lots of error output that we can't get rid of.

    MnSimplex simplex(theFCN, upar, highQuality);
    FunctionMinimum minS = simplex(0,massScale*0.000001);
    //const double etxAtMin = minS.UserState().Value("etx");
    //const double etyAtMin = minS.UserState().Value("ety");
    //MnMigrad migrad(theFCN, upar, highQuality);
    //FunctionMinimum minM = migrad(0,massScale*0.000001);
    //const double best = fmin(minS.Fval(), minM.Fval());
    const double best = minS.Fval();

    // DANG! Undoing our attempt to get rid of Minuit output:
    //if (DANG_log.is_open()) {
    //  std::cerr.rdbuf(DANG_save);
    //}

    return best;    
  }
  
  double  Basic_Mt2_332_CalculatorMod::mT2Fcn::operator()(const std::vector<double>& par) const {
    
    double qT1x = par[0];
    double qT1y = par[1];
    
    double qT2x = theExmiss - qT1x;
    double qT2y = theEymiss - qT1y;
    
    double ETj1 = sqrt(theMass1*theMass1 +  thePT1x*thePT1x + thePT1y*thePT1y);    // This is where (relative to SUSYPhys_Mt2_222_Calculator::mT2Fcn::operator()) we have removed the assumption that the visible particles are massless!
    double ETj2 = sqrt(theMass2*theMass2 +  thePT2x*thePT2x + thePT2y*thePT2y);     // This is where (relative to SUSYPhys_Mt2_222_Calculator::mT2Fcn::operator()) we have removed the assumption that the visible particles are massless! 
    
    double ETchi1 = sqrt(qT1x*qT1x + qT1y*qT1y + theMchi*theMchi);
    double ETchi2 = sqrt(qT2x*qT2x + qT2y*qT2y + theMchi*theMchi);
    
    double mTsq1 
      = theMass1*theMass1 // This is where (relative to SUSYPhys_Mt2_222_Calculator::mT2Fcn::operator()) we have removed the assumption that the visible particles are massless!
      + theMchi*theMchi + 2.0*(ETj1*ETchi1 - (thePT1x*qT1x + thePT1y*qT1y));
    double mTsq2 
      = theMass2*theMass2  // This is where (relative to SUSYPhys_Mt2_222_Calculator::mT2Fcn::operator()) we have removed the assumption that the visible particles are massless!
      + theMchi*theMchi + 2.0*(ETj2*ETchi2 - (thePT2x*qT2x + thePT2y*qT2y));
   

    // std::cout << "MOO " << sqrt(fmax(mTsq1,mTsq2)) << "\t" << qT1x << "\t" << qT1y << "\t" << 0.5*(qT1x-qT2x) << "\t" << 0.5*(qT1y-qT2y ) << std::endl;
 
    return fmax(mTsq1,mTsq2);
    
  }

} // end of Mt2 Namespace
