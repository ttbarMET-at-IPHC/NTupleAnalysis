#ifndef Resolution_h
#define Resolution_h

// IPHC headers
#include "NTFormat/interface/NTEvent.h"
#include "NTFormat/interface/NTGenParticle.h"

// ROOT headers
#include <TTree.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TMath.h>
#include <TF1.h>
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnSimplex.h"
#include "Minuit2/MnScan.h"
#include "Minuit2/FunctionMinimum.h"
#include "TFitterMinuit.h"

// STL headers
#include <iostream>
#include <memory>
#include <vector>
#include <assert.h>
#include <math.h>

class Resolution
{
    
  // -------------------------------------------------------------
  //                       method members
  // -------------------------------------------------------------
 
 public :

    //Constructor
    Resolution();

    //Destructor
    ~Resolution();

    static double GetSigmaJet(const TLorentzVector& jet);
    static double GetSigmaScaleFactor(double eta);
    double GetChi2(const std::vector<IPHCTree::NTJet>& jets, bool runningOnMonteCarlo);

    static double fc2 (double c1, double m12, double m22, double m02);

    static double fchi2 (double c1, double pt1, double sigma1, double pt2,
                         double sigma2, double m12, double m22, double m02);

    static void minuitFunction(int&, double* , double &result, double* par, int);

    // Sigma of reconstructed top
    const double GetSigmaRecoTop()
    {   return smtop_;  }

    // Mass of reconstructed top
    const double GetMRecoTop()
    {   return mtop_;   }    
    
    // LorentzVector of reconstructed top
    const TLorentzVector GetTLVRecoTop()
    {   return hadT_;   }
    

 private :
    
    double smtop_;
    double mtop_;
    TLorentzVector hadT_;

};
#endif
