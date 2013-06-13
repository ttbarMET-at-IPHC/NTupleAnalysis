#ifndef MT2_332_CALCULATOR_H
#define MT2_332_CALCULATOR_H

#include "Mt2Calculator.h"

//ROOT HEADERS
#include "TLorentzVector.h"
#include "TVector2.h"

namespace Mt2 {
/**
   Class which knows how to calculate M_T2.
   Please take care when choosing the appropriate function to use for your analysis.
   @author Alan Barr & Chris Lester
   @date 9 Feb 2006 and onwards
*/

  class Mt2_332_Calculator : public Mt2Calculator {
  public:
    /** 
	mt2_332
	
	Calculate M_T2 knowing ptmiss not pvis-transverse-lorentz vec
	
	- in principle this method
	has less information available to it than the method below called
	"mt2_4441" but in practice there is very little difference
	between the performance of the two.
    */
    virtual double mt2_332(const TLorentzVector& visibleA, // 3 d.o.f. 
			   const TLorentzVector& visibleB, // 3 d.o.f.
			   const TVector2& ptmiss,                 // 2 d.o.f.
			   double mInvisible) = 0;
    virtual double mt2_332_Sq(const TLorentzVector& visibleA, // 3 d.o.f. 
			      const TLorentzVector& visibleB, // 3 d.o.f.
			      const TVector2& ptmiss,                 // 2 d.o.f.
			      double mInvisible) = 0;
    Mt2_332_Calculator(const std::string & algName) : Mt2Calculator(algName) {};
  };
 
}

#endif
