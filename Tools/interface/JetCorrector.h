#ifndef JetCorrector_h
#define JetCorrector_h

// system include files
#include <vector>
#include <set>
#include <string>
#include <TH1F.h>

#include "NTFormat/interface/NTEvent.h"

class JetCorrector {

   public:
      JetCorrector(){ h_JEC=0; }
      ~JetCorrector(){}
      
      // Read root file containing histo with corrections values
      void LoadCorrections();
      // Correct energy of jets of the event
      void ApplyCorrections(IPHCTree::NTEvent *event);

   private:
      // histo containing the corrections
      TH1F* h_JEC; 
};

#endif
