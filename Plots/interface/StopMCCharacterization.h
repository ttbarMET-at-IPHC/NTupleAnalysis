#ifndef StopMCCharacterization_h
#define StopMCCharacterization_h

#include "TFile.h"
#include "TH3F.h"

#include "Selection/interface/StopMCinfo.h"

class StopMCCharacterization
{

  public:
    StopMCCharacterization();
    ~StopMCCharacterization();
  
  void Initialize(int nbinx, float xmin, float xmax, int nbiny, float ymin, float ymax);
  void Fill(StopMCinfo* stopMCinfo);
  void Compute();
  void Write(TFile* fout);

  private:

    //Pt distribution
    TH3F* hPt_Stop;
    TH3F* hPt_Neutralino;
    TH3F* hPt_Top;
    TH3F* hPt_W;

};

#endif
