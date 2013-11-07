#ifndef Dataset_h
#define Dataset_h

#include <string>
#include <iostream>
#include <vector>
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TROOT.h"

using namespace std;

class Dataset 
{
public:

    // Constructors

    Dataset() : 
        Name_(string("")), isData_(false), DotIt_(false),
        Color_(1), LineStyle_(1), LineWidth_(2), 
        NormFactor_(1), Xsection_(-1), PreSelEfficiency_(1.), eventTree_(NULL), NofEvtsToRunOver_(-999) { };

    Dataset(string name, bool isData, bool doIt, int color, int lineStyle, int lineWidth, float normFactor, float xsection) : 
        Name_(name), isData_(isData), DotIt_(doIt), 
        Color_(color), LineStyle_(lineStyle),LineWidth_(lineWidth), 
        NormFactor_(normFactor), Xsection_(xsection), PreSelEfficiency_(1.), eventTree_(NULL), NofEvtsToRunOver_(-999){ };

    Dataset(string name, bool isData, bool doIt, 
            int color, int lineStyle, int lineWidth, 
            float normFactor, float xsection, 
            vector<string> filenames): 
        Name_(name), isData_(isData), DotIt_(doIt), 
        Color_(color), LineStyle_(lineStyle), LineWidth_(lineWidth), 
        NormFactor_(normFactor), Xsection_(xsection), PreSelEfficiency_(1.), NofEvtsToRunOver_(-999), Filenames_(filenames)
    {
        //Dataset constructor
        nSkimmedEvent = 0;
        numberOfEventsBeforeMTSkimmer = 0;
        eventTree_ = new TChain("MyModule/Event"); 
        for(unsigned int i=0;i<Filenames_.size();i++)
        {
          eventTree_->AddFile(Filenames_[i].c_str());
          TFile * theTempFile = new TFile(Filenames_[i].c_str());
          theTempFile->cd("MyModule");
          
          TH1F* normHisto = (TH1F*)gROOT->FindObject("theNormHisto");
          nSkimmedEvent += normHisto->GetEntries();

          TH1F* histoFromMTSkimmer = (TH1F*)gROOT->FindObject("numberOfEventsBeforeMTSkimmer");
          numberOfEventsBeforeMTSkimmer += histoFromMTSkimmer->GetBinContent(1);

        }
        if(nSkimmedEvent>0) PreSelEfficiency_ = (float) eventTree_->GetEntries()/nSkimmedEvent;
         NormFactor_ = Xsection_/nSkimmedEvent;
    };

    // Destructor
    ~Dataset(){ 
        Filenames_.clear();
        eventTree_ = 0;
        delete eventTree_; 
    };
    
    string Name()              const { return Name_; };
    void  SetName(string name) {Name_ = name;};
    bool  isData()             const { return isData_; };
    void  SetIsData(const bool& isData) {isData_ = isData;};
    bool  DoIt()               const { return DotIt_; };
    int   Color()              const { return Color_; };
    int   LineStyle()          const { return LineStyle_; };
    int   LineWidth()          const { return LineWidth_; };
    float NormFactor()         const { return NormFactor_; };
    float Xsection()           const { return Xsection_; };
    float PreSelEfficiency()   const { return PreSelEfficiency_; };
    int   NofEvtsToRunOver()   const 
    { 
        if (NofEvtsToRunOver_>0 && eventTree_ && eventTree_->GetEntries()>=NofEvtsToRunOver_ ) 
            return NofEvtsToRunOver_; 
        if(eventTree_) 
            return eventTree_->GetEntries(); 
        return 0;
    };

    TChain* eventTree()        const { return eventTree_; };
    vector<string> Filenames() const { return Filenames_; };

    void SetEquivalentLuminosity(float EqLumi) { if(EqLumi>0) NormFactor_ = 1/EqLumi;};/** will recompute NormFactor = 1/EqLumi*/
    void SetOriginalNumberOfEvents(int NofEvts) { if(NofEvts>0) NormFactor_ = (Xsection_*PreSelEfficiency_)/NofEvts; NofEvtsToRunOver_ = NofEvts;};/** will compute NormFactor = Xsection/TNofEvts */
    void SetPreselEffAndNumberOfPreselEvents(float PreselEff, int NofSEvts)
    { 
        /** will compute NormFactor = Xsection*PreselEff/NofPreselEvts */

        if(NofSEvts>0) NormFactor_ = Xsection_*PreselEff/NofSEvts; 
        NofEvtsToRunOver_ = NofSEvts; 
        PreSelEfficiency_ = PreselEff;
    };
    
    void SetCrossSection(float xsection) { Xsection_ = xsection;};
    void SetCrossSectionError(float xsErrorMinus, float xsErrorPlus){ xsErrorMinus_ = xsErrorMinus ; xsErrorPlus_ = xsErrorPlus;};
    std::pair<float,float> GetCrossSectionError() const {return std::pair<float,float> (xsErrorMinus_,xsErrorPlus_);};
    
    
    int getNSkimmedEvent() { return nSkimmedEvent; };
    int getNumberOfEventsBeforeMTSkimmer() { return numberOfEventsBeforeMTSkimmer; };
    
private:

    string  Name_;
    bool    isData_;
    bool    DotIt_;
    int     Color_;
    int     LineStyle_;
    int     LineWidth_;
    float   NormFactor_;
    float   Xsection_;
    float   PreSelEfficiency_;
    float   xsErrorMinus_; //relative error
    float   xsErrorPlus_; //relative error
    TChain* eventTree_;
    int     NofEvtsToRunOver_;
    vector<string> Filenames_;
    
    int nSkimmedEvent;
    int numberOfEventsBeforeMTSkimmer;
};




#endif


