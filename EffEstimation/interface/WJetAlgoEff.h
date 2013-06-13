#ifndef WJetAlgoEff_h
#define WJetAlgoEff_h

#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"

#include <string>
#include <vector>

#include "NTFormat/interface/NTEvent.h"
#include "Selection/interface/StopMCinfo.h"

class WJetAlgoEff{

  public:
	WJetAlgoEff();
	~WJetAlgoEff();

	// Creating histograms
	void InitializeWTag(std::vector<std::string> algoNames, int nbins = 10, float xmin = 0, float xmax = 500);
	void InitializeTopTag(std::vector<std::string> algoNames, int nbins = 10, float xmin = 0, float xmax = 500);
	
	// Loading event
	void LoadEvent(IPHCTree::NTEvent* event, StopMCinfo* mcInfo);

	// Computing eff & fake-rate
	void Compute();

	// Writing into TFile the histograms
	void Write(TFile* fout);

  private:
 
    // -----------------------
    // Consider W-tagging
    // -----------------------
    
    //list of algorithms
    std::vector<std::string> WalgoNames_;
    //Charateristics of the plots
    int nbinsW_;
    float xminW_;
    float xmaxW_;
    //Histograms of normalisation
    TH1F* WJetAlgoEffNorm;
    std::vector<TH1F*> WJetAlgoFakeNorm;
    //Histograms of eff and fake rate for all algo
    std::vector<TH1F*> WJetAlgoEffPlots_;
    std::vector<TH1F*> WJetAlgoFakePlots_;

    //For display
    TCanvas* cWJetAlgoEff_;
    TCanvas* cWJetAlgoFake_;
    TLegend* legWEff_;
    TLegend* legWFake_;
    
    // -----------------------
    // Consider Top-tagging
    // -----------------------
    
    //list of algorithms
    std::vector<std::string> TopalgoNames_;
    //Charateristics of the plots
    int nbinsTop_;
    float xminTop_;
    float xmaxTop_;
    //Histograms of normalisation
    TH1F* TopJetAlgoEffNorm;
    std::vector<TH1F*> TopJetAlgoFakeNorm;
    //Histograms of eff and fake rate for all algo
    std::vector<TH1F*> TopJetAlgoEffPlots_;
    std::vector<TH1F*> TopJetAlgoFakePlots_;

    //For display
    TCanvas* cTopJetAlgoEff_;
    TCanvas* cTopJetAlgoFake_;
    TLegend* legTopEff_;
    TLegend* legTopFake_;

};

#endif
