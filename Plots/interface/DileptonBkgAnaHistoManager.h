#ifndef DileptonBkgAnaHistoManager_h
#define DileptonBkgAnaHistoManager_h

#include "NTFormat/interface/NTElectron.h"
#include "Plots/interface/HistoManager.h"


using namespace IPHCTree;

class DileptonBkgAnaHistoManager: public HistoManager{

  private:
	std::vector<std::string> histoNames;
	std::vector<std::string> histo2DNames;

  public:  
	DileptonBkgAnaHistoManager();
	~DileptonBkgAnaHistoManager();

	//Initialisation methods
	void CreateHistos(); /** Create a bunch of standard histos */

	//Fill methods
	void Fill(std::vector<std::string> dataName, std::vector<float> dataValue, const int& maxSelStep, const int& iChannel, const int& iDataset, const float& weight);
	void FillSelStep(std::vector<std::string> dataName, std::vector<float> dataValue,  const int& iSelStep, const int& iChannel, const int& iDataset, const float& weight);

	void Fill2D(std::vector<std::string> dataName, std::vector<float> dataValueX, std::vector<float> dataValueY, const int& maxSelStep, const int& iChannel, const int& iDataset, const float& weight);
	void Fill2DSelStep(std::vector<std::string> dataName, std::vector<float> dataValueX, std::vector<float> dataValueY, const int& iSelStep, const int& iChannel, const int& iDataset, const float& weight);


	void AddHisto_(string name, string title, string xaxis, const int& nbins, const float& min, const float& max);
	void AddHisto2D_(string name, string title, string xaxis, const int& nxbins, const float& xmin, const float& xmax, string yaxis, const int& nybins, const float& ymin, const float& ymax);

};

#endif
