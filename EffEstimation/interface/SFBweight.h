// Revision du 13 avril 2012

#ifndef SFBweight_h
#define SFBweight_h

#include "NTFormat/interface/NTEvent.h"
#include "Tools/interface/FileExists.h"

// system include files
#include <memory>
#include <assert.h>
#include <typeinfo>
#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include "TROOT.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TFile.h"



using namespace std;
using namespace IPHCTree;

class SFBweight {

   public:
	SFBweight();
	SFBweight(const SFBweight& sfw);
	SFBweight(int, float, int );
        ~SFBweight();

	void SFBinit(int, float, int );
        void LoadInfo(string fileName);
        void LoadInfo2(string fileName);
	void LoadTTsf(string fileName);
	void InitAlgoAndWP(int algo, float wp);
  /**
   * What you get depends on the info variable (but you never get a weight...) :
   *   0: SF from jet-enriched muon events (same for b and c jets)
   *   1: error on SF from jet-enriched muon events (for c jets, you get the b jet error, so use 2x err)
   *   2: data efficiency for b or light jets from jet-enriched muon events (there is no c-jet eff provided)
   *   3: error on data efficiency for b or light jets from jet-enriched muon events (there is no c-jet eff provided)
   *   4: MC efficiency from tt events
   *   5: error on MC efficiency from tt events
   *   6: SF from tt events (only for b & c jets)
   *   7: error on SF from tt events (only for b & c jets)
   */
        float GetWeight(int info, int, float, float) const;
        vector<float> GetWeigth4BSel(int,  int, const std::vector<NTJet> &  selJets) const;
        vector<float> GetWeigth4BSel(int method_b,  int syst_b, const std::vector<NTJet> &  selJets, float sf_val_for_b, float sf_val_for_l) const;
        const TH2D* GetHistoSFB() const;
	const TH2D* GetHistoEffMCb() const;
	const TH2D* GetHistoEffMCc() const;
	const TH2D* GetHistoEffMCl() const;

  /**
   * returns the flavour of the jet.
   * The flavour will be in the variable quarkorigin
   * sectectedflavour is 1 if the quarkorigin could be determined, 0 if not
   */
        void flavour(const NTJet & jet, int & sectectedflavour, int & quarkorigin) const;

   /**
   *  Returns the b-tagging discriminant of a particular jet according to
   *  the algorithm specified in the config file.
   */
	double getBtagDiscr(const NTJet & jet, const int& algo) const;

  /**
   * returns the data efficiency of a jet
   */
	float efficiency(const NTJet & jet,
	  bool applySFlSyst = false, float factorSFlSyst = 0.,
	  bool applySFbcSyst = false, float factorSFbcSyst = 0.) const;

	int Test() const {return map_effmcc_.size();}

	float getTTsf() const {return tt_sf;} 
	float getTTsfErr() const {return tt_sf_err;}

   private:
        std::string method_origin1_;
        std::string method_origin2_;

        TH2D* histo_sfvalb_;
        TH2D* histo_sferrb_;
        TH2D* histo_sfvall_;
        TH2D* histo_sferrl_;

        TH2D* histo_effvalb_;
        TH2D* histo_efferrb_;
        TH2D* histo_effvall_;
        TH2D* histo_efferrl_;

	//access to the current histo for the algo/wp initialized

        TH2D* histo_effmcb_;
        TH2D* histo_effmcc_;
        TH2D* histo_effmcl_;

        TH2D* histo_errmcb_;
        TH2D* histo_errmcc_;
        TH2D* histo_errmcl_;
 		
	//map containing all the algo/wp

        map<string,TH2D*> map_effmcb_;
        map<string,TH2D*> map_effmcc_;
        map<string,TH2D*> map_effmcl_;

        map<string,TH2D*> map_errmcb_;
        map<string,TH2D*> map_errmcc_;
        map<string,TH2D*> map_errmcl_;
	

        int btag_algo_;
        float btag_discri_;
        int n_bjets_;
	std::string algoname, wp;
	float tt_sf, tt_sf_err;

};

#endif
