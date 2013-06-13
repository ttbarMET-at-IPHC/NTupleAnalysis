#ifndef Mt2bblCom_h
#define Mt2bblCom_h

#include <vector>

#include "EventReco/interface/Mt2Com.h"

//IPHC HEADERS
#include "NTFormat/interface/NTEvent.h"

//ROOT HEADERS
#include <TTree.h>
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TMath.h"

class Mt2bblCom : public Mt2Com
{
	public:

	Mt2bblCom(const std::vector<IPHCTree::NTJet>& VBJet,
						const std::vector<IPHCTree::NTJet>& VJet,
						const std::vector<IPHCTree::NTElectron>& VElectron,
						const std::vector<IPHCTree::NTMuon>& VMuon);

	~Mt2bblCom();

// ----------------------------------------------------------------------------
// Global functions
// ----------------------------------------------------------------------------

	void ComputeMt2bbl();

	double GetMt2bl();
	double GetMt2b();


	private:

	std::vector<IPHCTree::NTJet> VBJet;
	std::vector<IPHCTree::NTJet> VJet;
	std::vector<IPHCTree::NTElectron> VElectron;
	std::vector<IPHCTree::NTMuon> VMuon;

	std::vector<IPHCTree::NTJet> NTBJet;
	std::vector<IPHCTree::NTJet> NTJet;
	std::vector<IPHCTree::NTElectron> NTElectron;
	std::vector<IPHCTree::NTMuon> NTMuon;

	std::vector<double> Mt2blCalc;
	std::vector<double> Mt2bCalc;
	int ctr;
	int it;
};

#endif
