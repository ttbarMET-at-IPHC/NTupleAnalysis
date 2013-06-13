#include "Tools/interface/JetCorrector.h"
#include "Tools/interface/FileExists.h"
#include <TFile.h>

// Read root file containing histo with corrections values
void JetCorrector::LoadCorrections(){

  // open file
  std::string fileName(getenv( "CMSSW_BASE" )+std::string("/src/NTuple/NTupleAnalysis/macros/data/GR_R_42_V23_JEC_L2L3Residual.root"));
  std::cout<<"Reading the JEC_L2L3Residual file "<<fileName<<endl;
  fexists(fileName, true);
  TFile *f_JEC = new TFile(fileName.c_str());
  f_JEC->cd();
  
  // get histo
  h_JEC = dynamic_cast<TH1F*>(gROOT->FindObject("AK5PFchs_L2L3Residual")->Clone(""));
  if(!h_JEC) std::cerr<<">>> No AK5PFchs_L2L3Residual correction available <<<"<<std::endl;
  else h_JEC->SetDirectory(0);
  
  f_JEC->Close();
  delete f_JEC;

}

// Correct energy of jets of the event
void JetCorrector::ApplyCorrections(IPHCTree::NTEvent *event){

  std::set<std::string> names;
  event->jets.GetCollectionList(names);
  
  // Loop over jets collections
  for (std::set<std::string>::const_iterator it = names.begin();
       it != names.end(); it++)
  if((*it)=="pf")
  {
    // select the 'pf' collection
    event->jets.SelectLabel(*it);
    // Loop over jets
    for(unsigned int i=0; i<event->jets.size(); i++)
    {
      // get correction factor that depends of eta
      float eta = event->jets[i].p4.Eta();
      float scale = 1;
      if(h_JEC) scale = h_JEC->GetBinContent(h_JEC->FindBin(eta));
      
      // apply correction
      event->jets[i].p4.SetPxPyPzE(scale*event->jets[i].p4.Px(),
                                   scale*event->jets[i].p4.Py(),
                                   scale*event->jets[i].p4.Pz(),
                                   scale*event->jets[i].p4.E());
    } 
  }

}

