// STL headers
#include <iomanip>
#include <iostream>
#include <time.h>

// IPHC Framework headers
#include "../../../../IPHCDataFormat/NTFormat/interface/NTEvent.h"
#include "../../../../IPHCDataFormat/NTFormat/interface/NTTransient.h"

// ROOT headers
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TDirectory.h>


// -----------------------------------------------------------------------------
//                            SKIMMING FUNCTION
// -----------------------------------------------------------------------------
bool NLeptonSkim(IPHCTree::NTEvent * event, int Nlep)
{
  // Initializing counters
  unsigned short nElectrons = 0;
  unsigned short nMuons     = 0;
    
  // Get electron collection names
  std::set<std::string> labels;
  event->electrons.GetCollectionList(labels);

  // Loop over electron collections
  for (std::set<std::string>::const_iterator theLabel = labels.begin();
       theLabel != labels.end(); theLabel++)
  {
    // Selecting a given collection
    event->electrons.SelectLabel(*theLabel);

    // Saving the number of electrons of the most populated collection
    if(event->electrons.size()>nElectrons) nElectrons = event->electrons.size();      
  }

  // Get muon collection names
  labels.clear();
  event->muons.GetCollectionList(labels);

  // Loop over muon collections
  for (std::set<std::string>::const_iterator theLabel = labels.begin();
       theLabel != labels.end(); theLabel++)
  {
    // Selecting a given collection
    event->muons.SelectLabel(*theLabel);

    // Saving the number of muons of the most populated collection
    if(event->muons.size()>nMuons) nMuons = event->muons.size();      
  }

  // Cut on the number of (light) leptons
  if ((nElectrons+nMuons) < Nlep) return false;
  else return true;
}


// -----------------------------------------------------------------------------
//                               MAIN PROGRAM
// -----------------------------------------------------------------------------
int main (int argc, char *argv[])
{
  // Displaying header
  std::cout << "#########################" << std::endl;
  std::cout << "Beginning of the program"  << std::endl;
  std::cout << "#########################" << std::endl;
  
  // Checking number of arguments
  if(argc<2) 
  { 
    std::cout << " Syntax Error !!" << std::endl;  
    std::cout << "Syntax is : NTupleSkimmer filename"  << std::endl;
    return 1;
  }

  // Extracting the filename
  std::string filename = argv[1];

  // Declaring a pointer to the current event
  IPHCTree::NTEvent * event = 0;

  //----------------------------------------------------------------------------
  //                                  INPUT
  //----------------------------------------------------------------------------

  // Opening input file (read-only mode)
  TFile* f = TFile::Open(filename.c_str());
  if (!f->IsOpen()) 
  {
    std::cout << "Sorry but cannot find the file '" << filename << "'" << std::endl;
    return 0;
  }

  // Accessing the input Tree
  TTree * eventTree = dynamic_cast<TTree*>(f->Get("MyModule/Event"));
  if (eventTree==0) 
  {
    std::cout << "Sorry but cannot find the tree 'MyModule/Event'" << std::endl;
    return 0;
  }

  // Linking the Tree with the event
  eventTree->SetBranchAddress ("NTEvent", &event);

  // Accessing to the histo 'theNormHisto'
  TH1F* theNormHisto        = dynamic_cast<TH1F*>(f->Get("MyModule/theNormHisto"));
  if (theNormHisto==0) 
  {
    std::cout << "Sorry but cannot find the histo called 'MyModule/theNormHisto'" << std::endl;
    return 0;
  }

  // Accessing to the histo 'theNormHistoByTMEME'
  TH1F* theNormHistoByTMEME = dynamic_cast<TH1F*>(f->Get("MyModule/theNormHistoByTMEME"));
  if (theNormHistoByTMEME==0) 
  {
    std::cout << "Sorry but cannot find the histo called 'MyModule/theNormHistoByTMEME'" << std::endl;
    return 0;
  }

  //----------------------------------------------------------------------------
  //                                  OUTPUT
  //----------------------------------------------------------------------------

  // Opening input file (write-only mode)
  TFile* oFile = TFile::Open("skimming.root","RECREATE");
  if (!oFile->IsOpen()) 
  {
    std::cout << "Sorry but cannot create the file 'skimming.root'" << std::endl;
    return 0;
  }

  // Creating directory MyModule in the TFile
  oFile->cd();
  oFile->mkdir("MyModule");
  Int_t btest = gDirectory->cd("skimming.root:/MyModule");
  if (btest!=1)
  {
    std::cout << "Impossible to create the folder MyModule in the output" << std::endl;
    return 0;
  }
  TDirectory* writeDir = oFile->GetDirectory("MyModule");
  if (writeDir==0)
  {
    std::cout << "Impossible to create the folder MyModule in the output" << std::endl;
    return 0;
  }

  // Creating the output Tree
  TTree* oTree = new TTree("Event", "");
  oTree->Branch("NTEvent", "IPHCTree::NTEvent", &event, 32000, 3);
  oTree->SetDirectory(writeDir);
  oTree->SetAutoSave();

  // Saving input histos in the output TFile 
  theNormHisto -> SetDirectory(writeDir);
  theNormHisto -> Write();
  theNormHistoByTMEME -> SetDirectory(writeDir);
  theNormHistoByTMEME -> Write();  


  //----------------------------------------------------------------------------
  //                               LOOP
  //----------------------------------------------------------------------------

  // Get number of entries in the input Tree
  Long64_t nentries = eventTree->GetEntriesFast();
  std::cout << "N event = " << nentries << std::endl;

  // Initializing counter for filtered events
  Long64_t nselected = 0;

  // Loop over input events
  for (unsigned int ievt=0;ievt<nentries; ievt++)
  {
    // Loading event
    eventTree->LoadTree(ievt);
    eventTree->GetEntry(ievt);
    IPHCTree::NTTransient::InitializeAfterReading(event);

    // Progress bar
    if (ievt%1000==0) 
      std::cout << "ievt=" << ievt << std::endl;

    // Skimming
    if(!NLeptonSkim(event, 2)) continue;
    nselected++;

    // Saving event
    IPHCTree::NTTransient::InitializeBeforeWriting(event);
    oTree->Fill();
  }

  // Writing the non-saved part of the output Tree
  oTree->Write();
  
  // Skimming summary
  std::cout << "N passed event = " << nselected << std::endl;
  std::cout << "Fraction of passed events = " 
            << nselected/static_cast<float>(nentries) *100 
            << " %" << std::endl;

  // End
  std::cout << "#########################" << std::endl;
  std::cout << "    End of the program   " << std::endl;
  std::cout << "#########################" << std::endl;

  return 0;
}

