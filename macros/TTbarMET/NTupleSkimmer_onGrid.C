#include <iomanip>
#include <iostream>
#include <time.h>
#include "../../../../IPHCDataFormat/NTFormat/interface/NTEvent.h"
#include "../../../../IPHCDataFormat/NTFormat/interface/NTTransient.h"

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TDirectory.h>

bool NLeptonSkim(IPHCTree::NTEvent * event, int Nlep, int Nj)
{

    short nTriggers = 0;
    short nElectrons = 0;
    short nMuons = 0;
    short nJets = 0;
   
    // -------------------------------------
    // Trigger
    // -------------------------------------

    std::vector<IPHCTree::NTTriggerPathType> myPaths;
	event->trigger.GetSubTable("HLT_IsoMu24_v*",myPaths);
	for (unsigned int i=0;i<myPaths.size();i++) 
	    if (myPaths[i].fired==1) nTriggers++;
	
    event->trigger.GetSubTable("HLT_IsoMu24_eta2p1_v*",myPaths);
	for (unsigned int i=0;i<myPaths.size();i++) 
	    if (myPaths[i].fired==1) nTriggers++;
	
    event->trigger.GetSubTable("HLT_Ele27_WP80_v*",myPaths);
	for (unsigned int i=0;i<myPaths.size();i++) 
	    if (myPaths[i].fired==1) nTriggers++;
    
    if(nTriggers == 0) return false;

    // -------------------------------------
    // Electrons
    // -------------------------------------
    {
      event->electrons.SelectLabel("selectedPatElectrons");
      for(unsigned int i=0;i<event->electrons.size();i++)
      {
         if (event->electrons[i].p4.Pt()>25 && abs(event->electrons[i].p4.Eta()) < 2.1 ) nElectrons++;
      }

    }
    
    // -------------------------------------
    // Muons
    // -------------------------------------
    {
      event->muons.SelectLabel("selectedPatMuons");
      for(unsigned int i=0;i<event->muons.size();i++)
      {
         if (event->muons[i].p4.Pt()>25 && abs(event->muons[i].p4.Eta()) < 2.5 ) nMuons++;
      }

    }

    if (nElectrons+nMuons < Nlep) return false;

    //-------------------------------------
    // JETS
    //-------------------------------------
    {
      event->jets.SelectLabel("pf");
      for(unsigned int i=0;i<event->jets.size();i++)
      {
         if (event->jets[i].p4.Pt()>25 && abs(event->jets[i].p4.Eta()) < 3 ) nJets++;
      }

    }
    
    if (nJets < Nj) return false;

    if(nTriggers >= 1 && nElectrons+nMuons >= Nlep && nJets >= Nj) return true;
    else return false;

}


int main (int argc, char *argv[])
{
 
//  std::cout << "#########################" << std::endl;
//  std::cout << "Beginning of the program"  << std::endl;
//  std::cout << "#########################" << std::endl;
  
  if(argc<5) { std::cout << "Syntax is : NTupleSkimmer  outputFolder outputFile inputFolder inputFile1 inputFile2 ..."  << std::endl; return 1;}
  
  // Input file
  
  TString outputFolder = argv[1];
  TString outputFile = argv[2];
  TString inputFolder = argv[3];
  vector<TString> inputFile;
  for (int i = 4 ; i < argc ; i++)
      inputFile.push_back(string(argv[i]));

  // Output file
  TString outFileName = outputFolder+"/"+outputFile;

  // Declaring a pointer to the current event
  IPHCTree::NTEvent * event = 0;

  //----------------------------------------------------------------------------
  //                                  OUTPUT
  //----------------------------------------------------------------------------
  
  // Opening input file (write-only mode)
  TFile* oFile = TFile::Open(outFileName,"RECREATE");
  if (!oFile->IsOpen())
  {
    std::cout << "Sorry but cannot create the file " << outFileName << std::endl;
    return 0;
  }

  // Creating directory MyModule in the TFile
  oFile->cd();
  oFile->mkdir("MyModule");
  TString nametest = outFileName+":/MyModule";
  Int_t btest = gDirectory->cd(nametest);
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
 
  TH1F* theNormHisto_total;
  TH1F* theNormHistoByTMEME_total;

  for (int iInput = 0 ; iInput < inputFile.size() ; iInput++)
  {

      //----------------------------------------------------------------------------
      //                           Open Input file
      //----------------------------------------------------------------------------

      TString file = inputFolder+"/"+inputFile[iInput] ;

      TFile* f = TFile::Open( file );
      
      if (!f->IsOpen())
      {
        std::cout << "Sorry but cannot find the file '" << file << "'" << std::endl;
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

      //----------------------------------------------------------------------------
      //                         NormHisto management
      //----------------------------------------------------------------------------


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

      if (iInput == 0)
      {
          theNormHisto_total = (TH1F*) theNormHisto->Clone();
          theNormHistoByTMEME_total = (TH1F*) theNormHistoByTMEME->Clone();
      }
      else
      {
          theNormHisto_total->Add(theNormHisto);
          theNormHistoByTMEME_total->Add(theNormHistoByTMEME);
      }

      //----------------------------------------------------------------------------
      //                               LOOP
      //----------------------------------------------------------------------------

      // Get number of entries
      Long64_t nentries = eventTree->GetEntriesFast();
      //nentries = 10000;
      //std::cout << "N event = " << nentries << std::endl;

      // Initializing counter for filtered events
      Long64_t nselected = 0;

      for (unsigned int ievt=0;ievt< nentries; ievt++)
      {
        // loading event
        eventTree->LoadTree(ievt);
        eventTree->GetEntry(ievt);
        IPHCTree::NTTransient::InitializeAfterReading(event);

      //  if (ievt%100000==0) std::cout << "ievt=" << ievt << std::endl;

        if(!NLeptonSkim(event, 1, 3 )) continue; 
        nselected++;

        // Saving event
        IPHCTree::NTTransient::InitializeBeforeWriting(event);
        oTree->Fill();
      }

      // Get number of entries
      std::cout << "Skimmed " << inputFile[iInput] << " : "
                << nselected << "/" << nentries
                << " (" << nselected/static_cast<float>(nentries) *100 << "%) passed"
                << endl;
      
  }

  oFile->cd();

  // Saving input histos in the output TFile 
  theNormHisto_total -> SetDirectory(writeDir);
  theNormHisto_total -> Write();
  theNormHistoByTMEME_total -> SetDirectory(writeDir);
  theNormHistoByTMEME_total -> Write();

  // Writing the non-saved part of the output Tree
  oTree->Write();



//  std::cout << "#########################" << std::endl;
//  std::cout << "    End of the program   " << std::endl;
//  std::cout << "#########################" << std::endl;

  return 0;
}

