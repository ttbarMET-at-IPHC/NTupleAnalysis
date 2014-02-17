#ifndef combined1LeptonStopSelection_h
#define combined1LeptonStopSelection_h

// IPHC headers
#include "Selection/interface/Selection.h"
#include "Selection/interface/SelectionTable.h"
#include "Tools/interface/Dataset.h"

#include "EventReco/interface/StopAnaReco.h"
#include "EventReco/interface/Resolution.h"
#include "EventReco/interface/Mt2Com_bisect.h"
#include "BTagReshaping/interface/BTagReshaping.h"

#include "JetCorrections/interface/JetCorrectionUncertainty.h"

// system include files
#include <memory>
#include <vector>
#include <sstream>

class combined1LeptonStopSelection: public Selection
{

    // -------------------------------------------------------------
    //                       method members
    // -------------------------------------------------------------
    public:

        //! Constructor without argument
        combined1LeptonStopSelection();

        //! Copy constructor
        combined1LeptonStopSelection(const combined1LeptonStopSelection &);

        //! Destructor
        ~combined1LeptonStopSelection()
        { }

        //! LoadEvent 
        bool LoadEvent(const NTEvent * event) 
        { 
            return Event::LoadEvent(event); 
        }

        void setBTagReshapingInput(string fileName);     
        void setMCJetCorrectionsInput(string fileName);  
        void setDataJetCorrectionsInput(string fileName);
        void setPileUpReweightingInput(string fileName); 
        void loadCorrections();

        // #####################################
        // #   Accessors to selected objects   #
        // #####################################

        //! Accessor to good electrons
        const std::vector<IPHCTree::NTElectron>& GetGoodElectrons() const { return goodElectrons; }

        //! Accessor to good muons
        const std::vector<IPHCTree::NTMuon>& GetGoodMuons() const { return goodMuons; }
        
        //! Accessor to good jets
        const std::vector<IPHCTree::NTJet>& GetGoodJets() const { return goodJets; }

        //! Accessor to selected jets
        const std::vector<IPHCTree::NTJet>& GetJetsForAna() const { return selectedJets; }

        //! Accessor to selected b-jets
        const std::vector<IPHCTree::NTJet>& GetBJetsForAna() const { return selectedBJets; }


        // Accessors to selected leptons

        TLorentzVector getTheLeadingLepton()      { return theLeadingLepton; };
        int            getTheLeadingLeptonPDGId() { return theLeadingLepton_pdgid; };
        TLorentzVector getTheSecondLepton()       { return theSecondLepton; };
        int            getTheSecondLeptonPDGId()  { return theSecondLepton_pdgid; };


        int getTheNumberOfSelectedLeptons() { return numberOfSelectedLeptons; }

        // #################################
        // #   Selection-related methods   #
        // #################################

        // Apply object selection on event
        void doObjectSelection(bool runningOnData, short int JESvariation = 0);

        // Return true if the event pass the selection
        bool passEventSelection(bool runningOnData);

        // Apply JES +/- 1sigma variation on a given jet collection
        void ApplyJESVariation(bool runningOnData, bool upOrDown);

        void FillKinematicP4();

        // Isolated track veto
        bool GetSUSYstopIsolatedTrackVeto(TLorentzVector lepton_p, float lepton_charge) const;

        // Tau veto
        bool GetSUSYstopTauVeto(TLorentzVector lepton_p, float lepton_charge) const;

        // Good muons
        std::vector<IPHCTree::NTMuon> GetSUSYstopGoodMuons()     const;

        // Selected muons
        std::vector<IPHCTree::NTMuon> GetSUSYstopSelectedMuons() const;

        // Good electrons
        std::vector<IPHCTree::NTElectron> GetSUSYstopGoodElectrons    () const;

        // Selected electrons
        std::vector<IPHCTree::NTElectron> GetSUSYstopSelectedElectrons() const;

        // Good jets
        std::vector<IPHCTree::NTJet>  GetSUSYstopGoodJets(int DataType) const;
        
        // Selected jets
        std::vector<IPHCTree::NTJet> GetSUSYstopSelectedJets(int DataType) const;

        // Selected b-jets
        std::vector<IPHCTree::NTJet> GetSUSYstopSelectedBJets(int DataType) const; 

        // Type-1 corrected PF MET
        IPHCTree::NTMET GetSUSYstopType1MET(int DataType) const;

        // Type-1 + Phi corrected PF MET
        IPHCTree::NTMET GetSUSYstopType1PhiMET(int DataType) const;

        // Check if a given trigger path has been trigerred
        bool checkPathHasBeenFired(string path);

        // Check if a given set of trigger (single electron, double muon, ...) contains a path that has been fired
        bool passTrigger(string channel);

        // #######################################
        // #   High-level variables definition   # 
        // #######################################

        float M3() const             { return top_hadronic.M();}
        float M_topleptonic() const  { return top_leptonic.M();}
        float MT_wleptonic() const   { return sqrt( 2.* theLeadingLepton.Pt() * theMET.Pt() *(1. - cos( theLeadingLepton.Phi() - theMET.Phi()) )) ;}
        float PT_tophad() const      { return top_hadronic.Pt();}
        float PT_topleptonic() const { return top_leptonic.Pt();}
        float PT_wleptonic()  const  { return w_leptonic.Pt();}
        float Dphi_lmet() const      { return theLeadingLepton.DeltaPhi(theMET);}
        float rawMet() const         { return rawMET.p2.Mod(); }
        float Met() const            { return theMET.Pt(); }

        // HT
        float HT() const      
        { 
            float HT_ = 0;
            for (unsigned int i = 0 ; i < selectedJets.size() ; i++)
                HT_ += selectedJets[i].p4.Pt();

            return HT_;
        }

        // HT ratio
        float HT_ratio() const      
        { 
            float HT_onTheSideOfMET = 0;
            float HT_total = 0;
            for (unsigned int i = 0 ; i < selectedJets.size() ; i++)
            {
                if (abs(theMET.DeltaPhi(selectedJets[i].p4)) < 3.1415/2.0) 
                    HT_onTheSideOfMET += selectedJets[i].p4.Pt();

                HT_total += selectedJets[i].p4.Pt();
            }

            return HT_onTheSideOfMET / HT_total;
        }

        // Hadronic chi2
        float HadronicChi2(bool runningOnData) const
        {
            Resolution Chi2;
            float value = Chi2.GetChi2(GetJetsForAna(),runningOnData);
            return value;
        }

        // MT2W
        float MT2W() const
        {
            Mt2Com_bisect Mt2;
            float value = Mt2.calculateMT2w(GetJetsForAna(),
                    GetBJetsForAna(),
                    theLeadingLepton,
                    theMET.Vect().XYvector(),
                    "MT2w");
            return value;
        }

        // DeltaPhi(MET, two first leading jets)
        float DPhi_MET_leadingJets() const
        {
            TLorentzVector firstLeadingJet;
            TLorentzVector secondLeadingJet;

            for (unsigned int j = 0 ; j < selectedJets.size() ; j++)
            {
                if (selectedJets[j].p4.Pt() > firstLeadingJet.Pt())
                    firstLeadingJet = selectedJets[j].p4;

                else if (selectedJets[j].p4.Pt() > secondLeadingJet.Pt())
                    secondLeadingJet = selectedJets[j].p4;
            }

            if (abs(theMET.DeltaPhi(firstLeadingJet)) < abs(theMET.DeltaPhi(secondLeadingJet)))
                return abs(theMET.DeltaPhi(firstLeadingJet));
            else
                return abs(theMET.DeltaPhi(secondLeadingJet));
        }

        // M3b
        float M3b()
        {
            TLorentzVector sum;
            if(selectedJets.size() == 3)
            {
                sum = selectedJets[0].p4 + selectedJets[1].p4 + selectedJets[2].p4;
            }
            else
            {
                // check which jet is closest to lepton, then take other 3
                double dphimin = 99.;
                int index_closest_jet = -1;
                for(int i=0; i < 4; i++)
                {
                    double dphi = theLeadingLepton.DeltaPhi(selectedJets[i].p4);
                    if (dphi < dphimin)
                    {
                        dphimin = dphi;
                        index_closest_jet = i;
                    }
                }

                for(int i=0; i<4; i++)
                {
                    if(i!=index_closest_jet)
                        sum = sum + selectedJets[i].p4;
                }

            }    

            return sum.M();    

        }
        
        TLorentzVector leadingBJet()
        {
            TLorentzVector leadingBJet = 0.0;
            
            for (unsigned int i = 0 ; i < selectedBJets.size() ; i++)
            {
                if (leadingBJet.Pt() < selectedBJets[i].p4.Pt())
                    leadingBJet      = selectedBJets[i].p4;
            }

            if (leadingBJet.Pt() == 0.0) return leadingJet();

            return leadingBJet;
        }

        TLorentzVector leadingJet()
        {
            TLorentzVector leadingJet = 0.0;
            
            for (unsigned int i = 0 ; i < selectedJets.size() ; i++)
            {
                if (leadingJet.Pt() < selectedJets[i].p4.Pt())
                    leadingJet      = selectedJets[i].p4;
            }

            return leadingJet;
        }

        TLorentzVector leadingJetByCSV(bool runningOnData)
        {
            TLorentzVector leadingCSVJet = 0.0;
            float maxDiscr = -1.0;

            string discrName;
            if (runningOnData) discrName = "combinedSecondaryVertexBJetTags";
            else               discrName = "zz1combinedSecondaryVertexBJetTagsReshapeNominal";

            for (unsigned int i = 0 ; i < selectedJets.size() ; i++)
            {
                if (maxDiscr < selectedJets[i].bTag[discrName])
                {
                    maxDiscr = selectedJets[i].bTag[discrName];
                    leadingCSVJet = selectedJets[i].p4;
                }
            }
            
            return leadingCSVJet;
        }

        TLorentzVector bJetClosestToLeadingLepton()
        {
            float minDeltaRLeptonB = 9999.0;
            TLorentzVector closestBJet;

            for (unsigned int i = 0 ; i < selectedBJets.size() ; i++)
            if (minDeltaRLeptonB    > selectedBJets[i].p4.DeltaR(theLeadingLepton))
            {
                closestBJet         = selectedBJets[i].p4;
                minDeltaRLeptonB    = closestBJet.DeltaR(theLeadingLepton);
            }

            if (minDeltaRLeptonB == 9999.0) return leadingJet();
            /*
            for (unsigned int i = 0 ; i < selectedJets.size() ; i++)
            if (minDeltaRLeptonB    > selectedJets[i].p4.DeltaR(theLeadingLepton))
            {
                closestBJet         = selectedJets[i].p4;
                minDeltaRLeptonB    = closestBJet.DeltaR(theLeadingLepton);
            }
            */

            return closestBJet;
        }

        float getPileUpWeight( int nvtx, bool runningOnData )
        {
            // From Verena

            if( runningOnData ) return 1;

            if( nvtx > pileUpReweightingHisto->GetNbinsX() )
                nvtx = pileUpReweightingHisto->GetNbinsX();

            float weight = 0;
            weight = pileUpReweightingHisto->GetBinContent( pileUpReweightingHisto->FindBin(nvtx) );
        
            //we don't want to kill events bc they have no weight
            if( weight <= 0 ) weight = 1.;
            return weight;

        }


        // -------------------------------------------------------------
        //                        data members
        // -------------------------------------------------------------
    private:

        // Raw objects

        std::vector<IPHCTree::NTJet>      rawJets;
        NTMET                             rawMET;

        // Partially selected objects

        std::vector<IPHCTree::NTElectron> goodElectrons;
        std::vector<IPHCTree::NTMuon>     goodMuons;
        std::vector<IPHCTree::NTJet>      goodJets;

        // Objects for analysis

        int numberOfSelectedLeptons;
        TLorentzVector theLeadingLepton;
        TLorentzVector theSecondLepton;
        int theLeadingLepton_pdgid;
        int theSecondLepton_pdgid;

        std::vector<IPHCTree::NTElectron> selectedElectrons;
        std::vector<IPHCTree::NTMuon>     selectedMuons;
        std::vector<IPHCTree::NTJet>      selectedJets;
        std::vector<IPHCTree::NTJet>      selectedBJets;
        TLorentzVector theMET;

        // High-level objects

        TLorentzVector all_hadronic;
        TLorentzVector top_hadronic;
        TLorentzVector the_4thjet;  
        TLorentzVector w_leptonic;
        TLorentzVector top_leptonic;

        // Corrections managers

        string bTagReshapingInput;     
        string jetCorrectionsMCInput;
        string jetCorrectionsDataInput;
        string pileUpReweightingInput;

        TH1F* pileUpReweightingHisto;
        JetCorrectionUncertainty* JESUncertainty_MC;
        JetCorrectionUncertainty* JESUncertainty_Data;

        BTagShapeInterface* BTagReshape_nominal; 
        BTagShapeInterface* BTagReshape_upBC; 
        BTagShapeInterface* BTagReshape_downBC;
        BTagShapeInterface* BTagReshape_upLight;
        BTagShapeInterface* BTagReshape_downLight;

};

#endif
