

#include <iomanip>
#include <iostream>

#include <TFile.h>
#include <TH2F.h>
#include <TMath.h>

using namespace std;



void simulateWtagRequirement(vector<string> processList, vector<float> inputYields, float e_tag, float e_mis, float* newB, float* newS);

int main()
{

    TH1F* sig_functionOf_mistag      = new TH1F("sig_mistag","S = f(mistag) for tag = 0.8",101,-0.005,1.005);
    TH1F* sig_functionOf_tag         = new TH1F("sig_tag","S = f(tag) for mistag = 0.4",101,-0.005,1.005);
    TH2F* sig_functionOf_tagVsMistag = new TH2F("sig_tagVsMistag","S = f(tag,mistag)",101,-0.005,1.005,101,-0.005,1.005);
    
    TH1F* gS_functionOf_mistag      = new TH1F("gS_mistag","gS = f(mistag) for tag = 0.8",101,-0.005,1.005);
    TH1F* gS_functionOf_tag         = new TH1F("gS_tag","gS = f(tag) for mistag = 0.4",101,-0.005,1.005);
    TH2F* gS_functionOf_tagVsMistag = new TH2F("gS_tagVsMistag","gS = f(tag,mistag)",101,-0.005,1.005,101,-0.005,1.005);

    TH1F* pur_functionOf_mistag      = new TH1F("pur_mistag","P = f(mistag) for tag = 0.8",101,-0.005,1.005);
    TH1F* pur_functionOf_tag         = new TH1F("pur_tag","P = f(tag) for mistag = 0.4",101,-0.005,1.005);
    TH2F* pur_functionOf_tagVsMistag = new TH2F("pur_tagVsMistag","P = f(tag,mistag)",101,-0.005,1.005,101,-0.005,1.005);
    
    TH1F* gP_functionOf_mistag      = new TH1F("gP_mistag","gP = f(mistag) for tag = 0.8",101,-0.005,1.005);
    TH1F* gP_functionOf_tag         = new TH1F("gP_tag","gP = f(tag) for mistag = 0.4",101,-0.005,1.005);
    TH2F* gP_functionOf_tagVsMistag = new TH2F("gP_tagVsMistag","gP = f(tag,mistag)",101,-0.005,1.005,101,-0.005,1.005);

    vector<string> processList;                vector<float> yield3jets;       vector<float> yield4jets;
    processList.push_back("1lepTop");          yield3jets.push_back(9.34);     yield4jets.push_back(4.88);
    processList.push_back("ttbarFullLept");    yield3jets.push_back(20.31);    yield4jets.push_back(12.32);
    processList.push_back("W+jets");           yield3jets.push_back(9.44);     yield4jets.push_back(5.71);
    processList.push_back("rare");             yield3jets.push_back(16.63);    yield4jets.push_back(4.25);
    processList.push_back("signal");           yield3jets.push_back(94.44/9.0);    yield4jets.push_back(70.03/9.0);

    float b_ref = 0;
    float s_ref = 0;

    for (unsigned int p = 0 ; p < processList.size() ; p++)
    {
        if (processList[p] != "signal") b_ref += yield4jets[p];
        else                            s_ref += yield4jets[p];
    }
        
    float sig_ref = (s_ref/sqrt(b_ref));    float pur_ref = s_ref / (s_ref + b_ref);

    // ---------------

    for (float e_tag = 0.0 ; e_tag <= 1.0 ; e_tag += 0.01)
    for (float e_mis = 0.0 ; e_mis <= 1.0 ; e_mis += 0.01)
    {
        float b_new = 0;    float s_new = 0;

        simulateWtagRequirement(processList,yield3jets,e_tag,e_mis,&b_new,&s_new);

        float sig_new = (s_new/sqrt(b_new));    float pur_new = s_new / (s_new + b_new);
        float gsig = sig_new / sig_ref;         float gpur = pur_new / pur_ref;
       
        sig_functionOf_tagVsMistag->Fill(e_tag,e_mis,sig_new);
        gS_functionOf_tagVsMistag->Fill(e_tag,e_mis,gsig);
        pur_functionOf_tagVsMistag->Fill(e_tag,e_mis,pur_new);
        gP_functionOf_tagVsMistag->Fill(e_tag,e_mis,gpur);
    }

    {
        float e_tag = 0.8;
        for (float e_mis = 0.0 ; e_mis <= 1.0 ; e_mis += 0.01)
        {
            float b_new = 0;    float s_new = 0;

            simulateWtagRequirement(processList,yield3jets,e_tag,e_mis,&b_new,&s_new);

            float sig_new = (s_new/sqrt(b_new));    float pur_new = s_new / (s_new + b_new);
            float gsig = sig_new / sig_ref;         float gpur = pur_new / pur_ref;
           
            sig_functionOf_mistag->Fill(e_mis,sig_new);
            gS_functionOf_mistag->Fill(e_mis,gsig);
            pur_functionOf_mistag->Fill(e_mis,pur_new);
            gP_functionOf_mistag->Fill(e_mis,gpur);
        }
    }

    {
        float e_mis = 0.4;
        for (float e_tag = 0.0 ; e_tag <= 1.0 ; e_tag += 0.01)
        {
            float b_new = 0;    float s_new = 0;

            simulateWtagRequirement(processList,yield3jets,e_tag,e_mis,&b_new,&s_new);

            float sig_new = (s_new/sqrt(b_new));    float pur_new = s_new / (s_new + b_new);
            float gsig = sig_new / sig_ref;         float gpur = pur_new / pur_ref;
           
            sig_functionOf_tag->Fill(e_tag,sig_new);
            gS_functionOf_tag->Fill(e_tag,gsig);
            pur_functionOf_tag->Fill(e_tag,pur_new);
            gP_functionOf_tag->Fill(e_tag,gpur);
        }
    }

    TFile f("gainEstimation.root","RECREATE");

    sig_functionOf_mistag->Write();
    gS_functionOf_mistag->Write();
    pur_functionOf_mistag->Write();
    gP_functionOf_mistag->Write();
 
    sig_functionOf_tag->Write();
    gS_functionOf_tag->Write();
    pur_functionOf_tag->Write();
    gP_functionOf_tag->Write();
   
    sig_functionOf_tagVsMistag->Write();
    gS_functionOf_tagVsMistag->Write();
    pur_functionOf_tagVsMistag->Write();
    gP_functionOf_tagVsMistag->Write();

    f.Close();

}


void simulateWtagRequirement(vector<string> processList, vector<float> inputYields, float e_tag, float e_mis, float* newB, float* newS)
{

    float p10[5];
    float p11[5];
    float p20[5];
    float p21[5];

    cout << endl;
    cout << "(e_tag,e_mis) = (" << e_tag << "," << e_mis << ")" << endl;

    for (int n = 1 ; n <= 5 ; n++)
    {
        //          Binomial term          Matched jet is   Unmatched jets      Unmatched jets
        //                                 or isnt tagged  that are tagged    that aren't tagged
        
        p10[n-1] = TMath::Binomial( n ,1) *                    e_mis        * pow(1 - e_mis , n-1);

        p11[n-1] =                             e_tag                        * pow(1 - e_mis , n-1)
                 + TMath::Binomial(n-1,1) * (1 - e_tag)  *     e_mis        * pow(1 - e_mis , n-2);

        p20[n-1] = TMath::Binomial( n ,2) *                 e_mis*e_mis     * pow(1 - e_mis , n-2);

        p21[n-1] = TMath::Binomial(n-1,1) *    e_tag     *     e_mis        * pow(1 - e_mis , n-2)
                 + TMath::Binomial(n-1,2) * (1 - e_tag)  *  e_mis*e_mis     * pow(1 - e_mis , n-3);

    }

    TFile f("multiplicityAK8_highPt_vs_matched.root","UPDATE");
    TH2F* highPt_vs_matched = 0;
    for (unsigned int p = 0 ; p < processList.size() ; p++)
    {
        float selectionProbability = 0;
        // Get the histo from the file
        f.GetObject(processList[p].c_str(),highPt_vs_matched);
        
        // Normalize histo
        highPt_vs_matched->Scale(1.0/highPt_vs_matched->Integral());

        // Compute total probability
        for (int n = 1 ; n <= 5 ; n++)
        {
            selectionProbability += highPt_vs_matched->GetBinContent(highPt_vs_matched->FindBin(n,0)) * (p10[n-1] + p20[n-1])
                                 +  highPt_vs_matched->GetBinContent(highPt_vs_matched->FindBin(n,1)) * (p11[n-1] + p21[n-1]);
        }


        cout << " [" << processList[p] << "] | p = " << selectionProbability; 
        cout                           << "  | Base yield : " << inputYields[p];
        cout                           << "  | After W-tag req. : " << inputYields[p] * selectionProbability;
        cout << endl;
        inputYields[p] *= selectionProbability;
    }

    float s = 0;
    float b = 0;

    for (unsigned int p = 0 ; p < processList.size() ; p++)
    {
        if (processList[p] != "signal") b += inputYields[p];
        else                            s += inputYields[p];
    }

    *newB = b;
    *newS = s;
}
