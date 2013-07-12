

#include <iomanip>
#include <iostream>

#include <TFile.h>
#include <TH2F.h>
#include <TMath.h>

using namespace std;





int main()
{



    float e_tag = 0.65;
    float e_mis = 0.4;
    
    float p10[5];
    float p11[5];
    float p20[5];
    float p21[5];

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
          cout << "n   = " << n 
          << " ; p10 = " << p10[n-1]
          << " ; p11 = " << p11[n-1]
          << " ; p20 = " << p20[n-1]
          << " ; p21 = " << p21[n-1]
          << endl;
    }

    vector<string> processList;                         vector<float> yield3jets;       vector<float> yield4jets;
    processList.push_back("1lepTop");                   yield3jets.push_back(9.34);     yield4jets.push_back(4.88);
    processList.push_back("ttbarFullLept");             yield3jets.push_back(20.31);    yield4jets.push_back(12.32);
    processList.push_back("W+jets");                    yield3jets.push_back(16.63);    yield4jets.push_back(4.25);
    processList.push_back("rare");                      yield3jets.push_back(9.44);     yield4jets.push_back(5.71);
    processList.push_back("signal");                    yield3jets.push_back(94.44);    yield4jets.push_back(70.03);

    TFile f("multiplicityAK8_highPt_vs_matched.root","UPDATE");
    TH2F* highPt_vs_matched = 0;
    for (unsigned int p = 0 ; p < processList.size() ; p++)
    {
        float selectionProbability = 0;
        // Get the histo from the file
        f.GetObject(processList[p].c_str(),highPt_vs_matched);
        
        // Normalize histo
        highPt_vs_matched->Scale(1.0/highPt_vs_matched->Integral());

        //cout << endl << endl;
        //cout << "process : " << processList[p] << endl;

        // Compute total probability
        for (int n = 1 ; n <= 5 ; n++)
        {
        //    cout << "proba n = " << n << " ; matched = 0      " << highPt_vs_matched->GetBinContent(highPt_vs_matched->FindBin(n,0)) << endl;
        //    cout << "proba n = " << n << " ; matched = 1      " << highPt_vs_matched->GetBinContent(highPt_vs_matched->FindBin(n,1)) << endl;

            selectionProbability += highPt_vs_matched->GetBinContent(highPt_vs_matched->FindBin(n,0)) * (p10[n-1] + p20[n-1])
                                 +  highPt_vs_matched->GetBinContent(highPt_vs_matched->FindBin(n,1)) * (p11[n-1] + p21[n-1]);
        }

        cout << "process : " << processList[p] << " ; p = " << selectionProbability << " ; previous yield : " << yield3jets[p] << " ; new yield : " << yield3jets[p] * selectionProbability << endl;
        yield3jets[p] *= selectionProbability;

    }

    float b_new = 0;
    float s_new = 0;
    float b_ref = 0;
    float s_ref = 0;

    for (unsigned int p = 0 ; p < processList.size() ; p++)
    {
        if (processList[p] != "signal")
        {
            b_new += yield3jets[p];
            b_ref += yield4jets[p];
        }
        else
        {
            s_new += yield3jets[p];
            s_ref += yield4jets[p];
        }
    }

    cout << endl;
    cout << " g   =   s_new / s_ref   =  " << (s_new / sqrt(b_new))   /   (s_ref / sqrt(b_ref))    << endl;;
    cout << endl;


    return 0;
}
