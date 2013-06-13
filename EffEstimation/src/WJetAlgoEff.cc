#include "../interface/WJetAlgoEff.h"

using namespace std;

WJetAlgoEff::WJetAlgoEff()
{
    // W-tagging
    nbinsW_ = 10;
    xminW_ = 0;
    xmaxW_ = 500;
    WJetAlgoEffNorm = 0;
    cWJetAlgoEff_ = 0;
    cWJetAlgoFake_ = 0;
    legWEff_ = 0;
    legWFake_ = 0;

    // Top-tagging
    nbinsTop_ = 10;
    xminTop_ = 0;
    xmaxTop_ = 500;
    TopJetAlgoEffNorm = 0;
    cTopJetAlgoEff_ = 0;
    cTopJetAlgoFake_ = 0;
    legTopEff_ = 0;
    legTopFake_ = 0;
}

WJetAlgoEff::~WJetAlgoEff()
{
    // W-tagging

    delete WJetAlgoEffNorm;
    delete cWJetAlgoEff_;
    delete cWJetAlgoFake_;
    delete legWEff_;
    delete legWFake_;
    for(unsigned int i=0;i<WJetAlgoEffPlots_.size();i++)
    {
    	delete WJetAlgoEffPlots_[i];
    }
    WJetAlgoEffPlots_.clear();
    for(unsigned int i=0;i<WJetAlgoFakePlots_.size();i++)
    {
    	delete WJetAlgoFakePlots_[i];
    }
    WJetAlgoFakePlots_.clear();
    for(unsigned int i=0;i<WJetAlgoFakeNorm.size();i++)
    {
    	delete WJetAlgoFakeNorm[i];
    }
    WJetAlgoFakeNorm.clear();


    // Top-tagging

    delete TopJetAlgoEffNorm;
    delete cTopJetAlgoEff_;
    delete cTopJetAlgoFake_;
    delete legTopEff_;
    delete legTopFake_;
    for(unsigned int i=0;i<TopJetAlgoEffPlots_.size();i++)
    {
    	delete TopJetAlgoEffPlots_[i];
    }
    TopJetAlgoEffPlots_.clear();
    for(unsigned int i=0;i<TopJetAlgoFakePlots_.size();i++)
    {
    	delete TopJetAlgoFakePlots_[i];
    }
    TopJetAlgoFakePlots_.clear();
    for(unsigned int i=0;i<TopJetAlgoFakeNorm.size();i++)
    {
    	delete TopJetAlgoFakeNorm[i];
    }
    TopJetAlgoFakeNorm.clear();

}

void WJetAlgoEff::InitializeWTag(std::vector<std::string> algoNames, int nbins, float xmin, float xmax){
	WalgoNames_ = algoNames;
	nbinsW_ = nbins;
	xminW_ = xmin;
	xmaxW_ = xmax;
	for (unsigned int jc = 0; jc<WalgoNames_.size(); jc++) 
	{
		string name = WalgoNames_[jc]+"_eff";
		string title = "Efficency of algo: " + WalgoNames_[jc];
		WJetAlgoEffPlots_.push_back( new TH1F(name.c_str(),name.c_str(),nbins,xmin,xmax));
		name = WalgoNames_[jc]+"_fake";
		title = "Fake-rate of algo: " + WalgoNames_[jc];
		WJetAlgoFakePlots_.push_back( new TH1F(name.c_str(),name.c_str(),nbins,xmin,xmax));
		name = WalgoNames_[jc]+"_fake_norm";
		title = "[Norm] Fake-rate of algo: " + WalgoNames_[jc];
		WJetAlgoFakeNorm.push_back( new TH1F(name.c_str(),name.c_str(),nbins,xmin,xmax));
      	}
	// normalisation plots
    	WJetAlgoEffNorm = new TH1F("WJetAlgoEffNorm","WJetAlgoEffNorm",nbins,xmin,xmax);
}
	
void WJetAlgoEff::InitializeTopTag(std::vector<std::string> algoNames, int nbins, float xmin, float xmax){
	TopalgoNames_ = algoNames;
	nbinsTop_ = nbins;
	xminTop_ = xmin;
	xmaxTop_ = xmax;
	for (unsigned int jc = 0; jc<TopalgoNames_.size(); jc++) 
	{
		string name = TopalgoNames_[jc]+"_eff";
		string title = "Efficency of algo: " + TopalgoNames_[jc];
		TopJetAlgoEffPlots_.push_back( new TH1F(name.c_str(),name.c_str(),nbins,xmin,xmax));
		name = TopalgoNames_[jc]+"_fake";
		title = "Fake-rate of algo: " + TopalgoNames_[jc];
		TopJetAlgoFakePlots_.push_back( new TH1F(name.c_str(),name.c_str(),nbins,xmin,xmax));
		name = TopalgoNames_[jc]+"_fake_norm";
		title = "[Norm] Fake-rate of algo: " + TopalgoNames_[jc];
		TopJetAlgoFakeNorm.push_back( new TH1F(name.c_str(),name.c_str(),nbins,xmin,xmax));
      	}
	// normalisation plots
    	TopJetAlgoEffNorm = new TH1F("TopJetAlgoEffNorm","TopJetAlgoEffNorm",nbins,xmin,xmax);
}
	
	
void WJetAlgoEff::LoadEvent(IPHCTree::NTEvent* event, StopMCinfo* mcInfo){

	
	// ----------------
	//   W-tagging
	// ----------------

	// Normalization for efficiencies
	for(unsigned int i=0;i<mcInfo->GetHadronicWs().size();i++){
		WJetAlgoEffNorm->Fill(mcInfo->GetHadronicWs()[i]->p4.Pt());
	}

	// Looking for reconstructed W jets
	for (unsigned int  jc=0; jc<WalgoNames_.size();jc++)
	{
		//change collection link
		if(!event->jets.DoYouKnowCollection(WalgoNames_[jc])){
			cerr<<" The algorithm: "<<WalgoNames_[jc]<<" is not found !!"<<endl;
			continue;
		}
		event->jets.SelectLabel(WalgoNames_[jc]);
		for(unsigned int ij=0;ij<event->jets.size();ij++)
			{
				//if(event->jets[ij].p4.M()>60 && event->jets[ij].p4.M()<100)
				if(event->jets[ij].p4.M()>60 && event->jets[ij].p4.M()<100)
				{
					//cout<<"We've found a jet in a W mass range ! mass = " <<event->jets[ij].p4.M()<<" Pt = "<<event->jets[ij].p4.Pt()<<" ["<<*jc<<"]"<<endl;
					//matching
					bool matched = false;
					vector<IPHCTree::NTGenParticle*> hadronicW = mcInfo->GetHadronicWs();
					for(unsigned int iw=0;iw<hadronicW.size();iw++)
					{
						if(!matched && hadronicW[iw]->p4.DeltaR(event->jets[ij].p4)<0.2)
						{
							//cout<<"matched"<<endl;
							WJetAlgoEffPlots_[jc]->Fill(hadronicW[iw]->p4.Pt());
							//WJetAlgoEff[jc]->Fill(event->jets[ij].p4.Pt());
							matched = true;
							break;
						}
					}
					if(!matched) WJetAlgoFakePlots_[jc]->Fill(event->jets[ij].p4.Pt());
					// WARNING: Fake rate normalization must be computed once per collection
					WJetAlgoFakeNorm[jc]->Fill((event->jets[ij].p4.Pt()));
				}
			}
	}
	
	
	// ----------------
	//   Top-tagging
	// ----------------
 	
	/*

	// Normalization for efficiencies
	for(unsigned int i=0;i<mcInfo->GetHadronicTops().size();i++){
		TopJetAlgoEffNorm->Fill(mcInfo->GetHadronicTops()[i]->p4.Pt());
	}

    
	// Looking for reconstructed Top jets
	for (unsigned int  jc=0; jc<TopalgoNames_.size();jc++)
	{
		//change collection link
		if(!event->jets.DoYouKnowCollection(TopalgoNames_[jc])){
			cerr<<" The algorithm: "<<TopalgoNames_[jc]<<" is not found !!"<<endl;
			continue;
		}
		event->jets.SelectLabel(TopalgoNames_[jc]);
		for(unsigned int ij=0;ij<event->jets.size();ij++)
			{
				//if(event->jets[ij].p4.M()>60 && event->jets[ij].p4.M()<100)
				if(event->jets[ij].p4.M()>120 && event->jets[ij].p4.M()<220)
				{
					//cout<<"Tope've found a jet in a Top mass range ! mass = " <<event->jets[ij].p4.M()<<" Pt = "<<event->jets[ij].p4.Pt()<<" ["<<*jc<<"]"<<endl;
					//matching
					bool matched = false;
					vector<IPHCTree::NTGenParticle*> hadronicTop = mcInfo->GetHadronicTops();
					for(unsigned int iw=0;iw<hadronicTop.size();iw++)
					{
						if(!matched && hadronicTop[iw]->p4.DeltaR(event->jets[ij].p4)<0.2)
						{
							//cout<<"matched"<<endl;
							TopJetAlgoEffPlots_[jc]->Fill(hadronicTop[iw]->p4.Pt());
							//TopJetAlgoEff[jc]->Fill(event->jets[ij].p4.Pt());
							matched = true;
							break;
						}
					}
					if(!matched) TopJetAlgoFakePlots_[jc]->Fill(event->jets[ij].p4.Pt());
					// TopARNING: Fake rate normalization must be computed once per collection
					TopJetAlgoFakeNorm[jc]->Fill((event->jets[ij].p4.Pt()));
				}
			}
	}
	*/


}


void WJetAlgoEff::Compute(){
	  
	  //------------------------
	  // W efficiencies
	  //------------------------

	  //Superposing plots
	  cWJetAlgoEff_ = new TCanvas("cWJetAlgoEff");
	  cWJetAlgoEff_->cd();
	  legWEff_ = new TLegend(0.2,0.3,0.5,0.6);
	  for(unsigned int plotIndex = 0;plotIndex<WJetAlgoEffPlots_.size();plotIndex++){
	  	// normalization of plots
		WJetAlgoEffPlots_[plotIndex]->Divide(WJetAlgoEffNorm);
		WJetAlgoEffPlots_[plotIndex]->Write();
	  	// change colors
		WJetAlgoEffPlots_[plotIndex]->SetLineColor(2+plotIndex);
		WJetAlgoEffPlots_[plotIndex]->SetLineWidth(2);
		legWEff_->AddEntry(WJetAlgoEffPlots_[plotIndex],WalgoNames_[plotIndex].c_str(),"l");
	  	if(plotIndex==0) WJetAlgoEffPlots_[plotIndex]->Draw();
	  	else WJetAlgoEffPlots_[plotIndex]->Draw("same");
	  }
	  legWEff_->Draw("same");

	  //------------------------
	  // W fake rate
	  //------------------------

	  cWJetAlgoFake_ = new TCanvas("cWJetAlgoFake");
	  cWJetAlgoFake_->cd();
	  legWFake_ = new TLegend(0.2,0.3,0.5,0.6);
	  for(unsigned int plotIndex = 0;plotIndex<WJetAlgoFakePlots_.size();plotIndex++){
	  	WJetAlgoFakePlots_[plotIndex]->Divide(WJetAlgoFakeNorm[plotIndex]);
	  	WJetAlgoFakePlots_[plotIndex]->Write();
	  	WJetAlgoFakePlots_[plotIndex]->SetLineColor(2+plotIndex);
		WJetAlgoFakePlots_[plotIndex]->SetLineWidth(2);
		legWFake_->AddEntry(WJetAlgoFakePlots_[plotIndex],WalgoNames_[plotIndex].c_str(),"l");
	  	if(plotIndex==0) WJetAlgoFakePlots_[plotIndex]->Draw();
	  	else WJetAlgoFakePlots_[plotIndex]->Draw("same");
	  }
	  legWFake_->Draw("same");
	  
	  //------------------------
	  // Top efficiencies
	  //------------------------

	  //Superposing plots
	  cTopJetAlgoEff_ = new TCanvas("cTopJetAlgoEff");
	  cTopJetAlgoEff_->cd();
	  legTopEff_ = new TLegend(0.2,0.3,0.5,0.6);
	  for(unsigned int plotIndex = 0;plotIndex<TopJetAlgoEffPlots_.size();plotIndex++){
	  	// normalization of plots
		TopJetAlgoEffPlots_[plotIndex]->Divide(TopJetAlgoEffNorm);
		TopJetAlgoEffPlots_[plotIndex]->Write();
	  	// change colors
		TopJetAlgoEffPlots_[plotIndex]->SetLineColor(2+plotIndex);
		TopJetAlgoEffPlots_[plotIndex]->SetLineWidth(2);
		legTopEff_->AddEntry(TopJetAlgoEffPlots_[plotIndex],TopalgoNames_[plotIndex].c_str(),"l");
	  	if(plotIndex==0) TopJetAlgoEffPlots_[plotIndex]->Draw();
	  	else TopJetAlgoEffPlots_[plotIndex]->Draw("same");
	  }
	  legTopEff_->Draw("same");

	  //------------------------
	  // Top fake rate
	  //------------------------

	  cTopJetAlgoFake_ = new TCanvas("cTopJetAlgoFake");
	  cTopJetAlgoFake_->cd();
	  legTopFake_ = new TLegend(0.2,0.3,0.5,0.6);
	  for(unsigned int plotIndex = 0;plotIndex<TopJetAlgoFakePlots_.size();plotIndex++){
	  	TopJetAlgoFakePlots_[plotIndex]->Divide(TopJetAlgoFakeNorm[plotIndex]);
	  	TopJetAlgoFakePlots_[plotIndex]->Write();
	  	TopJetAlgoFakePlots_[plotIndex]->SetLineColor(2+plotIndex);
		TopJetAlgoFakePlots_[plotIndex]->SetLineWidth(2);
		legTopFake_->AddEntry(TopJetAlgoFakePlots_[plotIndex],TopalgoNames_[plotIndex].c_str(),"l");
	  	if(plotIndex==0) TopJetAlgoFakePlots_[plotIndex]->Draw();
	  	else TopJetAlgoFakePlots_[plotIndex]->Draw("same");
	  }
	  legTopFake_->Draw("same");


}

void WJetAlgoEff::Write(TFile* foutAlgo)
{
	if(!foutAlgo) return;
	foutAlgo->cd();
	foutAlgo->mkdir("WJetAlgoPerf");
	foutAlgo->cd("WJetAlgoPerf");
	cWJetAlgoEff_->Write();
	cWJetAlgoFake_->Write();
	cTopJetAlgoEff_->Write();
	cTopJetAlgoFake_->Write();
}
