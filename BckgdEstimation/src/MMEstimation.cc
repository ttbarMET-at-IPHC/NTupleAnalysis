#include "BckgdEstimation/interface/MMEstimation.h"
		
MMEstimation::MMEstimation(vector<Dataset> datasets, float isoLoose, float isoTight, unsigned int n_bins, float min, float max, string Channel){

      theNBins = n_bins;
      theMinBin = min;
      theMaxBin = max;

      theIsolations.iso1[0] = isoLoose;
      theIsolations.iso1[1] = isoTight;
      theIsolations.iso1[2] = isoTight;
      theIsolations.iso2[0] = isoLoose;
      theIsolations.iso2[1] = isoLoose;
      theIsolations.iso2[2] = isoTight;

      theChannel = Channel;
      //Reinitialized to a value, during RunTheMatrixMethod() call, for each bin and each isolation
      for(unsigned int iso_index=0; iso_index<3; iso_index++){
       locNSelected[iso_index] = 0;
       locNSelectedSignal[iso_index] = 0;
       locNSelectedW[iso_index] = 0;
       locNSelectedQCD[iso_index] = 0;
      }

      //Reinitialized to a value, during RunTheMatrixMethod() call, for each bin
      EpsilonFake = 0;
      EpsilonFakeErr = 0;
      EpsilonSignal = 0;
      EpsilonSignalErr = 0;

      //Reinitialized to a value, for each bin, at each IncrementNSetBin() call
      //Then passed bin by bin to RunTheMatrixMethod() call
      theNSelected.reserve(100000);
      for(unsigned int bin_index=0; bin_index<theNBins; bin_index++){
       struct NSelected theNSelElement;
       for(unsigned int iso_index=0; iso_index<3; iso_index++){
 	 theNSelElement.NSel[iso_index] = 0.;
       }
       theNSelected.push_back(theNSelElement);
      }

      IsoNames[0] = "Loose";
      IsoNames[1] = "Medium";
      IsoNames[2] = "Tight";

      //Reinitialized to a value, during RunTheMatrixMethod() call, for each bin, each isolation and for each iteration a new value is filled
      theDistributions.reserve(100000);
      for(unsigned int bin_index=0; bin_index<theNBins; bin_index++){
       std::stringstream ss;
       ss << bin_index;
       struct Distribution theDistributionElement;
       for(unsigned int iso_index=0; iso_index<3; iso_index++){
        theDistributionElement.NMMEstimatedQCDDistribution[iso_index] = new TH1F(("NMMEstimated"+IsoNames[iso_index]+theChannel+"QCDDistribution_"+(ss.str())).c_str(), ("NMMEstimated"+IsoNames[iso_index]+theChannel+"QCDDistribution_"+(ss.str())).c_str(),  2000, -1000, 1000);
        theDistributionElement.NMMEstimatedWJetsDistribution[iso_index] = new TH1F(("NMMEstimated"+IsoNames[iso_index]+theChannel+"WJetsDistribution_"+(ss.str())).c_str(), ("NMMEstimated"+IsoNames[iso_index]+theChannel+"WJetsDistribution_"+(ss.str())).c_str(),  2000, -1000, 1000);
        theDistributionElement.NMMEstimatedSignalDistribution[iso_index] = new TH1F(("NMMEstimated"+IsoNames[iso_index]+theChannel+"SignalDistribution_"+(ss.str())).c_str(), ("NMMEstimated"+IsoNames[iso_index]+theChannel+"SignalDistribution_"+(ss.str())).c_str(),  2000, -1000, 1000);
	//        theDistributionElement.NMMEstimatedQCDDistribution[iso_index]->SetBit(TH1::kCanRebin);
	//        theDistributionElement.NMMEstimatedWJetsDistribution[iso_index]->SetBit(TH1::kCanRebin);
	//        theDistributionElement.NMMEstimatedSignalDistribution[iso_index]->SetBit(TH1::kCanRebin);
       }
       theDistributions.push_back(theDistributionElement);
      }


      //Initialized only here, once! Then filled elsewhere
      theMMEstimatedValues.reserve(100000);
      for(unsigned int bin_index=0; bin_index<theNBins; bin_index++){
       struct MMEstimated theMMEstimatedValuesElement;        
       for(unsigned int iso_index=0; iso_index<3; iso_index++){
        theMMEstimatedValuesElement.NofMMEstimatedQCD[iso_index]= 0.;
        theMMEstimatedValuesElement.MMEstimatedQCDErr[iso_index]= 0.;
        theMMEstimatedValuesElement.NofMMEstimatedWJets[iso_index]= 0.;
        theMMEstimatedValuesElement.MMEstimatedWJetsErr[iso_index]= 0.;
        theMMEstimatedValuesElement.NofMMEstimatedSignal[iso_index]= 0.;
        theMMEstimatedValuesElement.MMEstimatedSignalErr[iso_index]= 0.;
       }
       theMMEstimatedValues.push_back(theMMEstimatedValuesElement);
      }

      //Initialized only here, once! Then filled elsewhere
      for(unsigned int iso_index=0; iso_index<3; iso_index++){
       theMMEstimatedPlots.MMEstimated_QCD[iso_index] = new TH1F(("MMEstimated_"+IsoNames[iso_index]+theChannel+"_QCD").c_str(), ("MMEstimated_"+IsoNames[iso_index]+theChannel+"_QCD").c_str(), theNBins, theMinBin, theMaxBin);
       theMMEstimatedPlots.MMEstimated_WJets[iso_index] = new TH1F(("MMEstimated_"+IsoNames[iso_index]+theChannel+"_WJets").c_str(), ("MMEstimated_"+IsoNames[iso_index]+theChannel+"_WJets").c_str(), theNBins, theMinBin, theMaxBin);
       theMMEstimatedPlots.MMEstimated_Signal[iso_index] = new TH1F(("MMEstimated_"+IsoNames[iso_index]+theChannel+"_Signal").c_str(), ("MMEstimated_"+IsoNames[iso_index]+theChannel+"_Signal").c_str(), theNBins, theMinBin, theMaxBin);
      }




      //Initialized only here, once! Then filled elsewhere
      theMMExpectedPlots.reserve(100000);
      for(unsigned int dataset_index=0; dataset_index < datasets.size(); dataset_index++){
       if(datasets[dataset_index].Name()=="TTbar") // Create histos for Dilept and Semilept parts
       {
         Dataset dilept; dilept.SetName("TTbarSignal");
         datasets.push_back(dilept);
         Dataset semilept; semilept.SetName("TTbarSemileptonic");
         datasets.push_back(semilept);
       }
       struct MMExpectedPlots tmpMMExpectedPlot;
       for(unsigned int iso_index=0; iso_index<3; iso_index++){
 	 tmpMMExpectedPlot.Name[iso_index] = datasets[dataset_index].Name();
         TH1F *h = new TH1F(("MMExpected_"+IsoNames[iso_index]+theChannel+"_"+datasets[dataset_index].Name()).c_str(),("MMExpected_"+IsoNames[iso_index]+theChannel+"_"+datasets[dataset_index].Name()).c_str(), theNBins, theMinBin, theMaxBin); 
         h->Sumw2();
         tmpMMExpectedPlot.MMExpected[iso_index] = h;
         h = 0;
         delete h;
       }
       theMMExpectedPlots.push_back(tmpMMExpectedPlot);
      }


      for(unsigned int iso_index=0; iso_index<3; iso_index++){
       histoN[iso_index] = new TH1F(("histoN"+IsoNames[iso_index]+theChannel).c_str(), ("histoN"+IsoNames[iso_index]+theChannel).c_str(), theNBins, theMinBin, theMaxBin);
       histoSignal[iso_index] = new TH1F(("histoSignal"+IsoNames[iso_index]+theChannel).c_str(), ("histoSignal"+IsoNames[iso_index]+theChannel).c_str(), theNBins, theMinBin, theMaxBin);
       histoW[iso_index] = new TH1F(("histoW"+IsoNames[iso_index]+theChannel).c_str(), ("histoW"+IsoNames[iso_index]+theChannel).c_str(), theNBins, theMinBin, theMaxBin);
       histoQCD[iso_index] = new TH1F(("histoQCD"+IsoNames[iso_index]+theChannel).c_str(), ("histoQCD"+IsoNames[iso_index]+theChannel).c_str(), theNBins, theMinBin, theMaxBin);

      }


}

MMEstimation::~MMEstimation(){

      for(unsigned int bin_index=0; bin_index<theNBins; bin_index++){
        for(unsigned int iso_index=0; iso_index<3; iso_index++){
	 delete theDistributions[bin_index].NMMEstimatedQCDDistribution[iso_index];
  	 delete theDistributions[bin_index].NMMEstimatedWJetsDistribution[iso_index];
 	 delete theDistributions[bin_index].NMMEstimatedSignalDistribution[iso_index];
       }
      }


      for(unsigned int iso_index=0; iso_index<3; iso_index++){
       delete theMMEstimatedPlots.MMEstimated_QCD[iso_index];
       delete theMMEstimatedPlots.MMEstimated_WJets[iso_index];
       delete theMMEstimatedPlots.MMEstimated_Signal[iso_index];
      }

      for(unsigned int iso_index=0; iso_index<3; iso_index++){
       for(unsigned int i=0;i< theMMExpectedPlots.size(); i++){
	 delete theMMExpectedPlots[i].MMExpected[iso_index];
       }
      }


      for(unsigned int iso_index=0; iso_index<3; iso_index++){
       delete histoN[iso_index];
       delete histoSignal[iso_index];
       delete histoW[iso_index];
       delete histoQCD[iso_index];
      }


}


struct Isolations MMEstimation::GetIsolations(){
  return theIsolations;
}


void MMEstimation::IncrementNSelBin(unsigned int iso_index, unsigned int bin_index, float weight){
	 theNSelected[bin_index].NSel[iso_index] += weight;
}

void MMEstimation::CountNSel(const SSDiLeptonSelection & sel, Dataset dataset, string CandTypeRequested , float weight, int selStepCut, const NTMonteCarlo* mc){

    int selStepMM[3]= {-1, -1, -1};
    string CandTypeForMM[3]= {"", "", ""};
    unsigned int njets[3] = {-1, -1, -1};
    bool passSelection[3] = {false, false, false};
    
    // Apply selections
    for (unsigned int iso_index = 0; iso_index < 3; iso_index++) {
        vector<NTJet> SelectedJetsForMM; bool LepPairForMM[3]; vector<NTElectron> candElecForMM; vector<NTMuon> candMuonForMM;
        SSDiLeptonSelection sel_aux(sel);
        selStepMM[iso_index] = sel_aux.doFullSelection(&(dataset), CandTypeRequested, false, true, GetIsolations().iso1[iso_index],GetIsolations().iso2[iso_index] );         
	// Be careful if inverted or not (at the moment pairing is applied before isolation)
        LepPairForMM[iso_index] = sel_aux.GetLeptonPairForMM( candMuonForMM, candElecForMM, CandTypeForMM[iso_index], GetIsolations().iso1[iso_index],GetIsolations().iso2[iso_index]);
        SelectedJetsForMM = sel_aux.GetSelectedJets(candMuonForMM, candElecForMM);
	njets[iso_index] = SelectedJetsForMM.size();//  if(iso_index == 2){njets[2]=SelectedJetsForMM.size(); njets[1]=SelectedJetsForMM.size(); njets[0]=SelectedJetsForMM.size();}
    }


   // Fill Histos and do checks
   for(unsigned int bin_index = 0; bin_index < theNBins; bin_index++){
      //Here, njets[] is the variable vs. whih the plots of Matrix Method will be done. It could be any other variable!     
     for(unsigned int iso_index = 0; iso_index < 3; iso_index++){
       if((CandTypeForMM[iso_index] == CandTypeRequested) && (selStepMM[iso_index] >= selStepCut)  
            && ((njets[iso_index] >= (((theMaxBin-theMinBin)/theNBins)*bin_index)) && (njets[iso_index] < (((theMaxBin-theMinBin)/theNBins)*(bin_index+1))))){
         IncrementNSelBin(iso_index, bin_index, weight);
         passSelection[iso_index]=true;
       }
     }

     if( passSelection[2] && (!passSelection[1] || !passSelection[0])) 
     std::cout << "WARNING in CountNSel: Passing Selection Tight, but not passing selection Medium or Loose!" << std::endl;

     if( passSelection[1] && !passSelection[0]) 
     std::cout << "WARNING in CountNSel: Passing Selection Medium, but not passing selection Loose!" << std::endl;
   }

   
   // Fill MC truth
   if(mc)
   if(mc->partonFlavor.first!=0 && mc->partonFlavor.second!=0){ // Only for MC. Partons set to 0 for data
   
     for(unsigned int iso_index = 0; iso_index < 3; iso_index++){
       //Here, njets[] is the variable vs. whih the plots of Matrix Method will be done. It could be any other variable!     

       if((CandTypeForMM[iso_index] == CandTypeRequested) && (selStepMM[iso_index] >= selStepCut)){
        
	 if(CandTypeRequested=="mumu_ss"){
	   if(dataset.Name() == "TTbar" && (mc->TMEME == 20 || mc->TMEME == 11010 || mc->TMEME == 22000 || // mumu
                   mc->TMEME == 10 || // l(=mu)jets
                   mc->TMEME == 11 || mc->TMEME == 11001 || mc->TMEME == 10110 || mc->TMEME == 21100 || // emu
                   mc->TMEME == 10010 || // l(=mu)tau(->had)
		   mc->TMEME == 11000) ) // tau(->mu)jets
	     FillMMExpectedPlot("TTbarSemileptonic",iso_index,njets[iso_index],weight);
	   else FillMMExpectedPlot(dataset.Name(),iso_index,njets[iso_index],weight);
	 }
	 else {
	
         if(((CandTypeRequested.find("ee")!=string::npos && (mc->TMEME == 2 || mc->TMEME == 10101 || mc->TMEME == 20200)) || // ee
	     (CandTypeRequested.find("mumu")!=string::npos && (mc->TMEME == 20 || mc->TMEME == 11010 || mc->TMEME == 22000)) // mumu
	    ) && dataset.Name() == "TTbar"){
	   FillMMExpectedPlot("TTbarSignal",iso_index,njets[iso_index],weight);
	 }else if(((CandTypeRequested.find("ee")!=string::npos && (mc->TMEME == 1 || // l(=e)jets
                   mc->TMEME == 11 || mc->TMEME == 11001 || mc->TMEME == 10110 || mc->TMEME == 21100 || // emu
                   mc->TMEME == 10001 || // l(=e)tau(->had)
		   mc->TMEME == 10100)) || // tau(->e)jets
                   (CandTypeRequested.find("mumu")!=string::npos && (mc->TMEME == 10 || // l(=mu)jets
                   mc->TMEME == 11 || mc->TMEME == 11001 || mc->TMEME == 10110 || mc->TMEME == 21100 || // emu
                   mc->TMEME == 10010 || // l(=mu)tau(->had)
		   mc->TMEME == 11000)) // tau(->mu)jets
		   ) && dataset.Name() == "TTbar"
	 ){
	   FillMMExpectedPlot("TTbarSemileptonic",iso_index,njets[iso_index],weight);
 	 }else{
	   FillMMExpectedPlot(dataset.Name(),iso_index,njets[iso_index],weight);
	 }
	 
	 }

       
       }  
     } // End of loop on iso_index

   } // End of lines applicables only on MC

}



void MMEstimation::IncludeSystematicError(){
 float EpsilonFake0 = EpsilonFake;
 float EpsilonSignal0  = EpsilonSignal;
 EpsilonFake = ranEpsFake.Gaus(EpsilonFake0, EpsilonFakeErr);
 EpsilonSignal  = ranEpsSignal.Gaus(EpsilonSignal0,  EpsilonSignalErr);
}


void MMEstimation::IncludeStatisticalError(){
 float N3 = locNSelected[2];
 float N2 = locNSelected[1] - N3;
 float N1 = locNSelected[0] - N2 - N3;
 N1 = ranN1.Poisson(N1);
 N2 = ranN2.Poisson(N2);
 N3 = ranN3.Poisson(N3);
 locNSelected[2] = N3;
 locNSelected[1] = N2+N3;
 locNSelected[0] = N1+N2+N3;
}


void MMEstimation::IncludeStatisticalErrorForEpsilonsTest(float weight){

 float N3Signal = locNSelectedSignal[2];
 float N2Signal = locNSelectedSignal[1] - N3Signal;
 float N1Signal = locNSelectedSignal[0] - N2Signal - N3Signal;
 N1Signal = weight * ranN1.Poisson(N1Signal);
 N2Signal = weight * ranN2.Poisson(N2Signal);
 N3Signal = weight * ranN3.Poisson(N3Signal);
 locNSelectedSignal[2] = N3Signal;
 locNSelectedSignal[1] = N2Signal+N3Signal;
 locNSelectedSignal[0] = N1Signal+N2Signal+N3Signal;


 float N3W = locNSelectedW[2];
 float N2W = locNSelectedW[1] - N3W;
 float N1W = locNSelectedW[0] - N2W - N3W;
 N1W = weight * ranN1.Poisson(N1W);
 N2W = weight * ranN2.Poisson(N2W);
 N3W = weight * ranN3.Poisson(N3W);
 locNSelectedW[2] = N3W;
 locNSelectedW[1] = N2W+N3W;
 locNSelectedW[0] = N1W+N2W+N3W;


 float N3QCD = locNSelectedQCD[2];
 float N2QCD = locNSelectedQCD[1] - N3QCD;
 float N1QCD = locNSelectedQCD[0] - N2QCD - N3QCD;
 N1QCD = weight * ranN1.Poisson(N1QCD);
 N2QCD = weight * ranN2.Poisson(N2QCD);
 N3QCD = weight * ranN3.Poisson(N3QCD);
 locNSelectedQCD[2] = N3QCD;
 locNSelectedQCD[1] = N2QCD+N3QCD;
 locNSelectedQCD[0] = N1QCD+N2QCD+N3QCD;




}




void MMEstimation::SolveTheSystem(bool doCorrections){

      for(unsigned int iso_index=0; iso_index<3; iso_index++){
	NMMEstimatedSignal[iso_index] = 0;
	NMMEstimatedWJets[iso_index] = 0;
	NMMEstimatedQCD[iso_index] = 0;
      }


      TMatrixF theEffmatrix(3,3); // A
      TMatrixF solutions(3,1);    // x
      TMatrixF mesured(3,1);      // b

      float eff[9] ;  // A


      float Eff_ltS   = EpsilonSignal*EpsilonSignal;
      float Eff_ltW   = EpsilonSignal*EpsilonFake;
      float Eff_ltQCD = EpsilonFake*EpsilonFake;
      
      float Eff_lmS   = 2.*EpsilonSignal - EpsilonSignal*EpsilonSignal;
      float Eff_lmW   = EpsilonSignal + EpsilonFake - EpsilonSignal*EpsilonFake ;
      float Eff_lmQCD = 2.*EpsilonFake - EpsilonFake*EpsilonFake;

      eff[0] = Eff_ltS  ;
      eff[1] = Eff_ltW  ;
      eff[2] = Eff_ltQCD ;
      //  eff[2] = EpsilonFake0*EpsilonFake0;
      
      eff[3] = Eff_lmS  ;
      eff[4] = Eff_lmW  ;
      eff[5] = Eff_lmQCD ;
      //  eff[5] = 2.*EpsilonFake0 - EpsilonFake0*EpsilonFake0;
      
      eff[6] = 1.;
      eff[7] = 1.;
      eff[8] = 1.; //1 
      
      theEffmatrix.SetMatrixArray(eff);
      theEffmatrix.Invert();

      float ArrayMes[3] ; // b
      ArrayMes[0] = locNSelected[2];
      ArrayMes[1] = locNSelected[1];
      ArrayMes[2] = locNSelected[0];

      mesured.SetMatrixArray(ArrayMes);
      solutions = theEffmatrix*mesured; 

      NMMEstimatedSignal[2] = solutions(0,0)*Eff_ltS;
      NMMEstimatedWJets[2] = solutions(1,0)*Eff_ltW;
      NMMEstimatedQCD[2] = solutions(2,0)*Eff_ltQCD;
     
      NMMEstimatedSignal[1] = solutions(0,0)*Eff_lmS;
      NMMEstimatedWJets[1] = solutions(1,0)*Eff_lmW;
      NMMEstimatedQCD[1] = solutions(2,0)*Eff_lmQCD;
      
      NMMEstimatedSignal[0] = solutions(0,0);
      NMMEstimatedWJets[0] = solutions(1,0);
      NMMEstimatedQCD[0] = solutions(2,0);      

      if(doCorrections){
        float EpsilonSignal_bar = 1 - EpsilonSignal;

        NMMEstimatedSignal[2] = NMMEstimatedSignal[2] 
                              + NMMEstimatedSignal[0] * (2 * EpsilonSignal_bar - 2 * EpsilonSignal_bar * EpsilonSignal_bar) * EpsilonFake
                              + NMMEstimatedSignal[0] * EpsilonSignal_bar * EpsilonSignal_bar * EpsilonFake * EpsilonFake;

        NMMEstimatedWJets[2] = NMMEstimatedWJets[2] 
                             - NMMEstimatedSignal[0] * (2 * EpsilonSignal_bar - 2 * EpsilonSignal_bar * EpsilonSignal_bar) * EpsilonFake
                             + NMMEstimatedWJets[0] * EpsilonSignal_bar * EpsilonFake * EpsilonFake;


        NMMEstimatedQCD[2] = NMMEstimatedQCD[2]
	                   - NMMEstimatedSignal[0] * EpsilonSignal_bar * EpsilonSignal_bar * EpsilonFake * EpsilonFake
                           - NMMEstimatedWJets[0] * EpsilonSignal_bar * EpsilonFake * EpsilonFake;
      }



}

void MMEstimation::FillDistributions(unsigned int bin_index){

    for(unsigned int iso_index=0; iso_index<3; iso_index++){
          theDistributions[bin_index].NMMEstimatedSignalDistribution[iso_index]->Fill(NMMEstimatedSignal[iso_index]);
          theDistributions[bin_index].NMMEstimatedWJetsDistribution[iso_index]->Fill(NMMEstimatedWJets[iso_index]);
          theDistributions[bin_index].NMMEstimatedQCDDistribution[iso_index]->Fill(NMMEstimatedQCD[iso_index]);

          theMMEstimatedValues[bin_index].NofMMEstimatedSignal[iso_index] += NMMEstimatedSignal[iso_index];
          theMMEstimatedValues[bin_index].NofMMEstimatedQCD[iso_index] += NMMEstimatedQCD[iso_index];
          theMMEstimatedValues[bin_index].NofMMEstimatedWJets[iso_index] += NMMEstimatedWJets[iso_index];

	  /*
          theMMEstimatedValues[bin_index].MMEstimatedSignalErr[iso_index] += NMMEstimatedSignal[iso_index]*NMMEstimatedSignal[iso_index];
          theMMEstimatedValues[bin_index].MMEstimatedQCDErr[iso_index] += NMMEstimatedQCD[iso_index]*NMMEstimatedQCD[iso_index];
          theMMEstimatedValues[bin_index].MMEstimatedWJetsErr[iso_index] += NMMEstimatedWJets[iso_index]*NMMEstimatedWJets[iso_index];
	  */
    }
}


void MMEstimation::CalculateRms(unsigned int bin_index, unsigned int NbIterations){

    for(unsigned int iso_index=0; iso_index<3; iso_index++){
          theMMEstimatedValues[bin_index].MMEstimatedSignalErr[iso_index] += ((theMMEstimatedValues[bin_index].NofMMEstimatedSignal[iso_index]/NbIterations)-NMMEstimatedSignal[iso_index])*((theMMEstimatedValues[bin_index].NofMMEstimatedSignal[iso_index]/NbIterations)-NMMEstimatedSignal[iso_index]);
          theMMEstimatedValues[bin_index].MMEstimatedQCDErr[iso_index] += ((theMMEstimatedValues[bin_index].NofMMEstimatedQCD[iso_index]/NbIterations)-NMMEstimatedQCD[iso_index])*((theMMEstimatedValues[bin_index].NofMMEstimatedQCD[iso_index]/NbIterations)-NMMEstimatedQCD[iso_index]);
          theMMEstimatedValues[bin_index].MMEstimatedWJetsErr[iso_index] += ((theMMEstimatedValues[bin_index].NofMMEstimatedWJets[iso_index]/NbIterations)-NMMEstimatedWJets[iso_index])*((theMMEstimatedValues[bin_index].NofMMEstimatedWJets[iso_index]/NbIterations)-NMMEstimatedWJets[iso_index]);

    }
}



void MMEstimation::SetMMEstimated(unsigned int bin_index, unsigned int NbIterations){

     for(unsigned int iso_index=0; iso_index<3; iso_index++){
       /*
      theMMEstimatedValues[bin_index].NofMMEstimatedQCD[iso_index] = theDistributions[bin_index].NMMEstimatedQCDDistribution[iso_index]->GetMean();
      theMMEstimatedValues[bin_index].MMEstimatedQCDErr[iso_index] = theDistributions[bin_index].NMMEstimatedQCDDistribution[iso_index]->GetRMS();
      theMMEstimatedValues[bin_index].NofMMEstimatedWJets[iso_index] = theDistributions[bin_index].NMMEstimatedWJetsDistribution[iso_index]->GetMean();
      theMMEstimatedValues[bin_index].MMEstimatedWJetsErr[iso_index] = theDistributions[bin_index].NMMEstimatedWJetsDistribution[iso_index]->GetRMS();
      theMMEstimatedValues[bin_index].NofMMEstimatedSignal[iso_index] = theDistributions[bin_index].NMMEstimatedSignalDistribution[iso_index]->GetMean();
      theMMEstimatedValues[bin_index].MMEstimatedSignalErr[iso_index] = theDistributions[bin_index].NMMEstimatedSignalDistribution[iso_index]->GetRMS();
      */


      theMMEstimatedValues[bin_index].NofMMEstimatedSignal[iso_index] = (theMMEstimatedValues[bin_index].NofMMEstimatedSignal[iso_index])/NbIterations;
      theMMEstimatedValues[bin_index].NofMMEstimatedQCD[iso_index] = (theMMEstimatedValues[bin_index].NofMMEstimatedQCD[iso_index])/NbIterations;
      theMMEstimatedValues[bin_index].NofMMEstimatedWJets[iso_index] = (theMMEstimatedValues[bin_index].NofMMEstimatedWJets[iso_index])/NbIterations;
      theMMEstimatedValues[bin_index].MMEstimatedSignalErr[iso_index] = sqrt((theMMEstimatedValues[bin_index].MMEstimatedSignalErr[iso_index])/(NbIterations-1));
      theMMEstimatedValues[bin_index].MMEstimatedQCDErr[iso_index] = sqrt((theMMEstimatedValues[bin_index].MMEstimatedQCDErr[iso_index])/(NbIterations-1));
      theMMEstimatedValues[bin_index].MMEstimatedWJetsErr[iso_index] = sqrt((theMMEstimatedValues[bin_index].MMEstimatedWJetsErr[iso_index])/(NbIterations-1));

     }

}


void MMEstimation::ReadMMFile(string fileName){

  TFile* file_MM   = new TFile((fileName).c_str());
  file_MM->cd();
  
  TH1F* tmpMMNSelected[3];
  for(unsigned int iso_index=0; iso_index<3; iso_index++){
   tmpMMNSelected[iso_index]  = (TH1F*)gDirectory->Get(("MMNSelected_"+IsoNames[iso_index]+theChannel).c_str());
  }
  for(unsigned int iso_index=0; iso_index<3; iso_index++){
   for(unsigned int bin_index = 0; bin_index < theNBins; bin_index++){
    theNSelected[bin_index].NSel[iso_index] = tmpMMNSelected[iso_index]->GetBinContent(bin_index+1) ;
   }
  }
  for(unsigned int iso_index=0; iso_index<3; iso_index++){
   tmpMMNSelected[iso_index] = 0;
   delete tmpMMNSelected[iso_index];
  }

  file_MM->Close();
  delete file_MM;

}



void MMEstimation::ReadMMFileForPullTest(string fileName, float epsilon_s, float epsilon_f){

  TFile* file_MM   = new TFile((fileName).c_str());
  file_MM->cd();

  TH1F* histoTmpSignal[3];
  TH1F* histoTmpW[3];
  TH1F* histoTmpQCD[3];
  
  for(unsigned int iso_index=0; iso_index<3; iso_index++){
   histoTmpSignal[iso_index] = (TH1F*)gDirectory->Get(("MMExpected_"+IsoNames[iso_index]+theChannel+"_TTbarSignal").c_str());
   histoTmpW[iso_index] = (TH1F*)gDirectory->Get(("MMExpected_"+IsoNames[iso_index]+theChannel+"_TTbarSemileptonic").c_str());
   histoTmpQCD[iso_index] = (TH1F*)gDirectory->Get(("MMExpected_"+IsoNames[iso_index]+theChannel+"_TTbar").c_str());
  }  

  for(unsigned int iso_index=0; iso_index<3; iso_index++){
   for(unsigned int bin_index = 0; bin_index < theNBins; bin_index++){
     histoSignal[iso_index]->SetBinContent(bin_index+1, histoTmpSignal[iso_index]->GetBinContent(bin_index+1));
     histoW[iso_index]->SetBinContent(bin_index+1, histoTmpW[iso_index]->GetBinContent(bin_index+1));
     histoQCD[iso_index]->SetBinContent(bin_index+1, histoTmpQCD[iso_index]->GetBinContent(bin_index+1));
   }
  }  

  histoSignal[1]->Add(histoSignal[0], histoSignal[0], 0, 2*epsilon_s - epsilon_s * epsilon_s);
  histoW[1]->Add(histoW[0], histoW[0], 0, epsilon_s + epsilon_f - epsilon_s * epsilon_f);
  histoQCD[1]->Add(histoQCD[0], histoQCD[0], 0, 2*epsilon_f - epsilon_f * epsilon_f);

  histoSignal[2]->Add(histoSignal[0], histoSignal[0], 0, epsilon_s * epsilon_s);
  histoW[2]->Add(histoW[0], histoW[0], 0, epsilon_s * epsilon_f);
  histoQCD[2]->Add(histoQCD[0], histoQCD[0], 0, epsilon_f * epsilon_f);

  histoN[0]->Add(histoW[0], histoSignal[0], 1, 1);
  histoN[0]->Add(histoN[0], histoQCD[0], 1, 1);

  histoN[1]->Add(histoW[0], histoSignal[0], epsilon_s + epsilon_f - epsilon_s * epsilon_f, 2*epsilon_s - epsilon_s * epsilon_s);
  histoN[1]->Add(histoN[1], histoQCD[0], 1, 2*epsilon_f - epsilon_f * epsilon_f);

  histoN[2]->Add(histoW[0], histoSignal[0], epsilon_s * epsilon_f, epsilon_s * epsilon_s);
  histoN[2]->Add(histoN[2], histoQCD[0], 1, epsilon_f * epsilon_f);


  for(unsigned int iso_index=0; iso_index<3; iso_index++){
   for(unsigned int bin_index = 0; bin_index < theNBins; bin_index++){
    theNSelected[bin_index].NSel[iso_index] = histoN[iso_index]->GetBinContent(bin_index+1) ;
   }
  }

 
  for(unsigned int iso_index=0; iso_index<3; iso_index++){
   histoTmpSignal[iso_index] = 0;
   delete histoTmpSignal[iso_index];
   histoTmpW[iso_index] = 0;
   delete histoTmpW[iso_index];
   histoTmpQCD[iso_index] = 0;
   delete histoTmpQCD[iso_index];
  }

  file_MM->Close();
  delete file_MM;

}





void MMEstimation::RunTheMatrixMethod(vector<struct MMEpsilons> valMMEpsilons, unsigned int NbIterations, bool doStatistical, bool doSystematic, bool doCorrections){


  for(unsigned int bin_index = 0; bin_index < theNBins; bin_index++){

   for(unsigned int iso_index=0; iso_index<3; iso_index++){
    for(unsigned int i=0; i< (unsigned int) theDistributions[bin_index].NMMEstimatedQCDDistribution[iso_index]->GetXaxis()->GetNbins(); i++){
     theDistributions[bin_index].NMMEstimatedQCDDistribution[iso_index]->SetBinContent(i+1, 0.);
    }
    for(unsigned int i=0; i< (unsigned int) theDistributions[bin_index].NMMEstimatedWJetsDistribution[iso_index]->GetXaxis()->GetNbins(); i++){
     theDistributions[bin_index].NMMEstimatedWJetsDistribution[iso_index]->SetBinContent(i+1, 0.);
    }
    for(unsigned int i=0; i< (unsigned int) theDistributions[bin_index].NMMEstimatedSignalDistribution[iso_index]->GetXaxis()->GetNbins(); i++){
     theDistributions[bin_index].NMMEstimatedSignalDistribution[iso_index]->SetBinContent(i+1, 0.);
    }
   }



   for (unsigned int i=0 ; i<NbIterations ; i++ ) {   // pseudoexperiments
      for(unsigned int iso_index=0; iso_index<3; iso_index++){
       locNSelected[iso_index] = theNSelected[bin_index].NSel[iso_index];
      }
      EpsilonFake = valMMEpsilons[bin_index].EpsilonFake;
      EpsilonFakeErr = valMMEpsilons[bin_index].EpsilonFakeErr;
      EpsilonSignal = valMMEpsilons[bin_index].EpsilonSignal;
      EpsilonSignalErr = valMMEpsilons[bin_index].EpsilonSignalErr;

      if(doStatistical) IncludeStatisticalError();
      if(doSystematic) IncludeSystematicError();
      SolveTheSystem(doCorrections);
      FillDistributions(bin_index);
   }


   for (unsigned int i=0 ; i<NbIterations ; i++ ) {   // pseudoexperiments
      for(unsigned int iso_index=0; iso_index<3; iso_index++){
       locNSelected[iso_index] = theNSelected[bin_index].NSel[iso_index];
      }
      EpsilonFake = valMMEpsilons[bin_index].EpsilonFake;
      EpsilonFakeErr = valMMEpsilons[bin_index].EpsilonFakeErr;
      EpsilonSignal = valMMEpsilons[bin_index].EpsilonSignal;
      EpsilonSignalErr = valMMEpsilons[bin_index].EpsilonSignalErr;

      if(doStatistical) IncludeStatisticalError();
      if(doSystematic) IncludeSystematicError();
      SolveTheSystem(doCorrections);
      CalculateRms(bin_index, NbIterations); 
   }

 
   SetMMEstimated(bin_index, NbIterations);

   for(unsigned int iso_index=0; iso_index<3; iso_index++){
     /*
    theMMEstimatedPlots.MMEstimated_QCD[iso_index]->SetBinContent(bin_index+1, theDistributions[bin_index].NMMEstimatedQCDDistribution[iso_index]->GetMean());
    theMMEstimatedPlots.MMEstimated_WJets[iso_index]->SetBinContent(bin_index+1, theDistributions[bin_index].NMMEstimatedWJetsDistribution[iso_index]->GetMean());
    theMMEstimatedPlots.MMEstimated_Signal[iso_index]->SetBinContent(bin_index+1, theDistributions[bin_index].NMMEstimatedSignalDistribution[iso_index]->GetMean());
    theMMEstimatedPlots.MMEstimated_QCD[iso_index]->SetBinError(bin_index+1, theDistributions[bin_index].NMMEstimatedQCDDistribution[iso_index]->GetRMS());
    theMMEstimatedPlots.MMEstimated_WJets[iso_index]->SetBinError(bin_index+1, theDistributions[bin_index].NMMEstimatedWJetsDistribution[iso_index]->GetRMS());
    theMMEstimatedPlots.MMEstimated_Signal[iso_index]->SetBinError(bin_index+1, theDistributions[bin_index].NMMEstimatedSignalDistribution[iso_index]->GetRMS());
     */

    theMMEstimatedPlots.MMEstimated_QCD[iso_index]->SetBinContent(bin_index+1, theMMEstimatedValues[bin_index].NofMMEstimatedQCD[iso_index] );
    theMMEstimatedPlots.MMEstimated_WJets[iso_index]->SetBinContent(bin_index+1, theMMEstimatedValues[bin_index].NofMMEstimatedWJets[iso_index] );
    theMMEstimatedPlots.MMEstimated_Signal[iso_index]->SetBinContent(bin_index+1, theMMEstimatedValues[bin_index].NofMMEstimatedSignal[iso_index] );
    theMMEstimatedPlots.MMEstimated_QCD[iso_index]->SetBinError(bin_index+1, theMMEstimatedValues[bin_index].MMEstimatedQCDErr[iso_index] );
    theMMEstimatedPlots.MMEstimated_WJets[iso_index]->SetBinError(bin_index+1, theMMEstimatedValues[bin_index].MMEstimatedWJetsErr[iso_index] );
    theMMEstimatedPlots.MMEstimated_Signal[iso_index]->SetBinError(bin_index+1, theMMEstimatedValues[bin_index].MMEstimatedSignalErr[iso_index] );


   }

  }



}









void MMEstimation::RunTheMatrixMethodForEpsilonsTest(vector<struct MMEpsilons> valMMEpsilons, unsigned int NbIterations, bool doStatistical, bool doSystematic, float weight){


  for(unsigned int bin_index = 0; bin_index < theNBins; bin_index++){

   for(unsigned int iso_index=0; iso_index<3; iso_index++){
    for(unsigned int i=0; i< (unsigned int) theDistributions[bin_index].NMMEstimatedQCDDistribution[iso_index]->GetXaxis()->GetNbins(); i++){
     theDistributions[bin_index].NMMEstimatedQCDDistribution[iso_index]->SetBinContent(i+1, 0.);
    }
    for(unsigned int i=0; i< (unsigned int) theDistributions[bin_index].NMMEstimatedWJetsDistribution[iso_index]->GetXaxis()->GetNbins(); i++){
     theDistributions[bin_index].NMMEstimatedWJetsDistribution[iso_index]->SetBinContent(i+1, 0.);
    }
    for(unsigned int i=0; i< (unsigned int) theDistributions[bin_index].NMMEstimatedSignalDistribution[iso_index]->GetXaxis()->GetNbins(); i++){
     theDistributions[bin_index].NMMEstimatedSignalDistribution[iso_index]->SetBinContent(i+1, 0.);
    }
   }



   for (unsigned int i=0 ; i<NbIterations ; i++ ) {   // pseudoexperiments
      for(unsigned int iso_index=0; iso_index<3; iso_index++){
	locNSelectedSignal[iso_index] = histoSignal[iso_index]->GetBinContent(bin_index+1);
	locNSelectedW[iso_index] = histoW[iso_index]->GetBinContent(bin_index+1);
	locNSelectedQCD[iso_index] = histoQCD[iso_index]->GetBinContent(bin_index+1);
      }

      if(doStatistical) IncludeStatisticalErrorForEpsilonsTest(weight);
      if(doSystematic) IncludeSystematicError();

      for(unsigned int iso_index=0; iso_index<3; iso_index++){
        locNSelected[iso_index] = locNSelectedSignal[iso_index]
                                + locNSelectedW[iso_index]
                                + locNSelectedQCD[iso_index];
      }

      if (locNSelectedSignal[0] != 0) {
	EpsilonSignal = sqrt(locNSelectedSignal[2]/locNSelectedSignal[0]);
      }else{
        EpsilonSignal = 0;
      }
      EpsilonSignalErr = 0;
      if (locNSelectedW[0] != 0 && EpsilonSignal != 0) {
       EpsilonFake = (locNSelectedW[2]/locNSelectedW[0])/(EpsilonSignal);
      }else{
        EpsilonFake = 0;
      }
      EpsilonFakeErr = 0;

      SolveTheSystem(false);
      FillDistributions(bin_index);
   }
      


   for (unsigned int i=0 ; i<NbIterations ; i++ ) {   // pseudoexperiments
      for(unsigned int iso_index=0; iso_index<3; iso_index++){
	locNSelectedSignal[iso_index] = histoSignal[iso_index]->GetBinContent(bin_index+1);
	locNSelectedW[iso_index] = histoW[iso_index]->GetBinContent(bin_index+1);
	locNSelectedQCD[iso_index] = histoQCD[iso_index]->GetBinContent(bin_index+1);
      }

      if(doStatistical) IncludeStatisticalErrorForEpsilonsTest(weight);
      if(doSystematic) IncludeSystematicError();

      for(unsigned int iso_index=0; iso_index<3; iso_index++){
        locNSelected[iso_index] = locNSelectedSignal[iso_index]
                                + locNSelectedW[iso_index]
                                + locNSelectedQCD[iso_index];
      }

      if (locNSelectedSignal[0] != 0) {
	EpsilonSignal = sqrt(locNSelectedSignal[2]/locNSelectedSignal[0]);
      }else{
        EpsilonSignal = 0;
      }
      EpsilonSignalErr = 0;
      if (locNSelectedW[0] != 0 && EpsilonSignal != 0) {
       EpsilonFake = (locNSelectedW[2]/locNSelectedW[0])/(EpsilonSignal);
      }else{
        EpsilonFake = 0;
      }
      EpsilonFakeErr = 0;

      SolveTheSystem(false);
      CalculateRms(bin_index, NbIterations); 
   }

 
   SetMMEstimated(bin_index, NbIterations);

   for(unsigned int iso_index=0; iso_index<3; iso_index++){
     /*
    theMMEstimatedPlots.MMEstimated_QCD[iso_index]->SetBinContent(bin_index+1, theDistributions[bin_index].NMMEstimatedQCDDistribution[iso_index]->GetMean());
    theMMEstimatedPlots.MMEstimated_WJets[iso_index]->SetBinContent(bin_index+1, theDistributions[bin_index].NMMEstimatedWJetsDistribution[iso_index]->GetMean());
    theMMEstimatedPlots.MMEstimated_Signal[iso_index]->SetBinContent(bin_index+1, theDistributions[bin_index].NMMEstimatedSignalDistribution[iso_index]->GetMean());
    theMMEstimatedPlots.MMEstimated_QCD[iso_index]->SetBinError(bin_index+1, theDistributions[bin_index].NMMEstimatedQCDDistribution[iso_index]->GetRMS());
    theMMEstimatedPlots.MMEstimated_WJets[iso_index]->SetBinError(bin_index+1, theDistributions[bin_index].NMMEstimatedWJetsDistribution[iso_index]->GetRMS());
    theMMEstimatedPlots.MMEstimated_Signal[iso_index]->SetBinError(bin_index+1, theDistributions[bin_index].NMMEstimatedSignalDistribution[iso_index]->GetRMS());
     */

    theMMEstimatedPlots.MMEstimated_QCD[iso_index]->SetBinContent(bin_index+1, theMMEstimatedValues[bin_index].NofMMEstimatedQCD[iso_index] );
    theMMEstimatedPlots.MMEstimated_WJets[iso_index]->SetBinContent(bin_index+1, theMMEstimatedValues[bin_index].NofMMEstimatedWJets[iso_index] );
    theMMEstimatedPlots.MMEstimated_Signal[iso_index]->SetBinContent(bin_index+1, theMMEstimatedValues[bin_index].NofMMEstimatedSignal[iso_index] );
    theMMEstimatedPlots.MMEstimated_QCD[iso_index]->SetBinError(bin_index+1, theMMEstimatedValues[bin_index].MMEstimatedQCDErr[iso_index] );
    theMMEstimatedPlots.MMEstimated_WJets[iso_index]->SetBinError(bin_index+1, theMMEstimatedValues[bin_index].MMEstimatedWJetsErr[iso_index] );
    theMMEstimatedPlots.MMEstimated_Signal[iso_index]->SetBinError(bin_index+1, theMMEstimatedValues[bin_index].MMEstimatedSignalErr[iso_index] );


   }

  }



}

















vector<struct MMEstimated> MMEstimation::GetMMEstimated(){
return theMMEstimatedValues;
}


void MMEstimation::PrintMMEstimated(){

/*
   for(unsigned int bin_index=0; bin_index<theNBins; bin_index++){
    for(unsigned int iso_index=0; iso_index<3; iso_index++){
      std::cout <<"theMMEstimatedValue of NofQDC for bin_index: "<<bin_index <<" and iso_index: "<<iso_index <<" is: "<<theMMEstimatedValues[bin_index].NofMMEstimatedQCD[iso_index] << std::endl;
      std::cout <<"theMMEstimatedValue of QCDErr for bin_index: "<<bin_index <<" and iso_index: "<<iso_index <<" is: "<<theMMEstimatedValues[bin_index].MMEstimatedQCDErr[iso_index] << std::endl;
      std::cout <<"theMMEstimatedValue of NofWJets for bin_index: "<<bin_index <<" and iso_index: "<<iso_index <<" is: "<<theMMEstimatedValues[bin_index].NofMMEstimatedWJets[iso_index] << std::endl;
      std::cout <<"theMMEstimatedValue of WJetsErr for bin_index: "<<bin_index <<" and iso_index: "<<iso_index <<" is: "<<theMMEstimatedValues[bin_index].MMEstimatedWJetsErr[iso_index] << std::endl;
      std::cout <<"theMMEstimatedValue of NofSignal for bin_index: "<<bin_index <<" and iso_index: "<<iso_index <<" is: "<<theMMEstimatedValues[bin_index].NofMMEstimatedSignal[iso_index] << std::endl;
      std::cout <<"theMMEstimatedValue of SignalErr for bin_index: "<<bin_index <<" and iso_index: "<<iso_index <<" is: "<<theMMEstimatedValues[bin_index].MMEstimatedSignalErr[iso_index] << std::endl;
    }
   }
*/


   std::cout << "\\begin{table}" << std::endl;
   std::cout << "\\begin{center}" << std::endl;
   std::cout << "\\begin{tabular}{|";
   for(unsigned int bin_index=2; bin_index<theNBins; bin_index++){
    std::cout << " c |";
   }
   std::cout << "}" << std::endl;
   std::cout << "\\hline" << std::endl;


   for(unsigned int bin_index=2; bin_index<theNBins; bin_index++){
    std::cout << " & Njets = " << bin_index;
   }
   std::cout << "\\\\ \\hline" << std::endl;


   std::cout << "N of Signal-like";
   for(unsigned int bin_index=2; bin_index<theNBins; bin_index++){
    std::cout << " & " << theMMEstimatedValues[bin_index].NofMMEstimatedSignal[2] << "$\\pm$" << theMMEstimatedValues[bin_index].MMEstimatedSignalErr[2];
   }
   std::cout << "\\\\ \\hline" << std::endl;


   std::cout << "N of W-like";
   for(unsigned int bin_index=2; bin_index<theNBins; bin_index++){
    std::cout << " & " << theMMEstimatedValues[bin_index].NofMMEstimatedWJets[2] << "$\\pm$" << theMMEstimatedValues[bin_index].MMEstimatedWJetsErr[2];
   }
   std::cout << "\\\\ \\hline" << std::endl;


   std::cout << "N of QCD-like";
   for(unsigned int bin_index=2; bin_index<theNBins; bin_index++){
    std::cout << " & " << theMMEstimatedValues[bin_index].NofMMEstimatedQCD[2] << "$\\pm$" << theMMEstimatedValues[bin_index].MMEstimatedQCDErr[2];
   }
   std::cout << "\\\\ \\hline" << std::endl;

   std::cout << "\\end{tabular}" << std::endl;
   std::cout << "\\caption{FIX ME} \\label{FIX ME}" << std::endl;
   std::cout << "\\end{center}" << std::endl;
   std::cout << "\\end{table}" << std::endl;




   float NofMMEstimatedSignalTOTAL = 0;
   float NofMMEstimatedWJetsTOTAL = 0;
   float NofMMEstimatedQCDTOTAL = 0;

   float MMEstimatedSignalErrTOTAL = 0;
   float MMEstimatedWJetsErrTOTAL = 0;
   float MMEstimatedQCDErrTOTAL = 0;

   for(unsigned int bin_index=0; bin_index<theNBins; bin_index++){

    NofMMEstimatedSignalTOTAL += theMMEstimatedValues[bin_index].NofMMEstimatedSignal[2];
    NofMMEstimatedQCDTOTAL += theMMEstimatedValues[bin_index].NofMMEstimatedQCD[2];
    NofMMEstimatedWJetsTOTAL += theMMEstimatedValues[bin_index].NofMMEstimatedWJets[2];
    MMEstimatedSignalErrTOTAL += theMMEstimatedValues[bin_index].MMEstimatedSignalErr[2]*theMMEstimatedValues[bin_index].MMEstimatedSignalErr[2];
    MMEstimatedQCDErrTOTAL += theMMEstimatedValues[bin_index].MMEstimatedQCDErr[2]*theMMEstimatedValues[bin_index].MMEstimatedQCDErr[2];
    MMEstimatedWJetsErrTOTAL += theMMEstimatedValues[bin_index].MMEstimatedWJetsErr[2]*theMMEstimatedValues[bin_index].MMEstimatedWJetsErr[2];

   }

   std::cout << "NofMMEstimatedSignalTOTAL: " << NofMMEstimatedSignalTOTAL << std::endl;
   std::cout << "MMEstimatedSignalErrTOTAL: " << sqrt(MMEstimatedSignalErrTOTAL) << std::endl;
   std::cout << "NofMMEstimatedQCDTOTAL: " << NofMMEstimatedQCDTOTAL << std::endl;
   std::cout << "MMEstimatedQCDErrTOTAL: " << sqrt(MMEstimatedQCDErrTOTAL) << std::endl;
   std::cout << "NofMMEstimatedWJetsTOTAL: " << NofMMEstimatedWJetsTOTAL << std::endl;
   std::cout << "MMEstimatedWJetsErrTOTAL: " << sqrt(MMEstimatedWJetsErrTOTAL) << std::endl;




}


struct MMEstimatedPlots MMEstimation::GetMMEstimatedPlots(){
 return theMMEstimatedPlots;
}


void MMEstimation::FillMMExpectedPlot(string datasetname, unsigned int iso_index, float val, float weight){

       for(unsigned int i=0;i< theMMExpectedPlots.size(); i++){
	 if(datasetname == theMMExpectedPlots[i].Name[iso_index]) (theMMExpectedPlots[i].MMExpected[iso_index])->Fill(val, weight);
       }

}

vector <struct MMExpectedPlots> MMEstimation::GetMMExpectedPlots(){
  return theMMExpectedPlots;
}



void MMEstimation::WriteForProof(){

    for(unsigned int iso_index=0; iso_index<3; iso_index++){
       for(unsigned int i=0;i< theMMExpectedPlots.size(); i++){
	 (theMMExpectedPlots[i].MMExpected[iso_index])->Write();
       }
    }

    // Writing histos contaning all information concerning NSelected
    TH1F* MMNSelected[3];

    for(unsigned int iso_index=0; iso_index<3; iso_index++){
      MMNSelected[iso_index] = new TH1F(("MMNSelected_"+IsoNames[iso_index]+theChannel).c_str(), ("MMNSelected_"+IsoNames[iso_index]+theChannel).c_str(),  theNBins, theMinBin, theMaxBin);
    }
    for(unsigned int iso_index=0; iso_index<3; iso_index++){
     for(unsigned int bin_index = 0; bin_index < theNBins; bin_index++){
       MMNSelected[iso_index]->SetBinContent(bin_index+1, theNSelected[bin_index].NSel[iso_index]);
     }
     MMNSelected[iso_index]->Write();
    }
    for(unsigned int iso_index=0; iso_index<3; iso_index++){
     delete MMNSelected[iso_index];
    }
}


void MMEstimation::WriteMMFile(string fileName){
    
    TFile* file_MM = new TFile((fileName).c_str(),"RECREATE");
    file_MM->cd();


    for(unsigned int bin_index=0; bin_index<theNBins; bin_index++){
      for(unsigned int iso_index=0; iso_index<3; iso_index++){
	theDistributions[bin_index].NMMEstimatedQCDDistribution[iso_index]->Write();
  	theDistributions[bin_index].NMMEstimatedWJetsDistribution[iso_index]->Write();
 	theDistributions[bin_index].NMMEstimatedSignalDistribution[iso_index]->Write();
      }
    }



    for(unsigned int iso_index=0; iso_index<3; iso_index++){
       for(unsigned int i=0;i< theMMExpectedPlots.size(); i++){
	 (theMMExpectedPlots[i].MMExpected[iso_index])->Write();
       }
    }

    for(unsigned int iso_index=0; iso_index<3; iso_index++){
      theMMEstimatedPlots.MMEstimated_Signal[iso_index]->Write();
      theMMEstimatedPlots.MMEstimated_QCD[iso_index]->Write();
      theMMEstimatedPlots.MMEstimated_WJets[iso_index]->Write();
    }

    // Writing histos contaning all information concerning NSelected
    TH1F* MMNSelected[3];

    for(unsigned int iso_index=0; iso_index<3; iso_index++){
      MMNSelected[iso_index] = new TH1F(("MMNSelected_"+IsoNames[iso_index]+theChannel).c_str(), ("MMNSelected_"+IsoNames[iso_index]+theChannel).c_str(),  theNBins, theMinBin, theMaxBin);
    }
    for(unsigned int iso_index=0; iso_index<3; iso_index++){
     for(unsigned int bin_index = 0; bin_index < theNBins; bin_index++){
       MMNSelected[iso_index]->SetBinContent(bin_index+1, theNSelected[bin_index].NSel[iso_index]);
     }
     MMNSelected[iso_index]->Write();
    }
    for(unsigned int iso_index=0; iso_index<3; iso_index++){
     MMNSelected[iso_index] = 0;
     delete MMNSelected[iso_index];
    }

 
    file_MM->Write();
    file_MM->Close();
    delete file_MM;


}


void MMEstimation::WriteMMFileFast(string fileName){

    TFile* file_MM = new TFile((fileName).c_str(),"RECREATE");
    file_MM->cd();

    for(unsigned int bin_index=0; bin_index<theNBins; bin_index++){
      for(unsigned int iso_index=0; iso_index<3; iso_index++){
	theDistributions[bin_index].NMMEstimatedQCDDistribution[iso_index]->Write();
  	theDistributions[bin_index].NMMEstimatedWJetsDistribution[iso_index]->Write();
 	theDistributions[bin_index].NMMEstimatedSignalDistribution[iso_index]->Write();
      }
    }

    for(unsigned int iso_index=0; iso_index<3; iso_index++){
      theMMEstimatedPlots.MMEstimated_Signal[iso_index]->Write();
      theMMEstimatedPlots.MMEstimated_QCD[iso_index]->Write();
      theMMEstimatedPlots.MMEstimated_WJets[iso_index]->Write();
    }

    // Writing histos contaning all information concerning NSelected
    TH1F* MMNSelected[3];
    for(unsigned int iso_index=0; iso_index<3; iso_index++){
      MMNSelected[iso_index] = new TH1F(("MMNSelected_"+IsoNames[iso_index]+theChannel).c_str(), ("MMNSelected_"+IsoNames[iso_index]+theChannel).c_str(),  theNBins, theMinBin, theMaxBin);
    }
    for(unsigned int iso_index=0; iso_index<3; iso_index++){
     for(unsigned int bin_index = 0; bin_index < theNBins; bin_index++){
       MMNSelected[iso_index]->SetBinContent(bin_index+1, theNSelected[bin_index].NSel[iso_index]);
     }
     MMNSelected[iso_index]->Write();
    }
    for(unsigned int iso_index=0; iso_index<3; iso_index++){
     MMNSelected[iso_index] = 0;
     delete MMNSelected[iso_index];
    }

    file_MM->Write();
    file_MM->Close();
    delete file_MM;


}





void MMEstimation::WriteMMFileFastForPullTest(string fileINPUT, string fileOUTPUT, float epsilon_s, float epsilon_f){

 
    TFile* file_MM_INPUT   = new TFile((fileINPUT).c_str());
    file_MM_INPUT->cd();
     
    TH1F* tmpSignal[3];
    TH1F* tmpW[3];
    TH1F* tmpQCD[3];

    for(unsigned int iso_index=0; iso_index<3; iso_index++){
     tmpSignal[iso_index] = (TH1F*)gDirectory->Get(("MMExpected_"+IsoNames[iso_index]+theChannel+"_TTbarSignal").c_str());
     tmpW[iso_index] = (TH1F*)gDirectory->Get(("MMExpected_"+IsoNames[iso_index]+theChannel+"_TTbarSemileptonic").c_str());
     tmpQCD[iso_index] = (TH1F*)gDirectory->Get(("MMExpected_"+IsoNames[iso_index]+theChannel+"_TTbar").c_str());
    }

  
    tmpSignal[1]->Add(tmpSignal[1], tmpSignal[0], 0, 2*epsilon_s - epsilon_s * epsilon_s);
    tmpW[1]->Add(tmpW[1], tmpW[0], 0, epsilon_s + epsilon_f - epsilon_s * epsilon_f);
    tmpQCD[1]->Add(tmpQCD[1], tmpQCD[0], 0, 2*epsilon_f - epsilon_f * epsilon_f);

    tmpSignal[2]->Add(tmpSignal[2], tmpSignal[0], 0, epsilon_s * epsilon_s);
    tmpW[2]->Add(tmpW[2], tmpW[0], 0, epsilon_s * epsilon_f);
    tmpQCD[2]->Add(tmpQCD[2], tmpQCD[0], 0, epsilon_f * epsilon_f);

    
   
    TFile* file_MM_OUTPUT = new TFile((fileOUTPUT).c_str(),"RECREATE");
    file_MM_OUTPUT->cd();

   
    for(unsigned int iso_index=0; iso_index<3; iso_index++){
     tmpSignal[iso_index]->Write();
     tmpW[iso_index]->Write();
     tmpQCD[iso_index]->Write();
    }
   
    for(unsigned int iso_index=0; iso_index<3; iso_index++){
     tmpSignal[iso_index] = 0;
     tmpW[iso_index] = 0;
     tmpQCD[iso_index] = 0;
     delete tmpSignal[iso_index];
     delete tmpW[iso_index];
     delete tmpQCD[iso_index];
    }
    
    file_MM_INPUT->Close();
    delete file_MM_INPUT;



    for(unsigned int bin_index=0; bin_index<theNBins; bin_index++){
      for(unsigned int iso_index=0; iso_index<3; iso_index++){
	theDistributions[bin_index].NMMEstimatedQCDDistribution[iso_index]->Write();
  	theDistributions[bin_index].NMMEstimatedWJetsDistribution[iso_index]->Write();
 	theDistributions[bin_index].NMMEstimatedSignalDistribution[iso_index]->Write();
      }
    }

    for(unsigned int iso_index=0; iso_index<3; iso_index++){
      theMMEstimatedPlots.MMEstimated_Signal[iso_index]->Write();
      theMMEstimatedPlots.MMEstimated_QCD[iso_index]->Write();
      theMMEstimatedPlots.MMEstimated_WJets[iso_index]->Write();
    }

    // Writing histos contaning all information concerning NSelected
    TH1F* MMNSelected[3];
    for(unsigned int iso_index=0; iso_index<3; iso_index++){
      MMNSelected[iso_index] = new TH1F(("MMNSelected_"+IsoNames[iso_index]+theChannel).c_str(), ("MMNSelected_"+IsoNames[iso_index]+theChannel).c_str(),  theNBins, theMinBin, theMaxBin);
    }
    for(unsigned int iso_index=0; iso_index<3; iso_index++){
     for(unsigned int bin_index = 0; bin_index < theNBins; bin_index++){
       MMNSelected[iso_index]->SetBinContent(bin_index+1, theNSelected[bin_index].NSel[iso_index]);
     }
     MMNSelected[iso_index]->Write();
    }
    for(unsigned int iso_index=0; iso_index<3; iso_index++){
     MMNSelected[iso_index] = 0;
     delete MMNSelected[iso_index];
    }




    file_MM_OUTPUT->Write();
    file_MM_OUTPUT->Close();
    delete file_MM_OUTPUT;


}
