// Revision du 13 avril 2012

#include "EffEstimation/interface/SFBweight.h"
#include "Tools/interface/FileExists.h"
#include <iostream>
#include <fstream>
#include <sstream>

SFBweight::SFBweight(){
  btag_algo_=-1;
  btag_discri_=-999.;
  n_bjets_= 0;
  method_origin1_="default";
  method_origin2_="default";

  histo_sfvalb_=0;
  histo_sferrb_=0;
  histo_sfvall_=0;
  histo_sferrl_=0;

  histo_effvalb_=0;
  histo_efferrb_=0;
  histo_effvall_=0;         
  histo_efferrl_=0;

  histo_effmcb_=0;
  histo_effmcc_=0;
  histo_effmcl_=0;
  histo_errmcb_=0;
  histo_errmcc_=0;
  histo_errmcl_=0;

  map_effmcb_.clear();
  map_effmcc_.clear();
  map_effmcl_.clear();
  map_errmcb_.clear();
  map_errmcc_.clear();
  map_errmcl_.clear();

  tt_sf = 0.;
  tt_sf_err = 0.;
}


SFBweight::SFBweight(const SFBweight& w){
  cout<<"Call the copy constructor of SFBweight"<<endl;
  btag_algo_= w.btag_algo_;
  btag_discri_= w.btag_discri_;
  n_bjets_= w.n_bjets_;
  method_origin1_= w.method_origin1_;
  method_origin2_= w.method_origin2_;

  if(histo_sfvalb_!=0) histo_sfvalb_ = (TH2D*)  w.histo_sfvalb_->Clone("");
  if(histo_sferrb_!=0) histo_sferrb_ = (TH2D*)  w.histo_sferrb_->Clone("");
  if(histo_sfvall_!=0) histo_sfvall_ = (TH2D*)  w.histo_sfvall_->Clone("");
  if(histo_sferrl_!=0) histo_sferrl_ = (TH2D*)  w.histo_sferrl_->Clone("");

  if(histo_effvalb_!=0) histo_effvalb_ = (TH2D*)  w.histo_effvalb_->Clone("");
  if(histo_efferrb_!=0) histo_efferrb_ = (TH2D*)  w.histo_efferrb_->Clone("");
  if(histo_effvall_!=0) histo_effvall_ = (TH2D*)  w.histo_effvall_->Clone("");         
  if(histo_efferrl_!=0) histo_efferrl_ = (TH2D*)  w.histo_efferrl_->Clone("");

  if(histo_effmcb_!=0) histo_effmcb_ = (TH2D*)  w.histo_effmcb_->Clone("");
  if(histo_effmcc_!=0) histo_effmcc_ = (TH2D*)  w.histo_effmcc_->Clone("");
  if(histo_effmcl_!=0) histo_effmcl_ = (TH2D*)  w.histo_effmcl_->Clone("");
  if(histo_errmcb_!=0) histo_errmcb_ = (TH2D*)  w.histo_errmcb_->Clone("");
  if(histo_errmcc_!=0) histo_errmcc_ = (TH2D*)  w.histo_errmcc_->Clone("");
  if(histo_errmcl_!=0) histo_errmcl_ = (TH2D*)  w.histo_errmcl_->Clone("");
 
  //clear maps
  map_effmcc_.clear(); 
  map_effmcb_.clear(); 
  map_effmcl_.clear(); 
  map_errmcc_.clear(); 
  map_errmcb_.clear(); 
  map_errmcl_.clear(); 


  //normaly a copy of the map is needed ....
  map<string,TH2D*>::iterator it;
  
  //for ( it=w.map_effmcc_.begin() ; it != w.map_effmcc_.end(); it++ ) cout<<"toto"<<endl;//map_effmcc_.insert(pair<string,TH2D*>(it->first,(TH2D*) it->second->Clone("")));
  /*
  for ( it=w.map_effmcc_.begin() ; it != w.map_effmcc_.end(); it++ ) map_effmcc_.insert(pair<string,TH2D*>(it->first,(TH2D*) it->second->Clone("")));
  for ( it=w.map_effmcb_.begin() ; it != w.map_effmcb_.end(); it++ ) map_effmcb_.insert(pair<string,TH2D*>(it->first,(TH2D*) it->second->Clone("")));
  for ( it=w.map_effmcl_.begin() ; it != w.map_effmcl_.end(); it++ ) map_effmcl_.insert(pair<string,TH2D*>(it->first,(TH2D*) it->second->Clone("")));
  for ( it=w.map_effmcc_.begin() ; it != w.map_effmcc_.end(); it++ ) map_effmcc_.insert(pair<string,TH2D*>(it->first,(TH2D*) it->second->Clone("")));
  for ( it=w.map_effmcb_.begin() ; it != w.map_effmcb_.end(); it++ ) map_effmcb_.insert(pair<string,TH2D*>(it->first,(TH2D*) it->second->Clone("")));
  for ( it=w.map_effmcl_.begin() ; it != w.map_effmcl_.end(); it++ ) map_effmcl_.insert(pair<string,TH2D*>(it->first,(TH2D*) it->second->Clone("")));
  */
   tt_sf = w.getTTsf();
   tt_sf_err = w.getTTsfErr();
       
}

SFBweight::SFBweight(int btag_algo, float btag_discri, int n_bjets){
         SFBinit(btag_algo,btag_discri,n_bjets);
}


SFBweight::~SFBweight(){
}

void SFBweight::SFBinit(int btag_algo, float btag_discri, int n_bjets){

         //std::cout << " initialisation of SFBweight " << std::endl;
         btag_algo_=btag_algo;
         btag_discri_=btag_discri;
         n_bjets_=n_bjets;
  	//is it usefull if methodb=2 ??
         //std ::cout << " ALGO " <<  btag_algo << " & DISCRI CUT " << btag_discri << std::endl;
         if (btag_algo==0) {
             // TrackCounting
             algoname="TCHE";
             if (btag_discri>1.69 &&  btag_discri<1.71) {
                    method_origin1_="BTAG"+algoname+"L";
                    method_origin2_="MISTAG"+algoname+"L";
		    wp = "L";
             }
             else if (btag_discri>3.29 && btag_discri<3.31) {
                    method_origin1_="BTAG"+algoname+"M";
                    method_origin2_="MISTAG"+algoname+"M";
		    wp = "M";
             }
             else if (btag_discri>10.19 && btag_discri<10.21) {
                    //cout << " version hackee par Caro car TCHET pas dans DB " << endl;
                    method_origin1_="BTAG"+algoname+"M";
                    method_origin2_="MISTAG"+algoname+"M";
             }
             else {
                    method_origin1_="default";
                    method_origin2_="default";
             }
         }
         else if (btag_algo==1) {
              // SecondaryVertex
              algoname="SSVHE";
              if (btag_discri>1.73 && btag_discri<1.75) {
                    method_origin1_="BTAG"+algoname+"M";
                    method_origin2_="MISTAG"+algoname+"M";
              }
              else if (btag_discri>3.04 && btag_discri<3.06) {
                    //cout << " version hackee par Caro car SSVHET pas dans DB " << endl;
                    method_origin1_="BTAG"+algoname+"M";
                    method_origin2_="MISTAG"+algoname+"M";
		    wp = "M";
              }
              else {
                    method_origin1_="default";
                    method_origin2_="default";
              }
         }
         else if (btag_algo==2) {
             // TrackCounting
             algoname="TCHP";
             if (btag_discri>1.92 &&  btag_discri<1.94) {
                    method_origin1_="BTAG"+algoname+"M";
                    method_origin2_="MISTAG"+algoname+"M";
		    wp = "M";
             }
             else if (btag_discri>3.40 && btag_discri<3.42) {
                    method_origin1_="BTAG"+algoname+"T";
                    method_origin2_="MISTAG"+algoname+"T";
		    wp = "T";
              }
             else if (btag_discri>1.18 && btag_discri<1.20) {
                    //cout << " version hackee par Caro car TCHPL pas dans DB " << endl;
                    method_origin1_="BTAG"+algoname+"M";
                    method_origin2_="MISTAG"+algoname+"M";
              }
              else {
                    method_origin1_="default";
                    method_origin2_="default";
              }
         }
         else if (btag_algo==3) {
             // SecondaryVertex
              algoname="SSVHP";
             if (btag_discri>1.99 && btag_discri<2.01) {
                    method_origin1_="BTAG"+algoname+"T";
                    method_origin2_="MISTAG"+algoname+"T";
		    wp = "T";
              }
              else {
                    method_origin1_="default";
                    method_origin2_="default";
              }
         }
         else if (btag_algo==4) {
             // JP
              algoname="JP";
             if (btag_discri>0.274 && btag_discri<0.276) {
                    method_origin1_="BTAG"+algoname+"L";
                    method_origin2_="MISTAG"+algoname+"L";
		    wp = "L";
             }
             else if (btag_discri>0.544 && btag_discri<0.546) {
                    method_origin1_="BTAG"+algoname+"M";
                    method_origin2_="MISTAG"+algoname+"M";
		    wp = "M";
             }      
             else if (btag_discri>0.78 && btag_discri<0.8) {
                    method_origin1_="BTAG"+algoname+"T";
                    method_origin2_="MISTAG"+algoname+"T";
		    wp = "T";
              }     
              else {
                    method_origin1_="default";
                    method_origin2_="default";
              }
         }   
         else if (btag_algo==5) { 
             // JBP 
              algoname="JBP";
             if (btag_discri>1.32 && btag_discri<1.34) {
                    method_origin1_="BTAG"+algoname+"L";
                    method_origin2_="MISTAG"+algoname+"L";
		    wp = "L";
             }      
             else if (btag_discri>2.54 && btag_discri<2.56) {
                    method_origin1_="BTAG"+algoname+"M";
                    method_origin2_="MISTAG"+algoname+"M";
		    wp = "M";
             }      
             else if (btag_discri>3.73 && btag_discri<3.75) {
                    method_origin1_="BTAG"+algoname+"T";
                    method_origin2_="MISTAG"+algoname+"T";
		    wp = "T";
              }
              else {
                    method_origin1_="default";
                    method_origin2_="default";
              }     
         }    
         else if (btag_algo==6) {
             // CSV 
              algoname="CSV";
             if (btag_discri>0.243 && btag_discri<0.245) {
                    method_origin1_="BTAG"+algoname+"L";
                    method_origin2_="MISTAG"+algoname+"L";
		    wp = "L";
             }
             else if (btag_discri>0.678 && btag_discri<0.68) {
                    method_origin1_="BTAG"+algoname+"M";
                    method_origin2_="MISTAG"+algoname+"M";
		    wp = "M";
             }
             else if (btag_discri>0.897 && btag_discri<0.899) {
                    method_origin1_="BTAG"+algoname+"T";
                    method_origin2_="MISTAG"+algoname+"T";
		    wp = "T";
              }
              else {
                    method_origin1_="default";
                    method_origin2_="default";
              }
         }
         else {
             std::cout << " WHICH B-TAG ALGO DO YOU WANT??? " << std::endl;
             method_origin1_="default";
             method_origin2_="default";
         }
         //std::cout << " --> ALGO : " << algoname << ", METHODS for SF&Eff : " << method_origin1_ << " and " << method_origin2_ << std::endl;


}

void SFBweight::LoadInfo(string fileName){


  string fullFileName(getenv( "CMSSW_BASE" )+string("/src/NTuple/NTupleAnalysis/macros/data/")+fileName);
  std::cout<<"Reading the b-tag info from file "<<fullFileName<<endl;
  fexists(fullFileName, true);
  TFile* f = new TFile(fullFileName.c_str());

        std::string histo_name1= "h_" + method_origin1_ + "_BTAG" + "B" +"EFFCORR";
        std::string histo_name2= "h_" + method_origin1_ + "_BTAG" + "B" +"ERRCORR";
//        std::string histo_name3= "h_" + method_origin1_ + "_BTAG" + "B" +"EFF";  // no Eff in DB for b-quark
//        std::string histo_name4= "h_" + method_origin1_ + "_BTAG" + "B" +"ERR";  // no Err(Eff) in DB for b-quark
        std::string histo_name3= "h_" + method_origin1_ + "_BTAG" + "B" +"EFFCORR";
        std::string histo_name4= "h_" + method_origin1_ + "_BTAG" + "B" +"ERRCORR";
        std::string histo_name5= "h_" + method_origin2_ + "_BTAG" + "L" +"EFFCORR";
        std::string histo_name6= "h_" + method_origin2_ + "_BTAG" + "L" +"ERRCORR";
        std::string histo_name7= "h_" + method_origin2_ + "_BTAG" + "L" +"EFF";
        std::string histo_name8= "h_" + method_origin2_ + "_BTAG" + "L" +"ERR";
 
        if (method_origin1_ !="default" && method_origin2_!="default") {
         histo_sfvalb_    = (TH2D*) gDirectory->Get( histo_name1.c_str() )->Clone("") ;
         histo_sferrb_    = (TH2D*) gDirectory->Get( histo_name2.c_str() )->Clone("") ;
         histo_effvalb_   = (TH2D*) gDirectory->Get( histo_name3.c_str() )->Clone("") ;
         histo_efferrb_   = (TH2D*) gDirectory->Get( histo_name4.c_str() )->Clone("") ;
         histo_sfvall_    = (TH2D*) gDirectory->Get( histo_name5.c_str() )->Clone("") ;
         histo_sferrl_    = (TH2D*) gDirectory->Get( histo_name6.c_str() )->Clone("") ;
         histo_effvall_   = (TH2D*) gDirectory->Get( histo_name7.c_str() )->Clone("") ;
         histo_efferrl_   = (TH2D*) gDirectory->Get( histo_name8.c_str() )->Clone("") ;



         histo_sfvalb_->SetDirectory(0); // 
         histo_sferrb_->SetDirectory(0); // 
         histo_effvalb_->SetDirectory(0); // 
         histo_efferrb_->SetDirectory(0); // 
         histo_sfvall_->SetDirectory(0); // 
         histo_sferrl_->SetDirectory(0); // 
         histo_effvall_->SetDirectory(0); // 
         histo_efferrl_->SetDirectory(0); // 

         std::cout << " loaded histo: " <<  histo_name1.c_str() << std::endl;
         std::cout << " loaded histo: " <<  histo_name2.c_str() << std::endl;
         std::cout << " loaded histo: " <<  histo_name3.c_str() << std::endl;
         std::cout << " loaded histo: " <<  histo_name4.c_str() << std::endl;
         std::cout << " loaded histo: " <<  histo_name5.c_str() << std::endl;
         std::cout << " loaded histo: " <<  histo_name6.c_str() << std::endl;
         std::cout << " loaded histo: " <<  histo_name7.c_str() << std::endl;
         std::cout << " loaded histo: " <<  histo_name8.c_str() << std::endl;
        }
        else {
         std::cout << " no histogram loaded : not the correct methods! " << std::endl;
        }

        f->Close();
}




void SFBweight::LoadInfo2(string fileName){

  bool debug = true;

  if(debug) cout<<"SFBweight::LoadInfo2"<<endl;
  
  //////////////////////////////////////
  // Reinitialize the maps
  //////////////////////////////////////
  map_effmcb_.clear();
  map_effmcc_.clear();
  map_effmcl_.clear();
  map_errmcb_.clear();
  map_errmcc_.clear();
  map_errmcl_.clear();

	string fullFileName(getenv( "CMSSW_BASE" )+string("/src/NTuple/NTupleAnalysis/macros/data/")+fileName);
	if(debug) cout<<"SFBweight::LoadInfo2:: Loading the file "<<fullFileName<<endl;
	fexists(fullFileName, true);
        TFile* f2 = TFile::Open(fullFileName.c_str());

	TDirectoryFile* newdir = 0;
	for(int i=0;i<f2->GetListOfKeys()->GetSize();i++){
		if(f2->Get(f2->GetListOfKeys()->At(i)->GetName())->ClassName()==string("TDirectoryFile")){
			newdir = (TDirectoryFile*) f2->Get(f2->GetListOfKeys()->At(i)->GetName());
			string keyname = string(newdir->GetName());
			if(debug) cout<<"Loading histograms from the directory "<<keyname<<endl;
			if(newdir!=0){
        			map_effmcb_.insert(pair<string,TH2D*>(keyname,(TH2D*) newdir->Get( "h_eff_bq" )->Clone("") ));
        			map_effmcc_.insert(pair<string,TH2D*>(keyname,(TH2D*) newdir->Get( "h_eff_cq" )->Clone("") ));
        			map_effmcl_.insert(pair<string,TH2D*>(keyname,(TH2D*) newdir->Get( "h_eff_lq" )->Clone("") ));
       				map_errmcb_.insert(pair<string,TH2D*>(keyname,(TH2D*) newdir->Get( "h_err_bq" )->Clone("") ));
       				map_errmcc_.insert(pair<string,TH2D*>(keyname,(TH2D*) newdir->Get( "h_err_cq" )->Clone("") ));
    			   	map_errmcl_.insert(pair<string,TH2D*>(keyname,(TH2D*) newdir->Get( "h_err_lq" )->Clone("") ));
			}
		}
	}

	//SetDirectory(0);
	//This is needed before closing the file

	map<string,TH2D*>::iterator it;
	for ( it=map_effmcb_.begin() ; it != map_effmcb_.end(); it++ ) it->second->SetDirectory(0);
	for ( it=map_effmcc_.begin() ; it != map_effmcc_.end(); it++ ) it->second->SetDirectory(0);
	for ( it=map_effmcl_.begin() ; it != map_effmcl_.end(); it++ ) it->second->SetDirectory(0);
	for ( it=map_errmcb_.begin() ; it != map_errmcb_.end(); it++ ) it->second->SetDirectory(0);
	for ( it=map_errmcc_.begin() ; it != map_errmcc_.end(); it++ ) it->second->SetDirectory(0);
	for ( it=map_errmcl_.begin() ; it != map_errmcl_.end(); it++ ) it->second->SetDirectory(0);
        
	f2->Close();
	delete f2;
}

void SFBweight::LoadTTsf(string fileName){
  string fullFileName(getenv( "CMSSW_BASE" )+string("/src/NTuple/NTupleAnalysis/macros/data/")+fileName);
  cout<<"SFBweight::LoadTTsf:: Loading the file "<<fullFileName<<endl;
  fexists(fullFileName, true);
  std::ifstream file(fullFileName.c_str());
  string algo;
  float sf, sf_err;
  cout <<algoname+wp<<endl;
  
  while (!file.eof()) {
    file >> algo >> sf >> sf_err;
    // cout<<algo << sf << sf_err<<endl;
    if ((algoname+wp)==algo) {
      tt_sf = sf;
      tt_sf_err = sf_err;
      break;
    }
  }

  
}

void SFBweight::InitAlgoAndWP(int algo, float wp){

	bool debug = false;

	if(debug) cout<<"SFBweight::InitAlgoAndWP"<<endl;
	
  	stringstream ss (stringstream::in | stringstream::out);
  	string str_algo;
  	string str_discri;
  	ss << btag_discri_;
	ss << " ";
  	ss << btag_algo_;
  	ss >> str_discri;
  	ss >> str_algo;
	//solve the problem of 2 which should be 2.0
  	if(str_discri.find(".")>str_discri.size()) str_discri+=".0";

	if(debug) cout<<"Algo: "<<str_algo<<" DiscriCut: "<<str_discri<<endl;
  	string keyname = string("algo_")+str_algo+string("_discri_")+str_discri;
	if(debug) cout<<"Try to access: "<<keyname;
        if (map_effmcb_.find(keyname) == map_effmcb_.end()) {
		if(debug) cout<<": not found "<<endl;
		keyname +="0";
		if(debug) cout<<"Try to access: "<<keyname;
        	if (map_effmcb_.find(keyname) == map_effmcb_.end()) {
			if(debug) cout<<": not found "<<endl;
			keyname +="0";
			if(debug) cout<<"Try to access: "<<keyname;
        		if (map_effmcb_.find(keyname) == map_effmcb_.end()){ 
				if(debug) cout<<": not found"<<endl;
			}
			else if(debug) cout<<" found !"<<endl;
		}
		else if(debug) cout<<" found !"<<endl;
	} 
	else if(debug) cout<<" found !"<<endl;
        if (map_effmcb_.find(keyname) == map_effmcb_.end()) {
          cerr<<"SFBweight::LoadInfo2:: No btag efficiency computed for the given b-tagging algo & W.P. "<< keyname <<" --> use TCHEL info by default"<<endl;
          keyname = map_effmcb_.begin()->first;
        }
        if(debug) cout << "We will access to " << keyname << endl;

	//////////////////////////////////////////////////////////////////////////////////
	//Change the pointers: now point to the indicated algo and working point 
	//////////////////////////////////////////////////////////////////////////////////

         if(map_effmcb_.find(keyname)!=map_effmcb_.end()) histo_effmcb_ = map_effmcb_.find(keyname)->second;
	 else cerr<<"SFweight::InitAlgoAndWP:: histo effmcb not loaded for "<<keyname<<endl;
         if(map_effmcc_.find(keyname)!=map_effmcb_.end()) histo_effmcc_ = map_effmcc_.find(keyname)->second;
	 else cerr<<"SFweight::InitAlgoAndWP:: histo effmcc not loaded for "<<keyname<<endl;
         if(map_effmcl_.find(keyname)!=map_effmcb_.end()) histo_effmcl_ = map_effmcl_.find(keyname)->second;
	 else cerr<<"SFweight::InitAlgoAndWP:: histo effmcl not loaded for "<<keyname<<endl;

         if(map_errmcb_.find(keyname)!=map_effmcb_.end()) histo_errmcb_ = map_errmcb_.find(keyname)->second;
	 else cerr<<"SFweight::InitAlgoAndWP:: histo errmcb not loaded for "<<keyname<<endl;
         if(map_errmcc_.find(keyname)!=map_effmcb_.end()) histo_errmcc_ = map_errmcc_.find(keyname)->second;
	 else cerr<<"SFweight::InitAlgoAndWP:: histo errmcc not loaded for "<<keyname<<endl;
         if(map_errmcl_.find(keyname)!=map_effmcb_.end()) histo_errmcl_ = map_errmcl_.find(keyname)->second;
	 else cerr<<"SFweight::InitAlgoAndWP:: histo errmcl not loaded for "<<keyname<<endl;
}

  /**
   * What you get depends on the info variable (but you never get a weight...) 
   *   0: SF from jet-enriched muon events (same for b and c jets)
   *   1: error on SF from jet-enriched muon events (for c jets, you get the b jet error, so use 2x err)
   *   2: data efficiency for b or light jets from jet-enriched muon events (there is no c-jet eff provided)
   *   3: error on data efficiency for b or light jets from jet-enriched muon events (there is no c-jet eff provided)
   *   4: MC efficiency from tt events
   *   5: error on MC efficiency from tt events
   *   6: SF from tt events (only for b & c jets)
   *   7: error on SF from tt events (only for b & c jets)
   */

float SFBweight::GetWeight(int info,  int quarkorigin, float pt, float eta) const{


//        std::cout << " Extraction of weight for info = " << info << ", pt = " << pt << ", eta = " << eta 
//                  << " (" << whichquark << "-jet) " << std::endl;

        if (info<4 && method_origin1_=="default" && method_origin2_=="default") {
            std::cout << " NOT THE CORRECT METHODS ==> return -1!" << std::endl;
            return -1.;
        }

        TH2D* histo_local = 0;

        if (quarkorigin==5 || quarkorigin==4) {
          if (info==0){ histo_local= histo_sfvalb_; }
          else if (info==1) { histo_local= histo_sferrb_; }
          else if (info==2)  {
               if (quarkorigin==5) histo_local= histo_effvalb_; 
               else  { 
	         std::cout << "Data efficiency for c jets from jet-enriched muon events does not exist"<< std::endl;
                 return -1.;
               }
          }
          else if (info==3) {
               if (quarkorigin==5) histo_local= histo_efferrb_; 
               else  { 
	         std::cout << "Data efficiency for c jets from jet-enriched muon events does not exist"<< std::endl;
                 return -1.;
               }
	  }
          else if (info==4) {
               if (quarkorigin==5) histo_local= histo_effmcb_; 
               else  histo_local= histo_effmcc_; 
          }
          else if (info==5) {
               if (quarkorigin==5) histo_local= histo_errmcb_; 
               else  histo_local= histo_errmcc_; 
          }
	  else if (info==6) {
               return tt_sf;
          }
	  else if (info==7) {
               return tt_sf_err;
          }
          else { std::cout << "This option for info doesn't exist ! Choose another info value : [0,7]"<< std::endl;
                return -1.;
          }
        }
        else {
          // case : quarkorigin==0
          if (info==0){ histo_local= histo_sfvall_; }   
          else if (info==1) { histo_local=  histo_sferrl_; }
          else if (info==2) { histo_local=  histo_effvall_; }
          else if (info==3) { histo_local=  histo_efferrl_; }
          else if (info==4) { histo_local= histo_effmcl_; }
          else if (info==5) { histo_local= histo_errmcl_; }
          else { std::cout << "This option for info doesn't exist ! Choose another info value : [0,5]"<< std::endl;
                return -1.;
          }
        }

        float aa = -1.;

        if (quarkorigin==5 && pt>=240) pt=239;  // because limit in DB
	if (quarkorigin==0 && pt>=520) pt=519;  // because limit in DB
        aa= histo_local->GetBinContent( histo_local->FindBin(pt,eta) );

//         histo_local= new TH2D("aa","",10,0,1,10,0,1);
//         histo_local->Delete();

        return aa;
}
vector<float> SFBweight::GetWeigth4BSel(int method_b,  int syst_b, const std::vector<NTJet> &  selJets) const{
      return GetWeigth4BSel(method_b,  syst_b, selJets, -1., -1.);
}

vector<float> SFBweight::GetWeigth4BSel(int method_b,  int syst_b, const std::vector<NTJet> &  selJets, float sf_val_for_b, float sf_val_for_l) const{
       std::vector<float> proba_jetb;
       vector<float> weight;
       weight.push_back(1); // weightb of the event
       weight.push_back(0); // probability to get 0 jet
       weight.push_back(0); // probability to get 1 jet
       weight.push_back(0); // probability to get 2 jet
       weight.push_back(0); // probability to get at least 3 jet
       for(unsigned int j=0;j<selJets.size();j++){
             // only quark from a defined quarkorigin_
             int sectectedflavour=0;
             int quarkorigin=-1;
	     flavour(selJets[j], sectectedflavour, quarkorigin);

             float discri_val=getBtagDiscr(selJets[j], btag_algo_);
 
              float pt_val_jet = selJets[j].p4.Pt();
              if (pt_val_jet>1000.) pt_val_jet=997.;
              float eta_val_jet = selJets[j].p4.Eta();
              if (eta_val_jet<0) eta_val_jet*=-1.;
              if (eta_val_jet>=2.4) eta_val_jet=2.399;

              float pt_val_jet_orig = selJets[j].p4.Pt();
              float eta_val_jet_orig = selJets[j].p4.Eta();

              int info_b=method_b;
              if (method_b>=2) info_b=0;

              // Compute the weight, depending of method_b 
              // method_b = 0 : SF
              // method_b = 1 : Eff_Data
              // method_b = 2 : SF*Eff_MC
              // method_b = 3 : SF*Eff_MC with error depending on SF and Eff

              float weight_jet=0;
              //impose value for SF from outside
              if      (info_b==0 && quarkorigin==5 && sf_val_for_b>-1.) weight_jet=sf_val_for_b;  //b
              else if (info_b==0 && quarkorigin==4 && sf_val_for_b>-1.) weight_jet=sf_val_for_b;  //c
              else if (info_b==0 && quarkorigin==0 && sf_val_for_l>-1.) weight_jet=sf_val_for_l;  //udsg
              // take value from root files
              else weight_jet=GetWeight(info_b*2,quarkorigin,pt_val_jet,eta_val_jet);

              float ptvaleff=pt_val_jet;
              if (pt_val_jet>200.) ptvaleff=199.;
              float eff_mc=GetWeight(4,quarkorigin,ptvaleff,eta_val_jet);
              float sf1=weight_jet;
              float ersf1=0;
              float eref1=0;
              float errw=0;
              if (method_b==3 && syst_b>0) {
                          ersf1=GetWeight(1,quarkorigin,pt_val_jet,eta_val_jet);
                          eref1=GetWeight(5,quarkorigin,ptvaleff,eta_val_jet);
                          errw=sqrt(ersf1*ersf1*eff_mc*eff_mc + sf1*sf1*eref1*eref1)/eff_mc;  // from SF and Eff
              }
              if (method_b==4 && syst_b>0) {
                          eref1=GetWeight(5,quarkorigin,ptvaleff,eta_val_jet);  // from Eff only
              }


              if (syst_b==1) { //+ all
                    if (method_b==3) {
                          weight_jet+=errw;
                    }
                    else if (method_b==4) {
                          eff_mc+=eref1;
                    }
                    else {
                       weight_jet+=GetWeight(info_b*2+1,quarkorigin,pt_val_jet,eta_val_jet);
                    }
              }
              else if (syst_b==2) { //- all
                    if (method_b==3) {
                          weight_jet-=errw;
                    }
                    else if (method_b==4) {
                          eff_mc-=eref1;
                    }
                    else {
                       weight_jet-=GetWeight(info_b*2+1,quarkorigin,pt_val_jet,eta_val_jet);
                    }
              }
              else if (syst_b==3) { //+ only b!
                  if (quarkorigin==4|| quarkorigin==5) {
                    if (method_b==3) {
                          weight_jet+=errw;
                    }
                    else if (method_b==4) {
                          eff_mc+=eref1;
                    }
                    else {
                         weight_jet+=GetWeight(info_b*2+1,quarkorigin,pt_val_jet,eta_val_jet);
                    }
                  }
              }
              else if (syst_b==4) { //- only b!
                  if (quarkorigin==4|| quarkorigin==5) {
                    if (method_b==3) {
                          weight_jet-=errw;
                    }
                    else if (method_b==4) {
                          eff_mc-=eref1;
                    }
                    else {
                      weight_jet-=GetWeight(info_b*2+1,quarkorigin,pt_val_jet,eta_val_jet);
                    }
                  }
              }
              else if (syst_b==5) { //+ only light!
                  if (quarkorigin==0) {
                    if (method_b==3) {
                          weight_jet+=errw;
                    }
                    else if (method_b==4) {
                          eff_mc+=eref1;
                    }
                    else {
                      weight_jet+=GetWeight(info_b*2+1,quarkorigin,pt_val_jet,eta_val_jet);
                    }
                  }
              }
              else if (syst_b==6) { //- only light!
                  if (quarkorigin==0) {
                    if (method_b==3) {
                          weight_jet-=errw;
                    }
                    else if (method_b==4) {
                          eff_mc-=eref1;
                    }
                    else {
                       weight_jet-=GetWeight(info_b*2+1,quarkorigin,pt_val_jet,eta_val_jet);
                    }
                  }
              }
              else if (syst_b==7) { //+ b and - light!
                  if (quarkorigin==4|| quarkorigin==5) {
                    if (method_b==3) {
                          weight_jet+=errw;
                    }
                    else if (method_b==4) {
                          eff_mc+=eref1;
                    }
                    else {
                          weight_jet+=GetWeight(info_b*2+1,quarkorigin,pt_val_jet,eta_val_jet);
                    }
                  }
                  else if (quarkorigin==0) {
                    if (method_b==3) {
                          weight_jet-=errw;
                    }
                    else if (method_b==4) {
                          eff_mc-=eref1;
                    }
                    else {
                         weight_jet-=GetWeight(info_b*2+1,quarkorigin,pt_val_jet,eta_val_jet);
                    }
                  }
              }
              else if (syst_b==8) { //- b and + light!
                  if (quarkorigin==4|| quarkorigin==5) {
                    if (method_b==3) {
                          weight_jet-=errw;
                    }
                    else if (method_b==4) {
                          eff_mc-=eref1;
                    }
                    else {
                         weight_jet-=GetWeight(info_b*2+1,quarkorigin,pt_val_jet,eta_val_jet);
                    }
                  }
                  else if (quarkorigin==0) {
                    if (method_b==3) {
                          weight_jet+=errw;
                    }    
                    else if (method_b==4) {
                          eff_mc+=eref1;
                    }
                    else {
                          weight_jet+=GetWeight(info_b*2+1,quarkorigin,pt_val_jet,eta_val_jet);
                    }
                  }
              }



              if (method_b==0)  {
                // APPLY THE SCALE FACTOR
                if (discri_val>=btag_discri_) {
                   weight[0]*=weight_jet;
                }
                else {
                   //SF (notag) = (1-Eff_data)/(1-Eff_MC) = (1- SF * Eff_MC)/(1-Eff_MC)
                   if (eff_mc<1) weight[0]*=(1 - weight_jet*eff_mc)/(1-eff_mc);
                   else weight[0]*=1.;
                }
              }
              else {
                 // APPLY THE EFFICIENCY or SF*EFF
                 // ATTENTION : in that case, the Selection should stop at Step 6
                 if (method_b>=2) {
                    weight_jet*=eff_mc;
                   if (weight_jet>1.) weight_jet=1.;
                   if (weight_jet<0.) weight_jet=0.;
                 }
                 proba_jetb.push_back(weight_jet);
                 if (weight_jet<0) std::cout << " weight_jet " << weight_jet << " quark "<< quarkorigin 
                                   <<" eta " << eta_val_jet_orig << " pt " << pt_val_jet_orig
                                   << " eff_mc " << eff_mc << std::endl;
                 if (weight_jet>1) std::cout << " weight_jet " << weight_jet << " eta " << eta_val_jet << " pt " << pt_val_jet << std::endl;

              }



       }
       //         btag_nbjet;

       if (method_b>=1) {

           // combinatory
           float proba_0jet = 1.;
           float proba_1jet = 0.;
           float proba_2jet = 0.;
           float proba_atleast3jet = 0.;
           for (unsigned int ijet=0; ijet<proba_jetb.size(); ijet++) {
               proba_0jet *= (1-proba_jetb[ijet]);
               float pp1 = 1;
               for (unsigned int kjet=0; kjet<proba_jetb.size(); kjet++) {
                  if (kjet!=ijet) {
                      pp1 *= (1-proba_jetb[kjet]);
                   }

                   if (kjet>ijet) {
                       float pp2 = 1;
                       for (unsigned int ljet=0; ljet<proba_jetb.size(); ljet++) {
                          if (ljet!=ijet && ljet!=kjet) {
                             pp2 *= (1-proba_jetb[ljet]);
                           }
                        }
                        proba_2jet+=proba_jetb[ijet]*proba_jetb[kjet]*pp2;
                   }
               }
               proba_1jet+=proba_jetb[ijet]*pp1;
            }
            if(proba_jetb.size()>2) proba_atleast3jet= 1- proba_0jet - proba_1jet - proba_2jet;
            if (n_bjets_==0) {
              weight[0]*=1.;
            }
            else if (n_bjets_==1) {
              weight[0]*=1 - proba_0jet;
            }
            else if (n_bjets_==2) {
              weight[0]*=1 - proba_0jet - proba_1jet;
            }
            if (weight[0]<0)  {
                 // not good! :(
                 if (n_bjets_>0 && (proba_0jet<0 || proba_0jet>1 )) std::cout << " P0 " << proba_0jet << std::endl;
                 if (n_bjets_>1 && (proba_1jet<0 || proba_1jet>1 )) {
                       std::cout << " P1 " << proba_1jet << std::endl;
                 }
            }
            weight[1]=proba_0jet;
            weight[2]=proba_1jet;
            weight[3]=proba_2jet;
            weight[4]=proba_atleast3jet;

            if (weight[3]>1.) {
                 // not good! :(
                 std::cout << " weight[3]>1!!!??? " <<std::endl;
            }

       }


       return weight;


}

void SFBweight::flavour(const NTJet & jet, int & sectectedflavour, int & quarkorigin) const
{
  sectectedflavour=0;
  quarkorigin=-1;
  // LIGHT QUARKS 
  if (abs(jet.partonFlavour)==1) sectectedflavour=1; //d quark
  if (abs(jet.partonFlavour)==2) sectectedflavour=1; //u quark
  if (abs(jet.partonFlavour)==3) sectectedflavour=1; //s quark
  if (abs(jet.partonFlavour)==21) sectectedflavour=1; //gluon 
  if (abs(jet.partonFlavour)==21) sectectedflavour=1; //gluon 
  if (sectectedflavour==1) quarkorigin=0;
  // B QUARKS
  if (abs(jet.partonFlavour)==5) {
       sectectedflavour=1;
       quarkorigin=5;
  }
  // C QUARKS 
  else if (abs(jet.partonFlavour)==4) {
       sectectedflavour=1;
       quarkorigin=4;

  }
  // NOT RECOGNIZED --> LIGHT QUARKS, for the moment
  else if (jet.partonFlavour ==0) {
       sectectedflavour=1;
       quarkorigin=0;
  }

  if (sectectedflavour==0) {
     //DEBUG
     std::cout << " partonFlavour " << jet.partonFlavour << std::endl;
  }
}


double SFBweight::getBtagDiscr(const NTJet & jet, const int& algo) const
{
 switch (algo) {
   case 0:
     return jet.bTag["trackCountingHighEffBJetTags"];
     break;
   case 1:
     return jet.bTag["simpleSecondaryVertexHighEffBJetTags"];
     break;
   case 2:
     return jet.bTag["trackCountingHighPurBJetTags"];
     break;
   case 3:
     return jet.bTag["simpleSecondaryVertexHighPurBJetTags"];
     break;
   case 4:
     return jet.bTag["jetProbabilityBJetTags"];
     break;
   case 5:
     return jet.bTag["jetBProbabilityBJetTags"];
     break;
   case 6:
     return jet.bTag["combinedSecondaryVertexBJetTags"];
     break;
   default:
     cerr << "btagAlgo doesn't exist !" << endl;
     return -10000000;
     break;
   }
}


float SFBweight::efficiency(const NTJet & jet,
	bool applySFlSyst, float factorSFlSyst,
	bool applySFbcSyst, float factorSFbcSyst ) const
{
  int sectectedflavour=0;
  int quarkorigin=-1;
  flavour(jet, sectectedflavour, quarkorigin);

//   float discri_val=getBtagDiscr(jet, btag_algo_);

  float pt_val_jet = min(997., jet.p4.Pt());
  float eta_val_jet = min (2.339, fabs(jet.p4.Eta()));
  if (quarkorigin==5 || quarkorigin==4){
	float sf = tt_sf;
	if (applySFbcSyst) {
	  if (factorSFbcSyst>1.) sf += tt_sf_err*(quarkorigin==4?2:1);
	  else  sf -= tt_sf_err*(quarkorigin==4?2:1);
	}
	//cout << "use SFbc " << sf<<" +- "<<tt_sf_err*(quarkorigin==4?2:1)<<endl;
	return (GetWeight(4, quarkorigin, pt_val_jet, eta_val_jet) * sf);
  } else {
	float sf = GetWeight(0, quarkorigin, pt_val_jet, eta_val_jet);
	if (applySFlSyst) {
	  if (factorSFlSyst>1.) sf += GetWeight(1, quarkorigin, pt_val_jet, eta_val_jet);
	  else  sf -= GetWeight(1, quarkorigin, pt_val_jet, eta_val_jet);
	}
	//cout << "use SFl  " << sf<<" +- "<< GetWeight(1, quarkorigin, pt_val_jet, eta_val_jet) <<endl;
	//cout <<pt_val_jet <<" "<<eta_val_jet<<endl;
	return (GetWeight(4, quarkorigin, pt_val_jet, eta_val_jet) * sf);
  }
}


const TH2D* SFBweight::GetHistoSFB() const{
    return histo_sfvalb_;
}

const TH2D* SFBweight::GetHistoEffMCb() const{
	return histo_effmcb_;
}

const TH2D* SFBweight::GetHistoEffMCc() const{
	return histo_effmcc_;
}

const TH2D* SFBweight::GetHistoEffMCl() const{
	return histo_effmcl_;
}
