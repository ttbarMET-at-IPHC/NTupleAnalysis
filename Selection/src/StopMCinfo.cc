#include "../interface/StopMCinfo.h"


StopMCinfo::StopMCinfo()
{
  Reset();
}

void StopMCinfo::Reset()
{
    isSUSYEvent_ = false;
    isStop2topChiEvent_ = false;
    isStop2bCharginoEvent_ = false;
    ttbarChannel_ = unknown;
    TMEME_ = 0;
    SCNTWMEME_ = 0 ;
    mStop_ = -9999;
    mNeutralino_ = -9999;
    mChargino_ = -9999;
    event_ = 0;
    stop_.clear();
    charginos_.clear();
    neutralinos_.clear();
    top_.clear();
    W_.clear();
    lepton_.clear();
    neutrino_.clear();
    hadronicTau_.clear();
    TauDecay_.clear();
}
 
bool StopMCinfo::MassSelector(const float& mStop, const float& mNeutralino){
    if((mStop_-mStop)<1E-5 && (mNeutralino_-mNeutralino)<1E-5) return true;
    return false;
}

bool StopMCinfo::MassRangeSelector(const float& mStopMin, const float& mStopMax, const float& mNeutralinoMin, const float& mNeutralinoMax)
{
    if( mStop_>mStopMin && mStop_<mStopMax  && mNeutralino_>mNeutralinoMin && mNeutralino_<mNeutralinoMax) return true;
    return false;
}

void StopMCinfo::LoadEvent(IPHCTree::NTEvent* evt)
{
    Reset();
    event_ = evt;
    // loop other variables
    for (unsigned int ip = 0; ip< event_->mc.genParticles.size(); ip++)
    {
        //loop for stops
        if( abs(event_->mc.genParticles[ip].id) == 1000006 && event_->mc.genParticles[ip].isStatus3)
        {
            stop_.push_back(&event_->mc.genParticles[ip]);
            mStop_ = event_->mc.genParticles[ip].p4.M();
        }
        if( abs(event_->mc.genParticles[ip].id) == 1000022 && event_->mc.genParticles[ip].isStatus3)
        {
            neutralinos_.push_back(&event_->mc.genParticles[ip]);
            mNeutralino_ = event_->mc.genParticles[ip].p4.M();
        }
        if( abs(event_->mc.genParticles[ip].id) == 1000024 && event_->mc.genParticles[ip].isStatus3)
        {
            charginos_.push_back(&event_->mc.genParticles[ip]);
            mChargino_ = event_->mc.genParticles[ip].p4.M();
        }
        if( abs(event_->mc.genParticles[ip].id) == 6 && event_->mc.genParticles[ip].isStatus3)
        {
            top_.push_back(&event_->mc.genParticles[ip]);
        }
        if(abs(event_->mc.genParticles[ip].id) == 24 && event_->mc.genParticles[ip].isStatus3 )
        {
            W_.push_back(&event_->mc.genParticles[ip]);
        }
        //looking for a lepton coming from a W
        if( (abs(event_->mc.genParticles[ip].id)==11 || abs(event_->mc.genParticles[ip].id)==13 || abs(event_->mc.genParticles[ip].id)==15) &&
            (event_->mc.genParticles[ip].motherIndex_ >=0 && event_->mc.genParticles[ip].isStatus3 && abs(event_->mc.genParticles[event_->mc.genParticles[ip].motherIndex_].id) == 24))
        {
            lepton_.push_back(&event_->mc.genParticles[ip]);
        }
        //looking for a lepton coming from a W
        if( (abs(event_->mc.genParticles[ip].id)==12 || abs(event_->mc.genParticles[ip].id)==14 || abs(event_->mc.genParticles[ip].id)==16) &&
            (event_->mc.genParticles[ip].motherIndex_ >=0 && event_->mc.genParticles[ip].isStatus3 && abs(event_->mc.genParticles[event_->mc.genParticles[ip].motherIndex_].id) == 24))
        {
            lepton_.push_back(&event_->mc.genParticles[ip]);
        }
        //looking for a TauDecay coming from a W
        if( (abs(event_->mc.genParticles[ip].id)==11 || abs(event_->mc.genParticles[ip].id)==13 || abs(event_->mc.genParticles[ip].id)==15) &&
            (event_->mc.genParticles[ip].motherIndex_ >=0 && event_->mc.genParticles[ip].isStatus3 && abs(event_->mc.genParticles[event_->mc.genParticles[ip].motherIndex_].id) == 15) &&
            (event_->mc.genParticles[event_->mc.genParticles[ip].motherIndex_].motherIndex_ >=0) &&
            abs(event_->mc.genParticles[event_->mc.genParticles[event_->mc.genParticles[ip].motherIndex_].motherIndex_].id) == 24)
        {
            TauDecay_.push_back(&event_->mc.genParticles[ip]);
        }
        if((abs(event_->mc.genParticles[ip].id) >= 1000000 && abs(event_->mc.genParticles[ip].id) < 1000099 ))
        {
           isSUSYEvent_ = true;
           if((abs(event_->mc.genParticles[ip].id)!=1000006 && abs(event_->mc.genParticles[ip].id)!=1000022 && abs(event_->mc.genParticles[ip].id)!=1000024))
           {
              cerr<<"There are other SUSY particles than stop1, chargino1, neutralino1"<<endl;
              cerr<<"pid: "<<event_->mc.genParticles[ip].id<<endl;
           }
        }
    }
    // End of loop over particles

        // Filling hadronicTau
    hadronicTau_ = lepton_;
    for(unsigned int i=0;i<TauDecay_.size();i++)
    {
        for(unsigned int j=hadronicTau_.size();j>0;j--)
        {
            if(&(event_->mc.genParticles[TauDecay_[i]->motherIndex_]) == hadronicTau_[j])
                hadronicTau_.erase(hadronicTau_.begin()+j);
        }
    }
    // procedure to be validated !!
    
    // Fill TMEME & SCNTWMEME 
    TMEME_ = 0;
    SCNTWMEME_ = 0;
    TMEME_        += 10000*top_.size();
    SCNTWMEME_ += 100000*top_.size();
    for(unsigned int i=0;i<lepton_.size();i++)
    {
        if(abs(lepton_[i]->id) == 11){
            TMEME_ += 1;
            SCNTWMEME_ += 1;
        }
        if(abs(lepton_[i]->id) == 13){
            TMEME_ += 10;
            SCNTWMEME_ += 10;
        }
    }
    for(unsigned int i=0;i<TauDecay_.size();i++)
    {
        if(abs(TauDecay_[i]->id) == 11){
            TMEME_ += 100;
            SCNTWMEME_ += 100;
        }
        if(abs(TauDecay_[i]->id) == 13){
            TMEME_ += 1000;
            SCNTWMEME_ += 1000;
        }
    }
    SCNTWMEME_ += 100000000*(Long_t)stop_.size();
    SCNTWMEME_ += 10000000*(Long_t)charginos_.size();
    SCNTWMEME_ += 1000000*(Long_t)neutralinos_.size();

    //Filling channels
    if(TMEME_ == 21010 || TMEME_ == 22000 || TMEME_ == 20020) ttbarChannel_ = mumu;
    if(TMEME_ == 20101 || TMEME_ == 20200 || TMEME_ == 20002) ttbarChannel_ = ee;
    if(TMEME_ == 21100 || TMEME_ == 20011 || TMEME_ == 21001 || TMEME_ == 20110) ttbarChannel_ = emu;
    if(stop_.size()==2 && lepton_.size()==2 && hadronicTau_.size()==2) ttbarChannel_ = tautau; 
    if((TMEME_ == 21000 || TMEME_ == 20010) && hadronicTau_.size()==1) ttbarChannel_ = mutau;
    if((TMEME_ == 20100 || TMEME_ == 2001 ) && hadronicTau_.size()==1) ttbarChannel_ = etau;
    if(TMEME_ == 21000 || TMEME_ == 20010 ) ttbarChannel_ = muJets;
    if(TMEME_ == 20100 || TMEME_ == 2001 ) ttbarChannel_ = eJets;
    if(stop_.size()==2 && lepton_.size()==1 && hadronicTau_.size()==1) ttbarChannel_ = tauJets; 
    if(top_.size()==2 && lepton_.size()==0) ttbarChannel_ = fullyHad;

    //Filling booleans
    if( stop_.size()==2 && neutralinos_.size()==2) isStop2topChiEvent_ = true;
    if( stop_.size()==2 && charginos_.size()==2) isStop2bCharginoEvent_ = true;
    
    //Looking for hadronic W
    HadronicW_ = W_;
    for(unsigned int i=0;i<lepton_.size();i++)
    {
        for(unsigned int j=0;j<HadronicW_.size();j++){
            if(&(event_->mc.genParticles[lepton_[i]->motherIndex_]) == HadronicW_[j])
            HadronicW_.erase(HadronicW_.begin()+j);
            break;
        }
            
    }

    //Looking for hadronic top
    for(unsigned int i=0;i<HadronicW_.size();i++)
    {
        if(abs(event_->mc.genParticles[HadronicW_[i]->motherIndex_].id) == 6)
            HadronicTop_.push_back(&(event_->mc.genParticles[HadronicW_[i]->motherIndex_]));
    }
}


std::string StopMCinfo::Bool2Text(const bool& value) const{
    if(value) return string(" True ");
    else return string(" False ");
}

std::string StopMCinfo::PrintChannel() const{
    std::string out;
    if(ttbarChannel_ == ee) out=std::string("ee");
    if(ttbarChannel_ == emu) out=std::string("emu");
    if(ttbarChannel_ == mumu) out=std::string("mumu");
    if(ttbarChannel_ == etau) out=std::string("etau");
    if(ttbarChannel_ == tautau) out=std::string("tautau");
    if(ttbarChannel_ == eJets) out=std::string("eJets");
    if(ttbarChannel_ == muJets) out=std::string("muJets");
    if(ttbarChannel_ == tauJets) out=std::string("tauJets");
    if(ttbarChannel_ == fullyHad) out=std::string("fullyHad");
    if(ttbarChannel_ == unknown) out=std::string("unknown");
    return out;
}

void StopMCinfo::Print() const{
    cout<<"#############################"<<endl;
    cout<<"#   StopMCinfo::Print       #"<<endl;
    cout<<"#############################"<<endl;
    cout<<"# isSUSY : "<< Bool2Text(isSUSYEvent_)<<endl;
    cout<<"# stop->top+chi : "<<Bool2Text(isStop2topChiEvent_)<<endl;
    cout<<"# stop->chargino+b : "<<Bool2Text(isStop2bCharginoEvent_)<<endl;
    cout<<"# TMEME: "<<TMEME_<<endl;
    cout<<"# SCNTWMEME: "<<SCNTWMEME_<<endl;
    cout<<"# m(stop) = "<<mStop_<<" m(neutralino) = "<<mNeutralino_<<endl;
    cout<<"# channel = "<<PrintChannel()<<endl;
    cout<<"#############################"<<endl;
}


