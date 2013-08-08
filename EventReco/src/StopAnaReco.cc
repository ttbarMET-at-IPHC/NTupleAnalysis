#include "EventReco/interface/StopAnaReco.h"


StopAnaReco::StopAnaReco()
{
	HTSSM=0;
	HTOSM=0;
	overLap=0;
	max=0;
}

StopAnaReco::StopAnaReco(TH1F* HS,TH1F* HB)
{
	HistoS=HS;
	HistoB=HB;
	overLap=0;
	max=0;
	x=0;
	Fsig=0;
	EpsS=0;
	EpsB=0;
}

StopAnaReco::StopAnaReco(const IPHCTree::NTMET& NTMET,
												 const std::vector<IPHCTree::NTJet>& NTJet,
												 const std::vector<IPHCTree::NTElectron>& NTElectron,
												 const std::vector<IPHCTree::NTMuon>& NTMuon)
{
	MET=NTMET;
	Jet=NTJet;
	Electron=NTElectron;
	Muon=NTMuon;
}

StopAnaReco::~StopAnaReco() {}

// ----------------------------------------------------------------------------
// Global functions
// ----------------------------------------------------------------------------

bool StopAnaReco::IsSSM(const IPHCTree::NTMET& MET,const TLorentzVector& obj) const
{
		return (((obj.Phi()) < (MET.p2.Phi_mpi_pi(MET.p2.Phi())+TMath::Pi()/2)) && 
						((obj.Phi()) > (MET.p2.Phi_mpi_pi(MET.p2.Phi())-TMath::Pi()/2)));
}

std::pair<float,float> StopAnaReco::HT() const
{
	float HTSSM(0),HTOSM(0);

	for(unsigned int i=0 ; i < Jet.size() ; i++)
	{
		if (IsSSM(MET,Jet[i].p4)) HTSSM+=Jet[i].p4.Pt();
		else HTOSM+=Jet[i].p4.Pt();
	}
	for(unsigned int i=0 ; i < Electron.size() ; i++)
	{
		if (IsSSM(MET,Electron[i].p4)) HTSSM+=Electron[i].p4.Pt();
		else HTOSM+=Electron[i].p4.Pt();
	}
	for(unsigned int i=0 ; i < Muon.size() ; i++)
	{
		if (IsSSM(MET,Muon[i].p4)) HTSSM+=Muon[i].p4.Pt();
		else HTOSM+=Muon[i].p4.Pt();
	}

	return std::make_pair(HTSSM,HTOSM);
}

float StopAnaReco::GetHTSSM() const {return StopAnaReco::HT().first;}

float StopAnaReco::GetHTOSM() const {return StopAnaReco::HT().second;}

float StopAnaReco::GetHT() const {return StopAnaReco::HT().first+StopAnaReco::HT().second;}

float StopAnaReco::GetMt()
{
	Mt=0;
	Mt+=MET.met();
	for(unsigned int i=0 ; i < Electron.size() ; i++)
	{
		Mt+=Electron[i].p4.Mt();
	}
	for(unsigned int i=0 ; i < Muon.size() ; i++)
	{
		Mt+=Muon[i].p4.Mt();
	}
	return Mt;
}


//Normalisation d'un histogramme
void StopAnaReco::Normalize(TH1F* h1)
{
	if(h1->Integral()!=0) h1->Scale(1/h1->Integral());
	else cout<<"StopAnaReco::Normalize : AVERTISSEMENT : L'integrale de l'histogramme "<<h1->GetName()<<" est nulle."<<endl;
}



//Overlap de deux histogrammes
double StopAnaReco::GetOverlap(TH1F* HS,TH1F* HB)
{
	overLap=0;
	HistoS=HS;
	HistoB=HB;
	StopAnaReco::Normalize(HistoS);
	StopAnaReco::Normalize(HistoB);

	if(HS->GetNbinsX()==HB->GetNbinsX())
		for(int i=0; i<=HB->GetNbinsX(); i++)
			{
				if(HS->GetBinContent(i) < HB->GetBinContent(i))
					overLap += HS->GetBinContent(i);
				else 
					overLap += HB->GetBinContent(i);
			}
	else cout<<"StopAnaReco::GetOverlap : AVERTISSEMENT : Les deux histogrammes n'ont pas le même nombre de bins."<<endl;
	return overLap;
}

std::vector<TGraph*> StopAnaReco::ComputeEff(TH1F* HistoS,TH1F* HistoB)
{
	if(HistoS->FindLastBinAbove(0,1)>HistoB->FindLastBinAbove(0,1)) 
		max=HistoS->FindLastBinAbove(0,1);
	else 
		max=HistoB->FindLastBinAbove(0,1);

	double x[max],Fsig[max],EpsB[max],EpsS[max];

	if(HistoS->GetBinWidth(0)==HistoB->GetBinWidth(0))
		for(int i=0;i<=max;i++)
		{
			x[i]=HistoS->GetBinWidth(0)*i;
			EpsB[i]=HistoB->Integral(x[i],max)/HistoB->Integral(0,max);
			EpsS[i]=HistoS->Integral(x[i],max)/HistoS->Integral(0,max);

			if(EpsB[i]!=0) Fsig[i]=EpsS[i]/TMath::Sqrt(EpsB[i]);
			else Fsig[i]=0;
		}
	else cout<<"StopAnaReco::ComputeEff : AVERTISSEMENT : Les deux histogrammes n'ont pas le même nombre de bins."<<endl;

	std::vector<TGraph*> GraphVect(2);

	//Graphique de EpsS/sqrt(EpsB) en fonction de la coupure
	Cut=new TGraph(max,x,Fsig);
	Cut->SetTitle("#varepsilon_{S}/#sqrt{#varepsilon_{B}} en fonction de la coupure.;Coupure;#varepsilon_{S}/#sqrt{#varepsilon_{B}}");
	GraphVect[0]=Cut;

	//Graphique de EpsB en fonction de EpsS
	Rap=new TGraph(max,EpsS,EpsB);
	Rap->SetTitle("#varepsilon_{B}(#varepsilon_{S});#varepsilon_{S};#varepsilon_{B}");
	GraphVect[1]=Rap;

	return GraphVect;
}

double StopAnaReco::GetMaxCutEff(TH1F* HistoS,TH1F* HistoB)
{
	max=0;
	for(int j=0;j<=StopAnaReco::ComputeEff(HistoS,HistoB)[0]->GetMaxSize();j++)
		if(StopAnaReco::ComputeEff(HistoS,HistoB)[0]->GetY()[j]>max) max=StopAnaReco::ComputeEff(HistoS,HistoB)[0]->GetY()[j];

	return max;
}

//Calcul du Kfactor
double StopAnaReco::GetK(const double& Sigma,
												 const double& Luminosity,
												 const int& NofEvents) const
{return Sigma*Luminosity/NofEvents;}

TH1F* StopAnaReco::GetHistoAdd(std::vector<TH1F*>& HistoVect)
{
	for (unsigned int k=0;k<HistoVect.size();k++)
	{
		if (k!=0) HistoVect[0]->Add(HistoVect[k]);
	}
	return HistoVect[0];
}

void StopAnaReco::NormalizeKfactor(TH1F* h1,const double Kfactor)
{
	if(h1->GetEntries()!=0)
	{
		if(Kfactor!=0) h1->Scale(Kfactor);
		else cout<<"StopAnaReco::NormalizeKfactor : AVERTISSEMENT : Le Kfactor est nul."<<endl;
	}	
	else cout<<"StopAnaReco::NormalizeKfactor : AVERTISSEMENT : L'histogramme "<<h1->GetName()<<" est vide."<<endl;
}

TH1F* StopAnaReco::GetHistoNorm(const std::vector<TH1F*>& HistoVect,
															  const std::vector<Double_t>& KVect)
{
	if(HistoVect.size()==KVect.size())
	{
		//cout<<"GetHistoNorm HistoVect.size()= "<<HistoVect.size()<<endl;
		for(unsigned int k=0;k<HistoVect.size();k++)
		{
			HistoVect[k]->Scale(KVect[k]);
			if(k!=0) HistoVect[0]->Add(HistoVect[k]);
		}
	}
	return HistoVect[0];
}

double StopAnaReco::GetDeltaPhi()
{	
	return min(min(abs(MET.p2.Phi()-Jet[0].p4.Phi()-TMath::Pi()),
										 2*TMath::Pi()-abs(MET.p2.Phi()-Jet[0].p4.Phi()-TMath::Pi())),
						 min(abs(MET.p2.Phi()-Jet[1].p4.Phi()-TMath::Pi()),
										 2*TMath::Pi()-abs(MET.p2.Phi()-Jet[1].p4.Phi()-TMath::Pi())));
	/*float dPhi=9999;
	for (int i=0; i<2;i++)
	{
	float dPhitemp=abs(Jet[i].p4.Phi()-MET.p2.Phi()+TMath::Pi());
	cout<<"dPhitemp= "<<dPhitemp<<endl;
	if (dPhitemp > TMath::Pi()) dPhitemp = 2*TMath::Pi()-dPhitemp;
	if (dPhitemp < dPhi) dPhi=dPhitemp;
	}
	return dPhi;*/
}
