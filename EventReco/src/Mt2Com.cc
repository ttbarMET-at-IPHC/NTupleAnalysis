//Fait le 25/02/13 par Thomas DESCHLER.

#include "EventReco/interface/Mt2Com.h"

Mt2Com::Mt2Com()
{
	Apx=0;
	Apy=0;
	Apz=0;
	Ae=0;
	Am=0;
	Bpx=0;
	Bpy=0;
	Bpz=0;
	Be=0;
	Bm=0;
	m_invis_mass=0;
}

Mt2Com::~Mt2Com() {}

// ----------------------------------------------------------------------------
// Global functions
// ----------------------------------------------------------------------------

void Mt2Com::SetVisA(const TLorentzVector& TLorVectA) {VisA=TLorVectA;}
void Mt2Com::SetVisB(const TLorentzVector& TLorVectB) {VisB=TLorVectB;}

void Mt2Com::SetVisA(const Double_t Apx,
										 const Double_t Apy,
										 const Double_t Apz,
								 		 const Double_t Am)
{
	Double_t Ae(TMath::Sqrt(Am*Am + Apx*Apx + Apy*Apy + Apz*Apz));
	VisA.TLorentzVector::SetPxPyPzE(Apx,Apy,Apz,Ae);
}

void Mt2Com::SetVisB(const Double_t Bpx,
										 const Double_t Bpy,
										 const Double_t Bpz,
										 const Double_t Bm)
{
	Double_t Be(TMath::Sqrt(Bm*Bm + Bpx*Bpx + Bpy*Bpy + Bpz*Bpz));
	VisB.TLorentzVector::SetPxPyPzE(Bpx,Bpy,Bpz,Be);
}

void Mt2Com::SetMInvisMass(const double invis_mass) {m_invis_mass=invis_mass;}

TLorentzVector Mt2Com::GetVisA() const {return VisA;}
TLorentzVector Mt2Com::GetVisB() const {return VisB;}

TVector2 Mt2Com::GetpTMiss() const
{
	TLorentzVector m_otherVisible4Mom;
	return Mt2Com::neg(Mt2Com::transverse(VisA+VisB+m_otherVisible4Mom));
}

double Mt2Com::GetMInvisMass() const {return m_invis_mass;}

double Mt2Com::GetMt2_332() const
{
	Mt2::Basic_Mt2_332_CalculatorMod mt2Calculator;
	return mt2Calculator.mt2_332(VisA, VisB, Mt2Com::GetpTMiss(), m_invis_mass);	
}

TVector2 Mt2Com::transverse(const TLorentzVector& v) const
{
   return TVector2(v.TLorentzVector::Px(), v.TLorentzVector::Py());
}

TVector2 Mt2Com::neg(const TVector2& w) const
{
   return TVector2(-w.TVector2::Px(), -w.TVector2::Py());
}
