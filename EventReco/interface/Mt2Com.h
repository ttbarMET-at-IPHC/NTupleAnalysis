#ifndef Mt2Com_h
#define Mt2Com_h

#include <memory>
#include <vector>

//MT2 HEADERS (Modifiés)
#include "EventReco/interface/Basic_Mt2_332_CalculatorMod.h"

//ROOT HEADERS
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TMath.h"

class Mt2Com
{
public:

Mt2Com();

~Mt2Com();

// ----------------------------------------------------------------------------
// Global functions
// ----------------------------------------------------------------------------

void SetVisA(const TLorentzVector& TLorVectA);

void SetVisB(const TLorentzVector& TLorVectB);

void SetVisA(const Double_t Apx,
						 const Double_t Apy,
						 const Double_t Apz,
						 const Double_t Am);

void SetVisB(const Double_t Bpx,
						 const Double_t Bpy,
						 const Double_t Bpz,
						 const Double_t Bm);

void SetMInvisMass(const double invis_mass);

TLorentzVector GetVisA() const;
TLorentzVector GetVisB() const;

TVector2 GetpTMiss() const;
double GetMInvisMass() const;
double GetMt2_332() const;

TVector2 transverse(const TLorentzVector& v) const;
TVector2 neg(const TVector2& w) const;


private:

Double_t Apx;
Double_t Apy;
Double_t Apz;
Double_t Ae;
Double_t Am;

Double_t Bpx;
Double_t Bpy;
Double_t Bpz;
Double_t Be;
Double_t Bm;

double invis_mass;
double m_invis_mass;

TLorentzVector TLorVectA;
TLorentzVector TLorVectB;
TLorentzVector VisA;
TLorentzVector VisB;
TLorentzVector m_otherVisible4Mom;
TLorentzVector v;

TVector2 w;
};

#endif
