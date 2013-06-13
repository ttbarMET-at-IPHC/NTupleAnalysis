#include "../interface/StopMCCharacterization.h"

StopMCCharacterization::StopMCCharacterization()
{
}

StopMCCharacterization::~StopMCCharacterization()
{
  delete hPt_Stop;
  delete hPt_Neutralino;
  delete hPt_Top;
  delete hPt_W;
}

void StopMCCharacterization::Initialize(int nbinsx, float xmin, float xmax, int nbinsy, float ymin, float ymax)
{
    hPt_Stop = new TH3F("hPt_Stop","Stop - Pt",nbinsx,xmin,xmax,nbinsy,ymin,ymax,20,0,500);
    hPt_Neutralino = new TH3F("hPt_Neutralino","Neutralino - Pt",nbinsx,xmin,xmax,nbinsy,ymin,ymax,20,0,500);
    hPt_Top = new TH3F("hPt_Top","top - Pt",nbinsx,xmin,xmax,nbinsy,ymin,ymax,20,0,500);
    hPt_W = new TH3F("hPt_W","W -Pt",nbinsx,xmin,xmax,nbinsy,ymin,ymax,20,0,500);
}

void StopMCCharacterization::Fill(StopMCinfo* stopMCinfo)
{
   if(!stopMCinfo) return;
   if(!stopMCinfo->isStop2topChi2Event()) return;

   float mStop = stopMCinfo->GetStopMass();
   float mNeutralino = stopMCinfo->GetNeutralinoMass();

   // fill stop info
   for(unsigned int i=0;i<stopMCinfo->GetStops().size();i++)
   {
   	hPt_Stop->Fill(mStop,mNeutralino,stopMCinfo->GetStops()[i]->p4.Pt());
   }
   
   // fill neutralino info
   for(unsigned int i=0;i<stopMCinfo->GetNeutralinos().size();i++)
   {
   	hPt_Neutralino->Fill(mStop,mNeutralino,stopMCinfo->GetNeutralinos()[i]->p4.Pt());
   }
   
   // fill top info
   for(unsigned int i=0;i<stopMCinfo->GetTops().size();i++)
   {
   	hPt_Top->Fill(mStop,mNeutralino,stopMCinfo->GetTops()[i]->p4.Pt());
   }
   
   // fill W info
   for(unsigned int i=0;i<stopMCinfo->GetWs().size();i++)
   {
   	hPt_W->Fill(mStop,mNeutralino,stopMCinfo->GetWs()[i]->p4.Pt());
   }
}

void StopMCCharacterization::Compute()
{
}

void StopMCCharacterization::Write(TFile* fout)
{
	if(!fout) return;
	fout->cd();
	fout->mkdir("StopMCCharaterization");
	fout->cd("StopMCCharaterization");
	hPt_Stop->Write();
	hPt_Neutralino->Write();
	hPt_Top->Write();
	hPt_W->Write();

}
