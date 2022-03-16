#define SkimerBoost_cxx
#include "SkimerBoost.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <TMath.h>
#include <string>
#include <sstream>
using namespace std;


void SkimerBoost::Loop(TString OutputFile)
{
    TH1F* diEle_OS = new TH1F("diEle_OS","diEle_OS",40, 0, 160);
    TH1F* diEle_SS = new TH1F("diEle_SS","diEle_SS",40, 0, 160);

    TFile* file = TFile::Open(OutputFile, "RECREATE");
    
    fChain->SetBranchStatus("*",1);
    
    if (fChain == 0) return;
    
    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    float MuMass= 0.10565837;
    float eleMass= 0.000511;
    int ssRun [300] = {};
    int ssRunSize = 0;
    int osRun [5000] = {};
    int osRunSize = 0;
    
    
    for (int jentry=0; jentry<nentries;jentry++) {
        
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        
        
        
        if(jentry % 10000 == 0) cout << "Processed " << jentry << " events out of " <<nentries<<endl;
      
        
        
        TLorentzVector LeadEle4Momentum, SubEle4Momentum, ZCandidate;
        
        auto numDiEle(0);
        for (int iele = 0; iele < nEle; ++iele){
	  bool eleMVAId= false;
	  if (fabs (eleSCEta->at(iele)) <= 0.8 && eleIDMVANoIso->at(iele) >    0.837   ) eleMVAId= true;
        else if (fabs (eleSCEta->at(iele)) >  0.8 &&fabs (eleSCEta->at(iele)) <=  1.5 && eleIDMVANoIso->at(iele) >    0.715   ) eleMVAId= true;
        else if ( fabs (eleSCEta->at(iele)) >=  1.5 && eleIDMVANoIso->at(iele) >   0.357   ) eleMVAId= true;
        else eleMVAId= false;
        
        if (!eleMVAId) continue;

        float IsoLep1Value=elePFChIso->at(iele)/elePt->at(iele);
        if ( (elePFNeuIso->at(iele) + elePFPhoIso->at(iele) - 0.5* elePFPUIso->at(iele) )  > 0.0)
            IsoLep1Value= ( elePFChIso->at(iele) + elePFNeuIso->at(iele) + elePFPhoIso->at(iele) - 0.5* elePFPUIso->at(iele))/elePt->at(iele);
            
	LeadEle4Momentum.SetPtEtaPhiM(elePt->at(iele),eleEta->at(iele),elePhi->at(iele),eleMass);
            
            for (int jele = iele+1; jele < nEle; ++jele){
	      bool eleMVAId2= false;
	  if (fabs (eleSCEta->at(jele)) <= 0.8 && eleIDMVANoIso->at(jele) >    0.837   ) eleMVAId2= true;
        else if (fabs (eleSCEta->at(jele)) >  0.8 &&fabs (eleSCEta->at(jele)) <=  1.5 && eleIDMVANoIso->at(jele) >    0.715   ) eleMVAId2= true;
        else if ( fabs (eleSCEta->at(jele)) >=  1.5 && eleIDMVANoIso->at(jele) >   0.357   ) eleMVAId2= true;
        else eleMVAId2= false;
        
        if (!eleMVAId2) continue;
        
        float IsoLep2Value=elePFChIso->at(jele)/elePt->at(jele);
        if ( (elePFNeuIso->at(iele) + elePFPhoIso->at(iele) - 0.5* elePFPUIso->at(iele) )  > 0.0)
            IsoLep2Value= ( elePFChIso->at(jele) + elePFNeuIso->at(jele) + elePFPhoIso->at(jele) - 0.5* elePFPUIso->at(jele))/elePt->at(jele);	       

                
                SubEle4Momentum.SetPtEtaPhiM(elePt->at(jele),eleEta->at(jele),elePhi->at(jele),eleMass);

                numDiEle++;
		if (!(IsoLep1Value && IsoLep2Value && eleMVAId && eleMVAId2)) continue;
		ZCandidate = SubEle4Momentum + LeadEle4Momentum;
		bool OS = eleCharge->at(iele) * eleCharge->at(jele) < 0;
		bool SS = eleCharge->at(iele) * eleCharge->at(jele) > 0;
		if (ZCandidate.M() >= 86 && ZCandidate.M() <= 96) {
		  if (OS) {
		    osRun[osRunSize] = run;
		    osRunSize++;
		  }
		  if (SS) {
		    ssRun[ssRunSize] = run;
		    ssRunSize++;
		  }
		}
	       	if (OS){
		  diEle_OS->Fill(ZCandidate.M());
		  break;
		}
		if (SS){
		  diEle_SS->Fill(ZCandidate.M());
		  break;
		}
		break;    
            }
        }
        
        if(numDiEle < 1) continue;
        
 
    }
    int i;
    TH1::SetDefaultSumw2();
    TH1F* osEventsVsRun = new TH1F("OSrun", "OSrun", 100, 315200, 325200);
    TH1F* ssEventsVsRun = new TH1F("SSrun", "SSrun", 100, 315200, 325200);
    for (i = 0; i < osRunSize; i++) {
      osEventsVsRun->Fill(osRun[i]);
    }
    for (i = 0; i < ssRunSize; i++) {
      ssEventsVsRun->Fill(ssRun[i]);
    }
    osEventsVsRun->Write();
    ssEventsVsRun->Write();
    diEle_OS->Write();
    diEle_SS->Write();
    file->Close();
}

int main(int argc, char* argv[]){
    
    string InputFile=argv[1];
    string OutputFile=argv[2];
    
    cout<< "\n===\n input is "<<InputFile  <<"  and output is "<<OutputFile<<"\n===\n";
    
    SkimerBoost t(InputFile);
    t.Loop(OutputFile);
    
    return 0;
}
