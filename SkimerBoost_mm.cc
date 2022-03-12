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
    TH1F* diMu_OS = new TH1F("diMu_OS","diMu_OS",30, 0, 160);
    TH1F* diMu_SS = new TH1F("diMu_SS","diMu_SS",30, 0, 160);      

    TFile* file = TFile::Open(OutputFile, "RECREATE");
    
    fChain->SetBranchStatus("*",1);
    

    
    if (fChain == 0) return;
    
    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    float MuMass= 0.10565837;
    float eleMass= 0.000511;
    int ssRun [1000] = {};
    int ssRunSize = 0;
    int osRun [5000] = {};
    int osRunSize = 0;
    
    
    
    for (int jentry=0; jentry<nentries;jentry++) {
        
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        
        
        
        if(jentry % 10000 == 0) cout << "Processed " << jentry << " events out of " <<nentries<<endl;
        

        
        
        TLorentzVector LeadMu4Momentum, SubMu4Momentum, ZCandidate;
        
        auto numDiMu(0);
        for (int imu = 0; imu < nMu; ++imu){
	  float IsoLep1Value=muPFChIso->at(imu)/muPt->at(imu);
        if ( (muPFNeuIso->at(imu) + muPFPhoIso->at(imu) - 0.5* muPFPUIso->at(imu) )  > 0.0)
            IsoLep1Value= ( muPFChIso->at(imu) + muPFNeuIso->at(imu) + muPFPhoIso->at(imu) - 0.5* muPFPUIso->at(imu))/muPt->at(imu);
            if (muPt->at(imu) < 30 || muEta->at(imu) > 2.4) continue;
             bool MuId=( (muIDbit->at(imu) >> 1 & 1)  && fabs(muD0->at(imu)) < 0.045 && fabs(muDz->at(imu)) < 0.2);
            LeadMu4Momentum.SetPtEtaPhiM(muPt->at(imu),muEta->at(imu),muPhi->at(imu),MuMass);
            
            
            for (int jmu = imu+1; jmu < nMu; ++jmu){
	       float IsoLep2Value=muPFChIso->at(jmu)/muPt->at(jmu);
        if ( (muPFNeuIso->at(jmu) + muPFPhoIso->at(jmu) - 0.5* muPFPUIso->at(jmu) )  > 0.0)
            IsoLep2Value= ( muPFChIso->at(jmu) + muPFNeuIso->at(jmu) + muPFPhoIso->at(jmu) - 0.5* muPFPUIso->at(jmu))/muPt->at(jmu);
            if (muPt->at(jmu) < 30 || muEta->at(jmu) > 2.4) continue;
             bool MuId2=( (muIDbit->at(jmu) >> 1 & 1)  && fabs(muD0->at(jmu)) < 0.045 && fabs(muDz->at(jmu)) < 0.2);

		
                
                SubMu4Momentum.SetPtEtaPhiM(muPt->at(jmu),muEta->at(jmu),muPhi->at(jmu),MuMass);

                numDiMu++;
		if (!(IsoLep1Value && IsoLep2Value && MuId && MuId2)) continue;
		ZCandidate = SubMu4Momentum + LeadMu4Momentum;
		bool OS = muCharge->at(imu) * muCharge->at(jmu) < 0;
		bool SS = muCharge->at(imu) * muCharge->at(jmu) > 0;
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
		  diMu_OS->Fill(ZCandidate.M());
		  break;
		}
		if (SS){
		  diMu_SS->Fill(ZCandidate.M());
		  break;
		}
		break; 
            }
        }
        
        if(numDiMu < 1) continue;

        
        
    }
    TH1::SetDefaultSumw2();
    int i;
    TH1F* osEventsVsRun = new TH1F("OSrunMu", "OSrunMu", 50, 315200, 325200);
    TH1F* ssEventsVsRun = new TH1F("SSrunMu", "SSrunMu", 50, 315200, 325200);
    for (i = 0; i < osRunSize; i++) {
      osEventsVsRun->Fill(osRun[i]);
    }
    for (i = 0; i < ssRunSize; i++) {
      ssEventsVsRun->Fill(ssRun[i]);
    } 
    osEventsVsRun->Write();
    ssEventsVsRun->Write();
    diMu_OS->Write();
    diMu_SS->Write();
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

