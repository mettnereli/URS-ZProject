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
            
    TH1F* hEvents = (TH1F*)gDirectory->Get("ggNtuplizer/hEvents");
    TH1F* hPU     = (TH1F*)gDirectory->Get("ggNtuplizer/hPU");
    TH1F* hPUTrue = (TH1F*)gDirectory->Get("ggNtuplizer/hPUTrue");
    TH1F* diMu_OS = new TH1F("diMu_OS","diMu_OS",30, 0, 160);
    TH1F* diMu_SS = new TH1F("diMu_SS","diMu_SS",30, 0, 160);
    
    TFile* file = TFile::Open(OutputFile, "RECREATE");
    TTree* MyNewTree = fChain->CloneTree(0);
    
    fChain->SetBranchStatus("*",1);
    
    TH1F* hcount = new TH1F("hcount", "", 10, 0, 10);
    
    if (fChain == 0) return;
    
    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    float MuMass= 0.10565837;
    float eleMass= 0.000511;
    
    
    
    for (int jentry=0; jentry<nentries;jentry++) {
        
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        
        
        
        if(jentry % 10000 == 0) cout << "Processed " << jentry << " events out of " <<nentries<<endl;
        
        hcount->Fill(1);
        if (!isData)
            hcount->Fill(2,genWeight);
        
        
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
        hcount->Fill(3);

        
        
        MyNewTree->Fill();
    }
    
    
    MyNewTree->AutoSave();
    hEvents->Write();
    hcount->Write();
    diMu_OS->Write();
    diMu_SS->Write();
    if (hPU) hPU->Write();
    if (hPUTrue) hPUTrue->Write();
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


