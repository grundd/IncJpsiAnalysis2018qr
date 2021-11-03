// ZNClasses.c
// David Grund, Nov 3, 2021

// root headers
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
// my headers
#include "AnalysisManager.h"

void PrepareData();
void CalculateClasses();
Bool_t EventPassedZN(Int_t iMassCut, Int_t iPtCut, Int_t iPtBin);

void ZNClasses(){

    //PrepareData();

    CalculateClasses();

    return;
}

void CalculateClasses(){

    TFile *fFileIn = TFile::Open("Trees/ZNClasses/ZNClasses.root", "read");
    if(fFileIn) Printf("Input data loaded.");

    TTree *fTreeIn = dynamic_cast<TTree*> (fFileIn->Get("tZNClasses"));
    if(fTreeIn) Printf("Input tree loaded.");

    fTreeIn->SetBranchAddress("fPt", &fPt);
    fTreeIn->SetBranchAddress("fM", &fM);
    fTreeIn->SetBranchAddress("fZNA_energy", &fZNA_energy);
    fTreeIn->SetBranchAddress("fZNC_energy", &fZNC_energy);
    fTreeIn->SetBranchAddress("fZNA_time", &fZNA_time);
    fTreeIn->SetBranchAddress("fZNC_time", &fZNC_time); 

    Printf("%lli entries found in the tree.", fTreeIn->GetEntries());
    Int_t nEntriesAnalysed = 0;     

    SetPtBinning();

    for(Int_t iEntry = 0; iEntry < fTreeIn->GetEntries(); iEntry++){
        fTreeIn->GetEntry(iEntry);

        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            if(EventPassedZN(0, 1, iBin+1));
        }
                
        if((iEntry+1) % 500 == 0){
            nEntriesAnalysed += 500;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }

    return;
}

void PrepareData(){

    TFile *fFileIn = TFile::Open("Trees/AnalysisData/AnalysisResultsLHC18qrMerged.root", "read");
    if(fFileIn) Printf("Input data loaded.");

    TTree *fTreeIn = dynamic_cast<TTree*> (fFileIn->Get("AnalysisOutput/fTreeJPsi"));
    if(fTreeIn) Printf("Input tree loaded.");  

    ConnectTreeVariables(fTreeIn);  

    TFile fFileOut("Trees/ZNClasses/ZNClasses.root","RECREATE");

    TTree *tZNClasses = new TTree("tZNClasses", "tZNClasses");
    tZNClasses->Branch("fPt", &fPt, "fPt/D");
    tZNClasses->Branch("fM", &fM, "fM/D");
    tZNClasses->Branch("fZNA_energy", &fZNA_energy, "fZNA_energy/D");
    tZNClasses->Branch("fZNC_energy", &fZNC_energy, "fZNC_energy/D");
    tZNClasses->Branch("fZNA_time", &fZNA_time[0], "fZNA_time[4]/D");
    tZNClasses->Branch("fZNC_time", &fZNC_time[0], "fZNC_time[4]/D");

    Printf("%lli entries found in the tree.", fTreeIn->GetEntries());
    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < fTreeIn->GetEntries(); iEntry++){
        fTreeIn->GetEntry(iEntry);
        // iMassCut == 0    => 2.2 < m < 4.5 GeV
        // iPtCut == 3      => 0.2 < pT < 1.0 GeV
        if(EventPassed(0,3)) tZNClasses->Fill();

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }

    fFileOut.Write("",TObject::kWriteDelete);

    return;
}

Bool_t EventPassedZN(Int_t iMassCut, Int_t iPtCut, Int_t iPtBin){

    // Inv. mass cut
    Bool_t bMassCut = kFALSE;
    switch(iMassCut){
        case 0:
            if(fM > 3.0 && fM < 3.2) bMassCut = kTRUE;
            break;
        case 1:
            if((fM > 2.5 && fM < 2.9) || (fM > 3.3 && fM < 4.0)) bMassCut = kTRUE;
            break;
    }
    if(!bMassCut) return kFALSE;

    // Transverse momentu cut
    Bool_t bPtCut = kFALSE;
    switch(iPtBin){
        case 0:
            if(fPt > 0.2 && fPt < 1.0) bPtCut = kTRUE;
            break;
        case 1:
            if(fPt > ptBoundaries[iPtBin-1] && fPt <= ptBoundaries[iPtBin]) bPtCut = kTRUE;
            break;
    }
    if(!bPtCut) return kFALSE;

    // 10) Check ZN times: at least one time has to fall within pm 2 ns
    Bool_t ZN_times = kFALSE;
    for(Int_t i = 0; i < 4; i++){
        if(TMath::Abs(fZNA_time[i]) < 2 || TMath::Abs(fZNC_time[i]) < 2) ZN_times = kTRUE;
    }
    if(!ZN_times) return kFALSE;

    // Event passed all the selections =>
    return kTRUE;
}