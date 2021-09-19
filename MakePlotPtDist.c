// MakePlotPtDist.c
// David Grund, Sep 9, 2021
// To plot a pt distribution of measured events 

#include "TFile.h"
#include "TH1.h"

#include "AnalysisManager.h"

void MakePlotPtDist(){

    TFile *fFileIn = TFile::Open("Trees/AnalysisData/AnalysisResultsLHC18qrMerged.root", "read");
    if(fFileIn) Printf("Input data loaded.");

    TTree *fTreeIn = dynamic_cast<TTree*> (fFileIn->Get("AnalysisOutput/fTreeJPsi"));
    if(fTreeIn) Printf("Input tree loaded.");

    ConnectTreeVariables(fTreeIn);

    // Create new data tree with applied cuts
    TFile fFileOut("Trees/MakePlotPtDist/MakePlotPtDist.root","RECREATE");

    TTree *tPtDist = new TTree("tPtDist", "tPtDist");
    tPtDist->Branch("fPt", &fPt, "fPt/D");
    tPtDist->Branch("fM", &fM, "fM/D");

    Int_t n_bins = 100;
    Double_t fPtMin = 0.2;
    Double_t fPtMax = 2.0;
    TH1D *hPtDist = new TH1D("hPtDist", "hPtDist", n_bins, fPtMin, fPtMax);

    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < fTreeIn->GetEntries(); iEntry++){
        fTreeIn->GetEntry(iEntry);
        if(EventPassed(0,-1)){
            hPtDist->Fill(fPt);
            tPtDist->Fill();
        }
        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }

    fFileOut.Write("",TObject::kWriteDelete);

    return;
}