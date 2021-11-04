// ZNClasses.c
// David Grund, Nov 3, 2021

// cpp headers
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
// root headers
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TString.h"
// my headers
#include "AnalysisManager.h"

void PrepareData();
void CalculateClasses(Int_t iMassCut);
Bool_t EventPassedZN(Int_t iMassCut, Int_t iPtCut, Int_t iPtBin);
Double_t CalculateErrorBayes(Double_t k, Double_t n);

Double_t nEvTotal[nPtBins] = { 0 };        // total number of events per bin
Double_t nEvFired[nPtBins] = { 0 };        // number of events per bin with at least one fired ZN
Double_t nEvFired_ZNA[nPtBins] = { 0 };    // number of events per bin with fired ZNA
Double_t nEvFired_ZNC[nPtBins] = { 0 };    // number of events per bin with fired ZNC
Double_t nEvFired_OnlyZNA[nPtBins] = { 0 };// number of events per bin with fired ZNA only
Double_t nEvFired_OnlyZNC[nPtBins] = { 0 };// number of events per bin with fired ZNC only
Double_t nEvFired_both[nPtBins] = { 0 };   // number of events per bin with fired both ZN
Double_t nEvFired_none[nPtBins] = { 0 };   // number of events per bin with no fired ZN
Double_t nSums[8] = { 0 };

// Percentages (relative to the total number of events in a bin)
Double_t PercFiredBoth_val[nPtBins] = { 0 }; 
Double_t PercFiredNone_val[nPtBins] = { 0 }; 
Double_t PercFiredZNAOnly_val[nPtBins] = { 0 }; 
Double_t PercFiredZNCOnly_val[nPtBins] = { 0 }; 
Double_t PercFiredBoth_err[nPtBins] = { 0 }; 
Double_t PercFiredNone_err[nPtBins] = { 0 }; 
Double_t PercFiredZNAOnly_err[nPtBins] = { 0 }; 
Double_t PercFiredZNCOnly_err[nPtBins] = { 0 }; 

void ZNClasses(){

    //PrepareData();

    CalculateClasses(0); // Around J/psi peak: 3.0 < m < 3.2 GeV

    CalculateClasses(1); // Side-bands: 2.5 < m < 2.9 GeV && 3.3 < m < 4.0 GeV

    return;
}

void CalculateClasses(Int_t iMassCut){

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
            if(EventPassedZN(iMassCut, 1, iBin+1)){
                Bool_t ZN_fired = kFALSE;
                Bool_t ZNA_fired = kFALSE;
                Bool_t ZNC_fired = kFALSE;
                for(Int_t i = 0; i < 4; i++){
                    // Check ZN times: at least one time has to fall within pm 2 ns
                    if(TMath::Abs(fZNA_time[i]) < 2 || TMath::Abs(fZNC_time[i]) < 2) ZN_fired = kTRUE;
                    if(TMath::Abs(fZNA_time[i]) < 2) ZNA_fired = kTRUE;
                    if(TMath::Abs(fZNC_time[i]) < 2) ZNC_fired = kTRUE;
                }
                nEvTotal[iBin]++;
                if(ZN_fired) nEvFired[iBin]++;
                if(ZNA_fired && ZNC_fired) nEvFired_both[iBin]++;
                if(ZNA_fired) nEvFired_ZNA[iBin]++;
                if(ZNC_fired) nEvFired_ZNC[iBin]++;
                if(ZNA_fired && !ZNC_fired) nEvFired_OnlyZNA[iBin]++;
                if(ZNC_fired && !ZNA_fired) nEvFired_OnlyZNC[iBin]++;
                if(!ZNA_fired && !ZNC_fired) nEvFired_none[iBin]++;
            }
        }

        if((iEntry+1) % 500 == 0){
            nEntriesAnalysed += 500;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }

    // Calculate percentages
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        PercFiredBoth_val[iBin] = nEvFired_both[iBin] / nEvTotal[iBin] * 100;
        PercFiredBoth_err[iBin] = CalculateErrorBayes(nEvFired_both[iBin], nEvTotal[iBin]) * 100;
        PercFiredNone_val[iBin] = nEvFired_none[iBin] / nEvTotal[iBin] * 100;
        PercFiredNone_err[iBin] = CalculateErrorBayes(nEvFired_none[iBin], nEvTotal[iBin]) * 100;
        PercFiredZNAOnly_val[iBin] = nEvFired_OnlyZNA[iBin] / nEvTotal[iBin] * 100;
        PercFiredZNAOnly_err[iBin] = CalculateErrorBayes(nEvFired_OnlyZNA[iBin], nEvTotal[iBin]) * 100;
        PercFiredZNCOnly_val[iBin] = nEvFired_OnlyZNC[iBin] / nEvTotal[iBin] * 100;
        PercFiredZNCOnly_err[iBin] = CalculateErrorBayes(nEvFired_OnlyZNC[iBin], nEvTotal[iBin]) * 100;
    }

    // Print the results to console
    // Numbers
    Printf("Numbers of events:");
    Printf("Bin \tTotal \tFired \tZNA \tZNC \tOnlyZNA\tOnlyZNC\tBoth \tNone");
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        nSums[0] += nEvTotal[iBin];
        nSums[1] += nEvFired[iBin];
        nSums[2] += nEvFired_ZNA[iBin];
        nSums[3] += nEvFired_ZNC[iBin];
        nSums[4] += nEvFired_OnlyZNA[iBin];
        nSums[5] += nEvFired_OnlyZNC[iBin];
        nSums[6] += nEvFired_both[iBin];
        nSums[7] += nEvFired_none[iBin];
        Printf("%i \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f", 
            iBin+1, 
            nEvTotal[iBin], 
            nEvFired[iBin], 
            nEvFired_ZNA[iBin], 
            nEvFired_ZNC[iBin], 
            nEvFired_OnlyZNA[iBin], 
            nEvFired_OnlyZNC[iBin], 
            nEvFired_both[iBin],
            nEvFired_none[iBin]);
    }
    Printf("sum \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f", nSums[0], nSums[1], nSums[2], nSums[3], nSums[4], nSums[5], nSums[6], nSums[7]);
    // Percentages
    Printf("Percentages:");
    Printf("Bin \tBoth[%%]\terr \tZNA[%%]\terr \tZNC[%%]\terr \tNone[%%]\terr");
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        Printf("%i \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f", 
            iBin+1,
            PercFiredBoth_val[iBin],
            PercFiredBoth_err[iBin],
            PercFiredZNAOnly_val[iBin],
            PercFiredZNAOnly_err[iBin],
            PercFiredZNCOnly_val[iBin],
            PercFiredZNCOnly_err[iBin],
            PercFiredNone_val[iBin],
            PercFiredNone_err[iBin]);
    }

    // Print numbers to txt files
    TString FilePath = Form("Results/ZNClasses/%ibins_MassInterval%i_output.txt", nPtBins, iMassCut);
    ofstream outfile(FilePath.Data());
    outfile << Form("Numbers of events:\n");
    outfile << Form("Bin \tTotal \tFired \tZNA \tZNC \tOnlyZNA\tOnlyZNC\tBoth \tNone\n");
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        outfile << Form("%i \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \n", 
            iBin+1, 
            nEvTotal[iBin], 
            nEvFired[iBin], 
            nEvFired_ZNA[iBin], 
            nEvFired_ZNC[iBin], 
            nEvFired_OnlyZNA[iBin], 
            nEvFired_OnlyZNC[iBin], 
            nEvFired_both[iBin],
            nEvFired_none[iBin]);
    }
    outfile << Form("sum \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \n", nSums[0], nSums[1], nSums[2], nSums[3], nSums[4], nSums[5], nSums[6], nSums[7]);
    outfile << Form("\n\n");
    outfile << Form("Percentages:\n");
    outfile << Form("Bin \tBoth[%%]\terr \tZNA[%%]\terr \tZNC[%%]\terr \tNone[%%]\terr \n");
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        outfile << Form("%i \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \t%.1f \n", 
            iBin+1,
            PercFiredBoth_val[iBin],
            PercFiredBoth_err[iBin],
            PercFiredZNAOnly_val[iBin],
            PercFiredZNAOnly_err[iBin],
            PercFiredZNCOnly_val[iBin],
            PercFiredZNCOnly_err[iBin],
            PercFiredNone_val[iBin],
            PercFiredNone_err[iBin]);
    }
    outfile.close();
    Printf("Results printed to %s.", FilePath.Data()); 

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
    switch(iPtCut){
        case 0:
            if(fPt > 0.2 && fPt < 1.0) bPtCut = kTRUE;
            break;
        case 1:
            if(fPt > ptBoundaries[iPtBin-1] && fPt <= ptBoundaries[iPtBin]) bPtCut = kTRUE;
            break;
    }
    if(!bPtCut) return kFALSE;

    // Event passed all the selections =>
    return kTRUE;
}

Double_t CalculateErrorBayes(Double_t k, Double_t n){ // k = nominator, n = denominator

    Double_t var = (k + 1) * (k + 2) / (n + 2) / (n + 3) - (k + 1) * (k + 1) / (n + 2) / (n + 2);
    Double_t err = TMath::Sqrt(var);

    return err;
}