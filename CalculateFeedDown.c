// CalculateFeedDown.c
// David Grund, Sep 21, 2021

// cpp headers
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
// root headers
#include "TFile.h"
#include "TMath.h"
// my headers
#include "AnalysisManager.h"


#include<iostream>
using namespace std;

void CalculateAxE_PtBins(Int_t iFD);
// iFD:
// == 1 => coh charged
// == 2 => inc charged
// == 2 => coh neutral
// == 3 => inc neutral
Double_t CalculateErrorBayes(Double_t k, Double_t n);

// Temporary variables when loading the data:
Int_t i_bin;
Double_t pt_low, pt_upp;
Double_t axe_val, axe_err;

void CalculateFeedDown(){

    SetPtBinning();

    for(Int_t i = 1; i <= 4; i++) CalculateAxE_PtBins(i);

    return;
}

void CalculateAxE_PtBins(Int_t iFD){

    TFile *fRec = NULL;
    switch(iFD){
        case 1:
            fRec = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohPsi2sToMuPi.root", "read");
            break;
        case 2:
            fRec = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohPsi2sToMuPi.root", "read");
            break;
        case 3:
            fRec = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohPsi2sToMuPi_neutral.root", "read");
            break;
        case 4:
            fRec = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohPsi2sToMuPi_neutral.root", "read");
            break;
    }
    if(fRec) Printf("MC rec file loaded.");

    TTree *tRec = dynamic_cast<TTree*> (fRec->Get("AnalysisOutput/fTreeJPsiMCRec"));
    if(tRec) Printf("MC rec tree loaded.");
    
    ConnectTreeVariablesMCRec(tRec);

    Printf("Now calculating FD for %s.", fRec->GetName());

    TString Datasets[4] = {"CohCharged","IncCharged","CohNeutral","IncNeutral"};
    TString FilePath = Form("Results/FeedDown/AxE_%s_%ibins", (Datasets[iFD-1]).Data(), nPtBins);
    // Check if already calculated
    ifstream FileIn;
    FileIn.open((FilePath + ".txt").Data());
    TString s[5] = {""};
    if(!(FileIn.fail())){
        // Read data from the file
        while(!FileIn.eof()){
            FileIn >> i_bin >> pt_low >> pt_upp >> axe_val >> axe_err;
        }
        FileIn.close(); 
        Printf("AxE in %ibins for %s has already been calculated.", nPtBins, fRec->GetName());

        return;
    }

    Printf("Tree %s has %lli entries.", tRec->GetName(), tRec->GetEntries());

    Double_t NRec[nPtBins] = { 0 };
    // Loop over pt bins
    for(Int_t iPtBin = 1; iPtBin <= nPtBins; iPtBin++){
        // Loop over tree entries
        Printf("Bin %i:", iPtBin);
        Int_t nEntriesAnalysed = 0;
        for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++){
            tRec->GetEntry(iEntry);
            if(EventPassedMCRec(0, 4, iPtBin)) NRec[iPtBin-1]++;
            if((iEntry+1) % 200000 == 0){
                nEntriesAnalysed += 200000;
                Printf("%i entries analysed.", nEntriesAnalysed);
            }
        }
    }
    Printf("Done.");
    
    TTree *tGen = dynamic_cast<TTree*> (fRec->Get("AnalysisOutput/fTreeJPsiMCGen"));
    if(tGen) Printf("MC gen tree loaded.");
    
    ConnectTreeVariablesMCGen(tGen);

    Printf("Tree %s has %lli entries.", tGen->GetName(), tGen->GetEntries());

    Double_t NGen[nPtBins] = { 0 };
    // Loop over pt bins
    for(Int_t iPtBin = 1; iPtBin <= nPtBins; iPtBin++){
        // Loop over tree entries
        Printf("Bin %i:", iPtBin);
        Int_t nEntriesAnalysed = 0;
        for(Int_t iEntry = 0; iEntry < tGen->GetEntries(); iEntry++){
            tGen->GetEntry(iEntry);
            if(EventPassedMCGen(4, iPtBin)) NGen[iPtBin-1]++;
            if((iEntry+1) % 500000 == 0){
                nEntriesAnalysed += 500000;
                Printf("%i entries analysed.", nEntriesAnalysed);
            }
        }
    }
    Printf("Done.");

    // Calculate AxE per bin
    Double_t AxE[nPtBins] = { 0 };
    Double_t AxE_err[nPtBins] = { 0 };
    for(Int_t iPtBin = 0; iPtBin < nPtBins; iPtBin++){
        AxE[iPtBin] = NRec[iPtBin] / NGen[iPtBin];
        AxE_err[iPtBin] = CalculateErrorBayes(NRec[iPtBin], NGen[iPtBin]);
        Printf("AxE = (%.4f pm %.4f)%%", AxE[iPtBin]*100, AxE_err[iPtBin]*100);
    }

    // Save the results to text file
    ofstream outfile((FilePath + ".txt").Data());
    outfile << std::fixed << std::setprecision(3);
    //outfile << Form("Bin \tPtLow \tPtUpp \tAxE [%%]\tAxE_err [%%] \n");
    for(Int_t i = 1; i <= nPtBins; i++){
        outfile << i << "\t" << ptBoundaries[i-1] << "\t" << ptBoundaries[i] << "\t" << AxE[i-1]*100 << "\t" << AxE_err[i-1]*100 << "\n";
    }
    outfile.close();
    Printf("Results printed to %s.", (FilePath + ".txt").Data());    

    return;
}

Double_t CalculateErrorBayes(Double_t k, Double_t n){ // k = NRec, n = NGen

    Double_t var = (k + 1) * (k + 2) / (n + 2) / (n + 3) - (k + 1) * (k + 1) / (n + 2) / (n + 2);
    Double_t err = TMath::Sqrt(var);

    return err;
}