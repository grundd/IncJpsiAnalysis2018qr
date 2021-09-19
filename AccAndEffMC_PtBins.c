// AccAndEffMC_PtBins.c
// David Grund, 15-09-2021
// To calculate the acceptance x efficiency from MC data in defined pt bins

// cpp headers
#include <fstream> // print output to txt file
//#include <iomanip> // std::setprecision()
// root headers
#include "TH1.h"
#include "TString.h"
#include "TFile.h"
// my headers
#include "AnalysisManager.h"

TH1D* hNRec = new TH1D("hNRec","N rec per bin",nPtBins,ptBoundaries);
TH1D* hNGen = new TH1D("hNGen","N gen per bin",nPtBins,ptBoundaries);
TH1D* hAxE = new TH1D("hAxE","AxE per bin",nPtBins,ptBoundaries);

void FillHistNRec();
void FillHistNGen();
void SaveToFile(TH1D* hist, TString name);

void AccAndEffMC_PtBins(){

    SetPtBinning();

    FillHistNRec();

    FillHistNGen();

    return;
}

void FillHistNRec(){
    // Check if the corresponding text file already exists
    TString file("Results/AccAndEffMC_PtBins/");
    if(nPtBins == 4) file.Append("NRec_4bins.txt");
    if(nPtBins == 5) file.Append("NRec_5bins.txt");

    ifstream inFile;
    inFile.open(file);
    if(!(inFile.fail())){
        // This configuration has already been calculated
        Printf("*** The file %s already exists. ***", file.Data());
        // Fill hNRec with data from the text file
        Int_t inBin;
        Double_t inValue;
        while(!inFile.eof()){
            inFile >> inBin >> inValue; // fist and second column
            hNRec->SetBinContent(inBin, inValue);
        }
        inFile.close(); 

        return;
    } else {
        // This configuration is yet to be calculated
        Printf("*** Calculating N rec per bin for %s... ***", file.Data());

        TFile *fRec = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohJpsiToMu.root", "read");
        if(fRec) Printf("MC rec file loaded.");

        TTree *tRec = dynamic_cast<TTree*> (fRec->Get("AnalysisOutput/fTreeJPsiMCRec"));
        if(tRec) Printf("MC rec tree loaded.");
        
        ConnectTreeVariablesMCRec(tRec);

        // Loop over all pt bins
        for(Int_t iPtBin = 1; iPtBin <= nPtBins; iPtBin++){
            Int_t NRec = 0;
            for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++){
                tRec->GetEntry(iEntry);
                if(EventPassedMCRec(0, 4, iPtBin)) NRec++;
            }
            hNRec->SetBinContent(iPtBin, NRec);
            Printf("*** Bin %i done. ***", iPtBin);
        }
        Printf("*** Finished! ***");
        
        SaveToFile(hNRec, file);

        return;
    }
}

void FillHistNGen(){
    // Check if the corresponding text file already exists
    TString file("Results/AccAndEffMC_PtBins/");
    if(nPtBins == 4) file.Append("NGen_4bins.txt");
    if(nPtBins == 5) file.Append("NGen_5bins.txt");

    ifstream inFile;
    inFile.open(file);
    if(!(inFile.fail())){
        // This configuration has already been calculated
        Printf("*** The file %s already exists. ***", file.Data());
        // Fill hNRec with data from the text file
        Int_t inBin;
        Double_t inValue;
        while(!inFile.eof()){
            inFile >> inBin >> inValue; // fist and second column
            hNRec->SetBinContent(inBin, inValue);
        }
        inFile.close(); 

        return;
    } else {
        // This configuration is yet to be calculated
        Printf("*** Calculating N gen per bin for %s... ***", file.Data());

        TFile *fRec = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohJpsiToMu.root", "read");
        if(fRec) Printf("MC rec file loaded.");

        TTree *tGen = dynamic_cast<TTree*> (fRec->Get("AnalysisOutput/fTreeJPsiMCGen"));
        if(tGen) Printf("MC rec tree loaded.");
        
        ConnectTreeVariablesMCGen(tGen);

        // Loop over all pt bins
        for(Int_t iPtBin = 1; iPtBin <= nPtBins; iPtBin++){
            Int_t NGen = 0;
            for(Int_t iEntry = 0; iEntry < tGen->GetEntries(); iEntry++){
                tGen->GetEntry(iEntry);
                if(EventPassedMCGen(0, iPtBin)) NGen++;
            }
            hNGen->SetBinContent(iPtBin, NGen);
            Printf("*** Bin %i done. ***", iPtBin);
        }
        Printf("*** Finished! ***");
        
        SaveToFile(hNGen, file);

        return;
    }
    return;
}

void SaveToFile(TH1D* hist, TString name){
    ofstream outfile (name.Data());
    for(Int_t iBin = 1; iBin <= hist->GetNbinsX(); iBin++){
        outfile << iBin << "\t" << hist->GetBinContent(iBin) << "\n";
    }
    outfile.close();
    Printf("*** File saved in %s.***", name.Data());
}