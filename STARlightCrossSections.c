// STARlightCrossSections.c
// David Grund, Sep 23, 2021
// To find STARlight cross section for various processes for |y| < 0.8 instead of |y| < 1.0

// (... still has to be finished ...)

// cpp headers
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
// root headers
#include "TFile.h"
#include "TMath.h"
#include "TH1.h"
// my headers
#include "AnalysisManager.h"

Int_t nBins = 200;

// For |y| < 1.0
Double_t SigSL_j_inc = 5.247; //mb
Double_t SigSL_j_coh = 12.504;//mb
Double_t SigSL_p_inc = 0.92;  //mb
Double_t SigSL_p_coh = 2.52;  //mb
// For |y| < 0.8
Double_t SigSL_j_inc2;   //mb
Double_t SigSL_j_coh2;   //mb
Double_t SigSL_p_inc2;   //mb
Double_t SigSL_p_coh2;   //mb

void CalculateNewCS(Int_t Dataset);

void STARlightCrossSections(){

    for(Int_t i = 1; i <= 4; i++){
        CalculateNewCS(i);
    }

    return;
}

void CalculateNewCS(Int_t Dataset){

    TFile *file = NULL;
    switch(Dataset){
        case 1:
            file = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohJpsiToMu.root", "read");
            break;
        case 2:
            file = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohJpsiToMu.root", "read");
            break;
        case 3:
            file = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohPsi2sToMuPi.root", "read");
            break;
        case 4:
            file = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohPsi2sToMuPi.root", "read");
            break;
    }
    if(file) Printf("File %s loaded.", file->GetName());

    TTree *tGen = dynamic_cast<TTree*> (file->Get("AnalysisOutput/fTreeJPsiMCGen"));
    if(tGen) Printf("MC gen tree loaded.");
    
    ConnectTreeVariablesMCGen(tGen);

    Printf("Tree %s has %lli entries.", tGen->GetName(), tGen->GetEntries());

    TH1D *hRapDist = new TH1D("hRapDist","hRapDist",nBins,-2,2);

    // Loop over tree entries
    Int_t N_tot = 0; // total number of gen events
    Int_t N_y10 = 0; // number of gen events with |y| < 1.0
    Int_t N_y08 = 0; // number of gen events with |y| < 0.8
    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < tGen->GetEntries(); iEntry++){
        tGen->GetEntry(iEntry);

        hRapDist->Fill(fYGen);
        N_tot++;
        
        if(TMath::Abs(fYGen) <= 1.0) N_y10++;
        //else Printf("Event outside |y| <= 1.0.");

        if(TMath::Abs(fYGen) < 0.8) N_y08++;

        if((iEntry+1) % 200000 == 0){
            nEntriesAnalysed += 200000;
            //Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }
    Printf("Done.");
    Printf("N tot: %i", N_tot); 
    Printf("N with |y| <= 1.0: %i", N_y10); 
    Printf("N with |y| < 0.8: %i", N_y08); 

    hRapDist->Draw();

    return;
}