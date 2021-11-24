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

Double_t N_gen_tot[4] = { 0 };  // total number of gen events
Double_t N_gen_rap1[4] = { 0 }; // number of gen events with |y| < 1.0
Double_t N_gen_rap2[4] = { 0 }; // number of gen events with |y| < 0.8
Double_t Ratio1[4] = { 0 };      // ratio of N_gen_rap2[i] / N_gen_tot[i]
Double_t Ratio2[4] = { 0 };      // ratio of N_gen_rap2[i] / N_gen_rap1[i]

// 1): for |y| < 1.0, 2): for |y| < 0.8
Double_t SigSL_j_inc[2] = {5.247, 0}; //mb
Double_t SigSL_j_coh[2] = {12.504, 0};//mb
Double_t SigSL_p_inc[2] = {0.923, 0}; //mb
Double_t SigSL_p_coh[2] = {2.526, 0}; //mb

void CalculateNewCS(Int_t iMC);

void STARlight_CrossSections(){

    CalculateNewCS(0);

    /*
    for(Int_t i = 0; i < 4; i++){
        CalculateNewCS(i);
    }
    */

    return;
}

void CalculateNewCS(Int_t iMC){

    TFile *file = NULL;
    switch(iMC){
        case 0:
            file = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohJpsiToMu.root", "read");
            break;
        case 1:
            file = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohJpsiToMu.root", "read");
            break;
        case 2:
            file = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohPsi2sToMuPi.root", "read");
            break;
        case 3:
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
    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < tGen->GetEntries(); iEntry++){
        tGen->GetEntry(iEntry);

        hRapDist->Fill(fYGen);
        N_gen_tot[iMC]++;
        
        if(TMath::Abs(fYGen) < 1.0) N_gen_rap1[iMC]++;

        if(TMath::Abs(fYGen) < 0.8) N_gen_rap2[iMC]++;

        if((iEntry+1) % 200000 == 0){
            nEntriesAnalysed += 200000;
            //Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }

    Ratio1[iMC] = N_gen_rap2[iMC] / N_gen_tot[iMC];
    Ratio2[iMC] = N_gen_rap2[iMC] / N_gen_rap1[iMC];

    Printf("Done.");
    Printf("N_tot: %.0f", N_gen_tot[iMC]); 
    Printf("N_rap1 (with |y| < 1.0): %.0f", N_gen_rap1[iMC]); 
    Printf("N_rap2 (with |y| < 0.8): %.0f", N_gen_rap2[iMC]); 
    Printf("Ratio N_rap2/N_tot: %.4f", Ratio1[iMC]);
    Printf("Ratio N_rap2/N_rap1: %.4f", Ratio2[iMC]);

    hRapDist->Draw();

    return;
}