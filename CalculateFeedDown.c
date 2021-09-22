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

void CalculateFeedDownInBins();
Bool_t AxE_PtBins_Jpsi();
void CalculateAxE_PtBins_Psi2s(Int_t iFD);
// iFD:
// == 1 => coh charged
// == 2 => inc charged
// == 3 => coh neutral
// == 4 => inc neutral
Double_t CalculateErrorBayes(Double_t k, Double_t n);

TString Datasets[4] = {"CohCharged","IncCharged","CohNeutral","IncNeutral"};
// Temporary variables when loading the data:
Int_t i_bin;
Double_t pt_low, pt_upp;
// Arrays of AxE in bins
Double_t AxE_Jpsi_val[nPtBins] = { 0 };
Double_t AxE_Jpsi_err[nPtBins] = { 0 };
Double_t AxE_Psi2s_val[nPtBins] = { 0 };
Double_t AxE_Psi2s_err[nPtBins] = { 0 };
Double_t CorrFD[4][nPtBins] = { 0 };
Double_t CorrFD_err[4][nPtBins] = { 0 };
// Branching ratios
Double_t BR;
Double_t BR_err;
Double_t BR_ch = 0.3468;
Double_t BR_ch_err = 0.0030;
Double_t BR_ne = 0.2538;
Double_t BR_ne_err = 0.0032;
// STARlight cross sections for |y| < 1
// Needed for |y| < 0.8 (!)
Double_t SigSL_p;
Double_t SigSL_j_inc = 5.247; //mb
Double_t SigSL_j_coh = 12.504;//mb
Double_t SigSL_p_inc = 0.92;  //mb
Double_t SigSL_p_coh = 2.52;  //mb

void CalculateFeedDown(){

    CalculateFeedDownInBins();

    return;
}

void CalculateFeedDownInBins(){

    SetPtBinning();

    // Calculate AxE for kIncohJpsiToMu
    Bool_t AxE_Jpsi = AxE_PtBins_Jpsi();
    if(!AxE_Jpsi){
        Printf("Values of AxE for kIncohJpsiToMu not found!");
        Printf("Create them using AccAndEffMC.c.");
        Printf("Terminating...");
        return;
    }
    // Calculate the FD correction in all bins
    for(Int_t iFD = 1; iFD <= 4; iFD++){
        // Calculate AxE for Psi2s files
        CalculateAxE_PtBins_Psi2s(iFD);
        // Go over pt bins
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            // Cross-check: print the values that will be used for calculations
            Printf("Bin %i:", iBin+1);
            Printf("AxE J/psi:\t(%.3f pm %.3f)%%", AxE_Jpsi_val[iBin], AxE_Jpsi_err[iBin]);
            Printf("AxE Psi2s:\t(%.3f pm %.3f)%%", AxE_Psi2s_val[iBin], AxE_Psi2s_err[iBin]);
            
            // Load BR and cross sections
            if(iFD == 1 || iFD == 2){
                // charged processes
                BR = BR_ch;
                BR_err = BR_ch_err;
            }
            if(iFD == 3 || iFD == 4){
                // neutral processes
                BR = BR_ne;
                BR_err = BR_ne_err;
            }
            if(iFD == 1 || iFD == 3){
                // coherent Psi2s
                SigSL_p = SigSL_p_coh;
            }
            if(iFD == 2 || iFD == 4){
                // incoherent Psi2s
                SigSL_p = SigSL_p_inc;
            }   
            // Calculate FD corrections         
            CorrFD[iFD-1][iBin] = SigSL_p / SigSL_j_inc * AxE_Psi2s_val[iBin] / AxE_Jpsi_val[iBin] * BR * 100; // in percent already
            CorrFD_err[iFD-1][iBin] = CorrFD[iFD-1][iBin] * TMath::Sqrt(TMath::Power((AxE_Psi2s_err[iBin]/AxE_Psi2s_val[iBin]),2) + TMath::Power((AxE_Jpsi_err[iBin]/AxE_Jpsi_val[iBin]),2) + TMath::Power((BR_err/BR),2));
            Printf("FD corr: \t(%.3f pm %.3f)%%", CorrFD[iFD-1][iBin], CorrFD_err[iFD-1][iBin]);
            Printf(" ");
        }  
    } 
    // Define an output text file
    TString OutputFile(Form("Results/FeedDown/FeedDownCorrections_%ibins.txt", nPtBins));
    ofstream outfile(OutputFile.Data());
    outfile << std::fixed << std::setprecision(3);
    // Print the results to the output text file
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        outfile << iBin + 1 << "\t"
                << CorrFD[0][iBin] << "\t"      // coh charged
                << CorrFD_err[0][iBin] << "\t"
                << CorrFD[1][iBin] << "\t"      // inc charged
                << CorrFD_err[1][iBin] << "\t"
                << CorrFD[2][iBin] << "\t"      // coh neutral
                << CorrFD_err[2][iBin] << "\t"
                << CorrFD[3][iBin] << "\t"      // inc neutral
                << CorrFD_err[3][iBin] << "\n";
    }
    // Close the output file
    outfile.close();
    Printf("Results printed to %s.", OutputFile.Data());   

    return;
}

Bool_t AxE_PtBins_Jpsi(){

    TString FilePath = Form("Results/AccAndEffMC/AxE_%ibins.txt", nPtBins);
    // Check if already calculated
    ifstream FileIn;
    FileIn.open(FilePath.Data());
    if(!(FileIn.fail())){
        // Read data from the file
        Int_t i_line = 0;
        while(!FileIn.eof()){
            FileIn >> i_bin >> AxE_Jpsi_val[i_line] >> AxE_Jpsi_err[i_line];
            i_line++;
        }
        FileIn.close(); 
        Printf("AxE in %ibins for kIncohJpsiToMu has already been calculated.\n", nPtBins);

        return kTRUE;
    } else {
        // AxE for kIncohJpsiToMu not found
        return kFALSE;
    }
}

void CalculateAxE_PtBins_Psi2s(Int_t iFD){

    TString FilePath = Form("Results/FeedDown/AxE_%s_%ibins", (Datasets[iFD-1]).Data(), nPtBins);
    // Check if already calculated
    ifstream FileIn;
    FileIn.open((FilePath + ".txt").Data());
    if(!(FileIn.fail())){
        // Read data from the file
        Int_t i_line = 0;
        while(!FileIn.eof()){
            FileIn >> i_bin >> pt_low >> pt_upp >> AxE_Psi2s_val[i_line] >> AxE_Psi2s_err[i_line];
            i_line++;
        }
        FileIn.close(); 
        Printf("############################");
        Printf("AxE in %ibins for %s has already been calculated.\n", nPtBins, (Datasets[iFD-1]).Data());

        return;
    }

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
    for(Int_t iPtBin = 0; iPtBin < nPtBins; iPtBin++){
        AxE_Psi2s_val[iPtBin] = NRec[iPtBin] / NGen[iPtBin];
        AxE_Psi2s_err[iPtBin] = CalculateErrorBayes(NRec[iPtBin], NGen[iPtBin]);
        Printf("AxE = (%.4f pm %.4f)%%", AxE_Psi2s_val[iPtBin]*100, AxE_Psi2s_err[iPtBin]*100);
    }

    // Save the results to text file
    ofstream outfile((FilePath + ".txt").Data());
    outfile << std::fixed << std::setprecision(3);
    //outfile << Form("Bin \tPtLow \tPtUpp \tAxE_Psi2s_val [%%]\tAxE_Psi2s_err [%%] \n");
    for(Int_t i = 1; i <= nPtBins; i++){
        outfile << i << "\t" << ptBoundaries[i-1] << "\t" << ptBoundaries[i] << "\t" << AxE_Psi2s_val[i-1]*100 << "\t" << AxE_Psi2s_err[i-1]*100 << "\n";
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