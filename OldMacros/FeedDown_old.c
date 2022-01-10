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

void CalculateFD_Total(Int_t iPtCut);
void CalculateFD_Total_AOD();
void CalculateFD_PtBins();
void CalculateFD_PtFit();
Bool_t AxE_IncJpsi_Total(Int_t iMC, Int_t iMassCut, Int_t iPtCut);
// iMC == 0 => kCohJpsiToMu
//     == 1 => kIncohJpsiToMu
Bool_t AxE_IncJpsi_Total_AOD();
Bool_t AxE_IncJpsi_PtBins();
void CalculateAxE_Psi2s_Total(Int_t iPtCut, Int_t iFD);
void CalculateAxE_Psi2s_Total_AOD(Int_t iFD);
void CalculateAxE_Psi2s_PtBins(Int_t iFD);
void CalculateAxE_Psi2s_PtFit(Int_t iFD);
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
// Total feed-down coefficient (for PtCut0, PtCut3 and PtFit)
Double_t AxE_Jpsi_tot_val = 0;
Double_t AxE_Jpsi_tot_err = 0;
Double_t AxE_Psi2s_tot_val = 0;
Double_t AxE_Psi2s_tot_err = 0;
Double_t CorrFD_tot_val[4] = { 0 };
Double_t CorrFD_tot_err[4] = { 0 };
// Arrays of AxE in bins
Double_t AxE_IncJpsi_val[nPtBins] = { 0 };
Double_t AxE_IncJpsi_err[nPtBins] = { 0 };
Double_t AxE_Psi2s_val[nPtBins] = { 0 };
Double_t AxE_Psi2s_err[nPtBins] = { 0 };
Double_t CorrFD_val[4][nPtBins] = { 0 };
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
Double_t SigSL_j;
Double_t SigSL_j_inc = 5.247; //mb
Double_t SigSL_j_coh = 12.504;//mb
Double_t SigSL_p_inc = 0.92;  //mb
Double_t SigSL_p_coh = 2.52;  //mb

void CalculateFeedDown(){

    //CalculateFD_Total(0);
    //CalculateFD_Total(3);

    //CalculateFD_Total_AOD();

    //CalculateFD_PtBins();

    CalculateFD_PtFit();

    return;
}

void CalculateFD_Total(Int_t iPtCut){

    // Load AxE value for kIncohJpsiToMu
    Bool_t AxE_IncJpsi = AxE_IncJpsi_Total(1, 0, iPtCut);
    if(!AxE_IncJpsi){
        Printf("Values of AxE for kIncohJpsiToMu with pt cut no. %i not found!", iPtCut);
        Printf("Create them using AccAndEffMC.c.");
        Printf("Terminating...");
        return;
    }
    // Calculate all FD corrections
    for(Int_t iFD = 1; iFD <= 4; iFD++){
        CalculateAxE_Psi2s_Total(iPtCut, iFD);
        // Cross-check: print the values that will be used for calculations
        Printf("AxE J/psi:\t(%.3f pm %.3f)%%", AxE_Jpsi_tot_val, AxE_Jpsi_tot_err);
        Printf("AxE Psi2s:\t(%.3f pm %.3f)%%", AxE_Psi2s_tot_val, AxE_Psi2s_tot_err);

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
        CorrFD_tot_val[iFD-1] = SigSL_p / SigSL_j_inc * AxE_Psi2s_tot_val / AxE_Jpsi_tot_val * BR * 100; // in percent already
        CorrFD_tot_err[iFD-1]= CorrFD_tot_val[iFD-1] * TMath::Sqrt(TMath::Power((AxE_Psi2s_tot_err/AxE_Psi2s_tot_val),2) + TMath::Power((AxE_Jpsi_tot_err/AxE_Jpsi_tot_val),2) + TMath::Power((BR_err/BR),2));
        Printf("FD corr: \t(%.3f pm %.3f)%%", CorrFD_tot_val[iFD-1], CorrFD_tot_err[iFD-1]);
        Printf(" ");
    }
    // Define an output text file
    TString OutputFile(Form("Results/FeedDown/FeedDownCorrections_tot_PtCut%i.txt", iPtCut));
    ofstream outfile(OutputFile.Data());
    outfile << std::fixed << std::setprecision(3);
    // Print the results to the output text file
    outfile << CorrFD_tot_val[0] << "\t"    // coh charged
            << CorrFD_tot_err[0] << "\t"
            << CorrFD_tot_val[1] << "\t"    // inc charged
            << CorrFD_tot_err[1] << "\t"
            << CorrFD_tot_val[2] << "\t"    // coh neutral
            << CorrFD_tot_err[2] << "\t"
            << CorrFD_tot_val[3] << "\t"    // inc neutral
            << CorrFD_tot_err[3] << "\n";
    // Close the output file
    outfile.close();
    Printf("Results printed to %s.", OutputFile.Data());   

    return;
}

void CalculateFD_Total_AOD(){

    // Load AxE value for kIncohJpsiToMu
    Bool_t AxE_IncJpsi = AxE_IncJpsi_Total_AOD();
    if(!AxE_IncJpsi){
        Printf("Values of AxE for kIncohJpsiToMu not found!");
        Printf("Create them using AccAndEffMC.c.");
        Printf("Terminating...");
        return;
    }
    // Calculate all FD corrections
    for(Int_t iFD = 1; iFD <= 4; iFD++){
        CalculateAxE_Psi2s_Total_AOD(iFD);
        // Cross-check: print the values that will be used for calculations
        Printf("AxE J/psi:\t(%.3f pm %.3f)%%", AxE_Jpsi_tot_val, AxE_Jpsi_tot_err);
        Printf("AxE Psi2s:\t(%.3f pm %.3f)%%", AxE_Psi2s_tot_val, AxE_Psi2s_tot_err);

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
        CorrFD_tot_val[iFD-1] = SigSL_p / SigSL_j_inc * AxE_Psi2s_tot_val / AxE_Jpsi_tot_val * BR * 100; // in percent already
        CorrFD_tot_err[iFD-1]= CorrFD_tot_val[iFD-1] * TMath::Sqrt(TMath::Power((AxE_Psi2s_tot_err/AxE_Psi2s_tot_val),2) + TMath::Power((AxE_Jpsi_tot_err/AxE_Jpsi_tot_val),2) + TMath::Power((BR_err/BR),2));
        Printf("FD corr: \t(%.3f pm %.3f)%%", CorrFD_tot_val[iFD-1], CorrFD_tot_err[iFD-1]);
        Printf(" ");
    }
    // Define an output text file
    TString OutputFile(Form("Results/FeedDown/AOD/FeedDownCorrections.txt"));
    ofstream outfile(OutputFile.Data());
    outfile << std::fixed << std::setprecision(3);
    // Print the results to the output text file
    outfile << CorrFD_tot_val[0] << "\t"    // coh charged
            << CorrFD_tot_err[0] << "\t"
            << CorrFD_tot_val[1] << "\t"    // inc charged
            << CorrFD_tot_err[1] << "\t"
            << CorrFD_tot_val[2] << "\t"    // coh neutral
            << CorrFD_tot_err[2] << "\t"
            << CorrFD_tot_val[3] << "\t"    // inc neutral
            << CorrFD_tot_err[3] << "\n";
    // Close the output file
    outfile.close();
    Printf("Results printed to %s.", OutputFile.Data());   

    return;
}

void CalculateFD_PtBins(){

    SetPtBinning();

    // Load AxE values for kIncohJpsiToMu
    Bool_t AxE_IncJpsi = AxE_IncJpsi_PtBins();
    if(!AxE_IncJpsi){
        Printf("Values of AxE for kIncohJpsiToMu not found!");
        Printf("Create them using AccAndEffMC.c.");
        Printf("Terminating...");
        return;
    }
    // Calculate the FD correction in all bins
    for(Int_t iFD = 1; iFD <= 4; iFD++){
        // Calculate AxE for Psi2s files
        CalculateAxE_Psi2s_PtBins(iFD);
        // Go over pt bins
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            // Cross-check: print the values that will be used for calculations
            Printf("Bin %i:", iBin+1);
            Printf("AxE J/psi:\t(%.3f pm %.3f)%%", AxE_IncJpsi_val[iBin], AxE_IncJpsi_err[iBin]);
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
            CorrFD_val[iFD-1][iBin] = SigSL_p / SigSL_j_inc * AxE_Psi2s_val[iBin] / AxE_IncJpsi_val[iBin] * BR * 100; // in percent already
            CorrFD_err[iFD-1][iBin] = CorrFD_val[iFD-1][iBin] * TMath::Sqrt(TMath::Power((AxE_Psi2s_err[iBin]/AxE_Psi2s_val[iBin]),2) + TMath::Power((AxE_IncJpsi_err[iBin]/AxE_IncJpsi_val[iBin]),2) + TMath::Power((BR_err/BR),2));
            Printf("FD corr: \t(%.3f pm %.3f)%%", CorrFD_val[iFD-1][iBin], CorrFD_err[iFD-1][iBin]);
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
                << CorrFD_val[0][iBin] << "\t"      // coh charged
                << CorrFD_err[0][iBin] << "\t"
                << CorrFD_val[1][iBin] << "\t"      // inc charged
                << CorrFD_err[1][iBin] << "\t"
                << CorrFD_val[2][iBin] << "\t"      // coh neutral
                << CorrFD_err[2][iBin] << "\t"
                << CorrFD_val[3][iBin] << "\t"      // inc neutral
                << CorrFD_err[3][iBin] << "\n";
    }
    // Close the output file
    outfile.close();
    Printf("Results printed to %s.", OutputFile.Data());   

    return;
}

void CalculateFD_PtFit(){

    for(Int_t iMC = 0; iMC < 2; iMC++){
        // Load AxE value for kCohJpsiToMu (iMC == 0) and kIncohJpsiToMu (iMC == 1)
        Bool_t AxE_Jpsi = AxE_IncJpsi_Total(iMC, 1, 2);
        if(!AxE_Jpsi){
            if(iMC == 0) Printf("Values of AxE for kCohJpsiToMu for PtFit not found!");
            if(iMC == 1) Printf("Values of AxE for kIncohJpsiToMu for PtFit not found!");
            Printf("Create them using AccAndEffMC.c.");
            Printf("Terminating...");
            return;
        }

        for(Int_t i = iMC+1; i <= 4; i+=2){
            CalculateAxE_Psi2s_PtFit(i);
            // Cross-check: print the values that will be used for calculations
            Printf("AxE J/psi:\t(%.3f pm %.3f)%%", AxE_Jpsi_tot_val, AxE_Jpsi_tot_err);
            Printf("AxE Psi2s:\t(%.3f pm %.3f)%%", AxE_Psi2s_tot_val, AxE_Psi2s_tot_err);
            // Load BR and cross sections
            if(i == 1 || i == 2){
                // charged processes
                BR = BR_ch;
                BR_err = BR_ch_err;
            }
            if(i == 3 || i == 4){
                // neutral processes
                BR = BR_ne;
                BR_err = BR_ne_err;
            }
            if(i == 1 || i == 3){
                // coherent Psi2s
                SigSL_p = SigSL_p_coh;
                SigSL_j = SigSL_j_coh;
            }
            if(i == 2 || i == 4){
                // incoherent Psi2s
                SigSL_p = SigSL_p_inc;
                SigSL_j = SigSL_j_inc;
            }
            // Calculate FD corrections         
            CorrFD_tot_val[i-1] = SigSL_p / SigSL_j * AxE_Psi2s_tot_val / AxE_Jpsi_tot_val * BR * 100; // in percent already
            CorrFD_tot_err[i-1]= CorrFD_tot_val[i-1] * TMath::Sqrt(TMath::Power((AxE_Psi2s_tot_err/AxE_Psi2s_tot_val),2) + TMath::Power((AxE_Jpsi_tot_err/AxE_Jpsi_tot_val),2) + TMath::Power((BR_err/BR),2));
            Printf("FD corr: \t(%.3f pm %.3f)%%", CorrFD_tot_val[i-1], CorrFD_tot_err[i-1]);
            Printf(" ");
        }    

    }
    // Define an output text file for fD_coh
    TString OutputFileCoh(Form("Results/FeedDown/ForPtFit/FDCorr_PtFit_Coh.txt"));
    ofstream outfile_coh(OutputFileCoh.Data());
    outfile_coh << std::fixed << std::setprecision(3);
    // Print the results to the output text file
    outfile_coh << CorrFD_tot_val[0] << "\t"    // coh charged
                << CorrFD_tot_err[0] << "\t"
                << CorrFD_tot_val[2] << "\t"    // coh neutral
                << CorrFD_tot_err[2] << "\n";
    // Close the output file
    outfile_coh.close();
    Printf("Results printed to %s.", OutputFileCoh.Data());   
    // Define an output text file for fD_inc
    TString OutputFileInc(Form("Results/FeedDown/ForPtFit/FDCorr_PtFit_Inc.txt"));
    ofstream outfile_inc(OutputFileInc.Data());
    outfile_inc << std::fixed << std::setprecision(3);
    // Print the results to the output text file
    outfile_inc << CorrFD_tot_val[1] << "\t"    // inc charged
                << CorrFD_tot_err[1] << "\t"
                << CorrFD_tot_val[3] << "\t"    // inc neutral
                << CorrFD_tot_err[3] << "\n";
    outfile_inc.close();
    Printf("Results printed to %s.", OutputFileInc.Data()); 

    return;
}

Bool_t AxE_IncJpsi_Total(Int_t iMC, Int_t iMassCut, Int_t iPtCut){

    TString FilePath;
    if(iMC == 0) FilePath = Form("Results/AccAndEffMC/AxE_tot_CohJ_MassCut%i_PtCut%i.txt", iMassCut, iPtCut);
    if(iMC == 1) FilePath = Form("Results/AccAndEffMC/AxE_tot_IncJ_MassCut%i_PtCut%i.txt", iMassCut, iPtCut);
    // Check if already calculated
    ifstream FileIn;
    FileIn.open(FilePath.Data());
    if(!(FileIn.fail())){
        // Read data from the file
        while(!FileIn.eof()){
            FileIn >> AxE_Jpsi_tot_val >> AxE_Jpsi_tot_err;
        }
        FileIn.close(); 
        Printf("AxE for kIncohJpsiToMu with pt cut no. %i has already been calculated.\n", iPtCut); 

        return kTRUE;
    } else {
        // AxE for kIncohJpsiToMu not found
        return kFALSE;
    }
}

Bool_t AxE_IncJpsi_Total_AOD(){

    TString FilePath = Form("Results/AccAndEffMC/AxE_AOD_tot_MassCut1_PtCut0.txt");
    // Check if already calculated
    ifstream FileIn;
    FileIn.open(FilePath.Data());
    if(!(FileIn.fail())){
        // Read data from the file
        while(!FileIn.eof()){
            FileIn >> AxE_Jpsi_tot_val >> AxE_Jpsi_tot_err;
        }
        FileIn.close(); 
        Printf("AxE for kIncohJpsiToMu has already been calculated.\n"); 

        return kTRUE;
    } else {
        // AxE for kIncohJpsiToMu not found
        return kFALSE;
    }
}

Bool_t AxE_IncJpsi_PtBins(){

    TString FilePath = Form("Results/AccAndEffMC/AxE_%ibins.txt", nPtBins);
    // Check if already calculated
    ifstream FileIn;
    FileIn.open(FilePath.Data());
    if(!(FileIn.fail())){
        // Read data from the file
        Int_t i_line = 0;
        while(!FileIn.eof()){
            FileIn >> i_bin >> AxE_IncJpsi_val[i_line] >> AxE_IncJpsi_err[i_line];
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

void CalculateAxE_Psi2s_Total(Int_t iPtCut, Int_t iFD){

    TString FilePath = Form("Results/FeedDown/AxE_%s_tot_PtCut%i", (Datasets[iFD-1]).Data(), iPtCut);
    // Check if already calculated
    ifstream FileIn;
    FileIn.open((FilePath + ".txt").Data());
    if(!(FileIn.fail())){
        // Read data from the file
        while(!FileIn.eof()){
            FileIn >> AxE_Psi2s_tot_val >> AxE_Psi2s_tot_err;
        }
        FileIn.close(); 
        Printf("############################");
        Printf("AxE for %s with pt cut no. %i has already been calculated.\n", (Datasets[iFD-1]).Data(), iPtCut);

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

    Double_t NRec = 0;
    // Loop over tree entries
    Int_t nEntriesAnalysed = 0;
    for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++){
        tRec->GetEntry(iEntry);
        if(EventPassedMCRec(0, iPtCut)) NRec++;
        if((iEntry+1) % 200000 == 0){
            nEntriesAnalysed += 200000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }
    Printf("Done.");
    
    TTree *tGen = dynamic_cast<TTree*> (fRec->Get("AnalysisOutput/fTreeJPsiMCGen"));
    if(tGen) Printf("MC gen tree loaded.");
    
    ConnectTreeVariablesMCGen(tGen);

    Printf("Tree %s has %lli entries.", tGen->GetName(), tGen->GetEntries());

    Double_t NGen = 0;
    // Loop over tree entries
    nEntriesAnalysed = 0;
    for(Int_t iEntry = 0; iEntry < tGen->GetEntries(); iEntry++){
        tGen->GetEntry(iEntry);
        if(EventPassedMCGen(iPtCut)) NGen++;
        if((iEntry+1) % 500000 == 0){
            nEntriesAnalysed += 500000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }
    Printf("Done.");

    // Calculate AxE per bin
    AxE_Psi2s_tot_val = NRec / NGen * 100;
    AxE_Psi2s_tot_err = CalculateErrorBayes(NRec, NGen) * 100;
    Printf("AxE = (%.4f pm %.4f)%%", AxE_Psi2s_tot_val, AxE_Psi2s_tot_err);

    // Save the results to text file
    ofstream outfile((FilePath + ".txt").Data());
    outfile << std::fixed << std::setprecision(3);
    //outfile << Form("AxE_Psi2s [%%]\tAxE_Psi2s_err [%%] \n");
    outfile << AxE_Psi2s_tot_val << "\t" << AxE_Psi2s_tot_err << "\n";
    outfile.close();
    Printf("Results printed to %s.", (FilePath + ".txt").Data());    

    return;
}

void CalculateAxE_Psi2s_Total_AOD(Int_t iFD){
    TString FilePath = Form("Results/FeedDown/AOD/AxE_%s", (Datasets[iFD-1]).Data());
    // Check if already calculated
    ifstream FileIn;
    FileIn.open((FilePath + ".txt").Data());
    if(!(FileIn.fail())){
        // Read data from the file
        while(!FileIn.eof()){
            FileIn >> AxE_Psi2s_tot_val >> AxE_Psi2s_tot_err;
        }
        FileIn.close(); 
        Printf("############################");
        Printf("AxE for %s has already been calculated.\n", (Datasets[iFD-1]).Data());

        return;
    }

    TFile *fRec = NULL;
    switch(iFD){
        case 1:
            fRec = TFile::Open("Trees/AnalysisDataAOD/MC_rec_18qr_kCohPsi2sToMuPi.root", "read");
            break;
        case 2:
            fRec = TFile::Open("Trees/AnalysisDataAOD/MC_rec_18qr_kIncohPsi2sToMuPi.root", "read");
            break;
        case 3:
            fRec = TFile::Open("Trees/AnalysisDataAOD/MC_rec_18qr_kCohPsi2sToMuPi_NeutPi.root", "read");
            break;
        case 4:
            fRec = TFile::Open("Trees/AnalysisDataAOD/MC_rec_18qr_kIncohPsi2sToMuPi_NeutPi.root", "read");
            break;
    }
    if(fRec) Printf("MC rec file loaded.");

    TTree *tRec = dynamic_cast<TTree*> (fRec->Get("analysisTree"));
    if(tRec) Printf("MC rec tree loaded.");
    
    ConnectTreeVariablesMCRec_AOD(tRec);

    Printf("Now calculating FD for %s.", fRec->GetName());

    Printf("Tree %s has %lli entries.", tRec->GetName(), tRec->GetEntries());

    Double_t NRec = 0;
    // Loop over tree entries
    Int_t nEntriesAnalysed = 0;
    for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++){
        tRec->GetEntry(iEntry);
        if(EventPassedMCRec_AOD(1, 0)) NRec++;
        if((iEntry+1) % 200000 == 0){
            nEntriesAnalysed += 200000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }
    Printf("Done.");

    TFile *fGen = NULL;
    switch(iFD){
        case 1:
            fGen = TFile::Open("Trees/AnalysisDataAOD/MC_gen_18qr_kCohPsi2sToMuPi.root", "read");
            break;
        case 2:
            fGen = TFile::Open("Trees/AnalysisDataAOD/MC_gen_18qr_kIncohPsi2sToMuPi.root", "read");
            break;
        case 3:
            fGen = TFile::Open("Trees/AnalysisDataAOD/MC_gen_18qr_kCohPsi2sToMuPi_NeutPi.root", "read");
            break;
        case 4:
            fGen = TFile::Open("Trees/AnalysisDataAOD/MC_gen_18qr_kIncohPsi2sToMuPi_NeutPi.root", "read");
            break;
    }
    if(fGen) Printf("MC gen file loaded.");

    TTree *tGen = dynamic_cast<TTree*> (fGen->Get("MCgenTree"));
    if(tGen) Printf("MC gen tree loaded.");
    
    ConnectTreeVariablesMCGen_AOD_old(tGen);

    Printf("Tree %s has %lli entries.", tGen->GetName(), tGen->GetEntries());

    Double_t NGen = 0;
    // Loop over tree entries
    nEntriesAnalysed = 0;
    for(Int_t iEntry = 0; iEntry < tGen->GetEntries(); iEntry++){
        tGen->GetEntry(iEntry);
        if(EventPassedMCGen(0)) NGen++;
        if((iEntry+1) % 500000 == 0){
            nEntriesAnalysed += 500000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }
    Printf("Done.");

    // Calculate AxE per bin
    Printf("N_rec = %.0f", NRec);
    Printf("N_gen = %.0f", NGen);
    AxE_Psi2s_tot_val = NRec / NGen * 100;
    AxE_Psi2s_tot_err = CalculateErrorBayes(NRec, NGen) * 100;
    Printf("AxE = (%.4f pm %.4f)%%", AxE_Psi2s_tot_val, AxE_Psi2s_tot_err);

    // Save the results to text file
    ofstream outfile((FilePath + ".txt").Data());
    outfile << std::fixed << std::setprecision(3);
    //outfile << Form("AxE_Psi2s [%%]\tAxE_Psi2s_err [%%] \n");
    outfile << AxE_Psi2s_tot_val << "\t" << AxE_Psi2s_tot_err << "\n";
    outfile.close();
    Printf("Results printed to %s.", (FilePath + ".txt").Data());    

    return;    
}

void CalculateAxE_Psi2s_PtBins(Int_t iFD){

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
        AxE_Psi2s_val[iPtBin] = NRec[iPtBin] / NGen[iPtBin] * 100;
        AxE_Psi2s_err[iPtBin] = CalculateErrorBayes(NRec[iPtBin], NGen[iPtBin]) * 100;
        Printf("AxE = (%.4f pm %.4f)%%", AxE_Psi2s_val[iPtBin], AxE_Psi2s_err[iPtBin]);
    }

    // Save the results to text file
    ofstream outfile((FilePath + ".txt").Data());
    outfile << std::fixed << std::setprecision(3);
    //outfile << Form("Bin \tPtLow \tPtUpp \tAxE_Psi2s [%%]\tAxE_Psi2s_err [%%] \n");
    for(Int_t i = 1; i <= nPtBins; i++){
        outfile << i << "\t" << ptBoundaries[i-1] << "\t" << ptBoundaries[i] << "\t" << AxE_Psi2s_val[i-1] << "\t" << AxE_Psi2s_err[i-1] << "\n";
    }
    outfile.close();
    Printf("Results printed to %s.", (FilePath + ".txt").Data());    

    return;
}

void CalculateAxE_Psi2s_PtFit(Int_t iFD){

    TString FilePath = Form("Results/FeedDown/ForPtFit/AxE_%s_PtFit", (Datasets[iFD-1]).Data());
    // Check if already calculated
    ifstream FileIn;
    FileIn.open((FilePath + ".txt").Data());
    if(!(FileIn.fail())){
        // Read data from the file
        while(!FileIn.eof()){
            FileIn >> AxE_Psi2s_tot_val >> AxE_Psi2s_tot_err;
        }
        FileIn.close(); 
        Printf("############################");
        Printf("AxE for %s_PtFit has already been calculated.\n", (Datasets[iFD-1]).Data());

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

    Double_t NRec = 0;
    // Loop over tree entries
    Int_t nEntriesAnalysed = 0;
    for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++){
        tRec->GetEntry(iEntry);
        // Here with cuts:
        // 3.0 < m < 3.2 GeV
        // 0.0 < pt < 2.0 GeV
        if(EventPassedMCRec(1, 2)) NRec++;
        if((iEntry+1) % 200000 == 0){
            nEntriesAnalysed += 200000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }
    Printf("Done.");
    
    TTree *tGen = dynamic_cast<TTree*> (fRec->Get("AnalysisOutput/fTreeJPsiMCGen"));
    if(tGen) Printf("MC gen tree loaded.");
    
    ConnectTreeVariablesMCGen(tGen);

    Printf("Tree %s has %lli entries.", tGen->GetName(), tGen->GetEntries());

    Double_t NGen = 0;
    // Loop over tree entries
    nEntriesAnalysed = 0;
    for(Int_t iEntry = 0; iEntry < tGen->GetEntries(); iEntry++){
        tGen->GetEntry(iEntry);
        if(EventPassedMCGen(-1)) NGen++;
        if((iEntry+1) % 500000 == 0){
            nEntriesAnalysed += 500000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }
    Printf("Done.");

    // Calculate AxE per bin
    AxE_Psi2s_tot_val = NRec / NGen * 100;
    AxE_Psi2s_tot_err = CalculateErrorBayes(NRec, NGen) * 100;
    Printf("AxE = (%.4f pm %.4f)%%", AxE_Psi2s_tot_val, AxE_Psi2s_tot_err);

    // Save the results to text file
    ofstream outfile((FilePath + ".txt").Data());
    outfile << std::fixed << std::setprecision(3);
    //outfile << Form("AxE_Psi2s [%%]\tAxE_Psi2s_err [%%] \n");
    outfile << AxE_Psi2s_tot_val << "\t" << AxE_Psi2s_tot_err << "\n";
    outfile.close();
    Printf("Results printed to %s.", (FilePath + ".txt").Data());    

    return;
}

Double_t CalculateErrorBayes(Double_t k, Double_t n){ // k = NRec, n = NGen

    Double_t var = (k + 1) * (k + 2) / (n + 2) / (n + 3) - (k + 1) * (k + 1) / (n + 2) / (n + 2);
    Double_t err = TMath::Sqrt(var);

    return err;
}