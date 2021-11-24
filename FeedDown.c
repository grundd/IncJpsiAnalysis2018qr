// FeedDown.c
// David Grund, Sep 24, 2021

// cpp headers
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
#include <string> // getline
// root headers
#include "TFile.h"
#include "TMath.h"
// my headers
#include "AnalysisManager.h"

TString path = "Results/FeedDown/";

void CalculateFD_Total(Bool_t bAOD = kFALSE);
void CalculateFD_PtBins();
void CalculateFD_PtFit(Bool_t bAOD = kFALSE, Bool_t bRatioMeasured = kFALSE);

// Temporary variables when loading the data:
Int_t iBin;
// Arrays of AxEs
Double_t AxE_CohJ_val = 0;
Double_t AxE_CohJ_err = 0;
Double_t AxE_IncJ_val[nPtBins] = { 0 };
Double_t AxE_IncJ_err[nPtBins] = { 0 };
Double_t AxE_Psi2s_val[4][nPtBins] = { 0 };
Double_t AxE_Psi2s_err[4][nPtBins] = { 0 };
// Results
Double_t fD_coh_ch_val[nPtBins] = { 0 };
Double_t fD_inc_ch_val[nPtBins] = { 0 };
Double_t fD_coh_ne_val[nPtBins] = { 0 };
Double_t fD_inc_ne_val[nPtBins] = { 0 };
Double_t fD_coh_ch_err[nPtBins] = { 0 };
Double_t fD_inc_ch_err[nPtBins] = { 0 };
Double_t fD_coh_ne_err[nPtBins] = { 0 };
Double_t fD_inc_ne_err[nPtBins] = { 0 };
// STARlight cross sections
Double_t sig_SL_j_inc = 5.247; //mb
Double_t sig_SL_j_coh = 12.504;//mb
Double_t sig_SL_p_inc = 0.92;  //mb
Double_t sig_SL_p_coh = 2.52;  //mb
Double_t ratio_coh = 0.18; // from Michal's paper
// Weighing of SL cross sections by NGen(bin)/NGen(tot)
// order: JInc 	PCohCh 	PIncCh 	PCohNe 	PIncNe 
Double_t NGen_tot[5] = { 0 };
Double_t NGen_bin[5][nPtBins] = { 0 };
// Branching ratios
Double_t BR_ch = 0.3468;
Double_t BR_ch_err = 0.0030;
Double_t BR_ne = 0.2538;
Double_t BR_ne_err = 0.0032;

void FeedDown(){

    // ESDs
    CalculateFD_Total(kFALSE);

    // AODs
    CalculateFD_Total(kTRUE);

    // Pt bins
    CalculateFD_PtBins();

    // For PtFit, with ratio from SL
    // ESDs
    CalculateFD_PtFit(kFALSE, kFALSE);
    // AODs
    CalculateFD_PtFit(kTRUE, kFALSE);

    // For PtFit, with measured ratio
    // ESDs
    CalculateFD_PtFit(kFALSE, kTRUE);
    // AODs
    CalculateFD_PtFit(kTRUE, kTRUE);

    return;
}


void CalculateFD_Total(Bool_t bAOD){

    // 1) Load values of AxE for all datasets
    TString str;
    if(bAOD == kFALSE) str = "Results/AccAndEffMC/AxE_FeedDown_Total.txt";
    else str = "Results/AccAndEffMC/AOD/AxE_AOD_FeedDown_Total.txt";
    ifstream in_file(str.Data());

    if(in_file.fail()){
        Printf("Missing input file...");
        return;
    } else {
        Int_t i = 0;
        char ch[8];
        std::string str;
        while(std::getline(in_file,str)){
            Printf("Reading line %i: %s", i, str.data());
            istringstream in_stream(str);
            // skip first line
            if(i == 1) in_stream >> ch >> AxE_IncJ_val[0] >> AxE_IncJ_err[0];
            if(i > 1)  in_stream >> ch >> AxE_Psi2s_val[i-2][0] >> AxE_Psi2s_err[i-2][0];
            i++;   
        }
    }
    in_file.close();
    Printf("Input file loaded...");

    // 2) Define output file
    if(bAOD == kFALSE) str = Form("%sFeedDown_Total.txt", path.Data());
    else str = Form("%sFeedDown_Total_AOD.txt", path.Data());
    ofstream outfile(str.Data());
    outfile << std::fixed << std::setprecision(4);
    outfile << Form("fD[%%] \tCohCh \tErr \tIncCh \tErr \tCohNe \tErr \tIncNe \tErr \n");

    // 3) Calculate fD correction
    for(Int_t i = 0; i < 1; i++){
        // all values of fD in percent already
        fD_coh_ch_val[i] = sig_SL_p_coh / sig_SL_j_inc * AxE_Psi2s_val[0][i] / AxE_IncJ_val[i] * BR_ch * 100;
        fD_inc_ch_val[i] = sig_SL_p_inc / sig_SL_j_inc * AxE_Psi2s_val[1][i] / AxE_IncJ_val[i] * BR_ch * 100;
        fD_coh_ne_val[i] = sig_SL_p_coh / sig_SL_j_inc * AxE_Psi2s_val[2][i] / AxE_IncJ_val[i] * BR_ne * 100;
        fD_inc_ne_val[i] = sig_SL_p_inc / sig_SL_j_inc * AxE_Psi2s_val[3][i] / AxE_IncJ_val[i] * BR_ne * 100;
        if(AxE_Psi2s_val[0][i] != 0){
            fD_coh_ch_err[i] = fD_coh_ch_val[i] * TMath::Sqrt(
                TMath::Power((AxE_Psi2s_err[0][i]/AxE_Psi2s_val[0][i]),2) + 
                TMath::Power((AxE_IncJ_err[i]/AxE_IncJ_val[i]),2) + 
                TMath::Power((BR_ch_err/BR_ch),2));
        } else fD_coh_ch_err[i] = 0;
        if(AxE_Psi2s_val[1][i] != 0){
            fD_inc_ch_err[i] = fD_inc_ch_val[i] * TMath::Sqrt(
                TMath::Power((AxE_Psi2s_err[1][i]/AxE_Psi2s_val[1][i]),2) + 
                TMath::Power((AxE_IncJ_err[i]/AxE_IncJ_val[i]),2) + 
                TMath::Power((BR_ch_err/BR_ch),2));
        } else fD_inc_ch_err[i] = 0;
        if(AxE_Psi2s_val[2][i] != 0){
            fD_coh_ne_err[i] = fD_coh_ne_val[i] * TMath::Sqrt(
                TMath::Power((AxE_Psi2s_err[2][i]/AxE_Psi2s_val[2][i]),2) + 
                TMath::Power((AxE_IncJ_err[i]/AxE_IncJ_val[i]),2) + 
                TMath::Power((BR_ne_err/BR_ne),2));
        } else fD_coh_ne_err[i] = 0;
        if(AxE_Psi2s_val[3][i] != 0){
            fD_inc_ne_err[i] = fD_inc_ne_val[i] * TMath::Sqrt(
                TMath::Power((AxE_Psi2s_err[3][i]/AxE_Psi2s_val[3][i]),2) + 
                TMath::Power((AxE_IncJ_err[i]/AxE_IncJ_val[i]),2) + 
                TMath::Power((BR_ne_err/BR_ne),2));
        } else fD_inc_ne_err[i] = 0; 
        // Print the results to text file
        outfile << "Total:";
        outfile << "\t" << fD_coh_ch_val[i] << "\t" << fD_coh_ch_err[i] 
                << "\t" << fD_inc_ch_val[i] << "\t" << fD_inc_ch_err[i] 
                << "\t" << fD_coh_ne_val[i] << "\t" << fD_coh_ne_err[i] 
                << "\t" << fD_inc_ne_val[i] << "\t" << fD_inc_ne_err[i] << "\n";
    }
    outfile.close();
    Printf("*** Results printed to %s.***", str.Data());

    return;
}

void CalculateFD_PtBins(){

    // 1) Load values of AxE for all datasets
    ifstream in_file(Form("Results/AccAndEffMC/AxE_FeedDown_%ibins.txt", nPtBins));

    if(in_file.fail()){
        Printf("Missing input file...");
        return;
    } else {
        Int_t i = 0;
        std::string str;
        while(std::getline(in_file,str)){
            Printf("Reading line %i: %s", i, str.data());
            istringstream in_stream(str);
            // skip first line
            if(i != 0) in_stream >> iBin >> AxE_IncJ_val[i-1] >> AxE_IncJ_err[i-1]
                                         >> AxE_Psi2s_val[0][i-1] >> AxE_Psi2s_err[0][i-1]  // PCohCh
                                         >> AxE_Psi2s_val[1][i-1] >> AxE_Psi2s_err[1][i-1]  // PIncCh
                                         >> AxE_Psi2s_val[2][i-1] >> AxE_Psi2s_err[2][i-1]  // PCohNe
                                         >> AxE_Psi2s_val[3][i-1] >> AxE_Psi2s_err[3][i-1]; // PIncNe
            i++;   
        }
    }
    in_file.close();
    Printf("Input file loaded...");

    // 2) Load NGen per bin and NGen tot to weigh SL cross sections
    ifstream in_file2(Form("Results/AccAndEffMC/AxE_FeedDown_%ibins_NGen_SL.txt", nPtBins));

    if(in_file2.fail()){
        Printf("Missing input file...");
        return;
    } else {
        Int_t i = 0;
        char ch[8];
        std::string str;
        while(std::getline(in_file2,str)){
            Printf("Reading line %i: %s", i, str.data());
            istringstream in_stream(str);
            // skip first line
            if(i == 1) in_stream >> ch >> NGen_tot[0] >> NGen_tot[1] >> NGen_tot[2] >> NGen_tot[3] >> NGen_tot[4];
            if(i > 1) in_stream >> iBin >> NGen_bin[0][i-2] >> NGen_bin[1][i-2] >> NGen_bin[2][i-2] >> NGen_bin[3][i-2] >> NGen_bin[4][i-2];
            i++;   
        }
    }
    in_file2.close();
    Printf("Input file loaded...");
    // Cross-check:
    //Printf("%.0f", NGen_tot[0]); 
    //Printf("%.0f", NGen_bin[4][4]);

    // 3a) Define the output file for fD
    TString str = Form("%sFeedDown_%ibins.txt", path.Data(), nPtBins);
    ofstream outfile(str.Data());
    outfile << std::fixed << std::setprecision(4);
    outfile << Form("fD[%%] \tCohCh \tErr \tIncCh \tErr \tCohNe \tErr \tIncNe \tErr \n");
    // 3b) Define the output file for weighed STARlight cross sections
    TString str2 = Form("%sSTARlight_WeighedCrossSections_%ibins.txt", path.Data(), nPtBins);
    ofstream outfile2(str2.Data());
    outfile2 << Form("Bin \tIncNGen\tCChNgen\tIChNGen\tCNeNgen\tINeNgen\tIncSig \tCChSig \tIChSig \tCNeSig \tINeSig \n"); 

    // 4) Loop over bins, calculate fD corrections
    for(Int_t i = 0; i < nPtBins; i++){
        // all values of fD in percent already
        fD_coh_ch_val[i] = sig_SL_p_coh * (NGen_bin[1][i]/NGen_tot[1]) / (sig_SL_j_inc * (NGen_bin[0][i]/NGen_tot[0])) * AxE_Psi2s_val[0][i] / AxE_IncJ_val[i] * BR_ch * 100;
        fD_inc_ch_val[i] = sig_SL_p_inc * (NGen_bin[2][i]/NGen_tot[2]) / (sig_SL_j_inc * (NGen_bin[0][i]/NGen_tot[0])) * AxE_Psi2s_val[1][i] / AxE_IncJ_val[i] * BR_ch * 100;
        fD_coh_ne_val[i] = sig_SL_p_coh * (NGen_bin[3][i]/NGen_tot[3]) / (sig_SL_j_inc * (NGen_bin[0][i]/NGen_tot[0])) * AxE_Psi2s_val[2][i] / AxE_IncJ_val[i] * BR_ne * 100;
        fD_inc_ne_val[i] = sig_SL_p_inc * (NGen_bin[4][i]/NGen_tot[4]) / (sig_SL_j_inc * (NGen_bin[0][i]/NGen_tot[0])) * AxE_Psi2s_val[3][i] / AxE_IncJ_val[i] * BR_ne * 100;
        if(AxE_Psi2s_val[0][i] != 0){
            fD_coh_ch_err[i] = fD_coh_ch_val[i] * TMath::Sqrt(
                TMath::Power((AxE_Psi2s_err[0][i]/AxE_Psi2s_val[0][i]),2) + 
                TMath::Power((AxE_IncJ_err[i]/AxE_IncJ_val[i]),2) + 
                TMath::Power((BR_ch_err/BR_ch),2));
        } else fD_coh_ch_err[i] = 0;
        if(AxE_Psi2s_val[1][i] != 0){
            fD_inc_ch_err[i] = fD_inc_ch_val[i] * TMath::Sqrt(
                TMath::Power((AxE_Psi2s_err[1][i]/AxE_Psi2s_val[1][i]),2) + 
                TMath::Power((AxE_IncJ_err[i]/AxE_IncJ_val[i]),2) + 
                TMath::Power((BR_ch_err/BR_ch),2));
        } else fD_inc_ch_err[i] = 0;
        if(AxE_Psi2s_val[2][i] != 0){
            fD_coh_ne_err[i] = fD_coh_ne_val[i] * TMath::Sqrt(
                TMath::Power((AxE_Psi2s_err[2][i]/AxE_Psi2s_val[2][i]),2) + 
                TMath::Power((AxE_IncJ_err[i]/AxE_IncJ_val[i]),2) + 
                TMath::Power((BR_ne_err/BR_ne),2));
        } else fD_coh_ne_err[i] = 0;
        if(AxE_Psi2s_val[3][i] != 0){
            fD_inc_ne_err[i] = fD_inc_ne_val[i] * TMath::Sqrt(
                TMath::Power((AxE_Psi2s_err[3][i]/AxE_Psi2s_val[3][i]),2) + 
                TMath::Power((AxE_IncJ_err[i]/AxE_IncJ_val[i]),2) + 
                TMath::Power((BR_ne_err/BR_ne),2));
        } else fD_inc_ne_err[i] = 0; 
        // Print the results to text file
        outfile << i + 1;
        outfile << "\t" << fD_coh_ch_val[i] << "\t" << fD_coh_ch_err[i] 
                << "\t" << fD_inc_ch_val[i] << "\t" << fD_inc_ch_err[i] 
                << "\t" << fD_coh_ne_val[i] << "\t" << fD_coh_ne_err[i] 
                << "\t" << fD_inc_ne_val[i] << "\t" << fD_inc_ne_err[i] << "\n";  
        outfile2 << i + 1 << std::fixed << std::setprecision(0);
        outfile2 << "\t" << NGen_bin[0][i] << "\t" << NGen_bin[1][i] << "\t" << NGen_bin[2][i] << "\t" << NGen_bin[3][i] << "\t" << NGen_bin[4][i];
        outfile2 << std::fixed << std::setprecision(3);
        outfile2 << "\t" << sig_SL_j_inc * (NGen_bin[0][i]/NGen_tot[0]) 
                 << "\t" << sig_SL_p_coh * (NGen_bin[1][i]/NGen_tot[1]) 
                 << "\t" << sig_SL_p_inc * (NGen_bin[2][i]/NGen_tot[2])
                 << "\t" << sig_SL_p_coh * (NGen_bin[3][i]/NGen_tot[3]) 
                 << "\t" << sig_SL_p_inc * (NGen_bin[4][i]/NGen_tot[4]) << "\n";
    }
    outfile.close();
    Printf("*** Results printed to %s.***", str.Data());
    outfile2 << std::fixed << std::setprecision(0);
    outfile2 << "Total" << "\t" << NGen_tot[0] << "\t" << NGen_tot[1] << "\t" << NGen_tot[2] << "\t" << NGen_tot[3] << "\t" << NGen_tot[4] << "\n";
    outfile2.close();
    Printf("*** Results printed to %s.***", str2.Data());
    
    return;
}

void CalculateFD_PtFit(Bool_t bAOD, Bool_t bRatioMeasured){

    // 1) Load values of AxE for all datasets
    TString str;
    if(bAOD == kFALSE) str = "Results/AccAndEffMC/AxE_PtFit.txt";
    else str = "Results/AccAndEffMC/AOD/AxE_AOD_PtFit.txt";
    ifstream in_file(str.Data());

    if(in_file.fail()){
        Printf("Missing input file...");
        return;
    } else {
        Int_t i = 0;
        char ch[8];
        std::string str;
        while(std::getline(in_file,str)){
            Printf("Reading line %i: %s", i, str.data());
            istringstream in_stream(str);
            // skip first line
            if(i == 1) in_stream >> ch >> AxE_CohJ_val >> AxE_CohJ_err;
            if(i == 2) in_stream >> ch >> AxE_IncJ_val[0] >> AxE_IncJ_err[0];
            if(i > 1)  in_stream >> ch >> AxE_Psi2s_val[i-3][0] >> AxE_Psi2s_err[i-3][0];
            i++;   
        }
    }
    in_file.close();
    Printf("Input file loaded...");

    // 2) Define output file
    str = Form("%sFeedDown_PtFit", path.Data());
    if(bAOD) str.Append("_AOD");
    if(bRatioMeasured) str.Append("_ratMeas");
    str.Append(".txt");
    ofstream outfile(str.Data());
    outfile << std::fixed << std::setprecision(4);
    outfile << Form("fD[%%] \tCohCh \tErr \tIncCh \tErr \tCohNe \tErr \tIncNe \tErr \n");

    // 3) Calculate fD correction
    for(Int_t i = 0; i < 1; i++){
        // all values of fD in percent already
        if(!bRatioMeasured) fD_coh_ch_val[i] = sig_SL_p_coh / sig_SL_j_coh * AxE_Psi2s_val[0][i] / AxE_CohJ_val * BR_ch * 100;
        else fD_coh_ch_val[i] = ratio_coh * AxE_Psi2s_val[0][i] / AxE_CohJ_val * BR_ch * 100;
        fD_inc_ch_val[i] = sig_SL_p_inc / sig_SL_j_inc * AxE_Psi2s_val[1][i] / AxE_IncJ_val[i] * BR_ch * 100;
        if(!bRatioMeasured) fD_coh_ne_val[i] = sig_SL_p_coh / sig_SL_j_coh * AxE_Psi2s_val[2][i] / AxE_CohJ_val * BR_ne * 100;
        else fD_coh_ne_val[i] = ratio_coh * AxE_Psi2s_val[2][i] / AxE_CohJ_val * BR_ne * 100;
        fD_inc_ne_val[i] = sig_SL_p_inc / sig_SL_j_inc * AxE_Psi2s_val[3][i] / AxE_IncJ_val[i] * BR_ne * 100;
        if(AxE_Psi2s_val[0][i] != 0){
            fD_coh_ch_err[i] = fD_coh_ch_val[i] * TMath::Sqrt(
                TMath::Power((AxE_Psi2s_err[0][i]/AxE_Psi2s_val[0][i]),2) + 
                TMath::Power((AxE_CohJ_err/AxE_CohJ_val),2) + 
                TMath::Power((BR_ch_err/BR_ch),2));
        } else fD_coh_ch_err[i] = 0;
        if(AxE_Psi2s_val[1][i] != 0){
            fD_inc_ch_err[i] = fD_inc_ch_val[i] * TMath::Sqrt(
                TMath::Power((AxE_Psi2s_err[1][i]/AxE_Psi2s_val[1][i]),2) + 
                TMath::Power((AxE_IncJ_err[i]/AxE_IncJ_val[i]),2) + 
                TMath::Power((BR_ch_err/BR_ch),2));
        } else fD_inc_ch_err[i] = 0;
        if(AxE_Psi2s_val[2][i] != 0){
            fD_coh_ne_err[i] = fD_coh_ne_val[i] * TMath::Sqrt(
                TMath::Power((AxE_Psi2s_err[2][i]/AxE_Psi2s_val[2][i]),2) + 
                TMath::Power((AxE_CohJ_err/AxE_CohJ_val),2) + 
                TMath::Power((BR_ne_err/BR_ne),2));
        } else fD_coh_ne_err[i] = 0;
        if(AxE_Psi2s_val[3][i] != 0){
            fD_inc_ne_err[i] = fD_inc_ne_val[i] * TMath::Sqrt(
                TMath::Power((AxE_Psi2s_err[3][i]/AxE_Psi2s_val[3][i]),2) + 
                TMath::Power((AxE_IncJ_err[i]/AxE_IncJ_val[i]),2) + 
                TMath::Power((BR_ne_err/BR_ne),2));
        } else fD_inc_ne_err[i] = 0; 
        // Print the results to text file
        outfile << "Total:";
        outfile << "\t" << fD_coh_ch_val[i] << "\t" << fD_coh_ch_err[i] 
                << "\t" << fD_inc_ch_val[i] << "\t" << fD_inc_ch_err[i] 
                << "\t" << fD_coh_ne_val[i] << "\t" << fD_coh_ne_err[i] 
                << "\t" << fD_inc_ne_val[i] << "\t" << fD_inc_ne_err[i] << "\n";
    }
    outfile.close();
    Printf("*** Results printed to %s.***", str.Data());

    return;
}