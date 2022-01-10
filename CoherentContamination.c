// CoherentContamination.c
// David Grund, Sep 30, 2021

// cpp headers
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
#include <string> // getline
// root headers
#include "TString.h"
#include "TMath.h"
// my headers
#include "AnalysisManager.h"

TString path = "Results/CoherentContamination/";

Double_t sig_SL_j_inc = 5.247; //mb
Double_t sig_SL_j_coh = 12.504;//mb

Double_t AxE_CohJ_val[nPtBins];
Double_t AxE_CohJ_err[nPtBins];
Double_t AxE_IncJ_val[nPtBins];
Double_t AxE_IncJ_err[nPtBins];
Double_t NGen_tot[2] = { 0 };
Double_t NGen_bin[2][nPtBins] = { 0 };
Double_t fC_val[nPtBins];
Double_t fC_err[nPtBins];

void CalculateFC_Total();
void CalculateFC_PtBins();

void CoherentContamination(){

    //CalculateFC_Total();

    CalculateFC_PtBins();

    return;
}

void CalculateFC_Total(){

    SetPtBinning();

    // Load values of AxE
    TString str = "Results/AccAndEffMC/AxE_CorrFC_Total.txt";
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
            if(i == 1) in_stream >> ch >> AxE_CohJ_val[i-1] >> AxE_CohJ_err[i-1] >> AxE_IncJ_val[i-1] >> AxE_IncJ_err[i-1];
            i++;   
        }
    }
    in_file.close();
    Printf("Input file loaded...");

    // Calculate FC corrections and print them into text file
    str = Form("%sCoherentContamination_Total.txt", path.Data());
    ofstream outfile(str.Data());
    outfile << std::fixed << std::setprecision(4);
    outfile << Form("[%%] \tfC \tErr \n");  
    for(Int_t i = 0; i < 1; i++){  
        fC_val[i] = sig_SL_j_coh / sig_SL_j_inc * AxE_CohJ_val[i] / AxE_IncJ_val[i] * 100; // in percent already
        if(AxE_CohJ_val[i] != 0){
            fC_err[i] = fC_val[i] * TMath::Sqrt(
                TMath::Power((AxE_CohJ_err[i]/AxE_CohJ_val[i]),2) + 
                TMath::Power((AxE_IncJ_err[i]/AxE_IncJ_val[i]),2));
        } else fC_err[i] = 0;
        outfile << "Total:";
        outfile << "\t" << fC_val[i] << "\t" << fC_err[i];
    }
    outfile.close();
    Printf("*** Results printed to %s.***", str.Data());

    return;
}

void CalculateFC_PtBins(){

    // Load values of AxE
    TString str = Form("Results/AccAndEffMC/AxE_CorrFC_%ibins.txt", nPtBins);
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
            if(i > 0) in_stream >> ch >> AxE_CohJ_val[i-1] >> AxE_CohJ_err[i-1] >> AxE_IncJ_val[i-1] >> AxE_IncJ_err[i-1];
            i++;   
        }
    }
    in_file.close();
    Printf("Input file loaded...");

    // Load NGen per bin and NGen tot to weigh SL cross sections
    ifstream in_file2(Form("Results/AccAndEffMC/AxE_CorrFC_%ibins_NGen_SL.txt", nPtBins));

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
            if(i == 1) in_stream >> ch >> NGen_tot[0] >> NGen_tot[1];
            if(i > 1) in_stream >> ch >> NGen_bin[0][i-2] >> NGen_bin[1][i-2];
            i++;   
        }
    }
    in_file2.close();
    Printf("Input file loaded...");
    // Cross-check:
    Printf("%.0f", NGen_tot[0]); 
    Printf("%.0f", NGen_bin[1][3]);

    // Calculate FC corrections and print them into text file
    str = Form("%sCoherentContamination_%ibins.txt", path.Data(), nPtBins);
    ofstream outfile(str.Data());
    outfile << Form("Bin \tCohNGen\tIncNGen\tSigCoh \tSigInc \tAxECoh\tErr \tAxEInc\tErr \tfC [%%] \tErr \n");  
    for(Int_t i = 0; i < nPtBins; i++){  
        fC_val[i] = (sig_SL_j_coh * (NGen_bin[0][i]/NGen_tot[0])) / (sig_SL_j_inc * (NGen_bin[1][i]/NGen_tot[1])) * (AxE_CohJ_val[i] / AxE_IncJ_val[i]) * 100; // in percent already
        if(AxE_CohJ_val[i] != 0){
            fC_err[i] = fC_val[i] * TMath::Sqrt(
                TMath::Power((AxE_CohJ_err[i]/AxE_CohJ_val[i]),2) + 
                TMath::Power((AxE_IncJ_err[i]/AxE_IncJ_val[i]),2));
        } else fC_err[i] = 0;
        outfile << std::fixed  << std::setprecision(0);
        outfile << i+1 << "\t" << NGen_bin[0][i] << "\t" 
                               << NGen_bin[1][i] << "\t";
        outfile << std::fixed  << std::setprecision(4)
                               << sig_SL_j_coh * (NGen_bin[0][i]/NGen_tot[0]) << "\t" 
                               << sig_SL_j_inc * (NGen_bin[1][i]/NGen_tot[1]) << "\t";
        outfile << std::fixed  << std::setprecision(2)
                               << AxE_CohJ_val[i] << "\t" 
                               << AxE_CohJ_err[i] << "\t" 
                               << AxE_IncJ_val[i] << "\t" 
                               << AxE_IncJ_err[i] << "\t"                  
                               << fC_val[i] << "\t" 
                               << fC_err[i] << "\n";
    }
    outfile << std::fixed << std::setprecision(0) << "Total" << "\t" << NGen_tot[0] << "\t" << NGen_tot[1] << "\t";
    outfile << std::fixed << std::setprecision(3) << sig_SL_j_coh << "\t" << sig_SL_j_inc;
    outfile.close();
    Printf("*** Results printed to %s.***", str.Data());

    return;
}