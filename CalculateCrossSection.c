// CalculateCrossSection.c
// David Grund, Sep 21, 2021

// cpp headers
#include <fstream>  // print output to txt file
#include <iomanip>  // std::setprecision()
#include <string>   // getline
// root headers
#include "TMath.h"
// my headers
#include "AnalysisManager.h"

Double_t Lumi18q = 90.114;  // 1/(mu barn)
Double_t Lumi18r = 137.812; // 1/(mu barn)
Double_t BR = 0.05961;
Double_t BR_err = 0.00033;
Double_t RapWidth = 1.6;
Double_t EffVetoes = 0.846; // (!) 
// Old values
Double_t Lumi18q_aod = 91.431;  // 1/(mu barn)
Double_t Lumi18r_aod = 147.292; // 1/(mu barn)
// Cross section in pt bins
Double_t NYield[nPtBins] = { 0 };
Double_t NYield_err[nPtBins] = { 0 };
Double_t Pt2Widths[nPtBins] = { 0 };
Double_t AxE[nPtBins] = { 0 };
Double_t AxE_err[nPtBins] = { 0 };
Double_t CorrFD[4][nPtBins] = { 0 };
Double_t CorrFD_err[4][nPtBins] = { 0 };
Double_t CorrFC[nPtBins] = { 0 };
Double_t CorrFC_err[nPtBins] = { 0 };
Double_t Sigma[nPtBins] = { 0 };
Double_t Sigma_err[nPtBins] = { 0 };
// Total cross section for pt > 0.2 GeV/c
Double_t NYield_tot = 0;
Double_t NYield_tot_err = 0;
Double_t AxE_tot = 0;
Double_t AxE_tot_err = 0;
Double_t CorrFD_tot[4] = { 0 };
Double_t CorrFD_tot_err[4] = { 0 };
Double_t CorrFC_tot = 0;
Double_t CorrFC_tot_err = 0;
Double_t Sigma_tot = 0;
Double_t Sigma_tot_err = 0;
// Other systematic uncertainties (in percent)
Double_t Syst_Lumi = 2.7;
Double_t Syst_V0AD = 3.0;
Double_t Syst_EMdi = 2.0;
Double_t Syst_Trck = 2.8;
Double_t Syst_TrEf = 1.3;

// For temporarily loading bin numbers when reading the text files
Int_t i_bin;

void CalculateCrossSectionTotal(Bool_t bAOD = kFALSE);
void CalculateCrossSectionBins();

void CalculateCrossSection(){

    // ESDs
    CalculateCrossSectionTotal(kFALSE);

    // AODs
    CalculateCrossSectionTotal(kTRUE);

    CalculateCrossSectionBins();

    return;
}

void CalculateCrossSectionTotal(Bool_t bAOD){
    // for pt > 0.2 GeV/c

    ifstream file_in;

    // 1) Load N_yield
    if(!bAOD){
        file_in.open("Results/InvMassFit/inc/inc_signal.txt");
        if(!(file_in.fail())){
            // Read data from the file
            while(!file_in.eof()){
                file_in >> NYield_tot >> NYield_tot_err;
            }
            file_in.close(); 
            Printf("1) NYield file loaded.");
        } else {
            Printf("1) ERROR: NYield file missing. Terminating...");
            return;
        }
    } else {
        NYield_tot = 643;
        NYield_tot_err = 31;         
    }

    // 2) Load AxE
    if(!bAOD) file_in.open("Results/AccAndEffMC/AxE_JInc_MassCut0_PtCut0.txt");
    else file_in.open("Results/AccAndEffMC/AOD/AxE_AOD_JInc_MassCut1_PtCut0.txt");
    if(!(file_in.fail())){
        // Read data from the file
        Int_t i = 0;
        std::string str;
        while(std::getline(file_in,str)){
            istringstream in_stream(str);
            // skip first line
            if(i == 1) in_stream >> AxE_tot >> AxE_tot_err;
            i++;   
        } 
        file_in.close();
        Printf("2) AxE file loaded.");
    } else {
        Printf("2) ERROR: AxE file missing. Terminating...");
        return;            
    }

    // 3) Load all FD corrections
    if(!bAOD) file_in.open("Results/FeedDown/FeedDown_Total.txt");
    else file_in.open("Results/FeedDown/FeedDown_Total_AOD.txt");
    if(!(file_in.fail())){
        // Read data from the file
        Int_t i = 0;
        char ch[8];
        std::string str;
        while(std::getline(file_in,str)){
            istringstream in_stream(str);
            // skip first line
            if(i == 1) in_stream >> ch >> CorrFD_tot[0] >> CorrFD_tot_err[0] 
                                    >> CorrFD_tot[1] >> CorrFD_tot_err[1]
                                    >> CorrFD_tot[2] >> CorrFD_tot_err[2]
                                    >> CorrFD_tot[3] >> CorrFD_tot_err[3];
            i++;   
        } 
        file_in.close(); 
        Printf("3) FD file loaded.");
    } else {
        Printf("3) ERROR: FD file missing. Terminating...");
        return;
    }
    // Calculate total FD
    Double_t CorrFD_sum = 0;
    Double_t CorrFD_sum_err = 0;
    Double_t SumOfSquares = 0;
    for(Int_t iFD = 0; iFD < 4; iFD++){
        CorrFD_sum += CorrFD_tot[iFD];
        SumOfSquares += TMath::Power(CorrFD_tot_err[iFD],2);
    }
    CorrFD_sum_err = TMath::Sqrt(SumOfSquares);

    // 4) Load FC corr
    // (...)

    // 5) Total integrated luminosity
    Double_t LumiAll;
    if(!bAOD) LumiAll = (Lumi18q + Lumi18r) * 1000;     // 1/(mili barn)
    else LumiAll = (Lumi18q_aod + Lumi18r_aod) * 1000;  // 1/(mili barn)
    Printf("5) Total integrated lumi calculated.");

    // Calculate the total cross section
    Double_t Factors;
    Double_t Factors_err;
    Factors = 1.0 + CorrFD_sum / 100 + CorrFC_tot / 100;
    Factors_err = TMath::Sqrt(TMath::Power(CorrFD_sum_err / 100,2) + TMath::Power(CorrFC_tot_err / 100,2));
    Sigma_tot = NYield_tot / (Factors * EffVetoes * AxE_tot/100 * BR * RapWidth * LumiAll);
    Sigma_tot_err = Sigma_tot * TMath::Sqrt(
        TMath::Power(NYield_tot_err / NYield_tot, 2) +
        TMath::Power(AxE_tot_err / AxE_tot, 2) +
        TMath::Power(Factors_err / Factors, 2) +
        TMath::Power(BR_err / BR, 2));

    // Define output text file to print results
    TString FilePath;
    if(!bAOD) FilePath = "Results/CrossSection/total_ESDs_output.txt";
    else FilePath = "Results/CrossSection/total_AODs_output.txt";
    ofstream outfile(FilePath.Data());
    outfile << std::fixed << std::setprecision(3);
    // Print results to the text file
    outfile << Form("Lumi = %.3f 1/(mili barn)\n", LumiAll);
    outfile << Form("BR(J/psi -> mu mu) = (%.3f pm %.3f)%%\n", BR*100, BR_err*100);
    outfile << Form("Delta y = %.1f\n", RapWidth);
    outfile << Form("EffVetoes = %.1f%%\n", EffVetoes*100);
    outfile << Form("N \tN_er \tAxE \tAxE_er\tFD [%%]\tFD_err \tFC [%%]\tFC_er \tf [%%]\tf_er \tsig \tsig_er \n");
    outfile << std::fixed << std::setprecision(2);
    outfile << NYield_tot << "\t";
    outfile << NYield_tot_err << "\t";
    outfile << std::fixed << std::setprecision(3);
    outfile << AxE_tot << "\t";
    outfile << AxE_tot_err << "\t";
    outfile << CorrFD_sum << "\t";
    outfile << CorrFD_sum_err << "\t";
    outfile << CorrFC_tot << "\t";
    outfile << CorrFC_tot_err << "\t";
    outfile << Factors << "\t";
    outfile << Factors_err << "\t";
    outfile << Sigma_tot << "\t";
    outfile << Sigma_tot_err << "\n";
    outfile.close();
    Printf("Results printed to %s.", FilePath.Data()); 

    return;
}

void CalculateCrossSectionBins(){

    SetPtBinning();

    ifstream file_in;

    // 1) Load N_yield per bin
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        file_in.open(Form("Results/InvMassFit/%ibins/bin%i_signal.txt", nPtBins, iBin+1));
        // Read data from the file
        if(!(file_in.fail())){
            while(!file_in.eof()){
                file_in >> NYield[iBin] >> NYield_err[iBin];
            }
        } else {
            Printf("1) ERROR: N_yield file missing. Terminating...");
            return;
        }
        file_in.close(); 
    } 
    Printf("1) N_yield for %ibins loaded.", nPtBins);

    // 2) Load AxE per bin
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        file_in.open(Form("Results/AccAndEffMC/AxE_%ibins/AxE_bin%i_JInc.txt", nPtBins, iBin+1));
        // Read data from the file
        if(!(file_in.fail())){
            Int_t i = 0;
            std::string str;
            while(std::getline(file_in,str)){
                istringstream in_stream(str);
                // skip first line
                if(i == 1) in_stream >> AxE[iBin] >> AxE_err[iBin];
                i++;   
            }
        } else {
            Printf("2) ERROR: AxE file missing. Terminating...");
            return;            
        }
        file_in.close();
    }
    Printf("2) AxE from kIncohJpsiToMu for %ibins loaded.", nPtBins);

    // 3) Load all FD corrections per bin
    file_in.open(Form("Results/FeedDown/FeedDown_%ibins.txt", nPtBins));
    // Read data from the file
    if(!(file_in.fail())){
        Int_t i = 0;
        std::string str;
        while(std::getline(file_in,str)){
            istringstream in_stream(str);
            // skip first line
            if(i > 0) in_stream >> i_bin >> CorrFD[0][i-1] >> CorrFD_err[0][i-1] 
                                >> CorrFD[1][i-1] >> CorrFD_err[1][i-1]
                                >> CorrFD[2][i-1] >> CorrFD_err[2][i-1]
                                >> CorrFD[3][i-1] >> CorrFD_err[3][i-1];
            // Cross-check:
            //Printf("%.4f", CorrFD[0][i-1]);
            i++;   
        }
        Printf("3) FD corrections for %ibins loaded.", nPtBins);
        file_in.close();
    } else {
        Printf("3) ERROR: FD file missing. Terminating...");
        return;
    }
    // Calculate total FD per bin
    Double_t CorrFD_sum[nPtBins] = { 0 };
    Double_t CorrFD_sum_err[nPtBins] = { 0 };
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        Double_t SumOfSquares = 0;
        for(Int_t iFD = 0; iFD < 4; iFD++){
            CorrFD_sum[iBin] += CorrFD[iFD][iBin];
            SumOfSquares += TMath::Power(CorrFD_err[iFD][iBin],2);
        }
        CorrFD_sum_err[iBin] = TMath::Sqrt(SumOfSquares);
    }

    // 4) Load FC corr per bin
    file_in.open(Form("Results/PtFitWithoutBkg/CohContamination_Binning2_%ibins.txt", nPtBins));
    // Read data from the file
    if(!(file_in.fail())){
        Int_t i = 0;
        std::string str;
        while(std::getline(file_in,str)){
            istringstream in_stream(str);
            // skip first line
            if(i > 0) in_stream >> i_bin >> CorrFC[i-1] >> CorrFC_err[i-1];
            // Cross-check:
            //Printf("%.4f", CorrFC[i-1]);
            i++;   
        }
        Printf("4) FC corrections for %ibins loaded.", nPtBins);
        file_in.close();
    } else {
        Printf("4) ERROR: FC file missing. Terminating...");
        return;
    }

    // 5) Total integrated luminosity
    Double_t LumiAll = (Lumi18q + Lumi18r) * 1000; // 1/(mili barn)
    Printf("5) Total integrated lumi calculated.");

    // 6) Widths of pt intervals [in GeV^2]
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        Pt2Widths[iBin] = TMath::Power(ptBoundaries[iBin+1], 2) - TMath::Power(ptBoundaries[iBin], 2);
    }
    Printf("6) pt^2 widths calculated.");

    // Calculate the cross section per bin
    Double_t Factors[nPtBins];
    Double_t Factors_err[nPtBins];
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        Factors[iBin] = 1.0 + CorrFD_sum[iBin]/100 + CorrFC[iBin]/100;
        Factors_err[iBin] = TMath::Sqrt(TMath::Power(CorrFD_sum_err[iBin]/100,2) + TMath::Power(CorrFC_err[iBin]/100,2));
        Sigma[iBin] = NYield[iBin] / (Factors[iBin] * EffVetoes * AxE[iBin]/100 * BR * RapWidth * Pt2Widths[iBin] * LumiAll);
        Sigma_err[iBin] = Sigma[iBin] * TMath::Sqrt(
            TMath::Power(NYield_err[iBin] / NYield[iBin], 2)
            //TMath::Power(AxE_err[iBin] / AxE[iBin], 2) +
            //TMath::Power(Factors_err[iBin] / Factors[iBin], 2) +
            //TMath::Power(BR_err / BR, 2) + 
            //TMath::Power(Syst_Lumi / 100., 2) + 
            //TMath::Power(Syst_V0AD / 100., 2) + 
            //TMath::Power(Syst_EMdi / 100., 2) + 
            //TMath::Power(Syst_Trck / 100., 2) + 
            //TMath::Power(Syst_TrEf / 100., 2)
            );
    }

    // Define output text file to print results
    TString FilePath = Form("Results/CrossSection/%ibins_output.txt", nPtBins);
    ofstream outfile(FilePath.Data());
    outfile << std::fixed << std::setprecision(3);
    // Print results to the text file
    outfile << Form("Lumi = %.3f 1/(mili barn)\n", LumiAll);
    outfile << Form("BR(J/psi -> mu mu) = (%.3f pm %.3f)%%\n", BR*100, BR_err*100);
    outfile << Form("Delta y = %.1f\n", RapWidth);
    outfile << Form("EffVetoes = %.1f%%\n", EffVetoes*100);
    outfile << "Per bins:\n";
    outfile << Form("Bin\tPtLow \tPtUpp \tPt^2_w \tN \tN_er \tAxE \tAxE_er\tFD [%%]\tFD_err \tFC [%%]\tFC_er \tf [%%]\tf_er \tsig \tsig_er stat.\n");
    for(Int_t i = 0; i < nPtBins; i++){
        outfile << i+1 << "\t";
        outfile << ptBoundaries[i] << "\t";
        outfile << ptBoundaries[i+1] << "\t";
        outfile << std::fixed << std::setprecision(4);
        outfile << Pt2Widths[i] << "\t";
        outfile << std::fixed << std::setprecision(2);
        outfile << NYield[i] << "\t";
        outfile << NYield_err[i] << "\t";
        outfile << std::fixed << std::setprecision(3);
        outfile << AxE[i] << "\t";
        outfile << AxE_err[i] << "\t";
        outfile << CorrFD_sum[i] << "\t";
        outfile << CorrFD_sum_err[i] << "\t";
        outfile << CorrFC[i] << "\t";
        outfile << CorrFC_err[i] << "\t";
        outfile << Factors[i] << "\t";
        outfile << Factors_err[i] << "\t";
        outfile << Sigma[i] << "\t";
        outfile << Sigma_err[i] << "\n";
    }
    outfile.close();
    Printf("Results printed to %s.", FilePath.Data()); 

    TString FilePath2 = Form("Results/CrossSection/%ibins_values_plot.txt", nPtBins);
    ofstream outfile2(FilePath2.Data());
    outfile2 << std::fixed << std::setprecision(4);
    outfile2 << "Bin \ttLow \ttUpp \tSig \tErr \n";
    for(Int_t i = 0; i < nPtBins; i++){
        outfile2 << i+1 << "\t" << ptBoundaries[i] * ptBoundaries[i] << "\t" 
                                << ptBoundaries[i+1] * ptBoundaries[i+1] << "\t" 
                                << Sigma[i] << "\t"
                                << Sigma_err[i] << "\n";
    }
    outfile2.close();
    Printf("Results printed to %s.", FilePath2.Data()); 

    return;
}