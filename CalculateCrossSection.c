// CalculateCrossSection.c
// David Grund, Sep 21, 2021

// cpp headers
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
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

// For temporarily loading bin numbers when reading the text files
Int_t i_bin;

void CalculateCrossSectionTotal();
void CalculateCrossSectionBins();

void CalculateCrossSection(){

    CalculateCrossSectionTotal();

    CalculateCrossSectionBins();

    return;
}

void CalculateCrossSectionTotal(){

    return;
}

void CalculateCrossSectionBins(){

    SetPtBinning();

    // 1) Load N_yield per bin
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        TString NYield_path = Form("Results/InvMassFit/%ibins/bin%i_signal.txt", nPtBins, iBin+1);
        ifstream file_in;
        file_in.open(NYield_path.Data());
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
    TString AxE_path = Form("Results/AccAndEffMC/AxE_%ibins.txt", nPtBins);
    ifstream file_in;
    file_in.open(AxE_path.Data());
    if(!(file_in.fail())){
        // Read data from the file
        Int_t i_line = 0;
        while(!file_in.eof()){
            file_in >> i_bin >> AxE[i_line] >> AxE_err[i_line];
            i_line++;
        }
        file_in.close(); 
        Printf("2) AxE from kIncohJpsiToMu for %ibins loaded.", nPtBins);
    } else {
        Printf("2) ERROR: AxE file missing. Terminating...");
        return;
    }

    // 3) Load all FD corrections per bin
    TString CorrFD_path = Form("Results/FeedDown/FeedDownCorrections_%ibins.txt", nPtBins);
    file_in.open(CorrFD_path.Data());
    if(!(file_in.fail())){
        // Read data from the file
        Int_t i_line = 0;
        while(!file_in.eof()){
            file_in >> i_bin    >> CorrFD[0][i_line] >> CorrFD_err[0][i_line] 
                                >> CorrFD[1][i_line] >> CorrFD_err[1][i_line]
                                >> CorrFD[2][i_line] >> CorrFD_err[2][i_line]
                                >> CorrFD[3][i_line] >> CorrFD_err[3][i_line];
            i_line++;
        }
        file_in.close(); 
        Printf("3) Feed-down corrections for %ibins loaded.", nPtBins);
    } else {
        Printf("3) ERROR: CorrFD file missing. Terminating...");
        return;
    }
    // Calculate total FD per bin
    Double_t CorrFD_total[nPtBins];
    Double_t CorrFD_total_err[nPtBins];
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        Double_t SumOfSquares = 0;
        for(Int_t iFD = 0; iFD < 4; iFD++){
            CorrFD_total[iBin] += CorrFD[iFD][iBin];
            SumOfSquares += TMath::Power(CorrFD_err[iFD][iBin],2);
        }
        CorrFD_total_err[iBin] = TMath::Sqrt(SumOfSquares);
    }

    // 4) Load FC corr per bin
    // (...)

    // 5) Total integrated luminosity
    Double_t LumiAll = (Lumi18q + Lumi18r) * 1000; // 1/(mili barn)
    Printf("5) Total integrated lumi calculated.");

    // 6) Widths of pt intervals [in GeV^2]
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        Pt2Widths[iBin] = TMath::Power(ptBoundaries[iBin+1] - ptBoundaries[iBin], 2);
    }
    Printf("6) pt^2 widths calculated.");

    // Calculate the cross section per bin
    Double_t Factors[nPtBins];
    Double_t Factors_err[nPtBins];
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        Factors[iBin] = 1.0 + CorrFD_total[iBin]/100 + CorrFC[iBin]/100;
        Factors_err[iBin] = TMath::Sqrt(TMath::Power(CorrFD_total_err[iBin]/100,2) + TMath::Power(CorrFC_err[iBin]/100,2));
        Sigma[iBin] = NYield[iBin] / (Factors[iBin] * EffVetoes * AxE[iBin]/100 * BR * RapWidth * Pt2Widths[iBin] * LumiAll);
        Sigma_err[iBin] = Sigma[iBin] * TMath::Sqrt(
            TMath::Power(NYield_err[iBin] / NYield[iBin], 2) +
            TMath::Power(AxE_err[iBin] / AxE[iBin], 2) +
            TMath::Power(Factors_err[iBin] / Factors[iBin], 2) +
            TMath::Power(BR_err / BR, 2));
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
    outfile << "Bin\tPtLow \tPtUpp \tPt^2_w \tN \tN_er \tAxE \tAxE_er\tFD [%%]\tFD_err \tFC [%%]\tFC_er \tf [%%]\tf_er \tsig \tsig_er \n";
    for(Int_t i = 0; i < nPtBins; i++){
        outfile << i << "\t";
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
        outfile << CorrFD_total[i] << "\t";
        outfile << CorrFD_total_err[i] << "\t";
        outfile << CorrFC[i] << "\t";
        outfile << CorrFC_err[i] << "\t";
        outfile << Factors[i] << "\t";
        outfile << Factors_err[i] << "\t";
        outfile << Sigma[i] << "\t";
        outfile << Sigma_err[i] << "\n";
    }
    outfile.close();
    Printf("Results printed to %s.", FilePath.Data()); 

    return;
}