// PlotCrossSection.c
// David Grund, Oct 06, 2021

// cpp headers
#include <fstream>
#include <string>   // getline
// root headers
#include "TH1.h"
#include "TCanvas.h"
// my headers
#include "AnalysisManager.h"

Int_t iBin;
Double_t tLow[nPtBins], tUpp[nPtBins], Sig[nPtBins], Sig_err[nPtBins];

void PlotCrossSection(){

    SetPtBinning();

    ifstream file_in;
    file_in.open(Form("Results/CrossSection/%ibins_values_plot.txt", nPtBins));
    // Read data from the file
    if(!(file_in.fail())){
        Int_t i = 0;
        std::string str;
        while(std::getline(file_in,str)){
            istringstream in_stream(str);
            // skip first line
            if(i > 0) in_stream >> iBin >> tLow[i-1] >> tUpp[i-1] >> Sig[i-1] >> Sig_err[i-1];
            i++;   
        }
        Printf("Values loaded.");
        file_in.close();
    } else {
        Printf("Input file missing. Terminating...");
        return;
    }

    TH1D *hSigma_Pt = new TH1D("hSigma_Pt", "hSigma_Pt", nPtBins, ptBoundaries);
    Double_t tBoundaries[nPtBins+1];
    for(Int_t i = 0; i < nPtBins+1; i++){
        tBoundaries[i] = ptBoundaries[i]*ptBoundaries[i];
    }
    TH1D *hSigma_T = new TH1D("hSigma_T", "hSigma_T", nPtBins, tBoundaries);
    // Fill histograms
    for(Int_t i = 0; i < nPtBins; i++){
        hSigma_Pt->SetBinContent(i+1, Sig[i]);
        hSigma_Pt->SetBinError(i+1, Sig_err[i]);
        hSigma_T->SetBinContent(i+1, Sig[i]);
        hSigma_T->SetBinError(i+1, Sig_err[i]);
    }

    TCanvas *c_Pt = new TCanvas("c_Pt", "c_Pt", 900, 600);
    c_Pt->SetLogy();
    //c_Pt->SetLogx();
    hSigma_Pt->Draw();

    TCanvas *c_T = new TCanvas("c_T", "c_T", 900, 600);
    c_T->SetLogy();
    //c_T->SetLogx();
    hSigma_T->Draw();

    return;
}