// AccAndEffMC.c
// David Grund, 15-09-2021
// To calculate the acceptance x efficiency from MC data in defined pt bins

// cpp headers
#include <fstream> // print output to txt file
//#include <iomanip> // std::setprecision()
// root headers
#include "TH1.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h" // gStyle
#include "TLegend.h"
#include "TMath.h"
// my headers
#include "AnalysisManager.h"

TH1D* hNRec = NULL; 
TH1D* hNGen = NULL; 
TH1D* hAxE = NULL;

void CalculateAxE(Int_t iPtCut);
void CalculateAxE_AOD(Int_t iPtCut);
void CalculateAxEPtBins();
void FillHistNRec();
void FillHistNGen();
void SaveToFile(TH1D* hist, TString name);
Double_t CalculateErrorBayes(Double_t k, Double_t n);

void AccAndEffMC(){

    Bool_t bTotalESD = kTRUE;
    if(bTotalESD){
        CalculateAxE(0); // pt > 0.2 GeV/c
        CalculateAxE(3); // 0.2 < pt < 1 GeV/c
    }

    Bool_t bTotalAOD = kFALSE;
    if(bTotalAOD){
        CalculateAxE_AOD(0); // pt > 0.2 GeV/c
        CalculateAxE_AOD(3); // 0.2 < pt < 1 GeV/c
    }
    
    CalculateAxEPtBins();

    return;
}

void CalculateAxE(Int_t iPtCut){

    TFile *fRec = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohJpsiToMu.root", "read");
    //TFile *fRec = TFile::Open("Trees/AnalysisDataMC/old (from GRID)/AnalysisResults_kIncohJpsiToMu.root", "read");
    if(fRec) Printf("MC rec file loaded.");

    TTree *tRec = dynamic_cast<TTree*> (fRec->Get("AnalysisOutput/fTreeJPsiMCRec"));
    if(tRec) Printf("MC rec tree loaded.");
    
    ConnectTreeVariablesMCRec(tRec);

    Double_t NRec = 0;
    for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++){
        tRec->GetEntry(iEntry);
        if(EventPassedMCRec(0, iPtCut)) NRec++;
    }

    TTree *tGen = dynamic_cast<TTree*> (fRec->Get("AnalysisOutput/fTreeJPsiMCGen"));
    if(tGen) Printf("MC gen tree loaded.");
    
    ConnectTreeVariablesMCGen(tGen);

    Double_t NGen = 0;
    for(Int_t iEntry = 0; iEntry < tGen->GetEntries(); iEntry++){
        tGen->GetEntry(iEntry);
        if(EventPassedMCGen(iPtCut)) NGen++;
    }

    Printf("N_rec = %.0f", NRec);
    Printf("N_gen = %.0f", NGen);

    Double_t AxE = NRec / NGen;
    Double_t AxE_err = CalculateErrorBayes(NRec, NGen);

    Printf("AxE = (%.4f pm %.4f)%%", AxE*100, AxE_err*100);

    return;
}

void CalculateAxE_AOD(Int_t iPtCut){
    // comment the SPD matching in EventPassedMCRec when working with AODs

    TFile *fRec = TFile::Open("Trees/AnalysisDataAOD/MC_rec_18qr_kIncohJpsiToMu_migr.root", "read");
    if(fRec) Printf("MC rec file loaded.");

    TTree *tRec = dynamic_cast<TTree*> (fRec->Get("analysisTree"));
    if(tRec) Printf("MC rec tree loaded.");
    
    ConnectTreeVariablesMCRec_AOD(tRec);

    Double_t NRec = 0;
    for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++){
        tRec->GetEntry(iEntry);
        if(EventPassedMCRec_AOD(0, iPtCut)) NRec++;
    }

    TFile *fGen = TFile::Open("Trees/AnalysisDataAOD/MC_gen_18qr_kIncohJpsiToMu_migr.root", "read");
    if(fGen) Printf("MC gen file loaded.");

    TTree *tGen = dynamic_cast<TTree*> (fGen->Get("MCgenTree"));
    if(tGen) Printf("MC gen tree loaded.");
    
    ConnectTreeVariablesMCGen_AOD(tGen);

    Double_t NGen = 0;
    for(Int_t iEntry = 0; iEntry < tGen->GetEntries(); iEntry++){
        tGen->GetEntry(iEntry);
        if(EventPassedMCGen(iPtCut)) NGen++;
    }

    Printf("N_rec = %.0f", NRec);
    Printf("N_gen = %.0f", NGen);

    Double_t AxE = NRec / NGen;
    Double_t AxE_err = CalculateErrorBayes(NRec, NGen);

    Printf("AxE = (%.4f pm %.4f)%%", AxE*100, AxE_err*100);

    return;
}

void CalculateAxEPtBins(){

    SetPtBinning();

    hNRec = new TH1D("hNRec","N rec per bin",nPtBins,ptBoundaries);
    hNGen = new TH1D("hNGen","N gen per bin",nPtBins,ptBoundaries);

    FillHistNRec();
    FillHistNGen();

    hAxE = (TH1D*)hNRec->Clone("hAxE");
    hAxE->SetTitle("AxE per bin");
    hAxE->Sumw2();
    hAxE->Divide(hNGen);

    // Draw the histogram:
    TCanvas *c = new TCanvas("c", "c", 900, 600);
    c->SetTopMargin(0.02);
    c->SetBottomMargin(0.14);
    c->SetRightMargin(0.03);
    c->SetLeftMargin(0.145);
    // gStyle
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");
    // Marker and line
    //hAxE->SetMarkerStyle(21);
    //hAxE->SetMarkerColor(kBlue);
    //hAxE->SetMarkerSize(1.0);
    hAxE->SetLineColor(kBlue);
    hAxE->SetLineWidth(1.0);
    // Vertical axis
    hAxE->GetYaxis()->SetTitle("#it{N}_{rec}/#it{N}_{gen}");
    hAxE->GetYaxis()->SetTitleSize(0.056);
    hAxE->GetYaxis()->SetTitleOffset(1.3);
    hAxE->GetYaxis()->SetLabelSize(0.056);
    hAxE->GetYaxis()->SetDecimals(3);
    hAxE->GetYaxis()->SetRangeUser(0.0,hAxE->GetBinContent(1)*1.1);
    // Horizontal axis
    hAxE->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hAxE->GetXaxis()->SetTitleSize(0.056);
    hAxE->GetXaxis()->SetTitleOffset(1.2);
    hAxE->GetXaxis()->SetLabelSize(0.056);
    hAxE->GetXaxis()->SetLabelOffset(0.015);
    hAxE->GetXaxis()->SetDecimals(1);
    // Eventually draw it
    hAxE->Draw("P E1");
    // Legend
    TLegend *l = new TLegend(0.52,0.77,0.85,0.97);
    l->AddEntry((TObject*)0,Form("ALICE Simulation"),""); 
    l->AddEntry((TObject*)0,Form("Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV"),"");
    l->AddEntry((TObject*)0,Form("inc J/#psi #rightarrow #mu^{+}#mu^{-}"),"");
    l->SetTextSize(0.056);
    l->SetBorderSize(0); // no border
    l->SetFillStyle(0);  // legend is transparent
    l->Draw();
    // Legend 2
    TLegend *l2 = new TLegend(0.15,0.17,0.35,0.32);
    l2->AddEntry((TObject*)0,Form("|#it{y}| < 0.8"),""); 
    l2->AddEntry((TObject*)0,Form("2.2 < #it{m} < 4.5 GeV/#it{c}^{2}"),"");
    l2->SetTextSize(0.056);
    l2->SetBorderSize(0); // no border
    l2->SetFillStyle(0);  // legend is transparent
    l2->Draw();

    // Save the figures and print the results to txt file
    TString str;
    if(nPtBins == 4) str = "Results/AccAndEffMC/AxE_4bins";
    if(nPtBins == 5) str = "Results/AccAndEffMC/AxE_5bins";
    c->Print((str + ".pdf").Data());
    c->Print((str + ".png").Data());
    ofstream outfile((str + ".txt").Data());
    outfile << std::fixed << std::setprecision(5);
    outfile << "Bin \tAxE [%%] \tAxE_err [%%] \n";
    for(Int_t i = 1; i <= nPtBins; i++){
        outfile << i << "\t" << hAxE->GetBinContent(i)*100 << "\t\t" << hAxE->GetBinError(i)*100 << "\n";
    }
    outfile.close();
    Printf("*** Results printed to %s.***", (str + ".txt").Data());

    // Compare errors that Root gives with CalculateErrorBayes
    Bool_t DebugErrors = kFALSE;
    if(DebugErrors){
        Double_t ErrRoot = 0;
        Double_t ErrBayes = 0;    
        for(Int_t i = 1; i <= nPtBins; i++){
            ErrRoot = hAxE->GetBinError(i);
            ErrBayes = CalculateErrorBayes(hNRec->GetBinContent(i),hNGen->GetBinContent(i));
            Printf("Root: %.5f, Bayes: %.5f", ErrRoot, ErrBayes);
        }
    }

    // Cross-check: calculate the total value of AxE
    Double_t NRecTot = 0;
    Double_t NGenTot = 0;
    for(Int_t i = 1; i <= nPtBins; i++){
        NRecTot += hNRec->GetBinContent(i);
        NGenTot += hNGen->GetBinContent(i);
    }
    Double_t AxETot = NRecTot / NGenTot;
    Double_t AxETot_err = CalculateErrorBayes(NRecTot, NGenTot);
    Printf("Total AxE = (%.4f pm %.4f)%%", AxETot*100, AxETot_err*100);

    return;
}

void FillHistNRec(){
    // Check if the corresponding text file already exists
    TString file("Results/AccAndEffMC/");
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
    TString file("Results/AccAndEffMC/");
    if(nPtBins == 4) file.Append("NGen_4bins.txt");
    if(nPtBins == 5) file.Append("NGen_5bins.txt");

    ifstream inFile;
    inFile.open(file);
    if(!(inFile.fail())){
        // This configuration has already been calculated
        Printf("*** The file %s already exists. ***", file.Data());
        // Fill hNGen with data from the text file
        Int_t inBin;
        Double_t inValue;
        while(!inFile.eof()){
            inFile >> inBin >> inValue; // fist and second column
            hNGen->SetBinContent(inBin, inValue);
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
                if(EventPassedMCGen(4, iPtBin)) NGen++;
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

Double_t CalculateErrorBayes(Double_t k, Double_t n){ // k = NRec, n = NGen

    Double_t var = (k + 1) * (k + 2) / (n + 2) / (n + 3) - (k + 1) * (k + 1) / (n + 2) / (n + 2);
    Double_t err = TMath::Sqrt(var);

    return err;
}