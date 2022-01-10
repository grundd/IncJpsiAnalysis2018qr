// STARlight_NewPtShape_All.c
// David Grund, Jan 7, 2022

// cpp headers
#include <stdio.h>
// root headers
#include "TH1.h"
#include "TStyle.h"
#include "TCanvas.h"
//#include "TLegend.h"
// my headers
#include "AnalysisManager.h"
#include "STARlight_Utilities.h"

Double_t fPtGenerated;
Double_t fPtGenerated_443;
Double_t pT_low; // GeV/c
Double_t pT_upp; // GeV/c
Int_t nBins;

TH1D *hRecOld = NULL;
TH1D *hRatios = NULL;
TH1D *hRecNew = NULL;

void FillHistGen(const char* opt);
void FillHistRec(const char* opt);
void CalculateAndPlotRatios(const char* opt);

TString sIn = "";
TString sOut_subfolder = "";

void STARlight_NewPtShape_all()
{
    char MC[8] = "IncJ";
    FillHistGen(MC);
    FillHistRec(MC);
    CalculateAndPlotRatios(MC);

    return;
}

void FillHistGen(const char* opt)
{
    //*************************
    // opt == "CohJ" => kCohJpsiToMu
    // opt == "CohP" => kCohPsi2sToMuPi
    // opt == "IncJ" => kIncJpsiToMu
    // opt == "IncP" => kIncPsi2sToMuPi
    //*************************

    Printf("\n");
    Printf("*****");
    Printf("Filling hGen for %s, nGen = 2.000.000.", opt);
    Printf("*****");

    // https://www.cplusplus.com/reference/cstring/strncmp/
    if(strncmp(opt,"CohJ",4) == 0){
        pT_low = 0.0; pT_upp = 0.4; nBins = 80;
    } else if(strncmp(opt,"CohP",4) == 0){
        pT_low = 0.0; pT_upp = 0.6; nBins = 120;
    } else if(strncmp(opt,"IncJ",4) == 0){
        pT_low = 0.0; pT_upp = 1.2; nBins = 240;
    } else if(strncmp(opt,"IncP",4) == 0){
        pT_low = 0.0; pT_upp = 0.8; nBins = 160;
    }

    // #############################################################################################
    // Load generated MC events with R_A = 6.624 fm
    TFile *fOld = NULL;
    fOld = TFile::Open(Form("SimAliDPG/%s_6.624/tPtGen.root", opt), "read");
    //fOld = TFile::Open(Form("SimAliDPG/%s_6.624_new/tPtGen.root", opt), "read");
    if(fOld) Printf("File %s loaded.", fOld->GetName());
    TTree *tOld = NULL;
    tOld = dynamic_cast<TTree*> (fOld->Get("tPtGen"));
    tOld->SetBranchAddress("fPtGen", &fPtGenerated);
    //tOld->SetBranchAddress("fPtGen_443", &fPtGenerated_443);
    if(tOld) Printf("Tree %s loaded.", tOld->GetName());
    // Load generated MC events with R_A = 7.350 fm
    TFile *fNew = NULL;
    fNew = TFile::Open(Form("SimAliDPG/%s_7.350/tPtGen.root", opt), "read");
    //fNew = TFile::Open(Form("SimAliDPG/%s_7.350_new/tPtGen.root", opt), "read");
    if(fNew) Printf("File %s loaded.", fNew->GetName());
    TTree *tNew = NULL;
    tNew = dynamic_cast<TTree*> (fNew->Get("tPtGen"));
    tNew->SetBranchAddress("fPtGen", &fPtGenerated);
    //tNew->SetBranchAddress("fPtGen_443", &fPtGenerated_443);
    if(tNew) Printf("Tree %s loaded.", tNew->GetName());
    
    // #############################################################################################
    // Check if output already created, if not, create it
    TString str_fOut = Form("Trees/STARlight/SimAliDPG/%s.root", opt);
    TFile *fOut = TFile::Open(str_fOut.Data(),"read");

    if(fOut){

        Printf("%s already created.", str_fOut.Data());

    } else {

        fOut = new TFile(str_fOut.Data(),"RECREATE");

        TH1D *hGenOld = new TH1D("hGenOld","hGenOld",nBins,pT_low,pT_upp);
        TH1D *hGenNew = new TH1D("hGenNew","hGenNew",nBins,pT_low,pT_upp);

        Printf("Filling hGenOld.");
        Printf("tOld contains %lli entries.", tOld->GetEntries());

        // Loop over entries in tOld
        Int_t nEntriesAnalysed = 0;
        Int_t nEntriesProgress = (Double_t)tOld->GetEntries() / 20.;
        Int_t nPercent = 0;

        for(Int_t iEntry = 0; iEntry < tOld->GetEntries(); iEntry++){
            tOld->GetEntry(iEntry);

            hGenOld->Fill(fPtGenerated);

            // Update progress bar
            if((iEntry+1) % nEntriesProgress == 0){
                nPercent += 5;
                nEntriesAnalysed += nEntriesProgress;
                Printf("[%i%%] %i entries analysed.", nPercent, nEntriesAnalysed);
            }
        }

        Printf("Filling hGenNew.");
        Printf("tNew contains %lli entries.", tNew->GetEntries());

        // Loop over entries in tOld
        nEntriesAnalysed = 0;
        nEntriesProgress = (Double_t)tNew->GetEntries() / 20.;
        nPercent = 0;

        for(Int_t iEntry = 0; iEntry < tNew->GetEntries(); iEntry++){
            tNew->GetEntry(iEntry);

            hGenNew->Fill(fPtGenerated);

            // Update progress bar
            if((iEntry+1) % nEntriesProgress == 0){
                nPercent += 5;
                nEntriesAnalysed += nEntriesProgress;
                Printf("[%i%%] %i entries analysed.", nPercent, nEntriesAnalysed);
            }
        }

        fOut->Write("",TObject::kWriteDelete);
    }
    // #############################################################################################

    Printf("*****");
    Printf("Done.");
    Printf("*****");

    return;
}

void FillHistRec(const char* opt)
{
    //*************************
    // opt == "CohJ" => kCohJpsiToMu
    // opt == "CohP" => kCohPsi2sToMuPi
    // opt == "IncJ" => kIncJpsiToMu
    // opt == "IncP" => kIncPsi2sToMuPi
    //*************************

    hRecOld = new TH1D("hRecOld","hRecOld",nBins,pT_low,pT_upp);

    Printf("\n");
    Printf("*****");
    Printf("Filling hRec for %s.", opt);
    Printf("*****");

    // Load reconstructed MC events 
    TFile *fRec = NULL;
    // https://www.cplusplus.com/reference/cstring/strncmp/
    if(strncmp(opt,"CohJ",4) == 0){
        fRec = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohJpsiToMu.root", "read");
    } else if(strncmp(opt,"CohP",4) == 0){
        fRec = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohPsi2sToMuPi.root", "read");
    } else if(strncmp(opt,"IncJ",4) == 0){
        fRec = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohJpsiToMu.root", "read");
    } else if(strncmp(opt,"IncP",4) == 0){
        fRec = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohPsi2sToMuPi.root", "read");
    }
    if(fRec) Printf("File %s loaded.", fRec->GetName());
    // Get MC_rec tree
    TTree *tRec = dynamic_cast<TTree*> (fRec->Get("AnalysisOutput/fTreeJPsiMCRec"));
    if(tRec) Printf("Tree %s loaded.", tRec->GetName());
    ConnectTreeVariablesMCRec(tRec);

    Int_t nEntriesAnalysed = 0;
    Int_t nEntriesProgress = (Double_t)tRec->GetEntries() / 20.;
    Int_t nPercent = 0;

    // Prepare reconstructed events
    for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++){
        tRec->GetEntry(iEntry);
        if(EventPassedMCRec(1, 2)){
            hRecOld->Fill(fPt);
        } 

        if((iEntry+1) % nEntriesProgress == 0){
            nPercent += 5;
            nEntriesAnalysed += nEntriesProgress;
            Printf("[%i%%] %i entries analysed.", nPercent, nEntriesAnalysed);
        }
    }

    Printf("*****");
    Printf("Done.");
    Printf("*****");

    return;
}

void CalculateAndPlotRatios(const char* opt)
{

    Printf("\n");
    Printf("*****");
    Printf("Calculating ratios for %s.", opt);
    Printf("*****");

    // Load fGen
    TFile *fGen = TFile::Open(Form("Trees/STARlight/SimAliDPG/%s.root", opt), "read");
    if(fGen) Printf("File %s loaded.", fGen->GetName());

    // Load histograms with generated events
    TH1D *hGenOld = (TH1D*)fGen->Get("hGenOld");
    if(hGenOld) Printf("Histogram %s loaded.", hGenOld->GetName());
    TH1D *hGenNew = (TH1D*)fGen->Get("hGenNew");
    if(hGenNew) Printf("Histogram %s loaded.", hGenNew->GetName());

    hRatios = (TH1D*)hGenNew->Clone("hRatios");
    hRatios->SetTitle("hRatios");
    hRatios->Sumw2();
    hRatios->Divide(hGenOld);

    hRecNew = (TH1D*)hRecOld->Clone("hRecNew");
    hRecNew->SetTitle("hRecNew");
    hRecNew->Multiply(hRatios);

    // Print the results to console
    Printf("*****");
    Printf("Output:");
    Printf("pT_low\tpT_upp\tnEvOld\tnEvNew\tRatio");
    for(Int_t iBin = 1; iBin <= nBins; iBin++){
        Printf("%.3f\t%.3f\t%.0f\t%.0f\t%.3f",
            hRatios->GetBinLowEdge(iBin), hRatios->GetBinLowEdge(iBin+1), 
            hGenOld->GetBinContent(iBin), hGenNew->GetBinContent(iBin), hRatios->GetBinContent(iBin));
    }
    Printf("*****");

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");

    // Plot the results
    TCanvas *cRatios = new TCanvas("cRatios", "cRatios", 900, 600);
    // Canvas settings
    cRatios->SetTopMargin(0.04);
    cRatios->SetBottomMargin(0.15);
    cRatios->SetRightMargin(0.04);
    cRatios->SetLeftMargin(0.1);
    // Vertical axis
    hRatios->GetYaxis()->SetTitle("#it{N}_{gen}(#it{R}_{A} = 7.350 fm) / #it{N}_{gen}(#it{R}_{A} = 6.624 fm)");
    hRatios->GetYaxis()->SetTitleSize(0.05);
    hRatios->GetYaxis()->SetTitleOffset(0.9);
    hRatios->GetYaxis()->SetLabelSize(0.05);
    hRatios->GetYaxis()->SetDecimals(1);
    hRatios->GetYaxis()->SetRangeUser(0.0,5.0);
    // Horizontal axis
    hRatios->GetXaxis()->SetTitle("#it{p}_{T}^{gen} (GeV/#it{c})");
    hRatios->GetXaxis()->SetTitleSize(0.05);
    hRatios->GetXaxis()->SetTitleOffset(1.2);
    hRatios->GetXaxis()->SetLabelSize(0.05);
    // Draw
    hRatios->SetLineColor(215);
    hRatios->Draw("E1");

    TCanvas *cRec = new TCanvas("cRec", "cRec", 900, 600);
    cRec->SetTopMargin(0.04);
    cRec->SetBottomMargin(0.15);
    cRec->SetRightMargin(0.04);
    cRec->SetLeftMargin(0.1);
    cRec->SetLogy();
    // Vertical axis
    hRecOld->GetYaxis()->SetTitle("#it{N}_{rec} (after selections)");
    hRecOld->GetYaxis()->SetTitleSize(0.05);
    hRecOld->GetYaxis()->SetTitleOffset(0.9);
    hRecOld->GetYaxis()->SetLabelSize(0.05);
    hRecOld->GetYaxis()->SetRangeUser(0.1,1e5);
    // Horizontal axis
    hRecOld->GetXaxis()->SetTitle("#it{p}_{T}^{gen} (GeV/#it{c})");
    hRecOld->GetXaxis()->SetTitleSize(0.05);
    hRecOld->GetXaxis()->SetTitleOffset(1.2);
    hRecOld->GetXaxis()->SetLabelSize(0.05);
    // Draw
    hRecOld->SetLineColor(210);
    hRecOld->Draw("E1");
    hRecNew->SetLineColor(215);
    hRecNew->Draw("E1 SAME");

    // Print plots
    cRatios->Print(Form("Results/STARlight/SimAliDPG/Ratios/%s_7.350.pdf", opt));
    cRatios->Print(Form("Results/STARlight/SimAliDPG/Ratios/%s_7.350.png", opt));
    cRec->Print(Form("Results/STARlight/SimAliDPG/RecSpectra/%s_7.350.pdf", opt));
    cRec->Print(Form("Results/STARlight/SimAliDPG/RecSpectra/%s_7.350.png", opt));

    delete hGenOld;
    delete hGenNew;

    Printf("*****");
    Printf("Done.");
    Printf("*****");

    return;
}