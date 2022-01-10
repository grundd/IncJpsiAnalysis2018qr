// STARlight_NewPtShape_inc.c
// David Grund, Nov 23, 2021
// To study the effect of changing R_A in STARlight code to 7.53 fm
// on the shape of coherent J/psi -> mu mu events
// instructions: see the email from Guillermo, Nov 15, 2021

// root headers
#include "TH1.h"
#include "TStyle.h"
#include "TCanvas.h"
//#include "TLegend.h"
// my headers
#include "AnalysisManager.h"
#include "STARlight_Utilities.h"

Double_t fPtGenerated;
Double_t pT_low = 0.0; // GeV/c
Double_t pT_upp = 1.6; // GeV/c
Int_t nBins = 320;

TH1D *hRecOld = new TH1D("hRecOld","hRecOld",nBins,pT_low,pT_upp);
TH1D *hRatios = NULL;
TH1D *hRecNew = NULL;

void FillHistRec();
void FillTreeGen();
void CalculateAndPlotRatios();

void STARlight_NewPtShape_inc()
{
    FillHistRec();
    FillTreeGen();
    CalculateAndPlotRatios();

    return;
}

void FillTreeGen()
{

    Printf("*****");
    Printf("Filling trees with nGen = 2.000.000.");
    Printf("*****");

    Double_t nEvOld = 0;
    Double_t nEvNew = 0;
    
    // #############################################################################################
    
    TFile *fSLOld = NULL;
    TTree *tSLOld = NULL;
    // Load generated coherent MC events with R_A = 6.624 fm
    fSLOld = TFile::Open("Trees/STARlight/inc_2000000_stdRA/trees_starlight.root", "read");
    if(fSLOld) Printf("File %s loaded.", fSLOld->GetName());
    // Get the SL tree
    tSLOld = dynamic_cast<TTree*> (fSLOld->Get("starlightTree"));
    if(tSLOld) Printf("Tree %s loaded.", tSLOld->GetName());
    ConnectTreeVariables_tSL(tSLOld);
    TFile *fSLNew = NULL;
    TTree *tSLNew = NULL;
    // NGen = 6.000.000
    // Load generated coherent MC events with R_A = 7.350 fm
    fSLNew = TFile::Open("Trees/STARlight/inc_2000000_stdRA/trees_starlight.root", "read");
    if(fSLNew) Printf("File %s loaded.", fSLNew->GetName());
    // Get the SL tree
    tSLNew = dynamic_cast<TTree*> (fSLNew->Get("starlightTree"));
    if(tSLNew) Printf("Tree %s loaded.", tSLNew->GetName());
    ConnectTreeVariables_tSL(tSLNew);

    // #############################################################################################
    // Check if output already created, if not, create it

    TFile *fGenOld = TFile::Open("Trees/STARlight/IncoherentShape/tGenOld_2000k_6.624.root","read");

    if(fGenOld){

        Printf("File already created.");

    } else {

        fGenOld = new TFile("Trees/STARlight/IncoherentShape/tGenOld_2000k_6.624.root","RECREATE");

        TH1D *hGenOld = new TH1D("hGenOld","hGenOld",nBins,pT_low,pT_upp);

        TTree *tGenOld = new TTree("tGenOld","tGenOld");
        tGenOld->Branch("fPtGen",&fPtGenerated,"fPtGen/D");

        Printf("Filling tGenOld and hGenOld.");
        Printf("Old tree contains %lli entries.", tSLOld->GetEntries());

        // Loop over entries in tGen
        Int_t nEntriesAnalysed = 0;
        Int_t nEntriesProgress = (Double_t)tSLOld->GetEntries() / 20.;
        Int_t nPercent = 0;

        for(Int_t iEntry = 0; iEntry < tSLOld->GetEntries(); iEntry++){
            tSLOld->GetEntry(iEntry);
            if(TMath::Abs(fYGen) < 1.0){
                nEvOld++;
                fPtGenerated = parent->Pt();
                hGenOld->Fill(fPtGenerated);
                tGenOld->Fill();
            }
            // Update progress bar
            if((iEntry+1) % nEntriesProgress == 0){
                nPercent += 5;
                nEntriesAnalysed += nEntriesProgress;
                Printf("[%i%%] %i entries analysed.", nPercent, nEntriesAnalysed);
            }
        }
        Printf("No. of events with |y| < 1.0: %.0f", nEvOld);

        fGenOld->Write("",TObject::kWriteDelete);
    }

    TFile *fGenNew = TFile::Open("Trees/STARlight/IncoherentShape/tGenNew_2000k_7.350.root","read");

    if(fGenNew){

        Printf("File already created.");

    } else {

        fGenNew = new TFile("Trees/STARlight/IncoherentShape/tGenNew_2000k_7.350.root","RECREATE");

        TH1D *hGenNew = new TH1D("hGenNew","hGenNew",nBins,pT_low,pT_upp);

        TTree *tGenNew = new TTree("tGenNew","tGenNew");
        tGenNew->Branch("fPtGen",&fPtGenerated,"fPtGen/D");

        Printf("Filling tGenNew and hGenNew.");
        Printf("New tree contains %lli entries.", tSLNew->GetEntries());

        // Loop over entries in new SL tree
        Int_t nEntriesAnalysed = 0;
        Int_t nEntriesProgress = (Double_t)tSLNew->GetEntries() / 20.;
        Int_t nPercent = 0;


        for(Int_t iEntry = 0; iEntry < tSLNew->GetEntries(); iEntry++){
            tSLNew->GetEntry(iEntry);
            if(TMath::Abs(parent->Rapidity()) < 1.0){
                nEvNew++;
                fPtGenerated = parent->Pt();
                hGenNew->Fill(fPtGenerated);
                tGenNew->Fill();
            }
            // Update progress bar
            if((iEntry+1) % nEntriesProgress == 0){
                nPercent += 5;
                nEntriesAnalysed += nEntriesProgress;
                Printf("[%i%%] %i entries analysed.", nPercent, nEntriesAnalysed);
            }
        }
        Printf("No. of events with |y| < 1.0: %.0f", nEvNew);

        fGenNew->Write("",TObject::kWriteDelete);
    }

    // #############################################################################################

    Printf("*****");
    Printf("Done.");
    Printf("*****");
    Printf("\n");

    return;
}

void CalculateAndPlotRatios()
{
    Printf("*****");
    Printf("Calculating ratios nGen = 2.000.000.");
    Printf("*****");

    // Load tGenOld
    TFile *fGenOld = TFile::Open("Trees/STARlight/IncoherentShape/tGenOld_2000k_6.624.root", "read");
    if(fGenOld) Printf("File %s loaded.", fGenOld->GetName());

    TTree *tGenOld = (TTree*) fGenOld->Get("tGenOld");
    if(tGenOld) Printf("Tree %s loaded.", tGenOld->GetName()); 
    tGenOld->SetBranchAddress("fPtGen", &fPtGenerated);

    // Load tGenNew
    TFile *fGenNew = TFile::Open("Trees/STARlight/IncoherentShape/tGenNew_2000k_7.350.root", "read");
    if(fGenNew) Printf("File %s loaded.", fGenNew->GetName());

    TTree *tGenNew = (TTree*) fGenNew->Get("tGenNew");
    if(tGenNew) Printf("Tree %s loaded.", tGenNew->GetName()); 
    tGenNew->SetBranchAddress("fPtGen", &fPtGenerated);

    // Load histograms with generated events
    TH1D *hGenOld = (TH1D*)fGenOld->Get("hGenOld");
    if(hGenOld) Printf("Histogram %s loaded.", hGenOld->GetName());
    TH1D *hGenNew = (TH1D*)fGenNew->Get("hGenNew");
    if(hGenNew) Printf("Histogram %s loaded.", hGenNew->GetName());

    hRatios = (TH1D*)hGenNew->Clone("hRatios");
    hRatios->SetTitle("hRatios");
    hRatios->Sumw2();
    hRatios->Divide(hGenOld);

    hRecNew = (TH1D*)hRecOld->Clone("hRecNew");
    hRecNew->SetTitle("hRecNew");
    hRecNew->Multiply(hRatios);

    // Print the results to console
    Printf("\n");
    Printf("*****");

    Printf("Output:");
    Printf("pT_low\tpT_upp\tnEvOld\tnEvNew\tRatio");
    for(Int_t iBin = 1; iBin <= nBins; iBin++){
        Printf("%.3f\t%.3f\t%.0f\t%.0f\t%.3f",
            hRatios->GetBinLowEdge(iBin), hRatios->GetBinLowEdge(iBin+1), 
            hGenOld->GetBinContent(iBin), hGenNew->GetBinContent(iBin), hRatios->GetBinContent(iBin));
    }

    Printf("*****");
    Printf("\n");  

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
    cRatios->Print("Results/STARlight/IncoherentShape/Ratios/coh_modRA_7.350.pdf");
    cRatios->Print("Results/STARlight/IncoherentShape/Ratios/coh_modRA_7.350.png");
    cRec->Print("Results/STARlight/IncoherentShape/RecSpectra/coh_modRA_7.350.pdf");
    cRec->Print("Results/STARlight/IncoherentShape/RecSpectra/coh_modRA_7.350.png");

    delete hGenOld;
    delete hGenNew;

    return;
}

void FillHistRec()
{
    Printf("*****");
    Printf("Calculating histogram of original reconstructed events.");
    Printf("*****");

    TFile *fSL = NULL;
    fSL = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohJpsiToMu.root", "read");
    if(fSL) Printf("File %s loaded.", fSL->GetName());

    // Get the MCRec tree
    TTree *tRec = dynamic_cast<TTree*> (fSL->Get("AnalysisOutput/fTreeJPsiMCRec"));
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
    Printf("\n");

    return;
}