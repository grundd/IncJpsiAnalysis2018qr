// MigrationPt.c
// David Grund, Sep 28, 2021

// root headers
#include "TFile.h"
#include "TH2.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
// my headers
#include "AnalysisManager.h"

TString Path = "Results/MigrationPt/";

void PlotHistMigration(Bool_t vsPt, Bool_t inPercent);

void MigrationPt(){

    PlotHistMigration(kTRUE, kTRUE);

    PlotHistMigration(kFALSE, kTRUE);

    return;
}

void PlotHistMigration(Bool_t vsPt, Bool_t inPercent){

    // Load data
    TFile *fFileIn = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohJpsiToMu.root", "read");
    if(fFileIn) Printf("Input data loaded.");

    TTree *fTreeIn = dynamic_cast<TTree*> (fFileIn->Get("AnalysisOutput/fTreeJPsiMCRec"));
    if(fTreeIn) Printf("Input tree loaded.");

    ConnectTreeVariablesMCRec(fTreeIn);

    Printf("%lli entries found in the tree.", fTreeIn->GetEntries());
    Int_t nEntriesAnalysed = 0;

    // Create 2D histogram
    SetPtBinning();

    // pt gen on the horizontal axis, pt rec on the vertical axis
    Double_t ptBoundariesNew[nPtBins+3] = { 0 };
    Double_t pt2BoundariesNew[nPtBins+3] = { 0 };
    ptBoundariesNew[0] = 0.0;
    for(Int_t i = 1; i <= nPtBins+1; i++) ptBoundariesNew[i] = ptBoundaries[i-1];
    ptBoundariesNew[nPtBins+2] = 1.2;
    for(Int_t i = 0; i <= nPtBins+2; i++) pt2BoundariesNew[i] = ptBoundariesNew[i]*ptBoundariesNew[i];
    TH2D *hMigration = NULL;
    TH1D *hScaleByTotalNRec = NULL;
    if(vsPt){
        hMigration = new TH2D("hMigration", "pt rec vs pt gen", nPtBins+2, ptBoundariesNew, nPtBins+2, ptBoundariesNew);
        hScaleByTotalNRec = new TH1D("hScaleByTotalNRec", "Total NRec per bin in pt gen", nPtBins+2, ptBoundariesNew);

    } else {
        hMigration = new TH2D("hMigration", "pt^{2} rec vs pt^{2} gen", nPtBins+2, pt2BoundariesNew, nPtBins+2, pt2BoundariesNew);
        hScaleByTotalNRec = new TH1D("hScaleByTotalNRec", "Total NRec per bin in pt gen", nPtBins+2, pt2BoundariesNew);
    } 
    // Fill the histogram
    for(Int_t iEntry = 0; iEntry < fTreeIn->GetEntries(); iEntry++){
        fTreeIn->GetEntry(iEntry);

        if(EventPassed(0,-1)){
            if(vsPt){
                hMigration->Fill(fPtGen,fPt);
                hScaleByTotalNRec->Fill(fPtGen);
            } else {
                hMigration->Fill(fPtGen*fPtGen,fPt*fPt);
                hScaleByTotalNRec->Fill(fPtGen*fPtGen);
            }
        }
        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }
    // Scale the histogram
    for(Int_t iBinX = 1; iBinX <= nPtBins+2; iBinX++){
        Double_t NRec_tot = hScaleByTotalNRec->GetBinContent(iBinX);
        Printf("BinX %i: %.0f", iBinX, NRec_tot);
        for(Int_t iBinY = 1; iBinY <= nPtBins+2; iBinY++){
            Double_t ValueScaled;
            if(inPercent) ValueScaled = hMigration->GetBinContent(iBinX,iBinY) / NRec_tot * 100;
            else ValueScaled = hMigration->GetBinContent(iBinX,iBinY) / NRec_tot;
            hMigration->SetBinContent(iBinX,iBinY,ValueScaled);
        }
    }

    Printf("Number of entries in the histogram: %.0f", hMigration->GetEntries());

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    if(inPercent) gStyle->SetPaintTextFormat("4.2f");
    else gStyle->SetPaintTextFormat("4.3f");

    TCanvas *c = new TCanvas("c","c",1200,600);
    if(!vsPt){
        c->SetLogx();
        c->SetLogy();
    }
    c->SetTopMargin(0.03);
    c->SetBottomMargin(0.145);
    c->SetRightMargin(0.1);
    c->SetLeftMargin(0.075);

    hMigration->SetMarkerSize(1.8);
    // horizontal axis
    if(vsPt){
        hMigration->GetXaxis()->SetTitle("#it{p}_{T}^{gen} (GeV/#it{c})");
        hMigration->GetXaxis()->SetLabelOffset(0.015);
    } else {
        hMigration->GetXaxis()->SetTitle("(#it{p}_{T}^{gen})^{2} (GeV^{2}/#it{c}^{2})");
    }   
    hMigration->GetXaxis()->SetTitleSize(0.05);
    hMigration->GetXaxis()->SetTitleOffset(1.3);
    hMigration->GetXaxis()->SetLabelSize(0.05);
    hMigration->GetXaxis()->SetDecimals(1);
    // vertical axis
    if(vsPt) hMigration->GetYaxis()->SetTitle("#it{p}_{T}^{rec} (GeV/#it{c})");
    else     hMigration->GetYaxis()->SetTitle("(#it{p}_{T}^{rec})^{2} (GeV^{2}/#it{c}^{2})");
    hMigration->GetYaxis()->SetTitleSize(0.05);
    hMigration->GetYaxis()->SetLabelSize(0.05);
    hMigration->GetYaxis()->SetTitleOffset(0.65);
    hMigration->GetYaxis()->SetDecimals(1);
    // Z-axis
    hMigration->GetZaxis()->SetLabelSize(0.05);
    hMigration->Draw("COLZ TEXT");
    // Set ranges
    hMigration->GetXaxis()->SetRangeUser(0.0,1.2);
    hMigration->GetYaxis()->SetRangeUser(0.0,1.2);

    TString *str = NULL;
    if(vsPt) str = new TString(Form("%sMigrPt_%ibins_pt", Path.Data(), nPtBins));
    else     str = new TString(Form("%sMigrPt_%ibins_pt2", Path.Data(), nPtBins));
    c->Print((*str + ".pdf").Data());
    c->Print((*str + ".png").Data());

    return;
}