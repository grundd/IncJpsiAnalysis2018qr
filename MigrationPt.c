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

void PlotHistMigration(Bool_t vsPt);

void MigrationPt(){

    PlotHistMigration(kTRUE);

    //PlotHistMigration(kFALSE);

    return;
}

void PlotHistMigration(Bool_t vsPt){

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
    Double_t pt2Boundaries[nPtBins+1] = { 0 };
    for(Int_t i = 0; i <= nPtBins; i++) pt2Boundaries[i] = ptBoundaries[i]*ptBoundaries[i];
    TH2D *hMigration = NULL;
    TH1D *hScaleByTotalNRec = NULL;
    if(vsPt){
        hMigration = new TH2D("hMigration", "pt rec vs pt gen", nPtBins, ptBoundaries, nPtBins, ptBoundaries);
        hScaleByTotalNRec = new TH1D("hScaleByTotalNRec", "Total NRec per bin in pt gen", nPtBins, ptBoundaries);

    } else {
        hMigration = new TH2D("hMigration", "pt^{2} rec vs pt^{2} gen", nPtBins, pt2Boundaries, nPtBins, pt2Boundaries);
        hScaleByTotalNRec = new TH1D("hScaleByTotalNRec", "Total NRec per bin in pt gen", nPtBins, pt2Boundaries);
    } 
    // Fill the histogram
    for(Int_t iEntry = 0; iEntry < fTreeIn->GetEntries(); iEntry++){
        fTreeIn->GetEntry(iEntry);

        hMigration->Fill(fPtGen,fPt);
        hScaleByTotalNRec->Fill(fPtGen);

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }
    // Scale the histogram
    for(Int_t iBinX = 1; iBinX <= nPtBins; iBinX++){
        Double_t NRec_tot = hScaleByTotalNRec->GetBinContent(iBinX);
        Printf("BinX %i: %.0f", iBinX, NRec_tot);
        for(Int_t iBinY = 1; iBinY <= nPtBins; iBinY++){
            Double_t ValueScaled = hMigration->GetBinContent(iBinX,iBinY) / NRec_tot;
            hMigration->SetBinContent(iBinX,iBinY,ValueScaled);
        }
    }

    Printf("Number of entries in the histogram: %.0f", hMigration->GetEntries());

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.3f");

    TCanvas *c = new TCanvas("canvas","canvas",900,600);
    c->SetTopMargin(0.03);
    c->SetBottomMargin(0.13);
    c->SetRightMargin(0.12);

    hMigration->SetMarkerSize(1.8);
    // horizontal axis
    if(vsPt) hMigration->GetXaxis()->SetTitle("#it{p}_{T}^{gen} (GeV/#it{c})");
    else     hMigration->GetXaxis()->SetTitle("#it{p}^{2}_{T,gen} (GeV^{2}/#it{c}^{2})");
    hMigration->GetXaxis()->SetTitleSize(0.05);
    hMigration->GetXaxis()->SetLabelSize(0.05);
    hMigration->GetXaxis()->SetTitleOffset(1.1);
    // vertical axis
    if(vsPt) hMigration->GetYaxis()->SetTitle("#it{p}_{T}^{rec} (GeV/#it{c})");
    else     hMigration->GetYaxis()->SetTitle("#it{p}^{2}_{T,rec} (GeV^{2}/#it{c}^{2})");
    hMigration->GetYaxis()->SetTitleSize(0.05);
    hMigration->GetYaxis()->SetLabelSize(0.05);
    hMigration->GetYaxis()->SetTitleOffset(0.9);
    // Z-axis
    hMigration->GetZaxis()->SetLabelSize(0.05);
    hMigration->Draw("COLZ TEXT");

    TString *str = NULL;
    if(vsPt) str = new TString(Form("%sMigrPt_pt_%ibins", Path.Data(), nPtBins));
    else     str = new TString(Form("%sMigrPt_pt2_%ibins", Path.Data(), nPtBins));
    c->Print((*str + ".pdf").Data());
    c->Print((*str + ".png").Data());

    return;
}