// Run3Predictions.c
// David Grund, Oct 6, 2021

// root headers
#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
// my headers
#include "AnalysisManager.h"

Double_t nEvRun2 = 2836.5;  // events Roman had, with pt < 0.11
Double_t nEvRun3 = 1100000; // with arbitraty pt
Double_t SigRun2_val[6] = {1290., 1035., 743., 465., 229., 51.}; 
Double_t SigRun2_err[6] = {108., 77., 50., 37., 20., 6.}
Double_t ptBoundariesRoman[9] = {0.000, 
                                 0.027, 
                                 0.040, 
                                 0.051, 
                                 0.063, 
                                 0.079, 
                                 0.110, // Roman binning ends here
                                 0.155,
                                 0.200
                                 };
Double_t tBoundariesRoman[7] = {0.00000, 
                                0.00072, 
                                0.0016, 
                                0.0026, 
                                0.0040, 
                                0.0062, 
                                0.0121
                                };
Double_t LumiRun3 = 13; // 1/nanobarn

void PrepareData();

void Run3Predictions(){

    //PrepareData();

    TFile *file_in = TFile::Open("Trees/Run3Predictions/HistPtCoh.root", "read");
    if(file_in) Printf("Input data loaded.");

    TList *list_in = dynamic_cast<TList*> (file_in->Get("HistList"));
    if(list_in) Printf("Input list loaded.");

    TH1D *hCoh = (TH1D*)list_in->FindObject("hCoh");
    if(hCoh) Printf("Input histogram loaded.");

    Int_t nBins = 60;

    TH1D *hCohRun2 = new TH1D("hCohRun2", "hCohRun2", nBins, 0.0, 0.3);
    TH1D *hCohRun3 = new TH1D("hCohRun3", "hCohRun3", nBins, 0.0, 0.3);

    Double_t fPtGenerated = 0;

    // Generate Run 2 events
    Int_t nEvRun2Generated = 0;
    while(nEvRun2Generated < nEvRun2){
        fPtGenerated = hCoh->GetRandom();
        hCohRun2->Fill(fPtGenerated);
        if(fPtGenerated < 0.11) nEvRun2Generated++;
    }
    Printf("Run 2 events generated.");
    Double_t nEv = hCohRun2->Integral(1, 22);
    Double_t nEvAll = hCohRun2->Integral(1, nBins+1);
    Printf("Histogram %s contains %.0f events with pt < 0.11.", hCohRun2->GetName(), nEv);
    Printf("Histogram %s contains %.0f events in total.", hCohRun2->GetName(), nEvAll);

    // Generate Run 3 events
    Int_t nEvRun3Generated = 0;
    while(nEvRun3Generated < nEvRun3){
        fPtGenerated = hCoh->GetRandom();
        hCohRun3->Fill(fPtGenerated);        
        nEvRun3Generated++;
    }
    Printf("Run 3 events generated.");
    nEv = hCohRun3->Integral(1, 22);
    nEvAll = hCohRun3->Integral(1, nBins+1);
    Printf("Histogram %s contains %.0f events with pt < 0.11.", hCohRun3->GetName(), nEv);
    Printf("Histogram %s contains %.0f events in total.", hCohRun3->GetName(), nEvAll);

    TCanvas *c = new TCanvas("c", "c", 900, 600);
    c->SetLogy();
    // Canvas margins
    c->SetTopMargin(0.03);
    c->SetBottomMargin(0.14);
    c->SetRightMargin(0.03);
    c->SetLeftMargin(0.10);
    // TStyle
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");
    // hCohRun3
    // Style
    hCohRun3->SetLineColor(225);
    hCohRun3->SetLineWidth(2.0);
    hCohRun3->SetFillColor(225);
    hCohRun3->SetFillStyle(3357);
    // Vertical axis
    hCohRun3->GetYaxis()->SetTitle("Counts per 5 MeV");
    hCohRun3->GetYaxis()->SetTitleSize(0.05);
    hCohRun3->GetYaxis()->SetTitleOffset(1.0);
    hCohRun3->GetYaxis()->SetLabelSize(0.05);
    hCohRun3->GetYaxis()->SetRangeUser(1.0,1e6);
    //hAxE->GetYaxis()->SetDecimals(3);
    // Horizontal axis
    hCohRun3->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hCohRun3->GetXaxis()->SetTitleSize(0.05);
    hCohRun3->GetXaxis()->SetTitleOffset(1.2);
    hCohRun3->GetXaxis()->SetLabelSize(0.05);
    hCohRun3->GetXaxis()->SetLabelOffset(0.015);
    hCohRun3->GetXaxis()->SetDecimals(2);
    // Draw it
    hCohRun3->Draw("HIST");
    // hCohRun2
    // Line style
    hCohRun2->SetLineColor(221);
    hCohRun2->SetLineWidth(2.0);
    hCohRun2->SetFillColor(221);
    hCohRun2->SetFillStyle(3114);
    // Draw it
    hCohRun2->Draw("HIST SAME");
    // Legend 1
    TLegend *l1 = new TLegend(0.1,0.90,0.4,0.95);
    l1->AddEntry((TObject*)0,"#it{p}_{T}-distribution of coherently photoproduced J/#psi measured with ALICE","");
    l1->SetTextSize(0.042);
    l1->SetBorderSize(0);
    l1->SetFillStyle(0);
    l1->Draw();
    // Legend 2
    TLegend *l2 = new TLegend(0.54,0.72,0.76,0.90);
    l2->AddEntry(hCohRun2,"Run 2 (STARlight simulations)","f");
    l2->AddEntry(hCohRun3,"Run 3 predictions (STARlight sim)","f");
    l2->SetTextSize(0.042);
    l2->SetBorderSize(0);
    l2->SetFillStyle(0);
    l2->Draw();

    c->Print("Results/Run3Predictions/CohJpsi_Run3Predictions.pdf");
    c->Print("Results/Run3Predictions/CohJpsi_Run3Predictions.png");

    return;
}

void PrepareData(){

    TFile *fCoh = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohJpsiToMu.root", "read");
    if(fCoh) Printf("File %s loaded.", fCoh->GetName());

    TTree *tCoh = dynamic_cast<TTree*> (fCoh->Get("AnalysisOutput/fTreeJPsiMCRec"));
    if(tCoh) Printf("tCoh loaded.");

    ConnectTreeVariablesMCRec(tCoh);

    TList *l = new TList();
    TH1D *hCoh = new TH1D("hCoh","hCoh",200,0.0,0.5); 
    l->Add(hCoh);

    for(Int_t iEntry = 0; iEntry < tCoh->GetEntries(); iEntry++){
        tCoh->GetEntry(iEntry);
        hCoh->Fill(fPt);
    }
    Printf("Done.");

    TFile *f = new TFile("Trees/Run3Predictions/HistPtCoh.root","RECREATE");
    l->Write("HistList", TObject::kSingleKey);
    f->ls();

    return;
}