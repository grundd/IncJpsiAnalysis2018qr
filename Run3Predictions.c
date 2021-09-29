// Run3Predictions.c
// David Grund, Sep 29, 2021

// root headers
#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
// my headers
#include "AnalysisManager.h"

Double_t N_Run3_tot = 1100000;
Double_t sig_SL_j_inc = 5.247; //mb
Double_t sig_SL_j_coh = 12.504;//mb

void PrepareData();

void Run3Predictions(){

    //PrepareData();

    gStyle->SetOptTitle(0); // suppress title
    gStyle->SetOptStat(0);  // the type of information printed in the histogram statistics box
                            // 0 = no information
    gStyle->SetPalette(1);  // set color map
    gStyle->SetPaintTextFormat("4.2f"); // precision if plotted with "TEXT"

    TFile *file = TFile::Open("Trees/Run3Predictions/tPtNew.root", "read");
    if(file) Printf("File %s loaded.", file->GetName());

    TList *list = (TList*) file->Get("HistList");
    if(list) Printf("List %s loaded.", list->GetName()); 

    TH1D *hAll = (TH1D*)list->FindObject("hAll");
    if(hAll) Printf("Histogram %s loaded.", hAll->GetName()); 

    TCanvas *c = new TCanvas("c", "c", 900, 600);
    c->SetLogy();
    c->SetTopMargin(0.04);
    c->SetRightMargin(0.04);
    c->SetLeftMargin(0.10);
    hAll->Draw();

    return;
}

void PrepareData(){

    TFile *fCoh = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohJpsiToMu.root", "read");
    if(fCoh) Printf("File %s loaded.", fCoh->GetName());

    TFile *fInc = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohJpsiToMu.root", "read");
    if(fInc) Printf("File %s loaded.", fInc->GetName());

    TTree *tCoh = dynamic_cast<TTree*> (fCoh->Get("AnalysisOutput/fTreeJPsiMCRec"));
    if(tCoh) Printf("tCoh loaded.");

    TTree *tInc = dynamic_cast<TTree*> (fInc->Get("AnalysisOutput/fTreeJPsiMCRec"));
    if(tInc) Printf("tInc loaded.");

    Double_t ratio = sig_SL_j_coh / (sig_SL_j_coh + sig_SL_j_inc);
    Double_t N_coh = N_Run3_tot * ratio;
    Printf("N_coh = %.0f", N_coh);
    Double_t N_inc = N_Run3_tot * (1 - ratio);
    Printf("N_inc = %.0f", N_inc);

    TList *l = new TList();

    TH1D *hCoh = new TH1D("hCoh","hCoh",200,0.0,2.0); 
    TH1D *hInc = new TH1D("hInc","hInc",200,0.0,2.0); 
    
    ConnectTreeVariablesMCRec(tCoh);

    for(Int_t iEntry = 0; iEntry < tCoh->GetEntries(); iEntry++){
        tCoh->GetEntry(iEntry);
        hCoh->Fill(fPt);
    }
  
    ConnectTreeVariablesMCRec(tInc);

    for(Int_t iEntry = 0; iEntry < tInc->GetEntries(); iEntry++){
        tInc->GetEntry(iEntry);
        hInc->Fill(fPt);
    }

    Double_t factorCoh = N_coh / hCoh->GetEntries();
    Double_t factorInc = N_inc / hInc->GetEntries();
    hCoh->Scale(factorCoh);
    hInc->Scale(factorInc);

    TH1D *hAll = (TH1D*)hCoh->Clone("hAll");
    hAll->Add(hInc);

    l->Add(hCoh);
    l->Add(hInc);
    l->Add(hAll);

    Printf("Integral of hCoh: %.0f", hCoh->Integral());
    Printf("Integral of hInc: %.0f", hInc->Integral());
    Printf("Integral of hAll: %.0f", hAll->Integral());

    TFile *f = new TFile("Trees/Run3Predictions/tPtNew.root","RECREATE");
    l->Write("HistList", TObject::kSingleKey);
    f->ls();

    return;
}