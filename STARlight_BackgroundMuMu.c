// STARlight_BackgroundMuMu.c
// David Grund, Nov 21, 2021

// root headers
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
// my headers
#include "AnalysisManager.h"
#include "STARlight_Utilities.h"

TString str_folder = "Trees/STARlight/bkg_sudakov_3/";
Int_t nBins = 200;
Double_t pT_low = 0.0;
Double_t pT_upp = 1.0;

TH1D *hRec = new TH1D("hRec","hRec",nBins,pT_low,pT_upp);
TH1D *hGen = new TH1D("hGen","hGen",nBins,pT_low,pT_upp);
TH1D *hGenNew = new TH1D("hGen","hGen",nBins,pT_low,pT_upp);

void CompareBkgShapesRecGen();
void NewBkgShape();

void STARlight_BackgroundMuMu(){

    CompareBkgShapesRecGen();

    return;
}

void CompareBkgShapesRecGen(){

    TFile *file = NULL;
    file = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kTwoGammaToMuMedium.root", "read");

    TTree *tRec = dynamic_cast<TTree*> (file->Get("AnalysisOutput/fTreeJPsiMCRec"));
    if(tRec) Printf("MC rec tree loaded.");
    
    ConnectTreeVariablesMCRec(tRec);

    TTree *tGen = dynamic_cast<TTree*> (file->Get("AnalysisOutput/fTreeJPsiMCGen"));
    if(tGen) Printf("MC rec tree loaded.");
    
    ConnectTreeVariablesMCGen(tGen);

    Printf("Tree %s has %lli entries.", tRec->GetName(), tRec->GetEntries());
    Printf("Tree %s has %lli entries.", tGen->GetName(), tGen->GetEntries());

    // Loop over tree entries
    Int_t nEntriesAnalysed = 0;
    Int_t nEvPassed = 0;

    for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++){
        tRec->GetEntry(iEntry);

        // iMassCut == 1 => 3.0 < m < 3.2 GeV
        // iPtCut == 2 => pt < 2.0 GeV
        if(EventPassedMCRec(1, 2)){
            nEvPassed++;
            hRec->Fill(fPt);
        } 

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }
    Printf("fTreeJPsiMc tree Done.");
    Printf("%i events passed the selections.", nEvPassed);

    nEntriesAnalysed = 0;
    nEvPassed = 0;

    for(Int_t iEntry = 0; iEntry < tGen->GetEntries(); iEntry++){
        tGen->GetEntry(iEntry);

        // Rapidity |y| < 0.8
        // Mass 3.0 < m < 3.2 GeV
        //if(EventPassedMCGen() && fMGen > 3.0 && fMGen < 3.2){
        if(EventPassedMCGen()){
            nEvPassed++;
            hGen->Fill(fPtGen);
        } 

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }
    Printf("fTreeJPsiMCGen tree Done.");
    Printf("%i events passed the selections.", nEvPassed);

    // Normalize boths histograms to unity
    hRec->Scale(1./hRec->Integral());
    hGen->Scale(1./hGen->Integral());

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");

    TCanvas *c = new TCanvas("c", "c", 900, 600);
    c->SetTopMargin(0.04);
    c->SetBottomMargin(0.15);
    c->SetRightMargin(0.04);
    c->SetLeftMargin(0.1);
    c->SetLogy();
    // Vertical axis
    hGen->GetYaxis()->SetTitle("#it{N}_{rec} or #it{N}_{gen}");
    hGen->GetYaxis()->SetTitleSize(0.05);
    hGen->GetYaxis()->SetTitleOffset(0.95);
    hGen->GetYaxis()->SetLabelSize(0.05);
    // Horizontal axis
    hGen->GetXaxis()->SetTitle("#it{p}_{T}^{gen} (GeV/#it{c})");
    hGen->GetXaxis()->SetTitleSize(0.05);
    hGen->GetXaxis()->SetTitleOffset(1.2);
    hGen->GetXaxis()->SetLabelSize(0.05);
    // Draw
    hGen->SetLineColor(kRed);
    hGen->Draw("E1"); 
    //hRec->SetLineColor(kBlue);
    //hRec->Draw("E1 SAME");

    // New bkg shape (with Sudakov corrections)
    NewBkgShape();
    hGenNew->Scale(1./hGenNew->Integral());
    hGenNew->SetLineColor(kGreen);
    hGenNew->Draw("E1 SAME");     

    return;
}

void NewBkgShape(){

    TFile *fSL = TFile::Open((str_folder + "trees_starlight.root").Data(), "read");
    if(fSL) Printf("Input file SL loaded.");

    TTree *tSL = dynamic_cast<TTree*> (fSL->Get("starlightTree"));
    if(tSL) Printf("Tree %s loaded.", tSL->GetName());
    ConnectTreeVariables_tSL(tSL);

    for(Int_t iEntry = 0; iEntry < tSL->GetEntries(); iEntry++){
        tSL->GetEntry(iEntry);

        //if(parent->Rapidity() < 0.8 && parent->M() > 3.0 && parent->M() < 3.2){
        if(parent->Rapidity() < 0.8){
            hGenNew->Fill(parent->Pt());
        }
    }

    return;
}