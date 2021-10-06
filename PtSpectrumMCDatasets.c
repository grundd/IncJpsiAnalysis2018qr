// PtSpectrumMCDatasets.c
// David Grund, Oct 6, 2021

// root headers
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
// my headers
#include "AnalysisManager.h"

void PlotPtSpectrum(Int_t opt);

void PtSpectrumMCDatasets(){

    PlotPtSpectrum(0);

    PlotPtSpectrum(1);

    return;
}

void PlotPtSpectrum(Int_t opt){

    TH1D *hCoh = new TH1D("hCoh","hCoh",50,0.,1.);
    TH1D *hInc = new TH1D("hInc","hInc",50,0.,1.);

    // Load the data
    TFile *file1 = NULL;
    TFile *file2 = NULL;
    switch(opt){
        case 0:
            file1 = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohJpsiToMu.root", "read");
            file2 = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohJpsiToMu.root", "read");
            break;
        case 1:
            file1 = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohPsi2sToMuPi.root", "read");
            file2 = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohPsi2sToMuPi.root", "read");
            break;
    }
    if(file1) Printf("File %s loaded.", file1->GetName());
    if(file2) Printf("File %s loaded.", file2->GetName());

    TTree *tRec1 = dynamic_cast<TTree*> (file1->Get("AnalysisOutput/fTreeJPsiMCRec"));
    if(tRec1) Printf("MC rec tree loaded.");
    TTree *tRec2 = dynamic_cast<TTree*> (file2->Get("AnalysisOutput/fTreeJPsiMCRec"));
    if(tRec2) Printf("MC rec tree loaded.");

    Printf("Tree %s has %lli entries.", tRec1->GetName(), tRec1->GetEntries());
    Printf("Tree %s has %lli entries.", tRec2->GetName(), tRec2->GetEntries());
    
    ConnectTreeVariablesMCRec(tRec1);    

    // Loop over tree entries
    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < tRec1->GetEntries(); iEntry++){
        tRec1->GetEntry(iEntry);

        // iMassCut == 0 => 2.2 < m < 4.5 GeV, iPtCut == 2 => pt < 2.0 GeV
        if(EventPassedMCRec(1, 2)) hCoh->Fill(fPt);

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }
    Printf("Done.");
    hCoh->Scale(1./hCoh->Integral());

    ConnectTreeVariablesMCRec(tRec2);    

    // Loop over tree entries
    nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < tRec2->GetEntries(); iEntry++){
        tRec2->GetEntry(iEntry);

        // iMassCut == 0 => 2.2 < m < 4.5 GeV, iPtCut == 2 => pt < 2.0 GeV
        if(EventPassedMCRec(1, 2)) hInc->Fill(fPt);
         
        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }
    Printf("Done.");
    hInc->Scale(1./hInc->Integral());

    TCanvas *c = new TCanvas("c","c",900,600);
    c->SetTopMargin(0.02);
    c->SetBottomMargin(0.14);
    c->SetRightMargin(0.03);
    c->SetLeftMargin(0.13);

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");

    hCoh->SetLineColor(215);
    hCoh->SetLineWidth(2); 
    // Vertical axis
    hCoh->GetYaxis()->SetTitle("Probability density");
    hCoh->GetYaxis()->SetTitleSize(0.06);
    hCoh->GetYaxis()->SetTitleOffset(1.08);
    hCoh->GetYaxis()->SetLabelSize(0.06);
    hCoh->GetYaxis()->SetLabelOffset(0.01);
    hCoh->GetYaxis()->SetDecimals(2);
    // Horizontal axis
    hCoh->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hCoh->GetXaxis()->SetTitleSize(0.06);
    hCoh->GetXaxis()->SetLabelSize(0.06);
    //hCoh->GetXaxis()->SetDecimals(1);
    hCoh->Draw("HIST");

    hInc->SetLineColor(kRed);   
    hInc->SetLineWidth(2); 
    hInc->Draw("HIST SAME");

    // Legend1
    TLegend *l1 = new TLegend(0.10,0.9,0.95,0.97);
    l1->AddEntry((TObject*)0,Form("ALICE Simulation, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV"),"");
    l1->SetTextSize(0.058);
    l1->SetBorderSize(0); // no border
    l1->SetFillStyle(0);  // legend is transparent
    l1->Draw();

    // Legend2
    Double_t x_low = 0;
    if(opt == 0) x_low = 0.60;
    else         x_low = 0.44;
    TLegend *l2 = new TLegend(x_low,0.68,x_low + 0.3,0.89);
    l2->AddEntry((TObject*)0,Form("#it{p}_{T} of dimuons from:"),"");
    if(opt == 0){
        l2->AddEntry(hCoh,"coh J/#psi #rightarrow #mu^{+}#mu^{-}");
        l2->AddEntry(hInc,"inc J/#psi #rightarrow #mu^{+}#mu^{-}");
    } else {
        l2->AddEntry(hCoh,"coh #psi(2s) #rightarrow J/#psi + #pi^{+}#pi^{-} #rightarrow #mu^{+}#mu^{-} + #pi^{+}#pi^{-}");
        l2->AddEntry(hInc,"inc #psi(2s) #rightarrow J/#psi + #pi^{+}#pi^{-} #rightarrow #mu^{+}#mu^{-} + #pi^{+}#pi^{-}");
    }
    l2->SetTextSize(0.045);
    l2->SetBorderSize(0); // no border
    l2->SetFillStyle(0);  // legend is transparent
    l2->Draw();  

    if(opt == 0){
        c->Print("Results/PtSpectrumMCDatasets/SpectrumJpsi.pdf");
        c->Print("Results/PtSpectrumMCDatasets/SpectrumJpsi.png");
    } else {
        c->Print("Results/PtSpectrumMCDatasets/SpectrumPsi2s.pdf");
        c->Print("Results/PtSpectrumMCDatasets/SpectrumPsi2s.png");
    }

    return;
}