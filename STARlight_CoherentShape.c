// STARlight_CoherentShape.c
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

Double_t nEv_old1 = 0;
Double_t nEv_new1 = 0;
Double_t nEv_old2 = 0;
Double_t nEv_new2 = 0;

Double_t pT_low = 0.0; // GeV/c
Double_t pT_upp = 0.4; // GeV/c
Int_t nBins1 = 40;
Int_t nBins2 = 40;
TH1D *hGenOld1 = new TH1D("hGenOld1","hGenOld1",nBins1,pT_low,pT_upp);
TH1D *hGenNew1 = new TH1D("hGenNew1","hGenNew1",nBins1,pT_low,pT_upp);
TH1D *hRatios1 = NULL;
TH1D *hGenOld2 = new TH1D("hGenOld2","hGenOld2",nBins2,pT_low,pT_upp);
TH1D *hGenNew2 = new TH1D("hGenNew2","hGenNew2",nBins2,pT_low,pT_upp);
TH1D *hRatios2 = NULL;
TH1D *hRecOld1 = new TH1D("hRecOld","hRecOld",nBins1,pT_low,pT_upp);
TH1D *hRecOld2 = new TH1D("hRecOld","hRecOld",nBins2,pT_low,pT_upp);
TH1D *hRecNew1 = NULL;
TH1D *hRecNew2 = NULL;

void CalculateShapes();

void STARlight_CoherentShape(){

    CalculateShapes();

    return;
}

void CalculateShapes(){

    // 1) ##########################################################################################
    // NGen = 4.880.000
    // 1a) Load official root file: kCohJpsiToMu
    TFile *fSL_old1 = NULL;
    fSL_old1 = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohJpsiToMu.root", "read");
    if(fSL_old1) Printf("File %s loaded.", fSL_old1->GetName());

    // Get the MCGen tree
    TTree *tGen = dynamic_cast<TTree*> (fSL_old1->Get("AnalysisOutput/fTreeJPsiMCGen"));
    if(tGen) Printf("Tree %s loaded.", tGen->GetName());
    ConnectTreeVariablesMCGen(tGen);

    // Get the MCRec tree
    TTree *tRec = dynamic_cast<TTree*> (fSL_old1->Get("AnalysisOutput/fTreeJPsiMCRec"));
    if(tRec) Printf("Tree %s loaded.", tRec->GetName());
    ConnectTreeVariablesMCRec(tRec);

    // 1b) Load generated coherent MC events (with R_A = 7.53 fm)
    TFile *fSL_new1 = TFile::Open("Trees/STARlight/coh_modRA/trees_starlight.root", "read");
    if(fSL_new1) Printf("File %s loaded.", fSL_new1->GetName());

    // Get the SL tree
    TTree *tSL_new1 = dynamic_cast<TTree*> (fSL_new1->Get("starlightTree"));
    if(tSL_new1) Printf("Tree %s loaded.", tSL_new1->GetName());
    ConnectTreeVariables_tSL(tSL_new1);

    // 2) ##########################################################################################
    // NGen = 6.000.000
    // 2a) Load generated coherent MC events (with R_A = 6.62 fm)
    TFile *fSL_old2 = TFile::Open("Trees/STARlight/coh_stdRA/trees_starlight.root", "read");
    if(fSL_old2) Printf("File %s loaded.", fSL_old2->GetName());

    // Get the SL tree
    TTree *tSL_old2 = dynamic_cast<TTree*> (fSL_old2->Get("starlightTree"));
    if(tSL_old2) Printf("Tree %s loaded.", tSL_old2->GetName());
    ConnectTreeVariables_tSL(tSL_old2);  

    // 2b) Load generated coherent MC events (with R_A = 7.53 fm)
    TFile *fSL_new2 = TFile::Open("Trees/STARlight/coh_modRA_2/trees_starlight.root", "read");
    if(fSL_new2) Printf("File %s loaded.", fSL_new2->GetName());

    // Get the SL tree
    TTree *tSL_new2 = dynamic_cast<TTree*> (fSL_new2->Get("starlightTree"));
    if(tSL_new2) Printf("Tree %s loaded.", tSL_new2->GetName());
    ConnectTreeVariables_tSL(tSL_new2);      

    // #############################################################################################
    // Define the output
    TFile *f = new TFile("Trees/STARlight/CoherentShape/trees_CoherentShape.root","RECREATE");
    TList *l = new TList();

    TTree *tGenOld1 = new TTree("tGenOld1","tGenOld1");
    tGenOld1->Branch("fPtGen",&fPtGenerated,"fPtGen/D");
    l->Add(tGenOld1);
    l->Add(hGenOld1);
    TTree *tGenOld2 = new TTree("tGenOld2","tGenOld2");
    tGenOld2->Branch("fPtGen",&fPtGenerated,"fPtGen/D");
    l->Add(tGenOld2);
    l->Add(hGenOld2);
    TTree *tGenNew1 = new TTree("tGenNew1","tGenNew1");
    tGenNew1->Branch("fPtGen",&fPtGenerated,"fPtGen/D");
    l->Add(tGenNew1);
    l->Add(hGenNew1);
    TTree *tGenNew2 = new TTree("tGenNew2","tGenNew2");
    tGenNew2->Branch("fPtGen",&fPtGenerated,"fPtGen/D");
    l->Add(tGenNew2);
    l->Add(hGenNew2);

    // #############################################################################################
    Printf("\n***");
    Printf("Analysis 1: nGen = 4.880.000");
    Printf("kCohJpsiToMu MCGen tree contains %lli entries.", tGen->GetEntries());
    for(Int_t iEntry = 0; iEntry < tGen->GetEntries(); iEntry++){
        tGen->GetEntry(iEntry);
        if(TMath::Abs(fYGen) < 1.0){
            nEv_old1++;
            fPtGenerated = fPtGen;
            hGenOld1->Fill(fPtGenerated);
            tGenOld1->Fill();
        }
    }
    Printf("No. of events with |y| < 1.0: %.0f", nEv_old1);
    Printf("STARlight tree (7.53 fm) contains %lli entries.", tSL_new1->GetEntries());
    for(Int_t iEntry = 0; iEntry < tSL_new1->GetEntries(); iEntry++){
        tSL_new1->GetEntry(iEntry);
        if(TMath::Abs(parent->Rapidity()) < 1.0){
            nEv_new1++;
            fPtGenerated = parent->Pt();
            hGenNew1->Fill(fPtGenerated);
            tGenNew1->Fill();
        }
    }
    Printf("No. of events with |y| < 1.0: %.0f", nEv_new1);
    // #############################################################################################
    Printf("\n***");
    Printf("Analysis 2: nGen = 6.000.000");
    Printf("STARlight tree (6.62 fm) contains %lli entries.", tSL_old2->GetEntries());
    for(Int_t iEntry = 0; iEntry < tSL_old2->GetEntries(); iEntry++){
        tSL_old2->GetEntry(iEntry);
        if(TMath::Abs(parent->Rapidity()) < 1.0){
            nEv_old2++;
            fPtGenerated = parent->Pt();
            hGenOld2->Fill(fPtGenerated);
            tGenOld2->Fill();
        }
    }
    Printf("No. of events with |y| < 1.0: %.0f", nEv_old2);
    Printf("STARlight tree (7.53 fm) contains %lli entries.", tSL_new2->GetEntries());
    for(Int_t iEntry = 0; iEntry < tSL_new2->GetEntries(); iEntry++){
        tSL_new2->GetEntry(iEntry);
        if(TMath::Abs(parent->Rapidity()) < 1.0){
            nEv_new2++;
            fPtGenerated = parent->Pt();
            hGenNew2->Fill(fPtGenerated);
            tGenNew2->Fill();
        }
    }
    Printf("No. of events with |y| < 1.0: %.0f", nEv_new2);
    // #############################################################################################

    // Prepare reconstructed events
    for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++){
        tRec->GetEntry(iEntry);
        if(EventPassedMCRec(1, 2)){
            hRecOld1->Fill(fPt);
            hRecOld2->Fill(fPt);
        } 
    }

    hRatios1 = (TH1D*)hGenNew1->Clone("hRatios1");
    hRatios1->SetTitle("hRatios1");
    hRatios1->Sumw2();
    hRatios1->Divide(hGenOld1);

    hRecNew1 = (TH1D*)hRecOld1->Clone("hRecNew1");
    hRecNew1->SetTitle("hRecNew1");
    hRecNew1->Multiply(hRatios1);

    hRatios2 = (TH1D*)hGenNew2->Clone("hRatios2");
    hRatios2->SetTitle("hRatios2");
    hRatios2->Sumw2();
    hRatios2->Divide(hGenOld2);

    hRecNew2 = (TH1D*)hRecOld2->Clone("hRecNew2");
    hRecNew2->SetTitle("hRecNew2");
    hRecNew2->Multiply(hRatios2);

    l->Add(hRatios1);
    l->Add(hRatios2);

    // Print the results to console
    Printf("\n***");
    Printf("*** Analysis 1: nGen = 4.880.000");    
    Printf("pT_low\tpT_upp\tnEv_old\tnEv_new\tratio\n");
    for(Int_t iBin = 1; iBin <= nBins1; iBin++){
        Printf("%.3f\t%.3f\t%.0f\t%.0f\t%.3f",
            hRatios1->GetBinLowEdge(iBin), hRatios1->GetBinLowEdge(iBin+1), 
            hGenOld1->GetBinContent(iBin), hGenNew1->GetBinContent(iBin), hRatios1->GetBinContent(iBin));
    }
    Printf("\n***");
    Printf("Analysis 2: nGen = 6.000.000");    
    Printf("pT_low\tpT_upp\tnEv_old\tnEv_new\tratio\n");
    for(Int_t iBin = 1; iBin <= nBins2; iBin++){
        Printf("%.3f\t%.3f\t%.0f\t%.0f\t%.3f",
            hRatios2->GetBinLowEdge(iBin), hRatios2->GetBinLowEdge(iBin+1), 
            hGenOld2->GetBinContent(iBin), hGenNew2->GetBinContent(iBin), hRatios2->GetBinContent(iBin));
    }

    // Save the results to the root file
    l->Write("OutputList", TObject::kSingleKey);
    f->ls();

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
    hRatios1->GetYaxis()->SetTitle("#it{N}_{gen}(#it{R}_{A} = 7.53 fm) / #it{N}_{gen}(#it{R}_{A} = 6.62 fm)");
    hRatios1->GetYaxis()->SetTitleSize(0.05);
    hRatios1->GetYaxis()->SetTitleOffset(0.9);
    hRatios1->GetYaxis()->SetLabelSize(0.05);
    hRatios1->GetYaxis()->SetDecimals(1);
    hRatios1->GetYaxis()->SetRangeUser(0.0,5.0);
    // Horizontal axis
    hRatios1->GetXaxis()->SetTitle("#it{p}_{T}^{gen} (GeV/#it{c})");
    hRatios1->GetXaxis()->SetTitleSize(0.05);
    hRatios1->GetXaxis()->SetTitleOffset(1.2);
    hRatios1->GetXaxis()->SetLabelSize(0.05);
    // Draw
    hRatios1->SetLineColor(kRed);
    hRatios1->Draw("E1");
    hRatios2->SetLineColor(215);
    hRatios2->Draw("E1 SAME");

    TCanvas *cRec = new TCanvas("cRec", "cRec", 900, 600);
    cRec->SetTopMargin(0.04);
    cRec->SetBottomMargin(0.15);
    cRec->SetRightMargin(0.04);
    cRec->SetLeftMargin(0.1);
    cRec->SetLogy();
    // Vertical axis
    hRecOld1->GetYaxis()->SetTitle("#it{N}_{rec} (after selections)");
    hRecOld1->GetYaxis()->SetTitleSize(0.05);
    hRecOld1->GetYaxis()->SetTitleOffset(0.9);
    hRecOld1->GetYaxis()->SetLabelSize(0.05);
    hRecOld1->GetYaxis()->SetRangeUser(0.1,1e5);
    // Horizontal axis
    hRecOld1->GetXaxis()->SetTitle("#it{p}_{T}^{gen} (GeV/#it{c})");
    hRecOld1->GetXaxis()->SetTitleSize(0.05);
    hRecOld1->GetXaxis()->SetTitleOffset(1.2);
    hRecOld1->GetXaxis()->SetLabelSize(0.05);
    // Draw
    hRecOld1->SetLineColor(210);
    hRecOld1->Draw("E1");
    hRecNew1->SetLineColor(kRed);
    hRecNew1->Draw("E1 SAME");
    hRecNew2->SetLineColor(215);
    hRecNew2->Draw("E1 SAME");

    // Print plots
    cRatios->Print("Results/STARlight/CoherentShape/ratios.png");
    cRatios->Print("Results/STARlight/CoherentShape/ratios.pdf");
    cRec->Print("Results/STARlight/CoherentShape/reconstructed_spectra.png");
    cRec->Print("Results/STARlight/CoherentShape/reconstructed_spectra.pdf");

    return;
}