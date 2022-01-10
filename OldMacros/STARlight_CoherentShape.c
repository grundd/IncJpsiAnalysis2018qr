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

Double_t pT_low = 0.0; // GeV/c
Double_t pT_upp = 0.4; // GeV/c

void CalculateShapes(TString str_in, TString str_out_subfolder, TString str_out_name, Double_t R_A, Bool_t bCompareWith1);

    TString str_in = "";
    TString str_out_subfolder = "";
    TString str_out_name = "";

void STARlight_CoherentShape(){

    /*
    str_in = "coh_6000000_modRA";
    str_out_subfolder = "CoherentShape";
    str_out_name = "7.53";
    CalculateShapes(str_in, str_out_subfolder, str_out_name, 7.53, kTRUE);
    */

    str_in = "OptimalRA/coh_modRA_0_6.624";
    str_out_subfolder = "OptimalRA";
    str_out_name = "0_6.624";
    CalculateShapes(str_in, str_out_subfolder, str_out_name, 6.624, kFALSE);

    Int_t iUpTo = 7;
    Double_t R_A = 0.;
    for(Int_t i = 0; i < iUpTo; i++){
        R_A = 6.60 + i * 0.10;
        str_in = Form("OptimalRA/coh_modRA_%i_%.2f", i+1, R_A);
        str_out_subfolder = "OptimalRA";
        str_out_name = Form("%i_%.2f", i+1, R_A);
        CalculateShapes(str_in, str_out_subfolder, str_out_name, R_A, kFALSE);
    }

    return;
}

void CalculateShapes(TString str_in, TString str_out_subfolder, TString str_out_name, Double_t R_A, Bool_t bCompareWith1){

    // NGen = 4.880.000
    Int_t nBins1 = 40;
    Double_t nEv_old1 = 0;
    Double_t nEv_new1 = 0;
    TH1D *hGenOld1 = new TH1D("hGenOld1","hGenOld1",nBins1,pT_low,pT_upp);
    TH1D *hGenNew1 = new TH1D("hGenNew1","hGenNew1",nBins1,pT_low,pT_upp);
    TH1D *hRatios1 = NULL;
    TH1D *hRecOld1 = new TH1D("hRecOld1","hRecOld1",nBins1,pT_low,pT_upp);
    TH1D *hRecNew1 = NULL;
    // NGen == 6.000.000
    Int_t nBins2 = 40;
    Double_t nEv_old2 = 0;
    Double_t nEv_new2 = 0;
    TH1D *hGenOld2 = new TH1D("hGenOld2","hGenOld2",nBins2,pT_low,pT_upp);
    TH1D *hGenNew2 = new TH1D("hGenNew2","hGenNew2",nBins2,pT_low,pT_upp);
    TH1D *hRatios2 = NULL;
    TH1D *hRecOld2 = new TH1D("hRecOld2","hRecOld2",nBins2,pT_low,pT_upp);
    TH1D *hRecNew2 = NULL;

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
    TFile *fSL_new1 = TFile::Open("Trees/STARlight/coh_4880000_modRA/trees_starlight.root", "read");
    if(fSL_new1) Printf("File %s loaded.", fSL_new1->GetName());

    // Get the SL tree
    TTree *tSL_new1 = dynamic_cast<TTree*> (fSL_new1->Get("starlightTree"));
    if(tSL_new1) Printf("Tree %s loaded.", tSL_new1->GetName());
    ConnectTreeVariables_tSL(tSL_new1);

    // 2) ##########################################################################################
    // NGen = 6.000.000
    // 2a) Load generated coherent MC events (with R_A = 6.62 fm)
    TFile *fSL_old2 = TFile::Open("Trees/STARlight/coh_6000000_stdRA/trees_starlight.root", "read");
    if(fSL_old2) Printf("File %s loaded.", fSL_old2->GetName());

    // Get the SL tree
    TTree *tSL_old2 = dynamic_cast<TTree*> (fSL_old2->Get("starlightTree"));
    if(tSL_old2) Printf("Tree %s loaded.", tSL_old2->GetName());
    ConnectTreeVariables_tSL(tSL_old2);  

    // 2b) Load generated coherent MC events (with different R_A)
    TFile *fSL_new2 = TFile::Open(Form("Trees/STARlight/%s/trees_starlight.root", str_in.Data()), "read");
    if(fSL_new2) Printf("File %s loaded.", fSL_new2->GetName());

    // Get the SL tree
    TTree *tSL_new2 = dynamic_cast<TTree*> (fSL_new2->Get("starlightTree"));
    if(tSL_new2) Printf("Tree %s loaded.", tSL_new2->GetName());
    ConnectTreeVariables_tSL(tSL_new2);      

    // #############################################################################################
    // Define the output
    TFile *f = new TFile(Form("Trees/STARlight/%s/trees_CohShape_%s.root", str_out_subfolder.Data(), str_out_name.Data()),"RECREATE");
    TList *l = new TList();

    TTree *tGenOld1 = new TTree("tGenOld1","tGenOld1");
    tGenOld1->Branch("fPtGen",&fPtGenerated,"fPtGen/D");
    if(bCompareWith1) l->Add(tGenOld1);
    if(bCompareWith1) l->Add(hGenOld1);
    TTree *tGenOld2 = new TTree("tGenOld2","tGenOld2");
    tGenOld2->Branch("fPtGen",&fPtGenerated,"fPtGen/D");
    l->Add(tGenOld2);
    l->Add(hGenOld2);
    TTree *tGenNew1 = new TTree("tGenNew1","tGenNew1");
    tGenNew1->Branch("fPtGen",&fPtGenerated,"fPtGen/D");
    if(bCompareWith1) l->Add(tGenNew1);
    if(bCompareWith1) l->Add(hGenNew1);
    TTree *tGenNew2 = new TTree("tGenNew2","tGenNew2");
    tGenNew2->Branch("fPtGen",&fPtGenerated,"fPtGen/D");
    l->Add(tGenNew2);
    l->Add(hGenNew2);

    // #############################################################################################
    if(bCompareWith1){
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
        Printf("STARlight tree (new R_A) contains %lli entries.", tSL_new1->GetEntries());
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
    }
    // #############################################################################################
    Printf("\n***");
    Printf("Analysis 2: nGen = 6.000.000");
    Printf("STARlight tree (6.624 fm) contains %lli entries.", tSL_old2->GetEntries());
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
    Printf("STARlight tree (new R_A) contains %lli entries.", tSL_new2->GetEntries());
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
            if(bCompareWith1) hRecOld1->Fill(fPt);
            hRecOld2->Fill(fPt);
        } 
    }

    if(bCompareWith1){
        hRatios1 = (TH1D*)hGenNew1->Clone("hRatios1");
        hRatios1->SetTitle("hRatios1");
        hRatios1->Sumw2();
        hRatios1->Divide(hGenOld1);

        hRecNew1 = (TH1D*)hRecOld1->Clone("hRecNew1");
        hRecNew1->SetTitle("hRecNew1");
        hRecNew1->Multiply(hRatios1);
    }

    hRatios2 = (TH1D*)hGenNew2->Clone("hRatios2");
    hRatios2->SetTitle("hRatios2");
    hRatios2->Sumw2();
    hRatios2->Divide(hGenOld2);

    hRecNew2 = (TH1D*)hRecOld2->Clone("hRecNew2");
    hRecNew2->SetTitle("hRecNew2");
    hRecNew2->Multiply(hRatios2);

    if(bCompareWith1) l->Add(hRatios1);
    l->Add(hRatios2);

    // Print the results to console
    if(bCompareWith1){
        Printf("\n***");
        Printf("*** Analysis 1: nGen = 4.880.000");    
        Printf("pT_low\tpT_upp\tnEv_old\tnEv_new\tratio\n");
        for(Int_t iBin = 1; iBin <= nBins1; iBin++){
            Printf("%.3f\t%.3f\t%.0f\t%.0f\t%.3f",
                hRatios1->GetBinLowEdge(iBin), hRatios1->GetBinLowEdge(iBin+1), 
                hGenOld1->GetBinContent(iBin), hGenNew1->GetBinContent(iBin), hRatios1->GetBinContent(iBin));
        }
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
    hRatios2->GetYaxis()->SetTitle(Form("#it{N}_{gen}(#it{R}_{A} = %.3f fm) / #it{N}_{gen}(#it{R}_{A} = 6.624 fm)", R_A));
    hRatios2->GetYaxis()->SetTitleSize(0.05);
    hRatios2->GetYaxis()->SetTitleOffset(0.9);
    hRatios2->GetYaxis()->SetLabelSize(0.05);
    hRatios2->GetYaxis()->SetDecimals(1);
    hRatios2->GetYaxis()->SetRangeUser(0.0,5.0);
    // Horizontal axis
    hRatios2->GetXaxis()->SetTitle("#it{p}_{T}^{gen} (GeV/#it{c})");
    hRatios2->GetXaxis()->SetTitleSize(0.05);
    hRatios2->GetXaxis()->SetTitleOffset(1.2);
    hRatios2->GetXaxis()->SetLabelSize(0.05);
    // Draw
    hRatios2->SetLineColor(215);
    hRatios2->Draw("E1");
    if(bCompareWith1){
        hRatios1->SetLineColor(kRed);
        hRatios1->Draw("E1 SAME");
    }

    TCanvas *cRec = new TCanvas("cRec", "cRec", 900, 600);
    cRec->SetTopMargin(0.04);
    cRec->SetBottomMargin(0.15);
    cRec->SetRightMargin(0.04);
    cRec->SetLeftMargin(0.1);
    cRec->SetLogy();
    // Vertical axis
    hRecOld2->GetYaxis()->SetTitle("#it{N}_{rec} (after selections)");
    hRecOld2->GetYaxis()->SetTitleSize(0.05);
    hRecOld2->GetYaxis()->SetTitleOffset(0.9);
    hRecOld2->GetYaxis()->SetLabelSize(0.05);
    hRecOld2->GetYaxis()->SetRangeUser(0.1,1e5);
    // Horizontal axis
    hRecOld2->GetXaxis()->SetTitle("#it{p}_{T}^{gen} (GeV/#it{c})");
    hRecOld2->GetXaxis()->SetTitleSize(0.05);
    hRecOld2->GetXaxis()->SetTitleOffset(1.2);
    hRecOld2->GetXaxis()->SetLabelSize(0.05);
    // Draw
    hRecOld2->SetLineColor(210);
    hRecOld2->Draw("E1");
    if(bCompareWith1){
        hRecNew1->SetLineColor(kRed);
        hRecNew1->Draw("E1 SAME");
    }
    hRecNew2->SetLineColor(215);
    hRecNew2->Draw("E1 SAME");


    // Print plots
    cRatios->Print(Form("Results/STARlight/%s/Ratios/coh_modRA_%s.pdf", str_out_subfolder.Data(), str_out_name.Data()));
    cRatios->Print(Form("Results/STARlight/%s/Ratios/coh_modRA_%s.png", str_out_subfolder.Data(), str_out_name.Data()));
    cRec->Print(Form("Results/STARlight/%s/RecSpectra/coh_modRA_%s.pdf", str_out_subfolder.Data(), str_out_name.Data()));
    cRec->Print(Form("Results/STARlight/%s/RecSpectra/coh_modRA_%s.png", str_out_subfolder.Data(), str_out_name.Data()));

    delete hGenOld1;
    delete hGenNew1;
    delete hRatios1;
    delete hRecOld1;
    delete hRecNew1;
    delete hGenOld2;
    delete hGenNew2;
    delete hRatios2;
    delete hRecOld2;
    delete hRecNew2;

    return;
}

void FillHistRec()
{
    

    return;
}