// PtFitUtilities.h
// David Grund, Sep 29, 2021

// cpp headers
#include <vector>
#include <stdio.h> // printf
#include <fstream> // print output to txt file
// root headers
#include "TFile.h"
#include "TList.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TStyle.h"
// roofit headers
#include "RooBinning.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooAbsReal.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"

using namespace RooFit;

TString NamesPDFs[6] = {"hCohJ","hIncJ","hCohP","hIncP","hBkg","hDiss"};
TString OutputPDFs = "Trees/PtFit/";

Double_t fPtLow = 0.0;
Double_t fPtUpp = 2.0;
Int_t nBins;

Int_t BinningOpt;
// == 0 => uniform binning with step of 10 MeV
// > 0  => variable binnings
Bool_t debugBinning = kFALSE;
RooBinning fPtBins(fPtLow, fPtUpp);

vector<Double_t> edges;
Double_t *ptEdges;

void SetPtBins(Int_t opt){

    // Set the global BinningOpt variable
    BinningOpt = opt;

    vector<Double_t> BinsSize;
    vector<Double_t> BinsMax;
    Int_t nBinTypes = 0;
    // Define binning
    if(opt == 0){
        // opt == 0 => uniform binning with the step of 10 MeV
        nBinTypes = 1;
        BinsSize = {0.010};     // GeV
        BinsMax = {0.00, 2.00}; // GeV
    } else if(opt == 1){
        // opt == 1 => variable binning (optimal for pt fit)
        nBinTypes = 3;
        BinsSize = {0.010, 0.025, 0.050};   // GeV
        BinsMax = {0.00, 0.40, 0.80, 2.00}; // GeV        
    } else if(opt == 2){
        // opt == 2 => variable binning optimal for PtFitWithoutBkg.c 
        nBinTypes = 4;
        BinsSize = {0.010, 0.020, 0.080, 0.400};  // GeV 
        BinsMax = {0.00, 0.20, 0.40, 1.20, 2.00}; // GeV 
    } else if(opt == 3){
        // opt == 3 => variable binning optimal for the SideBand PDF, see PtFitPrepareSideBandPDF.c
        nBinTypes = 4;
        BinsSize = {0.010, 0.025, 0.200, 0.400};  // GeV 
        BinsMax = {0.00, 0.10, 0.20, 1.20, 2.00}; // GeV 
    }
    // Create a vector containing the boundaries
    edges.push_back(0.);
    Double_t fPtNow = 0.0;
    for(Int_t iBinType = 0; iBinType < nBinTypes; iBinType++){
        Double_t nBinsThisType = (BinsMax[iBinType+1] - BinsMax[iBinType]) / BinsSize[iBinType];
        Printf("Bin type %i: width %.3f GeV, from %.3f GeV to %.3f GeV, nBins = %.1f", 
            iBinType+1, BinsSize[iBinType], BinsMax[iBinType], BinsMax[iBinType+1], nBinsThisType);
            Int_t iBin = 0;
        while(iBin < nBinsThisType){
            fPtNow += BinsSize[iBinType];
            edges.push_back(fPtNow);
            iBin++;
        }
    }
    nBins = edges.size() - 1;
    // Load the values from the vector to an array
    ptEdges = &edges[0];
    // Set RooBinning
    for(Int_t i = 0; i < nBins - 1; i++) fPtBins.addBoundary(edges[i + 1]);
    // For debugBinningging
    if(debugBinning){
        Printf("%i bins created with the following boundaries:", nBins);
        for(Int_t i = 0; i < nBins + 1; i++) Printf("%.2f", edges[i]);
        Printf("RooBinning:");
        for(Int_t i = 0; i < nBins; i++) Printf("%.2f", fPtBins.binLow(i));
        Printf("%.2f", fPtBins.binHigh(nBins - 1));
    }
    return;
}

void FillHistogramsMC(Int_t iMC, TH1D *hist){

    // Load the data
    TFile *file = NULL;
    switch(iMC){
        case 0:
            file = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohJpsiToMu.root", "read");
            break;
        case 1:
            file = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohJpsiToMu.root", "read");
            break;
        case 2:
            file = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kCohPsi2sToMuPi.root", "read");
            break;
        case 3:
            file = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohPsi2sToMuPi.root", "read");
            break;
        case 4:
            file = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kTwoGammaToMuMedium.root", "read");
            break;
    }
    if(file) Printf("File %s loaded.", file->GetName());

    TTree *tRec = dynamic_cast<TTree*> (file->Get("AnalysisOutput/fTreeJPsiMCRec"));
    if(tRec) Printf("MC rec tree loaded.");
    
    ConnectTreeVariablesMCRec(tRec);

    Printf("Tree %s has %lli entries.", tRec->GetName(), tRec->GetEntries());

    // Loop over tree entries
    Int_t nEntriesAnalysed = 0;
    Int_t nEvPassed = 0;

    for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++){
        tRec->GetEntry(iEntry);

        // iMassCut == 1 => 3.0 < m < 3.2 GeV
        // iPtCut == 2 => pt < 2.0 GeV
        if(EventPassedMCRec(1, 2)){
            nEvPassed++;
            hist->Fill(fPt);
        } 

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }
    Printf("Done.");
    Printf("%i events passed the selections.", nEvPassed);

    return;
}

void PreparePDFs_MC(){

    TFile *file = TFile::Open(Form("%sPDFs_MC_Binning%i.root", OutputPDFs.Data(), BinningOpt),"read");
    if(file){
        Printf("MC PDFs for this binning already created.");
        return;
    }

    // ***************************************************************
    // Go over MC data
    // Define output histograms with predefined binning to create PDFs
    TList *l = new TList();
    TH1D *HistPDFs[6] = { NULL };
    HistPDFs[0] = new TH1D(NamesPDFs[0].Data(), NamesPDFs[0].Data(), nBins, ptEdges);
    HistPDFs[1] = new TH1D(NamesPDFs[1].Data(), NamesPDFs[1].Data(), nBins, ptEdges);
    HistPDFs[2] = new TH1D(NamesPDFs[2].Data(), NamesPDFs[2].Data(), nBins, ptEdges);
    HistPDFs[3] = new TH1D(NamesPDFs[3].Data(), NamesPDFs[3].Data(), nBins, ptEdges);
    HistPDFs[4] = new TH1D(NamesPDFs[4].Data(), NamesPDFs[4].Data(), nBins, ptEdges);

    for(Int_t i = 0; i < 5; i++){
        FillHistogramsMC(i, HistPDFs[i]);
        l->Add(HistPDFs[i]);
    }
    // ***************************************************************
    // Create the dissociative PDF
    HistPDFs[5] = new TH1D(NamesPDFs[5].Data(), NamesPDFs[5].Data(), nBins, ptEdges);
    TF1 *fDissH1 = new TF1("fDissH1","x*pow((1 + x*x*[0]/[1]),-[1])",fPtLow,fPtUpp);
    Double_t b_pd = 1.79;
    Double_t n_pd = 3.58;
    fDissH1->SetParameter(0,b_pd);
    fDissH1->SetParameter(1,n_pd);
    Int_t i = 0;
    while(i < 1e6){
        HistPDFs[5]->Fill(fDissH1->GetRandom());
        i++;
    }
    //HistPDFs[5]->Draw();
    l->Add(HistPDFs[5]);

    // ***************************************************************
    // Save results to the output file
    // Create the output file
    TFile *f = new TFile(Form("%sPDFs_MC_Binning%i.root", OutputPDFs.Data(), BinningOpt),"RECREATE");
    l->Write("HistList", TObject::kSingleKey);
    f->ls();

    return;
}

void ConnectTreeVariablesPt(TTree *t){

    t->SetBranchAddress("fPt", &fPt);

    Printf("Variables from %s connected.", t->GetName());

    return;
}

void SetStyle(){

    gStyle->SetOptTitle(0); // suppress title
    gStyle->SetOptStat(0);  // the type of information printed in the histogram statistics box
                            // 0 = no information
    gStyle->SetPalette(1);  // set color map
    gStyle->SetPaintTextFormat("4.2f"); // precision if plotted with "TEXT"

    return;
}

void SetCanvas(TCanvas *c, Bool_t isLogScale){

    if(isLogScale == kTRUE) c->SetLogy();
    c->SetTopMargin(0.05);
    c->SetBottomMargin(0.12);
    c->SetRightMargin(0.02);

    return;
}

void DrawCorrelationMatrix(TCanvas *cCM, RooFitResult* ResFit){

    // Set margins
    cCM->SetTopMargin(0.03);
    cCM->SetBottomMargin(0.11);
    cCM->SetRightMargin(0.17);
    cCM->SetLeftMargin(0.15);
    // Get 2D corr hist
    TH2* hCorr = ResFit->correlationHist();
    // Set X axis
    hCorr->GetXaxis()->SetBinLabel(1,"#it{N}_{coh}");
    hCorr->GetXaxis()->SetBinLabel(2,"#it{N}_{diss}");
    hCorr->GetXaxis()->SetBinLabel(3,"#it{N}_{coh}");
    // Set Y axis
    hCorr->GetYaxis()->SetBinLabel(1,"#it{N}_{inc}");
    hCorr->GetYaxis()->SetBinLabel(2,"#it{N}_{diss}");
    hCorr->GetYaxis()->SetBinLabel(3,"#it{N}_{coh}");
    // Set corr hist and draw it
    hCorr->SetMarkerSize(3.6);
    hCorr->GetXaxis()->SetLabelSize(0.13);
    hCorr->GetYaxis()->SetLabelSize(0.13);
    hCorr->GetZaxis()->SetLabelSize(0.08);
    hCorr->Draw("colz,text");

    return;
}