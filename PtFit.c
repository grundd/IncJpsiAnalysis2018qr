// PtFit.c
// David Grund, Sep 28, 2021

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
#include "RooAbsData.h"
#include "RooAbsReal.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
// my headers
#include "AnalysisManager.h"

using namespace RooFit;

//#######################################
// Options to set:
Int_t BinningOpt = 0;
// == 0 => uniform binning with step of 10 MeV
// == 1 => variable binning
Bool_t isBkgData = kFALSE; // kTRUE if side-band data are used for bkg
Bool_t Debug = kFALSE;
TString OutputPDFs = "Trees/PtFit/";
TString OutputFigures = "Results/PtFit/";
//#######################################

TString NamesPDFs[6] = {"hCohJ","hIncJ","hCohP","hIncP","hBkg","hDiss"};
Double_t fPtLow = 0.0;
Double_t fPtUpp = 2.0;
Int_t nBins;
vector<Double_t> edges;
Double_t *ptEdges;
RooBinning fPtBins(fPtLow, fPtUpp);

void SetPtBins();
void PreparePDFs_MC();
void PreparePDFs_SideBand();
void FillHistogramsMC(Int_t iMC, TH1D *hist);
void DoPtFit(); 
void ConnectTreeVariablesPt(TTree *t);
void SetStyle();
void DrawCorrelationMatrix(TCanvas *cCM, RooFitResult* ResFit);

void PtFit(){

    SetPtBins();

    PreparePDFs_MC();

    //PreparePDFs_SideBand();

    DoPtFit();

    return;
}

void SetPtBins(){

    vector<Double_t> BinsSize;
    vector<Double_t> BinsMax;
    Int_t nBinTypes = 0;
    // Define binning
    if(BinningOpt == 0){
        // BinningOpt == 0 => uniform binning
        nBinTypes = 1;
        BinsSize = {0.010};     // GeV
        BinsMax = {0.00, 2.00}; // GeV
    } else if(BinningOpt == 1){
        // BinningOpt == 1 => variable binning
        nBinTypes = 3;
        BinsSize = {0.010, 0.025, 0.050};   // GeV
        BinsMax = {0.00, 0.40, 0.80, 2.00}; // GeV        
    } else if(BinningOpt == 2){
        // BinningOpt == 2 => variable binning
        nBinTypes = 3;
        BinsSize = {0.010, 0.050, 0.100};   // GeV
        BinsMax = {0.00, 0.15, 0.80, 2.00}; // GeV        
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
    // For debugging
    if(Debug){
        Printf("%i bins created with the following boundaries:", nBins);
        for(Int_t i = 0; i < nBins + 1; i++) Printf("%.2f", edges[i]);
        Printf("RooBinning:");
        for(Int_t i = 0; i < nBins; i++) Printf("%.2f", fPtBins.binLow(i));
        Printf("%.2f", fPtBins.binHigh(nBins - 1));
    }
    return;
}

void PreparePDFs_MC(){

    TFile *file = TFile::Open(Form("%sPDFs_MC_Binning%i.root", OutputPDFs.Data(), BinningOpt),"read");
    if(file){
        Printf("PDFs for these boundaries already created.");
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

void PreparePDFs_SideBand(){

    TList *l = new TList();

    TH1D *hSideBandPDF = new TH1D("hSideBandPDF", "hSideBandPDF", nBins, ptEdges);
    l->Add(hSideBandPDF);

    TFile *file = TFile::Open("Trees/AnalysisData/AnalysisResultsLHC18qrMerged.root", "read");
    if(file) Printf("File %s loaded.", file->GetName());

    TTree *fTreeIn = dynamic_cast<TTree*> (file->Get("AnalysisOutput/fTreeJPsi"));
    if(fTreeIn) Printf("Input tree loaded.");

    ConnectTreeVariables(fTreeIn);

    Printf("Tree %s has %lli entries.", fTreeIn->GetName(), fTreeIn->GetEntries());

    // Loop over tree entries
    Int_t nEntriesAnalysed = 0;
    Int_t nEvPassed = 0;

    for(Int_t iEntry = 0; iEntry < fTreeIn->GetEntries(); iEntry++){
        fTreeIn->GetEntry(iEntry);

        // iMassCut == 0 => 2.2 < m < 4.5 GeV
        // iPtCut   == 2 => pt < 2.0 GeV
        if(EventPassed(0, 2) && (fM > 3.3 && fM < 4.5)){ // 3.3 < m < 4.5 GeV
            nEvPassed++;
            hSideBandPDF->Fill(fPt);
        } 

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }
    Printf("Done.");
    Printf("%i events passed the selections.", nEvPassed);    

    TCanvas *c = new TCanvas("c","c",900,600);
    c->SetLogy();
    hSideBandPDF->Scale(1., "width");
    hSideBandPDF->Draw("L");

    TFile *f = new TFile(Form("%sPDFs_SideBand.root", OutputPDFs.Data()),"RECREATE");
    l->Write("HistList", TObject::kSingleKey);
    f->ls();

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

void DoPtFit(){

    SetPtBinning();

    // Load the data file
    TFile *file = TFile::Open(Form("%sPDFs_MC_Binning%i.root", OutputPDFs.Data(), BinningOpt),"read");
    if(file) Printf("Input file %s loaded.", file->GetName()); 

    TList *list = (TList*) file->Get("HistList");
    if(list) Printf("List %s loaded.", list->GetName()); 

    // Load histograms
    // 1) kCohJpsiToMu
    TH1D *hCohJ = (TH1D*)list->FindObject(NamesPDFs[0].Data());
    if(hCohJ) Printf("Histogram %s loaded.", hCohJ->GetName());

    // 2) kIncohJpsiToMu
    TH1D *hIncJ = (TH1D*)list->FindObject(NamesPDFs[1].Data());
    if(hIncJ) Printf("Histogram %s loaded.", hIncJ->GetName());

    // 3) kCohPsi2sToMuPi
    TH1D *hCohP = (TH1D*)list->FindObject(NamesPDFs[2].Data());
    if(hCohP) Printf("Histogram %s loaded.", hCohP->GetName());

    // 4) kincohPsi2sToMuPi
    TH1D *hIncP = (TH1D*)list->FindObject(NamesPDFs[3].Data());
    if(hIncP) Printf("Histogram %s loaded.", hIncP->GetName());

    // 6) Dissociative
    TH1D *hDiss = (TH1D*)list->FindObject(NamesPDFs[5].Data());
    if(hDiss) Printf("Histogram %s loaded.", hDiss->GetName());

    // 5) Background
    TH1D *hBkgr = NULL;
    if(!isBkgData){
        hBkgr = (TH1D*)list->FindObject(NamesPDFs[4].Data());
        file->Close();
    } else {
        file->Close();
        file = TFile::Open(Form("%sPDFs_SideBand.root", OutputPDFs.Data()),"read");
        if(file) Printf("Input file %s loaded.", file->GetName()); 

        TList *list = (TList*) file->Get("HistList");
        if(list) Printf("List %s loaded.", list->GetName()); 

        hBkgr = (TH1D*)list->FindObject("hSideBandPDF");
        file->Close();
    }
    if(hBkgr) Printf("Histogram %s loaded.", hBkgr->GetName());

    // Definition of roofit variables
    RooRealVar fPt("fPt", "fPt", fPtLow, fPtUpp);
    RooArgSet fSetOfVariables(fPt);

    // Create PDFs
    // 1) kCohJpsiToMu
    RooDataHist DHisCohJ("DHisCohJ","DHisCohJ",fPt,hCohJ);
    RooHistPdf  hPDFCohJ("hPDFCohJ","hPDFCohJ",fSetOfVariables,DHisCohJ,0);

    // 2) kIncohJpsiToMu
    RooDataHist DHisIncJ("DHisIncJ","DHisIncJ",fPt,hIncJ);
    RooHistPdf  hPDFIncJ("hPDFIncJ","hPDFIncJ",fSetOfVariables,DHisIncJ,0);

    // 3) kCohPsi2sToMuPi
    RooDataHist DHisCohP("DHisCohP","DHisCohP",fPt,hCohP);
    RooHistPdf  hPDFCohP("hPDFCohP","hPDFCohP",fSetOfVariables,DHisCohP,0);

    // 4) kincohPsi2sToMuPi
    RooDataHist DHisIncP("DHisIncP","DHisIncP",fPt,hIncP);
    RooHistPdf  hPDFIncP("hPDFIncP","hPDFIncP",fSetOfVariables,DHisIncP,0);

    // 5) Background
    RooDataHist DHisBkgr("DHisBkgr","DHisBkgr",fPt,hBkgr);
    RooHistPdf  hPDFBkgr("hPDFBkgr","hPDFBkgr",fSetOfVariables,DHisBkgr,0);

    // 6) Dissociative
    RooDataHist DHisDiss("DHisDiss","DHisDiss",fPt,hDiss);
    RooHistPdf  hPDFDiss("hPDFDiss","hPDFDiss",fSetOfVariables,DHisDiss,0);

    // Close the file with PDFs
    file->Close();

    // Get the dataset
    file = TFile::Open("Trees/PtFit/PtFitTrees.root", "read");
    if(file) Printf("Input file %s loaded.", file->GetName()); 

    TTree *tData = dynamic_cast<TTree*> (file->Get("fDataTree")); 
    if(tData) Printf("Tree %s loaded.", tData->GetName());

    ConnectTreeVariablesPt(tData);
    RooDataSet DSetData("DSetData","DSetData",fSetOfVariables,Import(*tData));
    Printf("Tree %s contains %lli entries.", tData->GetName(), tData->GetEntries());
    Printf("DataSet %s contains %i entries.", tData->GetName(), DSetData.numEntries());
    Printf("LHC18qr data loaded.");
    file->Close();
    
    // 8) Create the model for fitting
    // 8.1) Normalizations:
    RooRealVar NCohJ("NCohJ","Number of coh J/psi events", 0.90*DSetData.numEntries(),0.50*DSetData.numEntries(),1.0*DSetData.numEntries());
    RooRealVar NIncJ("NIncJ","Number of inc J/psi events", 0.05*DSetData.numEntries(),0.01*DSetData.numEntries(),0.3*DSetData.numEntries());
    RooRealVar NDiss("NDiss","Number of dis J/psi events", 0.05*DSetData.numEntries(),0.01*DSetData.numEntries(),0.3*DSetData.numEntries());
    // Load the number of bkg events
    Double_t N_bkg_val, N_bkg_err;
    ifstream file_in;
    file_in.open("Results/InvMassFit/all/all_bkg.txt");
    if(!(file_in.fail())){
        while(!file_in.eof()) file_in >> N_bkg_val >> N_bkg_err;
    } else {
        Printf("Results/InvMassFit/all/all_bkg.txt not found. Terminating...");
        return;
    }
    file_in.close();
    RooRealVar NBkgr("NBkgr","Number of background events",N_bkg_val,N_bkg_val-10,N_bkg_val+10);
    NBkgr.setVal(N_bkg_val);
    NBkgr.setConstant(kTRUE);
    Printf("Number of background events with 3.0 < m < 3.2 GeV: %.0f pm %.0f", N_bkg_val, N_bkg_err);

    // Load the values of fD coefficients
    Double_t fDCohCh, fDCohChErr, fDIncCh, fDIncChErr, fDCohNe, fDCohNeErr, fDIncNe, fDIncNeErr;
    char ch[8];
    file_in.open("Results/FeedDown/FeedDown_PtFit.txt");
    if(!(file_in.fail())){
        // Read data from the file
        Int_t i = 0;
        std::string str;
        while(std::getline(file_in,str)){
            istringstream in_stream(str);
            // skip first line
            if(i == 1) in_stream >> ch >> fDCohCh >> fDCohChErr >> fDIncCh >> fDIncChErr >> fDCohNe >> fDCohNeErr >> fDIncNe >> fDIncNeErr;
            i++;   
        } 
        file_in.close();
        Printf("fD coefficients loaded.");
    } else {
        Printf("fD coefficients missing. Terminating...");
        return;
    }

    Double_t fDCoh = (fDCohCh + fDCohNe) / 100;
    Double_t fDInc = (fDIncCh + fDIncNe) / 100;
    Printf("********************");
    Printf("fD_coh = %.4f", fDCoh);
    Printf("fD_inc = %.4f", fDInc);
    Printf("********************");
    
    RooGenericPdf NCohP("NCohP","Number of coh FD events",Form("NCohJ*%.4f", fDCoh),RooArgSet(NCohJ));
    RooGenericPdf NIncP("NIncP","Number of inc FD events",Form("NIncJ*%.4f", fDInc),RooArgSet(NIncJ));

    // 8.2) The model:
    RooAddPdf Mod("Mod","Sum of all PDFs",
        RooArgList(hPDFCohJ, hPDFIncJ, hPDFCohP, hPDFIncP, hPDFDiss, hPDFBkgr),
        RooArgList(NCohJ, NIncJ, NCohP, NIncP, NDiss, NBkgr)
    );

    // 9) Perform fitting
    RooFitResult* ResFit = Mod.fitTo(DSetData,Extended(kTRUE),Save(),Range(fPtLow,fPtUpp));

    // 10) Output to text file
    TString *str = NULL;
    if(!isBkgData) str = new TString(Form("%sPtFit_BkgMC_Binning%i", OutputFigures.Data(),BinningOpt));
    else str = new TString(Form("%sPtFit_BkgData_Binning%i", OutputFigures.Data(),BinningOpt));
    // Total number of entries in the dataset:
    Double_t N_all = DSetData.numEntries();
    // Integrals of the PDFs in the whole pt range 
    fPt.setRange("fPtAll",0.0,2.0);
    RooAbsReal *fN_CohJ_all = hPDFCohJ.createIntegral(fPt,NormSet(fPt),Range("fPtAll"));
    RooAbsReal *fN_IncJ_all = hPDFIncJ.createIntegral(fPt,NormSet(fPt),Range("fPtAll"));
    RooAbsReal *fN_Diss_all = hPDFDiss.createIntegral(fPt,NormSet(fPt),Range("fPtAll"));  
    RooAbsReal *fN_CohP_all = hPDFCohP.createIntegral(fPt,NormSet(fPt),Range("fPtAll")); 
    RooAbsReal *fN_IncP_all = hPDFIncP.createIntegral(fPt,NormSet(fPt),Range("fPtAll"));
    RooAbsReal *fN_Bkgr_all = hPDFBkgr.createIntegral(fPt,NormSet(fPt),Range("fPtAll"));   
    // Integrals of the PDFs in the incoherent-enriched sample (IES)
    fPt.setRange("fPtIES",0.2,2.0);
    RooAbsReal *fN_CohJ_ies = hPDFCohJ.createIntegral(fPt,NormSet(fPt),Range("fPtIES"));
    RooAbsReal *fN_IncJ_ies = hPDFIncJ.createIntegral(fPt,NormSet(fPt),Range("fPtIES"));
    RooAbsReal *fN_Diss_ies = hPDFDiss.createIntegral(fPt,NormSet(fPt),Range("fPtIES"));  
    RooAbsReal *fN_CohP_ies = hPDFCohP.createIntegral(fPt,NormSet(fPt),Range("fPtIES")); 
    RooAbsReal *fN_IncP_ies = hPDFIncP.createIntegral(fPt,NormSet(fPt),Range("fPtIES"));
    RooAbsReal *fN_Bkgr_ies = hPDFBkgr.createIntegral(fPt,NormSet(fPt),Range("fPtIES"));   
    // Number of events in the whole pt range
    Double_t N_CohJ_all = fN_CohJ_all->getVal()*NCohJ.getVal();
    Double_t N_IncJ_all = fN_IncJ_all->getVal()*NIncJ.getVal();
    Double_t N_Diss_all = fN_Diss_all->getVal()*NDiss.getVal();
    Double_t N_CohP_all = fN_CohP_all->getVal()*fDCoh*NCohJ.getVal();
    Double_t N_IncP_all = fN_IncP_all->getVal()*fDInc*NIncJ.getVal();
    Double_t N_Bkgr_all = fN_Bkgr_all->getVal()*NBkgr.getVal();
    // Number of events in the IES
    Double_t N_CohJ_ies = fN_CohJ_ies->getVal()*NCohJ.getVal();
    Double_t N_IncJ_ies = fN_IncJ_ies->getVal()*NIncJ.getVal();
    Double_t N_Diss_ies = fN_Diss_ies->getVal()*NDiss.getVal();
    Double_t N_CohP_ies = fN_CohP_ies->getVal()*fDCoh*NCohJ.getVal();
    Double_t N_IncP_ies = fN_IncP_ies->getVal()*fDInc*NIncJ.getVal();
    Double_t N_Bkgr_ies = fN_Bkgr_ies->getVal()*NBkgr.getVal();
    // Total fC correction (for the whole IES)
    Double_t fC = N_CohJ_ies / (N_IncJ_ies + N_Diss_ies);
    Double_t N_CohJ_ies_err = fN_CohJ_ies->getVal()*NCohJ.getError();
    Double_t N_IncJ_ies_err = fN_IncJ_ies->getVal()*NIncJ.getError();
    Double_t N_Diss_ies_err = fN_Diss_ies->getVal()*NDiss.getError();
    Double_t denominator_err = TMath::Sqrt(TMath::Power(N_IncJ_ies_err,2) + TMath::Power(N_Diss_ies_err,2));
    Double_t fC_err = fC * TMath::Sqrt(TMath::Power((N_CohJ_ies_err/N_CohJ_ies),2) + TMath::Power((denominator_err/(N_IncJ_ies + N_Diss_ies)),2));
    // Integrals of the PDFs in the pt bins
    RooAbsReal *fN_CohJ_bins[nPtBins] = { NULL };
    RooAbsReal *fN_IncJ_bins[nPtBins] = { NULL };
    RooAbsReal *fN_Diss_bins[nPtBins] = { NULL };
    RooAbsReal *fN_CohP_bins[nPtBins] = { NULL };
    RooAbsReal *fN_IncP_bins[nPtBins] = { NULL };
    RooAbsReal *fN_Bkgr_bins[nPtBins] = { NULL };
    Double_t N_CohJ_bins[nPtBins] = { 0 };
    Double_t N_IncJ_bins[nPtBins] = { 0 };
    Double_t N_Diss_bins[nPtBins] = { 0 };
    Double_t N_CohP_bins[nPtBins] = { 0 };
    Double_t N_IncP_bins[nPtBins] = { 0 };
    Double_t N_Bkgr_bins[nPtBins] = { 0 };
    Double_t fC_bins_val[nPtBins] = { 0 };
    Double_t fC_bins_err[nPtBins] = { 0 };
    for(Int_t i = 0; i < nPtBins; i++){
        fPt.setRange(Form("fPtBin%i",i+1), ptBoundaries[i], ptBoundaries[i+1]);
        Printf("Now calculating for bin %i, (%.3f, %.3f) GeV", i+1, ptBoundaries[i], ptBoundaries[i+1]);
        fN_CohJ_bins[i] = hPDFCohJ.createIntegral(fPt,NormSet(fPt),Range(Form("fPtBin%i",i+1)));
        fN_IncJ_bins[i] = hPDFIncJ.createIntegral(fPt,NormSet(fPt),Range(Form("fPtBin%i",i+1)));
        fN_Diss_bins[i] = hPDFDiss.createIntegral(fPt,NormSet(fPt),Range(Form("fPtBin%i",i+1)));
        fN_CohP_bins[i] = hPDFCohP.createIntegral(fPt,NormSet(fPt),Range(Form("fPtBin%i",i+1)));
        fN_IncP_bins[i] = hPDFIncP.createIntegral(fPt,NormSet(fPt),Range(Form("fPtBin%i",i+1)));
        fN_Bkgr_bins[i] = hPDFBkgr.createIntegral(fPt,NormSet(fPt),Range(Form("fPtBin%i",i+1)));  
        N_CohJ_bins[i] = fN_CohJ_bins[i]->getVal()*NCohJ.getVal();
        N_IncJ_bins[i] = fN_IncJ_bins[i]->getVal()*NIncJ.getVal();
        N_Diss_bins[i] = fN_Diss_bins[i]->getVal()*NDiss.getVal();
        N_CohP_bins[i] = fN_CohP_bins[i]->getVal()*fDCoh*NCohJ.getVal();
        N_IncP_bins[i] = fN_IncP_bins[i]->getVal()*fDInc*NIncJ.getVal();
        N_Bkgr_bins[i] = fN_Bkgr_bins[i]->getVal()*NBkgr.getVal();
        fC_bins_val[i] = N_CohJ_bins[i] / (N_IncJ_bins[i] + N_Diss_bins[i]);
        Double_t N_CohJ_bins_err = fN_CohJ_bins[i]->getVal()*NCohJ.getError();
        Double_t N_IncJ_bins_err = fN_IncJ_bins[i]->getVal()*NIncJ.getError();
        Double_t N_Diss_bins_err = fN_Diss_bins[i]->getVal()*NDiss.getError();
        Double_t denominator_err = TMath::Sqrt(TMath::Power(N_IncJ_bins_err,2) + TMath::Power(N_Diss_bins_err,2));
        fC_bins_err[i] = fC_bins_val[i] * TMath::Sqrt(TMath::Power((N_CohJ_bins_err/N_CohJ_bins[i]),2) + TMath::Power((denominator_err/(N_IncJ_bins[i] + N_Diss_bins[i])),2));    
    }
    // Print to text file
    ofstream outfile((*str + ".txt").Data());
    outfile << std::fixed << std::setprecision(2);
    outfile << Form("Dataset contains %.0f events.\n***\n", N_all);
    outfile << "In 0.0 < pt 2.0 GeV/c:\n";
    outfile << "NCohJ \tNIncJ \tNCohP \tNIncP \tNDiss \tNBkgr\n";
    outfile << N_CohJ_all << "\t" << N_IncJ_all << "\t" << N_CohP_all << "\t" << N_IncP_all << "\t" << N_Diss_all << "\t" << N_Bkgr_all << "\n***\n";
    outfile << "In 0.2 < pt 2.0 GeV/c:\n";
    outfile << "NCohJ \tNIncJ \tNCohP \tNIncP \tNDiss \tNBkgr \tfC \tfC_err \n";
    outfile << N_CohJ_ies << "\t" << N_IncJ_ies << "\t" << N_CohP_ies << "\t" << N_IncP_ies << "\t" << N_Diss_ies << "\t" << N_Bkgr_ies << "\t";
    outfile << std::fixed << std::setprecision(4);
    outfile << fC << "\t" << fC_err << "\n***\n";
    for(Int_t i = 0; i < nPtBins; i++){
        outfile << Form("In bin %i, (%.3f, %.3f) GeV/c:\n", i+1, ptBoundaries[i], ptBoundaries[i+1]);
        outfile << std::fixed << std::setprecision(2);
        outfile << "NCohJ \tNIncJ \tNCohP \tNIncP \tNDiss \tNBkgr \tfC \tfC_err \n";
        outfile << N_CohJ_bins[i] << "\t" << N_IncJ_bins[i] << "\t" << N_CohP_bins[i] << "\t" << N_IncP_bins[i] << "\t" << N_Diss_bins[i] << "\t" << N_Bkgr_bins[i] << "\t";
        outfile << std::fixed << std::setprecision(4);
        outfile << fC_bins_val[i] << "\t" << fC_bins_err[i] << "\n***\n";
    }
    outfile.close();
    Printf("*** Results printed to %s. ***", (*str + ".txt").Data());

    // 11) Plot the results
    SetStyle();

    // 11.1) Draw the Correlation Matrix
    TCanvas *cCM = new TCanvas("cCM","cCM",600,500);
    DrawCorrelationMatrix(cCM,ResFit);

    // 11.2) Draw the pt fit
    TCanvas *cPt = new TCanvas("cPt","cPt",900,600);
    cPt->SetLogy();
    cPt->SetTopMargin(0.05);
    cPt->SetBottomMargin(0.12);
    cPt->SetRightMargin(0.02);  

    RooPlot* PtFrame = fPt.frame(Title("Pt fit"));
    DSetData.plotOn(PtFrame,Name("DSetData"),MarkerStyle(20), MarkerSize(1.),Binning(fPtBins));
    Mod.plotOn(PtFrame,Name("Mod"),                           LineColor(215),   LineStyle(1),LineWidth(3),Normalization(DSetData.numEntries(),RooAbsReal::NumEvent));
    Mod.plotOn(PtFrame,Name("hPDFCohJ"),Components(hPDFCohJ), LineColor(222),   LineStyle(1),LineWidth(3),Normalization(DSetData.numEntries(),RooAbsReal::NumEvent));
    Mod.plotOn(PtFrame,Name("hPDFIncJ"),Components(hPDFIncJ), LineColor(kRed),  LineStyle(1),LineWidth(3),Normalization(DSetData.numEntries(),RooAbsReal::NumEvent));
    Mod.plotOn(PtFrame,Name("hPDFCohP"),Components(hPDFCohP), LineColor(222),   LineStyle(7),LineWidth(3),Normalization(DSetData.numEntries(),RooAbsReal::NumEvent));
    Mod.plotOn(PtFrame,Name("hPDFIncP"),Components(hPDFIncP), LineColor(kRed),  LineStyle(7),LineWidth(3),Normalization(DSetData.numEntries(),RooAbsReal::NumEvent));
    Mod.plotOn(PtFrame,Name("hPDFBkgr"),Components(hPDFBkgr), LineColor(kBlack),LineStyle(1),LineWidth(3),Normalization(DSetData.numEntries(),RooAbsReal::NumEvent));
    Mod.plotOn(PtFrame,Name("hPDFDiss"),Components(hPDFDiss), LineColor(15),    LineStyle(1),LineWidth(3),Normalization(DSetData.numEntries(),RooAbsReal::NumEvent));

    PtFrame->SetAxisRange(0,2,"X");
    // Set X axis
    PtFrame->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    PtFrame->GetXaxis()->SetTitleSize(0.05);
    PtFrame->GetXaxis()->SetLabelSize(0.05);
    // Set Y axis
    PtFrame->GetYaxis()->SetTitle(Form("Counts per bin"));
    PtFrame->GetYaxis()->SetTitleSize(0.05);
    PtFrame->GetYaxis()->SetTitleOffset(0.95);
    PtFrame->GetYaxis()->SetLabelSize(0.05);
    PtFrame->GetYaxis()->SetLabelOffset(0.01);
    PtFrame->Draw("][");

    // 12) Draw the legends
    // Legend 1
    TLegend *l1 = new TLegend(0.16,0.71,0.52,0.935);
    l1->SetHeader("ALICE, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV","r"); 
    l1->AddEntry((TObject*)0,Form("J/#psi #rightarrow #mu^{+}#mu^{-}"),"");
    l1->AddEntry((TObject*)0,Form("|#it{y}| < 0.8"),"");
    l1->AddEntry((TObject*)0,Form("#it{m}_{#mu#mu} #in (3.0,3.2) GeV/#it{c}^{2}"),"");
    l1->SetTextSize(0.05);
    l1->SetBorderSize(0); // no border
    l1->SetFillStyle(0);  // legend is transparent
    l1->Draw();

    // Legend 2
    TLegend *l2 = new TLegend(0.59,0.50,0.9,0.935);
    //leg2->SetTextSize(0.027);
    l2->AddEntry("DSetData","Data", "P");
    l2->AddEntry("Mod","sum","L");
    l2->AddEntry("hPDFCohJ","coherent J/#psi", "L");
    l2->AddEntry("hPDFIncJ","incoherent J/#psi", "L");
    l2->AddEntry("hPDFDiss","inc. J/#psi with nucl. diss.", "L");
    l2->AddEntry("hPDFCohP","J/#psi from coh. #psi(2#it{S}) decay", "L");
    l2->AddEntry("hPDFIncP","J/#psi from inc. #psi(2S) decay", "L");
    l2->AddEntry("hPDFBkg","#gamma#gamma to #mu#mu", "L");
    //l2->AddEntry((TObject*)0,Form("f_{C} = %.4f #pm %.4f", f_C, f_C_err),""); // (#it{p}_T #in (0.2,2.0) GeV/#it{c})
    l2->SetTextSize(0.043);
    l2->SetBorderSize(0);
    l2->SetFillStyle(0);
    l2->Draw();

    // 13) Print the results to pdf and png
    cCM->Print((*str + "_CM.pdf").Data());
    cCM->Print((*str + "_CM.png").Data());
    cPt->Print((*str + ".pdf").Data());
    cPt->Print((*str + ".png").Data());

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
