// PtFitWithoutBkg.c
// David Grund, Sep 27, 2021

// roofit headers
#include "RooCBShape.h"
// my headers
#include "AnalysisManager.h"
#include "PtFitUtilities.h"

using namespace RooFit;

// Main function
void SubtractBackground();
void DoPtFitNoBkg();
void DoInvMassFitMain(Double_t fPtCutLow, Double_t fPtCutUpp, Bool_t save = kFALSE, Int_t bin = -1);
//void DoInvMassFitMainMC(Double_t fPtCutLow, Double_t fPtCutUpp, Bool_t save = kFALSE, Int_t bin = -1);
void PrepareDataTree();

TString OutputPtFitWithoutBkg = "Results/PtFitWithoutBkg/";
TString OutputTrees = "Trees/PtFit/";

Double_t N_Jpsi_val = 0;
Double_t N_Jpsi_err = 0;

void PtFitWithoutBkg(){

    //PrepareDataTree();

    SetPtBins(2);

    SubtractBackground();

    PreparePDFs_MC();

    DoPtFitNoBkg();

    return;
}

void SubtractBackground(){

    TString path = Form("%sJpsiSignalNoBkg_Binning%i.root", OutputTrees.Data(), BinningOpt);
    TFile *file = TFile::Open(path.Data(),"read");
    if(file){
        Printf("Background already subtracted for this binning.");
        return;
    }

    TList *l = new TList();
    TH1D *hNJpsiBins = new TH1D("hNJpsiBins", "hNJpsiBins", nBins, ptEdges);
    l->Add(hNJpsiBins);

    Double_t fPt = fPtLow;
    Int_t iBin = 1;
    while(iBin <= nBins){
        DoInvMassFitMain(ptEdges[iBin-1], ptEdges[iBin], kTRUE, iBin);
        hNJpsiBins->SetBinContent(iBin, N_Jpsi_val);
        hNJpsiBins->SetBinError(iBin, N_Jpsi_err);
        iBin++;
    }

    TCanvas *cHist = new TCanvas("cHist", "cHist", 900, 600);
    cHist->SetLogy();    
    // Vertical axis
    hNJpsiBins->GetYaxis()->SetTitle("Counts per bin");
    // Horizontal axis
    hNJpsiBins->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    hNJpsiBins->GetXaxis()->SetDecimals(1);
    // Draw the histogram
    //hNJpsiBins->Scale(1., "width");
    hNJpsiBins->Draw("E0");

    // Save results to the output file
    // Create the output file
    file = new TFile(path.Data(),"RECREATE");
    l->Write("HistList", TObject::kSingleKey);
    file->ls();

    return;
}

void DoPtFitNoBkg(){

    SetPtBinning();

    // Load the file with PDFs
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

    // 6) Dissociative
    RooDataHist DHisDiss("DHisDiss","DHisDiss",fPt,hDiss);
    RooHistPdf  hPDFDiss("hPDFDiss","hPDFDiss",fSetOfVariables,DHisDiss,0);

    // Close the file with PDFs
    file->Close();


    // Get the binned dataset
    file = TFile::Open(Form("%sJpsiSignalNoBkg_Binning%i.root", OutputTrees.Data(), BinningOpt), "read");
    if(file) Printf("Input file %s loaded.", file->GetName()); 

    list = (TList*) file->Get("HistList");
    if(list) Printf("List %s loaded.", list->GetName()); 

    TH1D *hData = (TH1D*)list->FindObject("hNJpsiBins");
    if(hData) Printf("Histogram %s loaded.", hData->GetName());

    RooDataHist DHisData("DHisData","DHisData",fPt,hData);
    Printf("Binned data with background subtracted loaded.");
    // Calculate the number of entries
    Double_t N_all = 0;
    for(Int_t i = 1; i <= hData->GetNbinsX(); i++){
        N_all += hData->GetBinContent(i);
    }
    Printf("Data contain %.0f entries in %i bins.", N_all, DHisData.numEntries());
    file->Close();

    // 8) Create the model for fitting
    // 8.1) Normalizations:
    RooRealVar NCohJ("NCohJ","Number of coh J/psi events", 0.90*N_all,0.50*N_all,1.0*N_all);
    RooRealVar NIncJ("NIncJ","Number of inc J/psi events", 0.05*N_all,0.01*N_all,0.3*N_all);
    RooRealVar NDiss("NDiss","Number of dis J/psi events", 0.05*N_all,0.01*N_all,0.3*N_all);

    // Load the values of fD coefficients
    Double_t fDCohCh, fDCohChErr, fDIncCh, fDIncChErr, fDCohNe, fDCohNeErr, fDIncNe, fDIncNeErr;
    char ch[8];
    ifstream file_in;
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
        RooArgList(hPDFCohJ, hPDFIncJ, hPDFCohP, hPDFIncP, hPDFDiss),
        RooArgList(NCohJ, NIncJ, NCohP, NIncP, NDiss)
    );

    // 9) Perform fitting
    RooFitResult *ResFit = Mod.fitTo(DHisData,Extended(kTRUE),Save(),Range(fPtLow,fPtUpp));

    // ###############################################################################################################
    // ###############################################################################################################
    // 10) Output to text file
    TString *str = new TString(Form("%sPtFitNoBkg_Binning%i_%ibins", OutputPtFitWithoutBkg.Data(),BinningOpt, nPtBins));
    // Integrals of the PDFs in the whole pt range 
    fPt.setRange("fPtAll",0.0,2.0);
    RooAbsReal *fN_CohJ_all = hPDFCohJ.createIntegral(fPt,NormSet(fPt),Range("fPtAll"));
    RooAbsReal *fN_IncJ_all = hPDFIncJ.createIntegral(fPt,NormSet(fPt),Range("fPtAll"));
    RooAbsReal *fN_Diss_all = hPDFDiss.createIntegral(fPt,NormSet(fPt),Range("fPtAll"));  
    RooAbsReal *fN_CohP_all = hPDFCohP.createIntegral(fPt,NormSet(fPt),Range("fPtAll")); 
    RooAbsReal *fN_IncP_all = hPDFIncP.createIntegral(fPt,NormSet(fPt),Range("fPtAll"));
    // Integrals of the PDFs in the incoherent-enriched sample (IES)
    fPt.setRange("fPtIES",0.2,2.0);
    RooAbsReal *fN_CohJ_ies = hPDFCohJ.createIntegral(fPt,NormSet(fPt),Range("fPtIES"));
    RooAbsReal *fN_IncJ_ies = hPDFIncJ.createIntegral(fPt,NormSet(fPt),Range("fPtIES"));
    RooAbsReal *fN_Diss_ies = hPDFDiss.createIntegral(fPt,NormSet(fPt),Range("fPtIES"));  
    RooAbsReal *fN_CohP_ies = hPDFCohP.createIntegral(fPt,NormSet(fPt),Range("fPtIES")); 
    RooAbsReal *fN_IncP_ies = hPDFIncP.createIntegral(fPt,NormSet(fPt),Range("fPtIES"));  
    // Integrals of the PDFs with 0.2 < pt < 1.0 GeV/c
    fPt.setRange("fPtTo1",0.2,1.0);
    RooAbsReal *fN_CohJ_to1 = hPDFCohJ.createIntegral(fPt,NormSet(fPt),Range("fPtTo1"));
    RooAbsReal *fN_IncJ_to1 = hPDFIncJ.createIntegral(fPt,NormSet(fPt),Range("fPtTo1"));
    RooAbsReal *fN_Diss_to1 = hPDFDiss.createIntegral(fPt,NormSet(fPt),Range("fPtTo1"));  
    RooAbsReal *fN_CohP_to1 = hPDFCohP.createIntegral(fPt,NormSet(fPt),Range("fPtTo1")); 
    RooAbsReal *fN_IncP_to1 = hPDFIncP.createIntegral(fPt,NormSet(fPt),Range("fPtTo1"));  
    // Number of events in the whole pt range
    Double_t N_CohJ_all = fN_CohJ_all->getVal()*NCohJ.getVal();
    Double_t N_IncJ_all = fN_IncJ_all->getVal()*NIncJ.getVal();
    Double_t N_Diss_all = fN_Diss_all->getVal()*NDiss.getVal();
    Double_t N_CohP_all = fN_CohP_all->getVal()*fDCoh*NCohJ.getVal();
    Double_t N_IncP_all = fN_IncP_all->getVal()*fDInc*NIncJ.getVal();
    // Number of events in the IES
    Double_t N_CohJ_ies = fN_CohJ_ies->getVal()*NCohJ.getVal();
    Double_t N_IncJ_ies = fN_IncJ_ies->getVal()*NIncJ.getVal();
    Double_t N_Diss_ies = fN_Diss_ies->getVal()*NDiss.getVal();
    Double_t N_CohP_ies = fN_CohP_ies->getVal()*fDCoh*NCohJ.getVal();
    Double_t N_IncP_ies = fN_IncP_ies->getVal()*fDInc*NIncJ.getVal();
    // Number of events with 0.2 < pt < 1.0 GeV/c
    Double_t N_CohJ_to1 = fN_CohJ_to1->getVal()*NCohJ.getVal();
    Double_t N_IncJ_to1 = fN_IncJ_to1->getVal()*NIncJ.getVal();
    Double_t N_Diss_to1 = fN_Diss_to1->getVal()*NDiss.getVal();
    Double_t N_CohP_to1 = fN_CohP_to1->getVal()*fDCoh*NCohJ.getVal();
    Double_t N_IncP_to1 = fN_IncP_to1->getVal()*fDInc*NIncJ.getVal();
    // Total fC correction (for the whole IES)
    Double_t fC = N_CohJ_ies / (N_IncJ_ies + N_Diss_ies);
    Double_t N_CohJ_ies_err = fN_CohJ_ies->getVal()*NCohJ.getError();
    Double_t N_IncJ_ies_err = fN_IncJ_ies->getVal()*NIncJ.getError();
    Double_t N_Diss_ies_err = fN_Diss_ies->getVal()*NDiss.getError();
    Double_t denominator_err = TMath::Sqrt(TMath::Power(N_IncJ_ies_err,2) + TMath::Power(N_Diss_ies_err,2));
    Double_t fC_err = fC * TMath::Sqrt(TMath::Power((N_CohJ_ies_err/N_CohJ_ies),2) + TMath::Power((denominator_err/(N_IncJ_ies + N_Diss_ies)),2));
    // Total fD correction (for the whole IES)
    Double_t fDCoh_val = N_CohP_ies / (N_IncJ_ies + N_Diss_ies);
    Double_t fDInc_val = N_IncP_ies / (N_IncJ_ies + N_Diss_ies);
    Double_t N_CohP_ies_err = fN_CohP_ies->getVal()*fDCoh*NCohJ.getError();
    Double_t N_IncP_ies_err = fN_IncP_ies->getVal()*fDInc*NIncJ.getError();
    Double_t fDCoh_err = fDCoh_val * TMath::Sqrt(TMath::Power((N_CohP_ies_err/N_CohP_ies),2) + TMath::Power((denominator_err/(N_IncJ_ies + N_Diss_ies)),2));
    Double_t fDInc_err = fDInc_val * TMath::Sqrt(TMath::Power((N_IncP_ies_err/N_IncP_ies),2) + TMath::Power((denominator_err/(N_IncJ_ies + N_Diss_ies)),2));
    // Integrals of the PDFs in the pt bins
    RooAbsReal *fN_CohJ_bins[nPtBins] = { NULL };
    RooAbsReal *fN_IncJ_bins[nPtBins] = { NULL };
    RooAbsReal *fN_Diss_bins[nPtBins] = { NULL };
    RooAbsReal *fN_CohP_bins[nPtBins] = { NULL };
    RooAbsReal *fN_IncP_bins[nPtBins] = { NULL };
    Double_t N_CohJ_bins[nPtBins] = { 0 };
    Double_t N_IncJ_bins[nPtBins] = { 0 };
    Double_t N_Diss_bins[nPtBins] = { 0 };
    Double_t N_CohP_bins[nPtBins] = { 0 };
    Double_t N_IncP_bins[nPtBins] = { 0 };
    Double_t fC_bins_val[nPtBins] = { 0 };
    Double_t fC_bins_err[nPtBins] = { 0 };
    Double_t fDCoh_bins_val[nPtBins] = { 0 };
    Double_t fDCoh_bins_err[nPtBins] = { 0 };
    Double_t fDInc_bins_val[nPtBins] = { 0 };
    Double_t fDInc_bins_err[nPtBins] = { 0 };
    for(Int_t i = 0; i < nPtBins; i++){
        fPt.setRange(Form("fPtBin%i",i+1), ptBoundaries[i], ptBoundaries[i+1]);
        Printf("Now calculating for bin %i, (%.3f, %.3f) GeV", i+1, ptBoundaries[i], ptBoundaries[i+1]);
        fN_CohJ_bins[i] = hPDFCohJ.createIntegral(fPt,NormSet(fPt),Range(Form("fPtBin%i",i+1)));
        fN_IncJ_bins[i] = hPDFIncJ.createIntegral(fPt,NormSet(fPt),Range(Form("fPtBin%i",i+1)));
        fN_Diss_bins[i] = hPDFDiss.createIntegral(fPt,NormSet(fPt),Range(Form("fPtBin%i",i+1)));
        fN_CohP_bins[i] = hPDFCohP.createIntegral(fPt,NormSet(fPt),Range(Form("fPtBin%i",i+1)));
        fN_IncP_bins[i] = hPDFIncP.createIntegral(fPt,NormSet(fPt),Range(Form("fPtBin%i",i+1)));  
        N_CohJ_bins[i] = fN_CohJ_bins[i]->getVal()*NCohJ.getVal();
        N_IncJ_bins[i] = fN_IncJ_bins[i]->getVal()*NIncJ.getVal();
        N_Diss_bins[i] = fN_Diss_bins[i]->getVal()*NDiss.getVal();
        N_CohP_bins[i] = fN_CohP_bins[i]->getVal()*fDCoh*NCohJ.getVal();
        N_IncP_bins[i] = fN_IncP_bins[i]->getVal()*fDInc*NIncJ.getVal();
        fC_bins_val[i] = N_CohJ_bins[i] / (N_IncJ_bins[i] + N_Diss_bins[i]);
        Double_t N_CohJ_bins_err = fN_CohJ_bins[i]->getVal()*NCohJ.getError();
        Double_t N_IncJ_bins_err = fN_IncJ_bins[i]->getVal()*NIncJ.getError();
        Double_t N_Diss_bins_err = fN_Diss_bins[i]->getVal()*NDiss.getError();
        Double_t denominator_err = TMath::Sqrt(TMath::Power(N_IncJ_bins_err,2) + TMath::Power(N_Diss_bins_err,2));
        if(N_CohJ_bins[i] != 0) fC_bins_err[i] = fC_bins_val[i] * TMath::Sqrt(TMath::Power((N_CohJ_bins_err/N_CohJ_bins[i]),2) + TMath::Power((denominator_err/(N_IncJ_bins[i] + N_Diss_bins[i])),2));    
        else fC_bins_err[i] = 0.;
        Double_t N_CohP_bins_err = fN_CohP_bins[i]->getVal()*fDCoh*NCohJ.getError();
        Double_t N_IncP_bins_err = fN_IncP_bins[i]->getVal()*fDInc*NIncJ.getError();
        fDCoh_bins_val[i] = N_CohP_bins[i] / (N_IncJ_bins[i] + N_Diss_bins[i]);
        fDInc_bins_val[i] = N_IncP_bins[i] / (N_IncJ_bins[i] + N_Diss_bins[i]);
        fDCoh_bins_err[i] = fDCoh_bins_val[i] * TMath::Sqrt(TMath::Power((N_CohP_bins_err/N_CohP_bins[i]),2) + TMath::Power((denominator_err/(N_IncJ_bins[i] + N_Diss_bins[i])),2));
        fDInc_bins_err[i] = fDInc_bins_val[i] * TMath::Sqrt(TMath::Power((N_IncP_bins_err/N_IncP_bins[i]),2) + TMath::Power((denominator_err/(N_IncJ_bins[i] + N_Diss_bins[i])),2));
    }
    // Print to text file
    ofstream outfile((*str + ".txt").Data());
    outfile << std::fixed << std::setprecision(2);
    outfile << Form("Dataset contains %.0f events.\n***\n", N_all);
    outfile << "In 0.0 < pt 2.0 GeV/c:\n";
    Double_t sum_all = N_CohJ_all + N_IncJ_all + N_CohP_all + N_IncP_all + N_Diss_all;
    outfile << "NCohJ \tNIncJ \tNCohP \tNIncP \tNDiss \tSum\n";
    outfile << N_CohJ_all << "\t" << N_IncJ_all << "\t" << N_CohP_all << "\t" << N_IncP_all << "\t" << N_Diss_all << "\t" << sum_all << "\n***\n";
    outfile << "In 0.2 < pt 2.0 GeV/c:\n";
    Double_t sum_ies = N_CohJ_ies + N_IncJ_ies + N_CohP_ies + N_IncP_ies + N_Diss_ies;
    outfile << "NCohJ \tNIncJ \tNCohP \tNIncP \tNDiss \tSum \tfC \tfC_err \tfDCoh \tfDC_err\tfDInc \tfDI_err\n";
    outfile << N_CohJ_ies << "\t" << N_IncJ_ies << "\t" << N_CohP_ies << "\t" << N_IncP_ies << "\t" << N_Diss_ies << "\t" << sum_ies << "\t";
    outfile << std::fixed << std::setprecision(4);
    outfile << fC << "\t" << fC_err << "\t" << fDCoh_val << "\t" << fDCoh_err << "\t" << fDInc_val << "\t" << fDInc_err << "\n***\n";
    outfile << "In 0.2 < pt 1.0 GeV/c:\n";
    Double_t sum_to1 = N_CohJ_to1 + N_IncJ_to1 + N_CohP_to1 + N_IncP_to1 + N_Diss_to1;
    outfile << "NCohJ \tNIncJ \tNCohP \tNIncP \tNDiss \tSum \n";
    outfile << std::fixed << std::setprecision(2);
    outfile << N_CohJ_to1 << "\t" << N_IncJ_to1 << "\t" << N_CohP_to1 << "\t" << N_IncP_to1 << "\t" << N_Diss_to1 << "\t" << sum_to1 << "\n***\n";
    for(Int_t i = 0; i < nPtBins; i++){
        outfile << Form("In bin %i, (%.3f, %.3f) GeV/c:\n", i+1, ptBoundaries[i], ptBoundaries[i+1]);
        outfile << std::fixed << std::setprecision(2);
        outfile << "NCohJ \tNIncJ \tNCohP \tNIncP \tNDiss \tfC \tfC_err \tfDCoh \tfDC_err\tfDInc \tfDI_err\n";
        outfile << N_CohJ_bins[i] << "\t" << N_IncJ_bins[i] << "\t" 
                << N_CohP_bins[i] << "\t" << N_IncP_bins[i] << "\t" 
                << N_Diss_bins[i] << "\t";
        outfile << std::fixed << std::setprecision(4);
        outfile << fC_bins_val[i] << "\t" << fC_bins_err[i] << "\t"
                << fDCoh_bins_val[i] << "\t" << fDCoh_bins_err[i] << "\t"
                << fDInc_bins_val[i] << "\t" << fDInc_bins_err[i] << "\n***\n";
    }
    outfile << "Sum over bins:\n";
    outfile << "NCohJ \tNIncJ \tNCohP \tNIncP \tNDiss\n";
    outfile << std::fixed << std::setprecision(2);
    Double_t sum_bins[5] = { 0 };
    for(Int_t i = 0; i < nPtBins; i++){
        sum_bins[0] += N_CohJ_bins[i];
        sum_bins[1] += N_IncJ_bins[i];
        sum_bins[2] += N_CohP_bins[i];
        sum_bins[3] += N_IncP_bins[i];
        sum_bins[4] += N_Diss_bins[i];
    }
    for(Int_t i = 0; i < 4; i++){
        outfile << sum_bins[i] << "\t";
    }
    outfile << sum_bins[4] << "\n***\n";
    outfile.close();
    Printf("*** Results printed to %s. ***", (*str + ".txt").Data());
    // ###############################################################################################################
    // ###############################################################################################################

    // 11) Plot the results
    SetStyle();

    // 11.1) Draw the Correlation Matrix
    TCanvas *cCM = new TCanvas("cCM","cCM",600,500);
    DrawCorrelationMatrix(cCM,ResFit);

    // 11.2) Draw the pt fit
    TCanvas *cPt = new TCanvas("cPt","cPt",900,600);
    SetCanvas(cPt, kTRUE); 

    RooPlot* PtFrame = fPt.frame(Title("Pt fit"));
    DHisData.plotOn(PtFrame,Name("DSetData"),MarkerStyle(20), MarkerSize(1.),Binning(fPtBins));
    Mod.plotOn(PtFrame,Name("Mod"),                           LineColor(215),   LineStyle(1),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));
    Mod.plotOn(PtFrame,Name("hPDFCohJ"),Components(hPDFCohJ), LineColor(222),   LineStyle(1),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));
    Mod.plotOn(PtFrame,Name("hPDFIncJ"),Components(hPDFIncJ), LineColor(kRed),  LineStyle(1),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));
    Mod.plotOn(PtFrame,Name("hPDFCohP"),Components(hPDFCohP), LineColor(222),   LineStyle(7),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));
    Mod.plotOn(PtFrame,Name("hPDFIncP"),Components(hPDFIncP), LineColor(kRed),  LineStyle(7),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));
    Mod.plotOn(PtFrame,Name("hPDFDiss"),Components(hPDFDiss), LineColor(15),    LineStyle(1),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));

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
    TLegend *l2 = new TLegend(0.59,0.60,0.9,0.935);
    //leg2->SetTextSize(0.027);
    l2->AddEntry("DSetData","Data", "P");
    l2->AddEntry("Mod","sum","L");
    l2->AddEntry("hPDFCohJ","coherent J/#psi", "L");
    l2->AddEntry("hPDFIncJ","incoherent J/#psi", "L");
    l2->AddEntry("hPDFDiss","inc. J/#psi with nucl. diss.", "L");
    l2->AddEntry("hPDFCohP","J/#psi from coh. #psi(2#it{S}) decay", "L");
    l2->AddEntry("hPDFIncP","J/#psi from inc. #psi(2S) decay", "L");
    //l2->AddEntry((TObject*)0,Form("f_{C} = %.4f #pm %.4f", f_C, f_C_err),""); // (#it{p}_T #in (0.2,2.0) GeV/#it{c})
    l2->SetTextSize(0.043);
    l2->SetBorderSize(0);
    l2->SetFillStyle(0);
    l2->Draw();

    // 13) Print the results to pdf and png
    str = new TString(Form("%sPtFitNoBkg_Binning%i", OutputPtFitWithoutBkg.Data(),BinningOpt));
    cCM->Print((*str + "_CM.pdf").Data());
    cCM->Print((*str + "_CM.png").Data());
    cPt->Print((*str + ".pdf").Data());
    cPt->Print((*str + ".png").Data());

    return;
}


void DoInvMassFitMain(Double_t fPtCutLow, Double_t fPtCutUpp, Bool_t save, Int_t bin){
    // Fit the invariant mass distribution using Double-sided CB function
    // Fix the values of the tail parameters to MC values
    // Peak corresponding to psi(2s) excluded

    // Cuts:
    char fStrReduce[120];
    Double_t fYCut      = 0.80;
    Double_t fMCutLow   = 2.2;
    Double_t fMCutUpp   = 4.5;

    sprintf(fStrReduce,"abs(fY)<%f && fPt>%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCutLow,fPtCutUpp,fMCutLow,fMCutUpp);

    // Binning:
    Int_t nBins = 115; // so that each bin between 2.2 and 4.5 GeV is 20 MeV wide
    RooBinning binM(nBins,fMCutLow,fMCutUpp);
    Double_t BinSizeDouble = (fMCutUpp - fMCutLow) * 1000 / nBins; // in MeV
    BinSizeDouble = BinSizeDouble + 0.5;
    // https://stackoverflow.com/questions/9695329/c-how-to-round-a-double-to-an-int
    Int_t BinSize = (Int_t)BinSizeDouble;

    Printf("\n");
    Printf("*** Bin size (double): %.3f ***", BinSizeDouble);
    Printf("*** Bin size (int): %i ***\n", BinSize);
  
    // Roofit variables
    RooRealVar fM("fM","fM",fMCutLow,fMCutUpp);
    RooRealVar fPt("fPt","fPt",0,10.);
    RooRealVar fY("fY","fY",-0.8,0.8);

    // Get the data trees
    TFile *fFileIn = new TFile("Trees/PtFit/PtFitWithoutBkgTree.root"); 
    TTree *fTreeIn = NULL;
    fFileIn->GetObject("Tree",fTreeIn);
        
    RooDataSet *fDataIn = new RooDataSet("fDataIn", "fDataIn", RooArgSet(fM,fY,fPt), Import(*fTreeIn));
    RooAbsData* fDataSet = fDataIn->reduce(fStrReduce);

    // Print the number of entries in the dataset
    Int_t nEvents = fDataSet->numEntries();
    Printf("*** Number of events in the dataset: %i ***\n", nEvents);

    // Crystal Ball parameters from MC (to be fixed)
    Double_t fAlpha_L;
    Double_t fAlpha_R;
    Double_t fN_L;
    Double_t fN_R;

    char name[20];
    Double_t values[4];
    Double_t errors[4];

    TString* pathMC = NULL;
    //pathMC = new TString(Form("%sInvMassFitsBinsMC/Binning%i/bin%i.txt", OutputPtFitWithoutBkg.Data(), BinningOpt, bin));
    if(fPtCutLow < 0.20) pathMC = new TString("Results/InvMassFitMC/coh/coh.txt");
    else pathMC = new TString("Results/InvMassFitMC/inc/inc.txt");

    ifstream file_in;
    file_in.open(pathMC->Data());
    if(file_in.fail()){
        Printf("\n");
        Printf("*** Warning! ***");
        Printf("*** MC values for tail parameters not found. Terminating... *** \n");
        return;
    } else {
        Int_t i_line = 0;
        while(!file_in.eof()){
            file_in >> name >> values[i_line] >> errors[i_line];
            i_line++;
        }
        file_in.close();
    }

    /*
    if(file_in.fail()){
        Printf("\n");
        Printf("*** MC values for this bin not found. Will be calculated... *** \n");
        DoInvMassFitMainMC(fPtCutLow,fPtCutUpp,kTRUE,bin);
    } 

    file_in.open(pathMC->Data());
    Int_t i_line = 0;
    while(!file_in.eof()){
        file_in >> name >> values[i_line] >> errors[i_line];
        i_line++;
    }
    file_in.close();
    */

    fAlpha_L = values[0];
    fAlpha_R = values[1];
    fN_L = values[2];
    fN_R = values[3];

    // RooFit: definition of tail parameters
    // DSCB = Double-sided Crystal Ball function
    RooRealVar alpha_L("alpha_L","alpha_L from DSCB",fAlpha_L,0.,10.);
    RooRealVar alpha_R("alpha_R","alpha_R from DSCB",fAlpha_R,-10.,0.);
    RooRealVar n_L("n_L","n_L from DSCB",fN_L,0.,20.);
    RooRealVar n_R("n_R","n_R from DSCB",fN_R,0.,20.);
    alpha_L.setConstant(kTRUE);
    alpha_R.setConstant(kTRUE);
    n_L.setConstant(kTRUE);
    n_R.setConstant(kTRUE);

    // Crystal Ball for J/Psi
    RooRealVar mass_Jpsi("mass_Jpsi","J/psi mass",3.097,3.00,3.20); 
    //mass_Jpsi.setConstant(kTRUE);
    RooRealVar sigma_Jpsi("sigma_Jpsi","J/psi resolution",0.08,0.01,0.1);
    RooGenericPdf mean_R("mean_R","J/psi mass","mass_Jpsi",RooArgSet(mass_Jpsi));
    RooGenericPdf sigma_R("sigma_R","J/psi resolution","sigma_Jpsi",RooArgSet(sigma_Jpsi));
    RooRealVar N_Jpsi("N_Jpsi","number of J/psi events",0.4*nEvents,0,nEvents);

    // Background
    RooRealVar lambda("lambda","background exp",-1.2,-10.,0.);
    RooRealVar N_bkg("N_bkg","number of background events",0.6*nEvents,0,nEvents);

    // Functions for fitting
    // J/psi:
    RooCBShape CB_left("CB_left","CB_left",fM,mass_Jpsi,sigma_Jpsi,alpha_L,n_L);
    RooCBShape CB_right("CB_right","CB_right",fM,mean_R,sigma_R,alpha_R,n_R);
    RooRealVar frac("frac","fraction of CBs",0.5);
    RooAddPdf DoubleSidedCB("DoubleSidedCB","DoubleSidedCB",RooArgList(CB_left,CB_right),RooArgList(frac));
    // Background:
    RooGenericPdf BkgPdf("BkgPdf","exp(fM*lambda)",RooArgSet(fM,lambda));

    // Create Model
    RooAddPdf DSCBAndBkgPdf("DSCBAndBkgPdf","Double sided CB and background PDFs", RooArgList(DoubleSidedCB,BkgPdf), RooArgList(N_Jpsi,N_bkg));
    // Perform fit
    RooFitResult* fResFit = DSCBAndBkgPdf.fitTo(*fDataSet,Extended(kTRUE),Range(fMCutLow,fMCutUpp),Save());

    // Calculate the number of J/psi events
    Double_t N_Jpsi_out[2];
    fM.setRange("JpsiMassRange",3.0,3.2);
    RooAbsReal *intDSCB = DoubleSidedCB.createIntegral(fM,NormSet(fM),Range("JpsiMassRange"));
    // Integral of the normalized PDF, DSCB => will range from 0 to 1

    N_Jpsi_out[0] = intDSCB->getVal()*N_Jpsi.getVal();
    N_Jpsi_out[1] = intDSCB->getVal()*N_Jpsi.getError();
    N_Jpsi_val = N_Jpsi_out[0];
    N_Jpsi_err = N_Jpsi_out[1];

    // ##########################################################
    // Plot the results
    SetStyle();

    // Draw histogram with fit results
    TCanvas *cHist = new TCanvas("cHist","cHist",800,600);
    cHist->SetTopMargin(0.055);
    cHist->SetBottomMargin(0.12);
    cHist->SetRightMargin(0.03);
    cHist->SetLeftMargin(0.11);

    RooPlot* fFrameM = fM.frame(Title("Mass fit")); 
    fDataSet->plotOn(fFrameM,Name("fDataSet"),Binning(binM),MarkerStyle(20),MarkerSize(1.));
    DSCBAndBkgPdf.plotOn(fFrameM,Name("DoubleSidedCB"),Components(DoubleSidedCB),LineColor(kBlack),LineStyle(kDashed),LineWidth(3));
    DSCBAndBkgPdf.plotOn(fFrameM,Name("BkgPdf"),Components(BkgPdf),LineColor(kRed),LineStyle(kDashed),LineWidth(3));
    DSCBAndBkgPdf.plotOn(fFrameM,Name("DSCBAndBkgPdf"),LineColor(215),LineWidth(3));
    // Vertical axis
    fFrameM->GetYaxis()->SetTitle(Form("Counts per %i MeV/#it{c}^{2}", BinSize));
    fFrameM->GetYaxis()->SetTitleSize(0.05);
    fFrameM->GetYaxis()->SetTitleOffset(1.1);
    fFrameM->GetYaxis()->SetLabelSize(0.05);
    fFrameM->GetYaxis()->SetLabelOffset(0.01);
    fFrameM->GetYaxis()->SetMaxDigits(3);
    // Horizontal axis
    fFrameM->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    fFrameM->GetXaxis()->SetTitleSize(0.05);
    fFrameM->GetXaxis()->SetLabelSize(0.05);
    fFrameM->GetXaxis()->SetDecimals(1);
    fFrameM->Draw();

    // Get chi2 
    Double_t chi2 = fFrameM->chiSquare("DSCBAndBkgPdf","data",fResFit->floatParsFinal().getSize());

    // -------------------------------------------------------------------------------- 
    // Legend1
    TLegend *l1 = new TLegend(0.09,0.76,0.3,0.935);
    //l1->SetHeader("ALICE, PbPb #sqrt{#it{s}_{NN}} = 5.02 TeV","r"); 
    l1->AddEntry((TObject*)0,Form("J/#psi #rightarrow #mu^{+}#mu^{-}"),"");
    l1->AddEntry((TObject*)0,Form("|#it{y}| < %.1f", fYCut),"");
    // Print the pt cut
    l1->AddEntry((TObject*)0,Form("#it{p}_{T} #in (%.2f,%.2f) GeV/#it{c}", fPtCutLow,fPtCutUpp),"");
    l1->SetTextSize(0.042);
    l1->SetBorderSize(0); // no border
    l1->SetFillStyle(0);  // legend is transparent
    l1->Draw();

    TLegend *lTitle = new TLegend(0.325,0.88,0.95,0.935);
    lTitle->AddEntry((TObject*)0,"ALICE, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV","");
    lTitle->SetTextSize(0.05);
    lTitle->SetBorderSize(0);
    lTitle->SetFillStyle(0);
    lTitle->Draw();

    // Legend2
    TLegend *l2 = new TLegend(0.465,0.29,0.95,0.87);
    //l2->SetHeader("ALICE, PbPb #sqrt{#it{s}_{NN}} = 5.02 TeV","r"); 
    l2->AddEntry("DSCBAndBkgPdf","sum","L");
    //l2->AddEntry((TObject*)0,Form("#chi^{2}/NDF = %.3f",chi2),"");
    l2->AddEntry("DoubleSidedCB","J/#psi signal","L");
    l2->AddEntry((TObject*)0,Form("#it{N}_{J/#psi} = %.0f #pm %.0f",N_Jpsi_out[0],N_Jpsi_out[1]),"");
    l2->AddEntry((TObject*)0,Form("#it{M}_{J/#psi} = %.3f #pm %.3f GeV/#it{c}^{2}", mass_Jpsi.getVal(), mass_Jpsi.getError()),"");
    l2->AddEntry((TObject*)0,Form("#sigma = %.3f #pm %.3f GeV/#it{c}^{2}", sigma_Jpsi.getVal(), sigma_Jpsi.getError()),"");
    l2->AddEntry((TObject*)0,Form("#alpha_{L} = %.3f", alpha_L.getVal()),"");
    l2->AddEntry((TObject*)0,Form("#alpha_{R} = %.3f", (-1)*(alpha_R.getVal())),"");
    l2->AddEntry((TObject*)0,Form("#it{n}_{L} = %.2f", n_L.getVal()),"");
    l2->AddEntry((TObject*)0,Form("#it{n}_{R} = %.2f", n_R.getVal()),"");
    l2->AddEntry("BkgPdf","background","L");
    l2->AddEntry((TObject*)0,Form("#lambda = %.3f #pm %.3f GeV^{-1}#it{c}^{2}",lambda.getVal(), lambda.getError()),"");
    l2->SetTextSize(0.042);
    l2->SetBorderSize(0);
    l2->SetFillStyle(0);
    l2->Draw();

    if(save){
        // Prepare path
        TString *str = NULL;
        str = new TString(Form("%sInvMassFitsBins/Binning%i/bin%i", OutputPtFitWithoutBkg.Data(), BinningOpt, bin));
        // Print the plots
        cHist->Print((*str + ".png").Data()); 
    }

    return;
}

/*
void DoInvMassFitMainMC(Double_t fPtCutLow, Double_t fPtCutUpp, Bool_t save, Int_t bin){
    // Fit the invariant mass distribution using Double-sided CB function
    // Peak corresponding to psi(2s) excluded

    // Cuts:
    char fStrReduce[120];
    Double_t fYCut      = 0.80;
    Double_t fMCutLow   = 2.95;
    Double_t fMCutUpp   = 3.25;

    sprintf(fStrReduce,"abs(fY)<%f && fPt>%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCutLow,fPtCutUpp,fMCutLow,fMCutUpp);

    // Binning:
    Int_t nBins = 75; // so that each bin between 2.95 and 3.25 GeV is 4 MeV wide
    RooBinning binM(nBins,fMCutLow,fMCutUpp);
    Double_t BinSizeDouble = (fMCutUpp - fMCutLow) * 1000 / nBins; // in MeV
    BinSizeDouble = BinSizeDouble + 0.5;
    // https://stackoverflow.com/questions/9695329/c-how-to-round-a-double-to-an-int
    Int_t BinSize = (Int_t)BinSizeDouble;

    Printf("\n");
    Printf("*** Bin size (double): %.3f ***", BinSizeDouble);
    Printf("*** Bin size (int): %i ***\n", BinSize);

    // Roofit variables
    RooRealVar fM("fM","fM",fMCutLow,fMCutUpp);
    RooRealVar fPt("fPt","fPt",0,10.);
    RooRealVar fY("fY","fY",-0.8,0.8);

    //fM.setBinning(binM);

    // Get the data trees
    TFile *fFileIn = new TFile("Trees/InvMassFitMC/InvMassFitMC.root"); 
    TTree *fTreeIn = NULL;
    fFileIn->GetObject("tMixedSample",fTreeIn);

    RooDataSet *fDataIn = new RooDataSet("fDataIn", "fDataIn", RooArgSet(fM,fY,fPt), Import(*fTreeIn));
    RooAbsData* fDataSet = fDataIn->reduce(fStrReduce);

    // Print the number of entries in the dataset
    Int_t nEvents = fDataSet->numEntries();
    Printf("*** Number of events in the dataset: %i ***\n", nEvents);

    // Roofit variables
    RooRealVar norm_L("norm_L","N_{L}(J/#psi)",nEvents,0,1e06);
    RooRealVar norm_R("norm_R","N_{R}(J/#psi)",nEvents,0,1e06);
    RooRealVar N("N","N(J/#psi)",nEvents,0,1e06);

    RooRealVar mean_L("m","m_{J/#psi}",3.097,3.0,3.2);
    RooRealVar sigma_L("sig","#sigma_{J/#psi}",0.0186,0.01,0.2);
    RooRealVar alpha_L("#alpha_{L}","alpha_{L}",1.,0.0,20.0);
    RooRealVar n_L("n_{L}","n_{L}",1.,0,30);

    RooGenericPdf mean_R("mean_R","m_{J/#psi}","m",RooArgSet(mean_L));
    RooGenericPdf sigma_R("sigma_R","#sigma_{J/#psi}","sig",RooArgSet(sigma_L));
    RooRealVar alpha_R("#alpha_{R}","alpha_{R}",-1.,-20.0,0.0); 
    RooRealVar n_R("n_{R}","n_{R}",1.,0,30);

    RooCBShape CB_left("CB_left","CB_left",fM,mean_L,sigma_L,alpha_L,n_L);
    RooCBShape CB_right("CB_right","CB_right",fM,mean_R,sigma_R,alpha_R,n_R);
    RooRealVar frac("frac","fraction of CBs",0.5);
    RooAddPdf DoubleSidedCB("DoubleSidedCB","DoubleSidedCB",RooArgList(CB_left,CB_right),RooArgList(frac));

    // Create model
    RooExtendPdf DSCBExtended("DSCBExtended","Extended DSCB",DoubleSidedCB,N);
    // Perform fit
    RooFitResult* fResFit = DSCBExtended.fitTo(*fDataSet,Extended(kTRUE),Range(fMCutLow,fMCutUpp),Save());

    // ##########################################################
    // Plot the results
    // Draw histogram and fit
    TCanvas *cHist = new TCanvas("cHist","cHist",800,600);
    SetCanvas(cHist,kFALSE);

    RooPlot* frameM = fM.frame(Title("Mass fit")); 
    fDataSet->plotOn(frameM,Name("data"),Binning(binM),MarkerStyle(20),MarkerSize(1.));
    DoubleSidedCB.plotOn(frameM,Name("DoubleSidedCB"),LineColor(215),LineWidth(3),LineStyle(kDashed));
    // Y axis
    frameM->GetYaxis()->SetTitleSize(0.045);
    frameM->GetYaxis()->SetLabelSize(0.045);
    frameM->GetYaxis()->SetLabelOffset(0.01);
    frameM->GetYaxis()->SetTitle(Form("Counts per %i MeV/#it{c}^{2}", BinSize));
    frameM->GetYaxis()->SetTitleOffset(1);
    frameM->GetYaxis()->SetMaxDigits(3);
    // X axis
    frameM->GetXaxis()->SetTitleSize(0.045);
    frameM->GetXaxis()->SetLabelSize(0.045);
    frameM->GetXaxis()->SetLabelOffset(0.01);
    frameM->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    frameM->GetXaxis()->SetTitleOffset(1.1);
    frameM->Draw();

    // Get chi2 
    Double_t chi2 = frameM->chiSquare("CB_ext","data",fResFit->floatParsFinal().getSize());

    TLegend *leg = new TLegend(0.655,0.5,0.945,0.935);
    leg->AddEntry((TObject*)0,Form("#it{N} = %.f #pm %.f", N.getVal(), N.getError()),"");
    leg->AddEntry((TObject*)0,Form("#mu = %.4f GeV/#it{c}^{2}", mean_L.getVal()),""); // mean_L.getError()
    leg->AddEntry((TObject*)0,Form("#sigma = %.4f GeV/#it{c}^{2}", sigma_L.getVal()),""); // sigma_L.getError()
    leg->AddEntry((TObject*)0,Form("#alpha_{L} = %.3f #pm %.3f", alpha_L.getVal(), alpha_L.getError()),"");
    leg->AddEntry((TObject*)0,Form("#alpha_{R} = %.3f #pm %.3f", (-1)*(alpha_R.getVal()), alpha_R.getError()),"");
    leg->AddEntry((TObject*)0,Form("#it{n}_{L} = %.2f #pm %.2f", n_L.getVal(), n_L.getError()),"");
    leg->AddEntry((TObject*)0,Form("#it{n}_{R} = %.2f #pm %.2f", n_R.getVal(), n_R.getError()),"");
    leg->SetTextSize(0.042);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();

    TLegend *leg2 = new TLegend(0.10,0.7,0.3,0.935);
    //leg2->SetHeader("ALICE, PbPb #sqrt{#it{s}_{NN}} = 5.02 TeV","r"); 
    leg2->AddEntry((TObject*)0,Form("ALICE Simulation"),"");
    leg2->AddEntry((TObject*)0,Form("Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV"),"");
    leg2->AddEntry((TObject*)0,Form("MC rec: J/#psi #rightarrow #mu^{+}#mu^{-}"),"");
    leg2->AddEntry((TObject*)0,Form("#it{p}_{T} #in (%.2f,%.2f) GeV/#it{c}", fPtCutLow,fPtCutUpp),"");
    leg2->SetTextSize(0.042);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->Draw();

    // Draw Histogram with log scale
    TCanvas *cHistLog = new TCanvas("cHistLog","cHistLog",800,600);
    SetCanvas(cHistLog,kTRUE);
    frameM->Draw();
    leg->Draw();
    leg2->Draw();

    if(save){
        // Prepare path
        TString *str = NULL;
        str = new TString(Form("%sInvMassFitsBinsMC/Binning%i/bin%i", OutputPtFitWithoutBkg.Data(), BinningOpt, bin));
        // Print the plots
        cHist->Print((*str + ".png").Data());
        cHistLog->Print((*str + "_log.png").Data());
        // Print the values of alpha and n to txt output files
        ofstream outfile((*str + ".txt").Data());
        outfile << "alpha_L \t" << alpha_L.getVal() << "\t" << alpha_L.getError() << "\n";
        outfile << "alpha_R \t" << alpha_R.getVal() << "\t" << alpha_R.getError() << "\n";
        outfile << "n_L \t" << n_L.getVal() << "\t" << n_L.getError() << "\n";
        outfile << "n_R \t" << n_R.getVal() << "\t" << n_R.getError() << "\n";
        outfile.close();
        Printf("*** Results printed to %s.***", (*str + ".txt").Data());
    }

    return;
}
*/

void PrepareDataTree(){

    TFile *fFileIn = TFile::Open("Trees/AnalysisData/AnalysisResultsLHC18qrMerged.root", "read");
    if(fFileIn) Printf("Input data loaded.");

    TTree *fTreeIn = dynamic_cast<TTree*> (fFileIn->Get("AnalysisOutput/fTreeJPsi"));
    if(fTreeIn) Printf("Input tree loaded.");

    ConnectTreeVariables(fTreeIn);

    // Create new data tree with applied cuts
    TFile fFileOut("Trees/PtFit/PtFitWithoutBkgTree.root","RECREATE");

    TTree *Tree = new TTree("Tree", "Tree");
    Tree->Branch("fPt", &fPt, "fPt/D");
    Tree->Branch("fM", &fM, "fM/D");
    Tree->Branch("fY", &fY, "fY/D");

    Printf("%lli entries found in the tree.", fTreeIn->GetEntries());
    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < fTreeIn->GetEntries(); iEntry++){
        fTreeIn->GetEntry(iEntry);
        if(EventPassed(0,2)) Tree->Fill();

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }

    fFileOut.Write("",TObject::kWriteDelete);

    return;
}