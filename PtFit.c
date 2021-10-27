// PtFit.c
// David Grund, Sep 28, 2021

// my headers
#include "AnalysisManager.h"
#include "PtFitUtilities.h"

//#######################################
// Options to set:
Bool_t isBkgData = kFALSE; // kTRUE if side-band data are used for bkg
TString OutputPtFit = "Results/PtFit/";
//#######################################

void DoPtFit(); 

void PtFit(){

    SetPtBins(1);

    PreparePDFs_MC();

    DoPtFit();

    return;
}

void DoPtFit(){

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

    // 5) Background
    TH1D *hBkgr = NULL;
    if(!isBkgData){
        hBkgr = (TH1D*)list->FindObject(NamesPDFs[4].Data());
        file->Close();
    } else {
        file->Close();
        file = TFile::Open("Trees/PtFit/SideBandPDF/PDF_SideBand.root","read");
        if(file) Printf("Input file %s loaded.", file->GetName()); 
        else Printf("File %s missing. Create it using PtFitPrepareSideBandPDF.c. terminating...", file->GetName());

        TList *list = (TList*) file->Get("HistList");
        if(list) Printf("List %s loaded.", list->GetName()); 

        hBkgr = (TH1D*)list->FindObject("hSideBand");
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

    // ###############################################################################################################
    // ###############################################################################################################
    // 10) Output to text file
    TString *str = NULL;
    if(!isBkgData) str = new TString(Form("%sPtFit_BkgMC_Binning%i_%ibins", OutputPtFit.Data(),BinningOpt, nPtBins));
    else str = new TString(Form("%sPtFit_BkgData_Binning%i_%ibins", OutputPtFit.Data(),BinningOpt, nPtBins));
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
    // Integrals of the PDFs with 0.2 < pt < 2.0 GeV/c
    fPt.setRange("fPtIES",0.2,2.0);
    RooAbsReal *fN_CohJ_ies = hPDFCohJ.createIntegral(fPt,NormSet(fPt),Range("fPtIES"));
    RooAbsReal *fN_IncJ_ies = hPDFIncJ.createIntegral(fPt,NormSet(fPt),Range("fPtIES"));
    RooAbsReal *fN_Diss_ies = hPDFDiss.createIntegral(fPt,NormSet(fPt),Range("fPtIES"));  
    RooAbsReal *fN_CohP_ies = hPDFCohP.createIntegral(fPt,NormSet(fPt),Range("fPtIES")); 
    RooAbsReal *fN_IncP_ies = hPDFIncP.createIntegral(fPt,NormSet(fPt),Range("fPtIES"));
    RooAbsReal *fN_Bkgr_ies = hPDFBkgr.createIntegral(fPt,NormSet(fPt),Range("fPtIES"));  
    // Integrals of the PDFs with 0.2 < pt < 1.0 GeV/c
    fPt.setRange("fPtTo1",0.2,1.0);
    RooAbsReal *fN_CohJ_to1 = hPDFCohJ.createIntegral(fPt,NormSet(fPt),Range("fPtTo1"));
    RooAbsReal *fN_IncJ_to1 = hPDFIncJ.createIntegral(fPt,NormSet(fPt),Range("fPtTo1"));
    RooAbsReal *fN_Diss_to1 = hPDFDiss.createIntegral(fPt,NormSet(fPt),Range("fPtTo1"));  
    RooAbsReal *fN_CohP_to1 = hPDFCohP.createIntegral(fPt,NormSet(fPt),Range("fPtTo1")); 
    RooAbsReal *fN_IncP_to1 = hPDFIncP.createIntegral(fPt,NormSet(fPt),Range("fPtTo1"));  
    RooAbsReal *fN_Bkgr_to1 = hPDFBkgr.createIntegral(fPt,NormSet(fPt),Range("fPtTo1"));  
    // Number of events in the whole pt range
        // values
        Double_t N_CohJ_all_val = fN_CohJ_all->getVal()*NCohJ.getVal();
        Double_t N_IncJ_all_val = fN_IncJ_all->getVal()*NIncJ.getVal();
        Double_t N_Diss_all_val = fN_Diss_all->getVal()*NDiss.getVal();
        Double_t N_CohP_all_val = fN_CohP_all->getVal()*fDCoh*NCohJ.getVal();
        Double_t N_IncP_all_val = fN_IncP_all->getVal()*fDInc*NIncJ.getVal();
        Double_t N_Bkgr_all_val = fN_Bkgr_all->getVal()*NBkgr.getVal();
        // errors
        Double_t N_CohJ_all_err = fN_CohJ_all->getVal()*NCohJ.getError();
        Double_t N_IncJ_all_err = fN_IncJ_all->getVal()*NIncJ.getError();
        Double_t N_Diss_all_err = fN_Diss_all->getVal()*NDiss.getError();
        Double_t N_CohP_all_err = fN_CohP_all->getVal()*fDCoh*NCohJ.getError();
        Double_t N_IncP_all_err = fN_IncP_all->getVal()*fDInc*NIncJ.getError();
        Double_t N_Bkgr_all_err = fN_Bkgr_all->getVal()*NBkgr.getError();    
    // Number of events with 0.2 < pt < 2.0 GeV/c
        // values
        Double_t N_CohJ_ies_val = fN_CohJ_ies->getVal()*NCohJ.getVal();
        Double_t N_IncJ_ies_val = fN_IncJ_ies->getVal()*NIncJ.getVal();
        Double_t N_Diss_ies_val = fN_Diss_ies->getVal()*NDiss.getVal();
        Double_t N_CohP_ies_val = fN_CohP_ies->getVal()*fDCoh*NCohJ.getVal();
        Double_t N_IncP_ies_val = fN_IncP_ies->getVal()*fDInc*NIncJ.getVal();
        Double_t N_Bkgr_ies_val = fN_Bkgr_ies->getVal()*NBkgr.getVal();
        // errors
        Double_t N_CohJ_ies_err = fN_CohJ_ies->getVal()*NCohJ.getError();
        Double_t N_IncJ_ies_err = fN_IncJ_ies->getVal()*NIncJ.getError();
        Double_t N_Diss_ies_err = fN_Diss_ies->getVal()*NDiss.getError();
        Double_t N_CohP_ies_err = fN_CohP_ies->getVal()*fDCoh*NCohJ.getError();
        Double_t N_IncP_ies_err = fN_IncP_ies->getVal()*fDInc*NIncJ.getError();
        Double_t N_Bkgr_ies_err = fN_Bkgr_ies->getVal()*NBkgr.getError();    
    // Number of events with 0.2 < pt < 1.0 GeV/c
        // values
        Double_t N_CohJ_to1_val = fN_CohJ_to1->getVal()*NCohJ.getVal();
        Double_t N_IncJ_to1_val = fN_IncJ_to1->getVal()*NIncJ.getVal();
        Double_t N_Diss_to1_val = fN_Diss_to1->getVal()*NDiss.getVal();
        Double_t N_CohP_to1_val = fN_CohP_to1->getVal()*fDCoh*NCohJ.getVal();
        Double_t N_IncP_to1_val = fN_IncP_to1->getVal()*fDInc*NIncJ.getVal();
        Double_t N_Bkgr_to1_val = fN_Bkgr_to1->getVal()*NBkgr.getVal();
        // errors
        Double_t N_CohJ_to1_err = fN_CohJ_to1->getVal()*NCohJ.getError();
        Double_t N_IncJ_to1_err = fN_IncJ_to1->getVal()*NIncJ.getError();
        Double_t N_Diss_to1_err = fN_Diss_to1->getVal()*NDiss.getError();
        Double_t N_CohP_to1_err = fN_CohP_to1->getVal()*fDCoh*NCohJ.getError();
        Double_t N_IncP_to1_err = fN_IncP_to1->getVal()*fDInc*NIncJ.getError();
        Double_t N_Bkgr_to1_err = fN_Bkgr_to1->getVal()*NBkgr.getError();
    // Total fC correction (for the whole IES, with 0.2 < pt < 2.0 GeV/c)
    Double_t fC = N_CohJ_ies_val / (N_IncJ_ies_val + N_Diss_ies_val);
    Double_t denominator_err = TMath::Sqrt(TMath::Power(N_IncJ_ies_err,2) + TMath::Power(N_Diss_ies_err,2));
    Double_t fC_err = fC * TMath::Sqrt(TMath::Power((N_CohJ_ies_err/N_CohJ_ies_val),2) + TMath::Power((denominator_err/(N_IncJ_ies_val + N_Diss_ies_val)),2));
    // Total fD correction (for the whole IES, with 0.2 < pt < 2.0 GeV/c)
    Double_t fDCoh_val = N_CohP_ies_val / (N_IncJ_ies_val + N_Diss_ies_val);
    Double_t fDInc_val = N_IncP_ies_val / (N_IncJ_ies_val + N_Diss_ies_val);
    Double_t fDCoh_err = fDCoh_val * TMath::Sqrt(TMath::Power((N_CohP_ies_err/N_CohP_ies_val),2) + TMath::Power((denominator_err/(N_IncJ_ies_val + N_Diss_ies_val)),2));
    Double_t fDInc_err = fDInc_val * TMath::Sqrt(TMath::Power((N_IncP_ies_err/N_IncP_ies_val),2) + TMath::Power((denominator_err/(N_IncJ_ies_val + N_Diss_ies_val)),2));
    // Integrals of the PDFs in the pt bins
    RooAbsReal *fN_CohJ_bins[nPtBins] = { NULL };
    RooAbsReal *fN_IncJ_bins[nPtBins] = { NULL };
    RooAbsReal *fN_Diss_bins[nPtBins] = { NULL };
    RooAbsReal *fN_CohP_bins[nPtBins] = { NULL };
    RooAbsReal *fN_IncP_bins[nPtBins] = { NULL };
    RooAbsReal *fN_Bkgr_bins[nPtBins] = { NULL };
    // values
    Double_t N_CohJ_bins_val[nPtBins] = { 0 };
    Double_t N_IncJ_bins_val[nPtBins] = { 0 };
    Double_t N_Diss_bins_val[nPtBins] = { 0 };
    Double_t N_CohP_bins_val[nPtBins] = { 0 };
    Double_t N_IncP_bins_val[nPtBins] = { 0 };
    Double_t N_Bkgr_bins_val[nPtBins] = { 0 };
    // errors
    Double_t N_CohJ_bins_err[nPtBins] = { 0 };
    Double_t N_IncJ_bins_err[nPtBins] = { 0 };
    Double_t N_Diss_bins_err[nPtBins] = { 0 };
    Double_t N_CohP_bins_err[nPtBins] = { 0 };
    Double_t N_IncP_bins_err[nPtBins] = { 0 };
    Double_t N_Bkgr_bins_err[nPtBins] = { 0 };
    // coefficients
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
        fN_Bkgr_bins[i] = hPDFBkgr.createIntegral(fPt,NormSet(fPt),Range(Form("fPtBin%i",i+1)));
        // values  
        N_CohJ_bins_val[i] = fN_CohJ_bins[i]->getVal()*NCohJ.getVal();
        N_IncJ_bins_val[i] = fN_IncJ_bins[i]->getVal()*NIncJ.getVal();
        N_Diss_bins_val[i] = fN_Diss_bins[i]->getVal()*NDiss.getVal();
        N_CohP_bins_val[i] = fN_CohP_bins[i]->getVal()*fDCoh*NCohJ.getVal();
        N_IncP_bins_val[i] = fN_IncP_bins[i]->getVal()*fDInc*NIncJ.getVal();
        N_Bkgr_bins_val[i] = fN_Bkgr_bins[i]->getVal()*NBkgr.getVal();
        // errors
        N_CohJ_bins_err[i] = fN_CohJ_bins[i]->getVal()*NCohJ.getError();
        N_IncJ_bins_err[i] = fN_IncJ_bins[i]->getVal()*NIncJ.getError();
        N_Diss_bins_err[i] = fN_Diss_bins[i]->getVal()*NDiss.getError();
        N_CohP_bins_err[i] = fN_CohP_bins[i]->getVal()*fDCoh*NCohJ.getError();
        N_IncP_bins_err[i] = fN_IncP_bins[i]->getVal()*fDInc*NIncJ.getError();
        N_Bkgr_bins_err[i] = fN_Bkgr_bins[i]->getVal()*NBkgr.getError();
        // fC correction        
        fC_bins_val[i] = N_CohJ_bins_val[i] / (N_IncJ_bins_val[i] + N_Diss_bins_val[i]);
        Double_t denominator_err = TMath::Sqrt(TMath::Power(N_IncJ_bins_err[i],2) + TMath::Power(N_Diss_bins_err[i],2));
        if(N_CohJ_bins_val[i] != 0) fC_bins_err[i] = fC_bins_val[i] * TMath::Sqrt(TMath::Power((N_CohJ_bins_err[i]/N_CohJ_bins_val[i]),2) + TMath::Power((denominator_err/(N_IncJ_bins_val[i] + N_Diss_bins_val[i])),2));    
        else fC_bins_err[i] = 0.;
        fDCoh_bins_val[i] = N_CohP_bins_val[i] / (N_IncJ_bins_val[i] + N_Diss_bins_val[i]);
        fDInc_bins_val[i] = N_IncP_bins_val[i] / (N_IncJ_bins_val[i] + N_Diss_bins_val[i]);
        fDCoh_bins_err[i] = fDCoh_bins_val[i] * TMath::Sqrt(TMath::Power((N_CohP_bins_err[i]/N_CohP_bins_val[i]),2) + TMath::Power((denominator_err/(N_IncJ_bins_val[i] + N_Diss_bins_val[i])),2));
        fDInc_bins_err[i] = fDInc_bins_val[i] * TMath::Sqrt(TMath::Power((N_IncP_bins_err[i]/N_IncP_bins_val[i]),2) + TMath::Power((denominator_err/(N_IncJ_bins_val[i] + N_Diss_bins_val[i])),2));
    }
    // Print to text file
    ofstream outfile((*str + ".txt").Data());
    outfile << std::fixed << std::setprecision(2);
    outfile << Form("Dataset contains %.0f events.\n***\n", N_all);
    outfile << "In 0.0 < pt 2.0 GeV/c:\n";
    Double_t sum_all = N_CohJ_all_val + N_IncJ_all_val + N_CohP_all_val + N_IncP_all_val + N_Diss_all_val + N_Bkgr_all_val;
    outfile << "NCohJ \tNIncJ \tNCohP \tNIncP \tNDiss \tNBkgr \tSum\n";
    outfile << N_CohJ_all_val << "\t" 
            << N_IncJ_all_val << "\t" 
            << N_CohP_all_val << "\t" 
            << N_IncP_all_val << "\t" 
            << N_Diss_all_val << "\t" 
            << N_Bkgr_all_val << "\t" 
            << sum_all << "\n***\n";
    outfile << "In 0.2 < pt 2.0 GeV/c:\n";
    Double_t sum_ies = N_CohJ_ies_val + N_IncJ_ies_val + N_CohP_ies_val + N_IncP_ies_val + N_Diss_ies_val + N_Bkgr_ies_val;
    outfile << "NCohJ \tNIncJ \tNCohP \tNIncP \tNDiss \tNBkgr \tSum \tfC \tfC_err \tfDCoh \tfDC_err\tfDInc \tfDI_err\n";
    outfile << N_CohJ_ies_val << "\t" 
            << N_IncJ_ies_val << "\t" 
            << N_CohP_ies_val << "\t" 
            << N_IncP_ies_val << "\t" 
            << N_Diss_ies_val << "\t" 
            << N_Bkgr_ies_val << "\t" 
            << sum_ies << "\t";
    outfile << std::fixed << std::setprecision(4);
    outfile << fC << "\t" << fC_err << "\t" << fDCoh_val << "\t" << fDCoh_err << "\t" << fDInc_val << "\t" << fDInc_err << "\n***\n";
    outfile << "In 0.2 < pt 1.0 GeV/c:\n";
    Double_t sum_to1 = N_CohJ_to1_val + N_IncJ_to1_val + N_CohP_to1_val + N_IncP_to1_val + N_Diss_to1_val + N_Bkgr_to1_val;
    outfile << "NCohJ \tNIncJ \tNCohP \tNIncP \tNDiss \tNBkgr \tSum \n";
    outfile << std::fixed << std::setprecision(2);
    outfile << N_CohJ_to1_val << "\t" 
            << N_IncJ_to1_val << "\t" 
            << N_CohP_to1_val << "\t" 
            << N_IncP_to1_val << "\t" 
            << N_Diss_to1_val << "\t" 
            << N_Bkgr_to1_val << "\t" 
            << sum_to1 << "\n***\n";
    for(Int_t i = 0; i < nPtBins; i++){
        outfile << Form("In bin %i, (%.3f, %.3f) GeV/c:\n", i+1, ptBoundaries[i], ptBoundaries[i+1]);
        outfile << std::fixed << std::setprecision(2);
        outfile << "NCohJ \tNIncJ \tNCohP \tNIncP \tNDiss \tNBkgr \tfC \tfC_err \tfDCoh \tfDC_err\tfDInc \tfDI_err\n";
        outfile << N_CohJ_bins_val[i] << "\t" << N_IncJ_bins_val[i] << "\t" 
                << N_CohP_bins_val[i] << "\t" << N_IncP_bins_val[i] << "\t" 
                << N_Diss_bins_val[i] << "\t" << N_Bkgr_bins_val[i] << "\t";
        outfile << std::fixed << std::setprecision(4);
        outfile << fC_bins_val[i] << "\t" << fC_bins_err[i] << "\t"
                << fDCoh_bins_val[i] << "\t" << fDCoh_bins_err[i] << "\t"
                << fDInc_bins_val[i] << "\t" << fDInc_bins_err[i] << "\n***\n";
    }
    outfile << "Sum over bins:\n";
    outfile << "NCohJ \tNIncJ \tNCohP \tNIncP \tNDiss \tNBkgr\n";
    outfile << std::fixed << std::setprecision(2);
    Double_t sum_bins[6] = { 0 };
    for(Int_t i = 0; i < nPtBins; i++){
        sum_bins[0] += N_CohJ_bins_val[i];
        sum_bins[1] += N_IncJ_bins_val[i];
        sum_bins[2] += N_CohP_bins_val[i];
        sum_bins[3] += N_IncP_bins_val[i];
        sum_bins[4] += N_Diss_bins_val[i];
        sum_bins[5] += N_Bkgr_bins_val[i];
    }
    for(Int_t i = 0; i < 5; i++){
        outfile << sum_bins[i] << "\t";
    }
    outfile << sum_bins[5] << "\n***\n";
    outfile.close();
    Printf("*** Results printed to %s. ***", (*str + ".txt").Data());
    // Print the Latex table
    outfile.open((*str + "_TeX.txt").Data());
    outfile << std::fixed << std::setprecision(1);
    outfile << R"($p_\mathrm{T} \in (0.2,1.0)$~GeV/$c$ & $)"
            << N_CohJ_to1_val << R"( \pm )" << N_CohJ_to1_err << "$ \t& $"
            << N_CohP_to1_val << R"( \pm )" << N_CohP_to1_err << "$ \t& $"
            << N_IncJ_to1_val << R"( \pm )" << N_IncJ_to1_err << "$ \t& $"
            << N_IncP_to1_val << R"( \pm )" << N_IncP_to1_err << "$ \t& $"
            << N_Diss_to1_val << R"( \pm )" << N_Diss_to1_err << "$ \t& $";
    Double_t sum_to1_err = TMath::Sqrt(TMath::Power(N_CohJ_to1_err, 2)
                                        + TMath::Power(N_CohP_to1_err, 2)
                                        + TMath::Power(N_IncJ_to1_err, 2)
                                        + TMath::Power(N_IncP_to1_err, 2)
                                        + TMath::Power(N_Diss_to1_err, 2));
    outfile << sum_to1 << R"( \pm )" << sum_to1_err << "$ & $"
            << N_Bkgr_to1_val << R"($ & - \\)"
            << "\n";
    for(Int_t i = 0; i < nPtBins; i++){
        outfile << std::fixed << std::setprecision(1)
                << Form("Bin %i", i+1) << "\t& $"
                << N_CohJ_bins_val[i] << R"( \pm )" << N_CohJ_bins_err[i] << "$ \t& $"
                << N_CohP_bins_val[i] << R"( \pm )" << N_CohP_bins_err[i] << "$ \t& $"
                << N_IncJ_bins_val[i] << R"( \pm )" << N_IncJ_bins_err[i] << "$ \t& $"
                << N_IncP_bins_val[i] << R"( \pm )" << N_IncP_bins_err[i] << "$ \t& $"
                << N_Diss_bins_val[i] << R"( \pm )" << N_Diss_bins_err[i] << "$ \t& $";
        Double_t signal_sum = N_CohJ_bins_val[i] + N_CohP_bins_val[i] + N_IncJ_bins_val[i] + N_IncP_bins_val[i] + N_Diss_bins_val[i];
        Double_t signal_err = TMath::Sqrt(TMath::Power(N_CohJ_bins_err[i], 2)
                                        + TMath::Power(N_CohP_bins_err[i], 2)
                                        + TMath::Power(N_IncJ_bins_err[i], 2)
                                        + TMath::Power(N_IncP_bins_err[i], 2)
                                        + TMath::Power(N_Diss_bins_err[i], 2));
        outfile << signal_sum << R"( \pm )" << signal_err << "$ \t& $"
                << N_Bkgr_bins_val[i] << "$ & $"
                << std::fixed << std::setprecision(3)
                << fC_bins_val[i]*100 << R"( \pm )" << fC_bins_err[i]*100 << R"( $ - \\ )"
                << "\n";
    }
    outfile.close();
    Printf("*** Results printed to %s. ***", (*str + "_TeX.txt").Data());
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
    DSetData.plotOn(PtFrame,Name("DSetData"),MarkerStyle(20), MarkerSize(1.),Binning(fPtBins));
    Mod.plotOn(PtFrame,Name("Mod"),                           LineColor(215),   LineStyle(1),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));
    Mod.plotOn(PtFrame,Name("hPDFCohJ"),Components(hPDFCohJ), LineColor(222),   LineStyle(1),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));
    Mod.plotOn(PtFrame,Name("hPDFIncJ"),Components(hPDFIncJ), LineColor(kRed),  LineStyle(1),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));
    Mod.plotOn(PtFrame,Name("hPDFCohP"),Components(hPDFCohP), LineColor(222),   LineStyle(7),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));
    Mod.plotOn(PtFrame,Name("hPDFIncP"),Components(hPDFIncP), LineColor(kRed),  LineStyle(7),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));
    Mod.plotOn(PtFrame,Name("hPDFBkgr"),Components(hPDFBkgr), LineColor(kBlack),LineStyle(1),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));
    Mod.plotOn(PtFrame,Name("hPDFDiss"),Components(hPDFDiss), LineColor(15),    LineStyle(1),LineWidth(3),Normalization(sum_all,RooAbsReal::NumEvent));

    PtFrame->SetAxisRange(0,2,"X");
    // Set X axis
    PtFrame->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    PtFrame->GetXaxis()->SetTitleSize(0.05);
    PtFrame->GetXaxis()->SetLabelSize(0.05);
    // Set Y axis
    //PtFrame->GetYaxis()->SetTitle(Form("Counts per bin"));
    PtFrame->GetYaxis()->SetTitle(Form("d#it{N}/d#it{p}_{T}"));
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
    if(!isBkgData) str = new TString(Form("%sPtFit_BkgMC_Binning%i", OutputPtFit.Data(),BinningOpt));
    else str = new TString(Form("%sPtFit_BkgData_Binning%i", OutputPtFit.Data(),BinningOpt));
    cCM->Print((*str + "_CM.pdf").Data());
    cCM->Print((*str + "_CM.png").Data());
    cPt->Print((*str + ".pdf").Data());
    cPt->Print((*str + ".png").Data());

    return;
}