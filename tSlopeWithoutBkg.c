// tSlopeWithoutBkg.c
// David Grund, Sep 15, 2021

// cpp headers
#include <fstream> // print output to txt file
// root headers
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TStyle.h" // gStyle
#include "TMath.h"
// roofit headers
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooGenericPdf.h"
#include "RooBinning.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
// my headers
#include "AnalysisManager.h"

using namespace RooFit;

// Main function
void DoInvMassFitMain(Double_t fPtCutLow, Double_t fPtCutUpp, Bool_t save = kFALSE, Int_t bin = -1);
// Support functions
void SetCanvas(TCanvas *c, Bool_t bLogScale);

Double_t YieldJpsi = 0;

void tSlopeWithoutBkg(){

    // PtBinning "Method 4"
    // Using Guillermo's macro after subtracting the background first using
    // inv mass fits in 10 predefined pt bins by 0.08 GeV (pt_step)

    const Double_t pt_step = 0.08;
    const Double_t pt_min = 0.20;
    const Double_t pt_max = 1.00;
    Int_t n_bins = (pt_max - pt_min) / pt_step;
    Printf("Number of bins: %i", n_bins);
    Int_t iBin = 1;

    // Print the results
    TString name = "PtBinning/tSlopeWithoutBkg/output.txt";
    ofstream outfile(name.Data());

    TH1D *hYieldJpsi = new TH1D("hYieldJpsi", "hYieldJpsi", n_bins, pt_min, pt_max);

    Double_t pt = pt_min;
    while(iBin <= n_bins){
        DoInvMassFitMain(pt, pt + pt_step, kTRUE, iBin);
        outfile << Form("Bin %i: (%.3f, %.3f), %.3f\n", iBin, pt, pt + pt_step, YieldJpsi);
        hYieldJpsi->SetBinContent(iBin, YieldJpsi);
        pt += pt_step;
        iBin++;
    }

    outfile.close();
    Printf("*** Results printed to %s.***", name.Data());

    TCanvas *cHist = new TCanvas("cHist", "cHist", 900, 600);
    cHist->SetLogy();    
    // Vertical axis
    hYieldJpsi->GetYaxis()->SetTitle(Form("Counts per %.0f MeV/#it{c}^{2}", pt_step*1000));
    // Horizontal axis
    hYieldJpsi->GetXaxis()->SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    hYieldJpsi->GetXaxis()->SetDecimals(1);
    // Draw the histogram
    hYieldJpsi->Draw("E0");

    // Fit the histogram with an exponential
    RooRealVar fPt("fPt", "fPt", pt_min, pt_max);
    RooRealVar slope("slope","slope",2.,1.,10.); // positive (extra minus in the exp)
    // Define RooFit datasets
    RooDataHist data("data","data",fPt,hYieldJpsi);
    // Define RooFit PDF
    RooGenericPdf ExpPdf("ExpPdf","exp(-fPt*slope)",RooArgSet(fPt,slope));
    RooFitResult* fResFit = ExpPdf.fitTo(data,Extended(kTRUE),Range(pt_min,pt_max),Save());

    TCanvas *cFit = new TCanvas("cFit","cFit",900,600);
    cFit->SetTopMargin(0.055);
    cFit->SetBottomMargin(0.12);
    cFit->SetRightMargin(0.03);
    cFit->SetLeftMargin(0.11);    
    cFit->SetLogy();

    // Roo frame
    RooPlot* fFrame = fPt.frame(Title("Fit of |p_T| distribution")); 
    data.plotOn(fFrame,Name("data"),MarkerStyle(20),MarkerSize(1.));
    ExpPdf.plotOn(fFrame,Name("ExpPdf"),Components(ExpPdf),LineColor(kBlue),LineStyle(kDashed),LineWidth(3));
    // Vertical axis
    fFrame->GetYaxis()->SetTitle(Form("Counts per %.0f MeV/#it{c}",pt_step*1000));
    fFrame->GetYaxis()->SetTitleOffset(1.0);
    fFrame->GetYaxis()->SetTitleSize(0.05);
    fFrame->GetYaxis()->SetLabelSize(0.05);
    // Horizontal axis
    fFrame->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fFrame->GetXaxis()->SetTitleSize(0.05);
    fFrame->GetXaxis()->SetLabelSize(0.05);
    // Finally plot it
    fFrame->Draw();
    // Legend
    TLegend *l = new TLegend(0.14,0.18,0.5,0.38);  
    l->AddEntry("ExpPdf","fit: exp(-#it{b}#it{p}_{T})","L");
    l->AddEntry((TObject*)0,Form("#it{b} = %.3f GeV^{-1}#it{c}", slope.getVal()),"");  
    l->AddEntry("data","data","P");
    l->SetTextSize(0.048);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->Draw();

    // Print the canvas
    cFit->Print("PtBinning/tSlopeWithoutBkg/exp_fit.png");

    // Calculate the pt binning
    Double_t b = slope.getVal(); // slope
    Double_t ptBoundariesNew[nPtBins+1] = { 0 };
    ptBoundariesNew[0] = pt_min;
    
    // Initialise
    Double_t pt_1 = pt_min;
    
    // Fraction of the integral occupied by each bin
    Double_t f = 1.0/nPtBins;

    // Go over the bins
    for(Int_t i = 1; i <= nPtBins; i++){
        // Exponentials
        Double_t e_i = TMath::Exp(-b*pt_min);
        Double_t e_f = TMath::Exp(-b*pt_max);
        // Current pt boundary
        Double_t e_1 = TMath::Exp(-b*pt_1);

        // Integral in the full range
        Double_t e_int = e_i - e_f; // the exact integral has an extra 1/b
    
        // New value of t
        Double_t pt_n = -TMath::Log(e_1-f*e_int)/b;

        ptBoundariesNew[i] = pt_n;

        // Prepare for the next bin
        pt_1 = pt_n;
    }

    // Make inv mass fits in defined bins
    for(Int_t i = 1; i <= nPtBins; i++){
        DoInvMassFitMain(ptBoundariesNew[i-1], ptBoundariesNew[i], kTRUE, 100+i);
    }    

    // Print the values
    for(Int_t i = 1; i <= nPtBins; i++){
        Printf("Bin %i: (%.3f, %.3f)", i, ptBoundariesNew[i-1], ptBoundariesNew[i]);
    }

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

    //fM.setBinning(binM);

    // Get the data trees
    TFile *fFileIn = new TFile("Trees/InvMassFit/InvMassFit.root"); 
    TTree *fTreeIn = NULL;
    fFileIn->GetObject("tIncEnrSample",fTreeIn);
        
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

    TString* path = new TString("Results/InvMassFitMC/");
    path->Append("inc_doubleCB.txt");

    ifstream fTxtFileIn;
    fTxtFileIn.open(path->Data());
    if(fTxtFileIn.fail()){
        Printf("\n");
        Printf("*** Warning! ***");
        Printf("*** MC values for tail parameters not found. Terminating... *** \n");
        return;
    } else {
        Int_t i_line = 0;
        while(!fTxtFileIn.eof()){
            fTxtFileIn >> name >> values[i_line] >> errors[i_line];
            i_line++;
        }
        fTxtFileIn.close();
    }
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
    YieldJpsi = N_Jpsi_out[0];

    // ##########################################################
    // Plot the results
    // Draw histogram with fit results
    TCanvas *cHist = new TCanvas("cHist","cHist",800,600);
    SetCanvas(cHist,kFALSE);

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
        str = new TString(Form("PtBinning/tSlopeWithoutBkg/bin%i", bin));
        // Print the plots
        //cHist->Print((*str + ".pdf").Data());
        cHist->Print((*str + ".png").Data()); 
    }

    return;
}

void SetCanvas(TCanvas *c, Bool_t bLogScale){

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");

    if(bLogScale == kTRUE) c->SetLogy();
    c->SetTopMargin(0.055);
    c->SetBottomMargin(0.12);
    c->SetRightMargin(0.03);
    c->SetLeftMargin(0.11);

    return;
}