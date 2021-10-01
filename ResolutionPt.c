// ResolutionPt.c
// David Grund, Sep 28, 2021

// root headers
#include "TFile.h"
#include "TList.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
// roofit headers
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooExtendPdf.h"
#include "RooCBShape.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
// my headers
#include "AnalysisManager.h"

using namespace RooFit;

Double_t res;
Double_t resLow = -0.5;
Double_t resUpp = 0.5;
Int_t nBins = 50;
Double_t BinSize = (resUpp - resLow) / (Double_t)nBins;

void FitResInBins(TH1D *h, Int_t iBin, Int_t FitOpt);
void CalculateResPerBin();

void ResolutionPt(){

    SetPtBinning();

    //CalculateResPerBin();

    // Load data file
    TFile *file = TFile::Open(Form("Trees/ResolutionPt/Trees_%ibins.root", nPtBins), "read");
    if(file) Printf("File %s loaded.", file->GetName());    
    else {
        Printf("Input file missing. Terminating...");
        return;
    }

    TList *list = (TList*) file->Get("Output");
    if(list) Printf("List %s loaded.", list->GetName()); 
    Printf("\n***\n");

    TTree *tResBins[nPtBins] = { NULL };
    TH1D *hResBins[nPtBins] = { NULL };

    // Loop over bins
    for(Int_t iBin = 0; iBin < nPtBins; iBin++){
        tResBins[iBin] = (TTree*)list->FindObject(Form("tResBin%i", iBin+1));
        hResBins[iBin] = (TH1D*)list->FindObject(Form("hResBin%i", iBin+1));
        if(hResBins[iBin]){
            Printf("Tree %s found.", tResBins[iBin]->GetName());
            Printf("Tree %s contains %lli entries.", tResBins[iBin]->GetName(), tResBins[iBin]->GetEntries());
            Printf("Histogram %s found.", hResBins[iBin]->GetName());
            FitResInBins(hResBins[iBin], iBin+1, 2);
            Printf("\n***\n");
        } 
    }

    return;
}

void FitResInBins(TH1D *h, Int_t iBin, Int_t FitOpt){

    RooRealVar res("res","res",resLow,resUpp);

    RooDataHist DHisData("DHisData","DHisData",res,h);
    // Calculate the number of entries
    Double_t N_all = 0;
    for(Int_t i = 1; i <= h->GetNbinsX(); i++){
        N_all += h->GetBinContent(i);
    }
    Printf("Histogram contains %.0f entries in %i bins.", N_all, DHisData.numEntries());

    RooRealVar mean_L("mean_L","mean_L",0.,resLow,resUpp);
    RooRealVar sigma_L("sigma_L","sigma_L",0.015,0.001,1);
    RooRealVar N("N","number of signal events", N_all, N_all*0.9, N_all*1.1);
    // Gaussian fit:
    RooGaussian gauss("gauss","gaussian PDF",res,mean_L,sigma_L);
    RooExtendPdf gauss_ext("gauss_ext","extended gaussian PDF", gauss, N);
    // Single CB fit:
    RooRealVar a_L("a_L","a from CB_L",1.,0.,10.);
    RooRealVar n_L("n_L","n from CB_L",1.,0.,20.);
    RooCBShape CB("CB","CB PDF",res,mean_L,sigma_L,a_L,n_L);
    RooExtendPdf CB_ext("CB_ext","extended CB PDF", CB, N);
    // Double CB fit:
    RooRealVar a_R("a_R","a from CB_R",1.,-10.,0.);
    RooRealVar n_R("n_R","n from CB_R",1.,0.,20.);
    RooGenericPdf mean_R("mean_R","mean_R","mean_L",RooArgSet(mean_L));
    RooGenericPdf sigma_R("sigma_R","sigma_R","sigma_L",RooArgSet(sigma_L));
    RooCBShape CB_L("CB_L","CB_L",res,mean_L,sigma_L,a_L,n_L);
    RooCBShape CB_R("CB_R","CB_R",res,mean_R,sigma_R,a_R,n_R);
    RooRealVar frac("frac","fraction of CBs",0.5);
    RooAddPdf DoubleSidedCB("DoubleSidedCB","DoubleSidedCB",RooArgList(CB_L,CB_R),RooArgList(frac));
    RooExtendPdf DSCB_ext("DSCB_ext","Extended DSCB PDF", DoubleSidedCB, N);
    // Fit:
    RooFitResult *ResFit = NULL;
    if(FitOpt == 0) RooFitResult *ResFit = gauss_ext.fitTo(DHisData,Save(),Extended(kTRUE),Range(resLow,resUpp));
    if(FitOpt == 1) RooFitResult *ResFit = CB_ext.fitTo(DHisData,Save(),Extended(kTRUE),Range(resLow,resUpp));
    if(FitOpt == 2) RooFitResult *ResFit = DSCB_ext.fitTo(DHisData,Save(),Extended(kTRUE),Range(resLow,resUpp));

    TCanvas *c = new TCanvas("c","c",900,600);
    //c->SetLogy();

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");

    c->SetTopMargin(0.06);
    c->SetBottomMargin(0.14);
    c->SetRightMargin(0.03);
    c->SetLeftMargin(0.10);

    RooPlot *frame = res.frame(Title("Resolution fit"));
    DHisData.plotOn(frame, Name("DHisData"), MarkerStyle(20), MarkerSize(1.));
    if(FitOpt == 0) gauss_ext.plotOn(frame, Name("frame"), LineColor(215), LineStyle(1), LineWidth(3));
    if(FitOpt == 1) CB_ext.plotOn(frame, Name("frame"), LineColor(215), LineStyle(1), LineWidth(3));
    if(FitOpt == 2) DSCB_ext.plotOn(frame, Name("frame"), LineColor(215), LineStyle(1), LineWidth(3));
    // Vertical axis
    frame->GetYaxis()->SetTitle(Form("Counts per %.3f [-]", BinSize));
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleOffset(0.95);
    frame->GetYaxis()->SetLabelSize(0.05);
    frame->GetYaxis()->SetLabelOffset(0.01);
    frame->GetYaxis()->SetMaxDigits(3);
    // Horizontal axis
    frame->GetXaxis()->SetTitle("(#it{p}_{T}^{rec} #minus #it{p}_{T}^{gen})/#it{p}_{T}^{gen} [-]");
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetTitleOffset(1.2);
    frame->GetXaxis()->SetLabelSize(0.05);
    frame->GetXaxis()->SetLabelOffset(0.01);
    frame->GetXaxis()->SetDecimals(1);
    frame->Draw();

    TLegend *l1 = new TLegend(0.02,0.55,0.5,0.92);
    l1->AddEntry((TObject*)0,Form("ALICE Simulation"),""); 
    l1->AddEntry((TObject*)0,Form("Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV"),"");
    l1->AddEntry((TObject*)0,Form("inc J/#psi #rightarrow #mu^{+}#mu^{-}"),"");
    l1->AddEntry((TObject*)0,Form("#it{p}_{T} #in (%.3f, %.3f) GeV/c", ptBoundaries[iBin-1], ptBoundaries[iBin]),"");
    l1->AddEntry((TObject*)0,Form("#mu = %.4f #pm %.4f", mean_L.getVal(), mean_L.getError()),"");
    l1->AddEntry((TObject*)0,Form("#sigma = %.4f #pm %.4f", sigma_L.getVal(), sigma_L.getError()),"");
    l1->SetTextSize(0.05);
    l1->SetBorderSize(0);
    l1->SetFillStyle(0);
    l1->Draw();

    ///*
    TLegend *l2 = new TLegend(0.62,0.65,0.9,0.92);
    l2->AddEntry((TObject*)0,Form("#alpha_{L} = %.2f #pm %.2f", a_L.getVal(), a_L.getError()),"");
    l2->AddEntry((TObject*)0,Form("#alpha_{R} = %.2f #pm %.2f", (-1)*a_R.getVal(), a_R.getError()),"");
    l2->AddEntry((TObject*)0,Form("#it{n}_{L} = %.2f #pm %.2f", n_L.getVal(), n_L.getError()),"");
    l2->AddEntry((TObject*)0,Form("#it{n}_{R} = %.2f #pm %.2f", n_R.getVal(), n_R.getError()),"");
    l2->SetTextSize(0.05);
    l2->SetBorderSize(0);
    l2->SetFillStyle(0);
    l2->Draw();
    //*/
    
    TCanvas *cLog = new TCanvas("cLog","cLog",900,600);
    cLog->SetLogy();
    cLog->SetTopMargin(0.06);
    cLog->SetBottomMargin(0.14);
    cLog->SetRightMargin(0.03);
    cLog->SetLeftMargin(0.10);
    frame->Draw();
    l1->Draw();
    l2->Draw();

    TString *str = new TString(Form("Results/ResolutionPt/%ibins/bin%i", nPtBins, iBin));
    //TString *str = new TString("Results/ResolutionPt/gauss_fit");
    c->Print((*str + ".pdf").Data());
    c->Print((*str + ".png").Data());
    cLog->Print((*str + "_log.pdf").Data());
    cLog->Print((*str + "_log.png").Data());

    return;
}

void CalculateResPerBin(){

    TFile *file = TFile::Open("Trees/AnalysisDataMC/AnalysisResults_MC_kIncohJpsiToMu.root", "read");
    if(file) Printf("File %s loaded.", file->GetName());

    TTree *tRec = dynamic_cast<TTree*> (file->Get("AnalysisOutput/fTreeJPsiMCRec"));
    if(tRec) Printf("Input tree loaded.");

    ConnectTreeVariablesMCRec(tRec);

    TFile *f = new TFile(Form("Trees/ResolutionPt/Trees_%ibins.root", nPtBins),"RECREATE");
    TList *l = new TList();

    TH1D *hResBins[nPtBins] = { NULL };
    TTree *tResBins[nPtBins] = { NULL };

    for(Int_t i = 0; i < nPtBins; i++){
        hResBins[i] = new TH1D(Form("hResBin%i", i+1), Form("hResBin%i", i+1), nBins, resLow, resUpp);
        tResBins[i] = new TTree(Form("tResBin%i", i+1), Form("tResBin%i", i+1));
        tResBins[i]->Branch("res", &res, "res/D");
        l->Add(hResBins[i]);
        l->Add(tResBins[i]);
    }

    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < tRec->GetEntries(); iEntry++){
        tRec->GetEntry(iEntry);
        for(Int_t iBin = 0; iBin < nPtBins; iBin++){
            if(EventPassedMCRec(0,4,iBin+1) && fPtGen > ptBoundaries[iBin] && fPtGen < ptBoundaries[iBin+1]){
                res = (fPt - fPtGen) / fPtGen;
                hResBins[iBin]->Fill(res);
                tResBins[iBin]->Fill();
            } 
        }
        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }
    Printf("Done.");

    l->Write("Output", TObject::kSingleKey);
    f->ls();

    return;
}