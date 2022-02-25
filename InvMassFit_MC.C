// InvMassFit_MC.c
// David Grund, Sep 17, 2021
// To perform fit of the invariant mass distribution of MC data

// cpp headers
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
// root headers
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TStyle.h" // gStyle
// roofit headers
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooGenericPdf.h"
#include "RooBinning.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
// my headers
#include "AnalysisManager.h"

using namespace RooFit;

// Main functions
void DoInvMassFitMainMC(Int_t opt, Bool_t pass3);
// Support functions
void DrawCorrelationMatrix(TCanvas *cCorrMat, RooFitResult* fResFit);
void SetCanvas(TCanvas *c, Bool_t bLogScale);
void PrepareMCTree(Bool_t pass3);

void InvMassFit_MC(){

    //PrepareMCTree(kFALSE);
    //PrepareMCTree(kTRUE);

    Bool_t pass3 = kFALSE;

    Bool_t main_fits = kFALSE;
    if(main_fits){
        DoInvMassFitMainMC(0, pass3);
        DoInvMassFitMainMC(1, pass3);
        DoInvMassFitMainMC(2, pass3);
        DoInvMassFitMainMC(3, pass3);
    }
    // bins:
    Bool_t bins = kTRUE;
    if(bins){
        SetPtBinning(); // PtBinning method must be chosen in PtBinsManager.h
        DoInvMassFitMainMC(4, pass3);
        DoInvMassFitMainMC(5, pass3);
        DoInvMassFitMainMC(6, pass3);
        DoInvMassFitMainMC(7, pass3);
        if(nPtBins == 5){
            DoInvMassFitMainMC(8, pass3);
        }
    }
    /*
    // debug bin 3 (value of n_R)
    Bool_t debug = kFALSE;
    if(debug){
        SetPtBinning();
        DoInvMassFitMainMC(6);
    }
    */

    Printf("Done.");

    return;
}

void DoInvMassFitMainMC(Int_t opt, Bool_t pass3){
    // Fit the invariant mass distribution using Double-sided CB function
    // Peak corresponding to psi(2s) excluded

    // Cuts:
    char fStrReduce[120];
    Double_t fPtCut     = -999;
    Double_t fPtCutLow  = -999;
    Double_t fPtCutUpp  = -999;
    Double_t fYCut      = 0.80;
    Double_t fMCutLow   = 2.90;
    Double_t fMCutUpp   = 3.30;

    switch(opt){
        case 0: // Incoherent-enriched sample
            fPtCut = 0.20;
            sprintf(fStrReduce,"abs(fY)<%f && fPt>%f && fM>%f && fM<%f",fYCut,fPtCut,fMCutLow,fMCutUpp);
            break;
        case 1: // Coherent-enriched sample
            fPtCut = 0.20;
            sprintf(fStrReduce,"abs(fY)<%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCut,fMCutLow,fMCutUpp);
            break;
        case 2: // Total sample (pt < 2.0 GeV/c)
            fPtCut = 2.00;
            sprintf(fStrReduce,"abs(fY)<%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCut,fMCutLow,fMCutUpp);
            break;
        case 3: // Sample with pt from 0.2 to 1 GeV/c 
            fPtCutLow = 0.20;
            fPtCutUpp = 1.00;
            sprintf(fStrReduce,"abs(fY)<%f && fPt>%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCutLow,fPtCutUpp,fMCutLow,fMCutUpp);
            break;
        case 4: // pt bin 1
            fPtCutLow = ptBoundaries[0];
            fPtCutUpp = ptBoundaries[1];
            sprintf(fStrReduce,"abs(fY)<%f && fPt>%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCutLow,fPtCutUpp,fMCutLow,fMCutUpp);
            break;
        case 5: // pt bin 2
            fPtCutLow = ptBoundaries[1];
            fPtCutUpp = ptBoundaries[2];
            sprintf(fStrReduce,"abs(fY)<%f && fPt>%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCutLow,fPtCutUpp,fMCutLow,fMCutUpp);
            break;
        case 6: // pt bin 3
            fPtCutLow = ptBoundaries[2];
            fPtCutUpp = ptBoundaries[3];
            sprintf(fStrReduce,"abs(fY)<%f && fPt>%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCutLow,fPtCutUpp,fMCutLow,fMCutUpp);
            break;
        case 7: // pt bin 4
            fPtCutLow = ptBoundaries[3];
            fPtCutUpp = ptBoundaries[4];
            sprintf(fStrReduce,"abs(fY)<%f && fPt>%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCutLow,fPtCutUpp,fMCutLow,fMCutUpp);
            break;
        case 8: // pt bin 5
            fPtCutLow = ptBoundaries[4];
            fPtCutUpp = ptBoundaries[5];
            sprintf(fStrReduce,"abs(fY)<%f && fPt>%f && fPt<%f && fM>%f && fM<%f",fYCut,fPtCutLow,fPtCutUpp,fMCutLow,fMCutUpp);
            break;
    }

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
    TString str_file = "";
    if(!pass3)  str_file = "Trees/InvMassFit_MC/InvMassFit_MC_pass1.root";
    else        str_file = "Trees/InvMassFit_MC/InvMassFit_MC_pass3.root";
    TFile *f_in = new TFile(str_file.Data()); 
    TTree *t_in = NULL;
    if(opt == 0 || opt == 3 || opt == 4 || opt == 5 || opt == 6 || opt == 7 || opt == 8){
        f_in->GetObject("tIncEnrSample",t_in);
    } else if(opt == 1){
        f_in->GetObject("tCohEnrSample",t_in);
    } else if(opt == 2){
        f_in->GetObject("tMixedSample",t_in);
    }

    RooDataSet *fDataIn = new RooDataSet("fDataIn", "fDataIn", RooArgSet(fM,fY,fPt), Import(*t_in));
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
    RooRealVar n_R("n_{R}","n_{R}",8.,0,30);

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
    // Draw Correlation Matrix
    TCanvas *cCorrMat = new TCanvas("cCorrMat","cCorrMat",700,600);
    DrawCorrelationMatrix(cCorrMat, fResFit);

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
    //leg2->AddEntry((TObject*)0,"#bf{This thesis}","");
    // Print the pt cut
    if(opt == 0){
        leg2->AddEntry((TObject*)0,Form("#it{p}_{T} > %.2f GeV/#it{c}", fPtCut),"");
    } else if(opt == 1 || opt == 2){
        leg2->AddEntry((TObject*)0,Form("#it{p}_{T} < %.2f GeV/#it{c}", fPtCut),"");
    } else if(opt == 3 || opt == 4 || opt == 5 || opt == 6 || opt == 7 || opt == 8){
        leg2->AddEntry((TObject*)0,Form("#it{p}_{T} #in (%.2f,%.2f) GeV/#it{c}", fPtCutLow,fPtCutUpp),"");
    }
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

    // Prepare path
    TString str = "";
    if(!pass3)  str = "Results/InvMassFit_MC/pass1/";
    else        str = "Results/InvMassFit_MC/pass3/";

    switch(opt){
        case 0:
            str = str + "inc/inc";
            break;
        case 1:
            str = str + "coh/coh";
            break;
        case 2:
            str = str + "all/all";
            break;
        case 3:
            str = str + "allbins/allbins";
            break;
        case 4:
            if(nPtBins == 4) str = str + "4bins/bin1"; 
            if(nPtBins == 5) str = str + "5bins/bin1"; 
            break;
        case 5:
            if(nPtBins == 4) str = str + "4bins/bin2"; 
            if(nPtBins == 5) str = str + "5bins/bin2"; 
            break;
        case 6:
            if(nPtBins == 4) str = str + "4bins/bin3"; 
            if(nPtBins == 5) str = str + "5bins/bin3"; 
            break;
        case 7:
            if(nPtBins == 4) str = str + "4bins/bin4"; 
            if(nPtBins == 5) str = str + "5bins/bin4"; 
            break;
        case 8:
            if(nPtBins == 5) str = str + "5bins/bin5"; 
            break;
    }
    // Print the plots
    cHist->Print((str + ".pdf").Data());
    cHist->Print((str + ".png").Data());
    cHistLog->Print((str + "_log.pdf").Data());
    cHistLog->Print((str + "_log.png").Data());
    cCorrMat->Print((str + "_cm.pdf").Data());
    cCorrMat->Print((str + "_cm.png").Data());  
    
    // Print the values of alpha and n to txt output files
    ofstream outfile((str + ".txt").Data());
    outfile << "alpha_L \t" << alpha_L.getVal() << "\t" << alpha_L.getError() << "\n";
    outfile << "alpha_R \t" << alpha_R.getVal() << "\t" << alpha_R.getError() << "\n";
    outfile << "n_L \t" << n_L.getVal() << "\t" << n_L.getError() << "\n";
    outfile << "n_R \t" << n_R.getVal() << "\t" << n_R.getError() << "\n";
    outfile.close();
    Printf("*** Results printed to %s.***", (str + ".txt").Data());

    return;
}

void DrawCorrelationMatrix(TCanvas *cCorrMat, RooFitResult* fResFit){

    cCorrMat->SetTopMargin(0.05);
    cCorrMat->SetRightMargin(0.12);
    cCorrMat->SetLeftMargin(0.12);

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");

    TH2* hCorr = fResFit->correlationHist();

    hCorr->GetXaxis()->SetBinLabel(7,"#sigma");
    hCorr->GetYaxis()->SetBinLabel(7,"#sigma");

    hCorr->SetMarkerSize(2.0);
    hCorr->GetXaxis()->SetLabelSize(0.08); //0.049
    hCorr->GetYaxis()->SetLabelSize(0.08);
    hCorr->Draw("colz,text");

    return;
}

void SetCanvas(TCanvas *c, Bool_t bLogScale){

    if(bLogScale == kTRUE) c->SetLogy();
    c->SetTopMargin(0.055);
    c->SetBottomMargin(0.12);
    c->SetRightMargin(0.03);
    c->SetLeftMargin(0.11);

    return;
}

void PrepareMCTree(Bool_t pass3){

    TString str_f_coh = "";
    TString str_f_inc = "";
    TString str_tree = "";
    TString str_fout = "";
    if(!pass3){
        str_f_coh = "Trees/AnalysisDataMC_pass1/AnalysisResults_MC_kCohJpsiToMu.root"; 
        str_f_inc = "Trees/AnalysisDataMC_pass1/AnalysisResults_MC_kIncohJpsiToMu.root";
        str_tree = "AnalysisOutput/fTreeJPsiMCRec"; 
        str_fout = "Trees/InvMassFit_MC/InvMassFit_MC_pass1.root";
    } else {
        str_f_coh = "Trees/AnalysisDataMC_pass3/AnalysisResults_MC_kCohJpsiToMu.root"; 
        str_f_inc = "Trees/AnalysisDataMC_pass3/AnalysisResults_MC_kIncohJpsiToMu.root";
        str_tree = "AnalysisOutput/fTreeJpsi"; 
        str_fout = "Trees/InvMassFit_MC/InvMassFit_MC_pass3.root";
    }

    // kCohJpsiToMu
    TFile *f_in_coh = TFile::Open(str_f_coh.Data(), "read");
    if(f_in_coh) Printf("Input data loaded.");

    TTree *t_in_coh = dynamic_cast<TTree*> (f_in_coh->Get(str_tree.Data()));
    if(t_in_coh) Printf("Input tree loaded.");

    // kIncohJpsiToMu
    TFile *f_in_inc = TFile::Open(str_f_inc.Data(), "read");
    if(f_in_inc) Printf("Input data loaded.");

    TTree *t_in_inc = dynamic_cast<TTree*> (f_in_inc->Get(str_tree.Data()));
    if(t_in_inc) Printf("Input tree loaded.");

    // Create new MC trees with applied cuts
    TFile f_out(str_fout.Data(),"RECREATE");

    TTree *tIncEnrSample = new TTree("tIncEnrSample", "tIncEnrSample");
    tIncEnrSample->Branch("fPt", &fPt, "fPt/D");
    tIncEnrSample->Branch("fM", &fM, "fM/D");
    tIncEnrSample->Branch("fY", &fY, "fY/D");

    TTree *tCohEnrSample = new TTree("tCohEnrSample", "tCohEnrSample");
    tCohEnrSample->Branch("fPt", &fPt, "fPt/D");
    tCohEnrSample->Branch("fM", &fM, "fM/D");
    tCohEnrSample->Branch("fY", &fY, "fY/D");

    TTree *tMixedSample = new TTree("tMixedSample", "tMixedSample");
    tMixedSample->Branch("fPt", &fPt, "fPt/D");
    tMixedSample->Branch("fM", &fM, "fM/D");
    tMixedSample->Branch("fY", &fY, "fY/D");

    ConnectTreeVariablesMCRec(t_in_coh, pass3);

    Printf("%lli entries found in the tree.", t_in_coh->GetEntries());
    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < t_in_coh->GetEntries(); iEntry++){
        t_in_coh->GetEntry(iEntry);
        if(EventPassedMCRec(0, 2, -1, pass3)){
            tCohEnrSample->Fill();
            tMixedSample->Fill();
        } 

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }

    ConnectTreeVariablesMCRec(t_in_inc, pass3);

    Printf("%lli entries found in the tree.", t_in_inc->GetEntries());
    nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < t_in_inc->GetEntries(); iEntry++){
        t_in_inc->GetEntry(iEntry);
        if(EventPassedMCRec(0, 2, -1, pass3)){
            tIncEnrSample->Fill();
            tMixedSample->Fill();
        } 

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }

    f_out.Write("",TObject::kWriteDelete);

    return;
}