// tSlope.c
// David Grund, Sep 10, 2021
// To calculate the slope of the |t| distribution of the measured data
// in the pt range 0.2 to 1.0 GeV/c

// root headers
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TBrowser.h"
// roofit headers
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooGenericPdf.h"
#include "RooBinning.h"
// my headers
#include "fTreeJPsiManager.h"

using namespace RooFit;

Double_t t;

void PrepareData();
void CalculateSlope();

void tSlope(){

    //PrepareData();

    CalculateSlope();

    return;
}

void CalculateSlope(){

    TFile *fFileIn = TFile::Open("Trees/tSlope/tSlope.root", "read");
    if(fFileIn) Printf("Input data loaded.");

    TTree *fTreeIn = dynamic_cast<TTree*> (fFileIn->Get("tSlope"));
    if(fTreeIn) Printf("Input tree loaded.");

    fTreeIn->SetBranchAddress("t", &t);

    Int_t n_bins = 100;
    Double_t pt_i = 0.2;
    Double_t pt_f = 1.0;
    Double_t t_i = pt_i*pt_i;
    Double_t t_f = pt_f*pt_f;
    Double_t bin_size_double = (t_f - t_i) * 1000000 / n_bins;
    bin_size_double = bin_size_double + 0.5;
    Int_t bin_size = (Int_t)bin_size_double; // in MeV^2

    Printf("Low t boundary: %f", t_i);
    Printf("Upp t boundary: %f", t_f);

    TCanvas *cHist = new TCanvas("cHist", "cHist", 900, 600);
    cHist->SetLogy();

    TH1D *hist = new TH1D("hist", "hist in Mandelstam t", n_bins, t_i, t_f); 
    fTreeIn->Draw("t >> hist");

    // Fitting in RooFit
    char fStrReduce[120];
    sprintf(fStrReduce,"t > %f && t < %f",t_i,t_f);
    RooRealVar t("t","t",t_i,t_f);
    RooRealVar slope("slope","slope",2.,1.,10.); // positive (extra minus in the exp)
    // Define RooFit datasets
    RooDataSet *fDataIn = new RooDataSet("fDataIn","fDataIn", RooArgSet(t), Import(*fTreeIn));
    RooAbsData* fDataSet = fDataIn->reduce(fStrReduce);
    // Print the number of entries in the dataset
    Int_t nEvents = fDataSet->numEntries();
    Printf("\n");
    Printf("*** Number of events in the dataset: %i ***\n", nEvents);
    // Define RooFit PDF
    RooGenericPdf ExpPdf("ExpPdf","exp(-t*slope)",RooArgSet(t,slope));
    RooFitResult* fResFit = ExpPdf.fitTo(*fDataSet,Extended(kTRUE),Range(t_i,t_f),Save());
    // Set Roo Binning
    RooBinning binT(n_bins,t_i,t_f);
    // Plot the results
    TCanvas *cFit = new TCanvas("cFit","cFit",900,600);
    cFit->SetLogy();
    RooPlot* fFrame = t.frame(Title("Fit of |t| distribution")); 
    fDataSet->plotOn(fFrame,Name("fDataSet"),Binning(binT),MarkerStyle(20),MarkerSize(1.));
    ExpPdf.plotOn(fFrame,Name("ExpPdf"),Components(ExpPdf),LineColor(kBlue),LineStyle(kDashed),LineWidth(3));
    // Vertical axis
    fFrame->GetYaxis()->SetTitle(Form("Counts per %i MeV^{2}",bin_size));
    // Horizontal axis
    fFrame->GetXaxis()->SetTitle("|#it{t}| (GeV^{2})");
    // Finally plot it
    fFrame->Draw();
    // Save the result
    cFit->Print("PtBinning/data_slope.png");
    cFit->Print("PtBinning/data_slope.pdf");


    return;
}

void PrepareData(){

    TFile *fFileIn = TFile::Open("Trees/AnalysisData/AnalysisResultsLHC18qrMerged.root", "read");
    if(fFileIn) Printf("Input data loaded.");

    TTree *fTreeIn = dynamic_cast<TTree*> (fFileIn->Get("AnalysisOutput/fTreeJPsi"));
    if(fTreeIn) Printf("Input tree loaded.");

    ConnectTreeVariables(fTreeIn);

    // Create new data tree with applied cuts
    TFile fFileOut("Trees/tSlope/tSlope.root","RECREATE");

    TTree *tSlope = new TTree("tSlope", "tSlope");
    tSlope->Branch("t", &t, "t/D");

    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < fTreeIn->GetEntries(); iEntry++){
        fTreeIn->GetEntry(iEntry);
        if(EventPassed(0,-1) && fPt > 0.0 && fPt < 1.0){
            t = fPt*fPt;
            tSlope->Fill();
        }
        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }

    fFileOut.Write("",TObject::kWriteDelete);

    new TBrowser;

    return;
}