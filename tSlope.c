// tSlope.c
// David Grund, Sep 10, 2021
// To calculate the slope of the |t| distribution of the measured data
// in the pt range 0.2 to 1.0 GeV/c

// cpp headers
#include <fstream> // print output to txt file
#include <iomanip> // std::setprecision()
// root headers
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TBrowser.h"
#include "TLegend.h"
#include "TStyle.h" // gStyle
#include "TString.h"
#include "TMath.h"
// roofit headers
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooGenericPdf.h"
#include "RooBinning.h"
// my headers
#include "TreesManager.h"

using namespace RooFit;

const Int_t nPtBins = 4;

Double_t t; // Mandelstam t
Double_t m; // inv mass
// Boundaries:
Double_t pt_i = 0.2;
Double_t pt_f = 1.0;
Double_t t_i = pt_i*pt_i;
Double_t t_f = pt_f*pt_f;
Double_t m_i = 3.0;
Double_t m_f = 3.2;

void CalcSlope();
void CalcBinning();
Int_t SumHistEntries(TH1D *h, Int_t i_low, Int_t i_upp);
Int_t FindBinIndex(TH1D *h, Double_t x);
void CalcBinning2();
void PrepareData();

void tSlope(){

    PrepareData();

    CalcSlope();

    CalcBinning();

    CalcBinning2();

    return;
}

void CalcSlope(){

    // PtBinning "Method 1"
    // In this way, background is not subtracted!
    // Better way is implemented in tSlopeWithoutBkg.c

    TFile *fFileIn = TFile::Open("Trees/tSlope/tSlope.root", "read");
    if(fFileIn) Printf("Input data loaded.");

    TTree *fTreeIn = dynamic_cast<TTree*> (fFileIn->Get("tSlope"));
    if(fTreeIn) Printf("Input tree loaded.");

    Int_t n_bins = 100;
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
    sprintf(fStrReduce,"t > %f && t < %f && m > %f && m < %f",t_i,t_f,m_i,m_f);
    RooRealVar t("t","t",t_i,t_f);
    RooRealVar m("m","m",m_i,m_f);
    RooRealVar slope("slope","slope",2.,1.,10.); // positive (extra minus in the exp)
    // Define RooFit datasets
    RooDataSet *fDataIn = new RooDataSet("fDataIn","fDataIn", RooArgSet(t,m), Import(*fTreeIn));
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
    // Set the canvas
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");

    cFit->SetTopMargin(0.055);
    cFit->SetBottomMargin(0.12);
    cFit->SetRightMargin(0.03);
    cFit->SetLeftMargin(0.11);    
    cFit->SetLogy();
    // Roo frame
    RooPlot* fFrame = t.frame(Title("Fit of |t| distribution")); 
    fDataSet->plotOn(fFrame,Name("fDataSet"),Binning(binT),MarkerStyle(20),MarkerSize(1.));
    ExpPdf.plotOn(fFrame,Name("ExpPdf"),Components(ExpPdf),LineColor(kBlue),LineStyle(kDashed),LineWidth(3));
    // Vertical axis
    fFrame->GetYaxis()->SetTitle(Form("Counts per %i MeV^{2}",bin_size));
    fFrame->GetYaxis()->SetTitleOffset(1.0);
    fFrame->GetYaxis()->SetTitleSize(0.05);
    fFrame->GetYaxis()->SetLabelSize(0.05);
    // Horizontal axis
    fFrame->GetXaxis()->SetTitle("|#it{t}| (GeV^{2})");
    fFrame->GetXaxis()->SetTitleSize(0.05);
    fFrame->GetXaxis()->SetLabelSize(0.05);
    // Finally plot it
    fFrame->Draw();
    // Legend
    TLegend *l = new TLegend(0.14,0.18,0.5,0.48);  
    l->AddEntry("ExpPdf","fit: #it{e}^{-#lambda t}","L");
    l->AddEntry((TObject*)0,Form("#lambda = %.3f GeV^{-2}", slope.getVal()),"");  
    l->AddEntry("fDataSet","data","P");
    l->AddEntry((TObject*)0,Form("#it{m} #in (%.1f, %.1f) GeV/#it{c}^{2}", m_i, m_f),""); 
    l->AddEntry((TObject*)0,Form("#it{p}_{T} #in (%.1f, %.1f) GeV/#it{c}", pt_i, pt_f),"");
    l->SetTextSize(0.048);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->Draw();
    // Save the result
    cFit->Print("PtBinning/data_slope.png");
    cFit->Print("PtBinning/data_slope.pdf");
    // Print the slope to the txt file
    TString name = "PtBinning/slope.txt";
    ofstream outfile (name.Data());
    outfile << slope.getVal();
    outfile.close();
    Printf("*** Results printed to %s.***", name.Data());

    return;
}

void CalcBinning(){

    Double_t b; // slope
    Double_t tBoundaries[nPtBins+1] = { 0 };
    Double_t ptBoundaries[nPtBins+1] = { 0 };
    tBoundaries[0] = t_i;
    ptBoundaries[0] = pt_i;

    // Load the calculated slope
    ifstream infile;
    infile.open("PtBinning/slope.txt");   
    while(!infile.eof()){
        infile >> b;
    }
    infile.close();
    Printf("Slope loaded: %.3f", b);

    // Initialise
    Double_t t_1 = t_i;
    
    // Fraction of the integral occupied by each bin
    Double_t f = 1.0/nPtBins;

    // Go over the bins
    for(Int_t i = 1; i <= nPtBins; i++){
        // Exponentials
        Double_t e_i = TMath::Exp(-b*t_i);
        Double_t e_f = TMath::Exp(-b*t_f);
        // Current t boundary
        Double_t e_1 = TMath::Exp(-b*t_1);

        // Integral in the full range
        Double_t e_int = e_i - e_f; // the exact integral has an extra 1/b
    
        // New value of t
        Double_t t_n = -TMath::Log(e_1-f*e_int)/b;

        tBoundaries[i] = t_n;
        ptBoundaries[i] = TMath::Sqrt(t_n);

        // Prepare for the next bin
        t_1 = t_n;
    }

    TFile *fFileIn = TFile::Open("Trees/tSlope/tSlope.root", "read");
    if(fFileIn) Printf("Input data loaded.");

    TTree *fTreeIn = dynamic_cast<TTree*> (fFileIn->Get("tSlope"));
    if(fTreeIn) Printf("Input tree loaded.");

    fTreeIn->SetBranchAddress("t", &t);
    fTreeIn->SetBranchAddress("m", &m);

    Int_t n_bins = 100;
    TH1D *hist = new TH1D("hist", "hist in Mandelstam t", n_bins, t_i, t_f); 

    for(Int_t iEntry = 0; iEntry < fTreeIn->GetEntries(); iEntry++){
        fTreeIn->GetEntry(iEntry);
        if(t > t_i && t < t_f && m > m_i && m < m_f) hist->Fill(t);
    }

    // Plot the histogram
    TCanvas *cHist2 = new TCanvas("cHist2", "cHist2", 900, 600);
    cHist2->SetLogy();

    hist->Draw();

    // Number of entries in the histogram
    Printf("# of entries in the histogram: %.0f", hist->GetEntries());
        
    // Calculate number of events per computed t bins
    Int_t nEvPerBin[nPtBins] = { 0 };
    Int_t nEvTotal = 0;
    for(Int_t iPtBin = 0; iPtBin < nPtBins; iPtBin++){
        // Find the bin ranges for given pt bin
        Int_t bin_i = -1;
        if(iPtBin == 0){
            bin_i = FindBinIndex(hist,tBoundaries[iPtBin]);
        } else {
            bin_i = FindBinIndex(hist,tBoundaries[iPtBin]) + 1;
        }
        Printf("bin_i = %i", bin_i);
        Int_t bin_f = FindBinIndex(hist,tBoundaries[iPtBin+1]);
        Printf("bin_f = %i", bin_f);
        // Sum the number of events in given pt bin
        nEvPerBin[iPtBin] = SumHistEntries(hist,bin_i,bin_f);
        Printf("*** Bin %i: %i", (iPtBin+1), nEvPerBin[iPtBin]);
        nEvTotal += nEvPerBin[iPtBin];
    }
    Printf("Total # of events = %i", nEvTotal);

    // Print the boundaries to the txt file
    TString name = "PtBinning/bin_boundaries1.txt";
    ofstream outfile (name.Data());
    outfile << std::fixed << std::setprecision(3) << "t:\t"; // Set the precision to 3 dec places
    for(Int_t i = 0; i <= nPtBins; i++){
        if(i != nPtBins){
            outfile <<  tBoundaries[i] << "\t";
        } else {
            outfile << tBoundaries[i] << std::endl;
        }
    }
    outfile << "pt:\t";
    for(Int_t i = 0; i <= nPtBins; i++){
        if(i != nPtBins){
            outfile << ptBoundaries[i] << "\t";
        } else {
            outfile << ptBoundaries[i] << std::endl;
        }
    }
    outfile << "nEv:\t";
    for(Int_t i = 0; i < nPtBins; i++){
        if(i != nPtBins){
            outfile << nEvPerBin[i] << "\t";
        } else {
            outfile << nEvPerBin[i] << std::endl;
        }
    }
    outfile.close();
    Printf("*** Results printed to %s.***", name.Data());

    return;
}

Int_t SumHistEntries(TH1D *h, Int_t i_low, Int_t i_upp){
    Int_t Int = 0;
    for(Int_t i = i_low; i <= i_upp; i++){
        Int += h->GetBinContent(i);
    }
    return Int;
}

Int_t FindBinIndex(TH1D *h, Double_t x){
    Int_t i = 1;
    while((h->GetBinLowEdge(i) + h->GetBinWidth(i)) < x) i++; 
    return i;
}

void CalcBinning2(){

    // PtBinning "Method 2"
    // Finer binning (n_bins = 200)
    // Simply add up the number of events until optimal value is reached
    // Then go to next bin
    // Note: in this method, background is also not subtracted

    TFile *fFileIn = TFile::Open("Trees/tSlope/tSlope.root", "read");
    if(fFileIn) Printf("Input data loaded.");

    TTree *fTreeIn = dynamic_cast<TTree*> (fFileIn->Get("tSlope"));
    if(fTreeIn) Printf("Input tree loaded.");

    fTreeIn->SetBranchAddress("t", &t);
    fTreeIn->SetBranchAddress("m", &m);

    Int_t n_bins = 200;
    TH1D *hist = new TH1D("hist", "hist in Mandelstam t", n_bins, t_i, t_f); 

    TCanvas *cHist3 = new TCanvas("cHist3", "cHist3", 900, 600);
    cHist3->SetLogy();
    hist->Draw();

    for(Int_t iEntry = 0; iEntry < fTreeIn->GetEntries(); iEntry++){
        fTreeIn->GetEntry(iEntry);
        if(t > t_i && t < t_f && m > m_i && m < m_f) hist->Fill(t);
    }
    // Number of entries in the histogram
    Printf("# of entries in the histogram: %.0f", hist->GetEntries()); 

    Double_t nOptEvPerBin = hist->GetEntries() / nPtBins;
    Printf("Optimal # of entries per bin: %.0f", nOptEvPerBin); 

    Double_t tBoundaries[nPtBins+1];
    tBoundaries[0] = hist->GetBinLowEdge(1);

    // Go over bins in hist and find pt bins
    Int_t nEvPerBin[nPtBins] = { 0 };
    Int_t iPtBin = 0;
    Int_t iHistBin = 1;
    Int_t SumEv = 0;
    while(kTRUE){
        SumEv += hist->GetBinContent(iHistBin);
        if(SumEv >= nOptEvPerBin){
            tBoundaries[iPtBin+1] = hist->GetBinLowEdge(iHistBin+1);
            nEvPerBin[iPtBin] = SumEv;
            iPtBin++;
            SumEv = 0;
        }
        if(iHistBin == hist->GetNbinsX()){
            tBoundaries[iPtBin+1] = hist->GetBinLowEdge(iHistBin) + hist->GetBinWidth(iHistBin);
            nEvPerBin[iPtBin] = SumEv;
            break;
        } 
        iHistBin++;
    }
    // Print the results to the console
    Int_t nEvTotal = 0;
    for(Int_t iPtBin = 0; iPtBin < nPtBins; iPtBin++){
        Printf("Bin %i: t_low = %.3f, t_upp = %.3f, entries: %i", iPtBin+1, tBoundaries[iPtBin], tBoundaries[iPtBin+1], nEvPerBin[iPtBin]);
        nEvTotal += nEvPerBin[iPtBin];
    }
    Printf("Total # of events = %i", nEvTotal);
    // Print the boundaries to the txt file
    TString name = "PtBinning/bin_boundaries2.txt";
    ofstream outfile (name.Data());
    outfile << std::fixed << std::setprecision(3) << "t:\t"; // Set the precision to 3 dec places
    for(Int_t i = 0; i <= nPtBins; i++){
        if(i != nPtBins){
            outfile <<  tBoundaries[i] << "\t";
        } else {
            outfile << tBoundaries[i] << std::endl;
        }
    }
    outfile << "pt:\t";
    for(Int_t i = 0; i <= nPtBins; i++){
        if(i != nPtBins){
            outfile << TMath::Sqrt(tBoundaries[i]) << "\t";
        } else {
            outfile << TMath::Sqrt(tBoundaries[i]) << std::endl;
        }
    }
    outfile << "nEv:\t";
    for(Int_t i = 0; i < nPtBins; i++){
        if(i != nPtBins){
            outfile << nEvPerBin[i] << "\t";
        } else {
            outfile << nEvPerBin[i] << std::endl;
        }
    }
    outfile.close();
    Printf("*** Results printed to %s.***", name.Data());

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
    tSlope->Branch("m", &m, "m/D"); 

    Int_t nEntriesAnalysed = 0;

    for(Int_t iEntry = 0; iEntry < fTreeIn->GetEntries(); iEntry++){
        fTreeIn->GetEntry(iEntry);
        if(EventPassed(-1,-1) && fPt > 0.0 && fPt < 2.0){
            t = fPt*fPt;
            m = fM;
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