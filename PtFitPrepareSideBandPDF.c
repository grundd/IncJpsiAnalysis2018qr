// PtFitPrepareSideBandPDF.c
// David Grund, Sep 29, 2021

// my headers
#include "AnalysisManager.h"
#include "PtFitUtilities.h"

void PrepareTree();
void PreparePDF();
void PlotPDF();
void CompareWithBkgFromInvMassFits();

void PtFitPrepareSideBandPDF(){

    //PrepareTree();

    //PreparePDF();

    //PlotPDF();

    CompareWithBkgFromInvMassFits();

    return;
}

void CompareWithBkgFromInvMassFits(){

    // Fill the SideBand histogram
    TFile *file = TFile::Open("Trees/PtFit/SideBandPDF/TreeData.root","read");
    if(file) Printf("Input file %s loaded.", file->GetName());     

    TTree *tree = dynamic_cast<TTree*> (file->Get("fDataTree"));
    if(tree) Printf("Input tree loaded.");

    tree->SetBranchAddress("fPt",&fPt);

    SetPtBins(3);

    TH1D *hSideBand = new TH1D("hSideBand", "hSideBand", nBins, ptEdges);
    for(Int_t iEntry = 0; iEntry < tree->GetEntries(); iEntry++){
        tree->GetEntry(iEntry);
        hSideBand->Fill(fPt);
    }
    // Normalize the histogram to bin widths
    hSideBand->Scale(1./hSideBand->Integral(),"width");

    // Load the bkgr histogram from the invariant mass fits
    file = TFile::Open("Trees/PtFit/JpsiSignalNoBkg_Binning2.root","read");
    if(file) Printf("Input file %s loaded.", file->GetName()); 

    TList *list = (TList*) file->Get("HistList");
    if(list) Printf("List %s loaded.", list->GetName()); 

    TH1D *hNBkgrBins = (TH1D*)list->FindObject("hNBkgrBins");
    // How many background events are there with 0.2 < pT < 1.0 GeV/c?
    Int_t bin1 = 1;
    Double_t PtEdgeLow = hNBkgrBins->GetBinLowEdge(bin1);
    while(PtEdgeLow < 0.2){
        bin1++;
        PtEdgeLow = hNBkgrBins->GetBinLowEdge(bin1);
    }
    Printf("Found low bin edge of %.3f GeV for bin number %i.", hNBkgrBins->GetBinLowEdge(bin1), bin1);
    Int_t bin2 = 1;
    PtEdgeLow = hNBkgrBins->GetBinLowEdge(bin2);
    while(PtEdgeLow < 1.0){
        bin2++;
        PtEdgeLow = hNBkgrBins->GetBinLowEdge(bin2);
    }
    Printf("Found low bin edge of %.3f GeV for bin number %i.", hNBkgrBins->GetBinLowEdge(bin2), bin2);

    Double_t nBkgEv = 0;
    for(Int_t iBin = bin1; iBin < bin2; iBin++){
        Printf("Current low bin edge: %.3f, bin width: %.3f", hNBkgrBins->GetBinLowEdge(iBin), hNBkgrBins->GetBinWidth(iBin));
        Double_t add;
        if(iBin < bin2-1){
            add = hNBkgrBins->GetBinContent(iBin);
        } else add = hNBkgrBins->GetBinContent(iBin) / 2;
        Printf("Bin content: %.3f, adding %.3f to nBkgEv.", hNBkgrBins->GetBinContent(iBin), add);
        nBkgEv += add;
    }
    Printf("There are %.2f bkgr events with 0.2 < pT < 1.0 GeV/c", nBkgEv);

    hNBkgrBins->Scale(1./hNBkgrBins->Integral(),"width");

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");

    // Plot the result (normalized to bin widths)
    TCanvas *c = new TCanvas("c","c",900,600);
    c->cd();
    // Canvas settings
    c->SetLogy();
    c->SetTopMargin(0.02);
    c->SetBottomMargin(0.14);
    c->SetRightMargin(0.03);
    c->SetLeftMargin(0.1);
    // Vertical axis
    hSideBand->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (normalized to bin width)");
    hSideBand->GetYaxis()->SetTitleSize(0.05);
    hSideBand->GetYaxis()->SetTitleOffset(0.9);
    hSideBand->GetYaxis()->SetLabelSize(0.05);
    // Horizontal axis
    hSideBand->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hSideBand->GetXaxis()->SetTitleSize(0.05);
    hSideBand->GetXaxis()->SetTitleOffset(1.2);
    hSideBand->GetXaxis()->SetLabelSize(0.05);
    hSideBand->GetXaxis()->SetDecimals(1);    
    // Draw it
    hSideBand->SetMarkerStyle(21);
    hSideBand->SetMarkerColor(kBlue);
    hSideBand->SetMarkerSize(1.0);
    hSideBand->SetLineColor(kBlack);
    hSideBand->SetLineWidth(1.0);
    hSideBand->Draw("P E1");  
    // Draw the other histogram
    hNBkgrBins->SetMarkerStyle(21);
    hNBkgrBins->SetMarkerColor(kRed);
    hNBkgrBins->SetMarkerSize(1.0);
    hNBkgrBins->SetLineColor(kBlack);
    hNBkgrBins->SetLineWidth(1.0);      
    hNBkgrBins->Draw("P E1 SAME");  

    c->Print("Results/PtFitWithoutBkg/BkgComparison/bkg_comparison.png");
    c->Print("Results/PtFitWithoutBkg/BkgComparison/bkg_comparison.pdf");

    return;
}

void PreparePDF(){

    TFile *file = TFile::Open("Trees/PtFit/SideBandPDF/TreeData.root","read");
    if(file) Printf("Input file %s loaded.", file->GetName());     

    TTree *tree = dynamic_cast<TTree*> (file->Get("fDataTree"));
    if(tree) Printf("Input tree loaded.");

    SetPtBins(3);

    TH1D *hSideBand = new TH1D("hSideBand", "hSideBand", nBins, ptEdges);
    tree->Draw("fPt >> hSideBand");

    TList *l = new TList();
    l->Add(hSideBand);

    TFile *f = new TFile(Form("%sSideBandPDF/PDF_SideBand.root", OutputPDFs.Data()),"RECREATE");
    l->Write("HistList", TObject::kSingleKey);
    f->ls();

    return;
}

void PlotPDF(){

    TFile *file = TFile::Open("Trees/PtFit/SideBandPDF/TreeData.root","read");
    if(file) Printf("Input file %s loaded.", file->GetName());     

    TTree *tree = dynamic_cast<TTree*> (file->Get("fDataTree"));
    if(tree) Printf("Input tree loaded.");

    tree->SetBranchAddress("fPt",&fPt);

    SetPtBins(3);

    // With variable binning
    TH1D *hSideBand = new TH1D("hSideBand", "hSideBand", nBins, ptEdges);
    // With uniform binning
    TH1D *hSideBand2 = new TH1D("hSideBand2", "hSideBand2", 150, 0.0, 2.0);
    // Fill both histograms
    for(Int_t iEntry = 0; iEntry < tree->GetEntries(); iEntry++){
        tree->GetEntry(iEntry);
        hSideBand->Fill(fPt);
        hSideBand2->Fill(fPt);
    }

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetPaintTextFormat("4.2f");

    // Plot the result (normalized to bin widths)
    TCanvas *c1 = new TCanvas("c1","c1",900,600);
    c1->cd();
    // Canvas settings
    c1->SetLogy();
    c1->SetTopMargin(0.02);
    c1->SetBottomMargin(0.14);
    c1->SetRightMargin(0.03);
    c1->SetLeftMargin(0.1);
    // Normalize the histogram to bin widths
    hSideBand->Scale(1.,"width");
    // Vertical axis
    hSideBand->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T} (normalized to bin width)");
    hSideBand->GetYaxis()->SetTitleSize(0.05);
    hSideBand->GetYaxis()->SetTitleOffset(0.9);
    hSideBand->GetYaxis()->SetLabelSize(0.05);
    // Horizontal axis
    hSideBand->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hSideBand->GetXaxis()->SetTitleSize(0.05);
    hSideBand->GetXaxis()->SetTitleOffset(1.2);
    hSideBand->GetXaxis()->SetLabelSize(0.05);
    hSideBand->GetXaxis()->SetDecimals(1);    
    // Draw it
    hSideBand->SetMarkerStyle(21);
    hSideBand->SetMarkerColor(kBlue);
    hSideBand->SetMarkerSize(1.0);
    hSideBand->SetLineColor(kBlack);
    hSideBand->SetLineWidth(1.0);
    hSideBand->Draw("P E1");

    c1->Print("Results/PtFit/SideBandPDF/hSideBand.png");
    c1->Print("Results/PtFit/SideBandPDF/hSideBand.pdf");

    // Plot the result with uniform binning
    TCanvas *c2 = new TCanvas("c2","c2",900,600);
    c2->cd();
    // Canvas settings
    c2->SetLogy();
    c2->SetTopMargin(0.02);
    c2->SetBottomMargin(0.14);
    c2->SetRightMargin(0.03);
    c2->SetLeftMargin(0.1);
    // Normalize the histogram to bin widths
    //hSideBand->Scale(1.,"width");
    // Vertical axis
    hSideBand2->GetYaxis()->SetTitle("Counts per bin");
    hSideBand2->GetYaxis()->SetTitleSize(0.05);
    hSideBand2->GetYaxis()->SetTitleOffset(0.9);
    hSideBand2->GetYaxis()->SetLabelSize(0.05);
    // Horizontal axis
    hSideBand2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hSideBand2->GetXaxis()->SetTitleSize(0.05);
    hSideBand2->GetXaxis()->SetTitleOffset(1.2);
    hSideBand2->GetXaxis()->SetLabelSize(0.05);
    hSideBand2->GetXaxis()->SetDecimals(1);    
    // Draw it
    hSideBand2->SetMarkerStyle(21);
    hSideBand2->SetMarkerColor(kBlue);
    hSideBand2->SetMarkerSize(1.0);
    hSideBand2->SetLineColor(kBlack);
    hSideBand2->SetLineWidth(1.0);
    hSideBand2->Draw("P E1");

    c2->Print("Results/PtFit/SideBandPDF/hSideBand2.png");
    c2->Print("Results/PtFit/SideBandPDF/hSideBand2.pdf");

    return;
}

void PrepareTree(){

    TFile *file = TFile::Open("Trees/AnalysisData/AnalysisResultsLHC18qrMerged.root", "read");
    if(file) Printf("File %s loaded.", file->GetName());

    TTree *fTreeIn = dynamic_cast<TTree*> (file->Get("AnalysisOutput/fTreeJPsi"));
    if(fTreeIn) Printf("Input tree loaded.");

    ConnectTreeVariables(fTreeIn);

    Printf("Tree %s has %lli entries.", fTreeIn->GetName(), fTreeIn->GetEntries());

    TFile *f = new TFile("Trees/PtFit/SideBandPDF/TreeData.root","RECREATE");

    TTree *TreeData = new TTree("fDataTree", "fDataTree");
    TreeData->Branch("fPt", &fPt, "fPt/D");

    // Loop over tree entries
    Int_t nEntriesAnalysed = 0;
    Int_t nEvPassed = 0;

    for(Int_t iEntry = 0; iEntry < fTreeIn->GetEntries(); iEntry++){
        fTreeIn->GetEntry(iEntry);

        // iMassCut == 0 => 2.2 < m < 4.5 GeV
        // iPtCut   == 2 => pt < 2.0 GeV
        if(EventPassed(0, 2) && (fM > 3.3 && fM < 4.5)){ // 3.3 < m < 4.5 GeV
            nEvPassed++;
            TreeData->Fill();
        } 

        if((iEntry+1) % 100000 == 0){
            nEntriesAnalysed += 100000;
            Printf("%i entries analysed.", nEntriesAnalysed);
        }
    }
    Printf("Done.");
    Printf("%i events passed the selections.", nEvPassed);    

    f->Write("",TObject::kWriteDelete);

    return;
}