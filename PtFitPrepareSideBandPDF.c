// PtFitPrepareSideBandPDF.c
// David Grund, Sep 29, 2021

// my headers
#include "AnalysisManager.h"
#include "PtFitUtilities.h"

void PrepareTree();

void PtFitPrepareSideBandPDF(){

    //PrepareTree();

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

    // Plot the result (normalized to bin widths)
    TCanvas *c = new TCanvas("c","c",900,600);
    c->SetLogy();
    hSideBand->Scale(1.,"width");
    hSideBand->Draw();

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