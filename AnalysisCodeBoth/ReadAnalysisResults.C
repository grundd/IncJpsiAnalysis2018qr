// ReadAnalysisResults.c

// root headers
#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TH1.h"

Double_t fMGen;

void ReadAnalysisResults()
{
    TFile *f = TFile::Open("AnalysisResults_loc_MC01.root", "read");
    if(f) Printf("Input file loaded.");

    // dE/dx histograms
    TList *l = dynamic_cast<TList*> (f->Get("AnalysisOutput/fOutputList"));
    if(l) Printf("Input list loaded.");

    TH2F *hTPCdEdx = (TH2F*)l->FindObject("hTPCdEdx");
    if(hTPCdEdx) Printf("Input hist loaded.");

    TH2F *hTPCdEdxMuon = (TH2F*)l->FindObject("hTPCdEdxMuon");
    if(hTPCdEdxMuon) Printf("Input hist loaded.");

    TH2F *hTPCdEdxElectron = (TH2F*)l->FindObject("hTPCdEdxElectron");
    if(hTPCdEdxElectron) Printf("Input hist loaded.");

    TCanvas *c1 = new TCanvas("c1","c1",900,600);
    hTPCdEdx->Draw("COLZ");
    TCanvas *c2 = new TCanvas("c2","c2",900,600);
    hTPCdEdxMuon->Draw("COLZ");
    TCanvas *c3 = new TCanvas("c3","c3",900,600);
    hTPCdEdxElectron->Draw("COLZ");

    // histogram with fMGen
    TTree *tGen = dynamic_cast<TTree*> (f->Get("AnalysisOutput/fTreeJpsiMCGen"));
    if(tGen) Printf("Input tree loaded.");
    tGen->SetBranchAddress("fMGen", &fMGen);
    TH1D *hMGen = new TH1D("hMGen","hMGen",100,3.0965,3.0975); // around J/psi mass peak

    Printf("Tree has %lli entries.", tGen->GetEntries());

    for(Int_t iEntry = 0; iEntry < tGen->GetEntries(); iEntry++){
        tGen->GetEntry(iEntry);
        hMGen->Fill(fMGen);
    }

    TCanvas *c4 = new TCanvas("c4","c4",900,600);
    hMGen->Draw();

    return;
}