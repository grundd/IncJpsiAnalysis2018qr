// ReadAnalysisResults.c

// root headers
#include "TFile.h"
#include "TList.h"
#include "TCanvas.h"
#include "TH2.h"

void ReadAnalysisResults()
{
    TFile *f = TFile::Open("AnalysisResults.root", "read");
    if(f) Printf("Input file loaded.");

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

    return;
}